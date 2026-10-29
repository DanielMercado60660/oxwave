// Copyright (c) 2008, Colorado School of Mines,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or
// without modification, are permitted provided that the following
// conditions are met:
//
//     *  Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     *  Redistributions in binary form must reproduce the above
//        copyright notice, this list of conditions and the following
//        disclaimer in the documentation and/or other materials provided
//        with the distribution.
//     *  Neither the name of the Colorado School of Mines nor the names of
//        its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
//
// Warranty Disclaimer:
// THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Rust port for OxWave, a Rust port of CWP/SU.

//! Byte metadata reader from `su/lib/readkfile.c`, used by `su/lib/gridread.c`.
//! Credits: Copyright (c) Colorado School of Mines, 2021. All rights reserved.
//! The C source gives no named author.

use oxwave_core::cio::{fgetc, fseeko, CFile, EOF, SEEK_SET};
use oxwave_k_file::getkvalscsv;

// fgets without the core wrapper's lossy UTF-8 conversion. Retain stale bytes
// after its terminating NUL, as C's fixed buffers do.
fn read_line(fp: &CFile, raw: &mut [u8; 10001]) -> bool {
    let mut n = 0;
    while n < 10000 {
        let c = fgetc(fp);
        if c == EOF {
            break;
        }
        raw[n] = c as u8;
        n += 1;
        if c == b'\n' as i32 {
            break;
        }
    }
    if n == 0 {
        return false;
    }
    raw[n] = 0;
    true
}

// Byte equivalent of tparse: C metadata has no encoding requirement.
fn parse_into(text: &[u8], fields: &mut [Vec<u8>]) -> usize {
    let count = text.split(|&b| b == b',').count();
    if count > fields.len() {
        // C would overrun gridread's fixed 999-element output arrays.
        oxwave_core::err!("gridread: too many K-file metadata fields.");
    }
    for (dst, field) in fields.iter_mut().zip(text.split(|&b| b == b',')) {
        *dst = if field.is_empty() { b"null".to_vec() } else { field.to_vec() };
    }
    count
}

// Private byte-preserving readkfile path for gridread. Keep the metadata
// control flow and error precedence from C; reuse the shared K-file crate's
// numeric record parser. This can be removed when that crate accepts bytes.
pub(super) fn read_grid_kfile(
    fp_r: &CFile,
    names: &mut [Vec<u8>],
    forms: &mut [Vec<u8>],
    dfield: &mut [f64],
    numcasesout: &mut i32,
    errwarn: &mut i32,
) {
    *errwarn = 0;
    let mut textraw = [0u8; 10001];
    let mut textbeg = [0u8; 10001];
    let mut num_names = 0;
    let mut num_forms = 0;
    let mut num_c_su_names = 0;
    let mut num_c_su_forms = 0;
    let mut num_c_su_setid = 0;
    let mut read_names = 0;
    let mut read_forms = 0;
    while read_line(fp_r, &mut textraw) {
        // C consumes another line before noticing both lists are complete.
        if read_names == -1 && read_forms == -1 {
            break;
        }
        let mut tsize = 0;
        for &c in &textraw {
            if c == 0 || c == b'\n' || c == b'\r' {
                break;
            }
            if c != b' ' && c != b'\t' {
                textbeg[tsize] = c;
                tsize += 1;
            }
        }
        // No new NUL until a list is parsed: short lines see stale bytes.
        for b in &mut textbeg[..10001] {
            *b = b.to_ascii_lowercase();
        }
        if textbeg.starts_with(b"c_su_setid") {
            num_c_su_setid += 1;
        }
        if textbeg.starts_with(b"c_su_names") {
            num_c_su_names += 1;
        }
        if textbeg.starts_with(b"c_su_forms") {
            num_c_su_forms += 1;
        }
        if read_names > 0 {
            textbeg[tsize] = 0;
            num_names = parse_into(&textbeg[..tsize], names);
            read_names = -1;
        }
        if textbeg.starts_with(b"c_su_names") {
            read_names = 1;
        }
        if read_forms > 0 {
            textbeg[tsize] = 0;
            num_forms = parse_into(&textbeg[..tsize], forms);
            read_forms = -1;
        }
        if textbeg.starts_with(b"c_su_forms") {
            read_forms = 1;
        }
    }
    fseeko(fp_r, 0, SEEK_SET);
    for (count, duplicate, missing) in [
        (num_c_su_names, 1, 2),
        (num_c_su_forms, 3, 4),
        (num_c_su_setid, 5, 6),
    ] {
        if count > 1 {
            *errwarn = duplicate;
            return;
        } else if count == 0 {
            *errwarn = missing;
            return;
        }
    }
    if num_forms != num_names {
        *errwarn = 7;
        return;
    }
    if num_names == 0 {
        // C would inspect an uninitialized names[0].
        oxwave_core::err!("gridread: missing K-file metadata lists.");
    }
    if names[0].starts_with(b"c_su_id") {
        names[0] = b"null".to_vec();
    }
    for n in 0..num_names {
        if !names[n].starts_with(b"null") {
            for m in n + 1..num_names {
                if names[n] == names[m] {
                    *errwarn = 9;
                    return;
                }
            }
        }
    }
    let mut nspot = Vec::with_capacity(num_names);
    let mut numcases = 0;
    for i in 0..num_names {
        if !names[i].starts_with(b"null") {
            names[numcases] = names[i].clone();
            forms[numcases] = forms[i].clone();
            nspot.push(i as i32);
            numcases += 1;
        }
    }
    if numcases > dfield.len() {
        oxwave_core::err!("gridread: too many K-file numeric fields.");
    }
    textraw.fill(0);
    textbeg.fill(0);
    let (mut ncount, mut comerr, mut morerr, mut numerr, mut nblank) = (0, 0, 0, 0, 0);
    let mut nextrow = false;
    while read_line(fp_r, &mut textraw) {
        ncount += 1;
        let mut textfront = [0u8; 10];
        for n in 0..10 {
            textfront[n] = textraw[n].to_ascii_lowercase();
        }
        if textfront.starts_with(b"c_su") || nextrow {
            nextrow = textfront.starts_with(b"c_su_names") || textfront.starts_with(b"c_su_forms");
        } else if textraw[0] == b'K' {
            getkvalscsv(
                &textraw,
                &mut textbeg,
                b',',
                &mut dfield[..numcases],
                &nspot,
                ncount,
                &mut comerr,
                &mut morerr,
                &mut numerr,
                &mut nblank,
            );
            break;
        }
    }
    for (count, error) in [(numerr, 10), (morerr, 11), (comerr, 12)] {
        if count > 0 {
            *errwarn = error;
            return;
        }
    }
    if nblank > 0 {
        *errwarn = -1;
    }
    *numcasesout = numcases as i32;
}

