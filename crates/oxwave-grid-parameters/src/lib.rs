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

//! Port of `su/lib/gridread.c` and `su/include/gridread.h`.
//! Credits: Copyright (c) Colorado School of Mines, 2021. All rights reserved;
//! the C source gives no named author.

use oxwave_core::cio::CFile;
use oxwave_core::par::getpardouble;

mod kfile;
use kfile::read_grid_kfile;

const NAMES: [&str; 18] = [
    "", "grid_lf", "grid_xa", "grid_ya", "grid_xb", "grid_yb",
    "grid_xc", "grid_yc", "grid_xd", "grid_yd", "grid_wb", "grid_wc",
    "grid_nb", "grid_nc", "grid_fp", "grid_lp", "grid_sb", "grid_cb",
];

/// C: `void gridcommand(int *enough);`
///
/// Reports whether all (1), none (0), or only some (-1) of the eight
/// required grid parameters are present in the current command line.
pub fn gridcommand(enough: &mut i32) {
    let mut nfound = 0;
    for i in 2..12 {
        if i == 8 || i == 9 { continue; }
        let mut gval = 0.0;
        if getpardouble(NAMES[i], &mut gval) != 0 { nfound += 1; }
    }
    *enough = 0;
    if nfound == 8 { *enough = 1; }
    else if nfound > 0 { *enough = -1; }
}

/// C: `void gridread(FILE *fpR, double *gvals, int *errwarn2);`
///
/// Reads grid parameters from an optional K-file and the current command
/// line, which takes precedence. `gvals` must have at least 20 elements.
/// Set `errwarn2` to zero to suppress diagnostics; on return positive values
/// report errors and -1 reports blank K-file fields.
/// Metadata is compared as bytes, including field names and formats that
/// are not UTF-8, just as in C.
pub fn gridread(fp_r: Option<&CFile>, gvals: &mut [f64], errwarn2: &mut i32) {
    if gvals.len() < 20 {
        oxwave_core::err!("gridread: gvals must contain at least 20 elements.");
    }
    let iwarn = *errwarn2;
    *errwarn2 = 0;
    let mut names = vec![Vec::new(); 999];
    let mut forms = vec![Vec::new(); 999];
    let mut dfield = vec![0.0; 999];
    let mut numcases = 0;
    let mut errwarn = 0;
    if let Some(fp) = fp_r {
        read_grid_kfile(fp, &mut names, &mut forms, &mut dfield, &mut numcases, &mut errwarn);
    }
    if iwarn != 0 {
        match errwarn {
            1 => oxwave_core::warn!("K-file read error: more than one C_SU_NAMES parameter record."),
            2 => oxwave_core::warn!("K-file read error: no C_SU_NAMES parameter record."),
            3 => oxwave_core::warn!("K-file read error: more than one C_SU_FORMS parameter record."),
            4 => oxwave_core::warn!("K-file read error: no C_SU_FORMS parameter record."),
            5 => oxwave_core::warn!("K-file read error: more than one C_SU_SETID record."),
            6 => oxwave_core::warn!("K-file read error: no C_SU_SETID record."),
            7 => oxwave_core::warn!("K-file read error: different number of values on C_SU_NAMES and C_SU_FORMS."),
            8 => oxwave_core::warn!("K-file read error: unable to allocate memory."),
            9 => oxwave_core::warn!("K-file read error: name exists at least twice in C_SU_NAMES list."),
            10 => oxwave_core::warn!("K-file read error: at least 1 field-unreadable as a number."),
            11 => oxwave_core::warn!("K-file read error: at least 1 field containing 2 numbers."),
            12 => oxwave_core::warn!("K-file read error: not-enough-commas to get all values for C_SU_NAMES list."),
            n if n > 0 => oxwave_core::warn!("K-file read error: returned with some unrecognized error code."),
            -1 => oxwave_core::warn!("K-file read warning: at least 1 all-blank field, assumed zero for all."),
            _ => {}
        }
    }
    *errwarn2 = errwarn;
    if errwarn > 0 { return; }
    for name in &mut names[..numcases as usize] { name.make_ascii_lowercase(); }
    for val in &mut gvals[1..18] { *val = -1.1e308; }
    for i in 2..12 {
        if i == 8 || i == 9 { continue; }
        if getpardouble(NAMES[i], &mut gvals[i]) == 0 {
            for j in 0..numcases as usize {
                if names[j] == NAMES[i].as_bytes() { gvals[i] = dfield[j]; }
            }
            if gvals[i] < -1.0e308 {
                gvals[i] = (i + 100) as f64;
                if iwarn != 0 { oxwave_core::warn!("grid parameter error: %s not found.", NAMES[i]); }
                *errwarn2 = 20;
                return;
            }
        }
    }
}
