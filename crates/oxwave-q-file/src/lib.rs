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

//! Port of `su/lib/qdefine.c`: q-file and command parameter readers.
//! Credits: the original source dates the q-file standard to December 2021
//! and refers to subinqcsv for examples; it gives no named author credit.

use oxwave_core::cio::{fgetc, fgets, fseeko, CFile, EOF, SEEK_SET};
use oxwave_core::par::{countnparval, countparname, countparval, getnpardouble, getpardouble};

/// Owned C `struct QInfo`. Readers leave kinf for the caller.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct QInfo {
    pub kinf: Vec<i32>,
    pub nto: i32,
    /// Includes zero-initialized spare space requested through numdind.
    pub dlots: Vec<f64>,
}

/// High-level reader outputs. Errors return the C errwarn code.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct QFile {
    pub names: Vec<String>,
    /// C's reported count (includes an unused slot in the no-tuple case).
    pub numpname: i32,
    pub iztuple: i32,
    pub ktuple: i32,
    pub ifixd: i32,
    pub records: Vec<QInfo>,
    pub pindepa: Vec<f64>,
    pub ndims: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QStats {
    pub num_names: i32,
    pub num_ndims: i32,
    pub ibsize: i32,
    pub ivsize: i32,
    pub irtotal: i32,
    pub errwarn: i32,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct QHead {
    pub names: Vec<String>,
    pub ndims: Vec<String>,
    pub errwarn: i32,
}

fn rewind(fp: &CFile) {
    fseeko(fp, 0, SEEK_SET);
}
fn prefix(s: &[u8], p: &[u8]) -> bool {
    s.get(..p.len()).is_some_and(|b| b.eq_ignore_ascii_case(p))
}

// Like fgets, retain the original byte count before converting header fields
// to strings. Core's lossy fgets can expand invalid UTF-8 beyond n-1 bytes.
fn fgets_bytes(n: usize, fp: &CFile) -> Option<Vec<u8>> {
    let mut bytes = Vec::new();
    while bytes.len() + 1 < n {
        let c = fgetc(fp);
        if c == EOF {
            break;
        }
        bytes.push(c as u8);
        if c == b'\n' as i32 {
            break;
        }
    }
    if bytes.is_empty() {
        None
    } else {
        Some(bytes)
    }
}

// Use the frozen core's numeric-prefix conversion for C's %lf input item.
// A zero consumed length is a matching failure, not a numeric zero.
fn scan_double(s: &str) -> Option<f64> {
    let (value, consumed) = oxwave_core::par::strtod(s);
    (consumed > 0).then_some(value)
}

/// C: `void tparseq(char *tbuf, char d, char **fields, int *numfields);`
/// Empty fields become "null". Stops at NUL; delimiter must be ASCII.
/// Strings are owned and terminated safely, unlike C's unterminated allocations.
pub fn tparseq(tbuf: &str, d: u8) -> Vec<String> {
    assert!(d.is_ascii());
    let b = tbuf.as_bytes();
    let mut fields = Vec::new();
    let mut nbeg = 0;
    for n in 0..=b.len() {
        let c = b.get(n).copied().unwrap_or(0);
        if c == d || c == 0 {
            fields.push(if n > nbeg {
                tbuf[nbeg..n].to_owned()
            } else {
                "null".to_owned()
            });
            nbeg = n + 1;
        }
        if c == 0 {
            break;
        }
    }
    fields
}

/// C: `void getqvalscsv(char *textraw, char *textbeg, int maxtext, char rdel,
/// int *nspot, int numcases, double *dfield, int *izlast, int *notblank, int *errwarn);`
/// Scratch storage is internal. nspot is sorted and nonempty; dfield is at
/// least as long. Returns errwarn, preserving partial writes on errors.
/// textraw replaces the C buffer and maxtext, with an implicit final NUL.
pub fn getqvalscsv(
    textraw: &str,
    rdel: u8,
    nspot: &[i32],
    dfield: &mut [f64],
    izlast: &mut i32,
    notblank: &mut i32,
) -> i32 {
    assert!(!nspot.is_empty() && dfield.len() >= nspot.len());
    let b = textraw.as_bytes();
    let mut nbeg = 0;
    let mut nfield = 0;
    let mut ineed = 0;
    *notblank = -1;
    for n in 0..=b.len() {
        let c = b.get(n).copied().unwrap_or(0);
        if c == rdel || c == 0 || c == b'\n' || c == b'\r' {
            if nfield == nspot[ineed] {
                let text = &textraw[nbeg..n];
                let mut nb = false;
                let mut ib = false;
                for ch in text.bytes() {
                    if ch != b' ' {
                        nb = true;
                        if ib {
                            return 1;
                        }
                    }
                    if ch == b' ' && nb {
                        ib = true;
                    }
                }
                let dval = if nb {
                    let Some(v) = scan_double(text) else {
                        return 2;
                    };
                    *notblank = ineed as i32 + 1;
                    v
                } else {
                    0.0
                };
                dfield[ineed] = dval;
                ineed += 1;
                if ineed >= nspot.len() {
                    break;
                }
            }
            if c == 0 || c == b'\n' || c == b'\r' {
                break;
            }
            nbeg = n + 1;
            nfield += 1;
        }
    }
    *izlast = ineed as i32;
    0
}

/// C: `void getqinfo(FILE *fpR, int *num_names, int *num_ndims, int iqflag,
/// int *ibsize, int *ivsize, int *irtotal, int *errwarn);`
/// Rewinds the file, except for early error 3. Only iqflag == 1 collects sizes.
/// Records must have line terminators: otherwise C scans uninitialized bytes.
pub fn getqinfo(fp: &CFile, iqflag: i32) -> QStats {
    let mut out = QStats {
        num_names: 1,
        num_ndims: 1,
        ibsize: 0,
        ivsize: 0,
        irtotal: 0,
        errwarn: 0,
    };
    let mut num_names = 0;
    let mut num_ndims = 0;
    rewind(fp);
    let mut nextone = false;
    while let Some(s) = fgets_bytes(10001, fp) {
        if nextone {
            if !prefix(&s, b"c_su_id") {
                out.errwarn = 3;
                return out;
            }
            for &c in &s {
                if c == b'\n' || c == b'\r' {
                    break;
                }
                if c == b',' {
                    out.num_names += 1;
                }
            }
        }
        nextone = false;
        if prefix(&s, b"c_su_names") {
            num_names += 1;
            nextone = true;
        }
        if prefix(&s, b"q,") {
            break;
        }
    }
    rewind(fp);
    while let Some(s) = fgets_bytes(10001, fp) {
        if prefix(&s, b"c_su_ndims") {
            num_ndims += 1;
            for &c in &s {
                if c == b'\n' || c == b'\r' {
                    break;
                }
                if c == b',' {
                    out.num_ndims += 1;
                }
            }
        }
        if prefix(&s, b"q,") {
            break;
        }
    }
    rewind(fp);
    if iqflag == 1 {
        let mut isiz = 0;
        let mut icom = 0;
        while let Some(s) = fgets_bytes(10001, fp) {
            let mut idone = None;
            for (n, &c) in s.iter().enumerate() {
                if c == b',' {
                    icom += 1;
                }
                if c == b'\n' || c == b'\r' {
                    out.irtotal += 1;
                    idone = Some(n as i32 + 1);
                    break;
                }
            }
            if let Some(n) = idone {
                if isiz + n > out.ibsize {
                    out.ibsize = isiz + n;
                }
                isiz = 0;
                if icom > out.ivsize {
                    out.ivsize = icom;
                }
                icom = 0;
            } else {
                isiz += 10001;
            } // C counts the NUL in full fgets chunks.
        }
        out.ivsize += 1;
        rewind(fp);
    }
    if num_names > 1 {
        out.errwarn = 1;
    } else if num_ndims > 1 {
        out.errwarn = 2;
    }
    out
}

/// C: `void readqhead(FILE *fpR, char *textraw, char *textbeg, int maxtext,
/// cwp_String *names, cwp_String *ndims, int *janames, int *jandims, int *errwarn);`
/// Scratch storage is internal. Starts at the current position; rewinds on
/// normal completion, retaining C's early-return position for error 2.
/// Non-UTF-8 field payloads use replacement characters in the owned strings;
/// line sizes, prefix checks, and scratch indexing operate on original bytes.
pub fn readqhead(fp: &CFile, maxtext: usize) -> QHead {
    assert!(maxtext > 1);
    let mut out = QHead::default();
    let mut num_names = 0;
    let mut read_names = 0;
    let mut read_ndims = 0;
    // C does not terminate textbeg before its prefix comparisons. Preserve
    // the previous line's initialized bytes when a shorter line overwrites it.
    let mut textbeg = vec![0u8; maxtext.max(11)];
    while let Some(raw) = fgets_bytes(maxtext, fp) {
        if read_names == -1 && read_ndims == -1 {
            break;
        }
        let mut tsize = 0;
        for c in raw {
            if c == 0 || c == b'\n' || c == b'\r' {
                break;
            }
            if c != b' ' && c != b'\t' {
                textbeg[tsize] = c;
                tsize += 1;
            }
        }
        for c in textbeg.iter_mut().take(10) {
            *c = c.to_ascii_lowercase();
        }
        if textbeg.starts_with(b"c_su_names") {
            num_names += 1;
        }
        if read_names > 0 {
            if !textbeg.starts_with(b"c_su_id") {
                out.errwarn = 2;
                return out;
            }
            textbeg[tsize] = 0;
            out.names = tparseq(&String::from_utf8_lossy(&textbeg[..tsize]), b',');
            read_names = -1;
        }
        if textbeg.starts_with(b"c_su_names") {
            read_names = 1;
        }
        if textbeg.starts_with(b"c_su_ndims") {
            textbeg[tsize] = 0;
            out.ndims = tparseq(&String::from_utf8_lossy(&textbeg[..tsize]), b',');
            read_ndims = -1;
        }
        if textbeg.starts_with(b"q,") {
            break;
        }
    }
    rewind(fp);
    if num_names > 1 {
        out.errwarn = 1;
    }
    out
}

/// C: `void qelementout(int iztuple, int ktuple, int numdind, int *kindx, int *ksize);`
/// Output slices need iztuple + max(ktuple-1, 0) entries for iztuple >= 0.
pub fn qelementout(iztuple: i32, ktuple: i32, numdind: i32, kindx: &mut [i32], ksize: &mut [i32]) {
    for k in 0..iztuple {
        kindx[k as usize] = k;
        ksize[k as usize] = 1;
    }
    for k in 0..ktuple - 1 {
        kindx[(k + iztuple) as usize] = iztuple.wrapping_add(k.wrapping_mul(numdind));
        ksize[(k + iztuple) as usize] = numdind;
    }
}

/// C: `void qsplit(struct QInfo *RecInfo, int ncdp, int **ngi, int *ngi_tot,
/// int **ngc, int *ngc_tot, int *errwarn);`
/// Requires nonempty records with three kinf entries. Returns (ngi, ngc).
/// Preserves C's acceptance of an incomplete final crossline group.
pub fn qsplit(records: &[QInfo]) -> Result<(Vec<i32>, Vec<i32>), i32> {
    assert!(!records.is_empty());
    let ncdp = records.len();
    let mut mgi_tot = ncdp;
    for jcdp in 0..ncdp {
        if records[jcdp].kinf[2] != records[0].kinf[2] {
            mgi_tot = jcdp;
            break;
        }
    }
    let mut igc_set = records[0].kinf[2];
    let mut iset_tot = mgi_tot;
    for jcdp in mgi_tot..ncdp {
        if records[jcdp - mgi_tot].kinf[1] != records[jcdp].kinf[1] {
            return Err(jcdp as i32 + 1);
        }
        if igc_set == records[jcdp].kinf[2] {
            iset_tot += 1;
        } else {
            if iset_tot != mgi_tot {
                return Err(jcdp as i32 + 1);
            }
            igc_set = records[jcdp].kinf[2];
            iset_tot = 1;
        }
    }
    let mgc_tot = ncdp / mgi_tot;
    let mut mgi = vec![0; mgi_tot];
    let mut mgc = vec![0; mgc_tot];
    for k in 0..mgi_tot {
        mgi[k] = records[k].kinf[1];
    }
    for k in 0..mgc_tot {
        mgc[k] = records[k * mgi_tot].kinf[2];
    }
    Ok((mgi, mgc))
}

/// C: `void getviacommand(cwp_String **pnameA, int *numpnameA, int *iztupleA,
/// int numdind, int *ktuple, int *ifixd, struct QInfo **RecInfoA, int *ncdpA,
/// double **pindepaA, cwp_String **ndimsA, int *errwarn);`
/// Reads oxwave-core's parameter state. Removes the fixed independent dimension
/// from pname, as C does. C's original numpname is retained separately.
/// Short later non-tuple parameters retain the preceding parameter's tail,
/// as in C's reused dvals buffer. Uninitialized reads or buffer overruns in C
/// instead terminate with an err! diagnostic.
pub fn getviacommand(
    pname: &mut Vec<String>,
    iztuple: i32,
    mut numdind: i32,
) -> Result<QFile, i32> {
    let mut out = QFile {
        numpname: pname.len() as i32,
        iztuple,
        ..QFile::default()
    };
    let ncdp = pname.first().map_or(0, |s| countparval(s));
    if iztuple < 1 || pname.len() < iztuple as usize || ncdp < 1 {
        return Err(1);
    }
    let iz = iztuple as usize;
    for i in 1..iz {
        if ncdp != countnparval(1, &pname[i]) {
            return Err(2);
        }
    }
    let mut iqnto = 0;
    let jtuple;
    if iz == pname.len() {
        out.ifixd = 2;
        out.ktuple = 1;
        jtuple = 1;
        numdind = if numdind > 0 {
            iztuple - 1
        } else {
            numdind.wrapping_neg()
        };
        // C leaves andim[0] uninitialized here; Rust returns no dimension name.
    } else {
        out.ktuple = pname.len() as i32 - iztuple;
        numdind = if numdind > 0 {
            iztuple + numdind * (out.ktuple - 1)
        } else {
            numdind.wrapping_neg()
        };
        if countparname(&pname[iz]) == 1 {
            out.ifixd = 1;
            jtuple = out.ktuple - 1;
            iqnto = countnparval(1, &pname[iz]);
            if iqnto < 1 {
                return Err(3);
            }
            getnpardouble(1, &pname[iz], &mut out.pindepa);
            out.ndims.push(pname.remove(iz));
        } else {
            jtuple = out.ktuple;
            out.ndims.push("vary".to_owned());
        }
    }
    for icdp in 0..ncdp {
        if out.ifixd == 0 {
            iqnto = countnparval(icdp + 1, &pname[iz]);
            if iqnto < 1 {
                return Err(4);
            }
        }
        let used = iztuple + iqnto * jtuple;
        let mut rec = QInfo {
            nto: iqnto,
            dlots: vec![0.0; numdind.max(used) as usize],
            ..QInfo::default()
        };
        if out.ifixd != 2 {
            if out.ifixd == 0 {
                let mut vals = Vec::new();
                getnpardouble(icdp + 1, &pname[iz], &mut vals);
                let start = (iztuple + (jtuple - 1) * iqnto) as usize;
                rec.dlots[start..start + vals.len()].copy_from_slice(&vals);
            }
            for k in 0..jtuple {
                let i = jtuple - 1 - k;
                if iqnto != countnparval(icdp + 1, &pname[iz + i as usize]) {
                    return Err(5);
                }
                let mut vals = Vec::new();
                getnpardouble(icdp + 1, &pname[iz + i as usize], &mut vals);
                let start = (iztuple + k * iqnto) as usize;
                rec.dlots[start..start + vals.len()].copy_from_slice(&vals);
            }
        }
        out.records.push(rec);
    }
    let mut dvals = vec![0.0; ncdp as usize];
    let mut initialized = 0;
    for i in 0..iz {
        let mut vals = Vec::new();
        getpardouble(&pname[i], &mut vals);
        if vals.len() > dvals.len() {
            oxwave_core::err!(
                "qdefine.c: getviacommand: too many non-tuple values for %s",
                pname[i]
            );
        }
        dvals[..vals.len()].copy_from_slice(&vals);
        initialized = initialized.max(vals.len());
        if initialized < dvals.len() {
            oxwave_core::err!(
                "qdefine.c: getviacommand: missing non-tuple values for %s",
                pname[i]
            );
        }
        for icdp in 0..ncdp as usize {
            out.records[icdp].dlots[i] = dvals[icdp];
        }
    }
    out.names = pname.clone();
    Ok(out)
}

/// C: `void getviaqfile(FILE *fpP, cwp_String **pnameA, int *numpnameA,
/// int *iztupleA, int numdind, int *ktuple, int *ifixd, struct QInfo **RecInfoA,
/// int *ncdpA, double **pindepaA, cwp_String **ndimsA, int *errwarn);`
/// Positive numpname selects names; nonpositive excludes abs(numpname) names.
/// Input names must have at least abs(numpname) entries. Output names/ndims
/// contain initialized logical entries only, not C's uninitialized slots.
pub fn getviaqfile(
    fp: &CFile,
    pname: &[String],
    numpname: i32,
    mut numdind: i32,
) -> Result<QFile, i32> {
    assert!(pname.len() >= numpname.unsigned_abs() as usize);
    let stats = getqinfo(fp, 1);
    if stats.errwarn > 0 {
        return Err(stats.errwarn);
    }
    if stats.num_names < 2 {
        return Err(27);
    }
    let head = readqhead(fp, stats.ibsize as usize);
    if head.errwarn > 0 {
        return Err(head.errwarn + 10);
    }
    let janames = head.names.len();
    let jandims = head.ndims.len();
    let ivsize = stats.ivsize.max(stats.num_names).max(stats.num_ndims) as usize;
    let mut aname = Vec::new();
    let mut andim = head.ndims.clone();
    let mut nspot = Vec::new();
    for (i, original) in head.names.iter().enumerate() {
        let mut name = original.clone();
        if i == 0 && name.starts_with("c_su_id") {
            name = "null".to_owned();
        }
        let underscores: Vec<_> = name.match_indices('_').map(|(i, _)| i).collect();
        if underscores.len() >= 2 {
            let k = underscores[0];
            let n = *underscores.last().unwrap();
            if n - k > 1
                && name[..k].bytes().all(|c| c.is_ascii_digit())
                && name[n + 1..].bytes().all(|c| c.is_ascii_digit())
            {
                name = name[k + 1..n].to_owned();
            }
        }
        if name.starts_with("null") {
            continue;
        }
        let found = pname[..numpname.unsigned_abs() as usize].contains(&name);
        if (numpname > 0 && found) || (numpname <= 0 && !found) {
            if jandims == janames {
                andim[aname.len()] = head.ndims[i].clone();
            }
            aname.push(name);
            nspot.push(i as i32);
        }
    }
    let mut numcases = aname.len();
    let mut iztuple = None;
    let mut ktuple = 0;
    let mut kspace = 0;
    'find: for i in 0..numcases {
        for j in i + 1..numcases {
            if aname[j] == aname[i] {
                iztuple = Some(i);
                ktuple = j - i;
                kspace = nspot[j] - nspot[i];
                break 'find;
            }
        }
    }
    let mut out = QFile::default();
    let mut iqnto = 0;
    let jtuple;
    let iz;
    if let Some(first) = iztuple {
        iz = first;
        out.ifixd = if jandims < 2 || andim[0] == "vary" {
            0
        } else {
            1
        };
        if out.ifixd == 1 && janames != jandims {
            return Err(22);
        }
        for i in 0..numcases {
            if aname[i] == aname[iz] {
                iqnto += 1;
            }
        }
        for k in 0..ktuple {
            let mut kqnto = 0;
            let mut i = iz + k;
            while i < numcases {
                if aname[iz + k] == aname[i] {
                    kqnto += 1;
                }
                i += ktuple;
            }
            if kqnto != iqnto {
                return Err(23);
            }
        }
        // C's inner loop can extend past ivsize by ktuple-1; provide storage.
        nspot.resize(ivsize + ktuple, 0);
        let mut i = iz + ktuple;
        while i < ivsize {
            for k in 0..ktuple {
                nspot[i + k] = nspot[i + k - ktuple] + kspace;
                numcases = i + k;
            }
            i += ktuple;
        }
        if out.ifixd == 1 {
            ktuple += 1;
        }
        jtuple = if out.ifixd == 1 { ktuple - 1 } else { ktuple };
        numdind = if numdind > 0 {
            iz as i32 + numdind * (ktuple as i32 - 1)
        } else {
            numdind.wrapping_neg()
        };
        if out.ifixd == 1 {
            for j in 0..iqnto {
                let s = &andim[iz + j * jtuple];
                if s.is_empty() {
                    return Err(24);
                }
                let Some(v) = scan_double(s) else {
                    return Err(25);
                };
                out.pindepa.push(v);
                for k in 1..jtuple {
                    let s = &andim[iz + j * jtuple + k];
                    if !s.starts_with("null") && !s.is_empty() {
                        let Some(test) = scan_double(s) else {
                            return Err(25);
                        };
                        if v != test {
                            return Err(26);
                        }
                    }
                }
            }
        }
    } else {
        iz = numcases;
        out.ifixd = 2;
        jtuple = 1;
        ktuple = 1;
        numdind = if numdind > 0 {
            iz as i32 - 1
        } else {
            numdind.wrapping_neg()
        };
    }
    let mut dfield = vec![0.0; ivsize + ktuple];
    let mut izlast = 0;
    let mut notblank = 0;
    let mut irec = 99;
    rewind(fp);
    while let Some(raw) = fgets(stats.ibsize as usize, fp) {
        irec += 1;
        if raw.starts_with("Q,") || raw.starts_with("q,") {
            // All fields filtered out is undefined in C (uninitialized nspot).
            if numcases == 0 {
                return Err(irec);
            }
            if getqvalscsv(
                &raw,
                b',',
                &nspot[..numcases],
                &mut dfield,
                &mut izlast,
                &mut notblank,
            ) > 0
            {
                return Err(irec);
            }
            if out.ifixd > 0 {
                if izlast != (iqnto * jtuple + iz) as i32 {
                    return Err(irec);
                }
            } else {
                let n = (izlast - iz as i32) / jtuple as i32;
                if n * jtuple as i32 + (iz as i32) < notblank {
                    return Err(irec);
                }
                // Truncation toward zero matters for all-blank tuple tails.
                let n = (notblank - iz as i32) / jtuple as i32;
                if n < 0 {
                    return Err(irec);
                } // invalid C allocation/indexing
                iqnto = n as usize;
            }
            let used = iz + iqnto * jtuple;
            let mut rec = QInfo {
                nto: iqnto as i32,
                dlots: vec![0.0; numdind.max(used as i32) as usize],
                ..QInfo::default()
            };
            rec.dlots[..iz].copy_from_slice(&dfield[..iz]);
            for j in 0..iqnto {
                for k in 0..jtuple {
                    let i = jtuple - 1 - k;
                    rec.dlots[iz + i * iqnto + j] = dfield[iz + j * jtuple + k];
                }
            }
            out.records.push(rec);
        }
    }
    out.iztuple = iz as i32;
    out.ktuple = ktuple as i32;
    out.numpname = (iz + ktuple) as i32;
    aname.truncate((iz + ktuple).min(aname.len()));
    out.names = aname;
    out.ndims = andim;
    Ok(out)
}
