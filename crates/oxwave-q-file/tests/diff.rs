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

//! Differential tests against su/lib/qdefine.c. The source credits the q-file
//! standard to December 2021, without naming an author.

use oxwave_cref::{CArgv, XorShift};
use oxwave_q_file as rust;
use std::ffi::{c_char, c_int, c_void, CString};
use std::ptr::null_mut;
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

#[repr(C)]
struct CQInfo {
    kinf: *mut c_int,
    nto: c_int,
    dlots: *mut f64,
}

extern "C" {
    fn getqinfo(
        fp: *mut c_void,
        num_names: *mut c_int,
        num_ndims: *mut c_int,
        iqflag: c_int,
        ibsize: *mut c_int,
        ivsize: *mut c_int,
        irtotal: *mut c_int,
        errwarn: *mut c_int,
    );
    fn getqvalscsv(
        textraw: *mut c_char,
        textbeg: *mut c_char,
        maxtext: c_int,
        rdel: c_char,
        nspot: *mut c_int,
        numcases: c_int,
        dfield: *mut f64,
        izlast: *mut c_int,
        notblank: *mut c_int,
        errwarn: *mut c_int,
    );
    fn getviacommand(
        pname: *mut *mut *mut c_char,
        numpname: *mut c_int,
        iztuple: *mut c_int,
        numdind: c_int,
        ktuple: *mut c_int,
        ifixd: *mut c_int,
        records: *mut *mut CQInfo,
        ncdp: *mut c_int,
        pindepa: *mut *mut f64,
        ndims: *mut *mut *mut c_char,
        errwarn: *mut c_int,
    );
    fn getviaqfile(
        fp: *mut c_void,
        pname: *mut *mut *mut c_char,
        numpname: *mut c_int,
        iztuple: *mut c_int,
        numdind: c_int,
        ktuple: *mut c_int,
        ifixd: *mut c_int,
        records: *mut *mut CQInfo,
        ncdp: *mut c_int,
        pindepa: *mut *mut f64,
        ndims: *mut *mut *mut c_char,
        errwarn: *mut c_int,
    );
    fn qelementout(
        iztuple: c_int,
        ktuple: c_int,
        numdind: c_int,
        kindx: *mut c_int,
        ksize: *mut c_int,
    );
    fn qsplit(
        records: *mut CQInfo,
        ncdp: c_int,
        ngi: *mut *mut c_int,
        ngi_tot: *mut c_int,
        ngc: *mut *mut c_int,
        ngc_tot: *mut c_int,
        errwarn: *mut c_int,
    );
    fn readqhead(
        fp: *mut c_void,
        textraw: *mut c_char,
        textbeg: *mut c_char,
        maxtext: c_int,
        names: *mut *mut c_char,
        ndims: *mut *mut c_char,
        janames: *mut c_int,
        jandims: *mut c_int,
        errwarn: *mut c_int,
    );
    fn tparseq(tbuf: *mut c_char, d: c_char, fields: *mut *mut c_char, numfields: *mut c_int);
    fn tmpfile() -> *mut c_void;
    fn fwrite(data: *const c_void, size: usize, count: usize, fp: *mut c_void) -> usize;
    fn fseek(fp: *mut c_void, offset: i64, whence: c_int) -> c_int;
    fn ftell(fp: *mut c_void) -> i64;
    fn fclose(fp: *mut c_void) -> c_int;
    fn free(ptr: *mut c_void);
}

// Anonymous files need no filesystem paths and are usable within the sandbox.
struct Files {
    c: *mut c_void,
    r: oxwave_core::cio::CFile,
}
impl Files {
    fn new(text: &str) -> Self {
        Self::from_bytes(text.as_bytes())
    }
    fn from_bytes(text: &[u8]) -> Self {
        let r = oxwave_core::cio::etmpfile();
        oxwave_core::cio::efwrite(text, &r);
        oxwave_core::cio::erewind(&r);
        let c = unsafe { tmpfile() };
        assert!(!c.is_null());
        assert_eq!(
            unsafe { fwrite(text.as_ptr().cast(), 1, text.len(), c) },
            text.len()
        );
        assert_eq!(unsafe { fseek(c, 0, 0) }, 0);
        Self { c, r }
    }
    fn positions_match(&self) {
        assert_eq!(unsafe { ftell(self.c) }, oxwave_core::cio::eftello(&self.r));
    }
}
impl Drop for Files {
    fn drop(&mut self) {
        unsafe {
            fclose(self.c);
        }
        oxwave_core::cio::efclose(&self.r);
    }
}

fn same(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert!(
            a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
            "{a:?} != {b:?}"
        );
    }
}

// tparseq allocates exactly the payload length, without NUL. Never use CStr on
// its outputs: compare only the known payload bytes, then free each allocation.
unsafe fn fields_match(ptrs: &[*mut c_char], expected: &[String]) {
    for (p, s) in ptrs.iter().zip(expected) {
        assert!(!p.is_null());
        assert_eq!(
            unsafe { std::slice::from_raw_parts(p.cast::<u8>(), s.len()) },
            s.as_bytes()
        );
    }
}
unsafe fn free_fields(ptrs: &[*mut c_char]) {
    for &p in ptrs {
        unsafe {
            free(p.cast());
        }
    }
}

fn size(rng: &mut XorShift, case: usize) -> usize {
    let choices = [1, 2, 3, 7, 16, 31, 64, 127, 257, 1024];
    if case < choices.len() {
        choices[case]
    } else {
        choices[rng.next_u32() as usize % 8]
    }
}
fn number(rng: &mut XorShift) -> String {
    match rng.next_u32() % 14 {
        0 => "0".into(),
        1 => "-0".into(),
        2 => "-1".into(),
        3 => "1e-300".into(),
        4 => "-1e300".into(),
        _ => format!("{:.15e}", rng.uniform(-1e8, 1e8)),
    }
}

#[test]
fn tparseq_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(41);
    for case in 0..3000 {
        let n = size(&mut rng, case);
        let delim = [b',', b';', b'|', b'\t'][case % 4];
        let mut text = String::new();
        for i in 0..n {
            if i > 0 {
                text.push(delim as char);
            }
            if rng.next_u32() % 4 != 0 {
                text.push_str(&number(&mut rng));
            }
        }
        if case % 17 == 0 {
            text.push('\0');
            text.push_str("ignored,tail");
        }
        let expected = rust::tparseq(&text, delim);
        let mut input = text.into_bytes();
        input.push(0);
        let mut fields = vec![null_mut(); n + 5];
        let mut count = -1;
        unsafe {
            tparseq(
                input.as_mut_ptr().cast(),
                delim as c_char,
                fields.as_mut_ptr(),
                &mut count,
            );
        }
        assert_eq!(count as usize, expected.len());
        unsafe {
            fields_match(&fields[..count as usize], &expected);
            free_fields(&fields[..count as usize]);
        }
    }
}

#[test]
fn qelementout_matches_c() {
    let mut rng = XorShift::new(42);
    for case in 0..5000 {
        let iz = size(&mut rng, case) as i32 - 1;
        let kt = rng.range_i32(1, 33);
        let nd = if case % 2 == 0 {
            rng.range_i32(-1024, 1024)
        } else {
            rng.next_u32() as i32
        };
        let len = (iz + kt + 5) as usize;
        let mut ci = vec![-999; len];
        let mut cs = ci.clone();
        let mut ri = ci.clone();
        let mut rs = ci.clone();
        unsafe {
            qelementout(iz, kt, nd, ci.as_mut_ptr(), cs.as_mut_ptr());
        }
        rust::qelementout(iz, kt, nd, &mut ri, &mut rs);
        assert_eq!(ci, ri);
        assert_eq!(cs, rs);
    }
}

#[test]
fn getqvalscsv_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(43);
    for case in 0..5000 {
        let n = size(&mut rng, case);
        let delim = [b',', b';', b'|'][case % 3];
        let mut text = String::new();
        for i in 0..n {
            if i > 0 {
                text.push(delim as char);
            }
            let token = match rng.next_u32() % 25 {
                0 => "".into(),
                1 => "   ".into(),
                2 => " 2 3".into(),
                3 => "abc".into(),
                4 => "2e+".into(),
                5 => "2junk".into(),
                6 => "\t-4\t".into(),
                7 => "nan".into(),
                8 => "-inf".into(),
                9 => "0x1.ap3".into(),
                _ => format!(" {} ", number(&mut rng)),
            };
            text.push_str(&token);
        }
        text.push_str(["\n", "\r\n", "", "\0ignored"][case % 4]);
        let mut spots: Vec<i32> = (0..n as i32).filter(|_| rng.next_u32() % 3 != 0).collect();
        if spots.is_empty() {
            spots.push(0);
        }
        if case % 5 == 0 {
            spots.push(n as i32 + 3);
        }
        let mut c = vec![1234.0; spots.len() + 2];
        let mut r = c.clone();
        let (mut cl, mut cn, mut ce) = (-27, -28, -29);
        let (mut rl, mut rn) = (cl, cn);
        let mut input = text.as_bytes().to_vec();
        input.push(0);
        let mut scratch = vec![0i8; input.len() + 20];
        unsafe {
            getqvalscsv(
                input.as_mut_ptr().cast(),
                scratch.as_mut_ptr(),
                input.len() as i32,
                delim as c_char,
                spots.as_mut_ptr(),
                spots.len() as i32,
                c.as_mut_ptr(),
                &mut cl,
                &mut cn,
                &mut ce,
            );
        }
        let re = rust::getqvalscsv(&text, delim, &spots, &mut r, &mut rl, &mut rn);
        assert_eq!((ce, cl, cn), (re, rl, rn), "case {case}: {text:?}");
        same(&c, &r);
    }
}

#[test]
fn getqvalscsv_numeric_prefixes_match_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(49);
    let edges = [
        "",
        " ",
        "\t",
        "+",
        "-",
        ".",
        "+.",
        "-.",
        "0",
        "-0",
        "+0.0",
        "1e",
        "1e+",
        "1e-",
        "-2E-junk",
        "0x",
        "0x.",
        "0xp1",
        "0x1p",
        "0x1p+",
        "0x1p-",
        "0x1.ap3",
        "-0X.8P-4",
        "inf",
        "infinity",
        "infi",
        "infin",
        "infinite",
        "-INFjunk",
        "nan",
        "nan(payload)",
        "nan(",
        "nan(!)",
        "-nan",
        "1e309",
        "-1e309",
        "1e-310",
        "-1e-999",
        "4.9406564584124654e-324",
        "2.2250738585072014e-308",
        "1.7976931348623157e308",
        "0x1p-1074",
        "0x1p1023",
        "1 2",
        " 1 ",
        "\t-4\t",
        "1\t2",
        "1junk",
        "1\0ignored",
        "\u{000b}-2",
        "\u{000c}+3",
    ];
    for case in 0..4000 {
        let token = if case < edges.len() {
            edges[case].to_owned()
        } else if case % 2 == 0 {
            // Exact binary fractions over a broad exponent range, including
            // overflow, subnormals, signed zeros, and ignored suffixes.
            format!(
                "{}0x{:x}.{:03x}p{:+}{}",
                if rng.next_u32() % 2 == 0 { "-" } else { "+" },
                rng.next_u32() % 16,
                rng.next_u32() % 4096,
                rng.range_i32(-1120, 1050),
                if case % 3 == 0 { "junk" } else { "" },
            )
        } else {
            format!(
                "{}.{:09}e{:+}",
                rng.next_u32() as i32,
                rng.next_u32() % 1_000_000_000,
                rng.range_i32(-340, 310),
            )
        };
        let mut input = token.as_bytes().to_vec();
        input.push(0);
        let mut scratch = vec![0i8; input.len() + 1];
        let mut spots = [0];
        let (mut c, mut r) = ([1234.0], [1234.0]);
        let (mut cl, mut cn, mut ce) = (-27, -28, -29);
        let (mut rl, mut rn) = (cl, cn);
        unsafe {
            getqvalscsv(
                input.as_mut_ptr().cast(),
                scratch.as_mut_ptr(),
                input.len() as i32,
                b',' as c_char,
                spots.as_mut_ptr(),
                1,
                c.as_mut_ptr(),
                &mut cl,
                &mut cn,
                &mut ce,
            );
        }
        let re = rust::getqvalscsv(&token, b',', &spots, &mut r, &mut rl, &mut rn);
        assert_eq!((ce, cl, cn), (re, rl, rn), "case {case}: {token:?}");
        assert!(
            c[0].to_bits() == r[0].to_bits() || (c[0].is_nan() && r[0].is_nan()),
            "case {case}: {token:?}: C {:?}, Rust {:?}",
            c[0],
            r[0],
        );
    }
}

fn check_stats(text: &str, flag: i32) {
    check_stats_bytes(text.as_bytes(), flag);
}

fn check_stats_bytes(text: &[u8], flag: i32) {
    let f = Files::from_bytes(text);
    let mut c = rust::QStats {
        num_names: -1,
        num_ndims: -1,
        ibsize: -1,
        ivsize: -1,
        irtotal: -1,
        errwarn: -1,
    };
    unsafe {
        getqinfo(
            f.c,
            &mut c.num_names,
            &mut c.num_ndims,
            flag,
            &mut c.ibsize,
            &mut c.ivsize,
            &mut c.irtotal,
            &mut c.errwarn,
        );
    }
    assert_eq!(c, rust::getqinfo(&f.r, flag), "{text:?}");
    f.positions_match();
}

#[test]
fn getqinfo_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(44);
    for case in 0..2500 {
        let n = size(&mut rng, case);
        let mut text = "C_SU_NAMES\nC_SU_ID".to_owned();
        for _ in 0..n {
            text.push_str(",field");
        }
        text.push('\n');
        if case % 3 == 0 {
            text.push_str("c_su_ndims,tims,0,1\n");
        }
        match case % 11 {
            0 => text.push_str("C_SU_NAMES\nC_SU_ID,cdp\n"),
            1 => text.push_str("C_SU_NDIMS,vary\nC_SU_NDIMS,vary\n"),
            2 => text.push_str("C_SU_NAMES\nwrong_id,cdp\n"),
            _ => (),
        }
        let count = rng.range_i32(1, 6);
        for _ in 0..count {
            text.push_str("Q");
            for _ in 0..n {
                text.push(',');
                text.push_str(&number(&mut rng));
            }
            text.push_str(if case % 7 == 0 { "\r\n" } else { "\n" });
        }
        check_stats(&text, [0, 1, -1, 2][case % 4]);
    }
    for text in [
        "",
        "# comment only\n",
        "q,1\nC_SU_NAMES\nwrong\n",
        "C_SU_NAMES\n",
    ] {
        check_stats(text, 1);
    }
    // Full 10,000-byte fgets chunks: preserve C's +1 size accounting.
    check_stats(&format!("Q,{}\n", "0,".repeat(16000)), 1);
}

#[test]
fn readqhead_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(45);
    for case in 0..2500 {
        let n = size(&mut rng, case);
        let mut text = "# preliminary comment\nC_SU_NAMES\nC_SU_ID".to_owned();
        for i in 0..n {
            text.push_str(&format!(", \tfield{}_{} ", i, rng.next_u32()));
        }
        text.push('\n');
        if case % 9 == 0 {
            text.push_str("C_SU_NAMES\nC_SU_ID,replaced\n");
        }
        if case % 9 == 1 {
            text.push_str("C_SU_NAMES\nbad_id,field\n");
        }
        if case % 3 != 0 {
            text.push_str(" C_SU_NDIMS , tims ,0,, 2 \n");
        }
        text.push_str("Q,1,2\nC_SU_NAMES\nC_SU_ID,late\n");
        let f = Files::new(&text);
        let maxtext = text.len() + 2;
        let mut raw = vec![0i8; maxtext];
        let mut beg = raw.clone();
        let mut names = vec![null_mut(); n + 20];
        let mut dims = names.clone();
        let (mut nn, mut nd, mut err) = (-1, -1, -1);
        unsafe {
            readqhead(
                f.c,
                raw.as_mut_ptr(),
                beg.as_mut_ptr(),
                maxtext as i32,
                names.as_mut_ptr(),
                dims.as_mut_ptr(),
                &mut nn,
                &mut nd,
                &mut err,
            );
        }
        let r = rust::readqhead(&f.r, maxtext);
        assert_eq!(err, r.errwarn);
        // C does not assign counts on early error 2. Previously written fields
        // are still present and can be compared with their known Rust lengths.
        if err == 0 || err == 1 {
            assert_eq!((nn as usize, nd as usize), (r.names.len(), r.ndims.len()));
        }
        unsafe {
            fields_match(&names[..r.names.len()], &r.names);
            fields_match(&dims[..r.ndims.len()], &r.ndims);
            free_fields(&names[..r.names.len()]);
            free_fields(&dims[..r.ndims.len()]);
        }
        f.positions_match();
    }
}

#[test]
fn readqhead_non_utf8_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(50);
    for case in 0..2500 {
        let n = size(&mut rng, case);
        let payloads: [&[u8]; 6] = [
            b"\xff",
            b"\x80",
            b"\xc0\xaf",
            b"\xe2\x82",
            b"\xff\xfe",
            "é".as_bytes(),
        ];
        let mut field = b"field".to_vec();
        for _ in 0..n {
            field.extend_from_slice(payloads[rng.next_u32() as usize % payloads.len()]);
        }
        let mut text = b"C_SU_NAMES\nC_SU_ID,".to_vec();
        text.extend_from_slice(&field);
        text.extend_from_slice(b"\nC_SU_NDIMS,");
        text.extend_from_slice(&field);
        text.extend_from_slice(b"\nQ,1\n");
        // A line fills (or nearly fills) the fgets buffer. Lossy decoding
        // expands beyond this capacity, including with a single invalid byte.
        let maxtext = text
            .split(|&b| b == b'\n')
            .map(|s| s.len() + 2)
            .max()
            .unwrap()
            + case % 2;
        check_stats_bytes(&text, 1);
        let f = Files::from_bytes(&text);
        let mut raw = vec![0i8; maxtext];
        let mut beg = raw.clone();
        let mut names = [null_mut(); 2];
        let mut dims = [null_mut(); 2];
        let (mut nn, mut nd, mut err) = (-1, -1, -1);
        unsafe {
            readqhead(
                f.c,
                raw.as_mut_ptr(),
                beg.as_mut_ptr(),
                maxtext as i32,
                names.as_mut_ptr(),
                dims.as_mut_ptr(),
                &mut nn,
                &mut nd,
                &mut err,
            );
        }
        let r = rust::readqhead(&f.r, maxtext);
        assert_eq!((err, nn, nd), (0, 2, 2));
        assert_eq!(r.errwarn, 0);
        assert_eq!(
            r.names,
            [
                "c_su_id".to_owned(),
                String::from_utf8_lossy(&field).into_owned()
            ]
        );
        assert_eq!(
            r.ndims,
            [
                "c_su_ndims".to_owned(),
                String::from_utf8_lossy(&field).into_owned()
            ]
        );
        // C returns raw, unterminated payloads. Compare their known byte
        // lengths before checking the Rust API's lossy string representation.
        unsafe {
            for (ptrs, first) in [(&names, &b"c_su_id"[..]), (&dims, &b"c_su_ndims"[..])] {
                assert_eq!(
                    std::slice::from_raw_parts(ptrs[0].cast::<u8>(), first.len()),
                    first
                );
                assert_eq!(
                    std::slice::from_raw_parts(ptrs[1].cast::<u8>(), field.len()),
                    field
                );
                free_fields(ptrs);
            }
        }
        f.positions_match();
    }
    // Invalid bytes on the required ID line still yield C's header error 2.
    let mut text = b"C_SU_NAMES\n".to_vec();
    text.extend_from_slice(&[0xff; 32]);
    text.push(b'\n');
    let f = Files::from_bytes(&text);
    let mut raw = [0i8; 34];
    let mut beg = raw;
    let mut names = [null_mut(); 1];
    let mut dims = names;
    let (mut nn, mut nd, mut err) = (-1, -1, -1);
    unsafe {
        readqhead(
            f.c,
            raw.as_mut_ptr(),
            beg.as_mut_ptr(),
            34,
            names.as_mut_ptr(),
            dims.as_mut_ptr(),
            &mut nn,
            &mut nd,
            &mut err,
        );
    }
    assert_eq!(err, 2);
    assert_eq!(rust::readqhead(&f.r, 34).errwarn, err);
    f.positions_match();
}

#[test]
fn qsplit_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(46);
    for case in 0..4000 {
        let ni = size(&mut rng, case);
        let nc = rng.range_i32(1, 8) as usize;
        let start_i = rng.range_i32(-10000, 10000);
        let start_c = rng.range_i32(-10000, 10000);
        let mut records = Vec::new();
        for c in 0..nc {
            for i in 0..ni {
                records.push(rust::QInfo {
                    kinf: vec![0, start_i + i as i32 * 3, start_c + c as i32 * 7],
                    ..Default::default()
                });
            }
        }
        match case % 4 {
            0 if records.len() > 1 => {
                let i = rng.next_u32() as usize % records.len();
                records[i].kinf[1] += 1;
            }
            1 if records.len() > 1 => {
                let i = rng.next_u32() as usize % records.len();
                records[i].kinf[2] += 1;
            }
            2 if records.len() > 1 => {
                records.pop();
            }
            _ => (),
        }
        let mut cr: Vec<_> = records
            .iter_mut()
            .map(|r| CQInfo {
                kinf: r.kinf.as_mut_ptr(),
                nto: 0,
                dlots: null_mut(),
            })
            .collect();
        let (mut gi, mut gc) = (null_mut(), null_mut());
        let (mut ni, mut nc, mut err) = (-1, -1, -1);
        unsafe {
            qsplit(
                cr.as_mut_ptr(),
                cr.len() as i32,
                &mut gi,
                &mut ni,
                &mut gc,
                &mut nc,
                &mut err,
            );
        }
        let r = rust::qsplit(&records);
        if err > 0 {
            assert_eq!(r, Err(err));
        } else {
            let (ri, rc) = r.unwrap();
            assert_eq!(unsafe { std::slice::from_raw_parts(gi, ni as usize) }, ri);
            assert_eq!(unsafe { std::slice::from_raw_parts(gc, nc as usize) }, rc);
            unsafe {
                free(gi.cast());
                free(gc.cast());
            }
        }
    }
}

struct CRead {
    names: *mut *mut c_char,
    count: i32,
    iz: i32,
    kt: i32,
    fixed: i32,
    records: *mut CQInfo,
    ncdp: i32,
    indep: *mut f64,
    dims: *mut *mut c_char,
    err: i32,
}
impl CRead {
    fn new(names: *mut *mut c_char, count: i32, iz: i32) -> Self {
        Self {
            names,
            count,
            iz,
            kt: -1,
            fixed: -1,
            records: null_mut(),
            ncdp: -1,
            indep: null_mut(),
            dims: null_mut(),
            err: -1,
        }
    }
    fn compare(&self, result: &Result<rust::QFile, i32>) {
        if self.err > 0 {
            assert_eq!(result.as_ref().unwrap_err(), &self.err);
            return;
        }
        let r = result.as_ref().unwrap();
        assert_eq!(
            (self.count, self.iz, self.kt, self.fixed, self.ncdp),
            (
                r.numpname,
                r.iztuple,
                r.ktuple,
                r.ifixd,
                r.records.len() as i32
            )
        );
        unsafe {
            let names = std::slice::from_raw_parts(self.names, r.names.len());
            fields_match(names, &r.names);
            if !r.ndims.is_empty() {
                fields_match(
                    std::slice::from_raw_parts(self.dims, r.ndims.len()),
                    &r.ndims,
                );
            }
            for (i, rec) in r.records.iter().enumerate() {
                let c = &*self.records.add(i);
                assert_eq!(c.nto, rec.nto);
                let jt = if r.ifixd == 0 { r.ktuple } else { r.ktuple - 1 };
                let used = (r.iztuple + rec.nto * jt) as usize;
                same(
                    &rec.dlots[..used],
                    std::slice::from_raw_parts(c.dlots, used),
                );
            }
            if !r.pindepa.is_empty() {
                same(
                    &r.pindepa,
                    std::slice::from_raw_parts(self.indep, r.pindepa.len()),
                );
            }
        }
    }
    fn free_numeric(&self) {
        // Header strings in getviaqfile are aliased/overwritten and not all
        // allocations remain reachable. Keep those C allocations alive; free
        // the uniquely owned numeric arrays. This also avoids reading any
        // unterminated C string to discover its length.
        if self.err == 0 {
            unsafe {
                for i in 0..self.ncdp as usize {
                    free((*self.records.add(i)).dlots.cast());
                }
                free(self.records.cast());
                free(self.indep.cast());
                free(self.dims.cast());
            }
        }
    }
}

#[test]
fn getviacommand_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(47);
    for case in 0..2500 {
        let n = size(&mut rng, case).min(if case < 10 { 1024 } else { 16 });
        let mode = case % 3;
        let tuple = 2 + (rng.next_u32() % 3) as usize;
        let mut names = vec!["cdp".to_owned(), "elev".to_owned()];
        if mode != 2 {
            for k in 0..tuple {
                names.push(format!("tuple{k}"));
            }
        }
        let mut args = vec!["diff".to_owned()];
        for name in &names[..2] {
            let vals: Vec<_> = (0..n).map(|_| number(&mut rng)).collect();
            args.push(format!("{name}={}", vals.join(",")));
        }
        let fixed_count = size(&mut rng, case).min(if case < 10 { 1024 } else { 31 });
        if mode != 2 {
            if mode == 1 {
                args.push(format!(
                    "tuple0={}",
                    (0..fixed_count)
                        .map(|_| number(&mut rng))
                        .collect::<Vec<_>>()
                        .join(",")
                ));
            }
            for i in 0..n {
                let count = if mode == 1 {
                    fixed_count
                } else {
                    1 + (rng.next_u32() % 15) as usize
                };
                for k in if mode == 1 { 1 } else { 0 }..tuple {
                    args.push(format!(
                        "tuple{k}={}",
                        (0..count)
                            .map(|_| number(&mut rng))
                            .collect::<Vec<_>>()
                            .join(",")
                    ));
                }
                let _ = i;
            }
        }
        // Cover errors 1, 2, 4, 5 as well as all three successful layouts.
        let iz = if case % 29 == 0 { 0 } else { 2 };
        if case % 29 == 1 {
            args[2] = "elev=1,2,3,4,5".into();
        }
        if case % 29 == 2 && mode == 0 {
            args.retain(|s| !s.starts_with("tuple0="));
        }
        if case % 29 == 3 && mode != 2 {
            let prefix = format!("tuple{}=", tuple - 1);
            let arg = args.iter_mut().find(|s| s.starts_with(&prefix)).unwrap();
            *arg = format!("{prefix}1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16");
        }
        let refs: Vec<_> = args.iter().map(String::as_str).collect();
        let mut argv = CArgv::new(&refs);
        unsafe {
            oxwave_cref::initargs(argv.argc(), argv.argv());
        }
        oxwave_core::par::initargs(&args);
        let cnames: Vec<_> = names
            .iter()
            .map(|s| CString::new(s.as_str()).unwrap())
            .collect();
        let mut pointers: Vec<_> = cnames.iter().map(|s| s.as_ptr() as *mut c_char).collect();
        let mut c = CRead::new(pointers.as_mut_ptr(), names.len() as i32, iz);
        let nd = if case % 2 == 0 {
            rng.range_i32(-100, 0)
        } else {
            rng.range_i32(1, 130)
        };
        unsafe {
            getviacommand(
                &mut c.names,
                &mut c.count,
                &mut c.iz,
                nd,
                &mut c.kt,
                &mut c.fixed,
                &mut c.records,
                &mut c.ncdp,
                &mut c.indep,
                &mut c.dims,
                &mut c.err,
            );
        }
        let result = rust::getviacommand(&mut names, iz, nd);
        c.compare(&result);
        if let Ok(r) = &result {
            assert_eq!(names, r.names);
            if !r.ndims.is_empty() {
                unsafe {
                    fields_match(std::slice::from_raw_parts(c.dims, 1), &r.ndims);
                }
            }
        }
        c.free_numeric();
    }
}

#[test]
fn getviacommand_short_later_parameters_match_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(51);
    for case in 0..2500 {
        let n = size(&mut rng, case).max(2);
        let a: Vec<_> = (0..n).map(|_| number(&mut rng)).collect();
        let short = rng.next_u32() as usize % n;
        let b: Vec<_> = (0..short).map(|_| number(&mut rng)).collect();
        let mut args = vec!["diff".to_owned(), format!("a={}", a.join(","))];
        match case % 3 {
            0 => {
                // Validation sees the first b, while getpardouble reads the last.
                args.push(format!("b={}", a.join(",")));
                args.push(format!("b={}", b.join(",")));
            }
            1 => {
                // A trailing comma is counted but is not converted. Its tail
                // remains initialized by the preceding non-tuple parameter.
                args.push(format!("b={},", a[..n - 1].join(",")));
            }
            _ => {
                args.push(format!("b={}", a.join(",")));
                args.push(format!("b={},", b.join(",")));
            }
        }
        // A second short read must retain b's prefix and a's remaining tail.
        args.push(format!("c={}", a.join(",")));
        args.push("c=-0".to_owned());
        let refs: Vec<_> = args.iter().map(String::as_str).collect();
        let mut argv = CArgv::new(&refs);
        unsafe {
            oxwave_cref::initargs(argv.argc(), argv.argv());
        }
        oxwave_core::par::initargs(&args);
        let cnames: Vec<_> = ["a", "b", "c"].map(|s| CString::new(s).unwrap()).into();
        let mut pointers: Vec<_> = cnames.iter().map(|s| s.as_ptr() as *mut c_char).collect();
        let mut c = CRead::new(pointers.as_mut_ptr(), 3, 3);
        unsafe {
            getviacommand(
                &mut c.names,
                &mut c.count,
                &mut c.iz,
                0,
                &mut c.kt,
                &mut c.fixed,
                &mut c.records,
                &mut c.ncdp,
                &mut c.indep,
                &mut c.dims,
                &mut c.err,
            );
        }
        let mut names = vec!["a".to_owned(), "b".to_owned(), "c".to_owned()];
        let r = rust::getviacommand(&mut names, 3, 0);
        assert_eq!(c.err, 0, "case {case}");
        c.compare(&r);
        c.free_numeric();
    }
}

#[test]
fn getviacommand_rejects_undefined_reads_and_overruns() {
    const CHILD_ENV: &str = "OXWAVE_Q_FILE_BAD_COMMAND";
    if let Ok(case) = std::env::var(CHILD_ENV) {
        let values = match case.as_str() {
            "empty" => vec!["a="],
            "trailing" => vec!["a=1,"],
            "overrun" => vec!["a=1,2", "b=1,2", "b=3,4,5"],
            _ => unreachable!(),
        };
        let args: Vec<_> = std::iter::once("diff")
            .chain(values)
            .map(str::to_owned)
            .collect();
        oxwave_core::par::initargs(&args);
        let mut names = if case == "overrun" {
            vec!["a".to_owned(), "b".to_owned()]
        } else {
            vec!["a".to_owned()]
        };
        let iz = names.len() as i32;
        let _ = rust::getviacommand(&mut names, iz, 0);
        panic!("undefined C input must exit through err!");
    }
    // The C oracle would read uninitialized doubles or overwrite its heap;
    // exercise only Rust, in subprocesses to verify the err! exit and message.
    for (case, message) in [
        ("empty", "missing non-tuple values for a"),
        ("trailing", "missing non-tuple values for a"),
        ("overrun", "too many non-tuple values for b"),
    ] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "getviacommand_rejects_undefined_reads_and_overruns",
                "--nocapture",
            ])
            .env(CHILD_ENV, case)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{case}: {output:?}");
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            format!("\ndiff: qdefine.c: getviacommand: {message}\n")
        );
    }
}

fn check_qfile(text: &str, filters: &[String], count: i32, nd: i32) {
    let f = Files::new(text);
    let cnames: Vec<_> = filters
        .iter()
        .map(|s| CString::new(s.as_str()).unwrap())
        .collect();
    let mut pointers: Vec<_> = cnames.iter().map(|s| s.as_ptr() as *mut c_char).collect();
    let mut c = CRead::new(pointers.as_mut_ptr(), count, -1);
    unsafe {
        getviaqfile(
            f.c,
            &mut c.names,
            &mut c.count,
            &mut c.iz,
            nd,
            &mut c.kt,
            &mut c.fixed,
            &mut c.records,
            &mut c.ncdp,
            &mut c.indep,
            &mut c.dims,
            &mut c.err,
        );
    }
    let result = rust::getviaqfile(&f.r, filters, count, nd);
    c.compare(&result);
    f.positions_match();
    c.free_numeric();
}

#[test]
fn getviaqfile_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(48);
    for case in 0..2500 {
        let mode = case % 3; // varying, fixed, no tuples
        let nt = size(&mut rng, case)
            .max(2)
            .min(if case < 10 { 1024 } else { 31 });
        let jt = if mode == 2 { 0 } else { 2 };
        let padding = if jt > 0 && case % 5 > 1 { 1 } else { 0 };
        let nr = if case < 10 {
            size(&mut rng, case).min(31)
        } else {
            rng.range_i32(1, 5) as usize
        };
        let mut text = "C_SU_NAMES\nC_SU_ID,cdp,elev".to_owned();
        for j in 0..nt {
            for k in 0..jt {
                text.push_str(&format!(",{:02}_val{k}_{j}", j));
            }
        }
        // Padding null fields avoids the C nspot allocation overrun: ivsize-iz
        // is a multiple of the stored tuple width, including the record ID.
        for _ in 0..padding {
            text.push_str(",null");
        }
        text.push('\n');
        // Varying input has no C_SU_NDIMS record. A literal `vary` in its
        // second field still enters C's fixed path and triggers error 22.
        if mode == 0 && case % 11 == 0 {
            text.push_str("C_SU_NDIMS,vary\n");
        }
        if mode == 1 {
            text.push_str("C_SU_NDIMS,tims,null");
            for j in 0..nt {
                text.push_str(&format!(",{j},{j}"));
            }
            for _ in 0..padding {
                text.push_str(",null");
            }
            text.push('\n');
        }
        for i in 0..nr {
            text.push_str(&format!("Q,{},{}", number(&mut rng), number(&mut rng)));
            let used = if mode == 0 {
                1 + (rng.next_u32() as usize % nt)
            } else {
                nt
            };
            for _ in 0..used {
                for _ in 0..jt {
                    text.push(',');
                    // Blank selected values become zero; a partially blank last
                    // tuple is dropped by the C varying-dimension reader.
                    match rng.next_u32() % 12 {
                        0 => (),
                        1 => text.push_str("   "),
                        2 => text.push_str(&format!(" {} ", number(&mut rng))),
                        _ => text.push_str(&number(&mut rng)),
                    }
                }
            }
            for _ in 0..padding {
                text.push(',');
            }
            text.push('\n');
            let _ = i;
        }
        let (filters, count) = match case % 5 {
            0 => (vec!["cdp".into(), "val0".into(), "val1".into()], 3),
            1 => (vec!["elev".into()], -1),
            2 => (vec!["nonexistent".into()], -1),
            _ => (Vec::new(), 0),
        };
        let nd = if case % 2 == 0 {
            rng.range_i32(-100, 0)
        } else {
            rng.range_i32(1, 100)
        };
        check_qfile(&text, &filters, count, nd);
    }
    for text in [
        "",
        "Q,1\n",
        "C_SU_NAMES\nwrong,cdp\n",
        "C_SU_NAMES\nC_SU_ID,cdp\nQ,abc\n",
        "C_SU_NAMES\nC_SU_ID,cdp\nQ,1 2\n",
        "C_SU_NAMES\nC_SU_ID,cdp,elev\nQ,1\n",
    ] {
        check_qfile(text, &[], 0, 10);
    }
    // Header errors are deterministic even though the C returns no records.
    for tail in [
        "C_SU_ID,cdp,t,v,t,other\n", // 23: inconsistent tuple members
        "C_SU_ID,cdp,v,w,v,w\nC_SU_NDIMS,tims,0,0\n", // 22: wrong dimension count
        "C_SU_ID,cdp,v,w,v,w\nC_SU_NDIMS,tims,abc,0,1,1\n", // 25: nonnumeric dimension
        "C_SU_ID,cdp,v,w,v,w\nC_SU_NDIMS,tims,0,3,1,1\n", // 26: unequal dimensions
    ] {
        check_qfile(&format!("C_SU_NAMES\n{tail}"), &[], 0, 8);
    }
}
