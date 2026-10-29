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

//! Differential tests for su/lib/readkfile.c.
//! Credits: Copyright (c) Colorado School of Mines, 2021. All rights reserved.
//! The C source gives no named author credit.

use oxwave_core::cio;
use oxwave_cref::XorShift;
use oxwave_k_file as rust;
use std::ffi::{c_char, c_int, c_void, CString};
use std::ptr::null_mut;
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

extern "C" {
    fn getkvalscsv(
        textraw: *mut c_char,
        textbeg: *mut c_char,
        maxtext: c_int,
        rdel: c_char,
        dfield: *mut f64,
        nspot: *mut c_int,
        numcases: c_int,
        ncount: c_int,
        comerr: *mut c_int,
        morerr: *mut c_int,
        numerr: *mut c_int,
        nblank: *mut c_int,
    );
    fn readkfile(
        fp_r: *mut c_void,
        names: *mut *mut c_char,
        forms: *mut *mut c_char,
        dfield: *mut f64,
        numcasesout: *mut c_int,
        errwarn: *mut c_int,
    );
    fn tparse(tbuf: *mut c_char, d: c_char, fields: *mut *mut c_char, numfields: *mut c_int);
    fn writekfile(
        fp_w: *mut c_void,
        names: *mut *mut c_char,
        forms: *mut *mut c_char,
        dfield: *mut f64,
        numcasesout: c_int,
        errwarn: *mut c_int,
    );
    fn tmpfile() -> *mut c_void;
    fn fwrite(data: *const c_void, size: usize, n: usize, fp: *mut c_void) -> usize;
    fn fread(data: *mut c_void, size: usize, n: usize, fp: *mut c_void) -> usize;
    fn fseek(fp: *mut c_void, offset: i64, whence: c_int) -> c_int;
    fn ftell(fp: *mut c_void) -> i64;
    fn fclose(fp: *mut c_void) -> c_int;
    fn free(ptr: *mut c_void);
}

struct Files {
    c: *mut c_void,
    r: cio::CFile,
}
impl Files {
    fn new(text: &str) -> Self {
        Self::from_bytes(text.as_bytes())
    }
    fn from_bytes(text: &[u8]) -> Self {
        let r = cio::etmpfile();
        cio::efwrite(text, &r);
        cio::erewind(&r);
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
        assert_eq!(unsafe { ftell(self.c) }, cio::eftello(&self.r));
    }
    fn contents(&self) -> (Vec<u8>, Vec<u8>) {
        self.positions_match();
        let n = unsafe { ftell(self.c) } as usize;
        assert_eq!(unsafe { fseek(self.c, 0, 0) }, 0);
        let mut c = vec![0; n];
        assert_eq!(unsafe { fread(c.as_mut_ptr().cast(), 1, n, self.c) }, n);
        cio::erewind(&self.r);
        let mut r = vec![0; n + 1];
        let count = cio::efread(&mut r, &self.r);
        r.truncate(count);
        (c, r)
    }
}
impl Drop for Files {
    fn drop(&mut self) {
        unsafe {
            fclose(self.c);
        }
        cio::efclose(&self.r);
    }
}

fn same(a: &[f64], b: &[f64], context: &str) {
    assert_eq!(a.len(), b.len());
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        assert!(
            a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
            "{context}: field {i}: {a:?} != {b:?}"
        );
    }
}

fn size(rng: &mut XorShift, case: usize) -> usize {
    let choices = [1, 2, 3, 7, 16, 31, 64, 127, 257];
    if case < choices.len() {
        choices[case]
    } else {
        choices[rng.next_u32() as usize % 7]
    }
}

fn number(rng: &mut XorShift) -> String {
    match rng.next_u32() % 24 {
        0 => "0".into(),
        1 => "-0".into(),
        2 => "-1".into(),
        3 => "1e-300".into(),
        4 => "-1e300".into(),
        5 => "1e309".into(),
        6 => "-1e-400".into(),
        7 => "0x1.234p-15".into(),
        8 => "inf".into(),
        9 => "-infinity".into(),
        10 => "nan".into(),
        11 => "nan(payload)".into(),
        12 => "1.25trailing".into(),
        13 => "1e+".into(),
        14 => "bad".into(),
        15 => "".into(),
        16 => "   ".into(),
        17 => "1 2".into(),
        18 => "  -1.5   ".into(),
        19 => "\t-3\t4".into(),
        20 => "\t".into(),
        21 => "1\t2".into(),
        22 => "-".into(),
        _ => format!("{:.17e}", rng.uniform(-1e12, 1e12)),
    }
}

// The original allocates exactly the payload length, WITHOUT a terminator.
// Read only known payload bytes, never CStr::from_ptr.
unsafe fn fields_match(ptrs: &[*mut c_char], expected: &[String]) {
    for (&p, s) in ptrs.iter().zip(expected) {
        assert!(!p.is_null());
        let actual = unsafe { std::slice::from_raw_parts(p.cast::<u8>(), s.len()) };
        assert_eq!(actual, s.as_bytes());
    }
}

#[test]
fn tparse_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b01);
    for case in 0..4000 {
        let n = if case == 9 {
            4096
        } else {
            size(&mut rng, case)
        };
        let delimiter = [b',', b';', b'|', b':', 0][case % 5];
        let mut text = String::new();
        for i in 0..n {
            if i > 0 {
                text.push(delimiter as char);
            }
            if rng.next_u32() % 4 != 0 {
                text.push_str(&format!("-{:x}.{}", rng.next_u32(), rng.next_u32()));
            }
        }
        if case % 17 == 0 {
            text.push_str("\0ignored,tail");
        }
        let expected = rust::tparse(&text, delimiter);
        let mut input = text.into_bytes();
        input.push(0);
        let mut fields = vec![null_mut(); n + 5];
        let mut count = -5;
        unsafe {
            tparse(
                input.as_mut_ptr().cast(),
                delimiter as c_char,
                fields.as_mut_ptr(),
                &mut count,
            );
        }
        assert_eq!(count as usize, expected.len(), "case {case}");
        unsafe {
            fields_match(&fields[..count as usize], &expected);
        }
        for p in &fields[..count as usize] {
            unsafe {
                free((*p).cast());
            }
        }
    }
}

fn compare_csv(raw: &mut [u8], spots: &mut [i32], initial: [i32; 4], ncount: i32) {
    let mut cv = vec![-987.5; spots.len()];
    let mut rv = cv.clone();
    let mut cb = vec![b'!'; raw.len() + 1];
    let mut rb = cb.clone();
    let [mut cc, mut cm, mut ce, mut cn] = initial;
    let [mut rc, mut rm, mut re, mut rn] = initial;
    unsafe {
        getkvalscsv(
            raw.as_mut_ptr().cast(),
            cb.as_mut_ptr().cast(),
            raw.len() as i32,
            b',' as c_char,
            cv.as_mut_ptr(),
            spots.as_mut_ptr(),
            spots.len() as i32,
            ncount,
            &mut cc,
            &mut cm,
            &mut ce,
            &mut cn,
        );
    }
    rust::getkvalscsv(
        raw, &mut rb, b',', &mut rv, spots, ncount, &mut rc, &mut rm, &mut re, &mut rn,
    );
    assert_eq!(
        [cc, cm, ce, cn],
        [rc, rm, re, rn],
        "raw {raw:?}, spots {spots:?}"
    );
    assert_eq!(cb, rb, "scratch buffer, raw {raw:?}");
    same(&cv, &rv, &format!("CSV {raw:?}"));
}

#[test]
fn getkvalscsv_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b02);
    for case in 0..6000 {
        let n = size(&mut rng, case);
        let mut text = String::new();
        let mut spots = Vec::new();
        for i in 0..n {
            if i > 0 {
                text.push(',');
            }
            text.push_str(&number(&mut rng));
            if rng.next_u32() % 3 != 0 {
                spots.push(i as i32);
            }
        }
        if spots.is_empty() {
            spots.push((n - 1) as i32);
        }
        if case % 3 == 0 {
            spots.push(n as i32 + rng.range_i32(1, 9));
        }
        text.push([0, b'\n', b'\r'][case % 3] as char);
        // Full buffer scans, truncation before a delimiter, and nonzero stale
        // bytes beyond the first NUL are all deliberate reference behaviors.
        let mut raw = text.into_bytes();
        let len = raw.len();
        raw.resize(len + rng.range_i32(1, 40) as usize, 0);
        if case % 13 == 0 {
            raw[len..].fill(b'x');
        }
        if case % 19 == 0 {
            raw.truncate((len / 2).max(1));
        }
        let initial = [
            rng.range_i32(-10, 10),
            rng.range_i32(-10, 10),
            rng.range_i32(-10, 10),
            rng.range_i32(-10, 10),
        ];
        let ncount = rng.range_i32(-100, 10000);
        compare_csv(&mut raw, &mut spots, initial, ncount);
    }
    // Near C's full fixed buffer, including one very long field.
    for n in [1, 2, 7, 31, 1024, 10000, 10001] {
        let mut raw = vec![b'0'; n];
        raw[n - 1] = 0;
        compare_csv(&mut raw, &mut [0], [0; 4], -1);
    }
}

#[test]
fn getkvalscsv_other_delimiters_match_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b03);
    for case in 0..2000 {
        let delim = [b';', b'|', b'\t', 0][case % 4];
        let n = size(&mut rng, case);
        let mut raw = Vec::new();
        for i in 0..n {
            if i > 0 {
                raw.push(delim);
            }
            raw.extend_from_slice(number(&mut rng).as_bytes());
        }
        raw.push(0);
        raw.resize(raw.len() + n + 8, 0);
        let mut spots: Vec<i32> = (0..n as i32).collect();
        let mut cv = vec![0.; n];
        let mut rv = cv.clone();
        let mut cb = vec![0; raw.len() + 1];
        let mut rb = cb.clone();
        let (mut cc, mut cm, mut ce, mut cn) = (0, 0, 0, 0);
        let (mut rc, mut rm, mut re, mut rn) = (0, 0, 0, 0);
        unsafe {
            getkvalscsv(
                raw.as_mut_ptr().cast(),
                cb.as_mut_ptr().cast(),
                raw.len() as i32,
                delim as c_char,
                cv.as_mut_ptr(),
                spots.as_mut_ptr(),
                n as i32,
                case as i32,
                &mut cc,
                &mut cm,
                &mut ce,
                &mut cn,
            );
        }
        rust::getkvalscsv(
            &raw,
            &mut rb,
            delim,
            &mut rv,
            &spots,
            case as i32,
            &mut rc,
            &mut rm,
            &mut re,
            &mut rn,
        );
        assert_eq!([cc, cm, ce, cn], [rc, rm, re, rn]);
        assert_eq!(cb, rb);
        same(&cv, &rv, "alternate delimiter");
    }
}

#[test]
fn getkvalscsv_numeric_prefixes_match_c() {
    let _guard = C_LOCK.lock().unwrap();
    for text in [
        "in",
        "infi",
        "infin",
        "infinity!",
        "+INF",
        "-NaN",
        "nan(",
        "nan(a-b)",
        "0x",
        "0x.",
        "0x.p3",
        "0x1p+",
        "0X.8P-1073",
        "1e+",
        "1e-xyz",
        "+.e",
        ".",
        "+",
        "1.2.3",
        "\u{000b}-2.5",
        "\u{000c}7",
        "-0x1p-1075",
        "0x1.00000000000018p-1022",
    ] {
        let mut raw = text.as_bytes().to_vec();
        raw.push(0);
        compare_csv(&mut raw, &mut [0], [0; 4], 0);
    }
    let mut rng = XorShift::new(0x4b06);
    for case in 0..5000 {
        let text = if case % 2 == 0 {
            // All bit patterns, including subnormals and nonfinite values.
            format!("{:.17e}tail", f64::from_bits(rng.next_u64()))
        } else {
            format!(
                "{}0x{:016x}.{:016x}p{:+}tail",
                if rng.next_u32() & 1 == 0 { "-" } else { "" },
                rng.next_u64(),
                rng.next_u64(),
                rng.range_i32(-1200, 1100)
            )
        };
        let mut raw = text.into_bytes();
        raw.push(0);
        compare_csv(&mut raw, &mut [0], [0; 4], case);
    }
}

fn compare_read(text: &str, capacity: usize) -> (i32, i32) {
    let files = Files::new(text);
    let mut cnames = vec![null_mut(); capacity];
    let mut cforms = vec![null_mut(); capacity];
    let mut rnames = vec![String::new(); capacity];
    let mut rforms = rnames.clone();
    let mut cv = vec![-6789.5; capacity];
    let mut rv = cv.clone();
    let (mut nc, mut ec) = (-123, 987);
    let (mut nr, mut er) = (nc, ec);
    unsafe {
        readkfile(
            files.c,
            cnames.as_mut_ptr(),
            cforms.as_mut_ptr(),
            cv.as_mut_ptr(),
            &mut nc,
            &mut ec,
        );
    }
    rust::readkfile(
        &files.r,
        &mut rnames,
        &mut rforms,
        &mut rv,
        &mut nr,
        &mut er,
    )
    .unwrap();
    assert_eq!((nc, ec), (nr, er), "input:\n{text}");
    same(&cv, &rv, text);
    files.positions_match();
    for (ptrs, expected) in [(&cnames, &rnames), (&cforms, &rforms)] {
        for (&p, value) in ptrs.iter().zip(expected) {
            if p.is_null() {
                assert!(value.is_empty());
            } else {
                unsafe {
                    fields_match(&[p], std::slice::from_ref(value));
                }
            }
        }
    }
    // readkfile itself leaks Rid, nspot, and overwritten string pointers.
    // Free accessible unique allocations, except its static "null" literal.
    let literal_null = if rnames[0] == "null" && ec >= 0 {
        cnames[0]
    } else {
        null_mut()
    };
    let mut owned: Vec<*mut c_char> = cnames.into_iter().chain(cforms).collect();
    owned.sort_unstable();
    owned.dedup();
    for p in owned {
        if !p.is_null() && p != literal_null {
            unsafe {
                free(p.cast());
            }
        }
    }
    (nr, er)
}

#[test]
fn readkfile_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b04);
    let mut seen = [0usize; 14];
    for case in 0..4000 {
        let n = size(&mut rng, case);
        let mut names = Vec::new();
        let mut forms = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            names.push(if i > 0 && rng.next_u32() % 4 == 0 {
                if case % 2 == 0 {
                    "null_ignored".into()
                } else {
                    "".into()
                }
            } else {
                // Distinct nonempty names differ inside allocated bytes:
                // strcmp need not read C's missing string terminators.
                format!("n{i:04}_{:08x}", rng.next_u32())
            });
            forms.push(["%.8e", "%12.6g", "%+.10f"][rng.next_u32() as usize % 3].to_owned());
            values.push(number(&mut rng));
        }
        let mode = case % 16;
        if mode == 7 {
            forms.push("%.3g".into());
        }
        let names_line = format!("C_SU_ID,{}\n", names.join(","));
        let forms_line = format!("C_SU_ID,{}\n", forms.join(","));
        let mut text = String::new();
        if mode != 6 {
            text.push_str("C_SU_SETID,K\n");
        }
        if mode == 5 {
            text.push_str("C_SU_SETID,K\n");
        }
        if mode != 2 {
            text.push_str("C_SU_NAMES\n");
            if mode == 1 {
                text.push_str("C_SU_NAMES\n");
            }
            text.push_str(&names_line);
        }
        if mode != 4 {
            text.push_str("C_SU_FORMS\n");
            if mode == 3 {
                text.push_str("C_SU_FORMS\n");
            }
            text.push_str(&forms_line);
        }
        // Metadata after both lists was intentionally not counted in C.
        if mode == 8 {
            text.push_str("C_SU_NAMES\nignored_header_payload\n");
        }
        text.push_str("X,other_record\n");
        if mode != 9 {
            text.push_str(if mode == 10 { "k," } else { "K," });
            if mode == 11 {
                text.push_str("bad,1 2");
            } else if mode == 12 {
                text.push_str("1 2");
            } else if mode == 13 {
                // Numeric data for the ID column only: insufficient commas.
                text.pop();
            } else if mode == 14 {
                text.push_str(&vec![""; n].join(","));
            } else {
                text.push_str(&values.join(","));
            }
            text.push('\n');
            text.push_str("K,9999\nunconsumed_tail\n");
        }
        // Case and whitespace are removed from metadata, never from K records.
        if case % 5 == 0 {
            let split = text.find("X,other_record").unwrap();
            let metadata = text[..split].replace("C_SU_", "c_sU_").replace(',', " ,\t");
            text = metadata + &text[split..];
        }
        if case % 7 == 0 {
            text = text.replace('\n', "\r\n");
        }
        let (_, error) = compare_read(&text, n + 16);
        seen[if error == -1 { 13 } else { error as usize }] += 1;
    }
    for code in [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13] {
        assert!(seen[code] > 0, "unexercised error code {code}");
    }
}

#[test]
fn readkfile_edge_records_match_c() {
    let _guard = C_LOCK.lock().unwrap();
    let base = "C_SU_SETID,K\nC_SU_FORMS\nC_SU_ID,%.9g\nC_SU_NAMES\nC_SU_ID,abcdefgh\n";
    for tail in [
        "",
        "K,-0",
        "K,1.25\r",
        "K,bad\n",
        " K,1\n",
        "k,2\n",
        "K,1\nK,2\n",
        "C_SU_COMMENT\nK,-4\n",
        "C_SU_NAMES\nK,unreadable\nK,7\n",
    ] {
        compare_read(&(base.to_owned() + tail), 16);
    }
    // Empty and zero-active-field files avoid C's undefined nspot[0] access
    // by having no K record. Missing data leaves dfield unchanged, as in C.
    compare_read("", 16);
    assert_eq!(compare_read(
        "C_SU_SETID,K\nC_SU_NAMES\nC_SU_ID,duplicate_value,duplicate_value\nC_SU_FORMS\nC_SU_ID,%.6g,%.6g\nK,1,2\n", 16), (-123, 9));
    assert_eq!(
        compare_read(
            "C_SU_SETID,K\nC_SU_NAMES\nC_SU_ID,nullskip\nC_SU_FORMS\nC_SU_ID,%.6g\n",
            16
        ),
        (0, 0)
    );
    // Exercise the fixed fgets boundary with long irrelevant records.
    for n in [9999, 10000, 10001, 20001] {
        compare_read(&(base.to_owned() + &"X".repeat(n) + "\nK,7\n"), 16);
    }
    compare_read(
        "C_SU_SETID,K\nC_SU_NAMES\nC_SU_ID,Aé中😀\nC_SU_FORMS\nC_SU_ID,%.9g\nK,-7\n",
        2,
    );
}

#[test]
fn readkfile_invalid_utf8_returns_error() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b08);
    let malformed: &[&[u8]] = &[
        b"\x80",
        b"\xff",
        b"\xc0\xaf",
        b"\xc2",
        b"\xe2\x82",
        b"\xed\xa0\x80",
        b"\xf0\x9f\x92",
        b"\xf4\x90\x80\x80",
    ];
    for case in 0..2048 {
        let output = if case % 2 == 0 { "names" } else { "forms" };
        let bad = if case < malformed.len() * 2 {
            malformed[case / 2].to_vec()
        } else {
            // Any isolated high byte followed by ASCII is invalid UTF-8.
            vec![0x80 | (rng.next_u32() as u8 & 0x7f)]
        };
        let prefix = format!("C_SU_ID,field_{:08x}", rng.next_u32());
        let mut bad_line = prefix.as_bytes().to_vec();
        bad_line.extend_from_slice(&bad);
        if case % 3 == 0 {
            bad_line.extend_from_slice(b"tail");
        }
        bad_line.push(b'\n');
        let mut text = b"C_SU_SETID,K\n".to_vec();
        for list in ["names", "forms"] {
            text.extend_from_slice(format!("C_SU_{}\n", list.to_ascii_uppercase()).as_bytes());
            text.extend_from_slice(if list == output {
                &bad_line
            } else if list == "names" {
                b"C_SU_ID,good_name\n"
            } else {
                b"C_SU_ID,%.9g\n"
            });
        }
        text.extend_from_slice(b"K,-1\n");
        let files = Files::from_bytes(&text);
        let mut names = vec!["unchanged".to_owned(); 2];
        let mut forms = names.clone();
        let mut values = [-6789.5];
        let (mut count, mut warning) = (-123, 987);
        let error = rust::readkfile(
            &files.r,
            &mut names,
            &mut forms,
            &mut values,
            &mut count,
            &mut warning,
        )
        .unwrap_err();
        match error {
            rust::ReadKFileError::InvalidUtf8 {
                output: actual,
                source,
            } => {
                assert_eq!(actual, output);
                assert_eq!(source.valid_up_to(), prefix.len());
            }
            other => panic!("unexpected error: {other}"),
        }
        assert_eq!((count, warning), (-123, 0));
        assert_eq!(values, [-6789.5]);
        assert_eq!(cio::eftello(&files.r), 0);
        let unchanged = if output == "names" { &names } else { &forms };
        assert!(unchanged.iter().all(|s| s == "unchanged"));
    }
}

#[test]
fn readkfile_capacity_errors_can_be_retried() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b09);
    for case in 0..1024 {
        let n = size(&mut rng, case);
        let names: Vec<_> = (0..n)
            .map(|i| {
                if i > 0 && rng.next_u32() % 3 == 0 {
                    "null_ignored".to_owned()
                } else {
                    format!("field_{i:04}")
                }
            })
            .collect();
        let values: Vec<_> = (0..n)
            .map(|i| {
                if i % 3 == 0 {
                    -0.0
                } else {
                    rng.range_i32(-10000, 10000) as f64
                }
            })
            .collect();
        let text = format!(
            "C_SU_SETID,K\nC_SU_NAMES\nC_SU_ID,{}\nC_SU_FORMS\nC_SU_ID,{}\nK,{}\n",
            names.join(","),
            vec!["%.9g"; n].join(","),
            values
                .iter()
                .map(f64::to_string)
                .collect::<Vec<_>>()
                .join(","),
        );
        let expected_names: Vec<_> = names
            .iter()
            .filter(|s| !s.starts_with("null"))
            .cloned()
            .collect();
        let expected_values: Vec<_> = names
            .iter()
            .zip(&values)
            .filter(|(s, _)| !s.starts_with("null"))
            .map(|(_, &value)| value)
            .collect();
        for output in ["names", "forms", "dfield"] {
            let files = Files::new(&text);
            let required = if output == "dfield" {
                expected_names.len()
            } else {
                n + 1 // The C_SU_ID and ignored null fields also need space.
            };
            let available = if case % 2 == 0 {
                0
            } else {
                rng.next_u32() as usize % required
            };
            let mut rnames = vec!["unchanged".to_owned(); n + 1];
            let mut rforms = rnames.clone();
            let mut rvalues = vec![-6789.5; expected_names.len()];
            match output {
                "names" => rnames.truncate(available),
                "forms" => rforms.truncate(available),
                _ => rvalues.truncate(available),
            }
            let (mut count, mut warning) = (-123, 987);
            assert_eq!(
                rust::readkfile(
                    &files.r,
                    &mut rnames,
                    &mut rforms,
                    &mut rvalues,
                    &mut count,
                    &mut warning,
                ),
                Err(rust::ReadKFileError::InsufficientCapacity {
                    output,
                    required,
                    available,
                }),
            );
            assert_eq!((count, warning), (-123, 0));
            assert!(rvalues.iter().all(|&value| value == -6789.5));
            assert_eq!(cio::eftello(&files.r), 0);
            match output {
                "names" => assert!(rnames.iter().all(|s| s == "unchanged")),
                "forms" => assert!(rforms.iter().all(|s| s == "unchanged")),
                _ => {}
            }
            // Resize only the output identified by the error, then retry
            // directly: readkfile has already rewound the stream.
            match output {
                "names" => rnames.resize(required, String::new()),
                "forms" => rforms.resize(required, String::new()),
                _ => rvalues.resize(required, -6789.5),
            }
            rust::readkfile(
                &files.r,
                &mut rnames,
                &mut rforms,
                &mut rvalues,
                &mut count,
                &mut warning,
            )
            .unwrap();
            assert_eq!((count, warning), (expected_names.len() as i32, 0));
            assert_eq!(&rnames[..count as usize], &expected_names);
            assert!(rforms[..count as usize].iter().all(|s| s == "%.9g"));
            same(&rvalues, &expected_values, &text);
        }
    }
}

fn compare_write(names: &[String], forms: &[String], values: &[f64]) {
    let files = Files::new("");
    let cnames: Vec<_> = names
        .iter()
        .map(|s| CString::new(s.split('\0').next().unwrap()).unwrap())
        .collect();
    let cforms: Vec<_> = forms
        .iter()
        .map(|s| CString::new(s.split('\0').next().unwrap()).unwrap())
        .collect();
    let mut np: Vec<_> = cnames.iter().map(|s| s.as_ptr() as *mut c_char).collect();
    let mut fp: Vec<_> = cforms.iter().map(|s| s.as_ptr() as *mut c_char).collect();
    let mut cv = values.to_vec();
    let (mut ce, mut re) = (-37, -37);
    unsafe {
        writekfile(
            files.c,
            np.as_mut_ptr(),
            fp.as_mut_ptr(),
            cv.as_mut_ptr(),
            values.len() as i32,
            &mut ce,
        );
    }
    rust::writekfile(&files.r, names, forms, values, &mut re);
    assert_eq!(ce, re);
    assert_eq!(ce, 0);
    let (c, r) = files.contents();
    assert_eq!(
        String::from_utf8(c).unwrap(),
        String::from_utf8(r).unwrap(),
        "names {names:?}, forms {forms:?}, values {values:?}"
    );
}

#[test]
fn writekfile_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b05);
    let specs = [
        "%.6g",
        "%+.8e",
        "% 15.7f",
        "%020.10E",
        "%-18.9G",
        "%#.0f",
        "%.17g",
        "%.3lf",
        "value=%.5e%%",
        "%.0f",
        "%#14.6G",
        "%.20e",
        "%a",
        "%+.5A",
        "%#20.0a",
        "%.13la",
        "literal%%",
        "%1$.5e/%1$.8g",
        "%'15.4f",
    ];
    let edges = [
        0.,
        -0.,
        1.,
        -1.,
        1e-300,
        -1e-300,
        1e307,
        -1e307,
        1e308,
        1.1e308,
        -1.1e308,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        f64::from_bits(1),
        f64::from_bits(1e308f64.to_bits() - 1),
    ];
    compare_write(&[], &[], &[]);
    compare_write(
        &vec!["x".into(); 9],
        &vec!["%.0a".into(); 9],
        &[1.49, 1.5, 1.51, 1.55, 1.6, 1.65, 1.75, 1.9, 1.99],
    );
    for x in [1.5, 1.75, 1.99, 3.5, 50382093.59389895] {
        compare_write(&["x".into()], &["%#.0a".into()], &[x]);
    }
    for x in [
        1.99,
        1.9999999999999998,
        0.9999999999999999,
        f64::from_bits(1),
        f64::MIN_POSITIVE,
    ] {
        for spec in ["%.1a", "%.5a", "%a", "%.13a", "%020.0a"] {
            compare_write(&["x".into()], &[spec.into()], &[x]);
        }
    }
    for case in 0..3000 {
        let n = size(&mut rng, case);
        let mut names = Vec::new();
        let mut forms = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            names.push(if case % 11 == 0 {
                "".to_owned()
            } else {
                format!("field_{i}")
            });
            let value = if i < edges.len() && case < 20 {
                edges[i]
            } else if rng.next_u32() % 3 == 0 {
                edges[rng.next_u32() as usize % edges.len()]
            } else {
                rng.uniform(-1e8, 1e8)
            };
            // Decimal fixed-point output of 257 huge negative doubles would
            // overrun C's 10001-byte buffer. Choose scientific format there.
            let form = if value.abs() > 1e100 {
                "%.17e"
            } else {
                specs[rng.next_u32() as usize % specs.len()]
            };
            forms.push(form.to_owned());
            values.push(value);
        }
        if case % 19 == 0 {
            names[0].push_str("\0ignored");
            forms[0].push_str("\0ignored");
        }
        compare_write(&names, &forms, &values);
    }
    for &x in &edges {
        for spec in specs {
            compare_write(&["x".into()], &[spec.into()], &[x]);
        }
    }
}

#[test]
fn writekfile_hex_rounding_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x4b07);
    for case in 0..4000 {
        let x = f64::from_bits(rng.next_u64());
        let form = format!(
            "%{}{}.{}{}",
            ["", "+", " ", "#", "0", "-", "+#0"][case % 7],
            [0, 1, 17, 32, 64][case % 5],
            case % 18,
            if case % 2 == 0 { "a" } else { "A" }
        );
        compare_write(&["x".into()], &[form], &[x]);
    }
    // Exact halfway cases and adjacent bit patterns exercise every rounding
    // position that has at least two omitted nibbles in the C formatter.
    for p in 0..12 {
        let shift = 52 - 4 * p;
        let halfway = (1u64 << 52) | (8u64 << (shift - 4));
        for bits in [halfway - 1, halfway, halfway + 1] {
            let x = f64::from_bits((1023u64 << 52) | (bits & ((1u64 << 52) - 1)));
            compare_write(&["x".into()], &[format!("%.{p}a")], &[x]);
        }
    }
}
