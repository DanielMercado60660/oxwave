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

//! Differential tests for su/lib/gridread.c.
//! Credits: Copyright (c) Colorado School of Mines, 2021; no named author.

use oxwave_cref::{CArgv, XorShift};
use std::ffi::{c_char, c_void, CString};
use std::sync::Mutex;

static LOCK: Mutex<()> = Mutex::new(());
const KEYS: [&str; 8] = ["grid_xa", "grid_ya", "grid_xb", "grid_yb",
    "grid_xc", "grid_yc", "grid_wb", "grid_wc"];

extern "C" {
    fn gridcommand(enough: *mut i32);
    fn gridread(fp: *mut c_void, gvals: *mut f64, errwarn: *mut i32);
    fn fopen(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn fclose(fp: *mut c_void) -> i32;
    fn ftello(fp: *mut c_void) -> i64;
}

fn setup(args: &[String]) -> CArgv {
    let refs: Vec<&str> = args.iter().map(String::as_str).collect();
    let mut c = CArgv::new(&refs);
    unsafe { oxwave_cref::initargs(c.argc(), c.argv()); }
    oxwave_core::par::initargs(args);
    c
}

#[test]
fn gridcommand_matches_c() {
    let _guard = LOCK.lock().unwrap();
    let mut rng = XorShift::new(123);
    for case in 0..4096 {
        let mask = if case < 256 { case } else { rng.next_u32() as usize & 255 };
        let mut args = vec!["gridtest".to_string(), "grid_xd=5".to_string()];
        for (i, key) in KEYS.iter().enumerate() {
            if mask & (1 << i) != 0 {
                let value = match case % 5 { 0 => 0.0, 1 => -0.0, _ => rng.uniform(-1e6, 1e6) };
                args.push(format!("{key}={value:.17e}"));
            }
        }
        let _argv = setup(&args);
        let (mut c, mut rust) = (77, 77);
        unsafe { gridcommand(&mut c); }
        oxwave_grid_parameters::gridcommand(&mut rust);
        assert_eq!(rust, c, "case {case} mask {mask}");
    }
}

#[test]
fn gridread_matches_c() {
    let _guard = LOCK.lock().unwrap();
    let mut rng = XorShift::new(8127);
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("gridread-diff-{}.k", std::process::id()));
    let cpath = CString::new(path.to_str().unwrap()).unwrap();
    let mode = CString::new("r").unwrap();
    for case in 0..3000 {
        let mask = if case < 256 { case } else { rng.next_u32() as usize & 255 };
        let mut args = vec!["gridtest".to_string()];
        for (i, key) in KEYS.iter().enumerate() {
            if mask & (1 << i) != 0 {
                let value = match case % 7 { 0 => 0.0, 1 => -0.0, 2 => -1.05e308,
                    _ => rng.uniform(-1e9, 1e9) };
                args.push(format!("{key}={value:.17e}"));
            }
        }
        let use_file = case % 3 != 0;
        if use_file {
            let count = match case % 8 { 0 => 1, 1 => 2, 2 => 3, 3 => 7, _ => 8 };
            let mut text = String::from("C_SU_SETID,K\nC_SU_FORMS\nC_SU_ID");
            for _ in 0..count { text.push_str(",%.17g"); }
            text.push_str("\nC_SU_NAMES\nC_SU_ID");
            for key in &KEYS[..count] { text.push(','); text.push_str(&key.to_ascii_uppercase()); }
            text.push_str("\nK");
            for i in 0..count {
                let value = if (case + i) % 11 == 0 { 0.0 } else { rng.uniform(-1e6, 1e6) };
                text.push_str(&format!(",{value:.17e}"));
            }
            text.push('\n');
            let header = "C_SU_SETID,K\nC_SU_FORMS\nC_SU_ID,%g,%g\nC_SU_NAMES\nC_SU_ID,grid_xa,grid_ya\n";
            match case % 31 {
                0 => text.clear(),
                1 => text = "C_SU_SETID,K\nC_SU_NAMES\nC_SU_ID,grid_xa\n".into(),
                2 => text = "C_SU_FORMS\nC_SU_ID,%g\nC_SU_NAMES\nC_SU_ID,grid_xa\nK,1\n".into(),
                3 => text = "C_SU_SETID,K\nC_SU_FORMS\nC_SU_ID,%g\nC_SU_NAMES\nC_SU_ID,grid_xa,grid_ya\nK,1,2\n".into(),
                4 => text = format!("{header}K, ,2\n"),
                5 => text = format!("{header}K,abc,2\n"),
                6 => text = format!("{header}K,1 2,2\n"),
                7 => text = format!("{header}K,1\n"),
                _ => {}
            }
            std::fs::write(&path, text).unwrap();
        }
        let _argv = setup(&args);
        let mut c = [42.25; 20];
        let mut rust = c;
        let (mut ce, mut re) = (0, 0);
        let cfp = if use_file { unsafe { fopen(cpath.as_ptr(), mode.as_ptr()) } } else { std::ptr::null_mut() };
        assert!(!use_file || !cfp.is_null());
        let rfp = use_file.then(|| oxwave_core::cio::fopen(path.to_str().unwrap(), "r").unwrap());
        unsafe { gridread(cfp, c.as_mut_ptr(), &mut ce); }
        oxwave_grid_parameters::gridread(rfp.as_ref(), &mut rust, &mut re);
        if use_file {
            unsafe { fclose(cfp); }
            oxwave_core::cio::fclose(rfp.as_ref().unwrap());
        }
        assert_eq!(re, ce, "error case {case}");
        for i in 0..20 {
            assert_eq!(rust[i].to_bits(), c[i].to_bits(), "case {case} index {i}");
        }
    }
    std::fs::remove_file(path).unwrap();
}

#[test]
fn gridread_byte_metadata_matches_c() {
    let _guard = LOCK.lock().unwrap();
    let mut rng = XorShift::new(0x85ff3127);
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("gridread-bytes-{}.k", std::process::id()));
    let cpath = CString::new(path.to_str().unwrap()).unwrap();
    let mode = CString::new("r").unwrap();
    for case in 0..2048 {
        let variant = case % 8;
        let mut extra = [1, 2, 3, 7, 8, 16, 31, 64, 127, 256, 509, 989, 990][case % 13];
        if variant == 3 || variant == 4 { extra = extra.max(2); }
        let mut names: Vec<Vec<u8>> = KEYS.iter().map(|key| key.as_bytes().to_vec()).collect();
        let mut forms = vec![b"%.17g".to_vec(); KEYS.len() + extra];
        for i in 0..extra {
            let mut name = format!("u{i:03x}").into_bytes();
            if variant != 0 {
                // Isolated high bytes cover every invalid UTF-8 byte value.
                name.push(0x80 | ((case + i) & 127) as u8);
            }
            names.push(name);
        }
        if variant != 1 {
            for form in &mut forms {
                form.push(0x80 | (rng.next_u32() & 127) as u8);
            }
        }
        let mut expected_error = 0;
        if variant == 3 {
            names[9] = names[8].clone();
            expected_error = 9;
        } else if variant == 4 {
            // Lossy decoding would collapse these distinct names and report 9.
            names[8] = b"u\x80".to_vec();
            names[9] = b"u\x81".to_vec();
        } else if variant == 5 {
            // C ignores every name beginning with "null", even non-UTF-8 ones.
            for name in &mut names[8..] { *name = b"null\xff".to_vec(); }
        } else if variant == 6 {
            names[0].push(0xff);
            expected_error = 20;
        }

        let mut args = vec!["gridtest".to_owned()];
        if variant == 7 {
            for key in KEYS {
                if rng.next_u32() & 1 != 0 {
                    args.push(format!("{key}={:.17e}", rng.uniform(-1e9, 1e9)));
                }
            }
        }
        let mut text = b"C_SU_SETID,K\n".to_vec();
        let lists = if case % 2 == 0 {
            [(b"C_SU_NAMES\n".as_slice(), &names), (b"C_SU_FORMS\n".as_slice(), &forms)]
        } else {
            [(b"C_SU_FORMS\n".as_slice(), &forms), (b"C_SU_NAMES\n".as_slice(), &names)]
        };
        for (heading, fields) in lists {
            text.extend_from_slice(heading);
            text.extend_from_slice(b"C_SU_ID");
            for field in fields {
                text.push(b',');
                text.extend_from_slice(field);
            }
            text.push(b'\n');
        }
        text.push(b'K');
        for i in 0..names.len() {
            text.push(b',');
            if i < KEYS.len() {
                let value = match (case + i) % 7 {
                    0 => 0.0, 1 => -0.0, _ => rng.uniform(-1e6, 1e6),
                };
                text.extend_from_slice(format!("{value:.17e}").as_bytes());
            } else if variant == 5 {
                text.extend_from_slice(b"bad");
            } else {
                text.extend_from_slice(b"-1");
            }
        }
        text.extend_from_slice(b"\nK,999\n");
        assert!(std::str::from_utf8(&text).is_err());
        std::fs::write(&path, text).unwrap();
        let _argv = setup(&args);
        let mut c = [42.25; 20];
        let mut rust = c;
        let (mut ce, mut re) = (0, 0);
        let cfp = unsafe { fopen(cpath.as_ptr(), mode.as_ptr()) };
        assert!(!cfp.is_null());
        let rfp = oxwave_core::cio::fopen(path.to_str().unwrap(), "r").unwrap();
        unsafe { gridread(cfp, c.as_mut_ptr(), &mut ce); }
        oxwave_grid_parameters::gridread(Some(&rfp), &mut rust, &mut re);
        let cpos = unsafe { ftello(cfp) };
        assert_eq!(oxwave_core::cio::ftello(&rfp), cpos, "position case {case}");
        unsafe { fclose(cfp); }
        oxwave_core::cio::fclose(&rfp);
        assert_eq!(ce, expected_error, "C error case {case}");
        assert_eq!(re, ce, "error case {case}");
        for i in 0..20 {
            assert_eq!(rust[i].to_bits(), c[i].to_bits(), "case {case} index {i}");
        }
    }
    std::fs::remove_file(path).unwrap();
}

#[test]
fn gridread_invalid_capacity_errors() {
    const CHILD: &str = "OXWAVE_GRIDREAD_CAPACITY_CHILD";
    if let Ok(kind) = std::env::var(CHILD) {
        oxwave_core::par::initargs(&["gridtest".to_owned()]);
        let mut values = vec![0.0; if kind == "short" { 19 } else { 20 }];
        let mut errwarn = 0;
        let file = std::env::var("OXWAVE_GRIDREAD_CAPACITY_FILE").ok()
            .map(|path| oxwave_core::cio::fopen(&path, "r").unwrap());
        oxwave_grid_parameters::gridread(file.as_ref(), &mut values, &mut errwarn);
        panic!("invalid gridread capacity should exit through err!");
    }
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("gridread-capacity-{}.k", std::process::id()));
    for kind in ["short", "names", "forms"] {
        let mut child = std::process::Command::new(std::env::current_exe().unwrap());
        child.args(["--exact", "gridread_invalid_capacity_errors", "--nocapture"])
            .env(CHILD, kind)
            .env_remove("OXWAVE_GRIDREAD_CAPACITY_FILE");
        let message = if kind == "short" {
            "gridread: gvals must contain at least 20 elements."
        } else {
            let header = if kind == "names" { "C_SU_NAMES" } else { "C_SU_FORMS" };
            let mut text = format!("C_SU_SETID,K\n{header}\nC_SU_ID");
            // 1,000 fields would overrun C's 999-element metadata arrays.
            for _ in 0..999 { text.push_str(",x"); }
            text.push('\n');
            std::fs::write(&path, text).unwrap();
            child.env("OXWAVE_GRIDREAD_CAPACITY_FILE", &path);
            "gridread: too many K-file metadata fields."
        };
        let output = child.output().unwrap();
        assert_eq!(output.status.code(), Some(1), "{kind}: {output:?}");
        assert_eq!(output.stderr, format!("\ngridtest: {message}\n").as_bytes());
    }
    std::fs::remove_file(path).unwrap();
}
