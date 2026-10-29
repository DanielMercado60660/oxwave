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

//! Differential tests of su/lib/fgetgthr.c and su/lib/fputgthr.c.
//! Credits: Potash Corporation Sascatchewan, Balasz Nemeth; given to CWP in 2008.

use oxwave_core::prelude::*;
use oxwave_cref::XorShift;
use oxwave_gather_io::{self as rust, GatherState};
use std::ffi::{c_char, c_void, CStr, CString};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::Mutex;

// C's segy rows have a 240-byte header followed by nt floats. These opaque
// pointers deliberately do not pretend the compact allocation is a full segy.
type Segy = c_void;
type File = c_void;
#[repr(C, align(8))]
struct CValue([u8; 8]);

extern "C" {
    // segy **fget_gather(FILE *, cwp_String *, cwp_String *, Value *,
    //                   int *, int *, float *, int *);
    fn fget_gather(
        fp: *mut File,
        key: *mut *mut c_char,
        type_: *mut *mut c_char,
        val: *mut CValue,
        nt: *mut i32,
        ntr: *mut i32,
        dt: *mut f32,
        first: *mut i32,
    ) -> *mut *mut Segy;
    // segy **get_gather(cwp_String *, cwp_String *, Value *,
    //                  int *, int *, float *, int *);
    fn get_gather(
        key: *mut *mut c_char,
        type_: *mut *mut c_char,
        val: *mut CValue,
        nt: *mut i32,
        ntr: *mut i32,
        dt: *mut f32,
        first: *mut i32,
    ) -> *mut *mut Segy;
    // segy **fput_gather(FILE *, segy **, int *, int *);
    fn fput_gather(
        fp: *mut File,
        rec: *mut *mut Segy,
        nt: *mut i32,
        ntr: *mut i32,
    ) -> *mut *mut Segy;
    // segy **put_gather(segy **, int *, int *);
    fn put_gather(rec: *mut *mut Segy, nt: *mut i32, ntr: *mut i32) -> *mut *mut Segy;
    fn tmpfile() -> *mut File;
    fn fclose(fp: *mut File) -> i32;
    fn rewind(fp: *mut File);
    fn fflush(fp: *mut File) -> i32;
    fn fread(p: *mut c_void, size: usize, n: usize, fp: *mut File) -> usize;
    fn fwrite(p: *const c_void, size: usize, n: usize, fp: *mut File) -> usize;
    fn malloc(n: usize) -> *mut c_void;
    fn free(p: *mut c_void);
}

static C_LOCK: Mutex<()> = Mutex::new(());
const SIZES: &[usize] = &[1, 2, 3, 7, 16, 17, 31, 64, 257, 1024, 4096, 65535];
const KEYS: &[&str] = &["cdp", "offset", "mark", "dt", "d1", "fldr"];
const ROUNDS: usize = 2200;

struct TempC(*mut File);
impl TempC {
    fn new() -> Self {
        let p = unsafe { tmpfile() };
        assert!(!p.is_null());
        Self(p)
    }
    fn bytes(&self) -> Vec<u8> {
        unsafe { rewind(self.0) };
        let mut result = Vec::new();
        let mut buf = [0u8; 8192];
        loop {
            let n = unsafe { fread(buf.as_mut_ptr().cast(), 1, buf.len(), self.0) };
            if n == 0 {
                break;
            }
            result.extend_from_slice(&buf[..n]);
        }
        result
    }
}
impl Drop for TempC {
    fn drop(&mut self) {
        assert_eq!(unsafe { fclose(self.0) }, 0);
    }
}

fn key_value(key: &str, group: usize) -> Value {
    let mut v = Value::default();
    // Include zero, negative keys, extremes, and a repeated (nonadjacent) key.
    let n = match group % 7 {
        0 => 0,
        1 => -1,
        2 => i32::MIN,
        3 => i32::MAX,
        _ => group as i32 * -71,
    };
    match hdtype(key) {
        "i" => v.set_i(n),
        "h" => v.set_h(n as i16),
        "u" => v.set_u(n as u16),
        "f" => v.set_f(match group % 9 {
            0 => -0.0,
            1 => 0.0,                        // Equal keys must merge.
            2 => f32::from_bits(0x7fc01234), // C valcmp treats NaN as equal.
            3 => f32::INFINITY,
            4 => f32::NEG_INFINITY,
            _ => n as f32 / 16.0,
        }),
        _ => unreachable!(),
    }
    v
}

fn trace(r: &mut XorShift, nt: usize, group: usize, serial: usize, key: &str) -> Trace {
    // Randomize every byte, including unassigned header words.
    let mut h = [0u8; HDRBYTES];
    for b in &mut h {
        *b = r.next_u32() as u8;
    }
    let mut hdr = SegyHeader::from_xdr(&h);
    hdr.ns = nt as u16;
    hdr.dt = r.next_u32() as u16;
    hdr.trid = 1;
    hdr.tracl = serial as i32;
    puthval(&mut hdr, getindex(key), &key_value(key, group));
    let mut data = vec![0.0; nt];
    for (i, x) in data.iter_mut().enumerate() {
        *x = match (i + serial) % 13 {
            0 => 0.0,
            1 => -0.0,
            2 => f32::MIN_POSITIVE,
            3 => f32::from_bits(1),
            4 => f32::MAX,
            5 => -f32::MAX,
            6 => f32::INFINITY,
            7 => f32::NEG_INFINITY,
            8 => f32::from_bits(0x7fc05678),
            _ => f32::from_bits(r.next_u32()),
        };
    }
    Trace { hdr, data }
}

fn encode_trace(out: &mut Vec<u8>, tr: &Trace) {
    out.extend_from_slice(oxwave_core::cio::as_bytes(std::slice::from_ref(&tr.hdr)));
    out.extend_from_slice(oxwave_core::cio::as_bytes(&tr.data));
}

// A corpus includes long and short gathers, lookahead transitions and a final
// singleton. Some datasets vary ns, testing zero-length later traces and the
// C reader's retention/truncation of samples to the first trace's nt.
fn corpus(nt: usize, key: &str, groups: usize) -> Vec<u8> {
    let mut r = XorShift::new(0x902d + nt as u64);
    let mut out = Vec::new();
    let mut serial = 0;
    for group in 0..groups {
        let ntr = if group + 1 == groups {
            1
        } else {
            match group {
                0 => 1,
                1 => 2,
                2 => 7,
                3 => 16,
                4 => 31,
                5 if nt == 1 => 4096,
                _ => 1 + r.next_u32() as usize % 5,
            }
        };
        for _ in 0..ntr {
            let ns = if serial > 0 && matches!(key, "offset" | "mark") {
                match serial % 5 {
                    0 => 0,
                    1 => 1,
                    2 => nt / 2,
                    3 => nt,
                    _ => (nt + 17).min(SU_NFLTS),
                }
            } else {
                nt
            };
            let tr = trace(&mut r, ns, group, serial, key);
            out.extend_from_slice(&tr.hdr.to_xdr());
            for x in tr.data {
                out.extend_from_slice(&x.to_bits().to_be_bytes());
            }
            serial += 1;
        }
    }
    out
}

fn read_cases(mut check: impl FnMut(&str, &[u8], &str)) {
    for (i, &nt) in SIZES.iter().enumerate() {
        let key = KEYS[i % KEYS.len()];
        check(
            key,
            &corpus(nt, key, if nt <= 64 { 301 } else { 3 }),
            &format!("nt={nt} key={key}"),
        );
    }
    // Distinguish the first singleton, a dropped final singleton, and a
    // returned final multi-trace gather. Cover the sampling-interval limits.
    for dt in [0, 1, u16::MAX] {
        for lengths in [&[1][..], &[2], &[1, 1], &[1, 2], &[2, 1], &[2, 2]] {
            let mut r = XorShift::new(0x51ab);
            let mut input = Vec::new();
            let mut serial = 0;
            for (group, &ntr) in lengths.iter().enumerate() {
                for _ in 0..ntr {
                    let mut tr = trace(&mut r, 2, group, serial, "cdp");
                    tr.dt = dt;
                    input.extend_from_slice(&tr.hdr.to_xdr());
                    for x in tr.data {
                        input.extend_from_slice(&x.to_bits().to_be_bytes());
                    }
                    serial += 1;
                }
            }
            check("cdp", &input, &format!("dt={dt} lengths={lengths:?}"));
        }
    }
}

// Include all out-parameters on every call, including the terminal None call.
fn metadata(
    out: &mut Vec<u8>,
    nt: i32,
    ntr: i32,
    dt: f32,
    first: i32,
    type_: &str,
    val: &[u8; 8],
    present: bool,
) {
    out.extend_from_slice(&nt.to_ne_bytes());
    out.extend_from_slice(&ntr.to_ne_bytes());
    out.extend_from_slice(&dt.to_bits().to_ne_bytes());
    out.extend_from_slice(&first.to_ne_bytes());
    out.extend_from_slice(type_.as_bytes());
    out.extend_from_slice(val);
    out.push(present as u8);
}

fn read_c(fp: Option<*mut File>, key: &str) -> Vec<u8> {
    let key = CString::new(key).unwrap();
    let mut kp = key.as_ptr() as *mut c_char;
    let mut tp = std::ptr::null_mut();
    let mut val = CValue([0xa5; 8]);
    let (mut nt, mut ntr, mut dt, mut first) = (-1, -1, -1.0, 0);
    let mut out = Vec::new();
    loop {
        let rec = unsafe {
            match fp {
                Some(fp) => fget_gather(
                    fp, &mut kp, &mut tp, &mut val, &mut nt, &mut ntr, &mut dt, &mut first,
                ),
                None => get_gather(
                    &mut kp, &mut tp, &mut val, &mut nt, &mut ntr, &mut dt, &mut first,
                ),
            }
        };
        let type_ = unsafe { CStr::from_ptr(tp) }.to_str().unwrap();
        metadata(&mut out, nt, ntr, dt, first, type_, &val.0, !rec.is_null());
        if rec.is_null() {
            break;
        }
        for i in 0..ntr as usize {
            unsafe {
                let row = *rec.add(i);
                out.extend_from_slice(std::slice::from_raw_parts(
                    row.cast::<u8>(),
                    HDRBYTES + nt as usize * 4,
                ));
                free(row);
            }
        }
        unsafe { free(rec.cast()) };
    }
    out
}

fn read_rust(fp: Option<&CFile>, key: &str) -> Vec<u8> {
    let mut state = GatherState::default();
    let mut val = Value { bytes: [0xa5; 8] };
    let mut type_ = "";
    let (mut nt, mut ntr, mut dt, mut first) = (-1, -1, -1.0, 0);
    let mut out = Vec::new();
    loop {
        let rec = match fp {
            Some(fp) => rust::fget_gather(
                fp, &mut state, key, &mut type_, &mut val, &mut nt, &mut ntr, &mut dt, &mut first,
            ),
            None => rust::get_gather(
                &mut state, key, &mut type_, &mut val, &mut nt, &mut ntr, &mut dt, &mut first,
            ),
        };
        metadata(
            &mut out,
            nt,
            ntr,
            dt,
            first,
            type_,
            &val.bytes,
            rec.is_some(),
        );
        let Some(rec) = rec else { break };
        assert_eq!(rec.len(), ntr as usize);
        for tr in rec {
            encode_trace(&mut out, &tr);
        }
    }
    out
}

fn exact(a: &[u8], b: &[u8], context: &str) {
    assert_eq!(a.len(), b.len(), "{context}: lengths");
    assert!(
        a == b,
        "{context}: first differing byte {:?}",
        a.iter().zip(b).position(|(a, b)| a != b)
    );
}

#[test]
fn fget_gather_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    read_cases(|key, input, context| {
        let cf = TempC::new();
        assert_eq!(
            unsafe { fwrite(input.as_ptr().cast(), 1, input.len(), cf.0) },
            input.len()
        );
        unsafe { rewind(cf.0) };
        let rf = etmpfile();
        efwrite(input, &rf);
        erewind(&rf);
        exact(
            &read_rust(Some(&rf), key),
            &read_c(Some(cf.0), key),
            &format!("fget_gather {context}"),
        );
    });
}

// Allocate with libc because the original writers free every row.
fn c_rows(rec: &[Trace], nt: usize) -> Vec<*mut Segy> {
    rec.iter()
        .map(|tr| {
            let mut bytes = Vec::new();
            encode_trace(&mut bytes, tr);
            assert_eq!(bytes.len(), HDRBYTES + nt * 4);
            unsafe {
                let p = malloc(bytes.len());
                assert!(!p.is_null());
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), p.cast(), bytes.len());
                p
            }
        })
        .collect()
}

fn write_cases(mut emit: impl FnMut(&[Trace], i32)) {
    let mut r = XorShift::new(0xaba347);
    emit(&[], 0);
    emit(&[trace(&mut r, 0, 0, 0, "cdp")], 0);
    for i in 0..ROUNDS {
        let nt = if i < SIZES.len() {
            SIZES[i]
        } else {
            SIZES[r.next_u32() as usize % 9]
        };
        let ntr = match i {
            0 => 0,
            1 => 1,
            2 => 2,
            3 => 7,
            4 => 16,
            5 => 31,
            6 => 257,
            7 => 1024,
            8..=11 => 1,
            _ => r.next_u32() as usize % 6,
        };
        let mut rec = Vec::new();
        for j in 0..ntr {
            let mut tr = trace(&mut r, nt, i, j, "cdp");
            // Variable ns, including zero; the copy length remains nt.
            tr.ns = if i % 3 == 0 {
                (r.next_u32() as usize % (nt + 1)) as u16
            } else {
                nt as u16
            };
            rec.push(tr);
        }
        emit(&rec, nt as i32);
    }
}

#[test]
fn fput_gather_matches_c() {
    let _lock = C_LOCK.lock().unwrap();
    let cf = TempC::new();
    let rf = etmpfile();
    write_cases(|rec, mut nt| {
        let saved = rec.to_vec();
        rust::fput_gather(&rf, rec, nt);
        let mut ptrs = c_rows(rec, nt as usize);
        let mut ntr = ptrs.len() as i32;
        assert!(unsafe { fput_gather(cf.0, ptrs.as_mut_ptr(), &mut nt, &mut ntr) }.is_null());
        assert_eq!(ntr, rec.len() as i32);
        // Borrowing must not mutate even NaNs or signed zero.
        let mut before = Vec::new();
        let mut after = Vec::new();
        for tr in &saved {
            encode_trace(&mut before, tr);
        }
        for tr in rec {
            encode_trace(&mut after, tr);
        }
        exact(&before, &after, "borrowed input");
    });
    erewind(&rf);
    let mut got = Vec::new();
    let mut buf = [0u8; 8192];
    loop {
        let n = oxwave_core::cio::fread(&mut buf, &rf);
        if n == 0 {
            break;
        }
        got.extend_from_slice(&buf[..n]);
    }
    exact(&got, &cf.bytes(), "fput_gather stream");
}

fn artifact(name: &str) -> PathBuf {
    PathBuf::from(std::env::var_os("CARGO_TARGET_DIR").unwrap_or_else(|| {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../.targets/lib-gather_io")
            .into_os_string()
    }))
    .join(format!("gather-diff-{}-{name}", std::process::id()))
}

fn child(test: &str, mode: &str, input: Option<&PathBuf>, output: &PathBuf, key: &str) {
    let mut cmd = Command::new(std::env::current_exe().unwrap());
    cmd.args(["--exact", test, "--nocapture"])
        .env("OXWAVE_GATHER_MODE", mode)
        .env("OXWAVE_GATHER_OUT", output)
        .env("OXWAVE_GATHER_KEY", key)
        .env("TMPDIR", artifact("tmp"))
        .stdout(Stdio::null())
        .stderr(Stdio::piped());
    cmd.stdin(match input {
        Some(p) => Stdio::from(std::fs::File::open(p).unwrap()),
        None => Stdio::null(),
    });
    let result = cmd.output().unwrap();
    assert!(
        result.status.success(),
        "{test} {mode}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn get_gather_matches_c() {
    if let Ok(mode) = std::env::var("OXWAVE_GATHER_MODE") {
        let key = std::env::var("OXWAVE_GATHER_KEY").unwrap();
        let bytes = if mode == "c" {
            read_c(None, &key)
        } else {
            read_rust(None, &key)
        };
        std::fs::write(std::env::var_os("OXWAVE_GATHER_OUT").unwrap(), bytes).unwrap();
        return;
    }
    let _lock = C_LOCK.lock().unwrap();
    let tmp = artifact("tmp");
    std::fs::create_dir_all(&tmp).unwrap();
    let input = artifact("input");
    let c_out = artifact("c-read");
    let r_out = artifact("r-read");
    read_cases(|key, bytes, context| {
        std::fs::write(&input, bytes).unwrap();
        child("get_gather_matches_c", "c", Some(&input), &c_out, key);
        child("get_gather_matches_c", "rust", Some(&input), &r_out, key);
        exact(
            &std::fs::read(&r_out).unwrap(),
            &std::fs::read(&c_out).unwrap(),
            &format!("get_gather {context}"),
        );
    });
    for p in [input, c_out, r_out] {
        std::fs::remove_file(p).unwrap();
    }
    std::fs::remove_dir_all(tmp).unwrap();
}

#[test]
fn readers_reject_invalid_input_like_c() {
    if let Ok(mode) = std::env::var("OXWAVE_GATHER_MODE") {
        // Both error packages prepend argv[0]; keep it identical in the child.
        let mut argv = oxwave_cref::CArgv::new(&["gather-diff"]);
        unsafe { oxwave_cref::initargs(argv.argc(), argv.argv()) };
        initargs(&["gather-diff".to_owned()]);
        let path = std::env::var_os("OXWAVE_GATHER_INPUT").unwrap();
        match mode.as_str() {
            "c-file" => {
                let bytes = std::fs::read(path).unwrap();
                let fp = TempC::new();
                unsafe {
                    assert_eq!(
                        fwrite(bytes.as_ptr().cast(), 1, bytes.len(), fp.0),
                        bytes.len()
                    );
                    rewind(fp.0);
                }
                read_c(Some(fp.0), "cdp");
            }
            "rust-file" => {
                let fp = efopen(path.to_str().unwrap(), "r");
                read_rust(Some(&fp), "cdp");
            }
            "c-stdin" => {
                read_c(None, "cdp");
            }
            "rust-stdin" => {
                read_rust(None, "cdp");
            }
            _ => panic!("unknown child mode"),
        }
        panic!("invalid input unexpectedly succeeded");
    }
    let _lock = C_LOCK.lock().unwrap();
    let input = artifact("invalid-input");
    let tmp = artifact("tmp");
    std::fs::create_dir_all(&tmp).unwrap();
    let mut hdr = SegyHeader::zeroed();
    hdr.trid = 1;
    hdr.ns = 2;
    let mut short_data = hdr.to_xdr().to_vec();
    short_data.extend_from_slice(&1.0f32.to_bits().to_be_bytes());
    let missing_data = hdr.to_xdr().to_vec();
    hdr.ns = 0;
    for bytes in [
        vec![],
        vec![0; HDRBYTES - 1],
        hdr.to_xdr().to_vec(),
        missing_data,
        short_data,
    ] {
        std::fs::write(&input, bytes).unwrap();
        for stream in ["file", "stdin"] {
            let run = |implementation| {
                Command::new(std::env::current_exe().unwrap())
                    .args([
                        "--exact",
                        "readers_reject_invalid_input_like_c",
                        "--nocapture",
                    ])
                    .env("OXWAVE_GATHER_MODE", format!("{implementation}-{stream}"))
                    .env("OXWAVE_GATHER_INPUT", &input)
                    .env("TMPDIR", &tmp)
                    .stdin(Stdio::from(std::fs::File::open(&input).unwrap()))
                    .output()
                    .unwrap()
            };
            let c = run("c");
            let rust = run("rust");
            assert_eq!(c.status.code(), Some(1));
            assert_eq!(rust.status.code(), c.status.code());
            exact(&rust.stderr, &c.stderr, "reader error message");
            exact(&rust.stdout, &c.stdout, "reader error stdout");
        }
    }
    std::fs::remove_file(input).unwrap();
    std::fs::remove_dir_all(tmp).unwrap();
}

#[test]
fn put_gather_matches_c() {
    if let Ok(mode) = std::env::var("OXWAVE_GATHER_MODE") {
        write_cases(|rec, mut nt| {
            if mode == "c" {
                let mut ptrs = c_rows(rec, nt as usize);
                let mut ntr = ptrs.len() as i32;
                assert!(unsafe { put_gather(ptrs.as_mut_ptr(), &mut nt, &mut ntr) }.is_null());
            } else {
                rust::put_gather(rec, nt);
            }
        });
        assert_eq!(unsafe { fflush(std::ptr::null_mut()) }, 0);
        efflush(&stdout());
        // Avoid mixing the test harness's result text with binary SU output.
        std::process::exit(0);
    }
    let _lock = C_LOCK.lock().unwrap();
    // Capture fd 1 directly, removing the harness's fixed initial banner.
    let tmp = artifact("tmp");
    std::fs::create_dir_all(&tmp).unwrap();
    let run = |mode: &str| {
        let output = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "put_gather_matches_c", "--nocapture", "--quiet"])
            .env("OXWAVE_GATHER_MODE", mode)
            .env("TMPDIR", &tmp)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let start = b"\nrunning 1 test\n";
        assert!(output.stdout.starts_with(start));
        output.stdout[start.len()..].to_vec()
    };
    exact(&run("rust"), &run("c"), "put_gather stdout");
    std::fs::remove_dir_all(tmp).unwrap();
}
