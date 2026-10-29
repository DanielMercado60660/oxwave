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

//! Differential tests for cwp/lib/tempname.c.
//! Credits: Andreas Klaedtke (2009); Jack K. Cohen, Colorado School of Mines (1995).

use oxwave_cref::XorShift;
use std::ffi::{c_char, c_int, c_void, CStr, CString};
use std::io::{Read, Seek, SeekFrom, Write};
use std::os::fd::{BorrowedFd, FromRawFd};
use std::os::unix::fs::MetadataExt;
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

extern "C" {
    fn temporary_filename(prefix: *mut c_char) -> *mut c_char;
    fn temporary_stream(prefix: *const c_char) -> *mut c_void;
    fn fwrite(ptr: *const c_void, size: usize, count: usize, stream: *mut c_void) -> usize;
    fn fread(ptr: *mut c_void, size: usize, count: usize, stream: *mut c_void) -> usize;
    fn fseek(stream: *mut c_void, offset: i64, whence: c_int) -> c_int;
    fn fclose(stream: *mut c_void) -> c_int;
}

#[test]
fn filename_matches_c_structure() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(345);
    let before = std::fs::read_dir("/dev/fd").unwrap().count();
    let mut rust_suffixes = std::collections::HashSet::new();
    for i in 0..2048 {
        let prefix = match i % 8 {
            0 => String::new(),
            1 => "/tmp/".to_string(),
            2 => format!("relative/{}", rng.next_u32()),
            3 => format!("-{}", rng.next_u32()),
            4 => "x".repeat((rng.next_u32() % 512) as usize),
            5 => "x".repeat([1, 2, 7, 127, 256, 509, 1000][i % 7]),
            6 => "directory\0ignored".to_owned(),
            _ => "séismique/地震".to_owned(),
        };
        let mut cprefix = prefix.as_bytes().to_vec();
        cprefix.push(0);
        let c = unsafe { CStr::from_ptr(temporary_filename(cprefix.as_mut_ptr().cast())) }
            .to_str()
            .unwrap()
            .to_owned();
        let rust = oxwave_tempfile::temporary_filename(&prefix);
        let prefix = prefix.split('\0').next().unwrap();
        for name in [&c, &rust] {
            let suffix = name.strip_prefix(prefix).unwrap();
            assert!(suffix.starts_with("/cmgui"));
            assert_eq!(suffix.len(), 12);
            assert!(suffix[6..].bytes().all(|b| b.is_ascii_alphanumeric()));
        }
        rust_suffixes.insert(rust.strip_prefix(prefix).unwrap().to_owned());

        // C leaks its mkstemp descriptor and file on success. Clean up only
        // that reference artifact, including on an unsandboxed test run.
        let path = format!("/tmp{}", c.strip_prefix(prefix).unwrap());
        if let Ok(metadata) = std::fs::metadata(&path) {
            let mut closed = false;
            for entry in std::fs::read_dir("/dev/fd").unwrap() {
                let entry = entry.unwrap();
                let fd = entry.file_name().to_str().unwrap().parse().unwrap();
                // Darwin's /dev/fd/N stat describes the fdescfs entry. Use
                // fstat through a duplicate to identify the actual open file.
                let borrowed = unsafe { BorrowedFd::borrow_raw(fd) };
                if let Ok(duplicate) = borrowed.try_clone_to_owned() {
                    let file = std::fs::File::from(duplicate);
                    let fd_metadata = file.metadata().unwrap();
                    if (fd_metadata.dev(), fd_metadata.ino()) == (metadata.dev(), metadata.ino()) {
                        drop(unsafe { std::fs::File::from_raw_fd(fd) });
                        closed = true;
                    }
                }
            }
            assert!(closed, "C's leaked descriptor must be found");
            std::fs::remove_file(path).unwrap();
        }
        assert_eq!(std::fs::read_dir("/dev/fd").unwrap().count(), before);
    }
    // Detect the previous fixed fallback without requiring random names to
    // be collision-free after their intermediate files have been removed.
    assert!(rust_suffixes.len() > 1);
}

#[test]
fn stream_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    let mut rng = XorShift::new(948);
    let dir = std::env::current_dir().unwrap();
    let dir = dir.to_str().unwrap();
    let before = std::fs::read_dir(dir).unwrap().count();
    for i in 0..2048 {
        let invalid = format!("{dir}/missing_{}", rng.next_u32());
        let prefix = match i % 8 {
            0 => invalid,
            1 => format!("{dir}/Cargo.toml"),
            2 => format!("{dir}/"),
            3 => format!("{dir}\0ignored/missing"),
            _ => dir.to_owned(),
        };
        let mut cp = prefix.as_bytes().to_vec();
        cp.push(0);
        let c = unsafe { temporary_stream(cp.as_ptr().cast()) };
        let rust = oxwave_tempfile::temporary_stream(Some(&prefix));
        assert_eq!(c.is_null(), rust.is_none());
        if let Some(mut rust) = rust {
            assert_eq!(rust.metadata().unwrap().mode() & 0o777, 0o600);
            assert_eq!(rust.metadata().unwrap().nlink(), 0);
            let len = [0, 1, 2, 7, 16, 127, 256, 1024, 8191, 65536][i % 10];
            let bytes: Vec<u8> = (0..len)
                .map(|_| if i % 11 == 0 { 0 } else { rng.next_u32() as u8 })
                .collect();
            unsafe {
                assert_eq!(fwrite(bytes.as_ptr().cast(), 1, len, c), len);
                assert_eq!(fseek(c, 0, 0), 0);
            }
            rust.write_all(&bytes).unwrap();
            rust.seek(SeekFrom::Start(0)).unwrap();
            let mut a = vec![0; len];
            let mut b = vec![0; len];
            assert_eq!(unsafe { fread(a.as_mut_ptr().cast(), 1, len, c) }, len);
            rust.read_exact(&mut b).unwrap();
            assert_eq!(a, b);
            assert_eq!(unsafe { fclose(c) }, 0);
        }
    }
    assert_eq!(std::fs::read_dir(dir).unwrap().count(), before);
    for prefix in [None, Some("")] {
        let cp = prefix.map(|s| CString::new(s).unwrap());
        let c = unsafe { temporary_stream(cp.as_ref().map_or(std::ptr::null(), |s| s.as_ptr())) };
        let r = oxwave_tempfile::temporary_stream(prefix);
        assert_eq!(c.is_null(), r.is_none());
        if !c.is_null() {
            assert_eq!(unsafe { fclose(c) }, 0);
        }
    }
}

#[test]
fn descriptor_exhaustion_matches_c() {
    let _guard = C_LOCK.lock().unwrap();
    const CHILD: &str = "OXWAVE_TEMPFILE_EXHAUST_FDS";
    if std::env::var_os(CHILD).is_none() {
        // Lower the limit only in a child, keeping the other tests and cargo
        // unaffected. Pass the executable as an argument, without interpolation.
        let output = std::process::Command::new("/bin/sh")
            .arg("-c")
            .arg("ulimit -n 64 && exec \"$1\" --exact descriptor_exhaustion_matches_c --nocapture")
            .arg("sh")
            .arg(std::env::current_exe().unwrap())
            .env(CHILD, "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "exhaustion child failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let directory = std::env::current_dir().unwrap();
    let directory = directory.to_str().unwrap();
    let c_directory = CString::new(directory).unwrap();
    let mut files = Vec::new();
    loop {
        match std::fs::File::open("/dev/null") {
            Ok(file) => files.push(file),
            Err(error) => {
                assert_eq!(error.raw_os_error(), Some(24)); // EMFILE on macOS/Linux.
                break;
            }
        }
    }

    let mut rng = XorShift::new(729);
    let mut c_suffixes = std::collections::HashSet::new();
    let mut rust_suffixes = std::collections::HashSet::new();
    for _ in 0..2048 {
        // A writable directory ensures this stream failure is caused by
        // descriptor exhaustion, even when the sandbox denies /tmp access.
        assert!(unsafe { temporary_stream(c_directory.as_ptr()) }.is_null());
        assert_eq!(std::io::Error::last_os_error().raw_os_error(), Some(24));
        assert!(oxwave_tempfile::temporary_stream(Some(directory)).is_none());

        let prefix = format!("relative/{}", rng.next_u32());
        let mut cprefix = prefix.as_bytes().to_vec();
        cprefix.push(0);
        let c = unsafe { CStr::from_ptr(temporary_filename(cprefix.as_mut_ptr().cast())) }
            .to_str()
            .unwrap()
            .to_owned();
        let rust = oxwave_tempfile::temporary_filename(&prefix);
        for name in [&c, &rust] {
            let suffix = name.strip_prefix(&prefix).unwrap();
            assert!(suffix.starts_with("/cmgui"));
            assert_eq!(suffix.len(), 12);
            assert!(suffix[6..].bytes().all(|b| b.is_ascii_alphanumeric()));
        }
        c_suffixes.insert(c.strip_prefix(&prefix).unwrap().to_owned());
        rust_suffixes.insert(rust.strip_prefix(&prefix).unwrap().to_owned());
    }
    // Follow libc's failure-template behavior, including whether failed calls
    // modify the template. Random suffix bytes cannot be compared directly.
    assert_eq!(c_suffixes.len() > 1, rust_suffixes.len() > 1);
    drop(files);
}
