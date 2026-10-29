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
//! Differential tests for par/lib/syscalls.c.
//! Authors: SEP: Rick Otolini, Ron, Jon Claerbout, Stew Levin;
//! CWP: Shuki Ronen, Jack Cohen.

use oxwave_cref::XorShift;
use oxwave_f_syscalls as rust;
use std::ffi::{c_char, c_void, CStr, CString};
use std::fs::{File, OpenOptions};
use std::os::fd::{AsRawFd, FromRawFd, IntoRawFd, OwnedFd};
use std::os::unix::fs::{MetadataExt, OpenOptionsExt};
use std::sync::atomic::{AtomicUsize, Ordering};

extern "C" {
    fn ecreat(path: *const c_char, perms: i32) -> i32;
    fn efork() -> i32;
    fn eopen(path: *const c_char, flags: i32, perms: i32) -> i32;
    fn eclose(fd: i32) -> i32;
    fn eunlink(path: *const c_char) -> i32;
    fn elseek(fd: i32, offset: i64, origin: i32) -> i64;
    fn epipe(fds: *mut i32) -> i32;
    fn eread(fd: i32, buf: *mut c_char, nbytes: usize) -> isize;
    fn ewrite(fd: i32, buf: *const c_char, nbytes: usize) -> isize;

    fn _exit(code: i32) -> !;
    fn waitpid(pid: i32, status: *mut i32, options: i32) -> i32;
    fn pipe(fds: *mut i32) -> i32;
    fn fcntl(fd: i32, cmd: i32, ...) -> i32;
    fn read(fd: i32, buf: *mut c_void, n: usize) -> isize;
    fn write(fd: i32, buf: *const c_void, n: usize) -> isize;
    fn initargs(argc: i32, argv: *mut *mut c_char);
    fn setrlimit(resource: i32, limit: *const Rlimit) -> i32;
}

// The syscall ABI and constants below are for 64-bit macOS and Linux.
#[repr(C)]
struct Rlimit {
    current: u64,
    maximum: u64,
}

#[cfg(target_os = "macos")]
const O_APPEND: i32 = 0x8;
#[cfg(target_os = "linux")]
const O_APPEND: i32 = 0x400;
#[cfg(target_os = "macos")]
const O_CREAT: i32 = 0x200;
#[cfg(target_os = "linux")]
const O_CREAT: i32 = 0x40;
#[cfg(target_os = "macos")]
const O_TRUNC: i32 = 0x400;
#[cfg(target_os = "linux")]
const O_TRUNC: i32 = 0x200;
#[cfg(target_os = "macos")]
const O_EXCL: i32 = 0x800;
#[cfg(target_os = "linux")]
const O_EXCL: i32 = 0x80;
#[cfg(target_os = "macos")]
const RLIMIT_NOFILE: i32 = 8;
#[cfg(target_os = "linux")]
const RLIMIT_NOFILE: i32 = 7;
const RLIMIT_FSIZE: i32 = 1;
const F_GETFL: i32 = 3;
const CASES: usize = 2048;
const SIZES: [usize; 13] = [
    0, 1, 2, 3, 7, 16, 127, 257, 1024, 4096, 65521, 65536, 131071,
];

static ID: AtomicUsize = AtomicUsize::new(0);

fn path() -> CString {
    CString::new(format!(
        "syscalls-test-{}-{}",
        std::process::id(),
        ID.fetch_add(1, Ordering::Relaxed)
    ))
    .unwrap()
}

fn own(fd: i32) -> OwnedFd {
    assert!(fd >= 0);
    // Each descriptor passed here is newly returned and has no other owner.
    unsafe { OwnedFd::from_raw_fd(fd) }
}

struct Paths {
    a: CString,
    b: CString,
}

impl Paths {
    fn new() -> Self {
        Self {
            a: path(),
            b: path(),
        }
    }

    fn fill(&self, bytes: &[u8]) {
        for path in [&self.a, &self.b] {
            std::fs::write(path.to_str().unwrap(), bytes).unwrap();
        }
    }

    fn open(&self) -> (OwnedFd, OwnedFd) {
        let open = |path: &CStr| {
            OpenOptions::new()
                .read(true)
                .write(true)
                .open(path.to_str().unwrap())
                .unwrap()
                .into()
        };
        (open(&self.a), open(&self.b))
    }
}

impl Drop for Paths {
    fn drop(&mut self) {
        for path in [&self.a, &self.b] {
            let _ = std::fs::remove_file(path.to_str().unwrap());
        }
    }
}

fn bytes(rng: &mut XorShift, n: usize) -> Vec<u8> {
    (0..n).map(|_| rng.next_u32() as u8).collect()
}

fn size(rng: &mut XorShift, i: usize) -> usize {
    if i < SIZES.len() {
        SIZES[i]
    } else {
        rng.next_u32() as usize % 4097
    }
}

fn raw_write(fd: i32, data: &[u8]) {
    assert_eq!(
        unsafe { write(fd, data.as_ptr().cast(), data.len()) },
        data.len() as isize
    );
}

fn raw_read(fd: i32, data: &mut [u8]) -> isize {
    unsafe { read(fd, data.as_mut_ptr().cast(), data.len()) }
}

#[test]
fn ecreat_matches() {
    let mut rng = XorShift::new(1);
    let modes = [
        -1,
        0,
        1,
        2,
        0o200,
        0o400,
        0o600,
        0o777,
        0o7777,
        i32::MIN,
        i32::MAX,
    ];
    for i in 0..CASES {
        let paths = Paths::new();
        // Alternate creating a new file and truncating an existing one. The
        // permissions argument only changes the mode when the file is new.
        if i % 2 == 1 {
            paths.fill(&bytes(&mut rng, i % 1031));
        }
        let mode = if i < modes.len() {
            modes[i]
        } else {
            rng.next_u32() as i32
        };
        let a = own(rust::ecreat(&paths.a, mode));
        let b = own(unsafe { ecreat(paths.b.as_ptr(), mode) });
        let am = std::fs::metadata(paths.a.to_str().unwrap()).unwrap();
        let bm = std::fs::metadata(paths.b.to_str().unwrap()).unwrap();
        assert_eq!(am.mode() & 0o7777, bm.mode() & 0o7777);
        assert_eq!((am.len(), bm.len()), (0, 0));
        // Even mode 000 files are writable through the returned descriptor.
        let data = bytes(&mut rng, i % 257);
        raw_write(a.as_raw_fd(), &data);
        raw_write(b.as_raw_fd(), &data);
        assert_eq!(unsafe { fcntl(a.as_raw_fd(), F_GETFL) }, unsafe {
            fcntl(b.as_raw_fd(), F_GETFL)
        });
        assert_eq!(
            std::fs::metadata(paths.a.to_str().unwrap()).unwrap().len(),
            data.len() as u64
        );
        assert_eq!(
            std::fs::metadata(paths.b.to_str().unwrap()).unwrap().len(),
            data.len() as u64
        );
    }
}

#[test]
fn eopen_matches() {
    let mut rng = XorShift::new(2);
    let flags = [
        0,
        1,
        2,
        1 | O_TRUNC,
        1 | O_APPEND,
        2 | O_CREAT,
        2 | O_CREAT | O_EXCL,
        2 | O_CREAT | O_APPEND,
    ];
    for i in 0..CASES {
        let paths = Paths::new();
        let flag = flags[if i < flags.len() {
            i
        } else {
            rng.next_u32() as usize % flags.len()
        }];
        if flag & O_CREAT == 0 {
            paths.fill(&bytes(&mut rng, i % 521));
        }
        let mode = (rng.next_u32() as i32 & 0o777) | 0o600;
        let a = own(rust::eopen(&paths.a, flag, mode));
        let b = own(unsafe { eopen(paths.b.as_ptr(), flag, mode) });
        assert_eq!(unsafe { fcntl(a.as_raw_fd(), F_GETFL) }, unsafe {
            fcntl(b.as_raw_fd(), F_GETFL)
        });
        assert_eq!(
            std::fs::metadata(paths.a.to_str().unwrap()).unwrap().mode() & 0o777,
            std::fs::metadata(paths.b.to_str().unwrap()).unwrap().mode() & 0o777
        );
        if flag & 3 != 0 {
            let data = bytes(&mut rng, i % 257);
            raw_write(a.as_raw_fd(), &data);
            raw_write(b.as_raw_fd(), &data);
        } else {
            let mut u = [0xa5; 1024];
            let mut v = u;
            assert_eq!(
                raw_read(a.as_raw_fd(), &mut u),
                raw_read(b.as_raw_fd(), &mut v)
            );
            assert_eq!(u, v);
        }
        assert_eq!(
            std::fs::read(paths.a.to_str().unwrap()).unwrap(),
            std::fs::read(paths.b.to_str().unwrap()).unwrap()
        );
    }
}

#[test]
fn eclose_matches() {
    let paths = Paths::new();
    paths.fill(&[]);
    let mut rng = XorShift::new(3);
    for _ in 0..CASES {
        if rng.next_u32() & 1 == 0 {
            let (a, b) = paths.open();
            let a = a.into_raw_fd();
            let b = b.into_raw_fd();
            assert_eq!(rust::eclose(a), unsafe { eclose(b) });
        } else {
            let mut a = [0; 2];
            let mut b = [0; 2];
            assert_eq!(unsafe { pipe(a.as_mut_ptr()) }, 0);
            assert_eq!(unsafe { pipe(b.as_mut_ptr()) }, 0);
            let first = rng.next_u32() as usize & 1;
            for i in [first, 1 - first] {
                assert_eq!(rust::eclose(a[i]), unsafe { eclose(b[i]) });
            }
        }
    }
}

#[test]
fn eunlink_matches() {
    let mut rng = XorShift::new(4);
    for i in 0..CASES {
        let paths = Paths::new();
        let n = size(&mut rng, i);
        let data = bytes(&mut rng, n);
        paths.fill(&data);
        let (a, b) = paths.open();
        assert_eq!(rust::eunlink(&paths.a), unsafe {
            eunlink(paths.b.as_ptr())
        });
        for path in [&paths.a, &paths.b] {
            assert_eq!(
                std::fs::metadata(path.to_str().unwrap())
                    .unwrap_err()
                    .kind(),
                std::io::ErrorKind::NotFound
            );
        }
        // Unlinking removes the names, but existing descriptors stay usable.
        let mut u = vec![0xa5; data.len() + 1];
        let mut v = u.clone();
        assert_eq!(raw_read(a.as_raw_fd(), &mut u), data.len() as isize);
        assert_eq!(raw_read(b.as_raw_fd(), &mut v), data.len() as isize);
        assert_eq!(u, v);
        assert_eq!(&u[..data.len()], data);
        assert_eq!(u[data.len()], 0xa5);
    }
}

#[test]
fn elseek_matches() {
    let mut rng = XorShift::new(5);
    let paths = Paths::new();
    let data = bytes(&mut rng, 4096);
    paths.fill(&data);
    let (a, b) = paths.open();
    let a = a.as_raw_fd();
    let b = b.as_raw_fd();
    let offsets = [0, 1, 2, 3, 127, 257, 4096, 65521, 1 << 31, 1 << 32, 1 << 40];
    for i in 0..5000 {
        let offset = if i < offsets.len() {
            offsets[i]
        } else {
            rng.next_u32() as i64
        };
        assert_eq!(rust::elseek(a, offset, 0), unsafe { elseek(b, offset, 0) });
        assert_eq!(rust::elseek(a, -(offset / 2), 1), unsafe {
            elseek(b, -(offset / 2), 1)
        });
        // Include valid negative offsets from a nonempty file's end.
        let offset = (rng.next_u32() % 8193) as i64 - 4096;
        let position = rust::elseek(a, offset, 2);
        assert_eq!(position, unsafe { elseek(b, offset, 2) });
        assert_eq!(position, 4096 + offset);
        let mut u = [0xa5; 32];
        let mut v = u;
        assert_eq!(raw_read(a, &mut u), raw_read(b, &mut v));
        assert_eq!(u, v);
    }
}

#[test]
fn epipe_matches() {
    let mut rng = XorShift::new(6);
    for i in 0..CASES {
        let mut a = [0; 2];
        let mut b = [0; 2];
        assert_eq!(rust::epipe(&mut a), unsafe { epipe(b.as_mut_ptr()) });
        let (ar, aw, br, bw) = (own(a[0]), own(a[1]), own(b[0]), own(b[1]));
        assert_ne!(ar.as_raw_fd(), aw.as_raw_fd());
        assert_ne!(br.as_raw_fd(), bw.as_raw_fd());
        // Stay below pipe capacity so these single-threaded writes cannot block.
        let n = if i < 8 {
            SIZES[i]
        } else {
            rng.next_u32() as usize % 1025
        };
        let data = bytes(&mut rng, n);
        raw_write(aw.as_raw_fd(), &data);
        raw_write(bw.as_raw_fd(), &data);
        let mut u = vec![0; n];
        let mut v = u.clone();
        assert_eq!(raw_read(ar.as_raw_fd(), &mut u), n as isize);
        assert_eq!(raw_read(br.as_raw_fd(), &mut v), n as isize);
        assert_eq!(u, v);
        assert_eq!(u, data);
        drop(aw);
        drop(bw);
        let mut eof = [0];
        assert_eq!(raw_read(ar.as_raw_fd(), &mut eof), 0);
        assert_eq!(raw_read(br.as_raw_fd(), &mut eof), 0);
    }
}

#[test]
fn eread_matches() {
    let paths = Paths::new();
    let mut rng = XorShift::new(7);
    let data = bytes(&mut rng, 65536);
    paths.fill(&data);
    let (a, b) = paths.open();
    let a = a.as_raw_fd();
    let b = b.as_raw_fd();
    for i in 0..3000 {
        let start = if i < SIZES.len() {
            0
        } else {
            rng.next_u32() as usize % (data.len() + 1024)
        };
        let n = size(&mut rng, i);
        rust::elseek(a, start as i64, 0);
        unsafe {
            elseek(b, start as i64, 0);
        }
        let mut u = vec![0xa5; n];
        let mut v = u.clone();
        let count = rust::eread(a, &mut u);
        assert_eq!(count, unsafe { eread(b, v.as_mut_ptr().cast(), n) });
        assert_eq!(count as usize, n.min(data.len().saturating_sub(start)));
        assert_eq!(u, v);
        if count > 0 {
            assert_eq!(&u[..count as usize], &data[start..start + count as usize]);
        }
        assert!(u[count as usize..].iter().all(|&v| v == 0xa5));
    }
}

#[test]
fn ewrite_matches() {
    let paths = Paths::new();
    paths.fill(&[]);
    let (a, b) = paths.open();
    let a = a.as_raw_fd();
    let b = b.as_raw_fd();
    let mut rng = XorShift::new(8);
    for i in 0..3000 {
        let n = size(&mut rng, i);
        let data = if i % 17 == 0 {
            vec![0; n]
        } else {
            bytes(&mut rng, n)
        };
        let start = (rng.next_u32() % 65537) as i64;
        rust::elseek(a, start, 0);
        unsafe {
            elseek(b, start, 0);
        }
        assert_eq!(rust::ewrite(a, &data), unsafe {
            ewrite(b, data.as_ptr().cast(), n)
        });
        assert_eq!(rust::elseek(a, 0, 1), start + n as i64);
        assert_eq!(unsafe { elseek(b, 0, 1) }, start + n as i64);
        rust::elseek(a, start, 0);
        unsafe {
            elseek(b, start, 0);
        }
        let mut u = vec![0xa5; n];
        let mut v = u.clone();
        assert_eq!(raw_read(a, &mut u), n as isize);
        assert_eq!(raw_read(b, &mut v), n as isize);
        assert_eq!(u, v);
        assert_eq!(u, data);
    }
    // Also compare unwritten holes and bytes left from overlapping writes.
    assert_eq!(
        std::fs::read(paths.a.to_str().unwrap()).unwrap(),
        std::fs::read(paths.b.to_str().unwrap()).unwrap()
    );
}

#[test]
fn efork_matches() {
    let mut rng = XorShift::new(9);
    for _ in 0..CASES {
        let code = (rng.next_u32() & 255) as i32;
        let data = rng.next_u32().to_ne_bytes();
        for c in [false, true] {
            let mut fds = [0; 2];
            assert_eq!(unsafe { pipe(fds.as_mut_ptr()) }, 0);
            let (reader, writer) = (own(fds[0]), own(fds[1]));
            let pid = if c { unsafe { efork() } } else { rust::efork() };
            if pid == 0 {
                // Only async-signal-safe libc calls after forking the threaded
                // test harness: no allocation, assertions, locks or destructors.
                unsafe {
                    write(fds[1], data.as_ptr().cast(), data.len());
                    _exit(code);
                }
            }
            assert!(pid > 0);
            drop(writer);
            let mut status = 0;
            assert_eq!(unsafe { waitpid(pid, &mut status, 0) }, pid);
            assert_eq!(status, code << 8);
            let mut got = [0; 4];
            assert_eq!(raw_read(reader.as_raw_fd(), &mut got), 4);
            assert_eq!(got, data);
            assert_eq!(raw_read(reader.as_raw_fd(), &mut got), 0);
        }
    }
}

#[test]
fn error_paths_match() {
    if let Ok(case) = std::env::var("OX_SYSCALL_ERROR_CASE") {
        let c = std::env::var_os("OX_SYSCALL_C").is_some();
        let arg = CString::new("syscalls-test").unwrap();
        let mut args = [arg.as_ptr() as *mut c_char, std::ptr::null_mut()];
        if c {
            unsafe {
                initargs(1, args.as_mut_ptr());
            }
        } else {
            oxwave_core::par::initargs(&["syscalls-test".into()]);
        }
        let p = CString::new("/dev/null/oxwave-missing").unwrap();
        let mut byte = [0u8; 1];
        unsafe {
            match (case.as_str(), c) {
                ("creat", true) => {
                    ecreat(p.as_ptr(), 0o600);
                }
                ("creat", false) => {
                    rust::ecreat(&p, 0o600);
                }
                ("open", true) => {
                    eopen(p.as_ptr(), 0, 0);
                }
                ("open", false) => {
                    rust::eopen(&p, 0, 0);
                }
                ("unlink", true) => {
                    eunlink(p.as_ptr());
                }
                ("unlink", false) => {
                    rust::eunlink(&p);
                }
                ("close", true) => {
                    eclose(-1);
                }
                ("close", false) => {
                    rust::eclose(-1);
                }
                ("close-twice", c) => {
                    // Isolate descriptor reuse checks from the other tests,
                    // which continually allocate descriptors in parallel.
                    let fd = File::open("/dev/null").unwrap().into_raw_fd();
                    if c {
                        assert_eq!(eclose(fd), 0);
                        eclose(fd);
                    } else {
                        assert_eq!(rust::eclose(fd), 0);
                        rust::eclose(fd);
                    }
                }
                ("seek", true) => {
                    elseek(-1, -1, -1);
                }
                ("seek", false) => {
                    rust::elseek(-1, -1, -1);
                }
                ("read", true) => {
                    eread(-1, byte.as_mut_ptr().cast(), 1);
                }
                ("read", false) => {
                    rust::eread(-1, &mut byte);
                }
                ("read-zero", true) => {
                    eread(i32::MIN, byte.as_mut_ptr().cast(), 0);
                }
                ("read-zero", false) => {
                    rust::eread(i32::MIN, &mut []);
                }
                ("write", true) => {
                    ewrite(-1, byte.as_ptr().cast(), 1);
                }
                ("write", false) => {
                    rust::ewrite(-1, &byte);
                }
                ("write-zero", true) => {
                    ewrite(i32::MAX, byte.as_ptr().cast(), 0);
                }
                ("write-zero", false) => {
                    rust::ewrite(i32::MAX, &[]);
                }
                ("seek-pipe", c) => {
                    let mut fds = [0; 2];
                    assert_eq!(pipe(fds.as_mut_ptr()), 0);
                    if c {
                        elseek(fds[0], 0, 0);
                    } else {
                        rust::elseek(fds[0], 0, 0);
                    }
                }
                ("seek-negative", c) | ("seek-origin", c) => {
                    let file = File::open("Cargo.toml").unwrap();
                    let (offset, origin) = if case == "seek-negative" {
                        (-1, 0)
                    } else {
                        (0, -1)
                    };
                    if c {
                        elseek(file.as_raw_fd(), offset, origin);
                    } else {
                        rust::elseek(file.as_raw_fd(), offset, origin);
                    }
                }
                ("read-directory", c) => {
                    let file = File::open(".").unwrap();
                    if c {
                        eread(file.as_raw_fd(), byte.as_mut_ptr().cast(), 1);
                    } else {
                        rust::eread(file.as_raw_fd(), &mut byte);
                    }
                }
                ("write-readonly", c) => {
                    let file = File::open("/dev/null").unwrap();
                    if c {
                        ewrite(file.as_raw_fd(), byte.as_ptr().cast(), 1);
                    } else {
                        rust::ewrite(file.as_raw_fd(), &byte);
                    }
                }
                ("pipe-exhausted", c) => {
                    // Only this subprocess's descriptor limit is changed.
                    let limit = Rlimit {
                        current: 16,
                        maximum: 16,
                    };
                    assert_eq!(setrlimit(RLIMIT_NOFILE, &limit), 0);
                    let mut files = Vec::with_capacity(16);
                    while let Ok(file) = File::open("/dev/null") {
                        files.push(file);
                    }
                    assert!(files.len() < 16);
                    let mut fds = [0; 2];
                    if c {
                        epipe(fds.as_mut_ptr());
                    } else {
                        rust::epipe(&mut fds);
                    }
                }
                ("write-partial", c) => {
                    let path = std::env::var("OX_SYSCALL_FILE").unwrap();
                    let file = OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .mode(0o600)
                        .open(path)
                        .unwrap();
                    let data = [0x5a; 2048];
                    let limit = Rlimit {
                        current: 1024,
                        maximum: 1024,
                    };
                    assert_eq!(setrlimit(RLIMIT_FSIZE, &limit), 0);
                    // The first write crosses the limit and transfers 1024
                    // bytes. Both wrappers must exit rather than return/retry.
                    if c {
                        ewrite(file.as_raw_fd(), data.as_ptr().cast(), data.len());
                    } else {
                        rust::ewrite(file.as_raw_fd(), &data);
                    }
                }
                _ => panic!("unknown case"),
            }
        }
        panic!("error wrapper returned");
    }
    for case in [
        "creat",
        "open",
        "unlink",
        "close",
        "close-twice",
        "seek",
        "read",
        "read-zero",
        "write",
        "write-zero",
        "seek-pipe",
        "seek-negative",
        "seek-origin",
        "read-directory",
        "write-readonly",
        "pipe-exhausted",
        "write-partial",
    ] {
        let paths = Paths::new();
        let run = |c| {
            let mut command = std::process::Command::new(std::env::current_exe().unwrap());
            command
                .args(["--exact", "error_paths_match", "--nocapture"])
                .env("OX_SYSCALL_ERROR_CASE", case)
                .env_remove("OX_SYSCALL_C")
                .env(
                    "OX_SYSCALL_FILE",
                    if c {
                        paths.b.to_str().unwrap()
                    } else {
                        paths.a.to_str().unwrap()
                    },
                );
            if c {
                command.env("OX_SYSCALL_C", "1");
            }
            command.output().unwrap()
        };
        let r = run(false);
        let c = run(true);
        assert_eq!(
            r.status.code(),
            Some(1),
            "{case}: {}",
            String::from_utf8_lossy(&r.stderr)
        );
        assert_eq!(
            r.status.code(),
            c.status.code(),
            "{case}: {}",
            String::from_utf8_lossy(&c.stderr)
        );
        assert_eq!(r.stderr, c.stderr, "{case}");
        assert!(!r.stderr.is_empty());
        if case == "write-partial" {
            let a = std::fs::read(paths.a.to_str().unwrap()).unwrap();
            let b = std::fs::read(paths.b.to_str().unwrap()).unwrap();
            assert_eq!(a, vec![0x5a; 1024]);
            assert_eq!(a, b);
        }
    }
}
