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
//! Port of par/lib/syscalls.c.
//! Authors: SEP: Rick Otolini, Ron, Jon Claerbout, Stew Levin;
//! CWP: Shuki Ronen, Jack Cohen.
//!
//! Each wrapper makes exactly one system call and exits through `syserr` on
//! failure. Reads may return fewer bytes than requested, including zero at EOF.
//! Writes must transfer the entire slice; a short write is fatal, as in C.

use oxwave_core::syserr;
use std::ffi::{c_char, c_int, c_void, CStr};

extern "C" {
    fn creat(path: *const c_char, mode: c_int) -> c_int;
    fn fork() -> c_int;
    fn open(path: *const c_char, flags: c_int, ...) -> c_int;
    fn close(fd: c_int) -> c_int;
    fn unlink(path: *const c_char) -> c_int;
    fn lseek(fd: c_int, offset: i64, whence: c_int) -> i64;
    fn pipe(fds: *mut c_int) -> c_int;
    fn read(fd: c_int, buf: *mut c_void, count: usize) -> isize;
    fn write(fd: c_int, buf: *const c_void, count: usize) -> isize;
    #[cfg(target_os = "macos")]
    fn __error() -> *mut c_int;
    #[cfg(not(target_os = "macos"))]
    fn __errno_location() -> *mut c_int;
}

/// C prototype: `int ecreat(char *path, int perms);`
pub fn ecreat(path: &CStr, perms: i32) -> i32 {
    let fd = unsafe { creat(path.as_ptr(), perms) };
    if fd == -1 {
        syserr!("%s: ecreat: creat failed", "syscalls.c");
    }
    fd
}

/// C prototype: `int efork(void);`
///
/// Returns the child PID in the parent and zero in the child.
///
/// After a fork in a multithreaded process, the child must follow POSIX
/// restrictions until it calls an exec function or exits.
pub fn efork() -> i32 {
    let fd = unsafe { fork() };
    if fd == -1 {
        syserr!("%s: efork: fork failed", "syscalls.c");
    }
    fd
}

/// C prototype: `int eopen(char *path, int flags, int perms);`
pub fn eopen(path: &CStr, flags: i32, perms: i32) -> i32 {
    let fd = unsafe { open(path.as_ptr(), flags, perms) };
    if fd == -1 {
        syserr!("%s: eopen: open failed", "syscalls.c");
    }
    fd
}

/// C prototype: `int eclose(int fd);`
pub fn eclose(fd: i32) -> i32 {
    let result = unsafe { close(fd) };
    if result == -1 {
        syserr!("%s: eclose: close failed", "syscalls.c");
    }
    result
}

/// C prototype: `int eunlink(char *path);`
pub fn eunlink(path: &CStr) -> i32 {
    let result = unsafe { unlink(path.as_ptr()) };
    if result == -1 {
        syserr!("%s: eunlink: unlink failed", "syscalls.c");
    }
    result
}

/// C prototype: `off_t elseek(int fd, off_t offset, int origin);`
pub fn elseek(fd: i32, offset: i64, origin: i32) -> i64 {
    let position = unsafe { lseek(fd, offset, origin) };
    if position == -1 {
        syserr!("%s: elseek: lseek failed", "syscalls.c");
    }
    position
}

/// C prototype: `int epipe(int fd[2]);`
pub fn epipe(fds: &mut [i32; 2]) -> i32 {
    let result = unsafe { pipe(fds.as_mut_ptr()) };
    if result == -1 {
        syserr!("%s: epipe: pipe failed", "syscalls.c");
    }
    result
}

/// C prototype: `ssize_t eread(int fd, char *buf, size_t nbytes);`
///
/// The slice length supplies `nbytes`; bytes beyond the returned count stay intact.
pub fn eread(fd: i32, buf: &mut [u8]) -> isize {
    let n = unsafe { read(fd, buf.as_mut_ptr().cast(), buf.len()) };
    if n == -1 {
        syserr!("%s: eread: read failed", "syscalls.c");
    }
    n
}

/// C prototype: `ssize_t ewrite(int fd, char *buf, size_t nbytes);`
///
/// The slice length supplies `nbytes`. Like C, this sets errno to `EFBIG`
/// on a short or failed write before reporting the error. It does not retry.
pub fn ewrite(fd: i32, buf: &[u8]) -> isize {
    let n = unsafe { write(fd, buf.as_ptr().cast(), buf.len()) };
    if n as usize != buf.len() {
        // The C compares nbytes, not nwrite, to -1 before setting errno.
        if buf.len() != usize::MAX {
            #[cfg(target_os = "macos")]
            unsafe {
                *__error() = 27; // EFBIG on macOS.
            }
            #[cfg(not(target_os = "macos"))]
            unsafe {
                *__errno_location() = 27; // EFBIG on Linux.
            }
        }
        syserr!("%s: ewrite: write failed", "syscalls.c");
    }
    n
}
