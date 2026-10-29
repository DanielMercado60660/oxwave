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

//! FCAT - fast cat with one read per file. Port of cwp/main/fcat.c.
//! Credits: Shuki.

use oxwave_core::prelude::*;
use std::ffi::c_void;
use std::fs::File;
use std::io::{self, Write};
use std::os::fd::AsRawFd;
use std::os::unix::ffi::OsStrExt;

const MAXSIZE: i32 = 8 * 1024 * 1024;

extern "C" {
    fn read(fd: i32, buf: *mut c_void, count: usize) -> isize;
    fn write(fd: i32, buf: *const c_void, count: usize) -> isize;
}

fn message(prefix: &[u8], name: &std::ffi::OsStr, suffix: &[u8]) {
    let mut stderr = io::stderr().lock();
    let _ = stderr.write_all(prefix);
    let _ = stderr.write_all(name.as_bytes());
    let _ = stderr.write_all(suffix);
}

fn main() {
    let argv = args();
    initargs(&argv);
    let av: Vec<_> = std::env::args_os().collect();
    if av.len() < 2 {
        let _ = io::stderr().write_all(b"usage: fcat file1 file2 ... > file3\n");
        exit(1);
    }

    let mut x = vec![0u8; MAXSIZE as usize];
    // A failed fstat leaves the C stat buffer unchanged. Before the first
    // success its contents are uninitialized and the result is undefined.
    let mut size = i32::MAX;
    for name in &av[1..] {
        let file = File::open(name).ok();
        if file.is_none() {
            message(b"can't open ", name, b"\n");
        }
        if let Some(ref f) = file {
            if let Ok(metadata) = f.metadata() {
                size = metadata.len() as i32;
            }
        }
        let n = size;
        if n > MAXSIZE {
            message(b"", name, b" is too big\n");
            exit(1);
        }
        let fd = file.as_ref().map_or(-1, |f| f.as_raw_fd());
        // A single syscall is intentional, including for negative sizes.
        if unsafe { read(fd, x.as_mut_ptr().cast(), n as usize) } != n as isize {
            message(b"read error from ", name, b"\n");
            exit(1);
        }
        drop(file);
        if unsafe { write(1, x.as_ptr().cast(), n as usize) } != n as isize {
            message(b"write error from ", name, b"\n");
            exit(1);
        }
    }
    exit(CWP_Exit());
}
