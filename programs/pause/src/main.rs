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

//! PAUSE - prompt and wait for user signal to continue.
//! Port of cwp/main/pause.c. Credits: CWP: Shuki.

use std::os::unix::ffi::OsStrExt;
use oxwave_core::cio::{exit, fgetc, fputs_bytes, restore_sigpipe, stdin, stdout};

fn main() {
    restore_sigpipe();
    let argv: Vec<_> = std::env::args_os().collect();
    let mut do_what = Vec::new();
    if argv.len() == 1 {
        do_what.extend_from_slice(b"continue");
    } else {
        do_what.extend_from_slice(argv[1].as_os_str().as_bytes());
        for arg in &argv[2..] {
            do_what.push(b' ');
            do_what.extend_from_slice(arg.as_os_str().as_bytes());
        }
    }
    fputs_bytes(b"\npress return key to ", &stdout());
    fputs_bytes(&do_what, &stdout());
    fputs_bytes(b"\n", &stdout());
    fgetc(&stdin());
    exit(0);
}
