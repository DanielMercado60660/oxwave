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

//! PERCENT - find 100 * ratio of arguments. Port of su/shell/percent.c
//! (Revision 1.1). Author: john.
//!
//! Usage: percent x y
//! This triviality is called by the rmaxdiff shell. Surely, there
//! is a better way, but ...
//!
//! The C source has no sdoc array (confirmed with tools/sdoc2rs.py),
//! parameter processing, or stdin processing.

use oxwave_core::prelude::*;

unsafe extern "C" {
    fn raise(sig: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
}

fn main() {
    let argv = args();
    // C starts with SIGPIPE at its default disposition; Rust ignores it.
    oxwave_core::cio::restore_sigpipe();
    // The original passes a null pointer to atof when an argument is absent.
    if argv.len() < 3 {
        // SAFETY: SIGSEGV is 11 and SIG_DFL is 0 on the target platform.
        // Reproduce the original's observed termination without a null read.
        unsafe {
            signal(11, 0);
            raise(11);
        }
        exit(1);
    }
    let x: f64 = atof(&argv[1]);
    let y: f64 = atof(&argv[2]);
    printf!("%g\n", if y == 0.0 { 0.0 } else { 100.0 * x / y });
    exit(CWP_Exit());
}
