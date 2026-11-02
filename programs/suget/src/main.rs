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


//! SUGET - Connect SU program to file descriptor for input stream.
//! Port of su/main/supromax/suget.c (Revision 1.7).
//! Credits (from the C source): John Anderson (visiting scholar from Mobil), July 1994.

use oxwave_core::prelude::*;
use std::fs::File;
use std::os::fd::FromRawFd;

unsafe extern "C" {
    fn fcntl(fd: i32, cmd: i32, ...) -> i32;
}

// Descriptor flags on macOS and Linux.
const F_GETFL: i32 = 3;
const O_ACCMODE: i32 = 3;
const O_WRONLY: i32 = 1;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUGET  - Connect SU program to file descriptor for input stream.\t",
    "\t\t\t\t\t\t\t\t\t",
    "    suget fd=$1 | next_su_module\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " This program is for interfacing \" outside processing systems \"\t",
    " with SU. Typically, an outside system would execute the su command file",
    " and a file descriptor would be passed by an outside system to\t\t",
    " the su command file so that output data from the outside system\t",
    " could be piped into the su programs executing inside the command file.",
    "\t\t\t\t\t\t\t\t\t",
    " Example:    suget fd=$1 | next_su_module\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "      fd=-1        file_descriptor_for_input_stream\t\t\t",
    "      verbose=0    minimal listing\t\t\t\t\t",
    "                   =1  asks for message with each trace processed.\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut fd = -1;
    let mut j = 0;
    let mut verbose = 0;

    initargs(&argv);
    requestdoc(1, SDOC);

    if getparint("fd", &mut fd) == 0 { fd = -1; }
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    checkpars();
    warn!("File descriptor passed to suget = %d", fd);

    // Match fdopen(fd, "r"): reject invalid/write-only descriptors, then
    // wrap the existing descriptor, preserving sockets, flags and position.
    // Darwin's FILE stores its descriptor in a signed short, so fdopen also
    // rejects valid descriptors above SHRT_MAX; File::from_raw_fd does not.
    // SAFETY: F_GETFL only inspects the descriptor and takes no third argument.
    let flags = unsafe { fcntl(fd, F_GETFL) };
    if flags == -1 || flags & O_ACCMODE == O_WRONLY
        || (cfg!(target_os = "macos") && fd > i16::MAX as i32)
    {
        err!("Bad file descriptor");
    }
    // SAFETY: F_GETFL validated fd. This stream owns it until process exit;
    // core's standard stream handles never close descriptors 0, 1 or 2.
    let file = unsafe { File::from_raw_fd(fd) };
    let fp = CFile::from_file(file);
    warn!("About to read first trace");
    if fgettr(&fp, &mut tr) == 0 { err!("Can't get first trace"); }

    loop {
        if verbose > 0 {
            warn!("read trace %d", j);
            j += 1;
        }
        puttr(&mut tr);
        if fgettr(&fp, &mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
