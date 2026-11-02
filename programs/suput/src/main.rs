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

//! SUPUT - Connect SU program to file descriptor for output stream.
//! Port of su/main/supromax/suput.c (Revision 1.7).
//! Credits: John Anderson (visiting scholar from Mobil), July 1994.

use oxwave_core::prelude::*;
use std::os::fd::FromRawFd;

unsafe extern "C" {
    fn dup(fd: std::os::raw::c_int) -> std::os::raw::c_int;
    fn fcntl(fd: std::os::raw::c_int, cmd: std::os::raw::c_int, ...) -> std::os::raw::c_int;
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUPUT - Connect SU program to file descriptor for output stream.\t",
    "\t\t\t\t\t\t\t\t\t",
    "       su_module | suput fp=$1\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " This program is for interfacing \" outside processing systems \"\t",
    " with SU. Typically, the outside system would execute the SU command file.",
    " The outside system provides the file descriptor it would like to read\t",
    " from to the command file to be an argument for suput.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example: su_module | suput fp=$1\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "       fd=-1       file_descriptor_for_output_stream_from_su\t\t",
    "       verbose=0   minimal listing\t\t\t\t\t",
    "                   =1  asks for message with each trace processed.\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut fd = -1;
    let mut nread = 0;
    let mut verbose = 0;

    initargs(&argv);
    requestdoc(1, SDOC);
    if getparint("fd", &mut fd) == 0 { fd = -1; }
    fprintf!(&stderr(), "File descriptor passed to suput = %d\n", fd);
    if getparint("verbose", &mut verbose) == 0 { verbose = 0; }
    checkpars();

    // fdopen(fd, "w") checks the descriptor's access mode without writing.
    // A zero-byte write can fail on a valid writable descriptor (for example,
    // an unconnected datagram socket), which fdopen must still accept.
    // F_GETFL and O_ACCMODE are 3, and O_RDONLY is 0, on macOS and Linux.
    // SAFETY: F_GETFL takes only an integer descriptor and no pointer argument.
    let flags = unsafe { fcntl(fd, 3) };
    if flags < 0 || flags & 3 == 0 { err!("Bad file descriptor \n"); }
    // dup accepts invalid descriptors and reports failure without creating a
    // Rust descriptor whose safety contract would require fd to be valid.
    // SAFETY: dup takes any integer descriptor and has no pointer arguments.
    let output_fd = unsafe { dup(fd) };
    if output_fd < 0 { err!("Bad file descriptor \n"); }
    // SAFETY: a successful dup returns a valid, newly owned descriptor.
    let file = unsafe { std::fs::File::from_raw_fd(output_fd) };
    let fp = CFile::from_file(file);

    if gettr(&mut tr) == 0 { err!("Can't get first trace \n"); }
    loop {
        if verbose > 0 {
            warn!("suput: read input traces %d", nread);
            nread += 1;
        }
        fputtr(&fp, &mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    oxwave_core::cio::fclose(&fp);
    exit(CWP_Exit());
}
