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

//! `err`: SU's `err`, `warn` and `syserr` (par/lib/errpkge.c).
//!
//! All three take a **C printf format** (see `cfmt`) and print
//! `"\n<argv0>: <message>\n"` on stderr after flushing stdout, where
//! `<argv0>` is `xargv[0]` exactly as invoked (set by `initargs`).
//! `err!` and `syserr!` then exit with status 1 (`EXIT_FAILURE`), flushing
//! all streams; `syserr!` appends `" (<strerror(errno)>)"`.
//!
//! | C                      | Rust                       |
//! |------------------------|----------------------------|
//! | `err("x=%d", x)`       | `err!("x=%d", x)` (type `!`) |
//! | `warn("x=%g", x)`      | `warn!("x=%g", x)`         |
//! | `syserr("can't %s", f)`| `syserr!("can't %s", f)`   |
//!
//! Library code that uses `__FILE__` in C keeps the literal C file name,
//! e.g. `err!("%s: bad first header", "fgettr.c")`.

use crate::cio::{fflush, fputs, stderr, stdout};

fn prologue(who: &str) {
    if fflush(&stdout()) == crate::cio::EOF {
        fputs(&format!("\n{}: fflush failed on stdout", who), &stderr());
    }
    let argv0 = crate::par::xargv().first().cloned().unwrap_or_default();
    fputs(&format!("\n{}: ", argv0), &stderr());
}

#[doc(hidden)]
pub fn err_impl(msg: &str) -> ! {
    prologue("err");
    fputs(msg, &stderr());
    fputs("\n", &stderr());
    crate::cio::exit(crate::cio::EXIT_FAILURE)
}

#[doc(hidden)]
pub fn warn_impl(msg: &str) {
    prologue("warn");
    fputs(msg, &stderr());
    fputs("\n", &stderr());
}

/// C `strerror(errnum)`.
pub fn strerror(errnum: i32) -> String {
    extern "C" {
        #[link_name = "strerror"]
        fn c_strerror(e: i32) -> *const std::ffi::c_char;
    }
    // SAFETY: strerror returns a valid NUL-terminated static string.
    unsafe { std::ffi::CStr::from_ptr(c_strerror(errnum)) }
        .to_string_lossy()
        .into_owned()
}

#[doc(hidden)]
pub fn syserr_impl(msg: &str) -> ! {
    let e = std::io::Error::last_os_error().raw_os_error().unwrap_or(0);
    prologue("syserr");
    fputs(msg, &stderr());
    fputs(&format!(" ({})\n", strerror(e)), &stderr());
    crate::cio::exit(crate::cio::EXIT_FAILURE)
}

/// `err!(c_fmt, args..)`: print "\nprog: msg\n" on stderr and exit(1).
#[macro_export]
macro_rules! err {
    ($fmt:expr $(, $arg:expr)* $(,)?) => {
        $crate::err::err_impl(&$crate::sprintf!($fmt $(, $arg)*))
    };
}

/// `warn!(c_fmt, args..)`: print "\nprog: msg\n" on stderr and continue.
#[macro_export]
macro_rules! warn {
    ($fmt:expr $(, $arg:expr)* $(,)?) => {
        $crate::err::warn_impl(&$crate::sprintf!($fmt $(, $arg)*))
    };
}

/// `syserr!(c_fmt, args..)`: like err! plus " (strerror(errno))".
#[macro_export]
macro_rules! syserr {
    ($fmt:expr $(, $arg:expr)* $(,)?) => {
        $crate::err::syserr_impl(&$crate::sprintf!($fmt $(, $arg)*))
    };
}
