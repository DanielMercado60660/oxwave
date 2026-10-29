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

//! `doc`: selfdoc (par/lib/docpkge.c).
//!
//! | C                   | Rust                          |
//! |---------------------|-------------------------------|
//! | `char *sdoc[] = {..., NULL};` | `const SDOC: &[&str] = &[...];` (no NULL; strings verbatim, incl. tabs/trailing blanks) |
//! | `requestdoc(flag)`  | `requestdoc(flag, SDOC)`      |
//! | `pagedoc()`         | `pagedoc(SDOC)` (never returns) |
//!
//! `requestdoc(1)`: selfdoc iff xargc==1 and stdin is a tty.
//! `requestdoc(0)`: iff xargc==1 and stdin and stdout are ttys.
//! `requestdoc(n>1)`: iff xargc <= n.
//! `pagedoc` flushes stdout, pipes each line + "\n" into `sh -c "$PAGER 1>&2"`
//! (`more 1>&2` when PAGER is unset), waits for it, then exits with status 1.
//! So selfdoc text ends up on **stderr** and the exit code is 1.

use crate::cio::{isatty, STDIN, STDOUT};

/// requestdoc: print selfdoc as directed by flag (see module docs).
pub fn requestdoc(flag: i32, sdoc: &[&str]) {
    let xargc = crate::par::xargc() as i32;
    match flag {
        1 => {
            if xargc == 1 && isatty(STDIN) {
                pagedoc(sdoc);
            }
        }
        0 => {
            if xargc == 1 && isatty(STDIN) && isatty(STDOUT) {
                pagedoc(sdoc);
            }
        }
        _ => {
            if xargc <= flag {
                pagedoc(sdoc);
            }
        }
    }
}

/// pagedoc: page the selfdoc to stderr through $PAGER (or more), exit(1).
pub fn pagedoc(sdoc: &[&str]) -> ! {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let pager = std::env::var("PAGER").ok();
    let cmd = match pager {
        Some(p) => format!("{} 1>&2", p),
        None => "more 1>&2".to_string(),
    };
    crate::cio::fflush(&crate::cio::stdout());
    if let Ok(mut child) = Command::new("/bin/sh").arg("-c").arg(&cmd).stdin(Stdio::piped()).spawn() {
        if let Some(mut si) = child.stdin.take() {
            let mut text = String::new();
            for line in sdoc {
                text.push_str(line);
                text.push('\n');
            }
            let _ = si.write_all(text.as_bytes());
        }
        let _ = child.wait();
    }
    crate::cio::exit(crate::cio::EXIT_FAILURE)
}
