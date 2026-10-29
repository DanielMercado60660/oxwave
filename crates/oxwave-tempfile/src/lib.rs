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

//! Port of cwp/lib/tempname.c.
//! Credits: Andreas Klaedtke (2009); Jack K. Cohen, Colorado School of Mines (1995).

use std::ffi::{c_char, c_int, CStr};
use std::fs::File;
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;

extern "C" {
    fn mkstemp(template: *mut c_char) -> c_int;
}

// All callers supply a writable, NUL-terminated template ending in XXXXXX.
// Retain libc's modified template even on failure, just as tempname.c does.
fn create_temp(template: &mut [u8]) -> Option<File> {
    let fd = unsafe { mkstemp(template.as_mut_ptr().cast()) };
    if fd == -1 {
        None
    } else {
        // mkstemp returns a new descriptor; File takes sole ownership.
        Some(unsafe { File::from_raw_fd(fd) })
    }
}

fn filename_from_template(prefix: &str, template: &[u8], file: Option<File>) -> String {
    let path = CStr::from_bytes_with_nul(template).unwrap().to_bytes();
    if let Some(file) = file {
        // Intentional repair of the C leak: this file only supplied a name.
        // Close before reporting an unlink failure, since err! exits.
        drop(file);
        let path = std::path::Path::new(std::ffi::OsStr::from_bytes(path));
        if let Err(error) = std::fs::remove_file(path) {
            oxwave_core::err!(
                "tempname.c: temporary_filename: cannot remove %s: %s",
                path.to_string_lossy().into_owned(),
                error.to_string()
            );
        }
    }
    // Templates are private UTF-8 paths; mkstemp replaces only ASCII Xs.
    let path = std::str::from_utf8(path).unwrap();
    let slash = path.rfind('/').unwrap();
    format!("{}{}", prefix.split('\0').next().unwrap(), &path[slash..])
}

/// Generate a name using the C routine's fixed /tmp template.
///
/// C prototype: `char *temporary_filename(char *prefix)`.
///
/// The result is owned, and embedded NUL bytes terminate the prefix.
/// Unlike C, a successfully created intermediate file is closed and removed.
/// Creation failure retains libc's modified template, as in C; the returned
/// name does not reserve a file. A cleanup failure exits with an `err!` message.
/// The C implementation ignores `mkstemp`'s return value, so creation failure
/// also returns a name here rather than reporting an error.
pub fn temporary_filename(prefix: &str) -> String {
    let mut template = *b"/tmp/cmguiXXXXXX\0";
    let file = create_temp(&mut template);
    filename_from_template(prefix, &template, file)
}

/// Create an unlinked read/write temporary file in the given directory.
///
/// C prototype: `FILE *temporary_stream(char const * const prefix)`.
///
/// A missing prefix denotes the root directory. Embedded NUL bytes terminate
/// the prefix. Returns `None` if file creation fails.
pub fn temporary_stream(prefix: Option<&str>) -> Option<File> {
    let prefix = prefix.unwrap_or("").split('\0').next().unwrap();
    let mut template = format!("{prefix}/fileXXXXXX\0").into_bytes();
    let file = create_temp(&mut template)?;
    let path = CStr::from_bytes_with_nul(&template).unwrap().to_bytes();
    let _ = std::fs::remove_file(std::ffi::OsStr::from_bytes(path));
    Some(file)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filename_closes_and_removes_intermediate_files() {
        // The sandbox forbids /tmp creation. Exercise the same private helper
        // with a writable template so the success and cleanup path is tested.
        let directory = std::env::current_dir().unwrap();
        let before = std::fs::read_dir("/dev/fd").unwrap().count();
        for _ in 0..2048 {
            let mut template = format!("{}/cmguiXXXXXX\0", directory.display()).into_bytes();
            let file = create_temp(&mut template).expect("writable template must succeed");
            let name = filename_from_template("elsewhere", &template, Some(file));
            let path = std::ffi::OsStr::from_bytes(&template[..template.len() - 1]);
            assert!(!std::path::Path::new(path).exists());
            assert!(name.starts_with("elsewhere/cmgui"));
            assert_eq!(std::fs::read_dir("/dev/fd").unwrap().count(), before);
        }
    }
}
