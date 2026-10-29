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

// Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.

//! DOWNFORT - change Fortran programs to lower case, preserving strings.
//! Port of cwp/main/downfort.c.
//! Credits: Brian Sumner, c. 1984.

use oxwave_core::prelude::*;

fn lower(c: i32) -> i32 {
    if (b'A' as i32..=b'Z' as i32).contains(&c) { c + 32 } else { c }
}

fn main() {
    let input = stdin();
    loop {
        let mut c = fgetc(&input);
        if c == -1 { break; }
        if c == b'\n' as i32 {
            putchar(c);
        } else if c == b'*' as i32 || c == b'C' as i32 {
            putchar(c);
            loop {
                c = fgetc(&input);
                if c == b'\n' as i32 { break; }
                putchar(c);
            }
            putchar(c);
        } else {
            putchar(c);
            loop {
                c = fgetc(&input);
                if c == b'\n' as i32 { break; }
                if c != b'\'' as i32 {
                    putchar(lower(c));
                } else {
                    putchar(c);
                    c = fgetc(&input);
                    if c != b'(' as i32 {
                        putchar(c);
                        loop {
                            c = fgetc(&input);
                            if c == b'\'' as i32 { break; }
                            putchar(c);
                        }
                        putchar(c);
                    } else {
                        putchar(c);
                        loop {
                            c = fgetc(&input);
                            if c == b'\'' as i32 { break; }
                            putchar(lower(c));
                        }
                        putchar(c);
                    }
                }
            }
            putchar(c);
        }
    }
    exit(CWP_Exit());
}
