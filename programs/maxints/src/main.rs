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

//! MAXINTS - compute maximum and minimum integer sizes.
//! Port of cwp/main/maxints.c.
//! Credits (from the C source): CWP: Jack K. Cohen.

use oxwave_core::prelude::*;

fn maxshort() -> i16 {
    let mut i: i16 = 1;
    let mut j: i16 = 0;
    while i > 0 { j = i; i = i.wrapping_shl(1); }
    j.wrapping_sub(1).wrapping_add(j)
}

fn maxlong() -> i64 {
    let mut i: i64 = 1;
    let mut j: i64 = 0;
    while i > 0 { j = i; i = i.wrapping_shl(1); }
    j.wrapping_sub(1).wrapping_add(j)
}

fn maxint() -> i32 {
    let mut i: i32 = 1;
    let mut j: i32 = 0;
    while i > 0 { j = i; i = i.wrapping_shl(1); }
    j.wrapping_sub(1).wrapping_add(j)
}

fn maxushort() -> u16 {
    let mut i: u16 = 1;
    let mut j: u16 = 0;
    while i > j { j = i; i = i.wrapping_shl(1); }
    j.wrapping_sub(1).wrapping_add(j)
}

fn maxulong() -> u64 {
    let mut i: u64 = 1;
    let mut j: u64 = 0;
    while i > j { j = i; i = i.wrapping_shl(1); }
    j.wrapping_sub(1).wrapping_add(j)
}

fn maxuint() -> u32 {
    let mut i: u32 = 1;
    let mut j: u32 = 0;
    while i > j { j = i; i = i.wrapping_shl(1); }
    j.wrapping_sub(1).wrapping_add(j)
}

fn main() {
    printf!("max short =  %d\n", maxshort() as i32);
    printf!("min short = %d\n", -(1 + maxshort() as i32));
    printf!("max long =  %ld\n", maxlong());
    printf!("min long = %ld\n", 1i64.wrapping_add(maxlong()).wrapping_neg());
    printf!("max int =  %d\n", maxint());
    printf!("min int = %d\n", 1i32.wrapping_add(maxint()).wrapping_neg());
    printf!("max unsigned short = %u\n", maxushort() as u32);
    printf!("max unsigned long = %lu\n", maxulong());
    printf!("max unsigned int = %u\n", maxuint());
    exit(CWP_Exit());
}
