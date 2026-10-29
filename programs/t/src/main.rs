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

//! T - time and date for non-military types. Port of cwp/main/t.c.
//! Credit: Jack.

use oxwave_core::prelude::*;
use std::ffi::{c_char, c_int, c_long};

#[repr(C)]
struct Tm {
    tm_sec: c_int,
    tm_min: c_int,
    tm_hour: c_int,
    tm_mday: c_int,
    tm_mon: c_int,
    tm_year: c_int,
    tm_wday: c_int,
    tm_yday: c_int,
    tm_isdst: c_int,
    tm_gmtoff: c_long,
    tm_zone: *const c_char,
}

unsafe extern "C" {
    fn time(timer: *mut i64) -> i64;
    fn localtime(timer: *const i64) -> *mut Tm;
}

fn main() {
    let wday = ["  Sunday", "  Monday", "  Tuesday", "  Wednesday",
        "  Thursday", "  Friday", "  Saturday"];
    let mon = [" January", " Febuary", " March", " April", " May", " June",
        " July", " August", " September", " October", " November", " December"];

    let seconds = unsafe { time(std::ptr::null_mut()) };
    let timer = unsafe { &mut *localtime(&seconds) };
    let am_pm = if timer.tm_hour >= 12 { " pm" } else { " am" };
    if timer.tm_hour > 12 { timer.tm_hour -= 12; }
    timer.tm_year += 1900;
    printf!("%d:%02d%s%s, %s %d, %d\n", timer.tm_hour, timer.tm_min,
        am_pm, wday[timer.tm_wday as usize], mon[timer.tm_mon as usize],
        timer.tm_mday, timer.tm_year);
    exit(CWP_Exit());
}

