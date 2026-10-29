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

//! Port of su/lib/headcase.c: trace header case lookup, reading, and writing.
//! Credits: the original SU source contains no named author credits.

use oxwave_core::header::SegyHeader;

/// Return the case number for a trace header key, or -1 for an unknown name.
///
/// C prototype: `int GetCase(char* cbuf);`
///
/// Names beginning with "null" return 0; names beginning with "numb" return
/// 1000 plus the suffix parsed as a C decimal integer. An embedded NUL ends
/// the name, as in C. Case numbers are those of the original source.
#[allow(non_snake_case)]
pub fn GetCase(cbuf: &str) -> i32 {
    let cbuf = cbuf.split('\0').next().unwrap();
    if cbuf.starts_with("null") {
        0
    } else {
        match cbuf {
            "tracl" => 1,
            "tracr" => 2,
            "fldr" => 3,
            "tracf" => 4,
            "ep" => 5,
            "cdp" => 6,
            "cdpt" => 7,
            "trid" => 8,
            "nvs" => 9,
            "nhs" => 10,
            "duse" => 11,
            "offset" => 12,
            "gelev" => 13,
            "selev" => 14,
            "sdepth" => 15,
            "gdel" => 16,
            "sdel" => 17,
            "swdep" => 18,
            "gwdep" => 19,
            "scalel" => 20,
            "scalco" => 21,
            "sx" => 22,
            "sy" => 23,
            "gx" => 24,
            "gy" => 25,
            "counit" => 26,
            "wevel" => 27,
            "swevel" => 28,
            "sut" => 29,
            "gut" => 30,
            "sstat" => 31,
            "gstat" => 32,
            "tstat" => 33,
            "laga" => 34,
            "lagb" => 35,
            "delrt" => 36,
            "muts" => 37,
            "mute" => 38,
            "ns" => 39,
            "dt" => 40,
            "gain" => 41,
            "igc" => 42,
            "igi" => 43,
            "corr" => 44,
            "sfs" => 45,
            "sfe" => 46,
            "slen" => 47,
            "styp" => 48,
            "stas" => 49,
            "stae" => 50,
            "tatyp" => 51,
            "afilf" => 52,
            "afils" => 53,
            "nofilf" => 54,
            "nofils" => 55,
            "lcf" => 56,
            "hcf" => 57,
            "lcs" => 58,
            "hcs" => 59,
            "year" => 60,
            "day" => 61,
            "hour" => 62,
            "minute" => 63,
            "sec" => 64,
            "timbas" => 65,
            "trwf" => 66,
            "grnors" => 67,
            "grnofr" => 68,
            "grnlof" => 69,
            "gaps" => 70,
            "otrav" => 71,
            "d1" => 72,
            "f1" => 73,
            "d2" => 74,
            "f2" => 75,
            "ungpow" => 76,
            "unscale" => 77,
            "ntr" => 78,
            "mark" => 79,
            "shortpad" => 80,
            _ if cbuf.starts_with("numb") => {
                1000i32.wrapping_add(oxwave_core::par::atoi(&cbuf[4..]))
            }
            _ => -1,
        }
    }
}

/// Read the raw header value selected by a case number from GetCase.
///
/// C prototype: `double fromhead(segy tr, int k);`
///
/// scalel and scalco are not applied. Accepts a header reference (also
/// obtainable by dereferencing an oxwave_core::trace::Trace), avoiding C's
/// copy of the entire trace and sample array.
///
/// # Panics
/// Panics unless k is in 1..=80. The C function returns an uninitialized
/// double outside that range and therefore has no defined result to port.
pub fn fromhead(tr: &SegyHeader, k: i32) -> f64 {
    match k {
        1 => tr.tracl as f64,
        2 => tr.tracr as f64,
        3 => tr.fldr as f64,
        4 => tr.tracf as f64,
        5 => tr.ep as f64,
        6 => tr.cdp as f64,
        7 => tr.cdpt as f64,
        8 => tr.trid as f64,
        9 => tr.nvs as f64,
        10 => tr.nhs as f64,
        11 => tr.duse as f64,
        12 => tr.offset as f64,
        13 => tr.gelev as f64,
        14 => tr.selev as f64,
        15 => tr.sdepth as f64,
        16 => tr.gdel as f64,
        17 => tr.sdel as f64,
        18 => tr.swdep as f64,
        19 => tr.gwdep as f64,
        20 => tr.scalel as f64,
        21 => tr.scalco as f64,
        22 => tr.sx as f64,
        23 => tr.sy as f64,
        24 => tr.gx as f64,
        25 => tr.gy as f64,
        26 => tr.counit as f64,
        27 => tr.wevel as f64,
        28 => tr.swevel as f64,
        29 => tr.sut as f64,
        30 => tr.gut as f64,
        31 => tr.sstat as f64,
        32 => tr.gstat as f64,
        33 => tr.tstat as f64,
        34 => tr.laga as f64,
        35 => tr.lagb as f64,
        36 => tr.delrt as f64,
        37 => tr.muts as f64,
        38 => tr.mute as f64,
        39 => tr.ns as f64,
        40 => tr.dt as f64,
        41 => tr.gain as f64,
        42 => tr.igc as f64,
        43 => tr.igi as f64,
        44 => tr.corr as f64,
        45 => tr.sfs as f64,
        46 => tr.sfe as f64,
        47 => tr.slen as f64,
        48 => tr.styp as f64,
        49 => tr.stas as f64,
        50 => tr.stae as f64,
        51 => tr.tatyp as f64,
        52 => tr.afilf as f64,
        53 => tr.afils as f64,
        54 => tr.nofilf as f64,
        55 => tr.nofils as f64,
        56 => tr.lcf as f64,
        57 => tr.hcf as f64,
        58 => tr.lcs as f64,
        59 => tr.hcs as f64,
        60 => tr.year as f64,
        61 => tr.day as f64,
        62 => tr.hour as f64,
        63 => tr.minute as f64,
        64 => tr.sec as f64,
        65 => tr.timbas as f64,
        66 => tr.trwf as f64,
        67 => tr.grnors as f64,
        68 => tr.grnofr as f64,
        69 => tr.grnlof as f64,
        70 => tr.gaps as f64,
        71 => tr.otrav as f64,
        72 => tr.d1 as f64,
        73 => tr.f1 as f64,
        74 => tr.d2 as f64,
        75 => tr.f2 as f64,
        76 => tr.ungpow as f64,
        77 => tr.unscale as f64,
        78 => tr.ntr as f64,
        79 => tr.mark as f64,
        80 => tr.shortpad as f64,
        _ => panic!("fromhead: case number must be in 1..=80"),
    }
}

/// Write a raw header value selected by a case number from GetCase.
///
/// C prototype: `void tohead(segy *tr, int k, double dvalue);`
///
/// Integer fields use nearest-integer rounding, with ties to even (the
/// default C floating-point rounding mode), followed by narrowing from C
/// long to the field type. Float fields narrow directly to float.
/// scalel and scalco are not applied. Invalid case numbers do nothing.
pub fn tohead(tr: &mut SegyHeader, k: i32, dvalue: f64) {
    match k {
        1 => tr.tracl = dvalue.round_ties_even() as i64 as i32,
        2 => tr.tracr = dvalue.round_ties_even() as i64 as i32,
        3 => tr.fldr = dvalue.round_ties_even() as i64 as i32,
        4 => tr.tracf = dvalue.round_ties_even() as i64 as i32,
        5 => tr.ep = dvalue.round_ties_even() as i64 as i32,
        6 => tr.cdp = dvalue.round_ties_even() as i64 as i32,
        7 => tr.cdpt = dvalue.round_ties_even() as i64 as i32,
        8 => tr.trid = dvalue.round_ties_even() as i64 as i16,
        9 => tr.nvs = dvalue.round_ties_even() as i64 as i16,
        10 => tr.nhs = dvalue.round_ties_even() as i64 as i16,
        11 => tr.duse = dvalue.round_ties_even() as i64 as i16,
        12 => tr.offset = dvalue.round_ties_even() as i64 as i32,
        13 => tr.gelev = dvalue.round_ties_even() as i64 as i32,
        14 => tr.selev = dvalue.round_ties_even() as i64 as i32,
        15 => tr.sdepth = dvalue.round_ties_even() as i64 as i32,
        16 => tr.gdel = dvalue.round_ties_even() as i64 as i32,
        17 => tr.sdel = dvalue.round_ties_even() as i64 as i32,
        18 => tr.swdep = dvalue.round_ties_even() as i64 as i32,
        19 => tr.gwdep = dvalue.round_ties_even() as i64 as i32,
        20 => tr.scalel = dvalue.round_ties_even() as i64 as i16,
        21 => tr.scalco = dvalue.round_ties_even() as i64 as i16,
        22 => tr.sx = dvalue.round_ties_even() as i64 as i32,
        23 => tr.sy = dvalue.round_ties_even() as i64 as i32,
        24 => tr.gx = dvalue.round_ties_even() as i64 as i32,
        25 => tr.gy = dvalue.round_ties_even() as i64 as i32,
        26 => tr.counit = dvalue.round_ties_even() as i64 as i16,
        27 => tr.wevel = dvalue.round_ties_even() as i64 as i16,
        28 => tr.swevel = dvalue.round_ties_even() as i64 as i16,
        29 => tr.sut = dvalue.round_ties_even() as i64 as i16,
        30 => tr.gut = dvalue.round_ties_even() as i64 as i16,
        31 => tr.sstat = dvalue.round_ties_even() as i64 as i16,
        32 => tr.gstat = dvalue.round_ties_even() as i64 as i16,
        33 => tr.tstat = dvalue.round_ties_even() as i64 as i16,
        34 => tr.laga = dvalue.round_ties_even() as i64 as i16,
        35 => tr.lagb = dvalue.round_ties_even() as i64 as i16,
        36 => tr.delrt = dvalue.round_ties_even() as i64 as i16,
        37 => tr.muts = dvalue.round_ties_even() as i64 as i16,
        38 => tr.mute = dvalue.round_ties_even() as i64 as i16,
        39 => tr.ns = dvalue.round_ties_even() as i64 as u16,
        40 => tr.dt = dvalue.round_ties_even() as i64 as u16,
        41 => tr.gain = dvalue.round_ties_even() as i64 as i16,
        42 => tr.igc = dvalue.round_ties_even() as i64 as i16,
        43 => tr.igi = dvalue.round_ties_even() as i64 as i16,
        44 => tr.corr = dvalue.round_ties_even() as i64 as i16,
        45 => tr.sfs = dvalue.round_ties_even() as i64 as i16,
        46 => tr.sfe = dvalue.round_ties_even() as i64 as i16,
        47 => tr.slen = dvalue.round_ties_even() as i64 as i16,
        48 => tr.styp = dvalue.round_ties_even() as i64 as i16,
        49 => tr.stas = dvalue.round_ties_even() as i64 as i16,
        50 => tr.stae = dvalue.round_ties_even() as i64 as i16,
        51 => tr.tatyp = dvalue.round_ties_even() as i64 as i16,
        52 => tr.afilf = dvalue.round_ties_even() as i64 as i16,
        53 => tr.afils = dvalue.round_ties_even() as i64 as i16,
        54 => tr.nofilf = dvalue.round_ties_even() as i64 as i16,
        55 => tr.nofils = dvalue.round_ties_even() as i64 as i16,
        56 => tr.lcf = dvalue.round_ties_even() as i64 as i16,
        57 => tr.hcf = dvalue.round_ties_even() as i64 as i16,
        58 => tr.lcs = dvalue.round_ties_even() as i64 as i16,
        59 => tr.hcs = dvalue.round_ties_even() as i64 as i16,
        60 => tr.year = dvalue.round_ties_even() as i64 as i16,
        61 => tr.day = dvalue.round_ties_even() as i64 as i16,
        62 => tr.hour = dvalue.round_ties_even() as i64 as i16,
        63 => tr.minute = dvalue.round_ties_even() as i64 as i16,
        64 => tr.sec = dvalue.round_ties_even() as i64 as i16,
        65 => tr.timbas = dvalue.round_ties_even() as i64 as i16,
        66 => tr.trwf = dvalue.round_ties_even() as i64 as i16,
        67 => tr.grnors = dvalue.round_ties_even() as i64 as i16,
        68 => tr.grnofr = dvalue.round_ties_even() as i64 as i16,
        69 => tr.grnlof = dvalue.round_ties_even() as i64 as i16,
        70 => tr.gaps = dvalue.round_ties_even() as i64 as i16,
        71 => tr.otrav = dvalue.round_ties_even() as i64 as i16,
        72 => tr.d1 = dvalue as f32,
        73 => tr.f1 = dvalue as f32,
        74 => tr.d2 = dvalue as f32,
        75 => tr.f2 = dvalue as f32,
        76 => tr.ungpow = dvalue as f32,
        77 => tr.unscale = dvalue as f32,
        78 => tr.ntr = dvalue.round_ties_even() as i64 as i32,
        79 => tr.mark = dvalue.round_ties_even() as i64 as i16,
        80 => tr.shortpad = dvalue.round_ties_even() as i64 as i16,
        _ => {}
    }
}
