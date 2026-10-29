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

//! Port of su/lib/elco_scalar.c, SEG-Y coordinate and elevation scalar utilities.
//! Credits: Colorado School of Mines (copyright above); no individual credits in the C source.

/// Selects a scalar for the maximum absolute coordinate.
///
/// C prototype: `short elco_scalar(int ncoords, double c[])`.
pub fn elco_scalar(c: &[f64]) -> i16 {
    let mut dtemp = 0.0;
    for &coordinate in c {
        let etemp = coordinate.abs();
        dtemp = if etemp > dtemp { etemp } else { dtemp };
    }
    if dtemp == 0.0 { return 1; }
    if dtemp >= 2147483647500.0 { return 10000; }
    if dtemp >= 214748364750.00 { return 1000; }
    if dtemp >= 21474836475.000 { return 100; }
    if dtemp >= 2147483647.5000 { return 10; }
    if dtemp >= 214748364.75000 { return 1; }
    if dtemp >= 21474836.475000 { return -10; }
    if dtemp >= 2147483.6475000 { return -100; }
    if dtemp >= 214748.36475000 { return -1000; }
    -10000
}

/// Returns the multiplier for retrieving coordinates from a SEG-Y header.
///
/// C prototype: `double from_segy_elco_multiplier(short scalar)`.
pub fn from_segy_elco_multiplier(scalar: i16) -> f64 {
    if scalar > 0 { return scalar as f64; }
    if scalar < 0 { return 1.0 / scalar as f64; }
    1.0
}

/// Returns the multiplier for storing coordinates in a SEG-Y header.
///
/// C prototype: `double to_segy_elco_multiplier(short scalar)`.
pub fn to_segy_elco_multiplier(scalar: i16) -> f64 {
    1.0 / from_segy_elco_multiplier(scalar)
}
