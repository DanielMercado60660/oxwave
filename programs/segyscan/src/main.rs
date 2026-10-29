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
//! SEGYSCAN - scan SEG-Y trace headers for minimum and maximum values.
//! Port of su/main/data_conversion/segyscan.c (Revision 1.3).
//! Credits (from the C source): Stew Levin, June 2013.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SEGYSCAN -- SCANs SEGY file trace headers for min-max in  several\t",
    "     possible formats.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   segyscan < segyfile\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The SEGY file trace headers are scanned assuming short, ushort, int,  ",
    " uint, float, and double and the results are printed as tables.\t",
    " \t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut inhdr = [0u8; 400];
    let mut s2min = [i16::MAX; 120];
    let mut s2max = [i16::MIN; 120];
    let mut u2min = [u16::MAX; 120];
    let mut u2max = [0u16; 120];
    let mut s4min = [i32::MAX; 60];
    let mut s4max = [i32::MIN; 60];
    let mut u4min = [u32::MAX; 60];
    let mut u4max = [0u32; 60];
    let mut f4min = [f32::MAX; 60];
    let mut f4max = [-f32::MAX; 60];
    let mut f8min = [f64::MAX; 30];
    let mut f8max = [-f64::MAX; 30];
    let input = stdin();
    for _ in 0..3200 { fgetc(&input); }
    fread(&mut inhdr, &input);
    for i in (0..26).step_by(2) { inhdr.swap(i, i + 1); }
    let mut trlen = u16::from_ne_bytes([inhdr[20], inhdr[21]]) as usize;
    let format = u16::from_ne_bytes([inhdr[24], inhdr[25]]);
    match format { 1 | 2 | 4 | 5 => trlen *= 4, 3 => trlen *= 2, _ => {} }
    loop {
        let nread = fread(&mut inhdr[..240], &input);
        for _ in 0..trlen { fgetc(&input); }
        if nread < 240 { break; }
        // IMLITTLEENDIAN in the C source tests an array address, always true.
        for i in (0..240).step_by(2) { inhdr.swap(i, i + 1); }
        for i in 0..120 {
            let u = u16::from_ne_bytes([inhdr[2*i], inhdr[2*i+1]]);
            let s = u as i16;
            if s2min[i] > s { s2min[i] = s; }
            if s2max[i] < s { s2max[i] = s; }
            if u2min[i] > u { u2min[i] = u; }
            if u2max[i] < u { u2max[i] = u; }
        }
        for i in (0..240).step_by(4) { inhdr.swap(i, i+2); inhdr.swap(i+1, i+3); }
        for i in 0..60 {
            let u = u32::from_ne_bytes(inhdr[4*i..4*i+4].try_into().unwrap());
            let s = u as i32;
            let f = f32::from_bits(u);
            if s4min[i] > s { s4min[i] = s; }
            if s4max[i] < s { s4max[i] = s; }
            if u4min[i] > u { u4min[i] = u; }
            if u4max[i] < u { u4max[i] = u; }
            if f4min[i] > f { f4min[i] = f; }
            if f4max[i] < f { f4max[i] = f; }
        }
        for i in (0..240).step_by(8) { for j in 0..4 { inhdr.swap(i+j, i+j+4); } }
        for i in 0..30 {
            let f = f64::from_ne_bytes(inhdr[8*i..8*i+8].try_into().unwrap());
            if f8min[i] > f { f8min[i] = f; }
            if f8max[i] < f { f8max[i] = f; }
        }
    }
    printf!("\nbyte#|     short      |     ushort     \n");
    printf!("---------------------------------------\n");
    for i in 0..120 {
        printf!("%4d | %6hd - %6hd | %6hu - %6hu\n", (1+i*2) as i32, s2min[i] as i32, s2max[i] as i32, u2min[i] as u32, u2max[i] as u32);
    }
    printf!("\nbyte#|           int            |           uint           |         float    \n");
    printf!("-------------------------------------------------------------------------------------\n");
    for i in 0..60 {
        printf!("%4d | %11d - %11d | %11u - %11u | %g - %g\n", (1+i*4) as i32, s4min[i], s4max[i], u4min[i], u4max[i], f4min[i] as f64, f4max[i] as f64);
    }
    printf!("\nbyte#|              double \n");
    printf!("------------------------------------------\n");
    for i in 0..30 { printf!("%4d | %g - %g\n", (1+i*8) as i32, f8min[i], f8max[i]); }
    exit(CWP_Exit());
}
