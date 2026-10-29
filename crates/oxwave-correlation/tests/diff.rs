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

//! Differential tests for cwp/lib/xcor.c.
//! Credits: Dave Hale, Colorado School of Mines, 11/23/91.

use oxwave_cref::XorShift;

extern "C" {
    fn xcor(lx: i32, ifx: i32, x: *const f32, ly: i32, ify: i32, y: *const f32,
            lz: i32, ifz: i32, z: *mut f32);
}

#[test]
fn xcor_matches_c() {
    let mut rng = XorShift::new(0x783c_019a);
    let sizes = [1usize, 2, 3, 7, 16, 29, 30, 31, 32, 53, 64, 127, 256, 513];
    for case in 0..12000 {
        let nx = sizes[(rng.next_u32() as usize) % sizes.len()];
        let ny = sizes[(rng.next_u32() as usize) % sizes.len()];
        let nz = sizes[(rng.next_u32() as usize) % sizes.len()];
        let ifx = rng.range_i32(-100, 100);
        let ify = rng.range_i32(-100, 100);
        let ifz = ify - ifx - nx as i32 + rng.range_i32(-40, (nx + ny + 40) as i32);
        let scale = [0.0, 1e-10, 1.0, 1e10][case % 4];
        let x: Vec<f32> = (0..nx).map(|i| if case % 7 == 0 && i % 3 == 0 { 0.0 } else { rng.uniform_f32(-scale, scale) }).collect();
        let y: Vec<f32> = (0..ny).map(|i| if case % 9 == 0 && i % 2 == 0 { -scale } else { rng.uniform_f32(-scale, scale) }).collect();
        let mut expected = vec![123.0; nz];
        let mut actual = vec![123.0; nz];
        unsafe { xcor(nx as i32, ifx, x.as_ptr(), ny as i32, ify, y.as_ptr(), nz as i32, ifz, expected.as_mut_ptr()) };
        oxwave_correlation::xcor(ifx, &x, ify, &y, ifz, &mut actual);
        let max = expected.iter().fold(0.0f32, |m, &v| m.max(v.abs()));
        for (i, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
            assert!((a - b).abs() <= 1e-5 * max, "case {case} index {i}: rust={a} C={b} max={max} nx={nx} ny={ny} nz={nz} ifx={ifx} ify={ify} ifz={ifz}");
        }
    }
}

