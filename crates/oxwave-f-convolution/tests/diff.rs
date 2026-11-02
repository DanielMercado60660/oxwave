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

//! Differential tests for cwp/lib/convolution.c, by Dave Hale.
use oxwave_cref::XorShift;
use oxwave_f_convolution::convolve_cwp as rust_convolve_cwp;

extern "C" {
    fn convolve_cwp(
        lx: i32,
        ifx: i32,
        x: *const f32,
        ly: i32,
        ify: i32,
        y: *const f32,
        lz: i32,
        ifz: i32,
        z: *mut f32,
    );
}

#[test]
fn convolve_matches_c() {
    let mut rng = XorShift::new(0x163826);
    let sizes = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 37, 64, 127, 256,
    ];
    for case in 0..8000 {
        let nx = sizes[(rng.next_u32() as usize) % sizes.len()];
        let ny = sizes[(rng.next_u32() as usize) % sizes.len()];
        let ifx = (rng.next_u32() % 81) as i32 - 40;
        let ify = (rng.next_u32() % 81) as i32 - 40;
        let ifz = ifx + ify + (rng.next_u32() % ((nx + ny + 81) as u32)) as i32 - 40;
        let nz = 1 + (rng.next_u32() % ((nx + ny + 40) as u32)) as usize;
        let mut sample = || ((rng.next_u32() % 20001) as i32 - 10000) as f32 / 137.0;
        let x: Vec<f32> = (0..nx)
            .map(|_| if case % 19 == 0 { 0.0 } else { sample() })
            .collect();
        let y: Vec<f32> = (0..ny)
            .map(|_| if case % 23 == 0 { 0.0 } else { sample() })
            .collect();
        let mut expected = vec![123.0; nz];
        let mut actual = expected.clone();
        // SU's convolve_cwp does dead loads of x[-1] (rolling on) and x[lx] (even-length middle) on the
        // shorter operand when it is longer than 30; they SIGBUS if a buffer abuts an unmapped page.
        // Give the C one guard element on each side of both inputs.
        let guard = |v: &[f32]| {
            let mut g = vec![f32::NAN; v.len() + 2];
            g[1..=v.len()].copy_from_slice(v);
            g
        };
        let (xg, yg) = (guard(&x), guard(&y));
        unsafe {
            convolve_cwp(
                nx as i32,
                ifx,
                xg.as_ptr().add(1),
                ny as i32,
                ify,
                yg.as_ptr().add(1),
                nz as i32,
                ifz,
                expected.as_mut_ptr(),
            );
        }
        rust_convolve_cwp(ifx, &x, ify, &y, ifz, &mut actual);
        let scale = expected.iter().fold(0.0f32, |a, &b| a.max(b.abs()));
        for (k, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
            assert!((a-b).abs() <= 1e-5*scale, "case {case}, nx={nx} ny={ny} ifx={ifx} ify={ify} ifz={ifz} k={k}: {a} != {b} scale={scale}");
        }
    }
}
