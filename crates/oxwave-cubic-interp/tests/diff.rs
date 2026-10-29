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

//! Differential tests for cwp/lib/intcub.c.
//! Author: Dave Hale, Colorado School of Mines, 06/02/89.

use oxwave_cref::XorShift;
use oxwave_cubic_interp::intcub;

mod c {
    extern "C" {
        // void intcub(int ideriv, int nin, float xin[], float ydin[][4],
        //             int nout, float xout[], float yout[]);
        pub fn intcub(
            ideriv: i32,
            nin: i32,
            xin: *const f32,
            ydin: *const [f32; 4],
            nout: i32,
            xout: *const f32,
            yout: *mut f32,
        );
    }
}

fn compare(ideriv: i32, xin: &[f32], ydin: &[[f32; 4]], xout: &[f32], case: usize) {
    // The trailing sentinel also checks that nout, rather than yout.len(),
    // determines how many elements are written.
    let mut expected = vec![123.0; xout.len() + 1];
    let mut actual = expected.clone();
    // SAFETY: both implementations receive matching lengths and valid buffers.
    // For empty xin, the cases below either have nout == 0 or an unsupported
    // derivative order, so the C implementation never dereferences xin/ydin.
    unsafe {
        c::intcub(
            ideriv,
            xin.len() as i32,
            xin.as_ptr(),
            ydin.as_ptr(),
            xout.len() as i32,
            xout.as_ptr(),
            expected.as_mut_ptr(),
        );
    }
    intcub(ideriv, xin, ydin, xout, &mut actual);
    for (i, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "n={} nout={} case={case} ideriv={ideriv} i={i} rust={a:?} c={b:?}",
            xin.len(),
            xout.len(),
        );
    }
}

#[test]
fn intcub_matches_c() {
    let mut rng = XorShift::new(0x4355424943);
    // 2,112 random tables, each tested for every supported derivative and
    // invalid orders. A single test keeps C's static search hint race-free.
    for &n in &[1usize, 2, 3, 7, 16, 31, 64, 257, 1024, 4093, 4096] {
        for case in 0..192 {
            let xscale = [0.0001, 1.0, 100.0][case % 3];
            let yscale = [1.0e-30, 1.0e-10, 1.0, 1.0e10, 1.0e15][case % 5];
            let mut xin = Vec::with_capacity(n);
            let mut x = rng.uniform_f32(-50.0, 0.0) * xscale;
            for _ in 0..n {
                x += rng.uniform_f32(0.05, 2.0) * xscale;
                xin.push(x);
            }
            if case % 2 != 0 {
                xin.reverse();
            }
            let mut ydin: Vec<[f32; 4]> = (0..n)
                .map(|_| {
                    [
                        rng.uniform_f32(-100.0, 100.0) * yscale,
                        rng.uniform_f32(-10.0, 10.0) * yscale,
                        rng.uniform_f32(-2.0, 2.0) * yscale,
                        rng.uniform_f32(-1.0, 1.0) * yscale,
                    ]
                })
                .collect();
            if case % 5 == 0 {
                ydin[0] = [0.0, -0.0, 0.0, 0.0];
            }
            if case % 19 == 0 {
                ydin.fill([0.0; 4]);
            }
            let lo = xin[0].min(xin[n - 1]) - 10.0 * xscale;
            let hi = xin[0].max(xin[n - 1]) + 10.0 * xscale;
            let nout = [0, 1, 2, 17, 64, 257][case % 6];
            let mut xout = rng.vec_f32(nout, lo, hi);
            xout.extend([xin[0], xin[n - 1], 0.0, lo, hi]);
            // Probe both neighboring representable values at an interior knot.
            let knot = xin[n / 2];
            xout.extend([knot.next_down(), knot, knot.next_up()]);
            match case % 3 {
                0 => xout.sort_by(f32::total_cmp),
                1 => xout.sort_by(|a, b| b.total_cmp(a)),
                _ => (),
            }
            for ideriv in [-2, -1, 0, 1, 2, 3, 4, i32::MIN, i32::MAX] {
                compare(ideriv, &xin, &ydin, &xout, case);
            }
            compare(0, &xin, &ydin, &[], case);
        }
    }
    // Signed zero and subnormal coefficients, without undefined C inputs.
    for &x in &[0.0, -0.0] {
        for y in [
            [0.0; 4],
            [-0.0; 4],
            [f32::from_bits(1), -f32::from_bits(1), 0.0, -0.0],
        ] {
            for ideriv in 0..=3 {
                compare(ideriv, &[x], &[y], &[0.0, -0.0, -1.0, 1.0], 0);
            }
        }
    }
    for ideriv in [-1, 4, 8, i32::MIN, i32::MAX] {
        compare(ideriv, &[], &[], &[0.0, -1.0], 0);
    }
    for ideriv in -1..=4 {
        compare(ideriv, &[], &[], &[], 0);
    }
}
