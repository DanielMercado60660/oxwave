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

//! Differential tests for cwp/lib/xindex.c.
//! Author: Dave Hale, Colorado School of Mines, 12/25/89.

use oxwave_cref::XorShift;

mod c {
    extern "C" {
        // void xindex(int nx, float ax[], float x, int *index);
        pub fn xindex(nx: i32, ax: *mut f32, x: f32, index: *mut i32);
    }
}

fn compare(ax: &[f32], x: f32, initial: i32) -> i32 {
    let mut expected = initial;
    let mut actual = initial;
    let mut c_ax = ax.to_vec();
    // C only reads ax; the copy keeps its non-const prototype exact and safe.
    unsafe { c::xindex(c_ax.len() as i32, c_ax.as_mut_ptr(), x, &mut expected) };
    oxwave_index_search::xindex(ax, x, &mut actual);
    assert_eq!(
        actual,
        expected,
        "n={}, x={x:?}, initial={initial}, ax={ax:?}",
        ax.len()
    );
    actual
}

#[test]
fn xindex_float_boundaries_match_c() {
    let mut rng = XorShift::new(0x9876_abcd);
    for n in [1usize, 2, 17, 128, 1021, 8192] {
        // Sorting arbitrary finite IEEE values covers subnormals and the full
        // exponent range without overflowing arithmetic used to build ax.
        let mut ax: Vec<f32> = (0..n)
            .map(|_| {
                let bits = rng.next_u32();
                f32::from_bits((bits & 0x8000_0000) | ((bits & 0x7fff_ffff) % 0x7f80_0000))
            })
            .collect();
        ax.sort_by(f32::total_cmp);
        for _ in 0..2 {
            let mut previous = 0;
            for k in 0..2000 {
                let value = ax[rng.next_u32() as usize % n];
                let x = match k % 6 {
                    0 => value,
                    1 => value.next_up(),
                    2 => value.next_down(),
                    3 => f32::from_bits(1),
                    4 => -f32::from_bits(1),
                    _ => rng.any_f32(),
                };
                let initial = if k % 2 == 0 {
                    previous
                } else {
                    rng.range_i32(-1, n as i32)
                };
                previous = compare(&ax, x, initial);
            }
            ax.reverse();
        }
    }

    // Duplicate runs, signed zeros, endpoints, and queries moving in both
    // directions exercise reuse of the output index as the next search hint.
    for ax in [
        vec![-4.0, -4.0, -1.0, 0.0, -0.0, 0.0, 3.0, 3.0],
        vec![3.0, 3.0, 0.0, -0.0, 0.0, -1.0, -4.0, -4.0],
        vec![0.0; 17],
        vec![
            f32::MIN,
            -f32::MIN_POSITIVE,
            0.0,
            f32::MIN_POSITIVE,
            f32::MAX,
        ],
    ] {
        for initial in -1..=ax.len() as i32 {
            let mut previous = initial;
            for &value in ax.iter().chain(ax.iter().rev()) {
                for x in [value.next_down(), value, value.next_up()] {
                    previous = compare(&ax, x, previous);
                }
            }
        }
    }
}

#[test]
fn xindex_matches_c() {
    let mut rng = XorShift::new(0x1234_5678);
    for n in [1usize, 2, 3, 7, 16, 31, 64, 127, 256, 1009, 4096] {
        for direction in 0..3 {
            for _ in 0..12 {
                let mut ax: Vec<f32> = (0..n)
                    .map(|_| rng.range_i32(-1000, 1000) as f32 * 0.25)
                    .collect();
                ax.sort_by(|a, b| a.partial_cmp(b).unwrap());
                if direction == 1 {
                    ax.reverse();
                }
                if direction == 2 {
                    ax.fill(rng.range_i32(-5, 5) as f32);
                }
                let mut previous = 0;
                for k in 0..160 {
                    let x = match k % 8 {
                        0 => 0.0,
                        1 => ax[rng.next_u32() as usize % n],
                        2 => f32::NEG_INFINITY,
                        3 => f32::INFINITY,
                        4 => f32::NAN,
                        _ => rng.uniform_f32(-300.0, 300.0),
                    };
                    let initial = match k % 7 {
                        0 => i32::MIN,
                        1 => i32::MAX,
                        2 => -1,
                        3 => n as i32,
                        4 => rng.range_i32(0, n as i32 - 1),
                        _ => previous,
                    };
                    previous = compare(&ax, x, initial);
                }
            }
        }
    }
    for ax in [
        &[0.0, -0.0, 0.0][..],
        &[-f32::INFINITY, 0.0, f32::INFINITY][..],
        &[f32::INFINITY, 0.0, -f32::INFINITY][..],
    ] {
        for x in [0.0, -0.0, f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            for initial in [-1, 0, 1, 2, 3] {
                compare(ax, x, initial);
            }
        }
    }
}
