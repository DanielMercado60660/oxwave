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

//! Port of cwp/lib/franuni.c.
//! Author: Dave Hale, Colorado School of Mines, 12/30/89.
//! Adapted from subroutine uni in Kahaner, Moler, and Nash (1988),
//! referencing unpublished notes by Marsaglia.

use std::sync::Mutex;

const CS: f64 = 362436.0 / 16777216.0;
const CD: f64 = 7654321.0 / 16777216.0;
const CM: f64 = 16777213.0 / 16777216.0;
const NBITS: usize = 24;

struct State {
    i: i32,
    j: i32,
    c: f32,
    u: [f32; 17],
}

static STATE: Mutex<State> = Mutex::new(State {
    i: 16,
    j: 4,
    c: CS as f32,
    u: [
        0.8668672834288,
        0.3697986366357,
        0.8008968294805,
        0.4173889774680,
        0.8254561579836,
        0.9640965269077,
        0.4508667414265,
        0.6451309529668,
        0.1645456024730,
        0.2787901807898,
        0.06761531340295,
        0.9663226330820,
        0.01963343943798,
        0.02947398211399,
        0.1636231515294,
        0.3976343250467,
        0.2631008574685,
    ],
});

/// Return the next pseudo-random float in [0, 1).
///
/// C prototype: `float franuni(void);`
///
/// The generator state is shared with [`sranuni`] across threads.
pub fn franuni() -> f32 {
    let mut state = STATE.lock().unwrap_or_else(|e| e.into_inner());

    // Basic generator is Fibonacci.
    let mut uni = state.u[state.i as usize] - state.u[state.j as usize];
    if uni < 0.0 {
        uni = (uni as f64 + 1.0) as f32;
    }
    let i = state.i as usize;
    state.u[i] = uni;
    state.i -= 1;
    if state.i < 0 {
        state.i = 16;
    }
    state.j -= 1;
    if state.j < 0 {
        state.j = 16;
    }

    // Second generator is congruential; CD and CM are double expressions in C.
    state.c = (state.c as f64 - CD) as f32;
    if state.c < 0.0 {
        state.c = (state.c as f64 + CM) as f32;
    }

    // Combination generator.
    uni -= state.c;
    if uni < 0.0 {
        uni = (uni as f64 + 1.0) as f32;
    }
    uni
}

/// Seed the shared pseudo-random generator.
///
/// C prototype: `void sranuni(int seed);`
///
/// Resets both generators and all 17 Fibonacci state values. Positive and
/// negative seeds with the same magnitude produce the same sequence.
/// `i32::MIN` uses wrapping absolute value to match the supplied C build;
/// the original C expression overflows for that seed.
pub fn sranuni(seed: i32) {
    let mut state = STATE.lock().unwrap_or_else(|e| e.into_inner());
    // ABS(INT_MIN) is undefined in C; wrapping matches the supplied C build.
    let magnitude = seed.wrapping_abs();
    let mut i1 = magnitude % 177 + 1;
    let mut j01 = magnitude % 167 + 1;
    let mut k1 = magnitude % 157 + 1;
    let mut l1 = magnitude % 147 + 1;
    for ii in 0..17 {
        let mut s = 0.0f32;
        let mut t = 0.5f32;
        for _ in 0..NBITS {
            let m1 = (((i1 * j01) % 179) * k1) % 179;
            i1 = j01;
            j01 = k1;
            k1 = m1;
            l1 = (53 * l1 + 1) % 169;
            if (l1 * m1) % 64 >= 32 {
                s += t;
            }
            t = (t as f64 * 0.5) as f32;
        }
        state.u[ii] = s;
    }
    state.i = 16;
    state.j = 4;
    state.c = CS as f32;
}
