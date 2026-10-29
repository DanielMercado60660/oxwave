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

//! Port of cwp/lib/frannor.c.
//! Author: Dave Hale, Colorado School of Mines, 01/21/89.
//! Adapted from subroutine rnor in Kahaner, Moler and Nash (1988),
//! which in turn was based on an algorithm by Marsaglia and Tsang (1984).
//!
//! The functions share one process-wide generator, as in C. Each call holds a
//! mutex; callers needing a reproducible sequence must avoid interleaving calls
//! from different threads. Before seeding, the C source's default state is used.

use std::sync::Mutex;

const V: [f32; 65] = [
    0.3409450, 0.4573146, 0.5397793, 0.6062427, 0.6631691, 0.7136975, 0.7596125, 0.8020356,
    0.8417227, 0.8792102, 0.9148948, 0.9490791, 0.9820005, 1.0138492, 1.0447810, 1.0749254,
    1.1043917, 1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635, 1.2982650,
    1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929, 1.4569915, 1.4834526, 1.5100121,
    1.5367061, 1.5635712, 1.5906454, 1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045,
    1.7598422, 1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166, 1.9839239,
    2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713, 2.2245175, 2.2725185, 2.3239338,
    2.3795007, 2.4402218, 2.5075117, 2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943,
    2.7769943,
];

struct State {
    i: usize,
    j: usize,
    u: [f32; 17],
}
static STATE: Mutex<State> = Mutex::new(State {
    i: 16,
    j: 4,
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

impl State {
    // The C UNI macro advances both indices before frannor chooses k.
    fn uni(&mut self) -> f32 {
        let mut uni = self.u[self.i] - self.u[self.j];
        if uni < 0.0 {
            uni = ((uni as f64) + 1.0) as f32;
        }
        self.u[self.i] = uni;
        self.i = if self.i == 0 { 16 } else { self.i - 1 };
        self.j = if self.j == 0 { 16 } else { self.j - 1 };
        uni
    }
}

fn abs(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

/// Return a normally distributed random float with zero mean and unit variance.
///
/// C prototype: `float frannor(void);`
///
/// Advances the shared generator initialized by [`srannor`].
pub fn frannor() -> f32 {
    let mut state = STATE.lock().unwrap_or_else(|e| e.into_inner());
    let uni = state.uni();
    let vni = ((uni + uni) as f64 - 1.0) as f32;
    // C uses the integer literal 128 here, so multiplication stays in float.
    let k = ((state.u[state.i] * 128.0f32) as i32 % 64) as usize;
    // Fast part.
    let rnor = vni * V[k + 1];
    if abs(rnor) <= V[k] {
        return rnor;
    }
    // Slow part. The unsuffixed C constants and math calls operate in double
    // precision, while x, y, s, bmbx, and xnmx are rounded to float on assignment.
    let mut x = (abs(rnor) - V[k]) / (V[k + 1] - V[k]);
    let mut y = state.uni();
    let s = x + y;
    if (s as f64) <= 1.301198 {
        if (s as f64) <= 0.9689279 {
            return rnor;
        }
        let bmbx = (0.4878992 - 0.4878992 * x as f64) as f32;
        if (y as f64) <= 12.67706 - 12.37586 * (-0.5 * bmbx as f64 * bmbx as f64).exp() {
            if (-0.5 * V[k + 1] as f64 * V[k + 1] as f64).exp()
                + y as f64 * 0.01958303 / V[k + 1] as f64
                <= (-0.5 * rnor as f64 * rnor as f64).exp()
            {
                return rnor;
            }
            loop {
                y = state.uni();
                x = (0.3601016 * (y as f64).ln()) as f32;
                y = state.uni();
                if !(-2.0 * (y as f64).ln() <= (x * x) as f64) {
                    break;
                }
            }
            let xnmx = (2.776994 - x as f64) as f32;
            return if rnor >= 0.0 { abs(xnmx) } else { -abs(xnmx) };
        }
    }
    let bmbx = (0.4878992 - 0.4878992 * x as f64) as f32;
    if rnor >= 0.0 {
        abs(bmbx)
    } else {
        -abs(bmbx)
    }
}

/// Seed the shared normal random number generator.
///
/// C prototype: `void srannor(int seed);`
///
/// Resets all generator state. Seeds with the same absolute value modulo 32707
/// produce the same sequence. For `i32::MIN`, reproduce this C build's wrapping
/// absolute value; the corresponding signed negation is undefined in C itself.
pub fn srannor(seed: i32) {
    let mut state = STATE.lock().unwrap_or_else(|e| e.into_inner());
    state.i = 16;
    state.j = 4;
    let mut ia = seed.wrapping_abs() % 32707;
    let mut ib = 1111;
    let mut ic = 1947;
    for ii in 0..17 {
        let mut s = 0.0f32;
        let mut t = 0.5f32;
        for _ in 0..64 {
            let mut id = ic - ia;
            if id < 0 {
                id += 32707;
                s += t;
            }
            ia = ib;
            ib = ic;
            ic = id;
            t = (t as f64 * 0.5) as f32;
        }
        state.u[ii] = s;
    }
}
