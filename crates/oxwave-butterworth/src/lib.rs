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

//! Port of cwp/lib/antialias.c and helpers from cwp/lib/butterworth.c.
//! Credits: Dave Hale, Colorado School of Mines, 06/06/90 and 06/02/89.

const PI: f64 = 3.141592653589793;

fn design(fpass: f32, apass: f32, fstop: f32, astop: f32) -> (i32, f32) {
    let wpass = (2.0 * (PI * fpass as f64).tan()) as f32;
    let wstop = (2.0 * (PI * fstop as f64).tan()) as f32;
    let num = (1.0 / (apass * apass) as f64 - 1.0) / (1.0 / (astop * astop) as f64 - 1.0);
    let (fnpoles, w3db) = if fstop > fpass {
        let f = (num.ln() / ((wpass / wstop) as f64).powf(2.0).ln()) as f32;
        (f, (wpass as f64 / (1.0 / (apass * apass) as f64 - 1.0).powf(0.5 / f as f64)) as f32)
    } else {
        let f = (num.ln() / ((wstop / wpass) as f64).powf(2.0).ln()) as f32;
        (f, (wpass as f64 * (1.0 / (apass * apass) as f64 - 1.0).powf(0.5 / f as f64)) as f32)
    };
    (1i32.wrapping_add(fnpoles as i32), ((0.5 * w3db as f64).atan() / PI) as f32)
}

fn lowpass(npoles: i32, f3db: f32, q: &mut [f32]) {
    let r = (2.0 * (PI * (f3db as f64).abs()).tan()) as f32;
    if npoles % 2 != 0 {
        let scale = (r as f64 + 2.0) as f32;
        let a = r / scale;
        let b1 = ((r as f64 - 2.0) / scale as f64) as f32;
        let mut pj = 0.0f32;
        let mut qjm1 = 0.0f32;
        for v in q.iter_mut() {
            let pjm1 = pj;
            pj = *v;
            *v = a * (pj + pjm1) - b1 * qjm1;
            qjm1 = *v;
        }
    }
    for jpair in 0..npoles / 2 {
        let theta = (PI * (2 * jpair + 1) as f64 / (2 * npoles) as f64) as f32;
        let scale = (4.0 + 4.0 * r as f64 * (theta as f64).sin() + (r * r) as f64) as f32;
        let a = r * r / scale;
        let b1 = ((2.0 * r as f64 * r as f64 - 8.0) / scale as f64) as f32;
        let b2 = ((4.0 - 4.0 * r as f64 * (theta as f64).sin() + (r * r) as f64) / scale as f64) as f32;
        let (mut pjm1, mut pj, mut qjm2, mut qjm1) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
        for v in q.iter_mut() {
            let pjm2 = pjm1;
            pjm1 = pj;
            pj = *v;
            *v = (a as f64 * (pj as f64 + 2.0 * pjm1 as f64 + pjm2 as f64)
                - (b1 * qjm1) as f64 - (b2 * qjm2) as f64) as f32;
            qjm2 = qjm1;
            qjm1 = *v;
        }
    }
}

/// Anti-alias filter before increasing the sampling interval.
///
/// C prototype: `void antialias(float frac, int phase, int n, float p[], float q[]);`
/// `q` must have at least as many elements as `p`. A nonzero `phase`
/// selects minimum phase; zero selects zero phase.
pub fn antialias(frac: f32, phase: i32, p: &[f32], q: &mut [f32]) {
    assert!(q.len() >= p.len());
    let q = &mut q[..p.len()];
    let abs_frac = if frac < 0.0 { -frac } else { frac };
    if abs_frac as f64 >= 1.0 {
        q.copy_from_slice(p);
        return;
    }
    let fnyq = (0.5 * abs_frac as f64) as f32;
    let fpass = (0.6 * fnyq as f64) as f32;
    let (npoles, f3db) = design(fpass, 0.99, fnyq, 0.01);
    if phase != 0 {
        q.copy_from_slice(p);
        lowpass(npoles.wrapping_mul(2), f3db, q);
    } else {
        let mut temp = vec![0.0; p.len() + 100];
        temp[..p.len()].copy_from_slice(p);
        lowpass(npoles, f3db, &mut temp);
        temp.reverse();
        lowpass(npoles, f3db, &mut temp);
        for i in 0..p.len() { q[i] = temp[temp.len() - 1 - i]; }
    }
}
