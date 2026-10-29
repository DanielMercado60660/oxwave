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

//! Quantile estimators, ported from cwp/lib/quest.c.
//! Author: Dave Hale, Colorado School of Mines, 05/07/89.
//! Algorithm: Jain, R. and Chlamtac, I., 1985, The PP algorithm for dynamic
//! calculation of quantiles and histograms without storing observations:
//! Comm. ACM, v. 28, n. 10.

/// Owned quantile estimator returned by [questalloc].
///
/// Its marker positions, heights, desired positions, and increments correspond
/// to the private QuestState in the C source. Drop it or call [questfree]
/// to release it.
#[derive(Debug)]
pub struct QuestState {
    m0: i32,
    m1: i32,
    m2: i32,
    m3: i32,
    m4: i32,
    q0: f32,
    q1: f32,
    q2: f32,
    q3: f32,
    q4: f32,
    f0: f32,
    f1: f32,
    f2: f32,
    f3: f32,
    f4: f32,
    d0: f32,
    d1: f32,
    d2: f32,
    d3: f32,
    d4: f32,
}

// Preserve insertion order for equal values, signed zeros, and NaNs.
fn initial_heights(x: &[f32]) -> [f32; 5] {
    assert!(
        x.len() >= 5,
        "quantile initialization requires at least five samples"
    );
    let mut q = [x[0], x[1], x[2], x[3], x[4]];
    for i in 1..5 {
        let mut j = i;
        while j > 0 && q[j - 1] > q[j] {
            let qtemp = q[j - 1];
            q[j - 1] = q[j];
            q[j] = qtemp;
            j -= 1;
        }
    }
    q
}

/// Estimate the quantile p of the samples x with the P-squared algorithm.
///
/// C prototype: `float quest(float p, int n, float x[]);`
///
/// As in C, p must be in [0, 1]. The slice replaces n and x.
/// Panics if fewer than five samples are provided (undefined behavior in C).
/// The result at exactly five samples is the middle sorted sample for every p.
pub fn quest(p: f32, x: &[f32]) -> f32 {
    let q = initial_heights(x);
    let (mut q0, mut q1, mut q2, mut q3, mut q4) = (q[0], q[1], q[2], q[3], q[4]);
    let (m0, mut m1, mut m2, mut m3, mut m4): (i32, i32, i32, i32, i32) = (0, 1, 2, 3, 4);
    let (mut f1, mut f2, mut f3) = (2.0 * p, 4.0 * p, 2.0 + 2.0 * p);
    let (d1, d2, d3) = (p / 2.0, p, (1.0 + p) / 2.0);
    let (mut mp, mut mm): (i32, i32);
    let (mut fm0, mut fm1, mut fm2, mut fm3, mut fm4, mut fmp, mut fmm): (
        f32,
        f32,
        f32,
        f32,
        f32,
        f32,
        f32,
    );
    let mut qtemp: f32;

    for &xtemp in &x[5..] {
        // increment marker locations and update min and max
        if xtemp < q0 {
            m1 = m1.wrapping_add(1);
            m2 = m2.wrapping_add(1);
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
            q0 = xtemp;
        } else if xtemp < q1 {
            m1 = m1.wrapping_add(1);
            m2 = m2.wrapping_add(1);
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
        } else if xtemp < q2 {
            m2 = m2.wrapping_add(1);
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
        } else if xtemp < q3 {
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
        } else if xtemp < q4 {
            m4 = m4.wrapping_add(1);
        } else {
            m4 = m4.wrapping_add(1);
            q4 = xtemp;
        }

        // increment desired marker positions
        f1 += d1;
        f2 += d2;
        f3 += d3;

        // adjust height and location of marker 1, if necessary
        mp = m1.wrapping_add(1);
        mm = m1.wrapping_sub(1);
        if f1 >= mp as f32 && m2 > mp {
            fmp = mp as f32;
            fm2 = m2 as f32;
            fm1 = m1 as f32;
            fm0 = m0 as f32;
            qtemp = q1
                + ((fmp - fm0) * (q2 - q1) / (fm2 - fm1) + (fm2 - fmp) * (q1 - q0) / (fm1 - fm0))
                    / (fm2 - fm0);
            if q0 < qtemp && qtemp < q2 {
                q1 = qtemp;
            } else {
                q1 = q1 + (q2 - q1) / (fm2 - fm1);
            }
            m1 = mp;
        } else if f1 <= mm as f32 && m0 < mm {
            fmm = mm as f32;
            fm2 = m2 as f32;
            fm1 = m1 as f32;
            fm0 = m0 as f32;
            qtemp = q1
                - ((fmm - fm0) * (q2 - q1) / (fm2 - fm1) + (fm2 - fmm) * (q1 - q0) / (fm1 - fm0))
                    / (fm2 - fm0);
            if q0 < qtemp && qtemp < q2 {
                q1 = qtemp;
            } else {
                q1 = q1 - (q0 - q1) / (fm0 - fm1);
            }
            m1 = mm;
        }

        // adjust height and location of marker 2, if necessary
        mp = m2.wrapping_add(1);
        mm = m2.wrapping_sub(1);
        if f2 >= mp as f32 && m3 > mp {
            fmp = mp as f32;
            fm3 = m3 as f32;
            fm2 = m2 as f32;
            fm1 = m1 as f32;
            qtemp = q2
                + ((fmp - fm1) * (q3 - q2) / (fm3 - fm2) + (fm3 - fmp) * (q2 - q1) / (fm2 - fm1))
                    / (fm3 - fm1);
            if q1 < qtemp && qtemp < q3 {
                q2 = qtemp;
            } else {
                q2 = q2 + (q3 - q2) / (m3.wrapping_sub(m2) as f32);
            }
            m2 = mp;
        } else if f2 <= mm as f32 && m1 < mm {
            fmm = mm as f32;
            fm3 = m3 as f32;
            fm2 = m2 as f32;
            fm1 = m1 as f32;
            qtemp = q2
                - ((fmm - fm1) * (q3 - q2) / (fm3 - fm2) + (fm3 - fmm) * (q2 - q1) / (fm2 - fm1))
                    / (fm3 - fm1);
            if q1 < qtemp && qtemp < q3 {
                q2 = qtemp;
            } else {
                q2 = q2 - (q1 - q2) / (fm1 - fm2);
            }
            m2 = mm;
        }

        // adjust height and location of marker 3, if necessary
        mp = m3.wrapping_add(1);
        mm = m3.wrapping_sub(1);
        if f3 >= mp as f32 && m4 > mp {
            fmp = mp as f32;
            fm4 = m4 as f32;
            fm3 = m3 as f32;
            fm2 = m2 as f32;
            qtemp = q3
                + ((fmp - fm2) * (q4 - q3) / (fm4 - fm3) + (fm4 - fmp) * (q3 - q2) / (fm3 - fm2))
                    / (fm4 - fm2);
            if q2 < qtemp && qtemp < q4 {
                q3 = qtemp;
            } else {
                // The C fallback intentionally uses q1 here.
                q3 = q3 + (q4 - q1) / (fm4 - fm3);
            }
            m3 = mp;
        } else if f3 <= mm as f32 && m2 < mm {
            fmm = mm as f32;
            fm4 = m4 as f32;
            fm3 = m3 as f32;
            fm2 = m2 as f32;
            qtemp = q3
                - ((fmm - fm2) * (q4 - q3) / (fm4 - fm3) + (fm4 - fmm) * (q3 - q2) / (fm3 - fm2))
                    / (fm4 - fm2);
            if q2 < qtemp && qtemp < q4 {
                q3 = qtemp;
            } else {
                q3 = q3 - (q2 - q3) / (fm2 - fm3);
            }
            m3 = mm;
        }
    }
    q2
}

/// Allocate and initialize an owned quantile estimator, processing all of x.
///
/// C prototype: `void *questalloc(float p, int n, float x[]);`
///
/// As in C, p must be in [0, 1]. The slice replaces n and x.
/// Panics if fewer than five samples are provided (undefined behavior in C).
/// Pass the returned state by mutable reference to [questupdate].
pub fn questalloc(p: f32, x: &[f32]) -> Box<QuestState> {
    let q = initial_heights(x);
    let mut s = Box::new(QuestState {
        m0: 0,
        m1: 1,
        m2: 2,
        m3: 3,
        m4: 4,
        q0: q[0],
        q1: q[1],
        q2: q[2],
        q3: q[3],
        q4: q[4],
        f0: 0.0,
        f1: 2.0 * p,
        f2: 4.0 * p,
        f3: 2.0 + 2.0 * p,
        f4: 4.0,
        d0: 0.0,
        d1: p / 2.0,
        d2: p,
        d3: (1.0 + p) / 2.0,
        d4: 1.0,
    });
    if x.len() > 5 {
        questupdate(&mut s, &x[5..]);
    }
    s
}

/// Update an estimator and return its current quantile estimate.
///
/// C prototype: `float questupdate(void *q, int n, float x[]);`
///
/// The slice replaces n and x. An empty slice returns the current estimate.
/// Integer differences are formed before conversion to f32, just as in C;
/// this differs from the float marker temporaries used by [quest].
pub fn questupdate(s: &mut QuestState, x: &[f32]) -> f32 {
    let (m0, mut m1, mut m2, mut m3, mut m4) = (s.m0, s.m1, s.m2, s.m3, s.m4);
    let (mut q0, mut q1, mut q2, mut q3, mut q4) = (s.q0, s.q1, s.q2, s.q3, s.q4);
    let (mut f0, mut f1, mut f2, mut f3, mut f4) = (s.f0, s.f1, s.f2, s.f3, s.f4);
    let (d0, d1, d2, d3, d4) = (s.d0, s.d1, s.d2, s.d3, s.d4);
    let (mut mp, mut mm): (i32, i32);
    let mut qtemp: f32;

    for &xtemp in x {
        // increment marker locations and update min and max
        if xtemp < q0 {
            m1 = m1.wrapping_add(1);
            m2 = m2.wrapping_add(1);
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
            q0 = xtemp;
        } else if xtemp < q1 {
            m1 = m1.wrapping_add(1);
            m2 = m2.wrapping_add(1);
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
        } else if xtemp < q2 {
            m2 = m2.wrapping_add(1);
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
        } else if xtemp < q3 {
            m3 = m3.wrapping_add(1);
            m4 = m4.wrapping_add(1);
        } else if xtemp < q4 {
            m4 = m4.wrapping_add(1);
        } else {
            m4 = m4.wrapping_add(1);
            q4 = xtemp;
        }

        // increment desired marker positions
        f0 += d0;
        f1 += d1;
        f2 += d2;
        f3 += d3;
        f4 += d4;

        // adjust height and location of marker 1, if necessary
        mp = m1.wrapping_add(1);
        mm = m1.wrapping_sub(1);
        if f1 >= mp as f32 && m2 > mp {
            qtemp = q1
                + ((mp.wrapping_sub(m0) as f32) * (q2 - q1) / (m2.wrapping_sub(m1) as f32)
                    + (m2.wrapping_sub(mp) as f32) * (q1 - q0) / (m1.wrapping_sub(m0) as f32))
                    / (m2.wrapping_sub(m0) as f32);
            if q0 < qtemp && qtemp < q2 {
                q1 = qtemp;
            } else {
                q1 = q1 + (q2 - q1) / (m2.wrapping_sub(m1) as f32);
            }
            m1 = mp;
        } else if f1 <= mm as f32 && m0 < mm {
            qtemp = q1
                - ((mm.wrapping_sub(m0) as f32) * (q2 - q1) / (m2.wrapping_sub(m1) as f32)
                    + (m2.wrapping_sub(mm) as f32) * (q1 - q0) / (m1.wrapping_sub(m0) as f32))
                    / (m2.wrapping_sub(m0) as f32);
            if q0 < qtemp && qtemp < q2 {
                q1 = qtemp;
            } else {
                q1 = q1 - (q0 - q1) / (m0.wrapping_sub(m1) as f32);
            }
            m1 = mm;
        }

        // adjust height and location of marker 2, if necessary
        mp = m2.wrapping_add(1);
        mm = m2.wrapping_sub(1);
        if f2 >= mp as f32 && m3 > mp {
            qtemp = q2
                + ((mp.wrapping_sub(m1) as f32) * (q3 - q2) / (m3.wrapping_sub(m2) as f32)
                    + (m3.wrapping_sub(mp) as f32) * (q2 - q1) / (m2.wrapping_sub(m1) as f32))
                    / (m3.wrapping_sub(m1) as f32);
            if q1 < qtemp && qtemp < q3 {
                q2 = qtemp;
            } else {
                q2 = q2 + (q3 - q2) / (m3.wrapping_sub(m2) as f32);
            }
            m2 = mp;
        } else if f2 <= mm as f32 && m1 < mm {
            qtemp = q2
                - ((mm.wrapping_sub(m1) as f32) * (q3 - q2) / (m3.wrapping_sub(m2) as f32)
                    + (m3.wrapping_sub(mm) as f32) * (q2 - q1) / (m2.wrapping_sub(m1) as f32))
                    / (m3.wrapping_sub(m1) as f32);
            if q1 < qtemp && qtemp < q3 {
                q2 = qtemp;
            } else {
                q2 = q2 - (q1 - q2) / (m1.wrapping_sub(m2) as f32);
            }
            m2 = mm;
        }

        // adjust height and location of marker 3, if necessary
        mp = m3.wrapping_add(1);
        mm = m3.wrapping_sub(1);
        if f3 >= mp as f32 && m4 > mp {
            qtemp = q3
                + ((mp.wrapping_sub(m2) as f32) * (q4 - q3) / (m4.wrapping_sub(m3) as f32)
                    + (m4.wrapping_sub(mp) as f32) * (q3 - q2) / (m3.wrapping_sub(m2) as f32))
                    / (m4.wrapping_sub(m2) as f32);
            if q2 < qtemp && qtemp < q4 {
                q3 = qtemp;
            } else {
                // The C fallback intentionally uses q1 here.
                q3 = q3 + (q4 - q1) / (m4.wrapping_sub(m3) as f32);
            }
            m3 = mp;
        } else if f3 <= mm as f32 && m2 < mm {
            qtemp = q3
                - ((mm.wrapping_sub(m2) as f32) * (q4 - q3) / (m4.wrapping_sub(m3) as f32)
                    + (m4.wrapping_sub(mm) as f32) * (q3 - q2) / (m3.wrapping_sub(m2) as f32))
                    / (m4.wrapping_sub(m2) as f32);
            if q2 < qtemp && qtemp < q4 {
                q3 = qtemp;
            } else {
                q3 = q3 - (q2 - q3) / (m2.wrapping_sub(m3) as f32);
            }
            m3 = mm;
        }
    }
    // Set state variables in the same order as C.
    s.m0 = m0;
    s.m1 = m1;
    s.m2 = m2;
    s.m3 = m3;
    s.m4 = m4;
    s.q0 = q0;
    s.q1 = q1;
    s.q2 = q2;
    s.q3 = q3;
    s.q4 = q4;
    s.f0 = f0;
    s.f1 = f1;
    s.f2 = f2;
    s.f3 = f3;
    s.f4 = f4;
    s.d0 = d0;
    s.d1 = d1;
    s.d2 = d2;
    s.d3 = d3;
    s.d4 = d4;
    q2
}

/// Release an estimator by consuming its ownership.
///
/// C prototype: `void questfree(void *q);`
///
/// Rust also releases the estimator automatically when its owner is dropped.
pub fn questfree(q: Box<QuestState>) {
    drop(q);
}
