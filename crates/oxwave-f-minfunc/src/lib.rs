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

//! Clean-room single-precision minimizers compatible with par/lib/minfunc.c.
//! Based on published bracketing, golden-section, Brent, and Powell algorithms;
//! the original C function bodies were not consulted.
//! Credits: Colorado School of Mines (2011). The original reports adaptation
//! from Numerical Recipes; this implementation uses no text from that source.

/// Evaluate `f(p + x * direction)` in single precision.
///
/// C: `float oneDfunction(float x);` (the C callback uses line-search globals).
/// Here the starting point and direction are explicit, so concurrent searches
/// do not share state.
#[allow(non_snake_case)]
pub fn oneDfunction<F: FnMut(&[f32]) -> f32>(
    x: f32,
    p: &[f32],
    direction: &[f32],
    mut f: F,
) -> f32 {
    assert_eq!(p.len(), direction.len());
    let point: Vec<f32> = p.iter().zip(direction).map(|(&a, &b)| a + x * b).collect();
    f(&point)
}

fn sign(a: f32, b: f32) -> f32 {
    if b > 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// Expand starting points `ax` and `bx` downhill until they bracket a minimum.
/// The six scalar arguments receive the three points and their function values.
/// The objective must have a finite minimum reachable by expansion.
///
/// C: `void bracket_minimum(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float));`
pub fn bracket_minimum<F: FnMut(f32) -> f32>(
    ax: &mut f32,
    bx: &mut f32,
    cx: &mut f32,
    fa: &mut f32,
    fb: &mut f32,
    fc: &mut f32,
    mut f: F,
) {
    *fa = f(*ax);
    *fb = f(*bx);
    if *fb > *fa {
        std::mem::swap(ax, bx);
        std::mem::swap(fa, fb);
    }
    *cx = *bx + 1.618034f32 * (*bx - *ax);
    *fc = f(*cx);
    while *fb > *fc {
        let r = (*bx - *ax) * (*fb - *fc);
        let q = (*bx - *cx) * (*fb - *fa);
        let mut u = (*bx as f64
            - ((*bx - *cx) * q - (*bx - *ax) * r) as f64
                / (2.0 * sign((q - r).abs().max(f32::EPSILON), q - r) as f64))
            as f32;
        let limit = (*bx as f64 + 100.0 * (*cx - *bx) as f64) as f32;
        let fu;
        if (*bx - u) * (u - *cx) > 0.0 {
            let v = f(u);
            if v < *fc {
                *ax = *bx;
                *bx = u;
                *fa = *fb;
                *fb = v;
                return;
            }
            if v > *fb {
                *cx = u;
                *fc = v;
                return;
            }
            u = *cx + 1.618034f32 * (*cx - *bx);
            fu = f(u);
        } else if (*cx - u) * (u - limit) > 0.0 {
            let v = f(u);
            if v < *fc {
                *bx = *cx;
                *cx = u;
                u = *cx + 1.618034f32 * (*cx - *bx);
                *fb = *fc;
                *fc = v;
                fu = f(u);
            } else {
                fu = v;
            }
        } else if (u - limit) * (limit - *cx) >= 0.0 {
            u = limit;
            fu = f(u);
        } else {
            u = *cx + 1.618034f32 * (*cx - *bx);
            fu = f(u);
        }
        *ax = *bx;
        *bx = *cx;
        *cx = u;
        *fa = *fb;
        *fb = *fc;
        *fc = fu;
    }
}

/// Minimize inside a bracket by golden-section interval contraction.
/// Returns the objective at `xmin`. `bx` must lie between the endpoints and
/// `tol` must be positive. As in C, termination is relative to the abscissae;
/// an exact zero minimizer can prevent convergence at subnormal spacing.
///
/// C: `float golden_bracket(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin);`
pub fn golden_bracket<F: FnMut(f32) -> f32>(
    ax: f32,
    bx: f32,
    cx: f32,
    mut f: F,
    tol: f32,
    xmin: &mut f32,
) -> f32 {
    let (mut x0, mut x3) = (ax, cx);
    let (mut x1, mut x2);
    if (cx - bx).abs() > (bx - ax).abs() {
        x1 = bx;
        x2 = bx + 0.3819660f32 * (cx - bx);
    } else {
        x2 = bx;
        x1 = bx - 0.3819660f32 * (bx - ax);
    }
    let (mut f1, mut f2) = (f(x1), f(x2));
    while (x3 - x0).abs() as f64 > tol as f64 * (x1.abs() as f64 + x2.abs() as f64) {
        if f2 < f1 {
            x0 = x1;
            x1 = x2;
            x2 = 0.61803399f32 * x1 + 0.3819660f32 * x3;
            f1 = f2;
            f2 = f(x2);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = 0.61803399f32 * x2 + 0.3819660f32 * x0;
            f2 = f1;
            f1 = f(x1);
        }
    }
    if f1 < f2 {
        *xmin = x1;
        f1
    } else {
        *xmin = x2;
        f2
    }
}

/// Minimize inside a bracket using parabolic interpolation and golden steps.
/// Returns the objective at `xmin`; `tol` is a positive relative tolerance.
///
/// C: `float brent_bracket(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin);`
pub fn brent_bracket<F: FnMut(f32) -> f32>(
    ax: f32,
    bx: f32,
    cx: f32,
    mut f: F,
    tol: f32,
    xmin: &mut f32,
) -> f32 {
    let (mut a, mut b) = if ax < cx { (ax, cx) } else { (cx, ax) };
    let (mut x, mut w, mut v) = (bx, bx, bx);
    let mut fx = f(x);
    let (mut fw, mut fv) = (fx, fx);
    let mut e = 0.0f32;
    for _ in 0..100 {
        let mid = (0.5 * (a + b) as f64) as f32;
        let t1 = (tol as f64 * x.abs() as f64 + f32::EPSILON as f64) as f32;
        let t2 = 2.0 * t1;
        if (x - mid).abs() as f64 <= t2 as f64 - 0.5 * (b - a) as f64 {
            *xmin = x;
            return fx;
        }
        let mut d;
        if e.abs() > t1 {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let old = e;
            if p.abs() as f64 >= (0.5 * q as f64 * old as f64).abs()
                || p <= q * (a - x)
                || p >= q * (b - x)
            {
                e = if x >= mid { a - x } else { b - x };
                d = 0.3819660f32 * e;
            } else {
                d = p / q;
                // SU follows an accepted interpolation with a golden step.
                e = 0.0;
                let u = x + d;
                if u - a < t2 || b - u < t2 {
                    d = sign(t1, mid - x);
                }
            }
        } else {
            e = if x >= mid { a - x } else { b - x };
            d = 0.3819660f32 * e;
        }
        let u = if d.abs() >= t1 {
            x + d
        } else {
            x + sign(t1, d)
        };
        let fu = f(u);
        if fu <= fx {
            if u >= x {
                a = x
            } else {
                b = x
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u
            } else {
                b = u
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    *xmin = x;
    fx
}

/// Minimize along `xi`, updating `p` and replacing `xi` with the displacement.
/// Both slices have length `n`; the objective value is written to `fret`.
///
/// C: `void linmin(float p[], float xi[], int n, float *fret, float (*func)());`
pub fn linmin<F: FnMut(&[f32]) -> f32>(p: &mut [f32], xi: &mut [f32], fret: &mut f32, mut f: F) {
    assert_eq!(p.len(), xi.len());
    let (mut ax, mut bx, mut cx, mut fa, mut fb, mut fc) = (0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    let mut line = |x| oneDfunction(x, p, xi, &mut f);
    bracket_minimum(
        &mut ax, &mut bx, &mut cx, &mut fa, &mut fb, &mut fc, &mut line,
    );
    let mut xmin = 0.0;
    *fret = brent_bracket(ax, bx, cx, &mut line, 2e-4, &mut xmin);
    for i in 0..p.len() {
        xi[i] *= xmin;
        p[i] += xi[i];
    }
}

/// Minimize using Powell direction searches, updating `p`, `xi`, and `fret`.
/// Directions are the columns of the square row-major matrix `xi`, with `n`
/// equal to `p.len()`. `iter` receives the number of completed sweeps.
/// Supply a nonempty point, spanning directions, and positive `ftol`.
///
/// C: `void powell_minimization(float p[], float **xi, int n, float ftol, int *iter, float *fret, float (*func)());`
pub fn powell_minimization<F: FnMut(&[f32]) -> f32>(
    p: &mut [f32],
    xi: &mut [Vec<f32>],
    ftol: f32,
    iter: &mut i32,
    fret: &mut f32,
    mut f: F,
) {
    let n = p.len();
    assert!(n > 0 && xi.len() == n && xi.iter().all(|r| r.len() == n));
    *fret = f(p);
    let mut pt = p.to_vec();
    let mut ptt = vec![0.0; n];
    let mut xit = vec![0.0; n];
    for it in 1..=200 {
        *iter = it;
        let fp = *fret;
        let mut ibig = 0;
        let mut del = 0.0;
        for i in 0..n {
            for j in 0..n {
                xit[j] = xi[j][i];
            }
            let previous = *fret;
            linmin(p, &mut xit, fret, &mut f);
            if previous - *fret > del {
                del = previous - *fret;
                ibig = i;
            }
        }
        if 2.0 * (fp - *fret).abs() as f64 <= ftol as f64 * (fp.abs() as f64 + fret.abs() as f64) {
            return;
        }
        for j in 0..n {
            ptt[j] = 2.0 * p[j] - pt[j];
            xit[j] = p[j] - pt[j];
            pt[j] = p[j];
        }
        let val = f(&ptt);
        if val < fp {
            let decrease = fp - *fret - del;
            let extrapolated = fp - val;
            let t =
                (2.0 * (fp as f64 - 2.0 * *fret as f64 + val as f64) * (decrease * decrease) as f64
                    - (del * (extrapolated * extrapolated)) as f64) as f32;
            if t < 0.0 {
                linmin(p, &mut xit, fret, &mut f);
                for j in 0..n {
                    xi[j][ibig] = xit[j];
                }
            }
        }
    }
}
