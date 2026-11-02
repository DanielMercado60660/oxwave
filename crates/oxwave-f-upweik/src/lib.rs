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

//! Upwind finite difference eikonal solver, ported from `par/lib/upweik.c`.
//! Authors: CWP: Zhenuye Liu, based on code by Dave Hale, 1992.
use oxwave_tridiag::tripp;

/// Extrapolate traveltimes and their radial/angular derivatives in polar coordinates.
///
/// C prototype: `void eikpex(int na, float da, float r, float dr, float sc[], float uc[], float wc[], float tc[], float sn[], float un[], float wn[], float tn[]);`
///
/// `sc` and `sn` are current/next slownesses. `uc`, `wc`, and `tc` are current
/// radial derivatives, angular derivatives, and times; `un`, `wn`, and `tn`
/// receive their extrapolated values. The angular count is `sc.len()`.
/// Arrays must have equal lengths; at least two samples, positive r and positive da
/// are required for positive dr.
/// Nonpositive dr copies the inputs unchanged, including for a single sample.
/// Excessive curvature prints the original SU diagnostic and exits with status 1.
/// Author: Dave Hale, Colorado School of Mines, 07/16/90.
#[allow(clippy::too_many_arguments)]
pub fn eikpex(
    da: f32,
    mut r: f32,
    mut dr: f32,
    sc: &[f32],
    uc: &[f32],
    wc: &[f32],
    tc: &[f32],
    sn: &[f32],
    un: &mut [f32],
    wn: &mut [f32],
    tn: &mut [f32],
) {
    let n = sc.len();
    assert!(
        (dr <= 0.0 || n >= 2)
            && [
                uc.len(),
                wc.len(),
                tc.len(),
                sn.len(),
                un.len(),
                wn.len(),
                tn.len()
            ]
            .iter()
            .all(|&x| x == n)
    );
    let mut wtemp = vec![0.; n];
    let mut s = sc.to_vec();
    let mut left = dr;
    let orig = dr;
    un.copy_from_slice(uc);
    wn.copy_from_slice(wc);
    tn.copy_from_slice(tc);
    let wrap =
        (((n as f32 * da) as f64 - 2.0 * std::f64::consts::PI).abs()) < 0.01 * (da.abs() as f64);
    while left > 0.0 {
        assert!(r > 0.0 && da > 0.0);
        let mut cmax = 0.001f32;
        for i in 0..n {
            if ((r * un[i].abs()) as f64) < 0.001 * wn[i].abs() as f64 {
                cmax = 1000.0;
            } else {
                let x = (wn[i] / (r * un[i])).abs();
                // MAX is a conditional expression in C (including for NaN).
                cmax = if cmax > x { cmax } else { x };
            }
        }
        let step = 0.98 / cmax as f64 * r as f64 * da as f64;
        dr = if (left as f64) < step {
            left
        } else {
            step as f32
        };
        if wrap {
            let calc = |i: usize, l: usize, h: usize| {
                let umaxl = if wn[l] > 0. { un[l] } else { s[i] };
                let (uminm, umaxm) = if wn[i] > 0. {
                    (s[i], un[i])
                } else {
                    (un[i], s[i])
                };
                let uminr = if wn[h] > 0. { s[i] } else { un[h] };
                let ueol = uminm + umaxl;
                let ueor = uminr + umaxm;
                wn[i] + dr * (ueor - ueol) / da
            };
            wtemp[0] = calc(0, n - 1, 1);
            wtemp[n - 1] = calc(n - 1, n - 2, 0);
        } else {
            wtemp[0] = if wn[0] <= 0. {
                wn[0] + dr * (un[1] - un[0]) / da
            } else {
                0.
            };
            wtemp[n - 1] = if wn[n - 1] >= 0. {
                wn[n - 1] + dr * (un[n - 1] - un[n - 2]) / da
            } else {
                0.
            };
        }
        for i in 1..n - 1 {
            let umaxl = if wn[i - 1] > 0. { un[i - 1] } else { s[i] };
            let (uminm, umaxm) = if wn[i] > 0. {
                (s[i], un[i])
            } else {
                (un[i], s[i])
            };
            let uminr = if wn[i + 1] > 0. { s[i] } else { un[i + 1] };
            let ueol = uminm + umaxl;
            let ueor = uminr + umaxm;
            wtemp[i] = wn[i] + dr * (ueor - ueol) / da;
        }
        left -= dr;
        r += dr;
        let inv = (1.0 / r as f64) as f32;
        let frac = left / orig;
        for i in 0..n {
            s[i] = ((frac * sc[i]) as f64 + (1.0 - frac as f64) * sn[i] as f64) as f32;
        }
        for i in 0..n {
            wn[i] = wtemp[i];
            let wor = wn[i] * inv;
            let uu = (s[i] - wor) * (s[i] + wor);
            if uu <= 0.0 {
                oxwave_core::err!(
                    "\tRaypath has a too large curvature!\n\t A smoother velocity is required. \n"
                );
            }
            let new = (uu as f64).sqrt() as f32;
            let old = un[i];
            un[i] = new;
            tn[i] = (tn[i] as f64 + 0.5 * dr as f64 * (new + old) as f64) as f32;
        }
    }
}
/// Extrapolate ray-theoretic sigma with a zero angular derivative at both boundaries.
///
/// C prototype: `void ray_theoretic_sigma(int na, float da, float r, float dr, float uc[], float wc[], float sc[], float un[], float wn[], float sn[]);`
///
/// `uc`/`un` and `wc`/`wn` are current/next radial and angular time derivatives.
/// `sc` is current sigma and `sn` receives next sigma; the count is `sc.len()`.
/// Arrays have equal lengths of at least three.
/// Author: Zhenyue Liu, Colorado School of Mines, 07/8/92.
/// For three samples, the Rust tripp solver avoids C's access before its workspace.
#[allow(clippy::too_many_arguments)]
pub fn ray_theoretic_sigma(
    da: f32,
    r: f32,
    dr: f32,
    uc: &[f32],
    wc: &[f32],
    sc: &[f32],
    un: &[f32],
    wn: &[f32],
    sn: &mut [f32],
) {
    let n = sc.len();
    assert!(
        n >= 3
            && [uc.len(), wc.len(), un.len(), wn.len(), sn.len()]
                .iter()
                .all(|&x| x == n)
    );
    let m = n - 2;
    let mut d = vec![0.; m];
    let mut e = d.clone();
    let mut c = d.clone();
    let mut b = d.clone();
    let r1 = r + dr;
    for i in 0..m {
        d[i] = ((uc[i + 1] + un[i + 1]) as f64 / (2.0 * dr as f64)) as f32;
        e[i] = ((wn[i + 1] / (r1 * r1) + wc[i + 1] / (r * r)) as f64 / (8.0 * da as f64)) as f32;
        b[i] = (1.0 - ((sc[i + 2] - sc[i]) * e[i]) as f64 + (d[i] * sc[i + 1]) as f64) as f32;
        c[i] = -e[i];
    }
    d[0] += c[0];
    d[m - 1] += e[m - 1];
    tripp(&mut d, &mut e, &mut c, &mut b);
    sn[1..n - 1].copy_from_slice(&b);
    sn[0] = sn[1];
    sn[n - 1] = sn[n - 2];
}
/// Extrapolate incident angle with an angular derivative of one at both boundaries.
///
/// C prototype: `void ray_theoretic_beta(int na, float da, float r, float dr, float uc[], float wc[], float bc[], float un[], float wn[], float bn[]);`
///
/// `uc`/`un` and `wc`/`wn` are current/next radial and angular time derivatives.
/// `bc` is current angle and `bn` receives next angle; the count is `bc.len()`.
/// Arrays have equal lengths of at least three.
/// Author: Zhenyue Liu, Colorado School of Mines, 07/8/92.
/// For three samples, the Rust tripp solver avoids C's access before its workspace.
#[allow(clippy::too_many_arguments)]
pub fn ray_theoretic_beta(
    da: f32,
    r: f32,
    dr: f32,
    uc: &[f32],
    wc: &[f32],
    bc: &[f32],
    un: &[f32],
    wn: &[f32],
    bn: &mut [f32],
) {
    let n = bc.len();
    assert!(
        n >= 3
            && [uc.len(), wc.len(), un.len(), wn.len(), bn.len()]
                .iter()
                .all(|&x| x == n)
    );
    let m = n - 2;
    let mut d = vec![0.; m];
    let mut e = d.clone();
    let mut c = d.clone();
    let mut b = d.clone();
    let r1 = r + dr;
    for i in 0..m {
        d[i] = uc[i + 1] * r * r + un[i + 1] * r1 * r1;
        e[i] = (((wn[i + 1] + wc[i + 1]) * dr) as f64 / (4.0 * da as f64)) as f32;
        b[i] = -(bc[i + 2] - bc[i]) * e[i] + d[i] * bc[i + 1];
        c[i] = -e[i];
    }
    d[0] += c[0];
    d[m - 1] += e[m - 1];
    b[0] += da * c[0];
    b[m - 1] -= da * e[m - 1];
    tripp(&mut d, &mut e, &mut c, &mut b);
    bn[1..n - 1].copy_from_slice(&b);
    bn[0] = bn[1] - da;
    bn[n - 1] = bn[n - 2] + da;
}

/// Compute first-arrival times, propagation angles, sigma, and incident angles.
///
/// C prototype: `void eiktam(float xs, float zs, int nz, float dz, float fz, int nx, float dx, float fx, float **vel, float **time, float **angle, float **sig, float **bet);`
///
/// Dimensions come from `vel`; all four output grids must have the same shape.
/// Grids are indexed [x][z], with both dimensions at least two and nonzero spacing.
/// The source must lie within the grid. Smooth positive velocities describe the
/// physical problem; singular values propagate IEEE infinities and NaNs as in C.
/// Revisor: Zhenyue Liu, Colorado School of Mines, 7/8/92.
#[allow(clippy::too_many_arguments)]
pub fn eiktam(
    xs: f32,
    zs: f32,
    dz: f32,
    mut fz: f32,
    dx: f32,
    mut fx: f32,
    vel: &[&[f32]],
    time: &mut [&mut [f32]],
    angle: &mut [&mut [f32]],
    sig: &mut [&mut [f32]],
    bet: &mut [&mut [f32]],
) {
    use oxwave_polar_grid::{polartorect, recttopolar};
    let nx = vel.len();
    assert!(nx >= 2);
    let nz = vel[0].len();
    assert!(nz >= 2 && dx != 0. && dz != 0.);
    assert!(vel.iter().all(|x| x.len() == nz));
    for grid in [&*time, &*angle, &*sig, &*bet] {
        assert!(grid.len() == nx && grid.iter().all(|x| x.len() == nz));
    }
    fx -= xs;
    fz -= zs;
    let ex = fx + (nx as i32 - 1) as f32 * dx;
    let ez = fz + (nz as i32 - 1) as f32 * dz;
    let mut rmaxs = fx * fx + fz * fz;
    for v in [fx * fx + ez * ez, ex * ex + ez * ez, ex * ex + fz * fz] {
        rmaxs = if rmaxs > v { rmaxs } else { v };
    }
    let rmax = (rmaxs as f64).sqrt() as f32;
    let mut dr = dx.abs().min(dz.abs());
    let nr = 1 + ((rmax / dr) as f64 + 0.5) as usize;
    assert!(nr >= 2);
    dr = rmax / (nr - 1) as f32;
    let pi = std::f64::consts::PI;
    let (fa, ea) = if fx == 0. && fz == 0. {
        (0., (pi / 2.) as f32)
    } else if fx < 0. && fz == 0. {
        ((-pi / 2.) as f32, (pi / 2.) as f32)
    } else if fx == 0. && fz < 0. {
        (0., pi as f32)
    } else {
        ((-pi) as f32, pi as f32)
    };
    let mut da = dr / rmax;
    let mut na = 1 + (((ea - fa) / da) as f64 + 0.5) as usize;
    da = (ea - fa) / (na - 1) as f32;
    // The C compares float values promoted to double with PI, so this is normally false.
    if fa as f64 == -pi && ea as f64 == pi {
        na -= 1;
    }
    let mut s = vec![vec![0.; nz]; nx];
    for i in 0..nx {
        for j in 0..nz {
            s[i][j] = (1.0 / vel[i][j] as f64) as f32;
        }
    }
    let mut sp = vec![vec![0.; na]; nr];
    let mut tp = sp.clone();
    let mut up = sp.clone();
    let mut wp = sp.clone();
    let mut ap = sp.clone();
    recttopolar(
        dz,
        fz,
        dx,
        fx,
        &s.iter().map(|x| x.as_slice()).collect::<Vec<_>>(),
        da,
        fa,
        dr,
        0.,
        &mut sp.iter_mut().map(|x| x.as_mut_slice()).collect::<Vec<_>>(),
    );
    let mut ss = 0f32;
    for row in sp.iter().take(2) {
        for &v in row {
            ss += v;
        }
    }
    ss /= (2 * na) as f32;
    let mut r = 0.;
    for ir in 0..2 {
        for ia in 0..na {
            up[ir][ia] = ss;
            wp[ir][ia] = 0.;
            tp[ir][ia] = r * ss;
        }
        r += dr;
    }
    r = dr;
    for ir in 1..nr - 1 {
        let (u0, u1) = up.split_at_mut(ir + 1);
        let (w0, w1) = wp.split_at_mut(ir + 1);
        let (t0, t1) = tp.split_at_mut(ir + 1);
        eikpex(
            da,
            r,
            dr,
            &sp[ir],
            &u0[ir],
            &w0[ir],
            &t0[ir],
            &sp[ir + 1],
            &mut u1[0],
            &mut w1[0],
            &mut t1[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        0.,
        &tp.iter().map(|x| x.as_slice()).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        time,
    );
    let mut a = fa;
    for ia in 0..na {
        ap[0][ia] = a;
        a += da;
    }
    r = dr;
    for ir in 1..nr {
        a = fa;
        for ia in 0..na {
            ap[ir][ia] = (a as f64 + ((wp[ir][ia] / (sp[ir][ia] * r)) as f64).asin()) as f32;
            a += da;
        }
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        0.,
        &ap.iter().map(|x| x.as_slice()).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        angle,
    );
    r = 0.;
    for row in tp.iter_mut().take(2) {
        for v in row {
            *v = r / ss;
        }
        r += dr;
    }
    r = dr;
    for ir in 1..nr - 1 {
        let (t0, t1) = tp.split_at_mut(ir + 1);
        ray_theoretic_sigma(
            da,
            r,
            dr,
            &up[ir],
            &wp[ir],
            &t0[ir],
            &up[ir + 1],
            &wp[ir + 1],
            &mut t1[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        0.,
        &tp.iter().map(|x| x.as_slice()).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        sig,
    );
    for row in tp.iter_mut().take(2) {
        a = fa;
        for v in row {
            *v = a;
            a += da;
        }
    }
    r = dr;
    for ir in 1..nr - 1 {
        let (t0, t1) = tp.split_at_mut(ir + 1);
        ray_theoretic_beta(
            da,
            r,
            dr,
            &up[ir],
            &wp[ir],
            &t0[ir],
            &up[ir + 1],
            &wp[ir + 1],
            &mut t1[0],
        );
        r += dr;
    }
    polartorect(
        da,
        fa,
        dr,
        0.,
        &tp.iter().map(|x| x.as_slice()).collect::<Vec<_>>(),
        dz,
        fz,
        dx,
        fx,
        bet,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eiktam_accepts_negative_depth_spacing() {
        let vel = vec![vec![2000.0; 9]; 11];
        let mut time = vec![vec![0.0; 9]; 11];
        let mut angle = time.clone();
        let mut sig = time.clone();
        let mut bet = time.clone();
        eiktam(
            500.0,
            0.0,
            -100.0,
            0.0,
            100.0,
            0.0,
            &vel.iter().map(Vec::as_slice).collect::<Vec<_>>(),
            &mut time.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
            &mut angle.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
            &mut sig.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
            &mut bet.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>(),
        );
        assert_eq!(time[5][0], 0.0);
        for grid in [&time, &angle, &sig, &bet] {
            assert!(grid.iter().flatten().all(|v| v.is_finite()));
        }
    }
}
