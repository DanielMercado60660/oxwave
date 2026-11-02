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

//! SUREMAC2D - Acoustic 2D Fourier modeling with REM time integration.
//! Port of su/main/synthetics_waveforms_testpatterns/suremac2d.c (Revision 1.7).
//! Credits: University of Hamburg: Ekkehart Tessmer, October 2012.
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_pfa_fft::{pfacc, Complex};
use oxwave_sinc_interp::ints8r;

pub fn stru(
    nx: usize,
    nz: usize,
    c2rho: &mut [f32],
    rhoinv: &mut [f32],
    opflag: i32,
    vel: &CFile,
    dens: Option<&CFile>,
) -> (f32, f32) {
    let mut bufz = vec![0.0f32; nz];
    for i in 0..nx {
        efread(&mut bufz, vel);
        for k in 0..nz {
            c2rho[k * nx + i] = bufz[k];
        }
    }
    match opflag {
        0 => {
            for i in 0..nx {
                efread(&mut bufz, dens.unwrap());
                for k in 0..nz {
                    rhoinv[k * nx + i] = (1.0 / bufz[k] as f64) as f32;
                    c2rho[k * nx + i] = c2rho[k * nx + i] * c2rho[k * nx + i] * bufz[k];
                }
            }
        }
        1 => {
            for k in 0..nz {
                for i in 0..nx {
                    c2rho[k * nx + i] = c2rho[k * nx + i] * c2rho[k * nx + i];
                }
            }
        }
        _ => {}
    }
    let mut vmax = 0.0f32;
    let mut vmin = 1.0e20f32;
    match opflag {
        0..=2 => {
            for k in 0..nz {
                for i in 0..nx {
                    let tmp = match opflag {
                        0 => ((c2rho[k * nx + i] * rhoinv[k * nx + i]) as f64).sqrt() as f32,
                        1 => (c2rho[k * nx + i] as f64).sqrt() as f32,
                        _ => c2rho[k * nx + i],
                    };
                    vmax = if tmp > vmax { tmp } else { vmax };
                    vmin = if tmp < vmin { tmp } else { vmin };
                }
            }
        }
        _ => {}
    }
    (vmin, vmax)
}

pub fn g2d(mut w: f32, nx: usize, nz: usize) -> Vec<f32> {
    let pi = (4.0 * 1.0f64.atan()) as f32;
    let i0 = nx as i32 / 2;
    let k0 = nz as i32 / 2;
    if (w as f64) < 0.1 {
        w = 0.1;
    }
    let sigma = ((-w * w) as f64 / (2.0 * 0.1f64.ln())).sqrt() as f32;
    let fac = (1.0 / (sigma as f64 * (2.0 * pi as f64).sqrt())) as f32;
    let mut sum = 0.0f32;
    let mut g = vec![0.0f32; nx * nz];
    for k in 0..nz {
        for i in 0..nx {
            let r2 = ((i as i32 - i0) * (i as i32 - i0) + (k as i32 - k0) * (k as i32 - k0)) as f32;
            let arg = (r2 as f64 / (2.0 * sigma as f64 * sigma as f64)) as f32;
            g[k * nx + i] = 0.0;
            if arg < 80.0 {
                g[k * nx + i] = (fac as f64 * (-(arg as f64)).exp()) as f32;
            }
            sum += g[k * nx + i];
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            g[k * nx + i] /= sum;
        }
    }
    g
}

pub fn rk(n: usize, d: f32, ind: i32) -> Vec<f32> {
    let pi = (4.0 * 1.0f64.atan()) as f32;
    let isign = -(-1.0f64).powf(ind as f32 as f64) as i32;
    let n2 = n / 2;
    let dn = (2.0 * pi as f64 / (n as f32 * d) as f64) as f32;
    let mut ak = vec![0.0f32; n];
    for i in 0..n {
        let mut c = dn * i as f32;
        if i > n2 {
            c = -dn * (n - i) as f32;
        }
        ak[i] = ((c as f64).powf(ind as f32 as f64) * isign as f64 / n as f64) as f32;
    }
    ak
}

// These work arrays replace the four independent sets of function-static
// arrays in difx, difx2, difz, and difz2 in the C program.
pub struct Derivative {
    a: Vec<Complex>,
    ak: Vec<f32>,
}
impl Derivative {
    pub fn new(n: usize, d: f32, ind: i32) -> Self {
        Self {
            a: vec![Complex::default(); n],
            ak: rk(n, d, ind),
        }
    }
}

pub fn difx(
    vin: Option<&[f32]>,
    vout: &mut [f32],
    bulk: &[f32],
    nx: usize,
    nz: usize,
    work: &mut Derivative,
    iadd: i32,
) {
    let ieo = nz % 2;
    let num = ((nz + 1) / 2) * 2;
    let a = &mut work.a;
    let rkx = &work.ak;
    for k in (0..num).step_by(2) {
        {
            // None is the C vin == vout case. Load the entire row pair
            // before overwriting it, using the existing FFT workspace.
            let vin = vin.unwrap_or(&*vout);
            for i in 0..nx {
                a[i].r = vin[k * nx + i];
                if ieo == 0 || k != num - 2 {
                    a[i].i = vin[(k + 1) * nx + i];
                } else {
                    a[i].i = 0.0;
                }
            }
        }
        pfacc(-1, a);
        for i in 0..nx {
            let tmp = a[i].r * rkx[i];
            a[i].r = -a[i].i * rkx[i];
            a[i].i = tmp;
        }
        pfacc(1, a);
        if iadd == 0 {
            for i in 0..nx {
                vout[k * nx + i] = a[i].r;
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] = a[i].i;
                }
            }
        } else if iadd == 1 {
            for i in 0..nx {
                vout[k * nx + i] += a[i].r;
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] += a[i].i;
                }
            }
        } else if iadd == 2 {
            for i in 0..nx {
                vout[k * nx + i] = a[i].r * bulk[k * nx + i];
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] = a[i].i * bulk[(k + 1) * nx + i];
                }
            }
        } else if iadd == 3 {
            for i in 0..nx {
                vout[k * nx + i] += a[i].r * bulk[k * nx + i];
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] += a[i].i * bulk[(k + 1) * nx + i];
                }
            }
        }
    }
}

pub fn difz(
    vin: &[f32],
    vout: &mut [f32],
    bulk: &[f32],
    nx: usize,
    nz: usize,
    nzpad: usize,
    work: &mut Derivative,
    iadd: i32,
    iload: i32,
) {
    let ieo = nx % 2;
    let num = ((nx + 1) / 2) * 2;
    let a = &mut work.a;
    let rkz = &work.ak;
    for i in (0..num).step_by(2) {
        let nzunl;
        if iload == 0 {
            for k in 0..nz {
                a[k].r = vin[k * nx + i];
                if ieo == 0 || i != num - 2 {
                    a[k].i = vin[k * nx + i + 1];
                } else {
                    a[k].i = 0.0;
                }
            }
            for k in nz..nzpad {
                a[k].r = 0.0;
                a[k].i = 0.0;
            }
            nzunl = nzpad;
        } else {
            for k in 0..nzpad {
                a[k].r = vin[k * nx + i];
                if ieo == 0 || i != num - 2 {
                    a[k].i = vin[k * nx + i + 1];
                } else {
                    a[k].i = 0.0;
                }
            }
            nzunl = nz;
        }
        pfacc(-1, a);
        for k in 0..nzpad {
            let tmp = a[k].r * rkz[k];
            a[k].r = -a[k].i * rkz[k];
            a[k].i = tmp;
        }
        pfacc(1, a);
        if iadd == 0 {
            for k in 0..nzunl {
                vout[k * nx + i] = a[k].r;
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] = a[k].i;
                }
            }
        } else if iadd == 1 {
            for k in 0..nzunl {
                vout[k * nx + i] += a[k].r;
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] += a[k].i;
                }
            }
        } else if iadd == 2 {
            for k in 0..nz {
                vout[k * nx + i] = a[k].r * bulk[k * nx + i];
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] = a[k].i * bulk[k * nx + i + 1];
                }
            }
            if nzunl == nzpad {
                for k in nz..nzpad {
                    vout[k * nx + i] = a[k].r * bulk[i];
                    if ieo == 0 || i != num - 2 {
                        vout[k * nx + i + 1] = a[k].i * bulk[i + 1];
                    }
                }
            }
        } else if iadd == 3 {
            for k in 0..nz {
                vout[k * nx + i] += a[k].r * bulk[k * nx + i];
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] += a[k].i * bulk[k * nx + i + 1];
                }
            }
            if nzunl == nzpad {
                for k in nz..nzpad {
                    vout[k * nx + i] += a[k].r * bulk[i];
                    if ieo == 0 || i != num - 2 {
                        vout[k * nx + i + 1] += a[k].i * bulk[i + 1];
                    }
                }
            }
        }
    }
}

pub fn bessel_jn(x: f64, n: usize, bes: &mut [f64]) {
    let eps = 1.0e-130;
    if x > n as f64 {
        err!("bessel_jn: n(%d) must be larger than x(%e)!\n", n, x);
    }
    if x.abs() > eps {
        bes[n - 1] = 0.0;
        bes[n - 2] = 1.0;
        let mut sum = 0.0;
        for i in (0..=n - 3).rev() {
            bes[i] = (2 * (i + 1)) as f64 / x * bes[i + 1] - bes[i + 2];
            if i % 2 == 0 {
                sum += 2.0 * bes[i];
                if sum.abs() > 1.0e30 {
                    for j in i..n {
                        bes[j] /= sum;
                    }
                    sum = 1.0;
                }
            }
        }
        sum = bes[0];
        for i in (2..n).step_by(2) {
            sum += 2.0 * bes[i];
        }
        for i in 0..n {
            bes[i] /= sum;
        }
    } else {
        bes[0] = 1.0;
        for i in 1..n {
            bes[i] = 0.0;
        }
    }
}

fn legpol(n: usize, x: f64) -> (f64, f64, f64) {
    let mut p0 = 1.0;
    let mut p1 = x;
    for i in 1..n {
        let tmp = ((2 * i + 1) as f64 * x * p1 - i as f64 * p0) / (i + 1) as f64;
        p0 = p1;
        p1 = tmp;
    }
    (p1, p0, n as f64 * (x * p1 - p0) / (x * x - 1.0))
}

fn null(n: usize, x0: &mut f64, eps: f64) -> (f64, f64) {
    let mut iter = 0;
    let mut d = 2.0 * eps;
    let mut p = 0.0;
    let mut pm = 0.0;
    while d > eps {
        iter += 1;
        if iter > 10 {
            fprintf!(stderr(), "stop in null: too many iterations!\n");
            exit(1);
        }
        let ps;
        (p, pm, ps) = legpol(n, *x0);
        let x1 = *x0 - p / ps;
        d = (x1 - *x0).abs();
        *x0 = x1;
    }
    (p, pm)
}

fn legcoeffs(n: usize, x: &mut [f64], a: &mut [f64]) {
    let eps = f64::EPSILON;
    let n2 = n / 2;
    let con =
        1.0 - 1.0 / (8.0 * n as f64 * n as f64) + 1.0 / (8.0 * n as f64 * n as f64 * n as f64);
    for i in 0..n2 {
        x[i] = con
            * (oxwave_cwp::cwp::PI * (4.0 * (i + 1) as f64 - 1.0) / (4.0 * n as f64 + 2.0)).cos();
    }
    for i in 0..n2 {
        let (_, pm) = null(n, &mut x[i], eps);
        a[i] = 2.0 * (1.0 - x[i] * x[i]) / ((n * n) as f64 * pm * pm);
    }
    if n % 2 != 0 {
        let i = n2;
        x[i] = 0.0;
        let (_, pm, _) = legpol(n, x[i]);
        a[i] = 2.0 / ((n * n) as f64 * pm * pm);
    }
    for i in 0..n2 {
        x[n - 1 - i] = -x[i];
        a[n - 1 - i] = a[i];
    }
}

pub fn lintgr(t: f32, r: f32, dt: f32, nt: usize, m: usize, wave: &[f32], sum: &mut [f64]) {
    let mut n = (t * r) as i32;
    n /= 2;
    if n < 20 {
        n = 20;
    }
    let n = n as usize;
    let m2 = m / 2;
    let mm = m + 50;
    let mut bes = vec![0.0f64; mm + 1];
    let mut x = vec![0.0f64; n];
    let mut a = vec![0.0f64; n];
    let mut u = vec![0.0f32; n];
    let mut xout = vec![0.0f32; n];
    let mut yout = vec![0.0f32; n];
    legcoeffs(n, &mut x, &mut a);
    for j in 0..m2 {
        sum[j] = 0.0;
    }
    for i in 0..n {
        u[i] = (0.5 * t as f64 * (x[i] + 1.0)) as f32;
        xout[i] = t - u[i];
    }
    let mut wmax = 0.0;
    for i in 0..nt {
        if (wave[i] as f64).abs() > wmax {
            wmax = (wave[i] as f64).abs();
        }
    }
    wmax *= 1.0e-6;
    ints8r(dt, dt, &wave[..nt], 0.0, 0.0, &xout, &mut yout);
    for i in 0..n {
        if (yout[i] as f64).abs() < wmax {
            continue;
        }
        let f = yout[i] as f64 * a[i];
        let arg = (u[i] * r) as f64;
        bessel_jn(arg, mm, &mut bes);
        for j in 0..m2 {
            sum[j] += bes[2 * j + 1] * f;
        }
    }
    for j in 0..m2 {
        sum[j] *= 0.5 * t as f64;
    }
}

pub fn fwave(t: f32, fmax: f32) -> f32 {
    let pi = (4.0 * 1.0f64.atan()) as f32;
    let pi2 = (2.0 * pi as f64) as f32;
    let agauss = (0.5 * fmax as f64) as f32;
    let tcut = (1.5 / agauss as f64) as f32;
    let s = (t - tcut) * agauss;
    if (s as f64).abs() < 4.0 {
        ((-2.0 * s as f64 * s as f64).exp() * ((pi2 * s) as f64).cos()) as f32
    } else {
        0.0
    }
}

pub fn damp(nbw: usize, abso: f32) -> Vec<f32> {
    let pi = (4.0 * 1.0f64.atan()) as f32;
    let delta = pi / nbw as f32;
    let mut bw = vec![0.0; nbw];
    for i in 0..nbw {
        bw[i] = (1.0 - abso as f64 * (1.0 + ((i as f32 * delta) as f64).cos()) * 0.5) as f32;
    }
    bw
}

pub fn tbc2d(a: &mut [f32], wbx: &[f32], wbz: &[f32], nwbx: i32, nwbz: i32, nx: i32, nz: i32) {
    // suremac2d.c:1853 writes before a when nwbx > nx; match fixed C.
    if nwbx > nx {
        err!("Error: nbwx must be <= nx!\n");
    }
    let mut ist = vec![0i32; nz as usize];
    let mut ng1 = nwbx;
    let mut ng2 = nwbz;
    let mut ihv = 0;
    if ng1 > ng2 {
        ng1 = nwbz;
        ng2 = nwbx;
        ihv = 1;
    }
    let mut ifl = -ng2;
    let mut istart = nwbx;
    let mut iz = 0;
    if ihv == 1 {
        ist[0] = nwbx;
    }
    for i in 0..ng2 {
        if ifl >= 0 {
            ifl -= ng2;
            iz += 1;
            if ihv == 1 {
                ist[iz] = nwbx - i;
            }
            istart -= 1;
        }
        ifl += ng1;
        if ihv == 0 {
            ist[i as usize] = istart;
        }
    }
    let nz0 = nz - 2 * nwbz;
    for k in 0..nz {
        let mut nwb = nwbx;
        if k >= nz0 + nwbz {
            nwb = ist[(nz - k - 1) as usize];
        }
        if k >= nz0 && k < nz0 + nwbz {
            nwb = ist[(k - nz0) as usize];
        }
        for i in 0..nwb {
            let wb = wbx[i as usize];
            let ind1 = i + k * nx;
            a[ind1 as usize] *= wb;
            let ind2 = nx - i - 1 + k * nx;
            a[ind2 as usize] *= wb;
        }
    }
    for k in 0..nwbz {
        let is = ist[k as usize];
        let wb = wbz[(nwbz - k - 1) as usize];
        for i in is..nx - is {
            let ind1 = i + (k + nz0) * nx;
            let ind2 = i + (nz - k - 1) * nx;
            a[ind1 as usize] *= wb;
            a[ind2 as usize] *= wb;
        }
    }
}

pub fn difx2(
    vin: &[f32],
    vout: &mut [f32],
    bulk: &[f32],
    nx: usize,
    nz: usize,
    work: &mut Derivative,
    iadd: i32,
) {
    let ieo = nz % 2;
    let num = ((nz + 1) / 2) * 2;
    let a = &mut work.a;
    let rkx = &work.ak;
    for k in (0..num).step_by(2) {
        for i in 0..nx {
            a[i].r = vin[k * nx + i];
            if ieo == 0 || k != num - 2 {
                a[i].i = vin[(k + 1) * nx + i];
            } else {
                a[i].i = 0.0;
            }
        }
        pfacc(-1, a);
        for i in 0..nx {
            a[i].r *= rkx[i];
            a[i].i *= rkx[i];
        }
        pfacc(1, a);
        if iadd == 0 {
            for i in 0..nx {
                vout[k * nx + i] = a[i].r;
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] = a[i].i;
                }
            }
        } else if iadd == 1 {
            for i in 0..nx {
                vout[k * nx + i] += a[i].r;
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] += a[i].i;
                }
            }
        } else if iadd == 2 {
            for i in 0..nx {
                vout[k * nx + i] = a[i].r * bulk[k * nx + i];
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] = a[i].i * bulk[(k + 1) * nx + i];
                }
            }
        } else if iadd == 3 {
            for i in 0..nx {
                vout[k * nx + i] += a[i].r * bulk[k * nx + i];
                if ieo == 0 || k != num - 2 {
                    vout[(k + 1) * nx + i] += a[i].i * bulk[(k + 1) * nx + i];
                }
            }
        }
    }
}

pub fn difz2(
    vin: &[f32],
    vout: &mut [f32],
    bulk: &[f32],
    nx: usize,
    nz: usize,
    nzpad: usize,
    work: &mut Derivative,
    iadd: i32,
) {
    let ieo = nx % 2;
    let num = ((nx + 1) / 2) * 2;
    let a = &mut work.a;
    let rkz = &work.ak;
    for i in (0..num).step_by(2) {
        let nzunl = nzpad;
        for k in 0..nz {
            a[k].r = vin[k * nx + i];
            if ieo == 0 || i != num - 2 {
                a[k].i = vin[k * nx + i + 1];
            } else {
                a[k].i = 0.0;
            }
        }
        for k in nz..nzpad {
            a[k].r = 0.0;
            a[k].i = 0.0;
        }
        pfacc(-1, a);
        for k in 0..nzpad {
            a[k].r *= rkz[k];
            a[k].i *= rkz[k];
        }
        pfacc(1, a);
        if iadd == 0 {
            for k in 0..nzunl {
                vout[k * nx + i] = a[k].r;
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] = a[k].i;
                }
            }
        } else if iadd == 1 {
            for k in 0..nzunl {
                vout[k * nx + i] += a[k].r;
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] += a[k].i;
                }
            }
        } else if iadd == 2 {
            for k in 0..nz {
                vout[k * nx + i] = a[k].r * bulk[k * nx + i];
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] = a[k].i * bulk[k * nx + i + 1];
                }
            }
        } else if iadd == 3 {
            for k in 0..nz {
                vout[k * nx + i] += a[k].r * bulk[k * nx + i];
                if ieo == 0 || i != num - 2 {
                    vout[k * nx + i + 1] += a[k].i * bulk[k * nx + i + 1];
                }
            }
        }
    }
}
