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

//! Port of par/lib/taup.c. Credits: John Anderson, Gabriel Alvarez,
//! Bjoern E. Rommel, and Dave Hale.

use oxwave_cwp::complex::{cadd, cmul, conjg, crmul, csub, Complex};

/// C: `float gofx(int igopt, float offset, float intercept_off, float refdepth)`
pub fn gofx(igopt: i32, offset: f32, intercept_off: f32, refdepth: f32) -> f32 {
    let offset = offset - intercept_off;
    match igopt {
        1 => offset * offset,
        2 => ((refdepth * refdepth + offset * offset) as f64).sqrt() as f32,
        3 => offset,
        4 => (offset as f64).abs() as f32,
        5 => (((refdepth * refdepth + offset * offset) as f64).sqrt() - refdepth as f64) as f32,
        _ => 0.0,
    }
}
/// C: `float freqweight(int j, float df, float f1, float f2)`
pub fn freqweight(j: i32, df: f32, f1: f32, f2: f32) -> f32 {
    let f = j as f32 * df;
    if f <= f1 {
        return 1.0;
    }
    if f >= f2 {
        return 0.0;
    }
    (f2 - f) / (f2 - f1)
}
/// C: `void compute_r(float w, int nx, float *g, int np, float dp, complex *r)`
pub fn compute_r(w: f32, g: &[f32], dp: f32, r: &mut [Complex]) {
    let fac = (1.0 / (g.len() * r.len()) as f64) as f32;
    for j in 0..r.len() {
        let mut rs = 0.0f32;
        let mut is = 0.0f32;
        for &x in g {
            let a = w * j as f32 * dp * x;
            rs = (rs as f64 + (a as f64).cos()) as f32;
            is = (is as f64 + (a as f64).sin()) as f32;
        }
        r[j] = Complex {
            r: fac * rs,
            i: fac * is,
        };
    }
}
/// C: `void compute_rhs(float w, int nx, float *g, complex *data, int np, float pmin, float dp, complex *rhs)`
pub fn compute_rhs(w: f32, g: &[f32], data: &[Complex], pmin: f32, dp: f32, rhs: &mut [Complex]) {
    let fac = (1.0 / g.len() as f64) as f32;
    for ip in 0..rhs.len() {
        let p = pmin + ip as f32 * dp;
        let mut rs = 0.0f32;
        let mut is = 0.0f32;
        for ix in 0..g.len() {
            let angle = (w * p * g[ix]) as f64;
            let tr = angle.cos() as f32;
            let ti = -angle.sin() as f32;
            let dr = data[ix].r;
            let di = data[ix].i;
            rs += tr * dr - ti * di;
            is += tr * di + ti * dr;
        }
        rhs[ip] = Complex {
            r: fac * rs,
            i: fac * is,
        };
    }
}
/// C: `float rcdot(int n, complex *a, complex *b)`
pub fn rcdot(a: &[Complex], b: &[Complex]) -> f32 {
    let mut sum = 0.0;
    for j in 0..a.len() {
        sum += a[j].r * b[j].r + a[j].i * b[j].i;
    }
    sum
}
/// C: `int ctoep(int n, complex *r, complex *a, complex *b, complex *f, complex *g)`
pub fn ctoep(
    r: &[Complex],
    a: &mut [Complex],
    b: &[Complex],
    f: &mut [Complex],
    g: &mut [Complex],
) -> i32 {
    let n = r.len();
    if r[0].r == 0.0 {
        return 0;
    }
    f[0] = Complex {
        r: (1.0 / r[0].r as f64) as f32,
        i: 0.0,
    };
    a[0] = Complex {
        r: b[0].r / r[0].r,
        i: b[0].i / r[0].r,
    };
    let (mut vr, mut vi, mut vsq) = (1.0f32, 0.0f32, 1.0f32);
    let mut j = 1;
    while j < n {
        f[j] = Complex::default();
        let (mut er, mut ei) = (0.0f32, 0.0f32);
        for k in 0..j {
            let t = j - k;
            er += r[t].r * f[k].r + r[t].i * f[k].i;
            ei += r[t].r * f[k].i - r[t].i * f[k].r;
        }
        let mut cr = (er * vr - ei * vi) / vsq;
        let mut ci = (er * vi + ei * vr) / vsq;
        vr = vr - (cr * er + ci * ei);
        vi = vi + (cr * ei - ci * er);
        vsq = vr * vr + vi * vi;
        if vsq <= 0.0 {
            break;
        }
        for k in 0..=j {
            let t = j - k;
            g[k].r = f[k].r - cr * f[t].r - ci * f[t].i;
            g[k].i = f[k].i + cr * f[t].i - ci * f[t].r;
        }
        f[..=j].copy_from_slice(&g[..=j]);
        a[j] = Complex::default();
        er = 0.0;
        ei = 0.0;
        for k in 0..j {
            let t = j - k;
            er += r[t].r * a[k].r + r[t].i * a[k].i;
            ei += r[t].r * a[k].i - r[t].i * a[k].r;
        }
        er -= b[j].r;
        ei -= b[j].i;
        cr = (er * vr - ei * vi) / vsq;
        ci = (er * vi + ei * vr) / vsq;
        for k in 0..=j {
            let t = j - k;
            a[k].r += -cr * f[t].r - ci * f[t].i;
            a[k].i += cr * f[t].i - ci * f[t].r;
        }
        j += 1;
    }
    vr = 0.0;
    vi = 0.0;
    for k in 0..j {
        vr += r[k].r * f[k].r - r[k].i * f[k].i;
        vi += r[k].r * f[k].i + r[k].i * f[k].r;
    }
    vsq = vr * vr + vi * vi;
    let er = vr / vsq;
    let ei = -vi / vsq;
    for k in 0..j {
        f[k].r = er * f[k].r - ei * f[k].i;
        f[k].i = er * f[k].i + ei * f[k].r;
    }
    j as i32
}
/// C: `void htmul(int n, complex *a, complex *x, complex *y)`
pub fn htmul(a: &[Complex], x: &[Complex], y: &mut [Complex]) {
    let n = y.len();
    for i in 0..n {
        y[i] = Complex::default();
        for j in 0..i {
            y[i] = cadd(cmul(conjg(a[i - j]), x[j]), y[i]);
        }
        for j in i..n {
            y[i] = cadd(cmul(a[j - i], x[j]), y[i]);
        }
    }
}
/// C: `void runav(int n, int len, float *a, float *b)`
pub fn runav(len: i32, a: &[f32], b: &mut [f32]) {
    let n = a.len();
    let lenh = (len / 2) as usize;
    let mut sum = 0.0f32;
    if len <= 1 {
        b[..n].copy_from_slice(a);
        return;
    }
    for j in 0..(len as usize).min(n) {
        sum += a[j];
    }
    for j in 0..lenh.min(n) {
        b[j] = sum;
    }
    for j in lenh..n.saturating_sub(lenh) {
        sum = sum + a[j + lenh] - a[j - lenh];
        b[j] = sum;
    }
    for j in n.saturating_sub(lenh)..n {
        b[j] = sum;
    }
    for v in &mut b[..n] {
        *v /= len as f32;
    }
}
/// C: `void taupmute(int ip, int ipa, int ipb, int nt, int ltap, float *rt)`
pub fn taupmute(ip: i32, ipa: i32, ipb: i32, ltap: i32, rt: &mut [f32]) {
    if ip <= ipa {
        return;
    }
    if ip >= ipb {
        rt.fill(0.0);
        return;
    }
    let mut w = (ipb - ipa).max(1) as f32;
    w = (ipb - ip) as f32 / w;
    let j = (w * rt.len() as f32) as usize;
    for k in 0..j {
        rt[k] = 0.0;
    }
    for k in j..rt.len().min(j + ltap.max(0) as usize) {
        rt[k] *= ((k - j) as i32 / ltap) as f32;
    }
}
/// C: `int ctoephcg(int niter, int n, complex *a, complex *x, complex *y, complex *s, complex *ss, complex *g, complex *rr)`
pub fn ctoephcg(
    niter: i32,
    a: &[Complex],
    x: &mut [Complex],
    y: &[Complex],
    s: &mut [Complex],
    ss: &mut [Complex],
    g: &mut [Complex],
    rr: &mut [Complex],
) -> i32 {
    let n = y.len();
    let eps = 1.0e-6f32;
    let mut rp = rcdot(y, y);
    let test = n as f32 * eps * eps * rp;
    x[..n].fill(Complex::default());
    rr[..n].copy_from_slice(y);
    htmul(a, rr, g);
    s[..n].copy_from_slice(&g[..n]);
    let mut gammam = rcdot(g, g);
    for iter in 0..niter {
        htmul(a, s, ss);
        let alpha = gammam / rcdot(ss, ss);
        for j in 0..n {
            x[j] = cadd(x[j], crmul(s[j], alpha));
            rr[j] = csub(rr[j], crmul(ss[j], alpha));
        }
        let rsq = rcdot(rr, rr);
        if iter > 0 && (rsq == rp || rsq < test) {
            return iter - 1;
        }
        rp = rsq;
        htmul(a, rr, g);
        let gamma = rcdot(g, g);
        if gamma < eps {
            return iter;
        }
        let beta = gamma / gammam;
        gammam = gamma;
        for j in 0..n {
            s[j] = cadd(g[j], crmul(s[j], beta));
        }
    }
    niter
}

/// C: `void fwd_tx_sstack(float dt, int nt, int nx, float xmin, float dx, int np, float pmin, float dp, float **traces, float **out_traces)`
#[allow(non_snake_case)]
pub fn fwd_tx_sstack(
    dt: f32,
    nt: usize,
    xmin: f32,
    dx: f32,
    pmin: f32,
    dp: f32,
    traces: &[Vec<f32>],
    out_traces: &mut [Vec<f32>],
) {
    for ip in 0..out_traces.len() {
        let p = pmin + ip as f32 * dp;
        out_traces[ip][..nt].fill(0.0);
        for ix in 0..traces.len() {
            let x = xmin + ix as f32 * dx;
            let delay = p * x / dt;
            let (id, fit, lit) = if delay >= 0.0 {
                (delay as i32, delay as i32 + 1, nt as i32 - 1)
            } else {
                (delay as i32 - 1, 1, nt as i32 + delay as i32 - 1)
            };
            let dfrac = delay - id as f32;
            for it in fit..lit {
                let it = it as usize;
                let v = traces[ix][it] + dfrac * (traces[ix][it + 1] - traces[ix][it]);
                let k = (it as i32 - id) as usize;
                out_traces[ip][k] =
                    (out_traces[ip][k] as f64 + (dx as f64).abs() * v as f64) as f32;
            }
        }
    }
}

/// C: `void rho_filter(int npoints, int nt, float dt, float *rho)`
pub fn rho_filter(nt: usize, dt: f32, rho: &mut [f32]) {
    let npoints = rho.len();
    let nph = npoints / 2;
    assert!(npoints % 2 == 1 && nph <= nt && nt <= 720720);
    let ntfft = oxwave_pfa_fft::npfa(nt as i32) as usize;
    let nfh = ntfft / 2;
    let df = (1.0 / (2.0f32 * dt * nfh as f32) as f64) as f32;
    let mut cx = vec![Complex::default(); ntfft];
    for iff in 1..=nfh {
        let z = Complex {
            r: iff as f32 * df,
            i: 0.0,
        };
        cx[iff] = z;
        cx[ntfft - iff] = z;
    }
    oxwave_pfa_fft::pfacc(-1, &mut cx);
    rho[nph] = 1.0;
    for it in 0..nph {
        let v = cx[it + 1].r / cx[0].r;
        rho[nph - it - 1] = v;
        rho[nph + it + 1] = v;
    }
}

/// C: `void inv_tx_sstack(float dt, int nt, int nx, int npoints, float xmin, float dx, int np, float pmin, float dp, float **traces, float **out_traces)`
pub fn inv_tx_sstack(
    dt: f32,
    nt: usize,
    npoints: usize,
    xmin: f32,
    dx: f32,
    pmin: f32,
    dp: f32,
    traces: &[Vec<f32>],
    out_traces: &mut [Vec<f32>],
) {
    let mut rho = vec![0.0; npoints];
    rho_filter(nt, dt, &mut rho);
    let mut filtered = vec![vec![0.0; nt]; traces.len()];
    for ip in 0..traces.len() {
        oxwave_f_convolution::convolve_cwp(
            -(npoints as i32 / 2),
            &traces[ip][..nt],
            0,
            &rho,
            0,
            &mut filtered[ip],
        );
    }
    fwd_tx_sstack(dt, nt, -pmin, -dp, xmin, dx, &filtered, out_traces);
}

/// C: `void inverse_p_transform(int nx, int nt, float dt, float pmax, float pmin, float dp, float depthref, float f1, float f2, float interoff, float offref, int igopt, float dx, float fx, float **in_traces, float **out_traces)`
pub fn inverse_p_transform(
    nt: usize,
    dt: f32,
    mut pmax: f32,
    mut pmin: f32,
    mut dp: f32,
    depthref: f32,
    f1: f32,
    f2: f32,
    interoff: f32,
    offref: f32,
    igopt: i32,
    dx: f32,
    fx: f32,
    in_traces: &[Vec<f32>],
    out_traces: &mut [Vec<f32>],
) {
    let nx = out_traces.len();
    let ntfft = oxwave_pfa_fft::npfar(nt as i32) as usize;
    let fac = (1000.0 * gofx(igopt, offref, interoff, depthref) as f64) as f32;
    pmin /= fac;
    pmax /= fac;
    dp /= fac;
    let np = (1.0 + ((pmax - pmin) / dp) as f64) as usize;
    oxwave_core::fprintf!(
        oxwave_core::cio::stderr(),
        "computing inverse radon transform: pmin=%g pmax=%g dp=%g np=%d\n",
        pmin,
        pmax,
        dp,
        np as i32
    );
    let nw = 1 + ntfft / 2;
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    let dw = (2.0 * std::f64::consts::PI * df as f64) as f32;
    let fac = (1.0 / ntfft as f64) as f32;
    let mut spec = vec![vec![Complex::default(); nw]; np];
    for ip in 0..np {
        let mut rt = vec![0.0; ntfft];
        for it in 0..nt {
            rt[it] = in_traces[ip][it] * fac;
        }
        oxwave_pfa_fft::pfarc(1, &rt, &mut spec[ip]);
    }
    let mut output = vec![vec![Complex::default(); nw]; nx];
    let fac = (1.0 / np as f64) as f32;
    for iw in 0..nw {
        let wa = freqweight(iw as i32, df, f1, f2);
        if wa > 0.0 {
            let w = iw as f32 * dw;
            for ix in 0..nx {
                let g = gofx(igopt, fx + ix as f32 * dx, interoff, depthref);
                let (mut rs, mut is) = (0.0f32, 0.0f32);
                for ip in 0..np {
                    let p = pmin + ip as f32 * dp;
                    let angle = w * p * g;
                    let tr = (angle as f64).cos() as f32;
                    let ti = (angle as f64).sin() as f32;
                    let z = if wa < 1.0 {
                        crmul(spec[ip][iw], wa)
                    } else {
                        spec[ip][iw]
                    };
                    rs += tr * z.r - ti * z.i;
                    is += tr * z.i + ti * z.r;
                }
                output[ix][iw] = Complex {
                    r: fac * rs,
                    i: fac * is,
                };
            }
        }
    }
    for ix in 0..nx {
        let mut rt = vec![0.0; ntfft];
        oxwave_pfa_fft::pfacr(-1, &output[ix], &mut rt);
        out_traces[ix][..nt].copy_from_slice(&rt[..nt]);
    }
}

/// C: `void inv_FK_sstack(float dt, int nt, int nx, float xmin, float dx, int np, float pmin, float dp, float fmin, float **traces, float **out_traces)`
#[allow(non_snake_case)]
pub fn inv_FK_sstack(
    dt: f32,
    nt: usize,
    xmin: f32,
    dx: f32,
    pmin: f32,
    dp: f32,
    fmin: f32,
    traces: &[Vec<f32>],
    out_traces: &mut [Vec<f32>],
) {
    let np = traces.len();
    let nx = out_traces.len();
    let pmax0 = pmin + (np - 1) as f32 * dp;
    let fp = (pmin as f64 + 0.5 * (pmax0 - pmin - (np - 1) as f32 * dp) as f64) as f32;
    let pmax = fp + (np - 1) as f32 * dp;
    let ntfft = oxwave_pfa_fft::npfar((2 * nt) as i32) as usize;
    let nkfft = oxwave_pfa_fft::npfa((4 * np) as i32) as usize;
    let nw = ntfft / 2 + 1;
    let fftscale = (1.0 / ntfft as f64) as f32;
    let dw = (2.0 * std::f64::consts::PI / (ntfft as f32 * dt) as f64) as f32;
    let dk = (2.0 * std::f64::consts::PI / (nkfft as f32 * dx) as f64) as f32;
    let zero = Complex {
        r: 0.0 * fmin,
        i: 0.0 * xmin,
    };
    let mut spec = vec![vec![Complex::default(); nw]; np];
    for ip in 0..np {
        let mut rt = vec![0.0; ntfft];
        rt[..nt].copy_from_slice(&traces[ip][..nt]);
        oxwave_pfa_fft::pfarc(1, &rt, &mut spec[ip]);
    }
    let mut out = vec![vec![Complex::default(); nw]; nx];
    let mut w = 0.0f32;
    for iw in 0..nw {
        let trp: Vec<_> = spec
            .iter()
            .map(|row| crmul(row[iw], dp * fftscale))
            .collect();
        let nk = (pmax * w / dk + 1.0) as usize;
        let mut pk = vec![0.0; nk];
        let mut k = 0.0f32;
        for v in &mut pk {
            *v = if w == 0.0 { 0.0 } else { k / w };
            k += dk;
        }
        let mut hk = vec![zero; nk.max(nkfft)];
        oxwave_sinc_interp::ints8c(dp, fp, &trp, zero, zero, &pk, &mut hk[..nk]);
        for ik in 0..nk {
            let phase = pk[ik] * w * xmin;
            let c = (phase as f64).cos() as f32;
            let s = (phase as f64).sin() as f32;
            let temp = hk[ik].r * c - hk[ik].i * s;
            hk[ik].i = hk[ik].r * s + hk[ik].i * c;
            hk[ik].r = temp;
        }
        oxwave_pfa_fft::pfacc(1, &mut hk[..nkfft]);
        for ix in 0..nx {
            out[ix][iw] = hk[ix];
        }
        w += dw;
    }
    for ix in 0..nx {
        let mut rt = vec![0.0; ntfft];
        oxwave_pfa_fft::pfacr(-1, &out[ix], &mut rt);
        out_traces[ix][..nt].copy_from_slice(&rt[..nt]);
    }
}

/// C: `void fwd_FK_sstack(float dt, int nt, int nx, float xmin, float dx, int np, float pmin, float dp, float fmin, float **traces, float **out_traces)`
#[allow(non_snake_case)]
pub fn fwd_FK_sstack(
    dt: f32,
    nt: usize,
    xmin: f32,
    dx: f32,
    pmin: f32,
    dp: f32,
    fmin: f32,
    traces: &[Vec<f32>],
    out_traces: &mut [Vec<f32>],
) {
    let nx = traces.len();
    let np = out_traces.len();
    let pmax = pmin + (np - 1) as f32 * dp;
    let fp = (pmin as f64 + 0.5 * (pmax - pmin - (np - 1) as f32 * dp) as f64) as f32;
    let pmax = if dp < 0.0 {
        fp
    } else {
        fp + (np - 1) as f32 * dp
    };
    let xmax = if dx < 0.0 {
        xmin
    } else {
        xmin + (nx - 1) as f32 * dx
    };
    let ntpad = (pmax * xmax).abs() / dt;
    let ntfft = oxwave_pfa_fft::npfar((nt + ntpad as usize) as i32) as usize;
    let nxfft = oxwave_pfa_fft::npfa((2 * (nx + 8)) as i32) as usize;
    let nw = ntfft / 2 + 1;
    let fftscale = (1.0 / ntfft as f64) as f32;
    let dw = (2.0 * std::f64::consts::PI / (ntfft as f32 * dt) as f64) as f32;
    let dk = (2.0 * std::f64::consts::PI / (nxfft as f32 * dx) as f64) as f32;
    let fk = (-std::f64::consts::PI / dx as f64) as f32;
    let fka = fk - 8.0 * dk;
    let zero = Complex {
        r: 0.0 * fmin,
        i: 0.0 * xmin,
    };
    let mut spec = vec![vec![Complex::default(); nw]; nx];
    for ix in 0..nx {
        let mut rt = vec![0.0; ntfft];
        rt[..nt].copy_from_slice(&traces[ix][..nt]);
        oxwave_pfa_fft::pfarc(1, &rt, &mut spec[ix]);
    }
    let mut out = vec![vec![Complex::default(); nw]; np];
    let mut w = 0.0f32;
    for iw in 0..nw {
        let mut ka = vec![zero; nxfft + 16];
        for ix in 0..nx {
            ka[ix + 8] = Complex {
                r: spec[ix][iw].r * dx * fftscale,
                i: spec[ix][iw].i * dx * fftscale,
            };
            if ix % 2 == 1 {
                ka[ix + 8].r = -ka[ix + 8].r;
                ka[ix + 8].i = -ka[ix + 8].i;
            }
        }
        oxwave_pfa_fft::pfacc(-1, &mut ka[8..8 + nxfft]);
        for ik in 0..8 {
            ka[ik] = ka[nxfft + ik];
            ka[8 + nxfft + ik] = ka[8 + ik];
        }
        let shift = (0.5 * (nx - 1) as f64 * dx as f64) as f32;
        let mut k = fka;
        for z in &mut ka {
            let phase = k * shift;
            let c = (phase as f64).cos() as f32;
            let s = (phase as f64).sin() as f32;
            let t = z.r * c - z.i * s;
            z.i = z.r * s + z.i * c;
            z.r = t;
            k += dk;
        }
        let mut kp = vec![0.0; np];
        let mut p = fp;
        for v in &mut kp {
            *v = w * p;
            p += dp;
        }
        let mut hp = vec![zero; np];
        oxwave_sinc_interp::ints8c(dk, fka, &ka, zero, zero, &kp, &mut hp);
        let shift = (-xmin as f64 - 0.5 * (nx - 1) as f64 * dx as f64) as f32;
        for ip in 0..np {
            let phase = kp[ip] * shift;
            let c = (phase as f64).cos() as f32;
            let s = (phase as f64).sin() as f32;
            out[ip][iw] = Complex {
                r: hp[ip].r * c - hp[ip].i * s,
                i: hp[ip].r * s + hp[ip].i * c,
            };
        }
        w += dw;
    }
    for ip in 0..np {
        let mut rt = vec![0.0; ntfft];
        oxwave_pfa_fft::pfacr(-1, &out[ip], &mut rt);
        out_traces[ip][..nt].copy_from_slice(&rt[..nt]);
    }
}

/// C: `void ga_xinterpolate(float **in_traces, float **out_traces, int ninterp, int nt, int nx, float freq1, float freq2, int lagc, int lent, int lenx, int xopt, float dt, int iopt)`
pub fn ga_xinterpolate(
    input: &[Vec<f32>],
    output: &mut [Vec<f32>],
    ninterp: usize,
    nt: usize,
    freq1: f32,
    freq2: f32,
    lagc: i32,
    lent: i32,
    lenx: i32,
    mut xopt: i32,
    dt: f32,
    iopt: i32,
) {
    let nx = input.len();
    if nx < 2 || (iopt == 0 && ninterp == 0) {
        for ix in 0..nx {
            output[ix][..nt].copy_from_slice(&input[ix][..nt]);
        }
        return;
    }
    let lent = 1 + 2 * (lent / 2);
    let lenx = 1 + 2 * (lenx / 2);
    let lagc = (1.0 + lagc as f64 * 0.001 / dt as f64) as i32;
    let ntfft = oxwave_pfa_fft::npfar(nt as i32) as usize;
    let nxfft = oxwave_pfa_fft::npfar(nx as i32) as usize;
    let nw = ntfft / 2 + 1;
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    let dff = (2.0 * std::f64::consts::PI / ntfft as f64) as f32;
    let mut a = vec![vec![0.0; nt]; nx];
    let mut b = a.clone();
    for ix in 0..nx {
        let mut rt = vec![0.0; ntfft];
        let av: Vec<_> = input[ix][..nt]
            .iter()
            .map(|v| (*v as f64).abs() as f32)
            .collect();
        let mut bv = vec![0.0; nt];
        let mut av2 = bv.clone();
        runav(lagc, &av, &mut bv);
        runav(lagc, &bv, &mut av2);
        for it in 0..nt {
            rt[it] = input[ix][it] / (av2[it] + 1.0e-30f32);
        }
        let mut c = vec![Complex::default(); nw];
        oxwave_pfa_fft::pfarc(1, &rt, &mut c);
        let mut cc = c.clone();
        for it in 0..nw {
            let f = it as f32 * df;
            let mut fcl = f / freq1;
            fcl = fcl * fcl * fcl * fcl;
            let mut fch = f / freq2;
            fch = fch * fch * fch * fch;
            let f = (fcl as f64 / ((1.0 + fcl as f64) * (1.0 + fch as f64))) as f32;
            c[it] = crmul(c[it], f);
            cc[it] = cmul(
                c[it],
                Complex {
                    r: 0.0,
                    i: -(it as f32) * dff,
                },
            );
        }
        oxwave_pfa_fft::pfacr(-1, &c, &mut rt);
        a[ix].copy_from_slice(&rt[..nt]);
        oxwave_pfa_fft::pfacr(-1, &cc, &mut rt);
        b[ix].copy_from_slice(&rt[..nt]);
    }
    if iopt == 1 {
        for ix in 0..nx {
            output[ix][..nt].copy_from_slice(&a[ix]);
        }
        return;
    }
    if nx < 16 {
        xopt = 1;
    }
    let dxx = (2.0 * std::f64::consts::PI / (nxfft * nxfft) as f64) as f32;
    for it in 0..nt {
        let mut rt = vec![0.0; nxfft];
        for ix in 0..nx {
            rt[ix] = a[ix][it];
        }
        if xopt != 0 {
            for j in 0..nx - 1 {
                rt[j] = rt[j + 1] - rt[j];
            }
            rt[nx - 1] = rt[nx - 2];
        } else {
            let mut c = vec![Complex::default(); nxfft / 2 + 1];
            oxwave_pfa_fft::pfarc(1, &rt, &mut c);
            for j in 0..c.len() {
                c[j] = cmul(
                    c[j],
                    Complex {
                        r: 0.0,
                        i: -(j as f32) * dxx,
                    },
                );
            }
            oxwave_pfa_fft::pfacr(-1, &c, &mut rt);
        }
        for ix in 0..nx {
            a[ix][it] = rt[ix];
        }
    }
    for ix in 0..nx {
        let aa: Vec<_> = (0..nt).map(|it| a[ix][it] * b[ix][it]).collect();
        let bb: Vec<_> = (0..nt).map(|it| b[ix][it] * b[ix][it]).collect();
        let mut tmp = vec![0.0; nt];
        runav(lent, &aa, &mut tmp);
        runav(lent, &tmp, &mut a[ix]);
        runav(lent, &bb, &mut tmp);
        runav(lent, &tmp, &mut b[ix]);
    }
    if lenx > 1 {
        for it in 0..nt {
            let aa: Vec<_> = (0..nx).map(|ix| a[ix][it]).collect();
            let bb: Vec<_> = (0..nx).map(|ix| b[ix][it]).collect();
            let mut tmp = vec![0.0; nx];
            let mut res = tmp.clone();
            runav(lenx, &aa, &mut tmp);
            runav(lenx, &tmp, &mut res);
            for ix in 0..nx {
                a[ix][it] = res[ix];
            }
            runav(lenx, &bb, &mut tmp);
            runav(lenx, &tmp, &mut res);
            for ix in 0..nx {
                b[ix][it] = res[ix];
            }
        }
    }
    let mut p = vec![0.0; nt];
    for ix in 1..nx {
        for it in 0..nt {
            p[it] = -a[ix - 1][it] / (b[ix - 1][it] + 1.0e-30f32);
        }
        if iopt == 2 {
            output[ix - 1][..nt].copy_from_slice(&p);
            continue;
        }
        output[(ix - 1) * (ninterp + 1)][..nt].copy_from_slice(&input[ix - 1][..nt]);
        for k in 0..ninterp {
            let wa = ((1.0 + k as f64) / (1 + ninterp) as f64) as f32;
            let wb = (1.0 - wa as f64) as f32;
            let time: Vec<_> = (0..nt).map(|it| it as f32 - p[it] * wa).collect();
            let mut aa = vec![0.0; nt];
            let mut bb = aa.clone();
            oxwave_sinc_interp::ints8r(1.0, 0.0, &input[ix - 1][..nt], 0.0, 0.0, &time, &mut aa);
            let time: Vec<_> = (0..nt).map(|it| it as f32 + p[it] * wb).collect();
            oxwave_sinc_interp::ints8r(1.0, 0.0, &input[ix][..nt], 0.0, 0.0, &time, &mut bb);
            for it in 0..nt {
                output[k + 1 + (ix - 1) * (ninterp + 1)][it] = wb * aa[it] + wa * bb[it];
            }
        }
    }
    if iopt == 2 {
        output[nx - 1][..nt].copy_from_slice(&p);
    } else {
        output[(nx - 1) * (ninterp + 1)][..nt].copy_from_slice(&input[nx - 1][..nt]);
    }
}

/// C: `void forward_p_transform(int nx, int nt, float dt, float pmax, float pmin, float dp, float depthref, float f1, float f2, float freq1, float freq2, int lagc, int lent, int lenx, int xopt, int ninterp, int nwin, float prewhite, float interoff, float offref, int igopt, float dx, float fx, float pmula, float pmulb, float **in_traces, float **out_traces)`
pub fn forward_p_transform(
    nt: usize,
    dt: f32,
    mut pmax: f32,
    mut pmin: f32,
    mut dp: f32,
    depthref: f32,
    f1: f32,
    f2: f32,
    freq1: f32,
    freq2: f32,
    lagc: i32,
    lent: i32,
    lenx: i32,
    xopt: i32,
    ninterp: usize,
    _nwin: i32,
    prewhite: f32,
    interoff: f32,
    offref: f32,
    igopt: i32,
    dx: f32,
    fx: f32,
    mut pmula: f32,
    mut pmulb: f32,
    input: &[Vec<f32>],
    output: &mut [Vec<f32>],
) {
    let nx = input.len();
    let ntfft = oxwave_pfa_fft::npfar(nt as i32) as usize;
    let fac = (1000.0 * gofx(igopt, offref, interoff, depthref) as f64) as f32;
    pmin /= fac;
    pmax /= fac;
    dp /= fac;
    pmula /= fac;
    pmulb /= fac;
    let ipa = (((pmula - pmin) / dp) as i32).max(0);
    let ipb = ((pmulb - pmin) / dp) as i32;
    let np = (1.0 + ((pmax - pmin) / dp) as f64) as usize;
    oxwave_core::fprintf!(
        oxwave_core::cio::stderr(),
        "computing forward radon transform: pmin=%g pmax=%g dp=%g np=%d\n",
        pmin,
        pmax,
        dp,
        np as i32
    );
    let nxi = (1 + ninterp) * (nx - 1) + 1;
    let mut interp = vec![vec![0.0; nt]; nxi];
    ga_xinterpolate(
        input,
        &mut interp,
        ninterp,
        nt,
        freq1,
        freq2,
        lagc,
        lent,
        lenx,
        xopt,
        dt,
        0,
    );
    let offset: Vec<_> = (0..nx).map(|ix| fx + ix as f32 * dx).collect();
    let xin: Vec<_> = (0..nx).map(|ix| ix as f32).collect();
    let d = (1.0 / (1 + ninterp) as f64) as f32;
    let rtr: Vec<_> = (0..nxi).map(|ix| ix as f32 * d).collect();
    let mut g = vec![0.0; nxi];
    oxwave_linear_interp::intlin(&xin, &offset, offset[0], offset[nx - 1], &rtr, &mut g);
    for x in &mut g {
        *x = gofx(igopt, *x, interoff, depthref);
    }
    let nw = 1 + ntfft / 2;
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    let dw = (2.0 * std::f64::consts::PI * df as f64) as f32;
    let fac = (1.0 / ntfft as f64) as f32;
    let mut spec = vec![vec![Complex::default(); nw]; nxi];
    for ix in 0..nxi {
        let mut rt = vec![0.0; ntfft];
        for it in 0..nt {
            rt[it] = interp[ix][it] * fac;
        }
        oxwave_pfa_fft::pfarc(1, &rt, &mut spec[ix]);
    }
    let mut out = vec![vec![Complex::default(); nw]; np];
    let mut solution = vec![Complex::default(); np];
    let mut r = solution.clone();
    let mut rhs = r.clone();
    let mut w1 = r.clone();
    let mut w2 = r.clone();
    let mut w3 = r.clone();
    let mut w4 = r.clone();
    for iw in 0..nw {
        let wa = freqweight(iw as i32, df, f1, f2);
        if wa > 0.0 {
            let w = iw as f32 * dw;
            let data: Vec<_> = spec
                .iter()
                .map(|v| if wa < 1.0 { crmul(v[iw], wa) } else { v[iw] })
                .collect();
            compute_rhs(w, &g, &data, pmin, dp, &mut rhs);
            compute_r(w, &g, dp, &mut r);
            r[0].r = (r[0].r as f64 * (1.0 + prewhite as f64)) as f32;
            let mut rsum = 0.0f32;
            for ip in 1..np {
                rsum =
                    (rsum as f64 + ((r[ip].r * r[ip].r + r[ip].i * r[ip].i) as f64).sqrt()) as f32;
            }
            rsum /= r[0].r;
            if rsum as f64 > 1.0 + (np / 5) as f64 {
                ctoephcg(
                    (np / 7) as i32,
                    &r,
                    &mut solution,
                    &rhs,
                    &mut w1,
                    &mut w2,
                    &mut w3,
                    &mut w4,
                );
            } else {
                ctoep(&r, &mut solution, &rhs, &mut w1, &mut w2);
            }
        } else {
            solution.fill(Complex::default());
        }
        for ip in 0..np {
            out[ip][iw] = solution[ip];
        }
    }
    for ip in 0..np {
        let mut rt = vec![0.0; ntfft];
        oxwave_pfa_fft::pfacr(-1, &out[ip], &mut rt);
        taupmute(ip as i32, ipa, ipb, 7, &mut rt);
        output[ip][..nt].copy_from_slice(&rt[..nt]);
    }
}
