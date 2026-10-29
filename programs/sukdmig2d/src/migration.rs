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

//! Numerical helpers from sukdmig2d.c, Revision 1.26.
//! Credits: Zhenyue Liu, Colorado School of Mines, 03/01/95.
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::*;

// Residual traveltime calculation based on reference time.
pub fn resit(
    nx: usize,
    fx: f32,
    dx: f32,
    nz: usize,
    nr: usize,
    dr: f32,
    tb: &[Vec<f32>],
    t: &mut [f32],
    x0: f32,
) {
    for ix in 0..nx {
        let xi = fx + ix as f32 * dx - x0;
        let ar = xi.abs() / dr;
        let mut jr = ar as usize;
        let sr = ar - jr as f32;
        let sr0 = (1.0 - sr as f64) as f32;
        if jr > nr - 2 {
            jr = nr - 2;
        }
        for iz in 0..nz {
            t[ix * nz + iz] -= sr0 * tb[jr][iz] + sr * tb[jr + 1][iz];
        }
    }
}

// Sum of two tables; separate destination permits the C in-place uses below.
pub fn sum2(a1: f32, a2: f32, t1: &[f32], t2: &[f32], t: &mut [f32]) {
    for i in 0..t.len() {
        t[i] = a1 * t1[i] + a2 * t2[i];
    }
}

// Compute reference traveltime and slowness.
#[allow(clippy::too_many_arguments)]
pub fn timeb(
    nr: usize,
    nz: usize,
    dr: f32,
    dz: f32,
    fz: f32,
    a: f32,
    v0: f32,
    t: &mut [Vec<f32>],
    p: &mut [Vec<f32>],
    cs0: &mut [Vec<f32>],
    ang: &mut [Vec<f32>],
) {
    if a == 0.0 {
        let mut r = 0.0f32;
        for ir in 0..nr {
            let mut z = fz;
            for iz in 0..nz {
                let mut rou = ((r * r + z * z) as f64).sqrt() as f32;
                if rou < dz {
                    rou = dz;
                }
                t[ir][iz] = rou / v0;
                p[ir][iz] = r / (rou * v0);
                cs0[ir][iz] = z / rou;
                ang[ir][iz] = ((r / rou) as f64).asin() as f32;
                z += dz;
            }
            r += dr;
        }
    } else {
        let oa = (1.0 / a as f64) as f32;
        let zc = v0 * oa;
        let mut r = 0.0f32;
        for ir in 0..nr {
            let mut z = fz + zc;
            for iz in 0..nz {
                let v = v0 + a * (z - zc);
                if ir == 0 {
                    t[ir][iz] = (((v / v0) as f64).ln() * oa as f64) as f32;
                    p[ir][iz] = 0.0;
                    ang[ir][iz] = 0.0;
                    cs0[ir][iz] = 1.0;
                } else {
                    let rc = ((r * r + z * z - zc * zc) as f64 / (2.0 * r as f64)) as f32;
                    let rou = ((zc * zc + rc * rc) as f64).sqrt() as f32;
                    t[ir][iz] =
                        (((v * (rou + rc) / (v0 * (rou + rc - r))) as f64).ln() * oa as f64) as f32;
                    p[ir][iz] = (((rou * rou - rc * rc) as f64).sqrt() / (rou * v0) as f64) as f32;
                    let mut temp = v * p[ir][iz];
                    if temp > 1.0 {
                        temp = 1.0;
                    }
                    ang[ir][iz] = (temp as f64).asin() as f32;
                    cs0[ir][iz] = rc / rou;
                }
                z += dz;
            }
            r += dr;
        }
    }
}

// The static state in filt(), owned by main and shared across input traces.
#[derive(Default)]
pub struct Filter {
    nfft: usize,
    itaper: i32,
    nw: i32,
    nwf: i32,
    taper: Vec<f32>,
    amp: Vec<f32>,
    ampi: Vec<f32>,
}
impl Filter {
    fn filt(
        &mut self,
        trace: &mut [f32],
        nt: usize,
        dt: f32,
        mut fmax: f32,
        ls: i32,
        m: usize,
        trf: &mut [f32],
    ) {
        fmax = (fmax as f64 * (2.0 * PI)) as f32;
        let ftaper = (0.1 * fmax as f64) as f32;
        let const2 = 2.0f64.sqrt() as f32;
        if self.nfft == 0 {
            self.nfft = npfaro((nt + m) as i32, (2 * (nt + m)) as i32) as usize;
            if self.nfft >= SU_NFLTS || self.nfft >= 720720 {
                err!("Padded nt=%d -- too big", self.nfft);
            }
            self.nw = (self.nfft / 2 + 1) as i32;
            let dw = (2.0 * PI / (self.nfft as f32 * dt) as f64) as f32;
            self.itaper = (0.5 + (ftaper / dw) as f64) as i32;
            self.taper = vec![0.0; (2 * self.itaper + 1) as usize];
            for iw in -self.itaper..=self.itaper {
                let temp = (iw as f32 as f64 / (1.0 + self.itaper as f64)) as f32;
                self.taper[(iw + self.itaper) as usize] =
                    (1.0 - temp) * (1.0 - temp) * (temp + 2.0) / 4.0;
            }
            self.nwf = (0.5 + (fmax / dw) as f64) as i32;
            if self.nwf > self.nw - self.itaper - 1 {
                self.nwf = self.nw - self.itaper - 1;
            }
            self.amp = vec![0.0; (self.nwf + self.itaper + 1) as usize];
            self.ampi = vec![0.0; self.amp.len()];
            for iw in 1..=self.nwf + self.itaper {
                self.amp[iw as usize] =
                    (((dw * iw as f32) as f64).sqrt() / self.nfft as f64) as f32;
                self.ampi[iw as usize] =
                    (0.5 / (1.0 - ((iw as f32 * dw * dt) as f64).cos())) as f32;
            }
        }
        let mut rt = vec![0.0; self.nfft];
        let mut ct = vec![Complex::default(); self.nw as usize];
        rt[..nt].copy_from_slice(&trace[..nt]);
        pfarc(1, &rt, &mut ct);
        for iw in self.nwf - self.itaper..=self.nwf + self.itaper {
            let itemp = iw - (self.nwf - self.itaper);
            ct[iw as usize].r = self.taper[itemp as usize] * ct[iw as usize].r;
            ct[iw as usize].i = self.taper[itemp as usize] * ct[iw as usize].i;
        }
        for iw in self.nwf + self.itaper + 1..self.nw {
            ct[iw as usize].r = 0.0;
            ct[iw as usize].i = 0.0;
        }
        for iw in 0..=self.nwf + self.itaper {
            let iw = iw as usize;
            if ls == 0 {
                let temp = (ct[iw].r - ct[iw].i) * self.amp[iw] * const2;
                ct[iw].i = (ct[iw].r + ct[iw].i) * self.amp[iw] * const2;
                ct[iw].r = temp;
            } else {
                ct[iw].i = ct[iw].i * self.amp[iw];
                ct[iw].r = ct[iw].r * self.amp[iw];
            }
        }
        pfacr(-1, &ct, &mut rt);
        for it in 0..nt {
            trace[it] = rt[it];
        }
        for iw in 0..=self.nwf + self.itaper {
            let iw = iw as usize;
            ct[iw].i = ct[iw].i * self.ampi[iw];
            ct[iw].r = ct[iw].r * self.ampi[iw];
        }
        pfacr(-1, &ct, &mut rt);
        for it in 0..m {
            trf[it] = rt[self.nfft - m + it];
        }
        for it in 0..nt + m {
            trf[it + m] = rt[it];
        }
    }
}

// Migrate one trace. Keep the arithmetic and interpolation order of the C.
#[allow(clippy::too_many_arguments)]
pub fn mig2d(
    trace: &mut [f32],
    nt: usize,
    ft: f32,
    dt: f32,
    sx: f32,
    gx: f32,
    mig: &mut [Vec<f32>],
    aperx: f32,
    nx: i32,
    fx: f32,
    dx: f32,
    nz: i32,
    fz: f32,
    dz: f32,
    ls: i32,
    mtmax: i32,
    dxm: f32,
    fmax: f32,
    angmax: f32,
    tb: &[Vec<f32>],
    pb: &[Vec<f32>],
    cs0b: &[Vec<f32>],
    angb: &[Vec<f32>],
    nr: i32,
    tsum: &[f32],
    nzt: i32,
    fzt: f32,
    dzt: f32,
    nxt: i32,
    fxt: f32,
    dxt: f32,
    npv: i32,
    cssum: &[f32],
    tvsum: &[f32],
    mig1: &mut [Vec<f32>],
    filter: &mut Filter,
) {
    let mut tmt = vec![vec![0.0f32; nzt as usize]; nxt as usize];
    let mut ampt = tmt.clone();
    let mut ampti = tmt.clone();
    let mut ampt1 = if npv != 0 { tmt.clone() } else { Vec::new() };
    let mut tm = vec![0.0f32; nzt as usize];
    let mut tzt = tm.clone();
    let mut amp = tm.clone();
    let mut ampi = tm.clone();
    let mut amp1 = if npv != 0 { tm.clone() } else { Vec::new() };
    let mut zpt = vec![0.0f32; nxt as usize];
    let mut trf = vec![0.0f32; nt + 2 * mtmax as usize];
    let z0 = (fz - fzt) / dzt;
    let rdz = dz / dzt;
    let pmin = (1.0 / (2.0 * dxm as f64 * fmax as f64)) as f32;
    let odt = (1.0 / dt as f64) as f32;
    filter.filt(trace, nt, dt, fmax, ls, mtmax as usize, &mut trf);
    let xm = (0.5 * (sx + gx) as f64) as f32;
    let rxz = if angmax == 90.0 {
        0.0
    } else {
        (1.0 / (angmax as f64 * PI / 180.0).tan()) as f32
    };
    let mut nxtf = ((xm - aperx - fxt) / dxt) as i32;
    if nxtf < 0 {
        nxtf = 0;
    }
    let mut nxte = ((xm + aperx - fxt) / dxt + 1.0) as i32;
    if nxte >= nxt {
        nxte = nxt - 1;
    }
    // Compute amplitudes and filter length.
    for ix in nxtf..=nxte {
        let x = fxt + ix as f32 * dxt;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = (((dis - dxt) * rxz - fzt) / dzt - 1.0) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let ar = if sx >= x {
            (sx - x) / dx
        } else {
            (x - sx) / dx
        };
        let mut jrs = ar as i32;
        if jrs > nr - 2 {
            jrs = nr - 2;
        }
        let srs = ar - jrs as f32;
        let srs0 = (1.0 - srs as f64) as f32;
        let ar = if gx >= x {
            (gx - x) / dx
        } else {
            (x - gx) / dx
        };
        let mut jrg = ar as i32;
        if jrg > nr - 2 {
            jrg = nr - 2;
        }
        let srg = ar - jrg as f32;
        let srg0 = (1.0 - srg as f64) as f32;
        let sigp = if (sx - x) * (gx - x) > 0.0 {
            1.0f32
        } else {
            -1.0f32
        };
        let (ix, jrs, jrg) = (ix as usize, jrs as usize, jrg as usize);
        zpt[ix] = fzt + (nzt - 1) as f32 * dzt;
        for iz in izt0 as usize..nzt as usize {
            let angs = srs0 * angb[jrs][iz] + srs * angb[jrs + 1][iz];
            let angg = srg0 * angb[jrg][iz] + srg * angb[jrg + 1][iz];
            let cs0s = srs0 * cs0b[jrs][iz] + srs * cs0b[jrs + 1][iz];
            let cs0g = srg0 * cs0b[jrg][iz] + srg * cs0b[jrg + 1][iz];
            let mut ampd = ((cs0s + cs0g) as f64 * (0.5 * (angs - angg) as f64).cos()) as f32;
            if ampd < 0.0 {
                ampd = -ampd;
            }
            ampt[ix][iz] = ampd;
            let mut pd = srs0 * pb[jrs][iz]
                + srs * pb[jrs + 1][iz]
                + sigp * (srg0 * pb[jrg][iz] + srg * pb[jrg + 1][iz]);
            if pd < 0.0 {
                pd = -pd;
            }
            let mut temp = pd * dxm * odt;
            if temp < 1.0 {
                temp = 1.0;
            }
            if temp > mtmax as f32 {
                temp = mtmax as f32;
            }
            ampti[ix][iz] = ampd / (temp * temp);
            tmt[ix][iz] = temp;
            if pd < pmin && zpt[ix] as f64 > fzt as f64 + (nzt as f64 - 1.1) * dzt as f64 {
                zpt[ix] = fzt + iz as f32 * dzt;
            }
            if npv != 0 {
                ampt1[ix][iz] = if cssum[ix * nzt as usize + iz] < 1.0 {
                    0.0
                } else {
                    tvsum[ix * nzt as usize + iz] / cssum[ix * nzt as usize + iz]
                };
            }
        }
    }
    let mut nxf = (((xm - aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    let mut nxe = (((xm + aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxe >= nx {
        nxe = nx - 1;
    }
    // Interpolate amplitudes and filter length along lateral.
    for ix in nxf..=nxe {
        let x = fx + ix as f32 * dx;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = ((dis * rxz - fzt) / dzt) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let mut iz0 = ((dis * rxz - fz) / dz) as i32;
        if iz0 < 0 {
            iz0 = 0;
        }
        if iz0 as f32 >= nz as f32 {
            iz0 = nz - 1;
        }
        let mut ax = (x - fxt) / dxt;
        let mut jx = ax as i32;
        ax = ax - jx as f32;
        if ax as f64 <= 0.01 {
            ax = 0.0;
        }
        if ax as f64 >= 0.99 {
            ax = 1.0;
        }
        let ax0 = (1.0 - ax as f64) as f32;
        if jx > nxte - 1 {
            jx = nxte - 1;
        }
        if jx < nxtf {
            jx = nxtf;
        }
        let ar = if sx >= x {
            (sx - x) / dx
        } else {
            (x - sx) / dx
        };
        let mut jrs = ar as i32;
        if jrs > nr - 2 {
            jrs = nr - 2;
        }
        let srs = ar - jrs as f32;
        let srs0 = (1.0 - srs as f64) as f32;
        let ar = if gx >= x {
            (gx - x) / dx
        } else {
            (x - gx) / dx
        };
        let mut jrg = ar as i32;
        if jrg > nr - 2 {
            jrg = nr - 2;
        }
        let srg = ar - jrg as f32;
        let srg0 = (1.0 - srg as f64) as f32;
        let (ix, jx, jrs, jrg) = (ix as usize, jx as usize, jrs as usize, jrg as usize);
        for iz in izt0 as usize..nzt as usize {
            tzt[iz] = ax0 * tsum[jx * nzt as usize + iz]
                + ax * tsum[(jx + 1) * nzt as usize + iz]
                + srs0 * tb[jrs][iz]
                + srs * tb[jrs + 1][iz]
                + srg0 * tb[jrg][iz]
                + srg * tb[jrg + 1][iz];
            amp[iz] = ax0 * ampt[jx][iz] + ax * ampt[jx + 1][iz];
            ampi[iz] = ax0 * ampti[jx][iz] + ax * ampti[jx + 1][iz];
            tm[iz] = ax0 * tmt[jx][iz] + ax * tmt[jx + 1][iz];
            if npv != 0 {
                amp1[iz] = ax0 * ampt1[jx][iz] + ax * ampt1[jx + 1][iz];
            }
        }
        let mut nzp = (((ax0 * zpt[jx] + ax * zpt[jx + 1] - fz) / dz) as f64 + 1.5) as i32;
        if nzp < iz0 {
            nzp = iz0;
        }
        if nzp > nz {
            nzp = nz;
        }
        // Interpolate along depth if operator aliasing.
        for iz in iz0..nzp {
            let az = z0 + iz as f32 * rdz;
            let mut jz = az as i32;
            if jz >= nzt - 1 {
                jz = nzt - 2;
            }
            let sz = az - jz as f32;
            let sz0 = (1.0 - sz as f64) as f32;
            let jz = jz as usize;
            let td = sz0 * tzt[jz] + sz * tzt[jz + 1];
            let at = (td - ft) * odt + mtmax as f32;
            let jt = at as i32;
            if jt > mtmax && jt < nt as i32 + mtmax - 1 {
                let ampd = sz0 * ampi[jz] + sz * ampi[jz + 1];
                let mt = (0.5 + (sz0 * tm[jz]) as f64 + (sz * tm[jz + 1]) as f64) as i32;
                let res = at - jt as f32;
                let res0 = (1.0 - res as f64) as f32;
                let temp = ((res0 as f64
                    * (-trf[(jt - mt) as usize] as f64 + 2.0 * trf[jt as usize] as f64
                        - trf[(jt + mt) as usize] as f64)
                    + res as f64
                        * (-trf[(jt - mt + 1) as usize] as f64
                            + 2.0 * trf[(jt + 1) as usize] as f64
                            - trf[(jt + mt + 1) as usize] as f64))
                    * ampd as f64) as f32;
                mig[ix][iz as usize] += temp;
                if npv != 0 {
                    mig1[ix][iz as usize] += temp * (sz0 * amp1[jz] + sz * amp1[jz + 1]);
                }
            }
        }
        // Interpolate along depth if not operator aliasing.
        for iz in nzp..nz {
            let az = z0 + iz as f32 * rdz;
            let mut jz = az as i32;
            if jz >= nzt - 1 {
                jz = nzt - 2;
            }
            let sz = az - jz as f32;
            let sz0 = (1.0 - sz as f64) as f32;
            let jz = jz as usize;
            let td = sz0 * tzt[jz] + sz * tzt[jz + 1];
            let at = (td - ft) * odt;
            let jt = at as i32;
            if jt > 0 && jt < nt as i32 - 1 {
                let ampd = sz0 * amp[jz] + sz * amp[jz + 1];
                let res = at - jt as f32;
                let res0 = (1.0 - res as f64) as f32;
                let temp = (res0 * trace[jt as usize] + res * trace[jt as usize + 1]) * ampd;
                mig[ix][iz as usize] += temp;
                if npv != 0 {
                    mig1[ix][iz as usize] += temp * (sz0 * amp1[jz] + sz * amp1[jz + 1]);
                }
            }
        }
    }
}
