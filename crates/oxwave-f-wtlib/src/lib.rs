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

//! Wavelet transforms from par/lib/wtlib.c (Revision 1.4).
//! CWP: Zhaobo Meng, Colorado School of Mines, Sept 1996.
//! Modified: Carlos E. Theodoro, Colorado School of Mines, Jun 97.
//! The modifications allow different resolution levels per dimension and
//! reconstruction of only the lowest resolution component.
#![allow(non_snake_case)]

/// Wavelet filter and generated scaling and mother wavelets. Coefficients are one-based.
#[derive(Clone, Debug, Default)]
pub struct WtFilter {
    pub order: i32,
    pub ioff: i32,
    pub joff: i32,
    pub cc: Vec<f32>,
    pub cr: Vec<f32>,
    pub phi: Vec<f32>,
    pub psi: Vec<f32>,
}
/// Wavelet dimensions and resolution levels; dimension arrays are one-based.
#[derive(Clone, Debug, Default)]
pub struct WtSizes {
    pub MaxLevel: i32,
    pub NDim: i32,
    pub order: i32,
    pub NPointsn: Vec<i32>,
    pub Mraleveln: Vec<i32>,
}
/// C enum ToCorD: inverse (0) or forward (1).
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ToCorD {
    ToC = 0,
    ToD = 1,
}
const C4: [f32; 5] = [
    0.0,
    0.4829629131445341,
    0.8365163037378079,
    0.2241438680420314,
    -0.1294095225512604,
];
const C6: [f32; 7] = [
    0.0,
    0.332670552950,
    0.806891509311,
    0.459877502118,
    -0.135011020010,
    -0.085441273882,
    0.035226291882,
];
const C8: [f32; 9] = [
    0.0,
    0.230377813309,
    0.714846570553,
    0.630880767930,
    -0.027983769417,
    -0.187034811719,
    0.030841381836,
    0.032883011667,
    -0.010597401785,
];
const C10: [f32; 11] = [
    0.0,
    0.160102397974,
    0.603829269797,
    0.724308528438,
    0.138428145901,
    -0.242294887066,
    -0.032244869585,
    0.077571493840,
    -0.006241490213,
    -0.012580751999,
    0.003335725285,
];
const C12: [f32; 13] = [
    0.0,
    0.111540743350,
    0.494623890398,
    0.751133908021,
    0.315250351709,
    -0.226264693965,
    -0.129766867567,
    0.097501605587,
    0.027522865530,
    -0.031582039318,
    0.000553842201,
    0.004777257511,
    -0.001077301085,
];
const C14: [f32; 15] = [
    0.0,
    0.077852054085,
    0.396539319482,
    0.729132090846,
    0.469782287405,
    -0.143906003929,
    -0.224036184994,
    0.071309219276,
    0.080612609151,
    -0.038029936935,
    -0.016574541631,
    0.012550998556,
    0.000429577973,
    -0.001801640704,
    0.000353713800,
];
const C16: [f32; 17] = [
    0.0,
    0.054415842243,
    0.312871590914,
    0.675930736297,
    0.585354683654,
    -0.015829105056,
    -0.284015542962,
    0.000472484574,
    0.128747426620,
    -0.017369301002,
    -0.044088253931,
    0.013981027917,
    0.008746094047,
    -0.004870352993,
    -0.000391740373,
    0.000675449406,
    -0.000117476784,
];
const C18: [f32; 19] = [
    0.0,
    0.038077947364,
    0.243834674613,
    0.604823123690,
    0.657288078051,
    0.133197385825,
    -0.293273783279,
    -0.096840783223,
    0.148540749338,
    0.030725681479,
    -0.067632829061,
    0.000250947115,
    0.022361662124,
    -0.004723204758,
    -0.004281503682,
    0.001847646883,
    0.000230385764,
    -0.000251963189,
    0.000039347320,
];
const C20: [f32; 21] = [
    0.0,
    0.026670057901,
    0.188176800078,
    0.527201188932,
    0.688459039454,
    0.281172343661,
    -0.249846424327,
    -0.195946274377,
    0.127369340336,
    0.093057364604,
    -0.071394147166,
    -0.029457536822,
    0.033212674059,
    0.003606553567,
    -0.010733175483,
    0.001395351747,
    0.001992405295,
    -0.000685856695,
    -0.000116466855,
    0.000093588670,
    -0.000013264203,
];

/// C: `void wto1dset(WtFilter *wfilt, WtSizes *wtsizes);`
/// Set up Daubechies filters of even order 4 through 20, then generate
/// `phi` and `psi` at `MaxLevel` (which must be nonnegative).
/// An unsupported order prints the C diagnostic and exits with status zero.
pub fn wto1dset(f: &mut WtFilter, s: &WtSizes) {
    f.order = s.order;
    f.cc = match s.order {
        4 => &C4[..],
        6 => &C6,
        8 => &C8,
        10 => &C10,
        12 => &C12,
        14 => &C14,
        16 => &C16,
        18 => &C18,
        20 => &C20,
        _ => {
            oxwave_core::fprintf!(
                oxwave_core::cio::stderr(),
                "unimplemented order in wto1dset"
            );
            oxwave_core::exit(0);
        }
    }
    .to_vec();
    f.cr = vec![0.0; f.cc.len()];
    let mut sig = -1.0f32;
    for k in 1..=s.order as usize {
        f.cr[s.order as usize + 1 - k] = sig * f.cc[k];
        sig = -sig;
    }
    f.ioff = -(s.order >> 1);
    f.joff = f.ioff;
    wt_cascade(f, s);
}
/// C: `void wt_cascade(WtFilter *wfilt, WtSizes *wtsizes);`
/// Generate zero-based scaling and mother wavelet arrays of length
/// `wtsizes->order << wtsizes->MaxLevel` using the one-based `cc` filter.
pub fn wt_cascade(f: &mut WtFilter, s: &WtSizes) {
    let interval = 1i32 << s.MaxLevel;
    let nphi = (s.order * interval) as usize;
    let minusn = interval - (nphi as i32) / 2;
    f.phi = vec![0.0; nphi];
    f.psi = vec![0.0; nphi];
    let mut plus = vec![0.0; nphi];
    for k in 0..s.order as usize {
        f.phi[k] = f.cc[k + 1];
    }
    let mut coe = 1.0f32;
    for _ in 0..s.MaxLevel {
        coe = (coe as f64 * 2.0f64.sqrt()) as f32;
    }
    for _ in 1..s.MaxLevel {
        for n in 0..nphi {
            plus[n] = 0.0;
            for k in 0..nphi {
                let ix = n as i32 - 2 * k as i32;
                if ix >= 0 && ix < s.order {
                    plus[n] += f.cc[ix as usize + 1] * f.phi[k];
                }
            }
        }
        f.phi.copy_from_slice(&plus);
    }
    let mut sum = 0.0f32;
    for k in 0..nphi {
        f.phi[k] *= coe;
        sum += f.phi[k];
    }
    sum /= interval as f32;
    let _ = sum; // The C routine computes this normalization check, then discards it.
    coe = 2.0f64.sqrt() as f32;
    for k in 0..nphi {
        let mut sig = -1i32;
        for n in 0..s.order {
            let ix = 2 * (k as i32 + minusn) + (n - 1) * interval;
            if ix >= 0 && ix < nphi as i32 {
                f.psi[k] += (sig as f32) * f.cc[n as usize + 1] * coe * f.phi[ix as usize];
            }
            sig = -sig;
        }
    }
}
/// C: `void wto1d(float *cd, int npoints, enum ToCorD tocord, WtFilter *wfilt);`
/// `cd` retains C's unused element zero and must contain at least `npoints + 1` elements.
/// The original `npoints - 1` bit mask is retained for every nonnegative length,
/// including non-powers of two. Custom filter offsets use the C build's wrapping
/// 32-bit index arithmetic before that mask is applied.
pub fn wto1d(cd: &mut [f32], npoints: i32, tocord: ToCorD, f: &WtFilter) {
    let mut temp = vec![0.0f32; npoints as usize + 1];
    let nmod = f.order.wrapping_mul(npoints);
    let n1 = npoints - 1;
    let nhalf = npoints >> 1;
    let mut icd = 1i32;
    let mut i = 1i32;
    while i <= npoints {
        // Custom offsets can overflow before the periodic index mask. Match
        // the original C build's 32-bit arithmetic rather than panicking.
        let ni = i.wrapping_add(nmod).wrapping_add(f.ioff);
        let nj = i.wrapping_add(nmod).wrapping_add(f.joff);
        if tocord == ToCorD::ToD {
            for k in 1..=f.order {
                let jf = n1 & ni.wrapping_add(k);
                let jr = n1 & nj.wrapping_add(k);
                temp[icd as usize] += f.cc[k as usize] * cd[(jf + 1) as usize];
                temp[(icd + nhalf) as usize] += f.cr[k as usize] * cd[(jr + 1) as usize];
            }
        } else {
            let ai = cd[icd as usize];
            let ai1 = cd[(icd + nhalf) as usize];
            for k in 1..=f.order {
                let jf = (n1 & ni.wrapping_add(k)) + 1;
                let jr = (n1 & nj.wrapping_add(k)) + 1;
                temp[jf as usize] += f.cc[k as usize] * ai;
                temp[jr as usize] += f.cr[k as usize] * ai1;
            }
        }
        i += 2;
        icd += 1;
    }
    for j in 1..=npoints {
        cd[j as usize] = temp[j as usize];
    }
}
/// C: `void wt1(float *cd, enum ToCorD tocord, int npoints, WtFilter *wfilt, WtSizes *wtsizes, int idim);`
/// Transform the one-based `cd[1..=npoints]`, with power-of-two length.
/// `idim` indexes the one-based resolution levels. The coarsest length
/// `npoints >> Mraleveln[idim]` must be at least four; otherwise the C
/// diagnostics are printed and the process exits with status zero.
pub fn wt1(cd: &mut [f32], tocord: ToCorD, npoints: i32, f: &WtFilter, s: &WtSizes, idim: i32) {
    let last = npoints >> s.Mraleveln[idim as usize];
    if last < 4 {
        oxwave_core::fprintf!(oxwave_core::cio::stderr(), "lenlastc=%d\n", last);
        oxwave_core::fprintf!(
            oxwave_core::cio::stderr(),
            "StLevel too large or npoints too small.\n"
        );
        oxwave_core::exit(0);
    }
    if tocord == ToCorD::ToD {
        let mut n = npoints;
        while n >= 2 * last {
            wto1d(cd, n, tocord, f);
            n >>= 1;
        }
    } else {
        let mut n = 2 * last;
        while n <= npoints {
            wto1d(cd, n, tocord, f);
            n <<= 1;
        }
    }
}
/// C: `void wtn(float *cd, enum ToCorD tocord, WtFilter *wfilt, WtSizes *wtsizes, int dconly);`
/// Transform one-based, contiguous n-dimensional data, first dimension fastest.
/// Each dimension has positive length and a nonnegative resolution level with
/// nonzero coarsest length. Unlike `wt1`, the C loops include that coarsest
/// length. With `dconly == 1` on a forward transform, retain only the first
/// two dimensions' DC rectangle, exactly as C does (even for higher dimensions).
pub fn wtn(cd: &mut [f32], tocord: ToCorD, f: &WtFilter, s: &WtSizes, dconly: i32) {
    let mut total = 1i64;
    for d in 1..=s.NDim as usize {
        total *= s.NPointsn[d] as i64;
    }
    let mut temp = vec![0.0; total as usize + 1];
    let mut prev = 1i64;
    for d in 1..=s.NDim as usize {
        let n = s.NPointsn[d];
        let last = n >> s.Mraleveln[d];
        let new = n as i64 * prev;
        if n > last {
            for i2 in (0..total).step_by(new as usize) {
                for i1 in 1..=prev {
                    for k in 1..=n {
                        temp[k as usize] = cd[(i1 + i2 + (k - 1) as i64 * prev) as usize];
                    }
                    if tocord == ToCorD::ToC {
                        let mut nt = last;
                        while nt <= n {
                            wto1d(&mut temp, nt, tocord, f);
                            nt <<= 1;
                        }
                    } else {
                        let mut nt = n;
                        while nt >= last {
                            wto1d(&mut temp, nt, tocord, f);
                            nt >>= 1;
                        }
                    }
                    for k in 1..=n {
                        cd[(i1 + i2 + (k - 1) as i64 * prev) as usize] = temp[k as usize];
                    }
                }
            }
        }
        prev = new;
    }
    if dconly == 1 && tocord == ToCorD::ToD {
        let old = cd[..=total as usize].to_vec();
        cd[1..=total as usize].fill(0.0);
        for iz in 1..=(s.NPointsn[1] >> s.Mraleveln[1]) as usize {
            for jx in 0..(s.NPointsn[2] >> s.Mraleveln[2]) as usize {
                let ix = iz + s.NPointsn[1] as usize * jx;
                cd[ix] = old[ix];
            }
        }
    }
}
/// C: `void fhierfromcd(float *f, float *cd, WtFilter *wfilt, WtSizes *wtsizes, int justDC);`
/// Calculate the two-dimensional MRA display, retaining the C routine's mixture
/// of zero-based display indices and one-based transforms. Buffers contain
/// `NPointsn[1] * NPointsn[2] + 1` elements. Dimensions must be powers of two,
/// levels must leave nonzero coarsest lengths, and at least one level is positive.
/// As in C, `justDC` is forwarded to inverse `wtn`, where it has no effect.
pub fn fhierfromcd(f: &mut [f32], cd: &[f32], wfilt: &WtFilter, wtsizes: &WtSizes, justDC: i32) {
    let n1 = wtsizes.NPointsn[1];
    let n2 = wtsizes.NPointsn[2];
    let ntot = (n2 * n1) as i64;
    let mut workf = vec![0.0f32; (ntot + 1) as usize];
    let mut workcd = vec![0.0f32; (ntot + 1) as usize];
    let StL1 = wtsizes.Mraleveln[1];
    let StL2 = wtsizes.Mraleveln[2];
    let StLevel = if StL1 <= StL2 { StL1 } else { StL2 };
    for i in 0..=ntot {
        f[i as usize] = 0.0;
    }
    for level in 0..StLevel {
        let sizez1 = n1 >> level;
        let sizex2 = n2 >> level;
        let leap = 2 << level;
        let leap2 = (leap * leap) as f32;
        for i in 0..ntot + 1 {
            workcd[i as usize] = 0.0;
        }
        for iz in sizez1 / 2..sizez1 {
            for jx in sizex2 / 2..sizex2 {
                workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
            }
        }
        wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
        for i in 0..ntot + 1 {
            workf[i as usize] = workcd[i as usize];
        }
        for iz in sizez1 / 2..sizez1 {
            for jx in sizex2 / 2..sizex2 {
                let mut sum = 0.0f32;
                for i in 0..leap {
                    for j in 0..leap {
                        sum += workf[(iz * leap + i - n1 + (jx * leap - n2 + j) * n1) as usize];
                    }
                }
                f[(iz + jx * n1) as usize] = sum / leap2;
            }
        }
        for i in 0..ntot + 1 {
            workcd[i as usize] = 0.0;
        }
        for iz in sizez1 / 2..sizez1 {
            for jx in 0..sizex2 / 2 {
                workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
            }
        }
        wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
        for i in 0..ntot + 1 {
            workf[i as usize] = workcd[i as usize];
        }
        for iz in sizez1 / 2..sizez1 {
            for jx in 0..sizex2 / 2 {
                let mut sum = 0.0f32;
                for i in 0..leap {
                    for j in 0..leap {
                        sum += workf[(iz * leap + i - n1 + (jx * leap + j) * n1) as usize];
                    }
                }
                f[(iz + jx * n1) as usize] = sum / leap2;
            }
        }
        for i in 0..ntot + 1 {
            workcd[i as usize] = 0.0;
        }
        for iz in 0..sizez1 / 2 {
            for jx in sizex2 / 2..sizex2 {
                workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
            }
        }
        wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
        for i in 0..ntot + 1 {
            workf[i as usize] = workcd[i as usize];
        }
        for iz in 0..sizez1 / 2 {
            for jx in sizex2 / 2..sizex2 {
                let mut sum = 0.0f32;
                for i in 0..leap {
                    for j in 0..leap {
                        sum += workf[(iz * leap + i + (jx * leap + j - n2) * n1) as usize];
                    }
                }
                f[(iz + jx * n1) as usize] = sum / leap2;
            }
        }
    }
    if StL1 == StL2 {
        assert!(
            StLevel > 0,
            "C shifts by a negative level when both levels are zero"
        );
        let level = StLevel - 1;
        let sizez1 = n1 >> level;
        let sizex2 = n2 >> level;
        let leap = 2 << level;
        let leap2 = (leap * leap) as f32;
        for i in 0..ntot + 1 {
            workcd[i as usize] = 0.0;
        }
        for iz in 0..sizez1 / 2 {
            for jx in 0..sizex2 / 2 {
                workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
            }
        }
        wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
        for i in 0..ntot + 1 {
            workf[i as usize] = workcd[i as usize];
        }
        for iz in 0..sizez1 / 2 {
            for jx in 0..sizex2 / 2 {
                let mut sum = 0.0f32;
                for i in 0..leap {
                    for j in 0..leap {
                        sum += workf[(iz * leap + i + (jx * leap + j) * n1) as usize];
                    }
                }
                f[(iz + jx * n1) as usize] = sum / leap2;
            }
        }
    }
    if StL1 > StL2 {
        for level in StLevel..StL1 {
            let sizez1 = n1 >> level;
            let sizex2 = n2 >> StL2;
            let leap = 2 << level;
            let _leap2 = (leap * leap) as f32; // Computed but unused in this C branch.
            for i in 0..ntot + 1 {
                workcd[i as usize] = 0.0;
            }
            for iz in sizez1 / 2..sizez1 {
                for jx in 0..sizex2 {
                    workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
                }
            }
            wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
            for i in 0..ntot + 1 {
                workf[i as usize] = workcd[i as usize];
            }
            for iz in sizez1 / 2..sizez1 {
                for jx in 0..sizex2 {
                    let mut sum = 0.0f32;
                    for i in 0..leap {
                        sum += workf[(iz * leap + i - n1 + jx * n1) as usize];
                    }
                    f[(iz + jx * n1) as usize] = sum / leap as f32;
                }
            }
            if level == (StL1 - 1) {
                for i in 0..ntot + 1 {
                    workcd[i as usize] = 0.0;
                }
                for iz in 0..sizez1 / 2 {
                    for jx in 0..sizex2 {
                        workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
                    }
                }
                wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
                for i in 0..ntot + 1 {
                    workf[i as usize] = workcd[i as usize];
                }
                for iz in 0..sizez1 / 2 {
                    for jx in 0..sizex2 {
                        let mut sum = 0.0f32;
                        for i in 0..leap {
                            sum += workf[(iz * leap + i + jx * n1) as usize];
                        }
                        f[(iz + jx * n1) as usize] = sum / leap as f32;
                    }
                }
            }
        }
    }
    if StL2 > StL1 {
        for level in StLevel..StL2 {
            let sizez1 = n1 >> StL1;
            let sizex2 = n2 >> level;
            let leap = 2 << level;
            let _leap2 = (leap * leap) as f32; // Computed but unused in this C branch.
            for i in 0..ntot + 1 {
                workcd[i as usize] = 0.0;
            }
            for iz in 0..sizez1 {
                for jx in sizex2 / 2..sizex2 {
                    workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
                }
            }
            wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
            for i in 0..ntot + 1 {
                workf[i as usize] = workcd[i as usize];
            }
            for iz in 0..sizez1 {
                for jx in sizex2 / 2..sizex2 {
                    let mut sum = 0.0f32;
                    for i in 0..leap {
                        sum += workf[(iz + (jx * leap + i - n2) * n1) as usize];
                    }
                    f[(iz + jx * n1) as usize] = sum / leap as f32;
                }
            }
            if level == StL2 - 1 {
                for i in 0..ntot + 1 {
                    workcd[i as usize] = 0.0;
                }
                for iz in 0..sizez1 {
                    for jx in 0..sizex2 / 2 {
                        workcd[(iz + jx * n1) as usize] = cd[(iz + jx * n1) as usize];
                    }
                }
                wtn(&mut workcd, ToCorD::ToC, wfilt, wtsizes, justDC);
                for i in 0..ntot + 1 {
                    workf[i as usize] = workcd[i as usize];
                }
                for iz in 0..sizez1 {
                    for jx in 0..sizex2 / 2 {
                        let mut sum = 0.0f32;
                        for i in 0..leap {
                            sum += workf[(iz + (jx * leap + i) * n1) as usize];
                        }
                        f[(iz + jx * n1) as usize] = sum / leap as f32;
                    }
                }
            }
        }
    }
}
