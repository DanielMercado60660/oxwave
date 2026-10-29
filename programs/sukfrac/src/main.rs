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

//! SUKFRAC - fractional powers of i|k|. Port of su/main/filters/sukfrac.c (Revision 1.15).
//! Credits: CWP: John Stockwell, June 1997, based on sufrac.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_pfa_fft::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUKFRAC - apply FRACtional powers of i|k| to data, with phase shift ",
    " \t\t\t\t\t\t\t\t\t",
    "     sukfrac <infile >outfile [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "  power=0\t\texponent of (i*sqrt(k1^2 + k2^2))^power\t\t",
    "\t\t\t=0 ===> phase shift only\t\t\t",
    "\t\t\t>0 ===> differentiation\t\t\t\t",
    "\t\t\t<0 ===> integration\t\t\t\t",
    "  sign=1\t\tsign on transform exponent       \t\t",
    "  d1=1.0\t\tx1 sampling interval\t\t\t\t",
    "  d2=1.0\t\tx2 sampling interval\t\t\t\t",
    "  phasefac=0\t\tphase shift by phase=phasefac*PI\t\t",
    "  ...directional derivative, active only if angle= is set ....\t\t",
    "  angle=\t\tif set applies operator directionally in the\t",
    "\t\t\tdirection specified by the value of angle,\t",
    "\t\t\t-360.0 <= angle < 359.99999 degress\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The relation: w = 2 pi F is well known for frequency, but there\t",
    " doesn't seem to be a commonly used letter corresponding to F for the\t",
    " spatial conjugate transform variables.  We use K1 and K2 for this.\t",
    " More specifically we assume a phase:\t\t\t\t\t",
    "\t\t-i(k1 x1 + k2 x2) = -2 pi i(K1 x1 + K2 x2).\t\t",
    " and K1, K2 define our respective wavenumbers.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Algorithms \t\t\t\t\t\t\t\t",
    " \tg(x1,x2)=Re[2DINVFFT{ ( (sign) i |k|)^power 2DFFT(f)}e^i(phase)]",
    " where: \t\t\t\t\t\t\t\t",
    "       |k| = sqrt[ (k1)^2 + (k2)^2 ] \t\t\t\t\t",
    " or when angle= is set \t\t\t\t\t\t",
    "       |k| = sqrt[ (k1  cos(angle) )^2 + (k2 sin(angle) )^2 ] \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " In the default mode a factor of (i|k|)^(power) is applied in the \t",
    " transform domain. For time data the time axis direction is taken to \t",
    " be the k1-direction.  The effect of this filter is to differentiate   ",
    " the input in the normal direction to any curvilinear features in\t",
    " in the data, and thus be a non-directional-specific edge enhancer.\t",
    " \t\t\t\t\t\t\t\t\t",
    " If angle= is set, then the intended effect is a derivative in the \t",
    " direction specified by the angle, with the k1-direction being angle=0,",
    " corresponding to curves whose normals lie in the x1-direction. \t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " Large amplitude errors will result of the data set has too few points.",
    " \t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " Edge sharpening: \t\t\t\t\t\t\t",
    " Laplacean : \t\t\t\t\t\t\t\t",
    "    sukfrac < image_data  power=2 phasefac=-1 | ... \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Image enhancement:\t\t\t\t\t\t\t",
    "   Derivative filter:\t\t\t\t\t\t\t",
    "    sukfrac < image_data  power=1 phasefac=-.5 | ... \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Image enhancement:\t\t\t\t\t\t\t",
    "   Half derivative (this one is the best for photographs): \t\t",
    "    sukfrac < image_data  power=.5 phasefac=-.25 | ... \t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    if tr.trid != 130 { warn!("tr.trid = %d", tr.trid as i32); }
    let nx1 = tr.ns as i32;
    let mut dx1 = 0.0f32;
    let mut dx2 = 0.0f32;
    if getparfloat("d1", &mut dx1) == 0 {
        if tr.d1 != 0.0 { dx1 = tr.d1; }
        else { dx1 = 1.0; warn!("tr.d1 not set, assuming d1=1.0"); }
    }
    if getparfloat("d2", &mut dx2) == 0 {
        if tr.d2 != 0.0 { dx2 = tr.d2; }
        else { dx2 = 1.0; warn!("tr.d2 not set, assuming d2=1.0"); }
    }
    let mut verbose = 1i32;
    let mut power = 0.0f32;
    let mut sign = 1.0f32;
    let mut phasefac = 0.0f32;
    let mut angle = 0.0f32;
    getparint("verbose", &mut verbose);
    getparfloat("power", &mut power);
    getparfloat("sign", &mut sign);
    getparfloat("phasefac", &mut phasefac);
    let phase = (phasefac as f64 * std::f64::consts::PI) as f32;
    let is_angle = getparfloat("angle", &mut angle) != 0;
    if is_angle { angle = (angle as f64 * std::f64::consts::PI / 180.0) as f32; }
    let arg = (sign as f64 * (0.5 * std::f64::consts::PI) * power as f64 + phase as f64) as f32;
    let exparg = cwp_cexp(crmul(cmplx(0.0, 1.0), arg));
    let mut headers = Vec::new();
    let tracefp = etmpfile();
    loop {
        headers.push(tr.hdr);
        efwrite(&tr.data[..nx1 as usize], &tracefp);
        if gettr(&mut tr) == 0 { break; }
    }
    let nx2 = headers.len() as i32;
    let nx1fft = npfaro(nx1, 2 * nx1);
    let nx2fft = npfa(nx2);
    if nx1fft >= SU_NFLTS as i32 || nx1fft >= 720720 { err!("Padded nx1=%d--too big", nx1fft); }
    if nx2fft >= SU_NFLTS as i32 || nx2fft >= 720720 { err!("Padded nx2=%d--too big", nx2fft); }
    let nk1 = nx1fft / 2 + 1;
    let nk2 = nx2fft / 2 + 1;
    let onk1 = (1.0f64 / nk1 as f64) as f32;
    let onk2 = (1.0f64 / nk2 as f64) as f32;
    let dk1 = (2.0 * std::f64::consts::PI * onk1 as f64 / dx1 as f64) as f32;
    let dk2 = (2.0 * std::f64::consts::PI * onk2 as f64 / dx2 as f64) as f32;
    let mut rt = vec![0.0f32; (nx1fft * nx2fft) as usize];
    let mut ct = vec![Complex::default(); (nk1 * nx2fft) as usize];
    let mut kfilter = vec![Complex::default(); (nx1fft * nx2fft) as usize];
    for ik2 in 0..nx2fft {
        for ik1 in 0..nk1 {
            let j = if ik2 < nk2 { ik2 } else { nx2fft - ik2 };
            let a = (ik1 * ik1) as f32 * dk1 * dk1;
            let b = (j * j) as f32 * dk2 * dk2;
            let mut k = if is_angle {
                let c = (angle as f64).cos();
                let s = (angle as f64).sin();
                (a as f64 * c * c + b as f64 * s * s).sqrt() as f32
            } else { ((a + b) as f64).sqrt() as f32 };
            if power < 0.0 && ik1 == 0 && ik2 == 0 { k = f32::MAX; }
            let amp = ((k as f64).powf(power as f64) * onk1 as f64 * onk2 as f64 / std::f64::consts::PI) as f32;
            kfilter[(ik2 * nx1fft + ik1) as usize] = crmul(exparg, amp);
        }
    }
    erewind(&tracefp);
    for ix2 in 0..nx2 as usize {
        let dest = ix2 * nx1fft as usize;
        efread(&mut rt[dest..dest + nx1 as usize], &tracefp);
    }
    efclose(&tracefp);
    pfa2rc(-1, 1, nx1fft, nx2, &rt, &mut ct);
    pfa2cc(-1, 2, nk1, nx2fft, &mut ct);
    for ik2 in 0..nx2fft {
        for ik1 in 0..nk1 {
            let i = (ik2 * nk1 + ik1) as usize;
            ct[i] = cmul(ct[i], kfilter[(ik2 * nx1fft + ik1) as usize]);
        }
    }
    pfa2cc(1, 2, nk1, nx2fft, &mut ct);
    pfa2cr(1, 1, nx1fft, nx2, &ct, &mut rt);
    for ix2 in 0..nx2 as usize {
        tr.hdr = headers[ix2];
        let start = ix2 * nx1fft as usize;
        tr.data[..nx1 as usize].copy_from_slice(&rt[start..start + nx1 as usize]);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
