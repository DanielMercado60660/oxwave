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
//! SUPWS - phase stack or phase-weighted stack of adjacent traces.
//! Port of su/main/stacking/supws.c.
//! Credits: Nils Maercklin, GeoForschungsZentrum Potsdam, 2001.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::nint;
use oxwave_hilbert::hilbert;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUPWS - Phase stack or phase-weighted stack (PWS) of adjacent traces\t",
    "\t having the same key header word\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " supws <stdin >stdout [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t \t",
    "\tkey=cdp\t   key header word to stack on\t\t\t\t",
    "\tpwr=1.0\t   raise phase stack to power pwr\t\t\t",
    "\tdt=(from header)  time sampling intervall in seconds\t\t",
    "\tsl=0.0\t\twindow length in seconds used for smoothing\t",
    "\t\t\tof the phase stack (weights)\t\t\t",
    "\tps=0\t\t0 = output is PWS, 1 = output is phase stack\t",
    "\tverbose=0\t 1 = echo additional information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note:\t\t\t\t\t\t\t\t \t",
    "\tPhase weighted stacking is a tool for efficient incoherent noise",
    "\treduction. An amplitude-unbiased coherency measure is designed\t",
    "\tbased on the instantaneous phase, which is used to weight the\t",
    "\tsamples of an ordinary, linear stack. The result is called the\t",
    "\tphase-weighted stack (PWS) and is cleaned from incoherent noise.",
    "\tPWS thus permits detection of weak but coherent arrivals.\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tThe phase-stack (coherency measure) has values between 0 and 1.\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tIf the stacking is over cdp and the PWS option is set, then the\t",
    "\toffset header field is set to zero. Otherwise, output traces get",
    "\ttheir headers from the first trace of each data ensemble to stack,",
    "\tincluding the offset field. Use \"sushw\" afterwards, if this is",
    "\tnot acceptable.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn do_smooth(data: &mut [f32], nt: usize, isl: i32) {
    let mut tmpdata = vec![0.0f32; nt];
    for it in 0..nt as i32 {
        let mut sval = 0.0f32;
        if it >= isl / 2 && it < nt as i32 - isl / 2 {
            for jt in it - isl / 2..it + isl / 2 { sval += data[jt as usize]; }
            tmpdata[it as usize] = sval / isl as f32;
        } else { tmpdata[it as usize] = data[it as usize]; }
    }
    data.copy_from_slice(&tmpdata);
}

fn main() {
    let argv = args();
    let mut intr = Trace::new();
    let mut outtr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut intr) == 0 { err!("can't get first trace"); }
    let nt = intr.ns as usize;
    let mut key = String::new();
    let mut pwr = 1.0f32;
    let mut dt = (intr.dt as f64 / 1000000.0) as f32;
    let mut sl = 0.0f32;
    let mut ps = 0i32;
    let mut verbose = 0i32;
    if getparstring("key", &mut key) == 0 { key = "cdp".to_string(); }
    getparfloat("pwr", &mut pwr);
    getparfloat("dt", &mut dt);
    getparfloat("sl", &mut sl);
    getparint("ps", &mut ps);
    getparint("verbose", &mut verbose);
    checkpars();
    let pws_and_cdp = key == "cdp" && ps == 0;
    if dt == 0.0 { dt = 0.004; warn!("dt not set, assuming dt=0.004"); }
    let isl = nint((sl as f64).abs() / dt as f64);
    if isl > nt as i32 { err!("sl=%g too long for trace", (sl as f64).abs()); }
    if verbose != 0 && isl != 0 { warn!("smoothing window = %d samples", isl); }
    let type_ = hdtype(&key);
    let index = getindex(&key);
    let mut val = Value::default();
    gethval(&intr, index, &mut val);
    let mut prevfval = vtof(type_, val);
    let mut fval = prevfval;
    let mut data = vec![0.0f32; nt];
    let mut hdata = vec![0.0f32; nt];
    let mut stdata = vec![0.0f32; nt];
    let mut psdata = vec![0.0f32; nt];
    let mut psct = vec![Complex::default(); nt];
    let mut ct = vec![Complex::default(); nt];
    let mut header = intr.hdr;
    let mut ntr = 0i32;
    let mut gottrace = true;
    loop {
        if gottrace {
            gethval(&intr, index, &mut val);
            fval = vtof(type_, val);
            data.copy_from_slice(&intr.data[..nt]);
        }
        hilbert(&data, &mut hdata);
        for i in 0..nt {
            ct[i] = cmplx(data[i], hdata[i]);
            let magnitude = rcabs(ct[i]);
            let a = if magnitude != 0.0 { (1.0f64 / rcabs(ct[i]) as f64) as f32 } else { 0.0 };
            ct[i] = crmul(ct[i], a);
        }
        if fval == prevfval && gottrace {
            ntr += 1;
            for i in 0..nt { stdata[i] += data[i]; psct[i] = cadd(psct[i], ct[i]); }
        }
        if fval != prevfval || !gottrace {
            if verbose != 0 { warn!("%s=%g, fold=%d\n", key.as_str(), prevfval, ntr); }
            for i in 0..nt {
                psdata[i] = rcabs(psct[i]) / ntr as f32;
                psdata[i] = (psdata[i] as f64).powf(pwr as f64) as f32;
            }
            if isl != 0 { do_smooth(&mut psdata, nt, isl); }
            if ps == 0 { for i in 0..nt { stdata[i] *= psdata[i] / ntr as f32; } }
            outtr.hdr = header;
            outtr.nhs = ntr as i16;
            if ps != 0 { outtr.data[..nt].copy_from_slice(&psdata); }
            else { outtr.data[..nt].copy_from_slice(&stdata); }
            if pws_and_cdp { outtr.offset = 0; }
            puttr(&mut outtr);
            if !gottrace { break; }
            prevfval = fval;
            ntr = 0;
            stdata.fill(0.0);
            psct.fill(Complex::default());
            if gottrace {
                ntr += 1;
                for i in 0..nt { stdata[i] += data[i]; psct[i] = cadd(psct[i], ct[i]); }
            }
            header = intr.hdr;
        }
        if gettr(&mut intr) == 0 { gottrace = false; }
    }
    exit(CWP_Exit());
}
