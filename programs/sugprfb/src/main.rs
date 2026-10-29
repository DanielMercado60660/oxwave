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

//! SUGPRFB - remove first breaks from GPR data.
//! Port of su/main/windowing_sorting_muting/sugprfb.c (Revision 1.3).

use oxwave_core::prelude::*;
use oxwave_array_extrema::max_index;
use oxwave_correlation::xcor;
use oxwave_qr::sqrst;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "SUGPRFB - SU program to remove First Breaks from GPR data\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  sugprfb < radar traces >outfile\t\t\t  \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " nx=51\t\tnumber of traces to sum to create pilot trace (odd)\t",
    " fbt=60\tlength of first break in number of samples\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The first fbt samples from nx traces are stacked to form a pilot\t",
    " first break trace, this is fitted to the actual traces by shifting\t",
    " and scaling.\t\t The nx traces long spatial window is\t\t",
    " slided along the section and a new pilot traces is formed for each\t",
    " position. The scalers in percent and the time shifts are stored in\t",
    " header words trwf and grnors.\t\t\t\t\t  \t",
    " \t\t\t\t\t\t\t\t\t",
];

fn remove_fb(yp: &mut [f32], ym: &[f32], n: usize) -> (i16, i16) {
    let mut y = vec![0.0; n];
    xcor(0, ym, 0, &yp[..n], -(n as i32)/2, &mut y);
    let b = -max_index(&y, 1) + n as i32 / 2;
    let ns = n - b.unsigned_abs() as usize;
    let mut x = ym.to_vec();
    if b < 0 {
        // The original reads before the pilot allocation; use zero padding.
        for i in 0..ns {
            let j = i as i32 + b;
            x[i] = if j < 0 { 0.0 } else { ym[j as usize] };
        }
    } else {
        for i in 0..ns { x[i + b as usize] = ym[i]; }
    }
    let mut a = [1.0];
    if ns > 0 {
        let mut res = vec![0.0; ns];
        let mut k = 0;
        sqrst(&mut [&mut x[..]], ns as i32, 1, &yp[..ns], 0.0,
              &mut a, &mut res, &mut k, &mut [0], &mut [0.0], &mut [0.0]);
    }
    let ramps = n - n / 5;
    let mut ramp = vec![1.0f32; n];
    for it in ramps..n {
        ramp[it] = (1.0f64 - ((it-ramps) as f32 / (n-ramps-1) as f32) as f64) as f32;
    }
    if b < 0 {
        for it in 0..ns { let j = (it as i32 - b) as usize; yp[j] -= ym[it] * a[0] * ramp[j]; }
    } else {
        for it in 0..ns {
            let j = it as i32 - b;
            let sample = if j < 0 { 0.0 } else { ym[j as usize] };
            yp[it] -= sample * a[0] * ramp[it + b as usize];
        }
    }
    let v = (1.0f64 - a[0] as f64) * 100.0;
    let scaler = if v >= 0.0 { (v + 0.5) as i32 } else { (v - 0.5) as i32 };
    (scaler as i16, b as i16)
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut nx = 51i32;
    getparint("nx", &mut nx);
    if nx % 2 == 0 { nx += 1; warn!(" nx has been changed to %d to be odd.\n", nx); }
    let mut fbt = 60i32;
    getparint("fbt", &mut fbt);
    checkpars();
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    if nx < 0 || fbt < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    if fbt == 0 { err!("fbt must be positive"); }
    let nt = tr.ns as usize;
    if fbt as usize > nt { err!("fbt exceeds number of trace samples"); }
    let nx = nx as usize;
    let fbt = fbt as usize;
    let mut stacked = vec![0.0f32; fbt];
    let mut traces = Vec::new();
    if traces.try_reserve_exact(nx).is_err() { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    traces.push(Trace { hdr: tr.hdr, data: tr.data[..nt].to_vec() });
    for _ in 1..nx { gettr(&mut tr); traces.push(Trace { hdr: tr.hdr, data: tr.data[..nt].to_vec() }); }
    for t in &traces { for it in 0..fbt { stacked[it] += t.data[it]; } }
    for v in &mut stacked { *v /= nx as f32; }
    let half = nx/2;
    for t in traces.iter().take(half+1) {
        tr.hdr = t.hdr;
        tr.data[..nt].copy_from_slice(&t.data);
        let (scaler, shft) = remove_fb(&mut tr.data, &stacked, fbt);
        tr.trwf = scaler; tr.grnors = shft; puttr(&mut tr);
    }
    gettr(&mut tr);
    loop {
        for it in 0..fbt { stacked[it] -= traces[0].data[it] / nx as f32; }
        traces.remove(0);
        traces.push(Trace { hdr: tr.hdr, data: tr.data[..nt].to_vec() });
        for it in 0..fbt { stacked[it] += traces[nx-1].data[it] / nx as f32; }
        tr.hdr = traces[half].hdr;
        tr.data[..nt].copy_from_slice(&traces[half].data);
        let (scaler, shft) = remove_fb(&mut tr.data, &stacked, fbt);
        tr.trwf = scaler; tr.grnors = shft; puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    for t in traces.iter().skip(half+1) {
        tr.hdr = t.hdr;
        tr.data[..nt].copy_from_slice(&t.data);
        let (scaler, shft) = remove_fb(&mut tr.data, &stacked, fbt);
        tr.trwf = scaler; tr.grnors = shft; puttr(&mut tr);
    }
    exit(CWP_Exit());
}
