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
//! SUDIPFILT - slope filter in the f-k domain, port of su/main/filters/sudipfilt.c.
//! Credits: CWP: Dave (algorithm), Jack (SU formatting).
use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, cmplx, cmul};
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfa, npfar, pfa2rc, pfa2cc, pfa2cr};
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;
unsafe extern "C" {
    fn access(path: *const i8, mode: i32) -> i32;
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUDIPFILT - DIP--or better--SLOPE Filter in f-k domain\t",
    " \t\t\t\t\t\t\t\t",
    " sudipfilt <infile >outfile [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t",
    " dt=(from header)\tif not set in header then mandatory\t",
    " dx=(from header, d1)\tif not set in header then mandatory\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " slopes=0.0\t\tmonotonically increasing slopes\t\t",
    " amps=1.0\t\tamplitudes corresponding to slopes\t",
    " bias=0.0\t\tslope made horizontal before filtering\t",
    "\t\t\t\t\t\t\t\t",
    " verbose=0\tverbose = 1 echoes information\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path",
    "\t\t prefix for storing temporary files; else if the",
    "\t         the CWP_TMPDIR environment variable is set use\t",
    "\t         its value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " d2 is an acceptable alias for dx in the getpar\t\t",
    "\t\t\t\t\t\t\t\t",
    " Slopes are defined by delta_t/delta_x, where delta\t\t",
    " means change. Units of delta_t and delta_x are the same\t",
    " as dt and dx. It is sometimes useful to fool the program\t",
    " with dx=1 dt=1, thus avoiding units and small slope values.\t",
    "\t\t\t\t\t\t\t\t",
    " Linear interpolation and constant extrapolation are used to\t",
    " determine amplitudes for slopes that are not specified.\t",
    " Linear moveout is removed before slope filtering to make\t",
    " slopes equal to bias appear horizontal.  This linear moveout\t",
    " is restored after filtering.  The bias parameter may be\t",
    " useful for spatially aliased data.  The slopes parameter is\t",
    " compensated for bias, so you need not change slopes when you\t",
    " change bias.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
];
fn main() {
    let argv = args(); let mut tr = Trace::new();
    initargs(&argv); requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut verbose = 0; getparint("verbose", &mut verbose);
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 { tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default(); }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_str()).unwrap();
        if unsafe { access(path.as_ptr(), 2) } != 0 { err!("you can't write in %s (or it doesn't exist)", tmpdir); }
    }
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (tr.dt as f64 / 1000000.0) as f32; }
    if dt == 0.0 { err!("dt field is zero and not getparred"); }
    let mut dx = 0.0f32;
    if getparfloat("dx", &mut dx) == 0 && getparfloat("d2", &mut dx) == 0 { dx = tr.d2; }
    if dx == 0.0 { err!("d2 field is zero and dx not getparred"); }
    let mut slopes = vec![0.0f32; countparval("slopes") as usize];
    if slopes.is_empty() { slopes.push(0.0); } else { getparfloat("slopes", &mut slopes); }
    let mut amps = vec![0.0f32; countparval("amps") as usize];
    if amps.is_empty() { amps.push(1.0); warn!("no amps given--doing no-op"); } else { getparfloat("amps", &mut amps); }
    let mut bias = 0.0f32; getparfloat("bias", &mut bias);
    if slopes.len() != amps.len() { err!("number of slopes (%d) must equal number of amps(%d)", slopes.len() as i32, amps.len() as i32); }
    for i in 1..slopes.len() { if slopes[i] <= slopes[i-1] { err!("slopes must be monotonically increasing"); } }
    let (tracefp, headerfp, tracefile, headerfile) = if tmpdir.is_empty() {
        let a = etmpfile(); let b = etmpfile();
        if verbose != 0 { warn!("using tmpfile() call"); }
        (a,b,String::new(),String::new())
    } else {
        let a = temporary_filename(&tmpdir); let mut b = temporary_filename(&tmpdir);
        // The library's fallback name is fixed when its /tmp creation fails.
        if a == b { b.push_str("_header"); }
        let af = efopen(&a, "w+"); let bf = efopen(&b, "w+");
        if verbose != 0 { warn!("putting temporary files in %s", tmpdir); }
        (af,bf,a,b)
    };
    let mut ntr = 0;
    loop {
        ntr += 1;
        let header = unsafe { std::slice::from_raw_parts((&tr.hdr as *const _) as *const u8, HDRBYTES) };
        efwrite(header, &headerfp);
        efwrite(&tr.data[..nt], &tracefp);
        if gettr(&mut tr) == 0 { break; }
    }
    slopefilter(&slopes, &amps, bias, nt, dt, ntr, dx, &tracefp);
    erewind(&headerfp); erewind(&tracefp);
    for _ in 0..ntr {
        let header = unsafe { std::slice::from_raw_parts_mut((&mut tr.hdr as *mut _) as *mut u8, HDRBYTES) };
        efread(header, &headerfp);
        efread(&mut tr.data[..nt], &tracefp); puttr(&mut tr);
    }
    efclose(&headerfp); if !headerfile.is_empty() { eremove(&headerfile); }
    efclose(&tracefp); if !tracefile.is_empty() { eremove(&tracefile); }
    exit(CWP_Exit());
}
fn slopefilter(slopes: &[f32], amps: &[f32], bias: f32, nt: usize, dt: f32, nx: usize, dx: f32, tracefp: &CFile) {
    let ntfft = npfar(nt as i32) as usize; let nxfft = npfa(nx as i32) as usize;
    let sfft = (1.0 / (ntfft * nxfft) as f64) as f32;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * std::f64::consts::PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = (0.000001 * dw as f64) as f32; let nk = nxfft;
    let dk = (2.0 * std::f64::consts::PI / (nxfft as f32 * dx) as f64) as f32;
    let mut cpfft = vec![Complex::default(); nw * nk]; let mut pfft = vec![0.0f32; ntfft * nxfft];
    erewind(tracefp);
    for ix in 0..nx { efread(&mut pfft[ix*ntfft..ix*ntfft+nt], tracefp); }
    pfa2rc(1, 1, ntfft as i32, nx as i32, &pfft, &mut cpfft);
    for ix in 0..nx {
        let mut w = 0.0f32;
        for iw in 0..nw {
            let phase = -(ix as f32) * dx * w * bias;
            let cshift = cmplx((phase as f64).cos() as f32, (phase as f64).sin() as f32);
            cpfft[ix*nw+iw] = cmul(cpfft[ix*nw+iw], cshift); w += dw;
        }
    }
    pfa2cc(-1, 2, nw as i32, nxfft as i32, &mut cpfft);
    for ik in 0..nk {
        let k = if ik <= nk/2 { ik as f32 * dk } else { (ik as i32 - nk as i32) as f32 * dk };
        let mut w = fw;
        for iw in 0..nw {
            let slope = k / w + bias; let mut amp = [0.0f32];
            intlin(slopes, amps, amps[0], amps[amps.len()-1], &[slope], &mut amp); amp[0] *= sfft;
            cpfft[ik*nw+iw].r *= amp[0]; cpfft[ik*nw+iw].i *= amp[0]; w += dw;
        }
    }
    pfa2cc(1, 2, nw as i32, nxfft as i32, &mut cpfft);
    for ix in 0..nx {
        let mut w = 0.0f32;
        for iw in 0..nw {
            let phase = ix as f32 * dx * w * bias;
            let cshift = cmplx((phase as f64).cos() as f32, (phase as f64).sin() as f32);
            cpfft[ix*nw+iw] = cmul(cpfft[ix*nw+iw], cshift); w += dw;
        }
    }
    pfa2cr(-1, 1, ntfft as i32, nx as i32, &cpfft, &mut pfft);
    erewind(tracefp); for ix in 0..nx { efwrite(&pfft[ix*ntfft..ix*ntfft+nt], tracefp); }
}
