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

//! UNISAM - uniformly sample a function. Port of par/main/resampling/unisam.c (Revision 1.15).
//! Credits: Dave Hale, Zhaobo Meng, John Stockwell, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::{cakima,cmonot,csplin};
use oxwave_f_smooth::smooth_1;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " UNISAM - UNIformly SAMple a function y(x) specified as x,y pairs\t",
    " \t\t\t\t\t\t\t\t\t",
    "   unisam xin= yin= nout= [optional parameters] >binaryfile\t\t",
    "    ... or ...\t\t\t\t\t\t\t\t",
    "   unisam xfile= yfile= npairs= nout= [optional parameters] >binaryfile",
    "    ... or ...\t\t\t\t\t\t\t\t",
    "   unisam xyfile= npairs= nout= [optional parameters] >binaryfile\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " xin=,,,\tarray of x values (number of xin = number of yin)\t",
    " yin=,,,\tarray of y values (number of yin = number of xin)\t",
    "  ... or\t\t\t\t\t\t\t\t",
    " xfile=\tbinary file of x values\t\t\t\t\t",
    " yfile=\tbinary file of y values\t\t\t\t\t",
    "  ... or\t\t\t\t\t\t\t\t",
    " xyfile=\tbinary file of x,y pairs\t\t\t\t",
    " npairs=\tnumber of pairs input (active only if xfile= and yfile=\t",
    " \t\tor xyfile= are set)\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " nout=\t\t number of y values output to binary file\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dxout=1.0\t output x sampling interval\t\t\t\t",
    " fxout=0.0\t output first x\t\t\t\t\t\t",
    " method=linear  =linear for linear interpolation (continuous y)\t",
    "\t\t =mono for monotonic cubic interpolation (continuous y')",
    "\t\t =akima for Akima's cubic interpolation (continuous y') ",
    "\t\t =spline for cubic spline interpolation (continuous y'')",
    " isint=,,,\t where these sine interpolations to apply\t\t",
    " amp=,,,\t amplitude of sine interpolations\t\t\t",
    " phase0=,,,\t starting phase (defaults: 0,0,0,...,0)\t\t\t",
    " totalphase=,,, total phase (default pi,pi,pi,...,pi.)\t\t\t",
    " nwidth=0       apply window smoothing if nwidth>0                     ",
    " sloth=0\t apply interpolation in input (velocities)\t\t",
    "\t\t =1 apply interpolation to 1/input (slowness),\t\t",
    " \t\t =2 apply interpolation to 1/input (sloth), and write\t",
    " \t\t out velocities in each case.\t\t\t\t",
    " smooth=0\t apply damped least squares smoothing to output\t\t",
    " r=10\t\t  ... damping coefficient, only active when smooth=1\t",
    " \t\t\t\t\t\t\t\t\t",
];
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut xfile, mut yfile, mut xyfile) = (String::new(), String::new(), String::new());
    getparstring("xfile", &mut xfile);
    getparstring("yfile", &mut yfile);
    getparstring("xyfile", &mut xyfile);
    let (mut sloth, mut smooth, mut r) = (0, 0, 10.0f32);
    getparint("sloth", &mut sloth);
    getparint("smooth", &mut smooth);
    getparfloat("r", &mut r);
    let (mut xin, mut yin) = (Vec::<f32>::new(), Vec::<f32>::new());
    let nin: usize;
    if xfile.is_empty() && yfile.is_empty() && xyfile.is_empty() {
        nin = getparfloat("xin", &mut xin) as usize;
        if nin == 0 { err!("Must specify xin!"); }
        if getparfloat("yin", &mut yin) != nin as i32 { err!("Number of yins must equal number of xins!"); }
    } else if !xfile.is_empty() && !yfile.is_empty() {
        let mut npairs = 0;
        mustgetparint("npairs", &mut npairs);
        nin = npairs as usize;
        let fp = efopen(&xfile, "r");
        xin = vec![0.0; nin];
        if fread(&mut xin, &fp) != nin { err!("error reading xfile=%s\n", xfile); }
        let fp = efopen(&yfile, "r");
        yin = vec![0.0; nin];
        if fread(&mut yin, &fp) != nin { err!("error reading yfile=%s\n", yfile); }
    } else if !xyfile.is_empty() {
        let mut npairs = 0;
        mustgetparint("npairs", &mut npairs);
        nin = npairs as usize;
        let fp = match fopen(&xyfile, "r") { Some(fp) => fp, None => err!("cannot open xyfile=%s\n", xyfile) };
        let mut xyin = vec![0.0f32; 2*nin];
        xin = vec![0.0; nin]; yin = vec![0.0; nin];
        if fread(&mut xyin, &fp) != 2*nin { err!("error reading xyfile=%s\n", xyfile); }
        for ix in 0..2*nin { if ix%2 == 0 { xin[ix/2]=xyin[ix]; } else { yin[ix/2]=xyin[ix]; } }
    } else { nin = 0; }
    if sloth != 0 {
        for y in &mut yin { if *y != 0.0 {
            let temp = *y;
            if sloth == 1 { *y = (1.0f64 / temp as f64) as f32; }
            else if sloth == 2 { *y = (1.0f64 / (temp*temp) as f64) as f32; }
        } }
    }
    let mut nout = 0;
    if getparint("nout", &mut nout) == 0 { err!("Must specify nout!"); }
    let (mut dxout, mut fxout) = (1.0f32, 0.0f32);
    getparfloat("dxout", &mut dxout); getparfloat("fxout", &mut fxout);
    let mut method = String::from("linear"); getparstring("method", &mut method);
    let mut xout = vec![0.0f32; nout as usize];
    let mut yout = vec![0.0f32; nout as usize];
    for i in 0..nout as usize { xout[i] = fxout + i as f32 * dxout; }
    let m = method.as_bytes().first().copied().unwrap_or(0);
    if m == b'l' || nin == 1 { intlin(&xin, &yin, yin[0], yin[nin-1], &xout, &mut yout); }
    else if m == b'm' || m == b'a' || m == b's' {
        let mut yind = vec![[0.0f32; 4]; nin];
        if m == b'm' { cmonot(&xin, &yin, &mut yind); }
        else if m == b'a' { cakima(&xin, &yin, &mut yind); }
        else { csplin(&xin, &yin, &mut yind); }
        intcub(0, &xin, &yind, &xout, &mut yout);
    } else { err!("%s is an unknown interpolation method!\n", method); }
    let nsint = countparval("isint") as usize;
    let mut isint = vec![0i32; nsint]; let mut amp = vec![0.0f32; nsint];
    let mut phase0 = vec![0.0f32; nsint]; let mut totalphase = vec![0.0f32; nsint];
    if nsint > 0 && nsint < nin - 1 {
        if getparint("isint", &mut isint) == 0 { err!("Must specify %d isint values", nsint as i32); }
        if getparfloat("amp", &mut amp) == 0 { err!("Must specity %d amp values", nsint as i32); }
        getparfloat("phase0", &mut phase0);
        if getparfloat("totalphase", &mut totalphase) == 0 { totalphase.fill(std::f64::consts::PI as f32); }
        for phase in &mut totalphase { if *phase == 0.0 { *phase = std::f64::consts::PI as f32; } }
    }
    for is in 0..nsint {
        if isint[is]+1 >= nin as i32-1 { break; }
        let j = isint[is] as usize;
        let ix0 = ((xin[j]-fxout)/dxout) as i32;
        let ix1 = ((xin[j+1]-fxout)/dxout) as i32;
        let denorm = (xin[j+1]-xin[j])/totalphase[is];
        for ix in ix0+1..ix1 {
            let arg = phase0[is]+(dxout*ix as f32-xin[j]+fxout)/denorm;
            yout[ix as usize] = (yout[ix as usize] as f64 + amp[is] as f64 * (arg as f64).sin()) as f32;
        }
    }
    let mut nwidth = 0i32; getparint("nwidth", &mut nwidth);
    nwidth = nwidth/2*2+1;
    checkpars();
    if nwidth > 1 {
        let mut yout0 = vec![0.0f32; nout as usize];
        for ix in nwidth/2..nout-nwidth/2 {
            let mut sum = 0.0f32;
            for ip in -nwidth/2..=nwidth/2 { sum += yout[(ix+ip) as usize]; }
            yout0[ix as usize] = sum / nwidth as f32;
        }
        for ix in nwidth/2..nout-nwidth/2 { yout[ix as usize] = yout0[ix as usize]; }
    }
    if sloth != 0 { for y in &mut yout { if *y != 0.0 {
        let temp = *y;
        if sloth == 1 { *y = (1.0f64/temp as f64) as f32; }
        else if sloth == 2 { *y = (1.0f64/temp as f64).sqrt() as f32; }
    } } }
    if smooth != 0 { smooth_1(&xout, &mut yout, r); }
    efwrite(&yout, &stdout());
    exit(CWP_Exit());
}
