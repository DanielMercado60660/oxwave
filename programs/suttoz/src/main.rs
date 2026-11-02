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

//! SUTTOZ - resample from time to depth; port of suttoz.c revision 1.20.
//! Credits: CWP: Dave Hale c. 1992.

use oxwave_core::prelude::*;
use oxwave_inverse_interp::yxtoxy;
use oxwave_linear_interp::intlin;
use oxwave_sinc_interp::ints8r;
use std::io::Read;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUTTOZ - resample from time to depth\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suttoz <stdin >stdout [optional parms]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nz=1+(nt-1)*dt*vmax/(2.0*dz)   number of depth samples in output\t",
    " dz=vmin*dt/2\t\tdepth sampling interval (defaults avoids aliasing)",
    " fz=v(ft)*ft/2\t\tfirst depth sample\t\t\t\t",
    " t=0.0,...\t\ttimes corresponding to interval velocities in v",
    " v=1500.0,...\t\tinterval velocities corresponding to times in v",
    " vfile=\t\t  binary (non-ascii) file containing velocities v(t)",
    " verbose=0\t\t>0 to print depth sampling information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Default value of nz set to avoid aliasing\t\t\t\t",
    " The t and v arrays specify an interval velocity function of time.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note that t and v are given  as arrays of floats separated by commas,  ",
    " for example:\t\t\t\t\t\t\t\t",
    " t=0.0,0.01,.2,... v=1500.0,1720.0,1833.5,... with the number of t values",
    " equaling the number of v values. The velocities are linearly interpolated",
    " to make a continuous, piecewise linear v(t) profile.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Linear interpolation and constant extrapolation is used to determine\t",
    " interval velocities at times not specified.  Values specified in t\t",
    " must increase monotonically.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Alternatively, interval velocities may be stored in a binary file\t",
    " containing one velocity for every time sample.  If vfile is specified,",
    " then the t and v arrays are ignored.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " see selfdoc of suztot  for depth to time conversion\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed:  ns, dt, and delrt\t\t\t",
    " Trace header fields modified:  trid, ns, d1, and f1\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut vt = vec![0.0f32; nt];
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        let mut ntpar = countparval("t");
        if ntpar == 0 { ntpar = 1; }
        let mut tpar = vec![0.0f32; ntpar as usize];
        if getparfloat("t", &mut tpar) == 0 { tpar[0] = 0.0; }
        let mut nvpar = countparval("v");
        if nvpar == 0 { nvpar = 1; }
        if nvpar != ntpar { err!("number of t and v values must be equal"); }
        let mut vpar = vec![0.0f32; nvpar as usize];
        if getparfloat("v", &mut vpar) == 0 { vpar[0] = 1500.0; }
        for i in 1..ntpar as usize {
            if tpar[i] <= tpar[i-1] { err!("tpar must increase monotonically"); }
        }
        let mut t = 0.0f32;
        for value in &mut vt {
            intlin(&tpar, &vpar, vpar[0], vpar[nvpar as usize-1], &[t], std::slice::from_mut(value));
            t += dt;
        }
    } else {
        let mut bytes = vec![0u8; nt * 4];
        let result = std::fs::File::open(&vfile).and_then(|mut f| f.read_exact(&mut bytes));
        if result.is_err() { err!("cannot read %d velocities from file %s", nt as i32, vfile); }
        for i in 0..nt { vt[i] = f32::from_ne_bytes(bytes[4*i..4*i+4].try_into().unwrap()); }
    }
    let mut vmin = vt[0];
    let mut vmax = vt[0];
    for &v in &vt[1..] { if v < vmin { vmin = v; } if v > vmax { vmax = v; } }
    let mut dz = 0.0f32;
    if getparfloat("dz", &mut dz) == 0 { dz = ((vmin * dt) as f64 / 2.0) as f32; }
    let mut fz = 0.0f32;
    if getparfloat("fz", &mut fz) == 0 { fz = ((vt[0] * ft) as f64 / 2.0) as f32; }
    let mut nz = 0i32;
    if getparint("nz", &mut nz) == 0 { nz = (1.0 + (((nt as i32 - 1) as f32 * dt * vmax) as f64 / (2.0 * dz as f64))) as i32; }
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);
    checkpars();
    CHECK_NT("nz", nz);
    if verbose != 0 {
        warn!("Input:");
        warn!("\tnumber of time samples = %d", nt as i32);
        warn!("\ttime sampling interval = %g", dt);
        warn!("\tfirst time sample = %g", ft);
        warn!("Output:\n");
        warn!("\tnumber of depth samples = %d", nz);
        warn!("\tdepth sampling interval = %g", dz);
        warn!("\tfirst depth sample = %g", fz);
    }
    let mut tz = vec![0.0f32; nz as usize];
    let mut temp = vec![0.0f32; nt];
    maketz(dt, ft, &vt, dz, fz, &mut tz);
    loop {
        tr.trid = 130;
        tr.ns = nz as u16;
        tr.d1 = dz;
        tr.f1 = fz;
        temp.copy_from_slice(&tr.data[..nt]);
        ints8r(dt, ft, &temp, 0.0, 0.0, &tz, &mut tr.data[..nz as usize]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}

fn maketz(dt: f32, ft: f32, v: &[f32], dz: f32, fz: f32, t: &mut [f32]) {
    let nt = v.len();
    let nz = t.len();
    let mut zt = vec![0.0f32; nt];
    zt[0] = (0.5 * ft as f64 * v[0] as f64) as f32;
    for it in 1..nt { zt[it] = (zt[it-1] as f64 + 0.5 * dt as f64 * v[it-1] as f64) as f32; }
    let lt = ft + (nt as i32 - 1) as f32 * dt;
    let lz = fz + (nz as i32 - 1) as f32 * dz;
    yxtoxy(dt, ft, &zt, dz, fz, 0.0, 0.0, t);
    let mut iz = 0usize;
    let mut z = fz;
    while iz < nz && z <= zt[0] {
        t[iz] = (2.0 * z as f64 / v[0] as f64) as f32;
        iz += 1; z += dz;
    }
    let mut iz = nz as isize - 1;
    let mut z = lz;
    while iz >= 0 && z >= zt[nt-1] {
        t[iz as usize] = (lt as f64 + 2.0 * (z - zt[nt-1]) as f64 / v[nt-1] as f64) as f32;
        iz -= 1; z -= dz;
    }
}
