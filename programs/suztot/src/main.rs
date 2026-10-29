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

//! SUZTOT - resample from depth to time, port of su/main/stretching_moveout_resamp/suztot.c (Revision 1.10).
//! Credits: CWP: John Stockwell, 2005; based on suttoz.c by Dave Hale, c. 1992.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_inverse_interp::yxtoxy;
use oxwave_linear_interp::intlin;
use oxwave_sinc_interp::ints8r;
use std::io::Read;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUZTOT - resample from depth to time\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suztot <stdin >stdout [optional parms]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nt=1+(nz-1)*2.0*dz/(vmax*dt)    number of time samples in output\t",
    " dt=2*dz/vmin\t\ttime sampling interval (defaults avoids aliasing)",
    " ft=2*fz/v(fz)\t\tfirst time sample\t\t\t\t",
    " z=0.0,...\t\tdepths corresponding to interval velocities in v",
    " v=1500.0,...\t\tinterval velocities corresponding to depths in v",
    " vfile=\t\tbinary (non-ascii) file containing velocities v(z)",
    " verbose=0\t\t>0 to print depth sampling information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Default value of nt set to avoid aliasing\t\t\t\t",
    " The z and v arrays specify an interval velocity function of depth.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note that z and v are given  as arrays of floats separated by commas,  ",
    " for example:\t\t\t\t\t\t\t\t",
    " z=0.0,100,200,... v=1500.0,1720.0,1833.5,... with the number of z values",
    " equaling the number of v values. The velocities are linearly interpolated",
    " to produce a piecewise linear v(z) profile. This fact must be taken into",
    " account when attempting to use this program as the inverse of suttoz.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Linear interpolation and constant extrapolation is used to determine\t",
    " interval velocities at times not specified.  Values specified in z\t",
    " must increase monotonically.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Alternatively, interval velocities may be stored in a binary file\t",
    " containing one velocity for every time sample.  If vfile is specified,",
    " then the z and v arrays are ignored.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " see the selfdoc of   suttoz  for time to depth conversion\t\t",
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
    let nz = tr.ns as usize;
    let mut dz = 0.0f32;
    if getparfloat("dz", &mut dz) == 0 { dz = tr.d1; }
    let mut vz = vec![0.0f32; nz];
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        let nzpar = countparval("z").max(1) as usize;
        let mut zpar = vec![0.0f32; nzpar];
        getparfloat("z", &mut zpar);
        let nvpar = countparval("v").max(1) as usize;
        if nvpar != nzpar { err!("number of t and v values must be equal"); }
        let mut vpar = vec![0.0f32; nvpar];
        if getparfloat("v", &mut vpar) == 0 { vpar[0] = 1500.0; }
        for izpar in 1..nzpar {
            if zpar[izpar] <= zpar[izpar-1] { err!("zpar must increase monotonically"); }
        }
        let mut z = 0.0f32;
        for value in &mut vz {
            intlin(&zpar, &vpar, vpar[0], vpar[nzpar-1], &[z], std::slice::from_mut(value));
            z += dz;
        }
    } else {
        let mut bytes = vec![0u8; nz * 4];
        let ok = std::fs::File::open(&vfile).and_then(|mut f| f.read_exact(&mut bytes)).is_ok();
        if !ok { err!("cannot read %d velocities from file %s", nz as i32, vfile); }
        for (v, b) in vz.iter_mut().zip(bytes.chunks_exact(4)) {
            *v = f32::from_ne_bytes([b[0], b[1], b[2], b[3]]);
        }
    }
    let mut vmin = vz[0];
    let mut vmax = vz[0];
    for &v in &vz[1..] { if v < vmin { vmin = v; } if v > vmax { vmax = v; } }
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 { dt = (2.0 * dz as f64 / vmin as f64) as f32; }
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 { ft = (2.0 * ft as f64 / vz[0] as f64) as f32; }
    let mut fz = 0.0f32;
    getparfloat("fz", &mut fz);
    let mut nt = 0i32;
    if getparint("nt", &mut nt) == 0 {
        nt = (1.0 + (((nz as i32 - 1) as f32 * dz) as f64 * 2.0 / (dt * vmax) as f64)) as i32;
    }
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);
    checkpars();
    CHECK_NT("nt", nt);
    if verbose != 0 {
        warn!("Input:");
        warn!("\tnumber of depth samples = %d", nz as i32);
        warn!("\tdepth sampling interval = %g", dz);
        warn!("\tfirst depth sample = %g", fz);
        warn!("Output:");
        warn!("\tnumber of time samples = %d", nt);
        warn!("\ttime sampling interval = %g", dt);
        warn!("\tfirst time sample = %g", ft);
    }
    if nt < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut zt = vec![0.0f32; nt as usize];
    makezt(dz, fz, &vz, dt, ft, &mut zt);
    let mut temp = vec![0.0f32; nz];
    loop {
        tr.trid = 1;
        tr.ns = nt as u16;
        tr.dt = (dt as f64 * 1000000.0) as u32 as u16;
        tr.f1 = ft;
        tr.d1 = 0.0;
        temp.copy_from_slice(&tr.data[..nz]);
        ints8r(dz, fz, &temp, 0.0, 0.0, &zt, &mut tr.data[..nt as usize]);
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}

fn makezt(dz: f32, fz: f32, v: &[f32], dt: f32, ft: f32, zt: &mut [f32]) {
    let nz = v.len();
    let nt = zt.len();
    let mut tz = vec![0.0f32; nz];
    tz[0] = (2.0 * fz as f64 / v[0] as f64) as f32;
    for iz in 1..nz {
        tz[iz] = (tz[iz-1] as f64 + 2.0 * dz as f64 / v[iz-1] as f64) as f32;
    }
    let lt = ft + (nt as i32 - 1) as f32 * dt;
    yxtoxy(dz, fz, &tz, dt, ft, 0.0, 0.0, zt);
    let mut it = 0usize;
    let mut t = ft;
    while it < nt && t <= tz[0] {
        zt[it] = ((v[0] * t) as f64 / 2.0) as f32;
        it += 1; t += dt;
    }
    let mut it = nt as i32 - 1;
    t = lt;
    while it >= 0 && t >= tz[nz-1] {
        zt[it as usize] = ((lt * v[nz-1]) as f64 / 2.0 + (v[nz-1] * (t - tz[nz-1])) as f64 / 2.0) as f32;
        it -= 1; t -= dt;
    }
}
