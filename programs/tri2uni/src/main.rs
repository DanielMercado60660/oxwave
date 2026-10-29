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

//! TRI2UNI - sample a triangulated sloth model on a uniform grid.
//! Port of tri/main/tri2uni.c (Revision 1.5).
//! Credits: Dave Hale, Colorado School of Mines, 04/23/91.

use oxwave_core::prelude::*;
use oxwave_lib_tri::*;
use std::ffi::c_void;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " TRI2UNI - convert a TRIangulated model to UNIformly sampled model\t",
    "\t\t\t\t\t\t\t\t\t",
    " tri2uni <triangfile >uniformfile n2= n1= [optional parameters]\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=                     number of samples in the first (fast) dimension",
    " n2=                     number of samples in the second dimension\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " d1=1.0                 sampling interval in first (fast) dimension\t",
    " d2=1.0                 sampling interval in second dimension\t\t",
    " f1=0.0                 first value in dimension 1 sampled\t\t",
    " f2=0.0                 first value in dimension 2 sampled\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note:\t\t\t\t\t\t\t\t\t",
    " The triangulated/uniformly-sampled quantity is assumed to be sloth=1/v^2",
    "\t\t\t\t\t\t\t\t\t",
];

extern "C" {
    fn fdopen(fd: i32, mode: *const i8) -> *mut c_void;
    fn raise(sig: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
}

#[repr(C)]
struct TriAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut n1, mut n2) = (0i32, 0i32);
    if getparint("n1", &mut n1) == 0 { err!("must specify n1!"); }
    if getparint("n2", &mut n2) == 0 { err!("must specify n2!"); }
    let (mut d1, mut d2, mut f1, mut f2) = (1.0f32, 1.0f32, 0.0f32, 0.0f32);
    getparfloat("d1", &mut d1);
    getparfloat("d2", &mut d2);
    getparfloat("f1", &mut f1);
    getparfloat("f2", &mut f2);
    checkpars();
    unsafe {
        let fp = fdopen(0, c"r".as_ptr());
        let m = readModel(fp);
        // The original dereferences the returned null model on empty input.
        if m.is_null() { signal(11, 0); raise(11); return; }
        let (xmin, xmax, zmin, zmax) = ((*m).ymin, (*m).ymax, (*m).xmin, (*m).xmax);
        if n1 < 0 || n2 < 0 {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        let count = (n1 as usize).checked_mul(n2 as usize)
            .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<f32>())
            .unwrap_or_else(|| err!("ealloc.c: malloc failed (Cannot allocate memory)"));
        let mut s = Vec::new();
        if s.try_reserve_exact(count).is_err() {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        s.resize(count, 0.0f32);
        let mut t = std::ptr::null_mut();
        let mut x = f2;
        for ix in 0..n2 {
            if ((x - xmin).abs() as f64) < 0.01 * d2 as f64 { x = xmin; }
            if ((x - xmax).abs() as f64) < 0.01 * d2 as f64 { x = xmax; }
            if x < xmin || x > xmax { err!("x=%g must be between xmin=%g and xmax=%g", x, xmin, xmax); }
            let mut z = f1;
            for iz in 0..n1 {
                if ((z - zmin).abs() as f64) < 0.01 * d1 as f64 { z = zmin; }
                if ((z - zmax).abs() as f64) < 0.01 * d1 as f64 { z = zmax; }
                if z < zmin || z > zmax { err!("z=%g must be between zmin=%g and zmax=%g", z, zmin, zmax); }
                t = insideTriInModel(m, t, z, x);
                let ta = &*((*t).fa.cast::<TriAttributes>());
                s[(ix * n1 + iz) as usize] = ta.s00 + x * ta.dsdx + z * ta.dsdz;
                z += d1;
            }
            x += d2;
        }
        fwrite(&s, &stdout());
    }
    exit(CWP_Exit());
}
