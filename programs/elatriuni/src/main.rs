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

//! ELATRIUNI - sample triangulated elastic models on a uniform grid.
//! Port of Trielas/main/elatriuni.c. Credits: Andreas Rueger and Dave Hale,
//! Colorado School of Mines, 1995 and 1991.

use oxwave_core::prelude::*;
use oxwave_lib_tri::*;
use std::ffi::c_void;
extern "C" {
    fn fdopen(fd: i32, mode: *const i8) -> *mut c_void;
    fn dup(fd: i32) -> i32;
    fn raise(sig: i32) -> i32;
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " ELATRIUNI - convert TRIangulated ELAstic models to UNIformly sampled models",
    "\t\t\t\t\t\t\t\t\t",
    "  elatriuni <modelfile nx= nz= [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Required Parameters:\t\t\t\t\t\t\t",
    "nx                     number of x samples\t\t\t\t",
    "nz                     number of z samples\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Optional Parameters:\t\t\t\t\t\t\t",
    "dx=1.0                 x sampling interval    \t\t\t    \t",
    "dz=1.0                 z sampling interval\t\t\t\t",
    "fx=0.0                 first x sampled\t\t\t\t\t",
    "fz=0.0                 first z sampled\t\t\t\t\t",
    "a1111file=a1111.bin    bin-file to store a1111 components \t\t",
    "a3333file=a3333.bin    bin-file to store a3333 components \t\t",
    "a1133file=a1133.bin    bin-file to store a1133 components \t\t",
    "a1313file=a1313.bin    bin-file to store a1313 components \t\t",
    "a1113file=a1113.bin    bin-file to store a1113 components \t\t",
    "a3313file=a3313.bin    bin-file to store a3313 components \t\t",
    "a1212file=a1212.bin    bin-file to store a1212 components \t\t",
    "a1223file=a1212.bin    bin-file to store a1223 components \t\t",
    "a2323file=a2323.bin    bin-file to store a2323 components \t\t",
    "rhofile=rho.bin        bin-file to store rho components \t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// Layout of FaceAttributes in Trielas/include/elastic.h.
#[repr(C)]
struct TriAttributes {
    a1111: f32, a3333: f32, a1133: f32, a1313: f32, a1113: f32,
    a3313: f32, a1212: f32, a2323: f32, a1223: f32, rho: f32, mindex: i32,
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut nx, mut nz) = (0, 0);
    if getparint("nx", &mut nx) == 0 { err!("must specify nx!"); }
    if getparint("nz", &mut nz) == 0 { err!("must specify nz!"); }
    let (mut dx, mut dz, mut fx, mut fz) = (1.0f32, 1.0f32, 0.0f32, 0.0f32);
    getparfloat("dx", &mut dx);
    getparfloat("dz", &mut dz);
    getparfloat("fx", &mut fx);
    getparfloat("fz", &mut fz);
    let names = ["a1111file", "a3333file", "a1133file", "a1313file", "a1113file", "a3313file", "a1212file", "a2323file", "a1223file", "rhofile"];
    let defaults = ["a1111.bin", "a3333.bin", "a1133.bin", "a1313.bin", "a1113.bin", "a3313.bin", "a1212.bin", "a2323.bin", "a1223.bin", "rho.bin"];
    let mut files = Vec::new();
    for i in 0..names.len() {
        let mut name = String::new();
        if getparstring(names[i], &mut name) == 0 { name = defaults[i].to_string(); }
        files.push(efopen(&name, "w"));
    }
    checkpars();
    unsafe {
        // readModel uses libc FILE*, while core owns the stdin descriptor.
        let fp = fdopen(dup(fileno(&stdin())), c"r".as_ptr());
        let m = readModel(fp);
        if m.is_null() { raise(11); exit(1); }
        let (xmin, xmax, zmin, zmax) = ((*m).ymin, (*m).ymax, (*m).xmin, (*m).xmax);
        // alloc2 allocates the row pointers before the contiguous float data.
        if nx < 0 || nz < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
        let n = (nx as usize).checked_mul(nz as usize)
            .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<f32>())
            .unwrap_or_else(|| err!("ealloc.c: malloc failed (Cannot allocate memory)"));
        let mut arrays = Vec::with_capacity(10);
        for _ in 0..10 {
            let mut a = Vec::new();
            if a.try_reserve_exact(n).is_err() { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
            a.resize(n, 0.0f32);
            arrays.push(a);
        }
        let mut t = std::ptr::null_mut();
        let mut x = fx;
        for ix in 0..nx {
            if x < xmin || x > xmax { err!("x=%g must be between xmin=%g and xmax=%g", x, xmin, xmax); }
            let mut z = fz;
            for iz in 0..nz {
                if z < zmin || z > zmax { err!("z=%g must be between zmin=%g and zmax=%g", z, zmin, zmax); }
                t = insideTriInModel(m, t, z, x);
                let ta = &*((*t).fa as *const TriAttributes);
                let values = [ta.a1111, ta.a3333, ta.a1133, ta.a1313, ta.a1113, ta.a3313, ta.a1212, ta.a2323, ta.a1223, ta.rho];
                for j in 0..10 { arrays[j][ix as usize * nz as usize + iz as usize] = values[j]; }
                z += dz;
            }
            x += dx;
        }
        for j in 0..10 { fwrite(&arrays[j], &files[j]); }
    }
    exit(1);
}
