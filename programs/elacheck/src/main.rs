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

//! ELACHECK - get elastic coefficients of model.
//! Port of Trielas/main/elacheck.c.
//! Credits: Andreas Rueger, Colorado School of Mines, 01/20/94.

use oxwave_core::prelude::*;
use oxwave_lib_tri::*;
use std::ffi::c_void;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " ELACHECK - get elastic coefficients of model  \t   \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " elacheck file= (required modelfile)\t\t\t\t\t",
    "                                                                       ",
    " ____ interactive program _____________                                ",
    " \t\t\t\t\t\t\t\t\t",
];

#[repr(C)]
struct FaceAttributes {
    a1111: f32, a3333: f32, a1133: f32, a1313: f32,
    a1113: f32, a3313: f32, a1212: f32, a2323: f32,
    a1223: f32, rho: f32, mindex: i32,
}

unsafe extern "C" {
    fn dup(fd: i32) -> i32;
    fn fdopen(fd: i32, mode: *const i8) -> *mut c_void;
    fn scanf(format: *const i8, ...) -> i32;
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut mfile = String::new();
    if getparstring("file", &mut mfile) == 0 {
        err!("ERROR: No modelfile defined");
    }
    let mfp = efopen(&mfile, "r");
    checkpars();

    let (mut a1111, mut a3333, mut a1133, mut a1313, mut a1113, mut a3313) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let (mut a1212, mut a1223, mut a2323, mut rho, mut v_p, mut v_s) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
    // readModel consumes a libc FILE pointer; duplicate the core stream's descriptor.
    let m = unsafe { readModel(fdopen(dup(fileno(&mfp)), c"r".as_ptr())) };
    let fp = stderr();
    fprintf!(fp, " **********************************************\n");
    fprintf!(fp, " *********check stiffness coefficients*********\n");
    fprintf!(fp, " **********************************************\n\n");
    let (mut x, mut z) = (0.0f32, 0.0f32);
    loop {
        let mut want = 1i32;
        fprintf!(fp, " Enter x-coordinate \n");
        unsafe { scanf(c"%f".as_ptr(), &mut x); }
        fprintf!(fp, " Enter z-coordinate \n");
        unsafe { scanf(c"%f".as_ptr(), &mut z); }
        unsafe {
            if x > (*m).ymax || x < (*m).ymin || z > (*m).xmax || z < (*m).xmin {
                fprintf!(fp, " Coordinate outside of model \n\n");
                continue;
            }
            let t = insideTriInModel(m, std::ptr::null_mut(), z, x);
            let fa = (*t).fa.cast::<FaceAttributes>();
            if !fa.is_null() {
                a1111 = (*fa).a1111; a3333 = (*fa).a3333; a1133 = (*fa).a1133;
                a1313 = (*fa).a1313; a1113 = (*fa).a1113; a3313 = (*fa).a3313;
                v_p = (a3333 as f64).sqrt() as f32;
                v_s = (a1313 as f64).sqrt() as f32;
                a1212 = (*fa).a1212; a1223 = (*fa).a1223; a2323 = (*fa).a2323; rho = (*fa).rho;
            }
        }
        fprintf!(fp, "\n ----Coordinates x=%g \t z=%g ----\n\n", x, z);
        fprintf!(fp, "\t\t a1111=%g \n", a1111);
        fprintf!(fp, "\t\t a3333=%g \n", a3333);
        fprintf!(fp, "\t\t a1133=%g \n", a1133);
        fprintf!(fp, "\t\t a1313=%g \n", a1313);
        fprintf!(fp, "\t\t a1113=%g \n", a1113);
        fprintf!(fp, "\t\t a3313=%g \n", a3313);
        fprintf!(fp, "\t\t a1212=%g \n", a1212);
        fprintf!(fp, "\t\t a1223=%g \n", a1223);
        fprintf!(fp, "\t\t a2323=%g \n", a2323);
        fprintf!(fp, "\t\t rho=%g \n", rho);
        fprintf!(fp, "\t\t   v_p=%g \n", v_p);
        fprintf!(fp, "\t\t   v_s=%g \n\n", v_s);
        fprintf!(fp, " More checks type: 1 if yes, 0 if no \n");
        unsafe { scanf(c"%i".as_ptr(), &mut want); }
        if want == 0 { break; }
    }
    exit(1);
}
