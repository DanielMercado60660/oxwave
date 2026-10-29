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

//! TRIRAY - dynamic ray tracing for a triangulated sloth model.
//! Port of tri/main/triray.c (Test Release 1.1, 2011/11/21).
//! Credits: Dave Hale, Colorado School of Mines, 02/16/91.
//! Modified by Andreas Rueger, Colorado School of Mines, 08/12/93:
//! Fresnel volumes, source-on-edge handling, reflection sequences, amplitude and attenuation.

#![allow(non_snake_case)]
use oxwave_core::prelude::*;
use oxwave_lib_tri::*;
use std::ffi::c_void;
use std::ptr::null_mut;
mod ray;
use ray::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t\t",
    " TRIRAY - dynamic RAY tracing for a TRIangulated sloth model\t\t",
    "\t\t\t\t\t\t\t\t\t\t",
    "  triray <modelfile >rayends [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " xs=(max-min)/2 x coordinate of source (default is halfway across model)",
    " zs=min\t z coordinate of source (default is at top of model)\t",
    " nangle=101\tnumber of takeoff angles\t\t\t\t",
    " fangle=-45\tfirst takeoff angle (in degrees)\t\t\t",
    " langle=45\tlast takeoff angle (in degrees)\t\t\t\t",
    " rayfile=\tfile of ray x,z coordinates of ray-edge intersections\t",
    " nxz=101\tnumber of (x,z) in optional rayfile (see notes below)\t",
    " wavefile=\tfile of ray x,z coordinates uniformly sampled in time\t",
    " nt=101\tnumber of (x,z) in optional wavefile (see notes below)\t",
    " infofile=\tASCII-file to store useful information \t\t",
    " fresnelfile=  used if you want to plot the fresnel volumes.\t\t",
    "\t\tdefault is <fresnelfile.bin> \t\t\t\t",
    " outparfile=\tcontains parameters for the plotting software.\t\t",
    "\t\tdefault is <outpar> \t\t\t\t\t",
    " krecord=\tif specified, only rays incident at interface with index",
    "\t\tkrecord are displayed and stored\t\t\t",
    " prim=\t   =1, only single-reflected rays are plotted\t\t\t",
    "\t\t=0, only direct hits are displayed  \t\t\t",
    " ffreq=-1\tFresnelVolume frequency \t\t\t\t",
    " refseq=1,0,0  index of reflector followed by sequence of reflection (1)",
    "\t\ttransmission(0) or ray stops(-1).\t\t\t",
    "\t\tThe default rayend is at the model boundary.\t\t",
    "\t\tNOTE:refseq must be defined for each reflector\t\t",
    " NOTES:\t\t\t\t\t\t\t\t",
    " The rayends file contains ray parameters for the locations at which\t",
    " the rays terminate.  \t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The rayfile is useful for making plots of ray paths.\t\t\t",
    " nxz should be larger than twice the number of triangles intersected\t",
    " by the rays.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The wavefile is useful for making plots of wavefronts.\t\t",
    " The time sampling interval in the wavefile is tmax/(nt-1),\t\t",
    " where tmax is the maximum time for all rays.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The infofile is useful for collecting information along the\t\t",
    " individual rays. The fresnelfile contains data used to plot \t\t",
    " the Fresnel volumes. The outparfile stores information used \t\t",
    " for the plotting software.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

#[repr(C)]
#[derive(Clone, Copy)]
struct RayStep {
    sigma: f32,
    x: f32,
    z: f32,
    px: f32,
    pz: f32,
    t: f32,
    q1: f32,
    p1: f32,
    q2: f32,
    p2: f32,
    kmah: i32,
    nref: i32,
    eu: *mut EdgeUse,
    atten: f32,
    ampli: f32,
    ampliphase: f32,
    f: *mut Face,
    rsNext: *mut RayStep,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct RefCheck {
    k: i32,
    nhits: i32,
    hitseq: *mut i32,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct FaceAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct EdgeAttributes {
    k: i32,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct EdgeUseAttributes {
    tx: f32,
    tz: f32,
    c: f32,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct RayEnd {
    sigma: f32,
    x: f32,
    z: f32,
    px: f32,
    pz: f32,
    t: f32,
    q1: f32,
    p1: f32,
    q2: f32,
    p2: f32,
    kmah: i32,
    nref: i32,
    sb: f32,
    se: f32,
    dsdxe: f32,
    dsdze: f32,
    ab: f32,
    kend: i32,
    ampli: f32,
    ampliphase: f32,
    atten: f32,
    dangle: f32,
}

// RayEnd has only 32-bit scalar fields and no padding: native C binary output.
unsafe impl oxwave_core::cio::Pod for RayEnd {}
const _: () = assert!(std::mem::size_of::<RayEnd>() == 88);

// The tri library reads the native SU model using a libc FILE. Core owns
// parameter, diagnostic and output streams; this duplicate only serves readModel.
extern "C" {
    fn fdopen(fd: i32, mode: *const i8) -> *mut c_void;
    fn dup(fd: i32) -> i32;
    fn fclose(fp: *mut c_void) -> i32;
    fn raise(signal: i32) -> i32;
    fn signal(signal: i32, handler: usize) -> usize;
}

// Restore the C signal disposition before reproducing a null-model fault:
// Rust's stack-overflow handler otherwise consumes a raised SIGSEGV once.
unsafe fn c_segv() -> ! {
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(139);
}

unsafe fn alloc_float(n: usize) -> *mut f32 {
    let mut values = Vec::new();
    if values.try_reserve_exact(n).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    values.resize(n, 0.0f32);
    Box::leak(values.into_boxed_slice()).as_mut_ptr()
}

fn open_optional(name: &str) -> *mut CFile {
    let mut path = String::new();
    if getparstring(name, &mut path) != 0 {
        Box::into_raw(Box::new(efopen(&path, "w")))
    } else {
        null_mut()
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    unsafe {
        run();
    }
    exit(CWP_Exit());
}

unsafe fn run() {
    unsafe {
        let input = fdopen(dup(0), c"r".as_ptr());
        let m = readModel(input);
        fclose(input);
        let (mut xs, mut zs, mut fangle, mut langle, mut ffreq) =
            (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
        let (mut nangle, mut nxz, mut nt, mut krecord, mut prim) = (0, 0, 0, 0, 0);
        if getparfloat("xs", &mut xs) == 0 {
            if m.is_null() {
                c_segv();
            }
            xs = (0.5 * ((*m).ymin + (*m).ymax) as f64) as f32;
        }
        if getparfloat("zs", &mut zs) == 0 {
            if m.is_null() {
                c_segv();
            }
            zs = (*m).xmin;
        }
        if getparint("nangle", &mut nangle) == 0 {
            nangle = 101;
        }
        if getparfloat("fangle", &mut fangle) == 0 {
            fangle = -45.0;
        }
        if getparfloat("langle", &mut langle) == 0 {
            langle = 45.0;
        }
        if getparint("nxz", &mut nxz) == 0 {
            nxz = 101;
        }
        if getparint("nt", &mut nt) == 0 {
            nt = 101;
        }
        if getparint("krecord", &mut krecord) == 0 {
            krecord = i32::MAX;
        }
        if getparint("prim", &mut prim) == 0 {
            prim = i32::MAX;
        }
        if getparfloat("ffreq", &mut ffreq) == 0 {
            ffreq = -1.0;
        }

        let rayfp = open_optional("rayfile");
        let wavefp = open_optional("wavefile");
        let infofp = open_optional("infofile");
        let mut fresnelfile = "fresnelfile.bin".to_string();
        let mut outparfile = "outpar".to_string();
        getparstring("fresnelfile", &mut fresnelfile);
        let fresnelfp = if ffreq > 0.0 {
            Box::into_raw(Box::new(efopen(&fresnelfile, "w")))
        } else {
            null_mut()
        };
        getparstring("outparfile", &mut outparfile);
        let outparfp = Box::into_raw(Box::new(efopen(&outparfile, "w")));

        let kk = countparname("refseq");
        if prim != i32::MAX && krecord == i32::MAX {
            err!("krecord= must be specified with prim=\n");
        }
        if prim != i32::MAX && kk == 0 {
            err!("must define target reflector(s) via refseq= \n");
        }

        let mut rc: Vec<RefCheck> = Vec::with_capacity(kk as usize);
        let mut hitseqa: Vec<Vec<i32>> = Vec::with_capacity(kk as usize);
        let mut nrefseq = vec![0i32; kk as usize];
        for i in 0..kk as usize {
            let nseq = countnparval(i as i32 + 1, "refseq");
            nrefseq[i] = nseq;
            let mut seq = vec![0i32; (nseq + 2) as usize];
            getnparint(i as i32 + 1, "refseq", &mut seq[..]);
            rc.push(RefCheck {
                k: seq[0],
                nhits: 0,
                hitseq: seq.as_mut_ptr().add(1),
            });
            hitseqa.push(seq);
        }
        if !infofp.is_null() {
            fprintf!(&*infofp,"THIS FILE CONTAINS RAY TRACING INFORMATION\n\n*********************************************\ndefined reflection/transmission sequences:\n*********************************************\n");
            for i in 0..kk as usize {
                fprintf!(&*infofp, "interface %d:", rc[i].k);
                for ii in 0..nrefseq[i] - 1 {
                    fprintf!(
                        &*infofp,
                        " %d%s",
                        *rc[i].hitseq.offset(ii as isize),
                        if ii < nrefseq[i] - 2 { "," } else { "" }
                    );
                }
                fprintf!(&*infofp, "\n");
            }
            fprintf!(
                &*infofp,
                "Interfaces without defined refseq are transmitting.\n"
            );
        }

        fangle = (fangle as f64 * (std::f64::consts::PI / 180.0)) as f32;
        langle = (langle as f64 * (std::f64::consts::PI / 180.0)) as f32;
        let dangle = if nangle > 1 {
            (langle - fangle) / (nangle - 1) as f32
        } else {
            0.0
        };
        let mut rsp: Vec<*mut RayStep> = Vec::new();
        if rsp.try_reserve_exact(nangle as usize).is_err() {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        rsp.resize(nangle as usize, null_mut());
        let mut re: Vec<RayEnd> = vec![std::mem::zeroed(); nangle as usize];
        shootRays(
            m,
            xs,
            zs,
            nangle,
            dangle,
            fangle,
            rc.as_mut_ptr(),
            kk,
            rsp.as_mut_ptr(),
            re.as_mut_ptr(),
            infofp,
        );
        if !fresnelfp.is_null() && !rayfp.is_null() {
            writeFresnel(
                m,
                rayfp,
                infofp,
                nxz,
                nangle,
                rsp.as_mut_ptr(),
                re.as_mut_ptr(),
                krecord,
                ffreq,
                prim,
                fresnelfp,
                outparfp,
            );
        } else if !rayfp.is_null() {
            writeRays(
                m,
                rayfp,
                nxz,
                nangle,
                rsp.as_mut_ptr(),
                re.as_mut_ptr(),
                krecord,
                prim,
                infofp,
                outparfp,
            );
        }
        if !wavefp.is_null() {
            writeWaves(wavefp, nt, nangle, rsp.as_mut_ptr());
        }
        if efwrite(&re, &stdout()) != nangle as usize {
            err!("Error writing ray ends to stdout!\n");
        }
    }
}
