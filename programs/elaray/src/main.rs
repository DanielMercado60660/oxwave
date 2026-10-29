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

//! ELARAY - ray tracing for elastic anisotropic models.
//! Port of Trielas/main/elaray.c (Test Release 1.2).
//! Author: Andreas Rueger, Colorado School of Mines, 01/02/95.
//! Based on gbray.c (Andreas Rueger, 08/12/93) and sdray.c
//! (Dave Hale, CSM, 02/26/91).
//! Original additional notices: Copyright (c) Colorado School of Mines,
//! 1999 and 2011. All rights reserved.

#![allow(
    non_snake_case,
    unused_mut,
    unused_variables,
    unused_assignments,
    unused_parens
)]
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
use oxwave_lib_tri::*;
use oxwave_lib_trielas::*;
mod helpers;
use helpers::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " ELARAY - ray tracing for elastic anisotropic models",
    "                                              \t\t\t\t",
    " elaray <modelfile >rayends [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " xs=(max-min)/2 x coordinate of source (default is halfway across model)",
    " zs=min         z coordinate of source (default is at top of model)\t",
    " nangle=101     number of takeoff angles\t\t\t\t",
    " fangle=-45     first takeoff angle (in degrees)\t\t\t",
    " langle=45      last takeoff angle (in degrees)\t\t\t",
    " nxz=101        number of (x,z) in optional rayfile (see notes below)\t",
    " mode=0         shoot P-rays\t\t\t\t\t\t",
    "\t      =1 shoot SV-rays\t\t\t\t\t\t",
    "\t      =2 shoot SH-rays\t\t\t\t\t\t",
    " prim        =1 only reflected rays are plotted \t\t",
    "             =0 only direct hits are displayed  \t\t\t",
    " refseq=1,0,0   index of reflector followed by sequence of:\t\t",
    "\t\t transmission(0)\t\t\t\t\t",
    "\t\t reflection (1)\t\t\t\t\t\t",
    "\t\t transmission with mode conversion (2)\t\t\t",
    "\t\t reflection with mode conversion (3)\t\t\t",
    "                ray stops(-1).\t\t\t\t\t\t",
    " krecord        if specified, only rays incident at interface with index",
    "                krecord are displayed and stored\t\t\t",
    " f0=1\t         force impact strenght\t\t\t\t\t",
    " fdip=0         force dip with respect to vertical\t\t\t",
    " fazi=0         force azimuth with respect to positive x-axis \t\t",
    " reftrans=0\t =1 include reflec/transm coeff(currently only for P)\t",
    " rayfile        file of ray x,z coordinates of ray-edge intersections\t",
    " wavefile       file of ray x,z coordinates uniformly sampled in time\t",
    " nt=\t\t number of (x,z) in optional wavefile (see notes below)\t",
    " tw=\t\t traveltime associated with wavefront (alternative to nt)",
    " infofile       ASCII-file to store useful information \t\t",
    " outparfile     contains parameters for the plotting software. \t",
    "                default is <outpar> \t\t\t\t\t",
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
    " where tmax is the maximum time for all rays. Alternatively, ",
    " one wavefront at time tw is plotted.\t",
    "\t\t\t\t\t\t\t\t\t",
    " The infofile is useful for collecting information along the\t\t",
    " individual rays. \t\t\t\t\t\t\t",
    " The outparfile stores information used for the plotting software\t",
    "\t\t\t\t\t\t\t\t\t",
];

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct FaceAttributes {
    a1111: f32,
    a3333: f32,
    a1133: f32,
    a1313: f32,
    a1113: f32,
    a3313: f32,
    a1212: f32,
    a2323: f32,
    a1223: f32,
    rho: f32,
    mindex: i32,
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct EdgeAttributes {
    k: i32,
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct EdgeUseAttributes {
    tx: f32,
    tz: f32,
    c: f32,
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
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
    vgx: f32,
    vgz: f32,
    g11: f32,
    g13: f32,
    g33: f32,
    kmah: i32,
    nref: i32,
    mode: i32,
    eu: *mut EdgeUse,
    atten: f32,
    ampli: f32,
    ampliphase: f32,
    f: *mut Face,
    rsNext: *mut RayStep,
    rho: f32,
}
#[repr(C)]
#[derive(Clone, Copy, Default)]
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
    vgx: f32,
    vgz: f32,
    ab: f32,
    kend: i32,
    ampli: f32,
    ampliphase: f32,
    dangle: f32,
    mode: i32,
    g11: f32,
    g13: f32,
    g33: f32,
    nameref: i32,
    num: i32,
    rhob: f32,
    rhoe: f32,
}
// repr(C), with only adjacent four-byte floats and integers and no padding.
unsafe impl oxwave_core::cio::Pod for RayEnd {}
#[derive(Clone, Copy, Default)]
struct RefCheck {
    k: i32,
    nhits: i32,
    hitseq: *mut i32,
}

// The topology library accepts libc FILE*, whereas core owns CFile streams.
// A duplicate descriptor connects stdin to readModel without competing buffers.
extern "C" {
    fn dup(fd: i32) -> i32;
    fn fdopen(fd: i32, mode: *const std::ffi::c_char) -> *mut std::ffi::c_void;
    fn fclose(fp: *mut std::ffi::c_void) -> i32;
    fn raise(sig: i32) -> i32;
    fn signal(sig: i32, handler: usize) -> usize;
    fn pipe(fds: *mut i32) -> i32;
}

unsafe fn segv() -> ! {
    signal(11, 0);
    raise(11);
    exit(139);
}

fn optional_file(name: &str) -> *mut CFile {
    let mut path = String::new();
    if getparstring(name, &mut path) != 0 {
        Box::into_raw(Box::new(efopen(&path, "w")))
    } else {
        std::ptr::null_mut()
    }
}

// The coefficient library's private copy of findqPqSV suppresses the C
// fprintf(NULL, ...) on root failure. Capture its diagnostic to reproduce
// that signal here, without duplicating the library's numerical routines.
unsafe fn rt_ani_real(
    gz: f32,
    gx: f32,
    a1111i: f32,
    a3333i: f32,
    a1133i: f32,
    a1313i: f32,
    a1113i: f32,
    a3313i: f32,
    rhoi: f32,
    a1111t: f32,
    a3333t: f32,
    a1133t: f32,
    a1313t: f32,
    a1113t: f32,
    a3313t: f32,
    rhot: f32,
    pxi: f32,
    pzi: f32,
    g11i: f32,
    g13i: f32,
    g33i: f32,
    pxg: f32,
    pzg: f32,
    g11g: f32,
    g13g: f32,
    g33g: f32,
    modei: i32,
    modeg: i32,
    rort: i32,
    coeff: &mut f32,
    ifp: Option<&CFile>,
) -> i32 {
    use std::io::Read;
    use std::os::fd::FromRawFd;
    // A coefficient call emits at most a few hundred diagnostic bytes.
    // A pipe captures them without adding a temporary-file dependency.
    let log = if ifp.is_none() {
        let mut fds = [0; 2];
        if pipe(fds.as_mut_ptr()) != 0 {
            syserr!("pipe failed");
        }
        Some((
            CFile::from_file(std::fs::File::from_raw_fd(fds[1])),
            std::fs::File::from_raw_fd(fds[0]),
        ))
    } else {
        None
    };
    let result = oxwave_lib_trielas::rt_ani_real(
        gz,
        gx,
        a1111i,
        a3333i,
        a1133i,
        a1313i,
        a1113i,
        a3313i,
        rhoi,
        a1111t,
        a3333t,
        a1133t,
        a1313t,
        a1113t,
        a3313t,
        rhot,
        pxi,
        pzi,
        g11i,
        g13i,
        g33i,
        pxg,
        pzg,
        g11g,
        g13g,
        g33g,
        modei,
        modeg,
        rort,
        coeff,
        ifp.or(log.as_ref().map(|(fp, _)| fp)),
    );
    if let Some((log, mut reader)) = log {
        efclose(&log);
        let mut bytes = Vec::new();
        let _ = reader.read_to_end(&mut bytes);
        if bytes
            .windows(b" rootfinder not sucessful !!! \n".len())
            .any(|w| w == b" rootfinder not sucessful !!! \n")
        {
            segv();
        }
    }
    result
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    unsafe {
        run();
    }
}

unsafe fn run() {
    let fp = fdopen(dup(0), c"r".as_ptr());
    let m = readModel(fp);
    fclose(fp);
    let mut mode = 0;
    getparint("mode", &mut mode);
    if mode != 0 && mode != 1 && mode != 2 {
        err!("ERROR Wrong raymode defined \n");
    }
    // C dereferences a null model for empty input, terminating with SIGSEGV.
    if m.is_null() {
        segv();
    }

    let mut xs = 0.0;
    if getparfloat("xs", &mut xs) == 0 {
        xs = (0.5 * ((*m).ymin + (*m).ymax) as f64) as f32;
    }
    let mut zs = 0.0;
    if getparfloat("zs", &mut zs) == 0 {
        zs = (*m).xmin;
    }
    let mut nangle = 101;
    getparint("nangle", &mut nangle);
    let mut reftrans = 0;
    getparint("reftrans", &mut reftrans);
    let mut fangle = -45.0f32;
    getparfloat("fangle", &mut fangle);
    let mut langle = 45.0f32;
    getparfloat("langle", &mut langle);
    let mut f0 = 1.0;
    getparfloat("f0", &mut f0);
    let mut fdip = 0.0;
    getparfloat("fdip", &mut fdip);
    let mut fazi = 0.0;
    getparfloat("fazi", &mut fazi);
    let rayfp = optional_file("rayfile");
    let ifp = optional_file("infofile");
    let mut nxz = 101;
    getparint("nxz", &mut nxz);
    let wavefp = optional_file("wavefile");
    let mut nt = i32::MAX;
    getparint("nt", &mut nt);
    let mut tw = f32::MAX;
    getparfloat("tw", &mut tw);
    let mut krecord = i32::MAX;
    getparint("krecord", &mut krecord);
    let mut prim = i32::MAX;
    getparint("prim", &mut prim);
    let mut outparfp = optional_file("outparfile");
    if outparfp.is_null() {
        outparfp = Box::into_raw(Box::new(efopen("outpar", "w")));
    }
    checkpars();

    let kk = countparname("refseq");
    if prim == 1 && kk == 0 {
        err!("\n first define target reflector \n ");
    }
    let mut rc = vec![RefCheck::default(); kk as usize];
    if reftrans != 0 {
        warn!("\n Refl/Transm Coeff. need more work and testing !!!!\n");
    }
    let mut hitseqa = vec![Vec::<i32>::new(); if kk > 0 { kk as usize } else { 1 }];
    let mut nrefseq = vec![0; hitseqa.len()];
    if !ifp.is_null() {
        fprintf!(&*ifp, "\n THIS FILE CONTAINS RAY TRACING INFORMATION \n");
        fprintf!(&*ifp, "\n *********************************************\n");
        fprintf!(&*ifp, " defined reflection/transmission sequences \n");
        fprintf!(&*ifp, " *********************************************\n");
    }
    for i in 0..kk as usize {
        let nseq = countnparval(i as i32 + 1, "refseq");
        nrefseq[i] = nseq - 1;
        // malloc's minimum allocation in the reference has room for four
        // integers, including the zero-filled tail read by an empty sequence.
        hitseqa[i] = vec![0; (nseq as usize).max(4)];
        getnparint(i as i32 + 1, "refseq", &mut hitseqa[i]);
        rc[i].k = hitseqa[i][0];
        rc[i].nhits = 0;
        rc[i].hitseq = hitseqa[i].as_mut_ptr().add(1);
        if !ifp.is_null() {
            fprintf!(&*ifp, "\n interface %i :", rc[i].k);
            for ii in 0..nseq - 1 {
                fprintf!(&*ifp, "%i ", *rc[i].hitseq.offset(ii as isize));
            }
        }
    }
    if !ifp.is_null() {
        fprintf!(
            &*ifp,
            "\n interfaces without defined refseq are transmitting \n"
        );
    }
    fangle = (fangle as f64 * (PI / 180.0)) as f32;
    langle = (langle as f64 * (PI / 180.0)) as f32;
    fdip = (fdip as f64 * (PI / 180.0)) as f32;
    fazi = (fazi as f64 * (PI / 180.0)) as f32;
    let dangle = if nangle == 1 {
        0.0
    } else {
        (langle - fangle) / nangle.wrapping_sub(1) as f32
    };
    // C alloc1 fails for a negative count, but both tracing loops still run
    // zero times. Its final efwrite then reports zero items of that count.
    let count = if nangle > 0 { nangle as usize } else { 0 };
    let mut rsp = vec![std::ptr::null_mut(); count];
    let mut re = vec![RayEnd::default(); count];
    shootRaysAni(
        m,
        xs,
        zs,
        nangle,
        reftrans,
        dangle,
        fangle,
        rc.as_mut_ptr(),
        kk,
        rsp.as_mut_ptr(),
        re.as_mut_ptr(),
        nrefseq.as_mut_ptr(),
        ifp,
        mode,
        fdip,
        f0,
        fazi,
    );
    if !rayfp.is_null() {
        writeRays(
            m,
            rayfp,
            nxz,
            nangle,
            rsp.as_mut_ptr(),
            re.as_mut_ptr(),
            krecord,
            prim,
            ifp,
            outparfp,
        );
    }
    if !wavefp.is_null() {
        if nt == i32::MAX && tw != f32::MAX {
            writeWaves2(wavefp, tw, nangle, rsp.as_mut_ptr());
        } else if nt != i32::MAX && tw == f32::MAX {
            writeWaves(wavefp, nt, nangle, rsp.as_mut_ptr());
        } else {
            fprintf!(stderr(), "\n define either nt or tw \n");
        }
    }
    if nangle < 0 {
        err!(
            "%s: efwrite: fwrite only %d items of %d",
            "subcalls.c",
            0,
            nangle
        );
    }
    if efwrite(&re, &stdout()) != count {
        err!("Error writing ray ends to stdout!\n");
    }
    exit(1);
}
