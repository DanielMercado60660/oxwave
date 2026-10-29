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

//! NORMRAY - dynamic ray tracing for normal incidence rays in a sloth model.
//! Port of tri/main/normray.c (2011/11/21).
//! Author: Dave Hale, Colorado School of Mines, 02/16/91.
//! Modified: Andreas Rueger, 08/12/93; Boyi Ou, 4/14/95.

#![allow(non_snake_case, unused, unused_parens)]
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, max};
use oxwave_lib_tri::*;
use std::ffi::c_void;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                              \t\t\t\t",
    " NORMRAY - dynamic ray tracing for normal incidence rays in a sloth model",
    "                                              \t\t\t\t",
    "    normray <modelfile >rayends [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " caustic\t 0: show all rays 1: show only caustic rays\t\t",
    " nonsurface\t 0: show rays which reach surface 1: show all rays      ",
    " surface\t 0: shot ray from subsurface 1: from surface               ",
    " nrays \t number of location to shoot rays                       ",
    " dangle \t increment of ray angle for one location                ",
    " nangle \t number of rays shot from one location                  ",
    " ashift \t shift first taking off angle                           ",
    " xs1 \t         x of shooting location                                 ",
    " zs1 \t         z of shooting location                                 ",
    " nangle=101     number of takeoff angles\t\t\t\t",
    " fangle=-45     first takeoff angle (in degrees)\t\t\t",
    " rayfile        file of ray x,z coordinates of ray-edge intersections\t",
    " nxz=101        number of (x,z) in optional rayfile (see notes below)\t",
    " wavefile       file of ray x,z coordinates uniformly sampled in time\t",
    " nt=101         number of (x,z) in optional wavefile (see notes below)\t",
    " infofile       ASCII-file to store useful information \t\t",
    " fresnelfile    used if you want to plot the fresnel volumes. \t\t",
    "                default is <fresnelfile.bin> \t\t\t\t",
    " outparfile     contains parameters for the plotting software. \t",
    "                default is <outpar> \t\t\t\t\t",
    " krecord        if specified, only rays incident at interface with index",
    "                krecord are displayed and stored\t\t\t",
    " prim           =1, only single-reflected rays are plotted \t\t",
    "                =0, only direct hits are displayed  \t\t\t",
    " ffreq=-1       FresnelVolume frequency \t\t\t\t",
    " refseq=1,0,0   index of reflector followed by sequence of reflection (1)",
    "                transmission(0) or ray stops(-1).\t\t\t",
    "                The default rayend is at the model boundary.\t\t",
    "                NOTE:refseq must be defined for each reflector\t\t",
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
struct FaceAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}

#[repr(C)]
struct EdgeUseAttributes {
    tx: f32,
    tz: f32,
    c: f32,
}

#[repr(C)]
struct EdgeAttributes {
    k: i32,
}

#[repr(C)]
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

#[repr(C)]
struct RefCheck {
    k: i32,
    nhits: i32,
    hitseq: *mut i32,
}

extern "C" {
    static mut __stdinp: *mut c_void;
    fn calloc(n: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn signal(sig: i32, handler: usize) -> usize;
    fn raise(sig: i32) -> i32;
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
    // readModel's public API takes a native FILE*. All other streams use core.
    let m = readModel(__stdinp);

    let mut nrays = 1;
    let mut caustic = 0;
    let mut surface = 0;
    let mut nonsurface = 0;
    let mut nangle = 1;
    let mut nxz = 101;
    let mut nt = 101;
    let mut krecord = i32::MAX;
    let mut prim = i32::MAX;
    let mut ffreq = -1.0f32;
    let mut dangle = 0.0f32;
    let mut ashift = 0.0f32;
    let mut xs1 = 0.0f32;
    let mut zs1 = 0.0f32;
    getparint("nrays", &mut nrays);
    getparint("caustic", &mut caustic);
    getparint("surface", &mut surface);
    getparint("nonsurface", &mut nonsurface);
    getparint("nangle", &mut nangle);
    getparint("nxz", &mut nxz);
    getparint("nt", &mut nt);
    getparint("krecord", &mut krecord);
    getparint("prim", &mut prim);
    getparfloat("ffreq", &mut ffreq);
    getparfloat("dangle", &mut dangle);
    getparfloat("ashift", &mut ashift);
    getparfloat("xs1", &mut xs1);
    getparfloat("zs1", &mut zs1);

    /* open files as necessary */
    let mut name = String::new();
    let rayfp = if getparstring("rayfile", &mut name) != 0 {
        Box::into_raw(Box::new(efopen(&name, "w")))
    } else {
        std::ptr::null_mut()
    };
    let wavefp = if getparstring("wavefile", &mut name) != 0 {
        Box::into_raw(Box::new(efopen(&name, "w")))
    } else {
        std::ptr::null_mut()
    };
    let infofp = if getparstring("infofile", &mut name) != 0 {
        Box::into_raw(Box::new(efopen(&name, "w")))
    } else {
        std::ptr::null_mut()
    };
    let mut fresnelfile = "fresnelfile.bin".to_string();
    getparstring("fresnelfile", &mut fresnelfile);
    let fresnelfp = if ffreq > 0.0 {
        Box::into_raw(Box::new(efopen(&fresnelfile, "w")))
    } else {
        std::ptr::null_mut()
    };
    let mut outparfile = "outpar".to_string();
    getparstring("outparfile", &mut outparfile);
    checkpars();
    let outparfp = efopen(&outparfile, "w");

    /* number of reflection sequences */
    let kk = countparname("refseq");
    if prim != i32::MAX && krecord == i32::MAX {
        err!("krecord= must be specified with prim=\n");
    }
    if prim != i32::MAX && kk == 0 {
        err!("must define target reflector(s) via refseq= \n");
    }

    /* allocate space for raychecks and hitseq */
    let mut rc: Vec<RefCheck> = (0..kk)
        .map(|_| RefCheck {
            k: 0,
            nhits: 0,
            hitseq: std::ptr::null_mut(),
        })
        .collect();
    let mut hitseqa = vec![Vec::<i32>::new(); kk as usize];
    let mut nrefseq = vec![0i32; kk as usize];
    let mut i = 0;
    while i < kk {
        let nseq = countnparval(i + 1, "refseq");
        nrefseq[i as usize] = nseq;
        // The C program reserves two extra entries after each sequence.
        hitseqa[i as usize] = vec![0i32; (nseq + 2) as usize];
        getnparint(i + 1, "refseq", &mut hitseqa[i as usize][..nseq as usize]);
        rc[i as usize].k = hitseqa[i as usize][0];
        rc[i as usize].nhits = 0;
        rc[i as usize].hitseq = hitseqa[i as usize].as_mut_ptr().add(1);
        i += 1;
    }

    /* print reflection sequence information if requested */
    if !infofp.is_null() {
        fprintf!(&*infofp, "THIS FILE CONTAINS RAY TRACING INFORMATION\n\n*********************************************\ndefined reflection/transmission sequences:\n*********************************************\n");
        i = 0;
        while i < kk {
            fprintf!(&*infofp, "interface %d:", rc[i as usize].k);
            let mut ii = 0;
            while ii < nrefseq[i as usize] - 1 {
                fprintf!(
                    &*infofp,
                    " %d%s",
                    *rc[i as usize].hitseq.offset(ii as isize),
                    if ii < nrefseq[i as usize] - 2 {
                        ","
                    } else {
                        ""
                    }
                );
                ii += 1;
            }
            fprintf!(&*infofp, "\n");
            i += 1;
        }
        fprintf!(
            &*infofp,
            "Interfaces without defined refseq are transmitting.\n"
        );
    }

    /* allocate space to xs, zs, as */
    if nrays < 0 {
        // C's unchecked alloc1float fails for a negative count, and fread
        // subsequently faults. Preserve SIGSEGV without a Rust capacity panic.
        signal(11, 0);
        raise(11);
    }
    let mut xs = vec![0.0f32; nrays as usize];
    let mut zs = vec![0.0f32; nrays as usize];
    // In C, the single-source assignment below writes as[i], where i==kk.
    // Reserve that otherwise out-of-bounds slot without changing as[0].
    let mut as_ = vec![0.0f32; (nrays as usize).max(kk as usize + 1)];
    /* open file of XInterface,ZInterface,AInterface */
    let xfp = fopen("XInterface", "r")
        .map(|f| Box::into_raw(Box::new(f)))
        .unwrap_or(std::ptr::null_mut());
    let zfp = fopen("ZInterface", "r")
        .map(|f| Box::into_raw(Box::new(f)))
        .unwrap_or(std::ptr::null_mut());
    let afp = if surface == 0 {
        fopen("AInterface", "r")
            .map(|f| Box::into_raw(Box::new(f)))
            .unwrap_or(std::ptr::null_mut())
    } else {
        std::ptr::null_mut()
    };
    fread(&mut xs, &*xfp);
    fread(&mut zs, &*zfp);
    if surface == 0 {
        fread(&mut as_[..nrays as usize], &*afp);
    } else {
        let mut ia = 0;
        while ia < nrays {
            as_[ia as usize] = 0.0;
            ia += 1;
        }
    }

    /* convert angles to radians and determine angle increment */
    dangle = ((dangle as f64) * (std::f64::consts::PI / 180.0)) as f32;
    ashift = ((ashift as f64) * (std::f64::consts::PI / 180.0)) as f32;
    if nrays == 1 {
        xs[0] = xs1;
        zs[0] = zs1;
        as_[i as usize] = 180.0;
    }

    let mut raycount = 0;
    let mut rsp: *mut *mut RayStep = std::ptr::null_mut();
    let mut re: *mut RayEnd = std::ptr::null_mut();
    i = 0;
    while i < nrays {
        /* allocate space for ray step lists and ray ends */
        rsp = calloc(nangle as usize, std::mem::size_of::<*mut RayStep>()) as *mut *mut RayStep;
        if rsp.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        re = calloc(nangle as usize, std::mem::size_of::<RayEnd>()) as *mut RayEnd;
        if re.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        let fangle = (as_[i as usize] as f64 * std::f64::consts::PI / 180.0
            - ((nangle / 2) as f32 * dangle) as f64
            - ashift as f64) as f32;
        shootRays(
            m,
            xs[i as usize],
            zs[i as usize],
            nangle,
            dangle,
            fangle,
            rc.as_mut_ptr(),
            kk,
            rsp,
            re,
            infofp,
        );
        if !fresnelfp.is_null() && !rayfp.is_null() {
            writeFresnel(
                m, rayfp, infofp, nxz, nangle, rsp, re, krecord, ffreq, prim, fresnelfp,
            );
        } else if !rayfp.is_null() {
            writeRays(
                rayfp,
                nxz,
                nangle,
                rsp,
                re,
                krecord,
                prim,
                infofp,
                caustic,
                nonsurface,
                &mut raycount,
            );
        }
        i += 1;
    }
    fprintf!(&outparfp, "%d\n", raycount);
    // The original retains only the last source's step lists and ray ends.
    if !wavefp.is_null() {
        writeWaves(wavefp, nt, nangle, rsp);
    }
    let bytes = if nangle == 0 {
        &[][..]
    } else {
        std::slice::from_raw_parts(
            re as *const u8,
            nangle as usize * std::mem::size_of::<RayEnd>(),
        )
    };
    if efwrite(bytes, &stdout()) / std::mem::size_of::<RayEnd>() != nangle as usize {
        err!("Error writing ray ends to stdout!\n");
    }
    exit(CWP_Exit());
}

unsafe fn shootRays(
    mut m: *mut Model,
    mut xs: f32,
    mut zs: f32,
    mut nangle: i32,
    mut dangle: f32,
    mut fangle: f32,
    mut rc: *mut RefCheck,
    mut kk: i32,
    mut rsp: *mut *mut RayStep,
    mut re: *mut RayEnd,
    mut infofp: *mut CFile,
) {
    let mut iangle: i32 = 0;
    let mut refl: i32 = 0;
    let mut stop: i32 = 0;
    let mut i: i32 = 0;
    let mut temp: i32 = 0;
    let mut test: i32 = 0;
    let mut angle: f32 = 0.0;
    let mut sb: f32 = 0.0;
    let mut zmax: f32 = 0.0;
    let mut xmax: f32 = 0.0;
    let mut eps: f32 = 0.0;
    let mut tris: *mut Tri = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    let mut ea: *mut EdgeAttributes = std::ptr::null_mut();
    let mut rs: *mut RayStep = std::ptr::null_mut();
    let mut rslast: *mut RayStep = std::ptr::null_mut();
    /* determine triangle in model containing source location */
    tris = insideTriInModel(m, std::ptr::null_mut(), zs, xs);
    /* if source sits on edge, move it by eps */
    if (({
        test = checkIfSourceOnEdge(tris, zs, xs);
        test
    }) != 0_i32)
    {
        zmax = (*(m)).xmax;
        eps = (*(m)).eps;
        if (test == 2_i32) {
            eps = ((5_i32 as f32) * eps);
        }
        zs = (if (zs < (zmax - eps)) {
            (zs + eps)
        } else {
            (zs - eps)
        });
        /* if we are still on an edge, move x value */
        if (checkIfSourceOnEdge(tris, zs, xs) != 0_i32) {
            xmax = (*(m)).ymax;
            xs = (if (xs < (xmax - eps)) {
                (xs + eps)
            } else {
                (xs - eps)
            });
        }
        /* find new triangle */
        tris = insideTriInModel(m, std::ptr::null_mut(), zs, xs);
    }
    /* determine sloth at source location (beginning of ray) */
    fa = ((*(tris)).fa as *mut FaceAttributes);
    sb = (((*(fa)).s00 + ((*(fa)).dsdx * xs)) + ((*(fa)).dsdz * zs));
    /* loop over takeoff angles */
    {
        iangle = 0_i32;
        angle = fangle;
        while (iangle < nangle) {
            /* initialize linked list of ray steps */
            rs = {
                (*(rsp).offset((iangle) as isize)) = ({
                    let p = calloc(
                        (1_i32 as usize) as usize,
                        (std::mem::size_of::<RayStep>()) as usize,
                    );
                    if p.is_null() {
                        syserr!("%s: malloc failed", "ealloc.c");
                    }
                    p
                } as *mut RayStep);
                (*(rsp).offset((iangle) as isize))
            };
            (*(rs)).sigma = (0.0_f64 as f32);
            (*(rs)).x = xs;
            (*(rs)).z = zs;
            (*(rs)).px = (((angle as f64).sin() * (sb as f64).sqrt()) as f32);
            (*(rs)).pz = (((angle as f64).cos() * (sb as f64).sqrt()) as f32);
            (*(rs)).q1 = (1.0_f64 as f32);
            (*(rs)).p1 = (0.0_f64 as f32);
            (*(rs)).q2 = (0.0_f64 as f32);
            (*(rs)).p2 = (1.0_f64 as f32);
            (*(rs)).t = (0.0_f64 as f32);
            (*(rs)).kmah = 0_i32;
            (*(rs)).nref = 0_i32;
            (*(rs)).eu = std::ptr::null_mut();
            (*(rs)).f = tris;
            (*(rs)).ampli = (1_i32 as f32);
            (*(rs)).ampliphase = (0_i32 as f32);
            (*(rs)).atten = (0_i32 as f32);
            (*(rs)).rsNext = std::ptr::null_mut();
            if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                fprintf!(&*(infofp), "\n****************************************\nRAY WITH TAKEOFF ANGLE: %g (degrees)\n****************************************\n" ,(((angle as f64) * 180.0_f64) / (3.141592653589793_f64)));
            }
            /* trace ray */
            loop {
                {
                    /* trace ray in triangle */
                    rs = traceRayInTri(rs);
                    /* remember last ray step */
                    rslast = rs;
                    /* edge attributes */
                    ea = ((*((*((*(rs)).eu)).e)).ea as *mut EdgeAttributes);
                    if (ea == ((std::ptr::null_mut()) as *mut EdgeAttributes)) {
                        /* null edges are transmitting */
                        refl = {
                            stop = {
                                temp = 0_i32;
                                temp
                            };
                            stop
                        };
                    } else {
                        /* check if edge is reflecting or stopping */
                        {
                            i = 0_i32;
                            temp = (*(ea)).k;
                            refl = {
                                stop = 0_i32;
                                stop
                            };
                            while (i < kk) {
                                if (temp == (*(rc).offset((i) as isize)).k) {
                                    if ((*((*(rc).offset((i) as isize)).hitseq)
                                        .offset((0_i32) as isize))
                                        == 1_i32)
                                    {
                                        refl = 1_i32;
                                    }
                                    if ((*((*(rc).offset((i) as isize)).hitseq)
                                        .offset((0_i32) as isize))
                                        == (-1_i32))
                                    {
                                        stop = 1_i32;
                                    }
                                    ((*(rc).offset((i) as isize)).nhits) =
                                        ((*(rc).offset((i) as isize)).nhits) + 1;
                                    ((*(rc).offset((i) as isize)).hitseq) =
                                        ((*(rc).offset((i) as isize)).hitseq).offset(1);
                                    break;
                                }
                                i = i + 1;
                            }
                        }
                    }
                    if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                        fprintf!(
                            &*(infofp),
                            "^^^\nInteraction with interface %d at (x=%g,z=%g).\n",
                            temp,
                            ((*(rslast)).x as f64),
                            ((*(rslast)).z as f64)
                        );
                    }
                    /* if ray at stopping edge, stop */
                    if (stop != 0) {
                        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                            fprintf!(&*(infofp), "Hit stopping edge.  Ray stopped.\n");
                        }
                        break;
                    } else {
                        if (refl != 0) {
                            rs = reflectRayFromEdge(rs, infofp);
                        } else {
                            rs = traceRayAcrossEdge(rs, infofp);
                        }
                    }
                }
                if !(rs != ((std::ptr::null_mut()) as *mut RayStep)) {
                    break;
                }
            }
            /* reinitialize RayChecks for new ray */
            {
                i = 0_i32;
                while (i < kk) {
                    (*(rc).offset((i) as isize)).hitseq = ((*(rc).offset((i) as isize)).hitseq)
                        .offset(-((*(rc).offset((i) as isize)).nhits) as isize);
                    (*(rc).offset((i) as isize)).nhits = 0_i32;
                    i = i + 1;
                }
            }
            /* save ray end parameters */
            fa = ((*((*(rslast)).f)).fa as *mut FaceAttributes);
            (*(re).offset((iangle) as isize)).sigma = (*(rslast)).sigma;
            (*(re).offset((iangle) as isize)).x = (*(rslast)).x;
            (*(re).offset((iangle) as isize)).z = (*(rslast)).z;
            (*(re).offset((iangle) as isize)).px = (*(rslast)).px;
            (*(re).offset((iangle) as isize)).pz = (*(rslast)).pz;
            (*(re).offset((iangle) as isize)).t = (*(rslast)).t;
            (*(re).offset((iangle) as isize)).q1 = (*(rslast)).q1;
            (*(re).offset((iangle) as isize)).p1 = (*(rslast)).p1;
            (*(re).offset((iangle) as isize)).q2 = (*(rslast)).q2;
            (*(re).offset((iangle) as isize)).p2 = (*(rslast)).p2;
            (*(re).offset((iangle) as isize)).kmah = (*(rslast)).kmah;
            (*(re).offset((iangle) as isize)).nref = (*(rslast)).nref;
            (*(re).offset((iangle) as isize)).sb = sb;
            (*(re).offset((iangle) as isize)).se =
                (((*(fa)).s00 + ((*(fa)).dsdx * (*(rslast)).x)) + ((*(fa)).dsdz * (*(rslast)).z));
            (*(re).offset((iangle) as isize)).dsdxe = (*(fa)).dsdx;
            (*(re).offset((iangle) as isize)).dsdze = (*(fa)).dsdz;
            (*(re).offset((iangle) as isize)).ab = angle;
            (*(re).offset((iangle) as isize)).kend = temp;
            (*(re).offset((iangle) as isize)).ampli = (*(rslast)).ampli;
            (*(re).offset((iangle) as isize)).ampliphase = (*(rslast)).ampliphase;
            (*(re).offset((iangle) as isize)).atten = (*(rslast)).atten;
            (*(re).offset((iangle) as isize)).dangle = dangle;
            /* optional output of rayend information */
            if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                fprintf!(&*(infofp), "---------------------------------------------------\n\nThe following information is stored at the rayend:\n\n");
                fprintf!(
                    &*(infofp),
                    "Ray stops at (%g,%g) at interface %d\n",
                    ((*(re).offset((iangle) as isize)).x as f64),
                    ((*(re).offset((iangle) as isize)).z as f64),
                    (*(re).offset((iangle) as isize)).kend
                );
                fprintf!(
                    &*(infofp),
                    "Takeoff angle=%g   Emergence angle=%g\n",
                    ((((*(re).offset((iangle) as isize)).ab as f64) * 180.0_f64)
                        / (3.141592653589793_f64)),
                    (((-((((*(re).offset((iangle) as isize)).px as f64) * 1.0_f64)
                        / ((*(re).offset((iangle) as isize)).se as f64).sqrt())
                    .asin())
                        * 180.0_f64)
                        / (3.141592653589793_f64))
                );
                fprintf!(
                    &*(infofp),
                    "Takeoff velocity=%g   Emergence velocity=%g\n",
                    (1.0_f64 / (sb as f64).sqrt()),
                    (1.0_f64 / ((*(re).offset((iangle) as isize)).se as f64).sqrt())
                );
                fprintf!(
                    &*(infofp),
                    "Slowness components px=%g   pz=%g\n",
                    ((*(re).offset((iangle) as isize)).px as f64),
                    ((*(re).offset((iangle) as isize)).pz as f64)
                );
                fprintf!(
                    &*(infofp),
                    "Traveltime=%g   Sigma=%g\n",
                    ((*(re).offset((iangle) as isize)).t as f64),
                    ((*(re).offset((iangle) as isize)).sigma as f64)
                );
                fprintf!(
                    &*(infofp),
                    "kmah index=%d   Number of reflections=%d\n",
                    (*(re).offset((iangle) as isize)).kmah,
                    (*(re).offset((iangle) as isize)).nref
                );
                fprintf!(
                    &*(infofp),
                    "Ray propagator q1=%g   q2=%g\n",
                    ((*(re).offset((iangle) as isize)).q1 as f64),
                    ((*(re).offset((iangle) as isize)).q2 as f64)
                );
                fprintf!(
                    &*(infofp),
                    "Ray propagator p1=%g   p2=%g\n",
                    ((*(re).offset((iangle) as isize)).p1 as f64),
                    ((*(re).offset((iangle) as isize)).p2 as f64)
                );
                fprintf!(
                    &*(infofp),
                    "Attenuation factor=%g   Angle increment=%g\n",
                    ((*(re).offset((iangle) as isize)).atten as f64),
                    (((dangle as f64) * 180.0_f64) / (3.141592653589793_f64))
                );
                fprintf!(
                    &*(infofp),
                    "Refl/Transm Amplitude=%g   Phase=%g\n",
                    ((*(re).offset((iangle) as isize)).ampli as f64),
                    ((*(re).offset((iangle) as isize)).ampliphase as f64)
                );
                fprintf!(
                    &*(infofp),
                    "---------------------------------------------------\n\n"
                );
            }
            iangle = iangle + 1;
            angle = ((angle as f32) + (dangle as f32)) as f32;
        }
    }
}
unsafe fn writeRays(
    mut rayfp: *mut CFile,
    mut nxz: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
    mut re: *mut RayEnd,
    mut krecord: i32,
    mut prim: i32,
    mut infofp: *mut CFile,
    caustic: i32,
    nonsurface: i32,
    raycount: &mut i32,
) {
    let mut iray: i32 = 0;
    let mut nrs: i32 = 0;
    let mut irs: i32 = 0;
    let mut nxzw: i32 = 0;
    let mut ns: i32 = 0;
    let mut is: i32 = 0;
    let mut icount: i32 = 0;
    let mut sigmamax: f32 = (0.0_f64 as f32);
    let mut sigma1: f32 = 0.0;
    let mut sigma2: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut dsdx: f32 = 0.0;
    let mut dsdz: f32 = 0.0;
    let mut ds: f32 = 0.0;
    let mut xs: f32 = 0.0;
    let mut zs: f32 = 0.0;
    let mut dsigma: f32 = 0.0;
    let mut rsi: *mut RayStep = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    /* initialize counter */
    icount = 0_i32;
    /* loop over rays */
    {
        iray = 0_i32;
        while (iray < nray) {
            if (caustic == 0_i32) {
                (*(re).offset((iray) as isize)).kmah = 1_i32;
            }
            if (nonsurface == 1_i32) {
                (*(re).offset((iray) as isize)).z = (0_i32 as f32);
            }
            if (((*(re).offset((iray) as isize)).kmah != 0_i32)
                && (((*(re).offset((iray) as isize)).z as f64) < 1e-06_f64))
            {
                if ((krecord == 2147483647_i32)
                    || ((krecord == (*(re).offset((iray) as isize)).kend)
                        && ((prim == 2147483647_i32)
                            || (prim == (*(re).offset((iray) as isize)).nref))))
                {
                    /* count rays */
                    icount = icount + 1;
                    (*raycount) = (*raycount) + 1;
                    /* count ray steps while determining maximum sigma */
                    {
                        nrs = 0_i32;
                        rsi = (*(rs).offset((iray) as isize));
                        while (rsi != ((std::ptr::null_mut()) as *mut RayStep)) {
                            sigmamax = (*(rsi)).sigma;
                            nrs = nrs + 1;
                            rsi = (*(rsi)).rsNext;
                        }
                    }
                    /* initialize number of (x,z) written for this ray */
                    nxzw = 0_i32;
                    /* loop over ray steps */
                    {
                        irs = 0_i32;
                        rsi = (*(rs).offset((iray) as isize));
                        while (irs < nrs) {
                            /* if sufficient number of (x,z) */
                            /* have been written */
                            if (nxzw >= nxz) {
                                break;
                            }
                            /* ray parameters at this step */
                            x = (*(rsi)).x;
                            z = (*(rsi)).z;
                            px = (*(rsi)).px;
                            pz = (*(rsi)).pz;
                            fa = ((*((*(rsi)).f)).fa as *mut FaceAttributes);
                            dsdx = (*(fa)).dsdx;
                            dsdz = (*(fa)).dsdz;
                            /* sigma at this step and next step */
                            sigma1 = (*(rsi)).sigma;
                            sigma2 = (if (irs < (nrs - 1_i32)) {
                                (*((*(rsi)).rsNext)).sigma
                            } else {
                                sigma1
                            });
                            /* number of sigma between this */
                            /* step and next */
                            ns = (((1_i32 as f64)
                                + ((((sigma2 - sigma1) * ((nxz - 1_i32) as f32)) as f64)
                                    / (2.0_f64 * (sigmamax as f64))))
                                as i32);
                            /* increment in sigma */
                            ds = ((sigma2 - sigma1) / (ns as f32));
                            /* loop over sigma */
                            {
                                is = 0_i32;
                                dsigma = (0.0_f64 as f32);
                                while (is < ns) {
                                    /* compute and write x and z */
                                    xs = (((x as f64)
                                        + ((dsigma as f64)
                                            * ((px as f64)
                                                + (((dsigma as f64) * 0.25_f64) * (dsdx as f64)))))
                                        as f32);
                                    zs = (((z as f64)
                                        + ((dsigma as f64)
                                            * ((pz as f64)
                                                + (((dsigma as f64) * 0.25_f64) * (dsdz as f64)))))
                                        as f32);
                                    if (efwrite(&[zs], &*(rayfp)) != (1_i32 as usize)) {
                                        err!("error writing rays!\n");
                                    }
                                    if (efwrite(&[xs], &*(rayfp)) != (1_i32 as usize)) {
                                        err!("error writing rays!\n");
                                    }
                                    /* if sufficient number written, */
                                    /* break */
                                    if ({
                                        let old = nxzw;
                                        nxzw = nxzw + 1;
                                        nxzw
                                    } >= nxz)
                                    {
                                        break;
                                    }
                                    is = is + 1;
                                    dsigma = ((dsigma as f32) + (ds as f32)) as f32;
                                }
                            }
                            irs = irs + 1;
                            rsi = (*(rsi)).rsNext;
                        }
                    }
                    /* if necessary, repeat last (x,z) */
                    while (nxzw < nxz) {
                        if (efwrite(&[zs], &*(rayfp)) != (1_i32 as usize)) {
                            err!("error writing rays!\n");
                        }
                        if (efwrite(&[xs], &*(rayfp)) != (1_i32 as usize)) {
                            err!("error writing rays!\n");
                        }
                        nxzw = nxzw + 1;
                    }
                } else {
                    /* manipulate rayend information */
                    (*(re).offset((iray) as isize)).kend = (-1_i32);
                }
            }
            iray = iray + 1;
        }
    }
    /* end loop over rays */
    if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
        fprintf!(&*(infofp), "\nWrote %d rays to ray file.\n", icount);
    }
}
unsafe fn writeWaves(
    mut wavefp: *mut CFile,
    mut nt: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
) {
    let errfp = Box::into_raw(Box::new(stderr()));
    let mut iray: i32 = 0;
    let mut it: i32 = 0;
    let mut tmax: f32 = 0.0;
    let mut dt: f32 = 0.0;
    let mut t1: f32 = 0.0;
    let mut t2: f32 = 0.0;
    let mut x1: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut px1: f32 = 0.0;
    let mut pz1: f32 = 0.0;
    let mut sigma1: f32 = 0.0;
    let mut sigma2: f32 = 0.0;
    let mut s00: f32 = 0.0;
    let mut dsdx: f32 = 0.0;
    let mut dsdz: f32 = 0.0;
    let mut ax: f32 = 0.0;
    let mut bx: f32 = 0.0;
    let mut az: f32 = 0.0;
    let mut bz: f32 = 0.0;
    let mut at: f32 = 0.0;
    let mut bt: f32 = 0.0;
    let mut ct: f32 = 0.0;
    let mut ti: f32 = 0.0;
    let mut sigmah: f32 = 0.0;
    let mut sigmam: f32 = 0.0;
    let mut tm: f32 = 0.0;
    let mut dsigma: f32 = 0.0;
    let mut sigma: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut xstart: f32 = 0.0;
    let mut zstart: f32 = 0.0;
    let mut rsi: *mut RayStep = std::ptr::null_mut();
    let mut rs1: *mut RayStep = std::ptr::null_mut();
    let mut rs2: *mut RayStep = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    /* determine maximum time for all rays */
    {
        iray = 0_i32;
        tmax = (0.0_f64 as f32);
        while (iray < nray) {
            {
                rsi = (*(rs).offset((iray) as isize));
                while (rsi != ((std::ptr::null_mut()) as *mut RayStep)) {
                    tmax = (max((tmax), ((*(rsi)).t)));
                    rsi = (*(rsi)).rsNext;
                }
            }
            iray = iray + 1;
        }
    }
    /* determine time sampling interval */
    dt = (tmax / ((if (nt > 1_i32) { (nt - 1_i32) } else { 1_i32 }) as f32));
    fprintf!(&*(errfp), "wavefront time increment = %g\n", (dt as f64));
    /* loop over rays */
    {
        iray = 0_i32;
        while (iray < nray) {
            /* initialize */
            rs1 = (*(rs).offset((iray) as isize));
            rs2 = (*(rs1)).rsNext;
            t1 = (*(rs1)).t;
            t2 = (*(rs2)).t;
            x1 = (*(rs1)).x;
            z1 = (*(rs1)).z;
            px1 = (*(rs1)).px;
            pz1 = (*(rs1)).pz;
            sigma1 = (*(rs1)).sigma;
            sigma2 = (*(rs2)).sigma;
            fa = ((*((*(rs1)).f)).fa as *mut FaceAttributes);
            s00 = (*(fa)).s00;
            dsdx = (*(fa)).dsdx;
            dsdz = (*(fa)).dsdz;
            ax = px1;
            bx = ((0.25_f64 * (dsdx as f64)) as f32);
            az = pz1;
            bz = ((0.25_f64 * (dsdz as f64)) as f32);
            at = ((s00 + (dsdx * x1)) + (dsdz * z1));
            bt = ((0.5_f64 * (((dsdx * px1) + (dsdz * pz1)) as f64)) as f32);
            ct = ((0.0833333_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64)) as f32);
            sigma = sigma1;
            t = t1;
            /* remember start x,z */
            xstart = x1;
            zstart = z1;
            /* loop over times */
            {
                it = 0_i32;
                ti = (0.0_f64 as f32);
                while (it < nt) {
                    /* if necessary, go to next ray step */
                    if (t2 < ti) {
                        loop {
                            {
                                rs1 = rs2;
                                rs2 = (*(rs1)).rsNext;
                                if (rs2 == ((std::ptr::null_mut()) as *mut RayStep)) {
                                    break;
                                }
                            }
                            if !((*(rs2)).t < ti) {
                                break;
                            }
                        }
                        if (rs2 == ((std::ptr::null_mut()) as *mut RayStep)) {
                            break;
                        }
                        t1 = (*(rs1)).t;
                        t2 = (*(rs2)).t;
                        x1 = (*(rs1)).x;
                        z1 = (*(rs1)).z;
                        px1 = (*(rs1)).px;
                        pz1 = (*(rs1)).pz;
                        sigma1 = (*(rs1)).sigma;
                        sigma2 = (*(rs2)).sigma;
                        fa = ((*((*(rs1)).f)).fa as *mut FaceAttributes);
                        s00 = (*(fa)).s00;
                        dsdx = (*(fa)).dsdx;
                        dsdz = (*(fa)).dsdz;
                        ax = px1;
                        bx = ((0.25_f64 * (dsdx as f64)) as f32);
                        az = pz1;
                        bz = ((0.25_f64 * (dsdz as f64)) as f32);
                        at = ((s00 + (dsdx * x1)) + (dsdz * z1));
                        bt = ((0.5_f64 * (((dsdx * px1) + (dsdz * pz1)) as f64)) as f32);
                        ct = ((0.0833333_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64)) as f32);
                        sigma = sigma1;
                        t = t1;
                    }
                    /* determine dsigma for this time by bisection */
                    sigmah = sigma2;
                    while (((t as f64) + (0.01_f64 * (dt as f64))) < (ti as f64)) {
                        sigmam = ((0.5_f64 * ((sigma + sigmah) as f64)) as f32);
                        dsigma = (sigmam - sigma1);
                        tm = (t1 + (dsigma * (at + (dsigma * (bt + (dsigma * ct))))));
                        if (tm < ti) {
                            t = tm;
                            sigma = sigmam;
                        } else {
                            sigmah = sigmam;
                        }
                    }
                    dsigma = (sigma - sigma1);
                    /* compute x and z */
                    x = (x1 + (dsigma * (ax + (dsigma * bx))));
                    z = (z1 + (dsigma * (az + (dsigma * bz))));
                    /* write x and z */
                    if (efwrite(&[z], &*(wavefp)) != (1_i32 as usize)) {
                        err!("error writing waves!\n");
                    }
                    if (efwrite(&[x], &*(wavefp)) != (1_i32 as usize)) {
                        err!("error writing waves!\n");
                    }
                    it = it + 1;
                    ti = ((ti as f32) + (dt as f32)) as f32;
                }
            }
            /* finish writing x and z */
            while (it < nt) {
                if (efwrite(&[zstart], &*(wavefp)) != (1_i32 as usize)) {
                    err!("error writing waves!\n");
                }
                if (efwrite(&[xstart], &*(wavefp)) != (1_i32 as usize)) {
                    err!("error writing waves!\n");
                }
                it = it + 1;
            }
            iray = iray + 1;
        }
    }
}
unsafe fn traceRayInTri(mut rs: *mut RayStep) -> *mut RayStep {
    let mut kmah: i32 = 0;
    let mut nref: i32 = 0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut ddt: f32 = 0.0;
    let mut s00: f32 = 0.0;
    let mut dsdx: f32 = 0.0;
    let mut dsdz: f32 = 0.0;
    let mut xa: f32 = 0.0;
    let mut za: f32 = 0.0;
    let mut xb: f32 = 0.0;
    let mut zb: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut a: f32 = 0.0;
    let mut b: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut ds: f32 = 0.0;
    let mut q: f32 = 0.0;
    let mut sigma: f32 = 0.0;
    let mut dsigma: f32 = 0.0;
    let mut dsigma1: f32 = 0.0;
    let mut dsigma2: f32 = 0.0;
    let mut small: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut ampli: f32 = 0.0;
    let mut ampliphase: f32 = 0.0;
    let mut atten: f32 = 0.0;
    let mut qfac: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    let mut eut: *mut EdgeUse = std::ptr::null_mut();
    let mut eusmall: *mut EdgeUse = std::ptr::null_mut();
    let mut euzero: *mut EdgeUse = std::ptr::null_mut();
    let mut f: *mut Face = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    /* get input parameters */
    sigma = (*(rs)).sigma;
    x = (*(rs)).x;
    z = (*(rs)).z;
    px = (*(rs)).px;
    pz = (*(rs)).pz;
    t = (*(rs)).t;
    q1 = (*(rs)).q1;
    p1 = (*(rs)).p1;
    q2 = (*(rs)).q2;
    p2 = (*(rs)).p2;
    kmah = (*(rs)).kmah;
    nref = (*(rs)).nref;
    eu = (*(rs)).eu;
    f = (*(rs)).f;
    ampli = (*(rs)).ampli;
    ampliphase = (*(rs)).ampliphase;
    atten = (*(rs)).atten;
    /* determine sloth function */
    fa = ((*(f)).fa as *mut FaceAttributes);
    s00 = (*(fa)).s00;
    dsdx = (*(fa)).dsdx;
    dsdz = (*(fa)).dsdz;
    /* determine Q-factor */
    qfac = (*(fa)).qfac;
    /* determine edge intersection with smallest positive dsigma */
    eusmall = std::ptr::null_mut();
    small = 3.40282347e+38_f32;
    eut = (*(f)).eu;
    loop {
        {
            /* edge endpoints */
            xa = (*((*((*(eut)).vu)).v)).y;
            za = (*((*((*(eut)).vu)).v)).x;
            xb = (*((*((*((*(eut)).euCW)).vu)).v)).y;
            zb = (*((*((*((*(eut)).euCW)).vu)).v)).x;
            /* coefficients b and c of equation to be solved for dsigma */
            dx = (xb - xa);
            dz = (zb - za);
            b = ((dx * pz) - (dz * px));
            c = ((if (eut == eu) {
                0.0_f64
            } else {
                (((dx * (z - za)) - (dz * (x - xa))) as f64)
            }) as f32);
            /* if sloth constant, solve linear equation */
            if (((dsdx as f64) == 0.0_f64) && ((dsdz as f64) == 0.0_f64)) {
                if ((b as f64) != 0.0_f64) {
                    dsigma1 = ((-c) / b);
                } else {
                    dsigma1 = 3.40282347e+38_f32;
                }
                dsigma2 = 3.40282347e+38_f32;
            } else {
                a = ((0.25_f64 * (((dx * dsdz) - (dz * dsdx)) as f64)) as f32);
                ds = ((((b * b) as f64) - ((4.0_f64 * (a as f64)) * (c as f64))) as f32);
                if ((ds as f64) < 0.0_f64) {
                    dsigma1 = {
                        dsigma2 = 3.40282347e+38_f32;
                        dsigma2
                    };
                } else {
                    q = (((-0.5_f64)
                        * (if ((b as f64) < 0.0_f64) {
                            ((b as f64) - (ds as f64).sqrt())
                        } else {
                            ((b as f64) + (ds as f64).sqrt())
                        })) as f32);
                    if ((a as f64) != 0.0_f64) {
                        dsigma1 = (q / a);
                    } else {
                        dsigma1 = 3.40282347e+38_f32;
                    }
                    if ((q as f64) != 0.0_f64) {
                        dsigma2 = (c / q);
                    } else {
                        dsigma2 = 3.40282347e+38_f32;
                    }
                }
            }
            /* remember edge with smallest positive dsigma */
            if ((0.0_f64 < (dsigma1 as f64)) && (dsigma1 < small)) {
                small = dsigma1;
                eusmall = eut;
            }
            if ((0.0_f64 < (dsigma2 as f64)) && (dsigma2 < small)) {
                small = dsigma2;
                eusmall = eut;
            }
            /* remember edge for which dsigma is zero */
            if ((dsigma1 as f64) == 0.0_f64) {
                euzero = eut;
            }
            if ((dsigma2 as f64) == 0.0_f64) {
                euzero = eut;
            }
            /* next edge use */
            eut = (*(eut)).euCW;
        }
        if !(eut != (*(f)).eu) {
            break;
        }
    }
    /* ray exits at edge with smallest positive dsigma */
    if (eusmall != ((std::ptr::null_mut()) as *mut EdgeUse)) {
        dsigma = small;
        eu = eusmall;
    } else {
        dsigma = (0.0_f64 as f32);
        eu = euzero;
    }
    /* update dynamic ray parameters */
    evaluateDynamic(
        dsigma,
        px,
        pz,
        dsdx,
        dsdz,
        (&raw mut q1),
        (&raw mut p1),
        (&raw mut q2),
        (&raw mut p2),
        (&raw mut kmah),
    );
    /* update ray parameters */
    sigma = ((sigma as f32) + (dsigma as f32)) as f32;
    ddt = (((dsigma as f64)
        * ((((s00 + (dsdx * x)) + (dsdz * z)) as f64)
            + ((dsigma as f64)
                * ((0.5_f64 * (((dsdx * px) + (dsdz * pz)) as f64))
                    + ((dsigma as f64)
                        * (0.0833333_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64)))))))
        as f32);
    t = ((t as f32) + (ddt as f32)) as f32;
    x = ((x as f64)
        + (((dsigma as f64) * ((px as f64) + ((0.25_f64 * (dsdx as f64)) * (dsigma as f64))))
            as f64)) as f32;
    z = ((z as f64)
        + (((dsigma as f64) * ((pz as f64) + ((0.25_f64 * (dsdz as f64)) * (dsigma as f64))))
            as f64)) as f32;
    px = ((px as f64) + (((0.5_f64 * (dsdx as f64)) * (dsigma as f64)) as f64)) as f32;
    pz = ((pz as f64) + (((0.5_f64 * (dsdz as f64)) * (dsigma as f64)) as f64)) as f32;
    /* don't let ray exit too close to a vertex */
    xa = (*((*((*(eu)).vu)).v)).y;
    dx = ((*((*((*((*(eu)).euCW)).vu)).v)).y - xa);
    za = (*((*((*(eu)).vu)).v)).x;
    dz = ((*((*((*((*(eu)).euCW)).vu)).v)).x - za);
    frac = (if ((abs((dx))) > (abs((dz)))) {
        ((x - xa) / dx)
    } else {
        ((z - za) / dz)
    });
    if ((frac as f64) < 0.0001_f64) {
        x = (((xa as f64) + (0.0001_f64 * (dx as f64))) as f32);
        z = (((za as f64) + (0.0001_f64 * (dz as f64))) as f32);
    } else {
        if ((frac as f64) > 0.9999_f64) {
            x = (((xa as f64) + (0.9999_f64 * (dx as f64))) as f32);
            z = (((za as f64) + (0.9999_f64 * (dz as f64))) as f32);
        }
    }
    /* compute contribution due to attenuation */
    if (qfac != 3.40282347e+38_f32) {
        atten = ((atten as f32) + ((ddt / qfac) as f32)) as f32;
    } else {
        atten = (0_i32 as f32);
    }
    /* return new raystep */
    (*(rs)).rsNext = ({
        let p = calloc(
            (1_i32 as usize) as usize,
            (std::mem::size_of::<RayStep>()) as usize,
        );
        if p.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        p
    } as *mut RayStep);
    rs = (*(rs)).rsNext;
    (*(rs)).sigma = sigma;
    (*(rs)).x = x;
    (*(rs)).z = z;
    (*(rs)).px = px;
    (*(rs)).pz = pz;
    (*(rs)).t = t;
    (*(rs)).q1 = q1;
    (*(rs)).p1 = p1;
    (*(rs)).q2 = q2;
    (*(rs)).p2 = p2;
    (*(rs)).kmah = kmah;
    (*(rs)).nref = nref;
    (*(rs)).eu = eu;
    (*(rs)).f = f;
    (*(rs)).rsNext = std::ptr::null_mut();
    (*(rs)).atten = atten;
    (*(rs)).ampli = ampli;
    (*(rs)).ampliphase = ampliphase;
    return rs;
}
unsafe fn traceRayAcrossEdge(mut rs: *mut RayStep, mut infofp: *mut CFile) -> *mut RayStep {
    let mut kmah: i32 = 0;
    let mut nref: i32 = 0;
    let mut sigma: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut temp1: f32 = 0.0;
    let mut temp2: f32 = 0.0;
    let mut s00: f32 = 0.0;
    let mut s1: f32 = 0.0;
    let mut ds1dx: f32 = 0.0;
    let mut ds1dz: f32 = 0.0;
    let mut s2: f32 = 0.0;
    let mut ds2dx: f32 = 0.0;
    let mut ds2dz: f32 = 0.0;
    let mut px1: f32 = 0.0;
    let mut pz1: f32 = 0.0;
    let mut px1r: f32 = 0.0;
    let mut pz1r: f32 = 0.0;
    let mut px2r: f32 = 0.0;
    let mut pz2r: f32 = 0.0;
    let mut pz2rs: f32 = 0.0;
    let mut c1ov1: f32 = 0.0;
    let mut c2ov2: f32 = 0.0;
    let mut oc2s: f32 = 0.0;
    let mut g: f32 = 0.0;
    let mut cterm: f32 = 0.0;
    let mut iterm: f32 = 0.0;
    let mut fac1: f32 = 0.0;
    let mut fac2: f32 = 0.0;
    let mut dv1dl: f32 = 0.0;
    let mut dv2dl: f32 = 0.0;
    let mut dv1dm: f32 = 0.0;
    let mut dv2dm: f32 = 0.0;
    let mut scale: f32 = 0.0;
    let mut gx: f32 = 0.0;
    let mut gz: f32 = 0.0;
    let mut hx: f32 = 0.0;
    let mut h_z: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut taper: f32 = 0.0;
    let mut ampli: f32 = 0.0;
    let mut atten: f32 = 0.0;
    let mut dens1: f32 = 0.0;
    let mut dens2: f32 = 0.0;
    let mut ampliphase: f32 = 0.0;
    let mut coeff: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    let mut eum: *mut EdgeUse = std::ptr::null_mut();
    let mut eua: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut euma: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut f: *mut Face = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    /* get input parameters */
    sigma = (*(rs)).sigma;
    x = (*(rs)).x;
    z = (*(rs)).z;
    px = (*(rs)).px;
    pz = (*(rs)).pz;
    t = (*(rs)).t;
    q1 = (*(rs)).q1;
    p1 = (*(rs)).p1;
    q2 = (*(rs)).q2;
    p2 = (*(rs)).p2;
    kmah = (*(rs)).kmah;
    nref = (*(rs)).nref;
    eu = (*(rs)).eu;
    f = (*(rs)).f;
    atten = (*(rs)).atten;
    ampli = (*(rs)).ampli;
    ampliphase = (*(rs)).ampliphase;
    /* check for boundary */
    if ((*((*(eu)).euMate)).f == ((std::ptr::null_mut()) as *mut Face)) {
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(&*(infofp), "Transmitted outside model.  Ray stopped.\n");
        }
        return std::ptr::null_mut();
    }
    /* determine sloth on this side of edge */
    fa = ((*(f)).fa as *mut FaceAttributes);
    s00 = (*(fa)).s00;
    ds1dx = (*(fa)).dsdx;
    ds1dz = (*(fa)).dsdz;
    s1 = ((s00 + (ds1dx * x)) + (ds1dz * z));
    /* determine density on this side */
    dens1 = (*(fa)).dens;
    /* determine sloth on other side of edge */
    eum = (*(eu)).euMate;
    f = (*(eum)).f;
    fa = ((*(f)).fa as *mut FaceAttributes);
    s00 = (*(fa)).s00;
    ds2dx = (*(fa)).dsdx;
    ds2dz = (*(fa)).dsdz;
    s2 = ((s00 + (ds2dx * x)) + (ds2dz * z));
    /* determine density on other side of the edge */
    dens2 = (*(fa)).dens;
    /* edge vector */
    dx = ((*((*((*(eum)).vu)).v)).y - (*((*((*(eu)).vu)).v)).y);
    dz = ((*((*((*(eum)).vu)).v)).x - (*((*((*(eu)).vu)).v)).x);
    /* fractional distance along edge */
    frac = (if ((abs((dx))) > (abs((dz)))) {
        ((x - (*((*((*(eu)).vu)).v)).y) / dx)
    } else {
        ((z - (*((*((*(eu)).vu)).v)).x) / dz)
    });
    /* here is speeding up potential */
    /* linearly interpolate unit vector g tangent to edge */
    eua = ((*(eu)).eua as *mut EdgeUseAttributes);
    euma = ((*(eum)).eua as *mut EdgeUseAttributes);
    if ((eua != ((std::ptr::null_mut()) as *mut EdgeUseAttributes))
        && (euma != ((std::ptr::null_mut()) as *mut EdgeUseAttributes)))
    {
        gx = ((((frac * (*(euma)).tx) as f64) - ((1.0_f64 - (frac as f64)) * ((*(eua)).tx as f64)))
            as f32);
        gz = ((((frac * (*(euma)).tz) as f64) - ((1.0_f64 - (frac as f64)) * ((*(eua)).tz as f64)))
            as f32);
    } else {
        gx = (-dx);
        gz = (-dz);
    }
    scale = ((1.0_f64 / (((gx * gx) + (gz * gz)) as f64).sqrt()) as f32);
    gx = ((gx as f32) * (scale as f32)) as f32;
    gz = ((gz as f32) * (scale as f32)) as f32;
    /* unit vector h normal to edge */
    hx = (-gz);
    h_z = gx;
    /* remember ray parameters on this side */
    px1 = px;
    pz1 = pz;
    /* rotated ray parameters on this side */
    px1r = ((px * h_z) - (pz * hx));
    pz1r = ((px * hx) + (pz * h_z));
    /* rotated ray parameters on other side */
    px2r = px1r;
    pz2rs = (s2 - (px2r * px2r));
    /* post-critical incidence*/
    if ((pz2rs as f64) <= 0.0_f64) {
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(&*(infofp), "Post-critical transmission.  Ray stopped.\n");
        }
        return std::ptr::null_mut();
    }
    pz2r = ((pz2rs as f64).sqrt() as f32);
    /* ray parameters on other side */
    px = ((px2r * h_z) + (pz2r * hx));
    pz = ((pz2r * h_z) - (px2r * hx));
    /* curvature term */
    c1ov1 = pz1r;
    c2ov2 = pz2r;
    oc2s = (s2 / pz2rs);
    if ((eua != ((std::ptr::null_mut()) as *mut EdgeUseAttributes))
        && (euma != ((std::ptr::null_mut()) as *mut EdgeUseAttributes)))
    {
        g = ((((frac * (*(euma)).c) as f64) - ((1.0_f64 - (frac as f64)) * ((*(eua)).c as f64)))
            as f32);
        cterm = ((g * oc2s) * (c1ov1 - c2ov2));
    } else {
        cterm = (0.0_f64 as f32);
    }
    /* update dynamic ray parameters */
    scale = ((((pz2r as f64) * (s1 as f64).sqrt()) / ((pz1r as f64) * (s2 as f64).sqrt())) as f32);
    q1 = (q1 * scale);
    p1 = ((p1 / scale) + (cterm * q1));
    q2 = (q2 * scale);
    p2 = ((p2 / scale) + (cterm * q2));
    /* velocity derivatives tangent and normal to ray */
    fac1 = (((-0.5_f64) / ((s1 * s1) as f64)) as f32);
    fac2 = (((-0.5_f64) / ((s2 * s2) as f64)) as f32);
    dv1dl = (fac1 * ((px1 * ds1dx) + (pz1 * ds1dz)));
    dv2dl = (fac2 * ((px * ds2dx) + (pz * ds2dz)));
    dv1dm = (fac1 * ((pz1 * ds1dx) - (px1 * ds1dz)));
    dv2dm = (fac2 * ((pz * ds2dx) - (px * ds2dz)));
    /* inhomogeneity term */
    iterm = (((((-px1r) * oc2s) as f64)
        * ((2.0_f64 * (((dv1dm * c1ov1) - (dv2dm * c2ov2)) as f64))
            + ((px1r * (dv1dl - dv2dl)) as f64))) as f32);
    /* update dynamic ray parameters */
    p1 = ((p1 as f32) + ((iterm * q1) as f32)) as f32;
    p2 = ((p2 as f32) + ((iterm * q2) as f32)) as f32;
    /* transmission effects on amplitudes */
    if ((((abs((s1 - s2))) as f64) > 0.001_f64) || (dens1 != dens2)) {
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(&*(infofp), "Transmitted WITH influence on amplitude.\n");
        }
        if (dens1 == 3.40282347e+38_f32) {
            dens1 = {
                dens2 = (1.0_f64 as f32);
                dens2
            };
        }
        temp1 = (dens2 / dens1);
        temp2 = (pz2r / pz1r);
        coeff = ((1.0_f64 + (((temp1 - temp2) / (temp1 + temp2)) as f64)) as f32);
        /* check if too close to critical */
        taper = ((px1r * px1r) / s2);
        taper = ((if ((taper as f64) > 0.9_f64) {
            ((((taper as f64) - 0.9_f64) * 5.0_f64) * (3.141592653589793_f64)).cos()
        } else {
            1.0_f64
        }) as f32);
        /* complete amplitude coeff */
        ampli = ((ampli as f64)
            * ((((coeff as f64) * ((pz2r / pz1r) as f64).sqrt()) * (taper as f64)) as f64))
            as f32;
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(
                &*(infofp),
                "  --incident side: v1=%g  dens1=%g  incidence angle=%g\n",
                ((1_i32 as f64) / (s1 as f64).sqrt()),
                (dens1 as f64),
                ((((pz1r as f64) / (s1 as f64).sqrt()).acos() * (180_i32 as f64))
                    / (3.141592653589793_f64))
            );
            fprintf!(
                &*(infofp),
                "  --opposite side: v2=%g  dens2=%g  refracted angle=%g\n",
                ((1_i32 as f64) / (s2 as f64).sqrt()),
                (dens2 as f64),
                ((((pz2r as f64) / (s2 as f64).sqrt()).acos() * (180_i32 as f64))
                    / (3.141592653589793_f64))
            );
            fprintf!(&*(infofp), "  --transmission coeff: %g\n", (coeff as f64));
            fprintf!(
                &*(infofp),
                "  --critical transm. taper: %g\n",
                (taper as f64)
            );
            fprintf!(&*(infofp), "  --total ampli. coeff: %g\n", (ampli as f64));
        }
        /* grazing incidence*/
        if (((pz1r * pz1r) as f64) <= (0.008_f64 * (s1 as f64))) {
            if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                fprintf!(&*(infofp), "Grazing incidence.  Ray stopped.\n");
            }
            return std::ptr::null_mut();
        }
    } else {
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(&*(infofp), "Transmitted WITHOUT influencing amplitude.\n");
        }
    }
    /* return new raystep */
    (*(rs)).rsNext = ({
        let p = calloc(
            (1_i32 as usize) as usize,
            (std::mem::size_of::<RayStep>()) as usize,
        );
        if p.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        p
    } as *mut RayStep);
    rs = (*(rs)).rsNext;
    (*(rs)).sigma = sigma;
    (*(rs)).x = x;
    (*(rs)).z = z;
    (*(rs)).px = px;
    (*(rs)).pz = pz;
    (*(rs)).t = t;
    (*(rs)).q1 = q1;
    (*(rs)).p1 = p1;
    (*(rs)).q2 = q2;
    (*(rs)).p2 = p2;
    (*(rs)).kmah = kmah;
    (*(rs)).nref = nref;
    (*(rs)).eu = eum;
    (*(rs)).f = f;
    (*(rs)).rsNext = std::ptr::null_mut();
    (*(rs)).ampli = ampli;
    (*(rs)).ampliphase = ampliphase;
    (*(rs)).atten = atten;
    return rs;
}
unsafe fn reflectRayFromEdge(mut rs: *mut RayStep, mut infofp: *mut CFile) -> *mut RayStep {
    let mut kmah: i32 = 0;
    let mut nref: i32 = 0;
    let mut sigma: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut s00: f32 = 0.0;
    let mut dsdx: f32 = 0.0;
    let mut dsdz: f32 = 0.0;
    let mut s: f32 = 0.0;
    let mut coeff: f32 = 0.0;
    let mut px1: f32 = 0.0;
    let mut pz1: f32 = 0.0;
    let mut pxr: f32 = 0.0;
    let mut pzr: f32 = 0.0;
    let mut c1ov1: f32 = 0.0;
    let mut c2ov2: f32 = 0.0;
    let mut oc2s: f32 = 0.0;
    let mut g: f32 = 0.0;
    let mut cterm: f32 = 0.0;
    let mut iterm: f32 = 0.0;
    let mut dv1dl: f32 = 0.0;
    let mut dv2dl: f32 = 0.0;
    let mut dv1dm: f32 = 0.0;
    let mut dv2dm: f32 = 0.0;
    let mut scale: f32 = 0.0;
    let mut gx: f32 = 0.0;
    let mut gz: f32 = 0.0;
    let mut hx: f32 = 0.0;
    let mut h_z: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut atten: f32 = 0.0;
    let mut ampli: f32 = 0.0;
    let mut ampliphase: f32 = 0.0;
    let mut dens1: f32 = 0.0;
    let mut dens2: f32 = 0.0;
    let mut s2: f32 = 0.0;
    let mut ds2dx: f32 = 0.0;
    let mut ds2dz: f32 = 0.0;
    let mut temp1: f32 = 0.0;
    let mut temp2: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    let mut eum: *mut EdgeUse = std::ptr::null_mut();
    let mut eua: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut euma: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut f: *mut Face = std::ptr::null_mut();
    let mut fn_: *mut Face = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    /* get input parameters */
    sigma = (*(rs)).sigma;
    x = (*(rs)).x;
    z = (*(rs)).z;
    px = (*(rs)).px;
    pz = (*(rs)).pz;
    t = (*(rs)).t;
    q1 = (*(rs)).q1;
    p1 = (*(rs)).p1;
    q2 = (*(rs)).q2;
    p2 = (*(rs)).p2;
    kmah = (*(rs)).kmah;
    nref = (*(rs)).nref;
    eu = (*(rs)).eu;
    f = (*(rs)).f;
    ampli = (*(rs)).ampli;
    ampliphase = (*(rs)).ampliphase;
    atten = (*(rs)).atten;
    /* determine sloth on incident side of edge */
    fa = ((*(f)).fa as *mut FaceAttributes);
    s00 = (*(fa)).s00;
    dsdx = (*(fa)).dsdx;
    dsdz = (*(fa)).dsdz;
    s = ((s00 + (dsdx * x)) + (dsdz * z));
    /* determine dens on incident side of edge */
    dens1 = (*(fa)).dens;
    /* edge vector */
    eum = (*(eu)).euMate;
    /* determine sloth on other side of edge */
    fn_ = (*(eum)).f;
    fa = ((*(fn_)).fa as *mut FaceAttributes);
    s00 = (*(fa)).s00;
    ds2dx = (*(fa)).dsdx;
    ds2dz = (*(fa)).dsdz;
    s2 = ((s00 + (ds2dx * x)) + (ds2dz * z));
    /* determine density on other side of the edge */
    dens2 = (*(fa)).dens;
    if (dens1 == 3.40282347e+38_f32) {
        dens1 = (1.0_f64 as f32);
    }
    if (dens2 == 3.40282347e+38_f32) {
        dens2 = (1.0_f64 as f32);
    }
    dx = ((*((*((*(eum)).vu)).v)).y - (*((*((*(eu)).vu)).v)).y);
    dz = ((*((*((*(eum)).vu)).v)).x - (*((*((*(eu)).vu)).v)).x);
    /* fractional distance along edge */
    frac = (if ((abs((dx))) > (abs((dz)))) {
        ((x - (*((*((*(eu)).vu)).v)).y) / dx)
    } else {
        ((z - (*((*((*(eu)).vu)).v)).x) / dz)
    });
    /* linearly interpolate unit vector g tangent to edge */
    eua = ((*(eu)).eua as *mut EdgeUseAttributes);
    euma = ((*(eum)).eua as *mut EdgeUseAttributes);
    gx = ((((frac * (*(euma)).tx) as f64) - ((1.0_f64 - (frac as f64)) * ((*(eua)).tx as f64)))
        as f32);
    gz = ((((frac * (*(euma)).tz) as f64) - ((1.0_f64 - (frac as f64)) * ((*(eua)).tz as f64)))
        as f32);
    scale = ((1.0_f64 / (((gx * gx) + (gz * gz)) as f64).sqrt()) as f32);
    gx = ((gx as f32) * (scale as f32)) as f32;
    gz = ((gz as f32) * (scale as f32)) as f32;
    /* unit vector h normal to edge */
    hx = (-gz);
    h_z = gx;
    /* remember incident ray parameters */
    px1 = px;
    pz1 = pz;
    /* rotated incident ray parameters */
    pxr = ((px * h_z) - (pz * hx));
    pzr = ((px * hx) + (pz * h_z));
    /* rotated reflected ray parameters */
    pxr = pxr;
    pzr = (-pzr);
    /* reflected ray parameters */
    px = ((pxr * h_z) + (pzr * hx));
    pz = ((pzr * h_z) - (pxr * hx));
    /* curvature term */
    c1ov1 = {
        c2ov2 = (-pzr);
        c2ov2
    };
    oc2s = (s / (pzr * pzr));
    g = ((((frac * (*(euma)).c) as f64) - ((1.0_f64 - (frac as f64)) * ((*(eua)).c as f64)))
        as f32);
    cterm = ((g * oc2s) * (c1ov1 + c2ov2));
    /* update dynamic ray parameters */
    scale = ((-c2ov2) / c1ov1);
    q1 = (q1 * scale);
    p1 = ((p1 / scale) + (cterm * q1));
    q2 = (q2 * scale);
    p2 = ((p2 / scale) + (cterm * q2));
    /* if sloth not constant */
    if (((dsdx as f64) != 0.0_f64) || ((dsdz as f64) != 0.0_f64)) {
        /* velocity derivatives tangent and normal to ray */
        scale = (((-0.5_f64) / ((s * s) as f64)) as f32);
        dv1dl = (scale * ((px1 * dsdx) + (pz1 * dsdz)));
        dv2dl = (scale * ((px * dsdx) + (pz * dsdz)));
        dv1dm = (scale * ((pz1 * dsdx) - (px1 * dsdz)));
        dv2dm = (scale * ((pz * dsdx) - (px * dsdz)));
        /* inhomogeneity term */
        iterm = (((((-pxr) * oc2s) as f64)
            * ((2.0_f64 * (((dv1dm * c1ov1) + (dv2dm * c2ov2)) as f64))
                + ((pxr * (dv1dl - dv2dl)) as f64))) as f32);
        /* update dynamic ray parameters */
        p1 = ((p1 as f32) + ((iterm * q1) as f32)) as f32;
        p2 = ((p2 as f32) + ((iterm * q2) as f32)) as f32;
    }
    /* increment number of reflections */
    nref = nref + 1;
    if ((s2 - (pxr * pxr)) >= (0_i32 as f32)) {
        /* pre-critical reflection */
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(
                &*(infofp),
                "Reflected pre-critically at %g degrees.\n",
                (((((-pzr) as f64) / (s as f64).sqrt()).acos() * (180_i32 as f64))
                    / (3.141592653589793_f64))
            );
        }
        temp1 = ((((((-pzr) * dens2) / dens1) as f64) / (s as f64).sqrt()) as f32);
        temp2 = ((((s2 / s) - ((pxr * pxr) / s)) as f64).sqrt() as f32);
        coeff = ((temp1 - temp2) / (temp1 + temp2));
        ampli = ((ampli as f32) * (coeff as f32)) as f32;
    } else {
        /* post-critical reflection */
        coeff = (1_i32 as f32);
        ampliphase = ((ampliphase as f64)
            - (((2_i32 as f64)
                * (((-(((pxr * pxr) - s2) as f64).sqrt()) * (dens1 as f64))
                    / ((dens2 * pzr) as f64))
                    .atan()) as f64)) as f32;
        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
            fprintf!(
                &*(infofp),
                "Reflected post-critically at %g degrees.\n",
                (((((-pzr) as f64) / (s as f64).sqrt()).acos() * (180_i32 as f64))
                    / (3.141592653589793_f64))
            );
        }
    }
    if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
        fprintf!(
            &*(infofp),
            "  --incident side: v1=%g  dens1=%g\n",
            ((1_i32 as f64) / (s as f64).sqrt()),
            (dens1 as f64)
        );
        fprintf!(
            &*(infofp),
            "  --opposite side: v2=%g  dens2=%g\n",
            ((1_i32 as f64) / (s2 as f64).sqrt()),
            (dens2 as f64)
        );
        fprintf!(&*(infofp), "  --reflection coeff: %g\n", (coeff as f64));
        fprintf!(&*(infofp), "  --total ampl coeff: %g\n", (ampli as f64));
        fprintf!(
            &*(infofp),
            "  --total phase shift: %g\n",
            (ampliphase as f64)
        );
    }
    /* return new raystep */
    (*(rs)).rsNext = ({
        let p = calloc(
            (1_i32 as usize) as usize,
            (std::mem::size_of::<RayStep>()) as usize,
        );
        if p.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        p
    } as *mut RayStep);
    rs = (*(rs)).rsNext;
    (*(rs)).sigma = sigma;
    (*(rs)).x = x;
    (*(rs)).z = z;
    (*(rs)).px = px;
    (*(rs)).pz = pz;
    (*(rs)).t = t;
    (*(rs)).q1 = q1;
    (*(rs)).p1 = p1;
    (*(rs)).q2 = q2;
    (*(rs)).p2 = p2;
    (*(rs)).kmah = kmah;
    (*(rs)).nref = nref;
    (*(rs)).eu = eu;
    (*(rs)).f = f;
    (*(rs)).rsNext = std::ptr::null_mut();
    (*(rs)).ampli = ampli;
    (*(rs)).ampliphase = ampliphase;
    (*(rs)).atten = atten;
    return rs;
}
unsafe fn evaluateDynamic(
    mut sigma: f32,
    mut px: f32,
    mut pz: f32,
    mut dsdx: f32,
    mut dsdz: f32,
    mut q1: *mut f32,
    mut p1: *mut f32,
    mut q2: *mut f32,
    mut p2: *mut f32,
    mut kmah: *mut i32,
) {
    let mut s0: f64 = 0.0;
    let mut s1: f64 = 0.0;
    let mut s2: f64 = 0.0;
    let mut ss: f64 = 0.0;
    let mut scale: f64 = 0.0;
    let mut a: f64 = 0.0;
    let mut b: f64 = 0.0;
    let mut c: f64 = 0.0;
    let mut d: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut q1i: f64 = 0.0;
    let mut p1i: f64 = 0.0;
    let mut q2i: f64 = 0.0;
    let mut p2i: f64 = 0.0;
    let mut q1o: f64 = 0.0;
    let mut p1o: f64 = 0.0;
    let mut q2o: f64 = 0.0;
    let mut p2o: f64 = 0.0;
    /* get input dynamic ray parameters */
    q1i = ((*(q1)) as f64);
    p1i = ((*(p1)) as f64);
    q2i = ((*(q2)) as f64);
    p2i = ((*(p2)) as f64);
    /* trivial case:  constant sloth or ray parallel to sloth gradient */
    if ((pz * dsdx) == (px * dsdz)) {
        q1o = (q1i + (p1i * (sigma as f64)));
        q2o = (q2i + (p2i * (sigma as f64)));
        p2o = p2i;
        p1o = p1i;
    } else {
        /* constants */
        s0 = (((px * px) + (pz * pz)) as f64);
        s1 = (((px * dsdx) + (pz * dsdz)) as f64);
        s2 = (0.25_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64));
        ss = ((((s2 * (sigma as f64)) + s1) * (sigma as f64)) + s0);
        scale = (1.0_f64 / (s0 * ss).sqrt());
        a = (s0 + ((0.5_f64 * s1) * (sigma as f64)));
        b = (((((0.25_f64 * s1) * s1) / s0) - s2) * (sigma as f64));
        c = (s0 + (s1 * (sigma as f64)));
        d = ((0.5_f64 * s1) + (2.0_f64 * b));
        b = ((b as f64) * ((sigma as f64) as f64)) as f64;
        e = (((0.5_f64 * s1) + (s2 * (sigma as f64))) / ss);
        /* update q1 and q2 */
        q1o = (scale * ((a * (q1i + (p1i * (sigma as f64)))) + (b * q1i)));
        q2o = (scale * ((a * (q2i + (p2i * (sigma as f64)))) + (b * q2i)));
        /* update p1 and p2 */
        p1o = ((scale * ((c * p1i) + (d * q1i))) - (e * q1o));
        p2o = ((scale * ((c * p2i) + (d * q2i))) - (e * q2o));
    }
    /* update kmah index */
    if ((((q2i * q2o) >= 0.0_f64) && ((p2i * p2o) < 0.0_f64)) && ((q2i * p2i) < 0.0_f64)) {
        (*(kmah)) = (((*(kmah)) as i32) + (2_i32 as i32)) as i32;
    } else {
        if (((q2o == 0.0_f64) || ((q2i * q2o) < 0.0_f64))
            || (((q2i == 0.0_f64) && ((p2i * p2o) < 0.0_f64)) && ((q2o * p2o) > 0.0_f64)))
        {
            (*(kmah)) = (((*(kmah)) as i32) + (1_i32 as i32)) as i32;
        }
    }
    /* return updated dynamic ray parameters */
    (*(q1)) = (q1o as f32);
    (*(p1)) = (p1o as f32);
    (*(q2)) = (q2o as f32);
    (*(p2)) = (p2o as f32);
}
unsafe fn evaluateDynamic2(
    mut sigma: f32,
    mut px: f32,
    mut pz: f32,
    mut dsdx: f32,
    mut dsdz: f32,
    mut q1: *mut f32,
    mut p1: *mut f32,
    mut q2: *mut f32,
    mut p2: *mut f32,
    mut s0: f32,
    mut s1: f32,
    mut s2: f32,
) {
    let mut ss: f64 = 0.0;
    let mut scale: f64 = 0.0;
    let mut a: f64 = 0.0;
    let mut b: f64 = 0.0;
    let mut c: f64 = 0.0;
    let mut d: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut q1i: f64 = 0.0;
    let mut p1i: f64 = 0.0;
    let mut q2i: f64 = 0.0;
    let mut p2i: f64 = 0.0;
    let mut q1o: f64 = 0.0;
    let mut p1o: f64 = 0.0;
    let mut q2o: f64 = 0.0;
    let mut p2o: f64 = 0.0;
    /* get input dynamic ray parameters */
    q1i = ((*(q1)) as f64);
    p1i = ((*(p1)) as f64);
    q2i = ((*(q2)) as f64);
    p2i = ((*(p2)) as f64);
    /* trivial case:  constant sloth or ray parallel to sloth gradient */
    if ((pz * dsdx) == (px * dsdz)) {
        q1o = (q1i + (p1i * (sigma as f64)));
        q2o = (q2i + (p2i * (sigma as f64)));
        p2o = p2i;
        p1o = p1i;
    } else {
        /* constants */
        ss = (((((s2 * sigma) + s1) * sigma) + s0) as f64);
        scale = (1.0_f64 / ((s0 as f64) * ss).sqrt());
        a = ((s0 as f64) + ((0.5_f64 * (s1 as f64)) * (sigma as f64)));
        b = (((((0.25_f64 * (s1 as f64)) * (s1 as f64)) / (s0 as f64)) - (s2 as f64))
            * (sigma as f64));
        c = ((s0 + (s1 * sigma)) as f64);
        d = ((0.5_f64 * (s1 as f64)) + (2.0_f64 * b));
        b = ((b as f64) * ((sigma as f64) as f64)) as f64;
        e = (((0.5_f64 * (s1 as f64)) + ((s2 * sigma) as f64)) / ss);
        /* update q1 and q2 */
        q1o = (scale * ((a * (q1i + (p1i * (sigma as f64)))) + (b * q1i)));
        q2o = (scale * ((a * (q2i + (p2i * (sigma as f64)))) + (b * q2i)));
        /* update p1 and p2 */
        p1o = ((scale * ((c * p1i) + (d * q1i))) - (e * q1o));
        p2o = ((scale * ((c * p2i) + (d * q2i))) - (e * q2o));
    }
    /* return updated dynamic ray parameters */
    (*(q1)) = (q1o as f32);
    (*(p1)) = (p1o as f32);
    (*(q2)) = (q2o as f32);
    (*(p2)) = (p2o as f32);
}
unsafe fn checkIfSourceOnEdge(mut tris: *mut Face, mut zs: f32, mut xs: f32) -> i32 {
    let mut x1: f32 = 0.0;
    let mut x2: f32 = 0.0;
    let mut x3: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut z2: f32 = 0.0;
    let mut z3: f32 = 0.0;
    let mut m12: f32 = 0.0;
    let mut m13: f32 = 0.0;
    let mut m23: f32 = 0.0;
    let mut b12: f32 = 0.0;
    let mut b13: f32 = 0.0;
    let mut b23: f32 = 0.0;
    let mut eps: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    eu = (*(tris)).eu;
    eps = (*((*(tris)).m)).eps;
    /* get vertices */
    x1 = (*((*((*(eu)).vu)).v)).y;
    z1 = (*((*((*(eu)).vu)).v)).x;
    x2 = (*((*((*((*(eu)).euCW)).vu)).v)).y;
    z2 = (*((*((*((*(eu)).euCW)).vu)).v)).x;
    x3 = (*((*((*((*(eu)).euCCW)).vu)).v)).y;
    z3 = (*((*((*((*(eu)).euCCW)).vu)).v)).x;
    /* source is sitting on vertex */
    if ((((((xs - x1) * (xs - x1)) + ((zs - z1) * (zs - z1))) < eps)
        || ((((xs - x2) * (xs - x2)) + ((zs - z2) * (zs - z2))) < eps))
        || ((((xs - x3) * (xs - x3)) + ((zs - z3) * (zs - z3))) < eps))
    {
        return 2_i32;
    }
    /* check special cases and compute slope of edges */
    if (((abs((x1 - x2))) as f64) < (0.1_f64 * (eps as f64))) {
        if (((abs((xs - x1))) as f64) < (0.1_f64 * (eps as f64))) {
            return 1_i32;
        } else {
            m12 = ((0.01_f64 * (3.40282347e+38_f32 as f64)) as f32);
        }
    } else {
        m12 = ((z1 - z2) / (x1 - x2));
    }
    if (((abs((x1 - x3))) as f64) < (0.1_f64 * (eps as f64))) {
        if (((abs((xs - x1))) as f64) < (0.1_f64 * (eps as f64))) {
            return 1_i32;
        } else {
            m13 = ((0.01_f64 * (3.40282347e+38_f32 as f64)) as f32);
        }
    } else {
        m13 = ((z1 - z3) / (x1 - x3));
    }
    if (((abs((x2 - x3))) as f64) < (0.1_f64 * (eps as f64))) {
        if (((abs((xs - x2))) as f64) < (0.1_f64 * (eps as f64))) {
            return 1_i32;
        } else {
            m23 = ((0.01_f64 * (3.40282347e+38_f32 as f64)) as f32);
        }
    } else {
        m23 = ((z2 - z3) / (x2 - x3));
    }
    b12 = (z1 - (m12 * x1));
    b13 = (z1 - (m13 * x1));
    b23 = (z2 - (m23 * x2));
    /* source on edge? */
    if (((((abs(((zs - (m12 * xs)) - b12))) as f64) < (0.1_f64 * (eps as f64)))
        || (((abs(((zs - (m13 * xs)) - b13))) as f64) < (0.1_f64 * (eps as f64))))
        || (((abs(((zs - (m23 * xs)) - b23))) as f64) < (0.1_f64 * (eps as f64))))
    {
        return 1_i32;
    }
    return 0_i32;
}
unsafe fn writeFresnel(
    mut m: *mut Model,
    mut rayfp: *mut CFile,
    mut infofp: *mut CFile,
    mut nxz: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
    mut re: *mut RayEnd,
    mut krecord: i32,
    mut freq: f32,
    mut prim: i32,
    mut fresnelfp: *mut CFile,
) {
    let mut iray: i32 = 0;
    let mut nrs: i32 = 0;
    let mut irs: i32 = 0;
    let mut nxzw: i32 = 0;
    let mut ns: i32 = 0;
    let mut is: i32 = 0;
    let mut icount: i32 = 0;
    let mut sigmamax: f32 = (0.0_f64 as f32);
    let mut sigma1: f32 = 0.0;
    let mut sigma2: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut dsdx: f32 = 0.0;
    let mut dsdz: f32 = 0.0;
    let mut dummy: f32 = 0.0;
    let mut ds: f32 = 0.0;
    let mut xs: f32 = 0.0;
    let mut zs: f32 = 0.0;
    let mut dsigma: f32 = 0.0;
    let mut mfa: f32 = 0.0;
    let mut mfb: f32 = 0.0;
    let mut scale: f32 = 0.0;
    let mut q1end: f32 = (0.0_f64 as f32);
    let mut q2end: f32 = (0.0_f64 as f32);
    let mut q1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut s0: f32 = 0.0;
    let mut s1: f32 = 0.0;
    let mut s2: f32 = 0.0;
    let mut q1i: f32 = 0.0;
    let mut q2i: f32 = 0.0;
    let mut p1i: f32 = 0.0;
    let mut p2i: f32 = 0.0;
    let mut xfresnel1: f32 = 0.0;
    let mut xfresnel2: f32 = 0.0;
    let mut zfresnel1: f32 = 0.0;
    let mut zfresnel2: f32 = 0.0;
    let mut xa: *mut f32 = std::ptr::null_mut();
    let mut za: *mut f32 = std::ptr::null_mut();
    let mut r: *mut f32 = std::ptr::null_mut();
    let mut rsi: *mut RayStep = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    /* allocate space for auxiliary arrays */
    xa = ({
        let p = calloc((nxz as usize) as usize, 4);
        if p.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        p
    } as *mut f32);
    za = ({
        let p = calloc((nxz as usize) as usize, 4);
        if p.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        p
    } as *mut f32);
    r = ({
        let p = calloc((nxz as usize) as usize, 4);
        if p.is_null() {
            syserr!("%s: malloc failed", "ealloc.c");
        }
        p
    } as *mut f32);
    /* initialize counter */
    icount = 0_i32;
    /* loop over rays */
    {
        iray = 0_i32;
        while (iray < nray) {
            if ((krecord == 2147483647_i32)
                || ((krecord == (*(re).offset((iray) as isize)).kend)
                    && ((prim == 2147483647_i32)
                        || (prim == (*(re).offset((iray) as isize)).nref))))
            {
                /* count rays */
                icount = icount + 1;
                /* count ray steps while determining maximum sigma */
                /* and the spreading values at the ray ends */
                {
                    nrs = 0_i32;
                    rsi = (*(rs).offset((iray) as isize));
                    while (rsi != ((std::ptr::null_mut()) as *mut RayStep)) {
                        sigmamax = (*(rsi)).sigma;
                        q1end = (*(rsi)).q1;
                        q2end = (*(rsi)).q2;
                        nrs = nrs + 1;
                        rsi = (*(rsi)).rsNext;
                    }
                }
                /* initialize number of (x,z) written for this ray */
                nxzw = 0_i32;
                /* loop over ray steps */
                {
                    irs = 0_i32;
                    rsi = (*(rs).offset((iray) as isize));
                    while (irs < nrs) {
                        /* if sufficient number of (x,z) */
                        /* have been written */
                        if (nxzw >= nxz) {
                            break;
                        }
                        /* ray parameters at this step */
                        x = (*(rsi)).x;
                        z = (*(rsi)).z;
                        px = (*(rsi)).px;
                        pz = (*(rsi)).pz;
                        fa = ((*((*(rsi)).f)).fa as *mut FaceAttributes);
                        dsdx = (*(fa)).dsdx;
                        dsdz = (*(fa)).dsdz;
                        /* constants for nontrivial case */
                        s0 = ((px * px) + (pz * pz));
                        s1 = ((px * dsdx) + (pz * dsdz));
                        s2 = ((0.25_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64)) as f32);
                        /* sigma at this step and next step */
                        sigma1 = (*(rsi)).sigma;
                        sigma2 = (if (irs < (nrs - 1_i32)) {
                            (*((*(rsi)).rsNext)).sigma
                        } else {
                            sigma1
                        });
                        /* q1, q2, p1, p2 at this step */
                        q1i = (*(rsi)).q1;
                        q2i = (*(rsi)).q2;
                        p1i = (*(rsi)).p1;
                        p2i = (*(rsi)).p2;
                        /* number of sigma between this */
                        /* step and next */
                        ns = (((1_i32 as f64)
                            + ((((sigma2 - sigma1) * ((nxz - 1_i32) as f32)) as f64)
                                / (2.0_f64 * (sigmamax as f64))))
                            as i32);
                        /* increment in sigma */
                        ds = ((sigma2 - sigma1) / (ns as f32));
                        /* loop over sigma */
                        {
                            is = 0_i32;
                            dsigma = (0.0_f64 as f32);
                            while (is < ns) {
                                /* compute and write x and z */
                                xs = (((x as f64)
                                    + ((dsigma as f64)
                                        * ((px as f64)
                                            + (((dsigma as f64) * 0.25_f64) * (dsdx as f64)))))
                                    as f32);
                                zs = (((z as f64)
                                    + ((dsigma as f64)
                                        * ((pz as f64)
                                            + (((dsigma as f64) * 0.25_f64) * (dsdz as f64)))))
                                    as f32);
                                (*(xa).offset((nxzw) as isize)) = xs;
                                (*(za).offset((nxzw) as isize)) = zs;
                                if (efwrite(&[zs], &*(rayfp)) != (1_i32 as usize)) {
                                    err!("error writing ray file!\n");
                                }
                                if (efwrite(&[xs], &*(rayfp)) != (1_i32 as usize)) {
                                    err!("error writing ray file!\n");
                                }
                                q1 = q1i;
                                q2 = q2i;
                                p1 = p1i;
                                p2 = p2i;
                                /* compute and write p1,p2,q1,q2 */
                                evaluateDynamic2(
                                    dsigma,
                                    px,
                                    pz,
                                    dsdx,
                                    dsdz,
                                    (&raw mut q1),
                                    (&raw mut p1),
                                    (&raw mut q2),
                                    (&raw mut p2),
                                    s0,
                                    s1,
                                    s2,
                                );
                                /* compute m(Of,A) */
                                if (q2 == (0_i32 as f32)) {
                                    mfa = 3.40282347e+38_f32;
                                } else {
                                    mfa = (p2 / q2);
                                }
                                /* compute m(Of,B) */
                                mfb = (((-q1) * q2end) + (q2 * q1end));
                                if (mfb == (0_i32 as f32)) {
                                    mfb = 3.40282347e+38_f32;
                                } else {
                                    mfb = (((p2 * q1end) - (p1 * q2end)) / mfb);
                                }
                                /* compute fresnel radius */
                                dummy = ((abs((1.0_f64 / (((mfa - mfb) * freq) as f64)))) as f32);
                                (*(r).offset((nxzw) as isize)) = ((dummy as f64).sqrt() as f32);
                                /* if sufficient number */
                                /* written, break */
                                if ({
                                    let old = nxzw;
                                    nxzw = nxzw + 1;
                                    nxzw
                                } >= nxz)
                                {
                                    break;
                                }
                                is = is + 1;
                                dsigma = ((dsigma as f32) + (ds as f32)) as f32;
                            }
                        }
                        irs = irs + 1;
                        rsi = (*(rsi)).rsNext;
                    }
                }
                /* if necessary, repeat last (x,z) and r */
                while (nxzw < nxz) {
                    (*(xa).offset((nxzw) as isize)) = xs;
                    (*(za).offset((nxzw) as isize)) = zs;
                    (*(r).offset((nxzw) as isize)) = (0_i32 as f32);
                    if (efwrite(&[zs], &*(rayfp)) != (1_i32 as usize)) {
                        err!("error writing ray file!\n");
                    }
                    if (efwrite(&[xs], &*(rayfp)) != (1_i32 as usize)) {
                        err!("error writing ray file!\n");
                    }
                    nxzw = nxzw + 1;
                }
                /* compute fresnel points */
                if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                    fprintf!(&*(infofp), "Fresnel points:\n");
                }
                {
                    is = 0_i32;
                    while (is < nxz) {
                        if ((is == 0_i32) || (is == (nxz - 1_i32))) {
                            xfresnel1 = {
                                xfresnel2 = (*(xa).offset((is) as isize));
                                xfresnel2
                            };
                            zfresnel1 = {
                                zfresnel2 = (*(za).offset((is) as isize));
                                zfresnel2
                            };
                        } else {
                            scale = ((((((*(za).offset((is - 1_i32) as isize))
                                - (*(za).offset((is + 1_i32) as isize)))
                                * ((*(za).offset((is - 1_i32) as isize))
                                    - (*(za).offset((is + 1_i32) as isize))))
                                + (((*(xa).offset((is + 1_i32) as isize))
                                    - (*(xa).offset((is - 1_i32) as isize)))
                                    * ((*(xa).offset((is + 1_i32) as isize))
                                        - (*(xa).offset((is - 1_i32) as isize)))))
                                as f64)
                                .sqrt() as f32);
                            if (((scale == (0_i32 as f32)) || (is == 0_i32))
                                || (is == (nxz - 1_i32)))
                            {
                                xfresnel1 = {
                                    xfresnel2 = (*(xa).offset((is) as isize));
                                    xfresnel2
                                };
                                zfresnel1 = {
                                    zfresnel2 = (*(za).offset((is) as isize));
                                    zfresnel2
                                };
                            } else {
                                dummy = ((((1.0_f64 / (scale as f64))
                                    * (((*(za).offset((is - 1_i32) as isize))
                                        - (*(za).offset((is + 1_i32) as isize)))
                                        as f64))
                                    * ((*(r).offset((is) as isize)) as f64))
                                    as f32);
                                xfresnel1 = ((*(xa).offset((is) as isize)) + dummy);
                                xfresnel2 = ((*(xa).offset((is) as isize)) - dummy);
                                dummy = ((((1.0_f64 / (scale as f64))
                                    * (((*(xa).offset((is + 1_i32) as isize))
                                        - (*(xa).offset((is - 1_i32) as isize)))
                                        as f64))
                                    * ((*(r).offset((is) as isize)) as f64))
                                    as f32);
                                zfresnel1 = ((*(za).offset((is) as isize)) + dummy);
                                zfresnel2 = ((*(za).offset((is) as isize)) - dummy);
                            }
                        }
                        if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
                            fprintf!(
                                &*(infofp),
                                "x1=%g  z1=%g  x2=%g  z2=%g  r=%g\n",
                                (xfresnel1 as f64),
                                (zfresnel1 as f64),
                                (xfresnel2 as f64),
                                (zfresnel2 as f64),
                                ((*(r).offset((is) as isize)) as f64)
                            );
                        }
                        if (efwrite(&[zfresnel1], &*(fresnelfp)) != (1_i32 as usize)) {
                            err!("error writing Fresnel file!\n");
                        }
                        if (efwrite(&[xfresnel1], &*(fresnelfp)) != (1_i32 as usize)) {
                            err!("error writing Fresnel file!\n");
                        }
                        if (efwrite(&[zfresnel2], &*(fresnelfp)) != (1_i32 as usize)) {
                            err!("error writing Fresnel file!\n");
                        }
                        if (efwrite(&[xfresnel2], &*(fresnelfp)) != (1_i32 as usize)) {
                            err!("error writing Fresnel file!\n");
                        }
                        is = is + 1;
                    }
                }
            } else {
                /* do not write rays and fresnel volumes; */
                /* manipulate rayend information */
                (*(re).offset((iray) as isize)).kend = (-1_i32);
            }
            iray = iray + 1;
        }
    }
    /* close loop over rays */
    if (infofp != ((std::ptr::null_mut()) as *mut CFile)) {
        fprintf!(&*(infofp), "\nWrote %d rays to ray file.\n", icount);
    }
    /* free space */
    free((xa) as *mut c_void);
    free((za) as *mut c_void);
    free((r) as *mut c_void);
}
