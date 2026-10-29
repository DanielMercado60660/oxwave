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

//! SUTETRARAY - 3-D tetrahedral wavefront construction raytracing.
//! Port of tetra/main/sutetraray.c (Revision 1.5).
//! Credits: CWP: Zhaobo Meng, 1996.
//! Reference: Zhaobo Meng and Norman Bleistein, CWP report 251, 1997.

#![allow(
    unused_parens,
    unused_mut,
    unused_variables,
    unused_assignments,
    dead_code,
    non_snake_case
)]
use oxwave_core::prelude::*;
use oxwave_lib_tetra::*;

// Explicit casts preserve C's intermediate float/double types. Raw pointers
// retain the original array layout, linked triangles, and aliased ray fields;
// their allocations and streams belong to State, owned by main.
// checked_div/checked_rem preserve arm64 C's zero result for integer division
// by zero, rather than Rust's panic.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUTETRARAY - 3-D TETRAhedral WaveFront construction RAYtracing\t",
    "\t\t\t\t\t\t\t\t\t",
    "   sutetray < tetrafile [parameters] > ttfile\t \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nxgd=\t\t\t\tnumber of x samples (fastest dimension)\t",
    " nygd=\t\t\t\tnumber of y samples (second dimension) \t",
    " fxgd=xmin(from tetrafile)\tfirst x sample\t\t\t  \t",
    " fygd=ymin(from tetrafile)\tfirst y sample\t\t\t\t",
    " shootupward=1\t\t\t=1 shooting upward, -1 downward\t\t",
    " dxgd=(xmax-xmin)/(nxgd-1)\tx sampling interval\t\t \t",
    " dygd=(ymax-ymin)/(nygd-1)\ty sampling interval\t\t\t",
    " takeoff=30.0\t\t\tmax takeoff angle for shooting aperture\t",
    " ntakeoff=5\t\t\tnumber of samples in takeoff angle\t",
    " nazimuth=15\t\t\tnumber of samples in azimuth\t\t",
    " maxntris=1000\t\t\tmaximum number of triangles allowed\t",
    " maxnrays=600\t\t\tmaximum number of rays allowed\t\t",
    " maxnsegs=100\t\t\tmaximum number of ray segments allowed \t",
    " tmax=10\t\t\tmaximum traveltime traced\t\t",
    " dtwf=0.05\t\t\twavefront step size for raytracing\t",
    " ntwf=(int)(tmax/dtwf)+1\tmaximum number of wavefronts\t  \t",
    " ntwf=MAX(ntwf,ifwf2dump+nwf2dump)\t\t\t\t     \t",
    " edgemax=4\t\t\tmaximum triangle length (in dxgd) not to",
    "\t\t\t\tbe split\t\t\t\t",
    " fsx=fxgd+nxgd/2*dxgd\t\tfirst source in x\t\t\t",
    " dsx=dxgd\t\t\tx increment in source\t\t\t",
    " nsx=1\t\t\t\tnumber of source in x\t\t\t",
    " fsy=fygd+nygd/2*dygd\t\tfirst source in y\t\t\t",
    " dsy=dygd\t\t\ty increment in source\t\t   \t",
    " nsy=1\t\t\t\tnumber of source in y\t\t   \t",
    " fsz=zmin+(zmax-zmin)/2\tfirst source in z\t\t\t",
    " dsy=(zmax-zmin)/2\t\ty increment in source\t\t   \t",
    " nsy=1\t\t\t\tnumber of source in y\t\t\t",
    " nxt=nxgd*2/2+1\tnumber of x-samples in ttable for sukdmig3d\t",
    " nyt=nygd*2/2+1\tnumber of y-samples in ttable for sukdmig3d\t",
    " irefseq=nhz-1,...     reflector sequence\t\t\t\t",
    " crfile=NULL\t   file for cosines and ray paths r (for sukdmig3d)\t",
    " sttfile=NULL\t  file for surface traveltimes (for visualization)\t",
    " verbose=0\t\t=1 print some useful information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tThe following two files are output for viewer3 to display   \t",
    " rayfile=NULL\t\tray path file\t\t    \t\t\t",
    " wffile=NULL\t\twavefront file (then must give ifwf2dump)\t",
    " ifwf2dump=20\t\tthe first wavefront to dump  \t\t\t",
    " nwf2dump=1\t\tthe number of wavefronts to dump\t\t",
    " jpfile=stderr\t\tjob print (default: stderr)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Disclaimer:\t\t\t\t\t\t\t\t",
    " This is a research code that will take considerable work to get into  ",
    " the form of a a production level 3D migration code. The code is       ",
    " offered as is, along with tetramod and sukdmig3d, to provide a starting",
    " point for researchers who wish to write their own 3D migration codes.",
    " \t\t\t\t\t\t\t\t\t",
];

#[repr(C)]
#[derive(Clone, Copy)]
struct POINT {
    x: [f32; 3],
    s: f32,
    n: [f32; 3],
}

#[repr(C)]
#[derive(Clone, Copy)]
struct FACET {
    ip: [i32; 3],
    cn: [f32; 3],
    itetra: [i32; 2],
    area: f32,
    ct: [f32; 10],
}

#[repr(C)]
#[derive(Clone, Copy)]
struct TETRA {
    ip: [i32; 4],
    ifacet: [i32; 4],
    ireg: i32,
    gs: [f32; 3],
    v: f32,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct WF {
    x: [f32; 3],
    p: [f32; 3],
}

#[repr(C)]
#[derive(Clone, Copy)]
struct TRI {
    i1: i32,
    i2: i32,
    i3: i32,
    next: *mut TRI,
}

#[repr(C)]
#[derive(Clone, Copy)]
struct RAY {
    icode: i32,
    itetra: i32,
    lastfacet: i32,
    ionce: i32,
    iref: i32,
    shootupward: i32,
    nseg: i32,
    t: f32,
    x: [f32; 3],
    ttotal: f32,
    r: f32,
    v: f32,
    xseg: *mut *mut f32,
    tseg: *mut f32,
    p: [f32; 3],
}

const NORMAL: i32 = 0;

const SPOILED: i32 = 1;

const HIT_SURFACE: i32 = 2;

const HIT_FACET_FIRST: i32 = 0;

const HIT_WF_FIRST: i32 = 1;

const HIT_SURFACE_FIRST: i32 = 2;

const OUT_OF_MODEL: i32 = 3;

const PLEASE: i32 = 0;

const PLEASE_DO_NOT: i32 = 1;

const ALREADY: i32 = 0;

const NOT_YET: i32 = 1;

const WF_YES: i32 = 0;

const WF_NO: i32 = 1;

const WF_TRUE: i32 = 0;

const WF_FALSE: i32 = 1;

// main (sutetraray.c:236).

unsafe fn run(state: *mut State) -> i32 {
    let mut shootupward: i32 = 0;
    /* flag: 1 upward, -1 downward */
    let mut npoint: i32 = 0;
    /* number of control points */
    let mut point: *mut POINT = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut POINT);
    /* control points */
    let mut nfacet: i32 = 0;
    /* number of facets */
    let mut facet: *mut FACET = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut FACET);
    /* facet information */
    let mut ntetra: i32 = 0;
    /* number of tetra */
    let mut tetra: *mut TETRA = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TETRA);
    /* tetra information */
    let mut nygd: i32 = 0;
    /* number of samples in y direction */
    let mut nxgd: i32 = 0;
    /* number of samples in x direction */
    /**********************************************************
    These four parameters are used for defining the ttable
    size for sukdmig3d
    **********************************************************/
    let mut nxt: i32 = 0;
    /* number of x-samples in ttable for sukdmig3d */
    let mut nyt: i32 = 0;
    /* number of y-samples in ttable for sukdmig3d */
    let mut ixtf: i32 = 0;
    /* number of samples in ttable on left of source */
    let mut iytf: i32 = 0;
    /* number of samples in ttable in front of source */
    let mut ntakeoff: i32 = 0;
    /* number of samples in take-off */
    let mut nazimuth: i32 = 0;
    /* number of samples in azimuth */
    let mut fygd: f32 = 0.0;
    /* first sample in y */
    let mut fxgd: f32 = 0.0;
    /* first sample in x */
    let mut fsx: f32 = 0.0;
    /* first source in x */
    let mut fsy: f32 = 0.0;
    /* first source in y */
    let mut fsz: f32 = 0.0;
    /* first source in z */
    let mut dsx: f32 = 0.0;
    /* x increment of sources */
    let mut dsy: f32 = 0.0;
    /* y increment of sources */
    let mut dsz: f32 = 0.0;
    /* z increment of sources */
    let mut nsx: i32 = 0;
    /* number of sources in x direction */
    let mut nsy: i32 = 0;
    /* number of sources in y direction */
    let mut nsz: i32 = 0;
    /* number of sources in z direction */
    let mut irefseq: *mut i32 = std::ptr::null_mut();
    /* reflector sequence */
    let mut dygd: f32 = 0.0;
    /* y sample interval */
    let mut dxgd: f32 = 0.0;
    /* x sample interval */
    let mut dtwf: f32 = 0.0;
    /* time step in raytracing */
    let mut ntwf: i32 = 0;
    /* number of time steps in raytracing */
    let mut tmax: f32 = 0.0;
    /* maximum time traced */
    let mut rayfile: *mut i8 = (concat!("", "\0").as_ptr() as *mut i8);
    /* file for ray information */
    let mut wffile: *mut i8 = (concat!("", "\0").as_ptr() as *mut i8);

    let mut crfile: *mut i8 = (concat!("", "\0").as_ptr() as *mut i8);

    let mut sttfile: *mut i8 = (concat!("", "\0").as_ptr() as *mut i8);

    let mut jpfile: *mut i8 = (concat!("", "\0").as_ptr() as *mut i8);
    /* file for job print */
    let mut rayfp: *mut CFile = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile);

    let mut rayfp0: *mut CFile = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile);

    let mut crfp: *mut CFile = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile);

    let mut sttfp: *mut CFile = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile);

    let mut wffp: *mut CFile = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile);

    let mut jpfp: *mut CFile = std::ptr::addr_of_mut!((*state).errors);

    // The C variables are uninitialized unless wffile is open. Zero gives
    // deterministic storage for the otherwise undefined no-wffile path.
    let mut ifwf2dump: i32 = 0;

    let mut nwf2dump: i32 = 0;

    let mut nrefseq: i32 = 0;

    let mut ns: i32 = 0;
    /* number of sources */
    let mut maxntris: i32 = 0;
    /* maximum number of triangles allowed */
    let mut maxnrays: i32 = 0;
    /* maximum number of rays allowed */
    let mut edgemax: f32 = 0.0;
    /* maximum edge length in dxgd */
    let mut verbose: i32 = 0;
    /* if =1 print some useful information */
    let mut source: [f32; 3] = [0.0; 3];
    /* source position */
    let mut isx: i32 = 0;
    /* source index for nsx */
    let mut isy: i32 = 0;
    /* source index for nsy */
    let mut isz: i32 = 0;
    /* source index for nsz */
    let mut ixsf: i32 = 0;
    /* first source in gd */
    let mut ixsr: i32 = 0;
    /* source spacing and gd spacing ratio */
    let mut iysf: i32 = 0;
    /* first source in gd */
    let mut iysr: i32 = 0;
    /* source spacing and gd spacing ratio */
    let mut xmin: f32 = 0.0;
    let mut xmax: f32 = 0.0;
    let mut ymin: f32 = 0.0;
    let mut ymax: f32 = 0.0;
    let mut zmin: f32 = 0.0;
    let mut zmax: f32 = 0.0;

    let mut takeoff: f32 = 0.0;
    /* max takeoff angle */
    let mut iref: i32 = 0;
    let mut maxnsegs: i32 = 0;
    /* hook up getpar */

    requestdoc(1i32, SDOC);
    /* get parameters */
    if (!(par_int("verbose", std::ptr::addr_of_mut!(verbose)) != 0)) {
        verbose = 0i32;
    }

    if (!(par_int("shootupward", std::ptr::addr_of_mut!(shootupward)) != 0)) {
        shootupward = 1i32;
    }

    if (par_string(state, "jpfile", std::ptr::addr_of_mut!(jpfile)) != 0) {
        if (({
            jpfp = open_file(state, jpfile, (concat!("w", "\0").as_ptr() as *mut i8));
            jpfp
        }) == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
        {
            jpfp = std::ptr::addr_of_mut!((*state).errors);
        }
    }

    if (!(getparfloat("takeoff", &mut *std::ptr::addr_of_mut!(takeoff)) != 0)) {
        takeoff = (30f64 as f32);
    }

    takeoff = ((((takeoff as f64) / 180f64) * (3.1415926535897931f64)) as f32);

    nrefseq = countparval("irefseq");

    irefseq = allocate(
        state,
        ((nrefseq as usize) as usize) * std::mem::size_of::<i32>(),
    )
    .cast::<i32>();

    par_int("irefseq", irefseq);

    if (verbose != 0) {
        {
            iref = 0i32;
            loop {
                if !(iref < nrefseq) {
                    break;
                }
                fprintf!(
                    &*jpfp,
                    "irefseq[%d]=%d\n",
                    iref,
                    (*(irefseq).offset((iref) as isize))
                );
                {
                    let old = iref;
                    iref = iref + 1;
                    old
                };
            }
        }
    }

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(xmin) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "=xmin\n") {
                break 'scan;
            }
        }
        scanned
    };

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(xmax) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "=xmax\n") {
                break 'scan;
            }
        }
        scanned
    };

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(ymin) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "=ymin\n") {
                break 'scan;
            }
        }
        scanned
    };

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(ymax) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "=ymax\n") {
                break 'scan;
            }
        }
        scanned
    };

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(zmin) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "=zmin\n") {
                break 'scan;
            }
        }
        scanned
    };

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(zmax) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "=zmax\n") {
                break 'scan;
            }
        }
        scanned
    };

    if (verbose != 0) {
        fprintf!(
            &*jpfp,
            "xmin=%f\nxmax=%f\nymin=%f\nymax=%f\nzmin=%f\nzmax=%f\n",
            (xmin as f64),
            (xmax as f64),
            (ymin as f64),
            (ymax as f64),
            (zmin as f64),
            (zmax as f64)
        );
    }

    if (!(par_int("nxgd", std::ptr::addr_of_mut!(nxgd)) != 0)) {
        nxgd = 21i32;
    }

    if (nxgd < 2i32) {
        err!("nxgd too small");
    }

    nxgd = (((nxgd).checked_div(2i32).unwrap_or(0)).wrapping_mul(2i32)).wrapping_add(1i32);

    if (!(getparfloat("fxgd", &mut *std::ptr::addr_of_mut!(fxgd)) != 0)) {
        fxgd = xmin;
    }

    if (fxgd < xmin) {
        err!("fxgd too small");
    }

    if (!(getparfloat("dxgd", &mut *std::ptr::addr_of_mut!(dxgd)) != 0)) {
        dxgd = ((xmax - xmin) / (((nxgd).wrapping_sub(1i32)) as f32));
    }

    if (((fxgd + ((((nxgd).wrapping_sub(1i32)) as f32) * dxgd)) as f64) > ((xmax as f64) * 1.01f64))
    {
        fprintf!(
            &*std::ptr::addr_of_mut!((*state).errors),
            "fxgd=%f nxgd=%d dxgd=%f xmax=%f\n",
            (fxgd as f64),
            nxgd,
            (dxgd as f64),
            (xmax as f64)
        );

        err!("fxgd+(nxgd-1)*dxgd too large");
    }

    if (!(par_int("nygd", std::ptr::addr_of_mut!(nygd)) != 0)) {
        nygd = 21i32;
    }

    if (nygd < 2i32) {
        err!("nygd too small");
    }

    nygd = (((nygd).checked_div(2i32).unwrap_or(0)).wrapping_mul(2i32)).wrapping_add(1i32);

    if (!(getparfloat("fygd", &mut *std::ptr::addr_of_mut!(fygd)) != 0)) {
        fygd = ymin;
    }

    if (fygd < ymin) {
        err!("fygd too small");
    }

    if (!(getparfloat("dygd", &mut *std::ptr::addr_of_mut!(dygd)) != 0)) {
        dygd = ((ymax - ymin) / (((nygd).wrapping_sub(1i32)) as f32));
    }

    if (((fygd + ((((nygd).wrapping_sub(1i32)) as f32) * dygd)) as f64) > ((ymax as f64) * 1.01f64))
    {
        err!("fygd+(nygd-1)*dygd too large");
    }

    if (verbose != 0) {
        fprintf!(&*jpfp, "shootupward=%d\n", shootupward);

        fprintf!(&*jpfp, "nxgd=%d\nnygd=%d\n", nxgd, nygd);

        fprintf!(&*jpfp, "dxgd=%f\ndygd=%f\n", (dxgd as f64), (dygd as f64));

        fprintf!(&*jpfp, "fxgd=%f\nfygd=%f\n", (fxgd as f64), (fygd as f64));
    }

    if (par_string(state, "rayfile", std::ptr::addr_of_mut!(rayfile)) != 0) {
        if (({
            rayfp = open_file(state, rayfile, (concat!("w", "\0").as_ptr() as *mut i8));
            rayfp
        }) == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
        {
            err!("Can not open rayfile to write");
        }
    }

    if (par_string(state, "wffile", std::ptr::addr_of_mut!(wffile)) != 0) {
        if (({
            wffp = open_file(state, wffile, (concat!("w", "\0").as_ptr() as *mut i8));
            wffp
        }) == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
        {
            err!("Can not open wffile to write");
        }
    }

    if (wffp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile)) {
        if (!(par_int("ifwf2dump", std::ptr::addr_of_mut!(ifwf2dump)) != 0)) {
            ifwf2dump = 20i32;
        }

        if (!(par_int("nwf2dump", std::ptr::addr_of_mut!(nwf2dump)) != 0)) {
            nwf2dump = 1i32;
        }
    }

    if (verbose != 0) {
        fprintf!(&*jpfp, "ifwf2dump=%d, nwf2dump=%d\n", ifwf2dump, nwf2dump);
    }

    if (par_string(state, "crfile", std::ptr::addr_of_mut!(crfile)) != 0) {
        if (({
            crfp = open_file(state, crfile, (concat!("w", "\0").as_ptr() as *mut i8));
            crfp
        }) == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
        {
            err!("Can not open crfile to write");
        }
    }

    if (par_string(state, "sttfile", std::ptr::addr_of_mut!(sttfile)) != 0) {
        if (({
            sttfp = open_file(state, sttfile, (concat!("w", "\0").as_ptr() as *mut i8));
            sttfp
        }) == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
        {
            err!("Can not open sttfile to write");
        }
    }

    if (!(par_int("maxntris", std::ptr::addr_of_mut!(maxntris)) != 0)) {
        maxntris = 1000i32;
    }

    if (!(par_int("maxnrays", std::ptr::addr_of_mut!(maxnrays)) != 0)) {
        maxnrays = 600i32;
    }

    if (!(par_int("maxnsegs", std::ptr::addr_of_mut!(maxnsegs)) != 0)) {
        maxnsegs = 100i32;
    }

    if (!(getparfloat("edgemax", &mut *std::ptr::addr_of_mut!(edgemax)) != 0)) {
        edgemax = (4i32 as f32);
    }

    if (!(getparfloat("fsx", &mut *std::ptr::addr_of_mut!(fsx)) != 0)) {
        fsx = (((fxgd as f64) + (((nxgd as f64) / 2f64) * (dxgd as f64))) as f32);
    }

    if (!(getparfloat("fsy", &mut *std::ptr::addr_of_mut!(fsy)) != 0)) {
        fsy = (((fygd as f64) + (((nygd as f64) / 2f64) * (dygd as f64))) as f32);
    }

    if (!(par_int("nsx", std::ptr::addr_of_mut!(nsx)) != 0)) {
        nsx = 1i32;
    }

    if (!(par_int("nsy", std::ptr::addr_of_mut!(nsy)) != 0)) {
        nsy = 1i32;
    }

    if (!(par_int("nsz", std::ptr::addr_of_mut!(nsz)) != 0)) {
        nsz = 1i32;
    }

    ns = (nsx).wrapping_mul(nsy);

    if (!(par_int("nxt", std::ptr::addr_of_mut!(nxt)) != 0)) {
        nxt = (((nxgd).checked_div(2i32).unwrap_or(0)).wrapping_mul(2i32)).wrapping_add(1i32);
    }

    if (!(par_int("nyt", std::ptr::addr_of_mut!(nyt)) != 0)) {
        nyt = (((nygd).checked_div(2i32).unwrap_or(0)).wrapping_mul(2i32)).wrapping_add(1i32);
    }

    nxt = (((nxt).checked_div(2i32).unwrap_or(0)).wrapping_mul(2i32)).wrapping_add(1i32);

    nyt = (((nyt).checked_div(2i32).unwrap_or(0)).wrapping_mul(2i32)).wrapping_add(1i32);

    if (!(getparfloat("dsx", &mut *std::ptr::addr_of_mut!(dsx)) != 0)) {
        dsx = dxgd;
    }

    if (!(getparfloat("dsy", &mut *std::ptr::addr_of_mut!(dsy)) != 0)) {
        dsy = dygd;
    }

    if (!(getparfloat("dsz", &mut *std::ptr::addr_of_mut!(dsz)) != 0)) {
        dsz = ((zmax - zmin) / (2i32 as f32));
    }

    ixsf = (((fsx - fxgd) / dxgd) as i32);

    if (ixsf < (nxt).checked_div(2i32).unwrap_or(0)) {
        err!(
            "ixsf=%d should >= nxt/2=%d. Use larger fsx or smaller nxt",
            ixsf,
            (nxt).checked_div(2i32).unwrap_or(0)
        );
    }

    ixsr = ((if (((1i32) as f32) > (dsx / dxgd)) {
        ((1i32) as f32)
    } else {
        (dsx / dxgd)
    }) as i32);

    if (((ixsf).wrapping_add((ixsr).wrapping_mul(nsx)))
        .wrapping_add((nxt).checked_div(2i32).unwrap_or(0))
        > nxgd)
    {
        err!(
            "ixsf+ixsr*nsx+nxt/2=%d should <= nxgd=%d",
            ((ixsf).wrapping_add((ixsr).wrapping_mul(nsx)))
                .wrapping_add((nxt).checked_div(2i32).unwrap_or(0)),
            nxgd
        );
    }

    iysf = (((fsy - fygd) / dygd) as i32);

    if (iysf < (nyt).checked_div(2i32).unwrap_or(0)) {
        err!(
            "iysf=%d should >= nyt/2=%d. Use larger fsy or smaller nyt",
            iysf,
            (nyt).checked_div(2i32).unwrap_or(0)
        );
    }

    iysr = ((if (((1i32) as f32) > (dsy / dygd)) {
        ((1i32) as f32)
    } else {
        (dsy / dygd)
    }) as i32);

    if (((iysf).wrapping_add((iysr).wrapping_mul(nsy)))
        .wrapping_add((nyt).checked_div(2i32).unwrap_or(0))
        > nygd)
    {
        err!(
            "iysf+iysr*nsy+nyt/2=%d should <= nygd=%d",
            ((iysf).wrapping_add((iysr).wrapping_mul(nsy)))
                .wrapping_add((nyt).checked_div(2i32).unwrap_or(0)),
            nygd
        );
    }

    fsx = (fxgd + ((ixsf as f32) * dxgd));

    fsy = (fygd + ((iysf as f32) * dygd));

    fsx = (((fsx as f64) + (0.01f64 * (dxgd as f64))) as f32);

    fsy = (((fsy as f64) + (0.01f64 * (dygd as f64))) as f32);

    dsx = ((ixsr as f32) * dxgd);

    dsy = ((iysr as f32) * dygd);

    if (!(getparfloat("fsz", &mut *std::ptr::addr_of_mut!(fsz)) != 0)) {
        fsz = (((zmin as f64) + (((zmax - zmin) as f64) / 2f64)) as f32);
    }

    if (!(par_int("ntakeoff", std::ptr::addr_of_mut!(ntakeoff)) != 0)) {
        ntakeoff = 5i32;
    }

    if (!(par_int("nazimuth", std::ptr::addr_of_mut!(nazimuth)) != 0)) {
        nazimuth = 15i32;
    }

    if (!(getparfloat("tmax", &mut *std::ptr::addr_of_mut!(tmax)) != 0)) {
        tmax = (10f64 as f32);
    }

    if (!(getparfloat("dtwf", &mut *std::ptr::addr_of_mut!(dtwf)) != 0)) {
        dtwf = (0.050000000000000003f64 as f32);
    }

    if (!(par_int("ntwf", std::ptr::addr_of_mut!(ntwf)) != 0)) {
        ntwf = ((tmax / dtwf) as i32).wrapping_add(1i32);

        ntwf = (if ((200i32) > (ntwf)) {
            (200i32)
        } else {
            (ntwf)
        });
    }

    checkpars();

    if (ntwf < (ifwf2dump).wrapping_add(nwf2dump)) {
        fprintf!(
            &*jpfp,
            "ntwf < ifwf2dump+nwf2dump, reset=%d\n",
            (ifwf2dump).wrapping_add(nwf2dump)
        );

        ntwf = (ifwf2dump).wrapping_add(nwf2dump);
    }

    if (ntwf > 2000i32) {
        err!("ntwf>2000, too large");
    }

    if (verbose != 0) {
        fprintf!(
            &*jpfp,
            "nxgd=%d\nnygd=%d\ndxgd=%g\n",
            nxgd,
            nygd,
            (dxgd as f64)
        );

        fprintf!(
            &*jpfp,
            "dygd=%g\nntakeoff=%d\nnazimuth=%d\n",
            (dygd as f64),
            ntakeoff,
            nazimuth
        );

        fprintf!(&*jpfp, "takeoff=%f\n", (takeoff as f64));

        fprintf!(
            &*jpfp,
            "fsx=%g,dsx=%g,nsx=%d\nfsy=%g,dsy=%g,nsy=%d\nfsz=%g,dsz=%g,nsz=%d\n",
            (fsx as f64),
            (dsx as f64),
            nsx,
            (fsy as f64),
            (dsy as f64),
            nsy,
            (fsz as f64),
            (dsz as f64),
            nsz
        );

        fprintf!(&*jpfp, "nxt=%d\nnyt=%d\n", nxt, nyt);

        fprintf!(
            &*jpfp,
            "maxntris=%d\nmaxnrays=%d\nmaxnsegs=%d\nedgemax=%g\n",
            maxntris,
            maxnrays,
            maxnsegs,
            (edgemax as f64)
        );

        fprintf!(
            &*jpfp,
            "ntwf=%d,dtwf=%g,tmax=%g\n",
            ntwf,
            (dtwf as f64),
            (tmax as f64)
        );
    }

    if ((nxt < 0i32) || (nxt > nxgd)) {
        err!("nxt must be >=0 <=nxgd");
    }

    if ((nyt < 0i32) || (nyt > nygd)) {
        err!("nyt must be >=0 <=nygd");
    }

    if (((((fsx + ((((nsx).wrapping_sub(1i32)) as f32) * dsx))
        > (fxgd + ((((nxgd).wrapping_sub(1i32)) as f32) * dxgd)))
        || ((fsy + ((((nsy).wrapping_sub(1i32)) as f32) * dsy))
            > (fygd + ((((nygd).wrapping_sub(1i32)) as f32) * dygd))))
        || (fsx < fxgd))
        || (fsy < fygd))
    {
        err!("Source out of medium");
    }
    /* number of sources */
    (*state).tro.hdr.nhs = (nsx as i16);

    (*state).tro.hdr.nvs = (nsy as i16);

    (*state).tro.hdr.duse = (nsz as i16);
    /* first source */
    (*state).tro.hdr.sdel = ixsf;

    (*state).tro.hdr.gdel = iysf;

    (*state).tro.hdr.sdepth = (((fsz as f64) * 1000f64) as i32);
    /* source spacing */
    (*state).tro.hdr.swdep = ixsr;

    (*state).tro.hdr.gwdep = iysr;

    (*state).tro.hdr.ep = (((dsz as f64) * 1000f64) as i32);
    /* first sample */
    (*state).tro.hdr.f1 = fxgd;

    (*state).tro.hdr.f2 = fygd;
    /* model size: nxgd x nygd */
    (*state).tro.hdr.gelev = nygd;

    (*state).tro.hdr.selev = nxgd;
    /* traveltime sample spacing */
    (*state).tro.hdr.d1 = dxgd;

    (*state).tro.hdr.d2 = dygd;
    /* number of traces and samples in each trace */
    (*state).tro.hdr.ns = (nxt as u16);

    (*state).tro.hdr.ntr = nyt;

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(npoint) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(
                &*std::ptr::addr_of_mut!((*state).input),
                "=npoint : x[3] : s : n[3]\n",
            ) {
                break 'scan;
            }
        }
        scanned
    };

    fprintf!(&*jpfp, "npoint=%d\n", npoint);

    point = (allocate(
        state,
        (npoint as u64).wrapping_mul((std::mem::size_of::<POINT>() as u64)) as usize,
    ) as *mut POINT);

    read_point(state, point, npoint);

    fprintf!(&*jpfp, "point %d read\n", npoint);

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(nfacet) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(
                &*std::ptr::addr_of_mut!((*state).input),
                "=nfacet : area : ct[10] : ip[3] : cn[3] : itetra[2]\n",
            ) {
                break 'scan;
            }
        }
        scanned
    };

    fprintf!(&*jpfp, "nfacet=%d\n", nfacet);

    facet = (allocate(
        state,
        (nfacet as u64).wrapping_mul((std::mem::size_of::<FACET>() as u64)) as usize,
    ) as *mut FACET);

    read_facet(state, facet, nfacet);

    fprintf!(&*jpfp, "facet %d read\n", nfacet);

    {
        let mut scanned = 0i32;
        'scan: {
            if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                *std::ptr::addr_of_mut!(ntetra) = v;
                scanned += 1;
            } else {
                break 'scan;
            }
            if !scan_literal(
                &*std::ptr::addr_of_mut!((*state).input),
                "=ntetra : v : ip[4] : ifacet[4] : ireg : gs[3]\n",
            ) {
                break 'scan;
            }
        }
        scanned
    };

    fprintf!(&*jpfp, "ntetra=%d\n", ntetra);

    tetra = (allocate(
        state,
        (ntetra as u64).wrapping_mul((std::mem::size_of::<TETRA>() as u64)) as usize,
    ) as *mut TETRA);

    read_tetra(state, tetra, ntetra);

    fprintf!(&*jpfp, "tetra %d read\n", ntetra);

    if (rayfp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile)) {
        fprintf!(&*rayfp, "%d =Number of shots\n", ns);
    }

    {
        isy = 0i32;
        loop {
            if !(isy < nsy) {
                break;
            }
            {
                (*(std::ptr::addr_of_mut!(source).cast::<f32>()).offset((1i32) as isize)) =
                    (fsy + ((isy as f32) * dsy));

                iytf = ((iysf).wrapping_add((iysr).wrapping_mul(isy)))
                    .wrapping_sub((nyt).checked_div(2i32).unwrap_or(0));

                {
                    isx = 0i32;
                    loop {
                        if !(isx < nsx) {
                            break;
                        }
                        {
                            (*(std::ptr::addr_of_mut!(source).cast::<f32>())
                                .offset((0i32) as isize)) = (fsx + ((isx as f32) * dsx));

                            ixtf = ((ixsf).wrapping_add((ixsr).wrapping_mul(isx)))
                                .wrapping_sub((nxt).checked_div(2i32).unwrap_or(0));

                            {
                                isz = 0i32;
                                loop {
                                    if !(isz < nsz) {
                                        break;
                                    }
                                    {
                                        (*(std::ptr::addr_of_mut!(source).cast::<f32>())
                                            .offset((2i32) as isize)) =
                                            (fsz + ((isz as f32) * dsz));

                                        if (isz < (nsz).wrapping_sub(1i32)) {
                                            rayfp0 = ((std::ptr::null_mut::<std::ffi::c_void>())
                                                as *mut CFile);
                                        } else {
                                            rayfp0 = rayfp;
                                        }
                                        /*******************************************************
                                                    Usually it is not every good to put the source right on
                                                    a facet, this is because the raytracer may not stably
                                                    find the tetrahetron that the source lies in. So addition
                                        of EPS will help.
                                        ********************************************************/
                                        /* shooting one shot */
                                        oneshot(
                                            state,
                                            shootupward,
                                            maxnsegs,
                                            npoint,
                                            point,
                                            nfacet,
                                            facet,
                                            ntetra,
                                            tetra,
                                            irefseq,
                                            nrefseq,
                                            jpfp,
                                            takeoff,
                                            dtwf,
                                            ntwf,
                                            maxntris,
                                            maxnrays,
                                            edgemax,
                                            nxt,
                                            nyt,
                                            ixtf,
                                            iytf,
                                            std::ptr::addr_of_mut!(source).cast::<f32>(),
                                            ntakeoff,
                                            nazimuth,
                                            fxgd,
                                            dxgd,
                                            nxgd,
                                            fygd,
                                            dygd,
                                            nygd,
                                            xmin,
                                            ymin,
                                            xmax,
                                            ymax,
                                            zmax,
                                            rayfp0,
                                            wffp,
                                            crfp,
                                            sttfp,
                                            ifwf2dump,
                                            nwf2dump,
                                        );

                                        fprintf!(
                                            &*std::ptr::addr_of_mut!((*state).errors),
                                            "shot (%d %d %d) is written\n",
                                            isx,
                                            isy,
                                            isz
                                        );

                                        if (wffp
                                            != ((std::ptr::null_mut::<std::ffi::c_void>())
                                                as *mut CFile))
                                        {
                                            fclose(&*wffp);
                                        }

                                        wffp = ((std::ptr::null_mut::<std::ffi::c_void>())
                                            as *mut CFile);
                                    }
                                    {
                                        let old = isz;
                                        isz = isz + 1;
                                        old
                                    };
                                }
                            }
                        }
                        {
                            let old = isx;
                            isx = isx + 1;
                            old
                        };
                    }
                }
            }
            {
                let old = isy;
                isy = isy + 1;
                old
            };
        }
    }

    if (rayfp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile)) {
        fclose(&*rayfp);
    }

    if (jpfp != std::ptr::addr_of_mut!((*state).errors)) {
        fclose(&*jpfp);
    }

    return 0i32;
}

// WF_normalize (sutetraray.c:649).

unsafe fn WF_normalize(mut n1: *mut f32, mut n2: *mut f32, mut n3: *mut f32) {
    let mut norm: f32 = 0.0;

    norm = ((((*n1) * (*n1)) + ((*n2) * (*n2))) + ((*n3) * (*n3)));

    norm = ((norm as f64).sqrt() as f32);

    (*n1) = ((*n1) / norm);

    (*n2) = ((*n2) / norm);

    (*n3) = ((*n3) / norm);
}

// oneshot (sutetraray.c:670).

unsafe fn oneshot(
    state: *mut State,
    mut shootupward: i32,
    mut maxnsegs: i32,
    mut npoint: i32,
    mut point: *mut POINT,
    mut nfacet: i32,
    mut facet: *mut FACET,
    mut ntetra: i32,
    mut tetra: *mut TETRA,
    mut irefseq: *mut i32,
    mut nrefseq: i32,
    mut jpfp: *mut CFile,
    mut takeoff: f32,
    mut dtwf: f32,
    mut ntwf: i32,
    mut maxntris: i32,
    mut maxnrays: i32,
    mut edgemax: f32,
    mut nxt: i32,
    mut nyt: i32,
    mut ixtf: i32,
    mut iytf: i32,
    mut source: *mut f32,
    mut ntakeoff: i32,
    mut nazimuth: i32,
    mut fxgd: f32,
    mut dxgd: f32,
    mut nxgd: i32,
    mut fygd: f32,
    mut dygd: f32,
    mut nygd: i32,
    mut xmin: f32,
    mut ymin: f32,
    mut xmax: f32,
    mut ymax: f32,
    mut zmax: f32,
    mut rayfp: *mut CFile,
    mut wffp: *mut CFile,
    mut crfp: *mut CFile,
    mut sttfp: *mut CFile,
    mut ifwf2dump: i32,
    mut nwf2dump: i32,
) -> i32 {
    let mut wf: *mut WF = std::ptr::null_mut();
    /*******************************************************
    The traveltime table start from ixgdl, end with ixgdl+nxt-1 in x
    and start with iygdf and end with iygdf+nyt-1
    *******************************************************/
    let mut azimuth: f32 = 0.0;
    /* azimuth angle */
    let mut sintakeoff: f32 = 0.0;
    /* sin(takeoff) */
    let mut cosazimuth: f32 = 0.0;
    /* cos(azimuth) */
    let mut costakeoff: f32 = 0.0;
    /* cos(takeoff) */
    let mut sinazimuth: f32 = 0.0;
    /* sin(azimuth) */
    let mut s: f32 = 0.0;
    /* slowness */
    let mut itakeoff: i32 = 0;
    /* index for azimuth */
    let mut iazimuth: i32 = 0;
    /* index for takeoff */
    let mut dazimuth: f32 = 0.0;
    /* azimuth increment */
    let mut dtakeoff: f32 = 0.0;
    /* takeoff increment */
    let mut ithis: i32 = 0;
    /* index for this ray */
    let mut inext: i32 = 0;

    let mut itwf: i32 = 0;
    /* index for t */
    let mut distmax: f32 = (0i32 as f32);
    /* the maximum edge of the triangle */
    let mut iedgemax: i32 = 0;
    /* the index for the largest triangle */
    let mut iedgemax0: i32 = 0;

    let mut iedge: i32 = 0;
    /* index for iedgemax */
    let mut enough_ntwf: i32 = WF_TRUE;

    let mut i1: i32 = 0;
    let mut i2: i32 = 0;
    let mut i3: i32 = 0;

    let mut ii1: i32 = 0;
    let mut ii2: i32 = 0;
    let mut ii3: i32 = 0;

    let mut iray: i32 = 0;
    /* index for rays */
    let mut head: TRI = std::mem::zeroed();
    /* first triangle */
    let mut headp: *mut TRI = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI);
    /* pointer to the first usable triangle */
    let mut trip: *mut TRI = std::ptr::null_mut();
    /* pointer to tris */
    let mut trip0: *mut TRI = std::ptr::null_mut();
    /* runing pointer */
    let mut itrip: *mut TRI = std::ptr::null_mut();
    /* index for trip */
    let mut temtrip: *mut TRI = std::ptr::null_mut();
    /* temporary pointer to triangles */
    let mut nrays: i32 = 0i32;
    /* number of rays */
    let mut ngoodrays: i32 = 0;

    let mut ntris: i32 = 0i32;
    /* number of triangles */
    let mut dist: [f32; 3] = [0.0; 3];

    let mut ixgd: i32 = 0;
    let mut iygd: i32 = 0;

    let mut x0: f32 = 0.0;
    let mut x1: f32 = 0.0;
    let mut x2: f32 = 0.0;

    let mut y0: f32 = 0.0;
    let mut y1: f32 = 0.0;
    let mut y2: f32 = 0.0;

    let mut z0: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut z2: f32 = 0.0;

    let mut xvertex0: f32 = 0.0;
    let mut xvertex1: f32 = 0.0;

    let mut yvertex0: f32 = 0.0;
    let mut yvertex1: f32 = 0.0;

    let mut zvertex0: f32 = 0.0;
    let mut zvertex1: f32 = 0.0;

    let mut initreg: i32 = 0;

    let mut ntr: i32 = 0;

    let mut initetra: i32 = 0;

    let mut tripi0: i32 = 0;
    let mut tripi1: i32 = 0;

    let mut ttable: *mut *mut f32 = std::ptr::null_mut();
    /* traveltime table */
    let mut ctable: *mut *mut f32 = std::ptr::null_mut();

    let mut rtable: *mut *mut f32 = std::ptr::null_mut();

    let mut check_tris: i32 = PLEASE;

    let mut result: i32 = 0;

    let mut done: i32 = 0;

    let mut ray: *mut RAY = std::ptr::null_mut();

    let mut temp: f32 = 0.0;
    let mut sqrts: f32 = 0.0;
    let mut tsum: f32 = 0.0;

    let mut iseg: i32 = 0;

    let mut adjust: f32 = 0.0;
    /* when interpolating, ajust this amount */
    let mut r: f32 = 0.0;
    let mut d: f32 = 0.0;
    let mut d2: f32 = 0.0;
    let mut sintheta: f32 = 0.0;
    let mut costheta: f32 = 0.0;

    ray = (allocate(
        state,
        (maxnrays as u64).wrapping_mul((std::mem::size_of::<RAY>() as u64)) as usize,
    ) as *mut RAY);

    fprintf!(
        &*jpfp,
        "source at %f %f %f\n",
        ((*(source).offset((0i32) as isize)) as f64),
        ((*(source).offset((1i32) as isize)) as f64),
        ((*(source).offset((2i32) as isize)) as f64)
    );

    ttable = array2(state, (nxgd as usize), (nygd as usize));

    ctable = array2(state, (nxgd as usize), (nygd as usize));

    rtable = array2(state, (nxgd as usize), (nygd as usize));

    {
        ixgd = 0i32;
        loop {
            if !(ixgd < nxgd) {
                break;
            }
            {
                iygd = 0i32;
                loop {
                    if !(iygd < nygd) {
                        break;
                    }
                    (*(*(ttable).offset((iygd) as isize)).offset((ixgd) as isize)) = (9f64 as f32);
                    {
                        let old = iygd;
                        iygd = iygd + 1;
                        old
                    };
                }
            }
            {
                let old = ixgd;
                ixgd = ixgd + 1;
                old
            };
        }
    }

    {
        std::ptr::write_bytes(
            ((*(ctable).offset((0i32) as isize)) as *mut std::ffi::c_void),
            0i32 as u8,
            ((nxgd).wrapping_mul(nygd) as u64).wrapping_mul((std::mem::size_of::<f32>() as u64))
                as usize,
        );
        ((*(ctable).offset((0i32) as isize)) as *mut std::ffi::c_void)
    };

    {
        std::ptr::write_bytes(
            ((*(rtable).offset((0i32) as isize)) as *mut std::ffi::c_void),
            0i32 as u8,
            ((nxgd).wrapping_mul(nygd) as u64).wrapping_mul((std::mem::size_of::<f32>() as u64))
                as usize,
        );
        ((*(rtable).offset((0i32) as isize)) as *mut std::ffi::c_void)
    };
    /* the ray information is bookkept on horizons */
    wf = (allocate(
        state,
        ((std::mem::size_of::<WF>() as u64).wrapping_mul((maxnrays as u64)))
            .wrapping_mul((ntwf as u64)) as usize,
    ) as *mut WF);

    initetra = in_which_tetra(xmin, ymin, xmax, ymax, zmax, source, point, tetra, ntetra);
    /* number of tetra */
    if (initetra < 0i32) {
        err!("Source out of model");
    }

    fprintf!(&*jpfp, "source in tetra %d\n", initetra);

    initreg = (*(tetra).offset((initetra) as isize)).ireg;

    fprintf!(&*jpfp, "initreg=%d\n", initreg);
    /*********************************************
    Initialize the ray position
    to be the source position.
    ********************************************* */
    {
        iray = 0i32;
        loop {
            if !(iray < maxnrays) {
                break;
            }
            {
                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).x).cast::<f32>())
                    .offset((0i32) as isize)) = (*(source).offset((0i32) as isize));

                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).x).cast::<f32>())
                    .offset((1i32) as isize)) = (*(source).offset((1i32) as isize));

                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).x).cast::<f32>())
                    .offset((2i32) as isize)) = (*(source).offset((2i32) as isize));

                (*(ray).offset((iray) as isize)).xseg =
                    array2(state, (3i32 as usize), (maxnsegs as usize));
                // A reference sweep of limits 1..=22 found invalid overflow
                // row pointers for the small allocations (1..=15), but a
                // usable neighbouring pointer from 16 onward. Preserve the
                // observed small-allocation signals without invalidating
                // row zero or every larger segment table's extra row.
                if maxnsegs > 0 && maxnsegs < 16 {
                    *(*(ray).offset(iray as isize)).xseg.add(maxnsegs as usize) =
                        std::ptr::null_mut();
                }

                (*(ray).offset((iray) as isize)).tseg = allocate(
                    state,
                    ((maxnsegs as usize) as usize + 1) * std::mem::size_of::<f32>(),
                )
                .cast::<f32>();
                if (*(ray).offset(iray as isize)).tseg.is_null() {
                    err!("ealloc.c: malloc failed (%s)", oxwave_core::err::strerror(12));
                }

                {
                    iseg = 0i32;
                    loop {
                        if !(iseg < maxnsegs) {
                            break;
                        }
                        (*((*(ray).offset((iray) as isize)).tseg).offset((iseg) as isize)) =
                            (0f64 as f32);
                        {
                            let old = iseg;
                            iseg = iseg + 1;
                            old
                        };
                    }
                }

                (*(ray).offset((iray) as isize)).icode = NORMAL;

                (*(ray).offset((iray) as isize)).shootupward = shootupward;

                (*(ray).offset((iray) as isize)).itetra = initetra;

                (*(ray).offset((iray) as isize)).lastfacet = (-1i32);

                (*(ray).offset((iray) as isize)).ionce = 0i32;

                (*(ray).offset((iray) as isize)).iref = 0i32;

                (*(ray).offset((iray) as isize)).nseg = 0i32;

                (*(ray).offset((iray) as isize)).r = (0f64 as f32);

                (*(ray).offset((iray) as isize)).ttotal = (0f64 as f32);
            }
            {
                let old = iray;
                iray = iray + 1;
                old
            };
        }
    }

    {
        iray = 0i32;
        loop {
            if !(iray < ((nazimuth).wrapping_mul(ntakeoff)).wrapping_add(1i32)) {
                break;
            }
            {
                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).x).cast::<f32>())
                    .offset((0i32) as isize)) = (*(source).offset((0i32) as isize));

                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).x).cast::<f32>())
                    .offset((1i32) as isize)) = (*(source).offset((1i32) as isize));

                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).x).cast::<f32>())
                    .offset((2i32) as isize)) = (*(source).offset((2i32) as isize));

                (*(*((*(ray).offset((iray) as isize)).xseg).offset((0i32) as isize))
                    .offset((0i32) as isize)) = (*(source).offset((0i32) as isize));

                (*(*((*(ray).offset((iray) as isize)).xseg).offset((0i32) as isize))
                    .offset((1i32) as isize)) = (*(source).offset((1i32) as isize));

                (*(*((*(ray).offset((iray) as isize)).xseg).offset((0i32) as isize))
                    .offset((2i32) as isize)) = (*(source).offset((2i32) as isize));
            }
            {
                let old = iray;
                iray = iray + 1;
                old
            };
        }
    }

    {
        itwf = 1i32;
        loop {
            if !(itwf < (ntwf).wrapping_sub(1i32)) {
                break;
            }
            {
                {
                    iray = 0i32;
                    loop {
                        if !(iray < maxnrays) {
                            break;
                        }
                        {
                            ithis = ((itwf).wrapping_mul(maxnrays)).wrapping_add(iray);

                            (*(std::ptr::addr_of_mut!((*(wf).offset((ithis) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize)) = ((-99999i32) as f32);

                            (*(std::ptr::addr_of_mut!((*(wf).offset((ithis) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize)) = ((-99999i32) as f32);

                            (*(std::ptr::addr_of_mut!((*(wf).offset((ithis) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize)) = ((-99999i32) as f32);
                        }
                        {
                            let old = iray;
                            iray = iray + 1;
                            old
                        };
                    }
                }
            }
            {
                let old = itwf;
                itwf = itwf + 1;
                old
            };
        }
    }

    azimuth = ((2f64 * 3.1415926000000001f64) as f32);

    dazimuth = (azimuth / (nazimuth as f32));

    dtakeoff = (takeoff / (ntakeoff as f32));
    /*********************************************************
    The maximum number of rays must not be more than maxnrays
    *********************************************************/
    if ((ntakeoff).wrapping_mul(nazimuth) >= maxnrays) {
        err!(
            "too many rays: (ntakeoff*nazimuth=)%d>=(maxnrays=%d",
            (ntakeoff).wrapping_mul(nazimuth),
            maxnrays
        );
    }

    done = NOT_YET;

    s = s_intp(
        point,
        std::ptr::addr_of_mut!((*(tetra).offset((initetra) as isize))),
        source,
    );
    /* the current position */
    sqrts = ((s as f64).sqrt() as f32);

    {
        itwf = 0i32;
        loop {
            if !(itwf < ntwf) {
                break;
            }
            {
                ithis = (itwf).wrapping_mul(maxnrays);

                inext = ((itwf).wrapping_add(1i32)).wrapping_mul(maxnrays);

                if ((done as u32) == (ALREADY as u32)) {
                    break;
                }

                if (itwf == 0i32) {
                    /*****************************************
                    Pointing to the first triangle
                    **************************************** */
                    headp = std::ptr::addr_of_mut!(head);
                    /********************************************
                    There are nazimuth*ntakeoff+1 rays starting for each shot
                    ********************************************/
                    {
                        iray = 0i32;
                        loop {
                            if !(iray <= (nazimuth).wrapping_mul(ntakeoff)) {
                                break;
                            }
                            {
                                /************************************
                                   Ray 0 is the center; takeoff is the
                                take-off angle; azimuth is the azimuth
                                ************************************/
                                if (iray == 0i32) {
                                    itakeoff = 0i32;

                                    sintakeoff = (0f64 as f32);

                                    costakeoff = (1f64 as f32);
                                } else {
                                    itakeoff = (1i32).wrapping_add(
                                        ((iray).wrapping_sub(1i32))
                                            .checked_div(nazimuth)
                                            .unwrap_or(0),
                                    );

                                    temp = (dtakeoff * (itakeoff as f32));

                                    sintakeoff = ((temp as f64).sin() as f32);

                                    costakeoff = ((temp as f64).cos() as f32);
                                }
                                /***************************************
                                Index for azimuth
                                ****************************************/
                                if (iray == 0i32) {
                                    iazimuth = 0i32;

                                    sinazimuth = (0f64 as f32);

                                    cosazimuth = (1f64 as f32);
                                } else {
                                    iazimuth = ((iray).wrapping_sub(1i32))
                                        .checked_rem(nazimuth)
                                        .unwrap_or(0);

                                    temp = (((((iazimuth as f64)
                                        + (0.5f64 * (((itakeoff).wrapping_sub(1i32)) as f64)))
                                        * (dazimuth as f64))
                                        + 0.10000000000000001f64)
                                        as f32);

                                    sinazimuth = ((temp as f64).sin() as f32);

                                    cosazimuth = ((temp as f64).cos() as f32);
                                }
                                /***************************************
                                the equation of the slowness vector is
                                px^2 + py^2 + pz^2 = s,
                                there are two freedoms.
                                ***************************************/
                                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).p)
                                    .cast::<f32>())
                                .offset((0i32) as isize)) = ((sintakeoff * cosazimuth) * sqrts);

                                (*(std::ptr::addr_of_mut!((*(wf).offset((iray) as isize)).p)
                                    .cast::<f32>())
                                .offset((1i32) as isize)) = ((sintakeoff * sinazimuth) * sqrts);
                                /* If the rays are shot upward, put -1 here;
                                otherwise, put +1 here */
                                if (shootupward == 1i32) {
                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset((iray) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize)) = ((-costakeoff) * sqrts);
                                } else {
                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset((iray) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize)) = (costakeoff * sqrts);
                                }
                            }
                            {
                                let old = iray;
                                iray = iray + 1;
                                old
                            };
                        }
                    }
                    /******************************************
                    At first step, we have nazimuth*ntakeoff+1 rays,
                    which form (2*ntakeoff-1)*nazimuth triangles.

                            4   3   2
                           \ 3|2 /
                            \ | /
                           4 \|/ 1
                            5-----0-----1
                           5 /|\ nazimuth-1
                            / | \
                           / 6|  \
                          6    nazimuth
                    ******************************************/
                    (head).i1 = 0i32;

                    (head).i2 = 1i32;

                    (head).i3 = 2i32;

                    trip = headp;

                    {
                        iray = 2i32;
                        loop {
                            if !(iray <= nazimuth) {
                                break;
                            }
                            {
                                (*trip).next =
                                    (allocate(state, (std::mem::size_of::<TRI>() as u64) as usize)
                                        as *mut TRI);

                                trip = (*trip).next;

                                (*trip).i1 = 0i32;

                                (*trip).i2 = iray;

                                (*trip).i3 =
                                    ((iray).checked_rem(nazimuth).unwrap_or(0)).wrapping_add(1i32);
                            }
                            {
                                let old = iray;
                                iray = iray + 1;
                                old
                            };
                        }
                    }

                    {
                        iazimuth = 0i32;
                        loop {
                            if !(iazimuth < (ntakeoff).wrapping_sub(1i32)) {
                                break;
                            }
                            {
                                {
                                    itakeoff = 0i32;
                                    loop {
                                        if !(itakeoff < nazimuth) {
                                            break;
                                        }
                                        {
                                            (*trip).next = (allocate(
                                                state,
                                                (std::mem::size_of::<TRI>() as u64) as usize,
                                            )
                                                as *mut TRI);

                                            trip = (*trip).next;

                                            (*trip).i1 = ((itakeoff).wrapping_add(1i32))
                                                .wrapping_add((iazimuth).wrapping_mul(nazimuth));

                                            (*trip).i2 = ((((itakeoff).wrapping_add(1i32))
                                                .checked_rem(nazimuth)
                                                .unwrap_or(0))
                                            .wrapping_add(1i32))
                                            .wrapping_add((iazimuth).wrapping_mul(nazimuth));

                                            (*trip).i3 = ((*trip).i1).wrapping_add(nazimuth);

                                            (*trip).next = (allocate(
                                                state,
                                                (std::mem::size_of::<TRI>() as u64) as usize,
                                            )
                                                as *mut TRI);

                                            trip = (*trip).next;

                                            (*trip).i1 = ((((itakeoff).wrapping_add(1i32))
                                                .checked_rem(nazimuth)
                                                .unwrap_or(0))
                                            .wrapping_add(1i32))
                                            .wrapping_add((iazimuth).wrapping_mul(nazimuth));

                                            (*trip).i2 = ((itakeoff).wrapping_add(1i32))
                                                .wrapping_add(
                                                    ((iazimuth).wrapping_add(1i32))
                                                        .wrapping_mul(nazimuth),
                                                );

                                            (*trip).i3 = ((*trip).i1).wrapping_add(nazimuth);
                                        }
                                        {
                                            let old = itakeoff;
                                            itakeoff = itakeoff + 1;
                                            old
                                        };
                                    }
                                }
                            }
                            {
                                let old = iazimuth;
                                iazimuth = iazimuth + 1;
                                old
                            };
                        }
                    }

                    (*trip).next = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI);

                    nrays = ((nazimuth).wrapping_mul(ntakeoff)).wrapping_add(1i32);

                    ntris =
                        (((2i32).wrapping_mul(ntakeoff)).wrapping_sub(1i32)).wrapping_mul(nazimuth);

                    fprintf!(&*jpfp, "Starting: ntrips=%d nrays=%d\n", ntris, nrays);
                } else {
                    /****************************************************
                    Test if the head is spoiled. If the first ray is spoiled,
                            the next ray will take the place; All triangles connected
                            with these rays need to be removed;
                    ****************************************************/
                    while ((((((*(ray).offset(((*headp).i1) as isize)).icode as u32)
                        == (SPOILED as u32))
                        || (((*(ray).offset(((*headp).i2) as isize)).icode as u32)
                            == (SPOILED as u32)))
                        || (((*(ray).offset(((*headp).i3) as isize)).icode as u32)
                            == (SPOILED as u32)))
                        && (ntris > 1i32))
                    {
                        headp = (*headp).next;

                        ntris = ntris.wrapping_sub(1);

                        fprintf!(&*jpfp, "head spoiled, cut off: ntris=%d\n", ntris);
                    }

                    if (ntris < 2i32) {
                        break;
                    }
                    /****************************************************
                      Need to test if all triangles are still good. Loop over
                      all triangle, to see is all three rays of which the
                      triangle is composed to be NORMAL.
                    ****************************************************/
                    {
                        trip = headp;
                        loop {
                            {
                                /* make sure that next triangle is not the end*/
                                if ((*trip).next
                                    == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI))
                                {
                                    break;
                                }
                                /* if next triangle is good, then continue to
                                test the next */
                                if (((((*(ray).offset(((*(*trip).next).i1) as isize)).icode
                                    as u32)
                                    != (SPOILED as u32))
                                    && (((*(ray).offset(((*(*trip).next).i2) as isize)).icode
                                        as u32)
                                        != (SPOILED as u32)))
                                    && (((*(ray).offset(((*(*trip).next).i3) as isize)).icode
                                        as u32)
                                        != (SPOILED as u32)))
                                {
                                    trip = (*trip).next;

                                    continue;
                                }
                                /************************************************
                                if the program execuates here, it means that
                                triangle trip is good; but trip->next is bad.
                                We need to cut off this triangle. What could
                                happen to the one or more rays of this this triangle
                                for this misfortune to happen?
                                ************************************************/
                                if (ntris < 1i32) {
                                    break;
                                }

                                ntris = ntris.wrapping_sub(1);
                                /***********************************************
                                trip is good; but trip->next is bad. We need to see
                                if trip0->next->next is also bad, or even further.
                                In these case, we need to cut off all the consecutive
                                bad triangles.
                                ************************************************/
                                {
                                    itrip = (*trip).next;
                                    loop {
                                        {
                                            /* bad to the very end */
                                            if ((*itrip).next
                                                == ((std::ptr::null_mut::<std::ffi::c_void>())
                                                    as *mut TRI))
                                            {
                                                break;
                                            }
                                            /* good. the bad chain has stopped */
                                            if (((((*(ray).offset(((*(*itrip).next).i1) as isize))
                                                .icode
                                                as u32)
                                                != (SPOILED as u32))
                                                && (((*(ray).offset(((*(*itrip).next).i2) as isize))
                                                    .icode
                                                    as u32)
                                                    != (SPOILED as u32)))
                                                && (((*(ray).offset(((*(*itrip).next).i3) as isize))
                                                    .icode
                                                    as u32)
                                                    != (SPOILED as u32)))
                                            {
                                                break;
                                            }

                                            ntris = ntris.wrapping_sub(1);

                                            retire_triangle((itrip as *mut std::ffi::c_void));
                                        }
                                        {
                                            itrip = (*itrip).next;
                                            itrip
                                        };
                                    }
                                }
                                /* now: both trip and trip->next are good ones */
                                (*trip).next = (*itrip).next;
                            }
                        }
                    }
                    /* from headp to the end, there are ntris triangles
                    they are all good ones. */
                    /* to dump the wavefront */
                    if ((wffp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
                        && (itwf == ifwf2dump))
                    {
                        fprintf!(&*wffp, "%d = nwf2dump\n", nwf2dump);
                    }
                    /****************************************************
                    Test if we need to add some triangles
                    ****************************************************/
                    {
                        trip = headp;
                        loop {
                            if !(trip != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI)) {
                                break;
                            }
                            {
                                if (nrays >= maxnrays) {
                                    break;
                                }

                                if (ntris >= maxntris) {
                                    fprintf!(
                                        &*jpfp,
                                        "ntris >= maxntris: wavefront after this step is spoiled\n"
                                    );

                                    break;
                                }

                                if (((((*(ray).offset(((*trip).i1) as isize)).icode as u32)
                                    == (SPOILED as u32))
                                    || (((*(ray).offset(((*trip).i2) as isize)).icode as u32)
                                        == (SPOILED as u32)))
                                    || (((*(ray).offset(((*trip).i3) as isize)).icode as u32)
                                        == (SPOILED as u32)))
                                {
                                    err!("Triangle wrong: Impossible! But that's life!");
                                }
                                /* ray one */
                                x0 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i1)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((0i32) as isize));

                                y0 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i1)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((1i32) as isize));

                                z0 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i1)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((2i32) as isize));
                                /* ray two */
                                x1 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i2)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((0i32) as isize));

                                y1 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i2)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((1i32) as isize));

                                z1 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i2)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((2i32) as isize));
                                /* ray three */
                                x2 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i3)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((0i32) as isize));

                                y2 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i3)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((1i32) as isize));

                                z2 = (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((ithis).wrapping_add((*trip).i3)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((2i32) as isize));
                                /************************************************
                                    Write the current wavefront to the wffile. This
                                can be visualized by viewer3
                                ************************************************/
                                if ((check_tris as u32) == (PLEASE as u32)) {
                                    (*(std::ptr::addr_of_mut!(dist).cast::<f32>())
                                        .offset((0i32) as isize)) = ((((x0 - x1) * (x0 - x1))
                                        + ((y0 - y1) * (y0 - y1)))
                                        + ((z0 - z1) * (z0 - z1)));

                                    (*(std::ptr::addr_of_mut!(dist).cast::<f32>())
                                        .offset((1i32) as isize)) = ((((x2 - x1) * (x2 - x1))
                                        + ((y2 - y1) * (y2 - y1)))
                                        + ((z2 - z1) * (z2 - z1)));

                                    (*(std::ptr::addr_of_mut!(dist).cast::<f32>())
                                        .offset((2i32) as isize)) = ((((x2 - x0) * (x2 - x0))
                                        + ((y2 - y0) * (y2 - y0)))
                                        + ((z2 - z0) * (z2 - z0)));

                                    distmax = (*(std::ptr::addr_of_mut!(dist).cast::<f32>())
                                        .offset((0i32) as isize));

                                    iedgemax = 0i32;

                                    {
                                        iedge = 1i32;
                                        loop {
                                            if !(iedge < 3i32) {
                                                break;
                                            }
                                            {
                                                if ((*(std::ptr::addr_of_mut!(dist).cast::<f32>())
                                                    .offset((iedge) as isize))
                                                    > distmax)
                                                {
                                                    iedgemax = iedge;

                                                    distmax = (*(std::ptr::addr_of_mut!(dist)
                                                        .cast::<f32>())
                                                    .offset((iedge) as isize));
                                                }
                                            }
                                            {
                                                let old = iedge;
                                                iedge = iedge + 1;
                                                old
                                            };
                                        }
                                    }
                                    /* this triangle is ok. all edges within limits,
                                    nothing needs to be done with this triangle */
                                    if (distmax < (edgemax * edgemax)) {
                                        {
                                            trip = (*trip).next;
                                            trip
                                        };
                                        continue;
                                    }
                                    /* this triangle is too big. needs to be splitted
                                    into two triangles */
                                    if (iedgemax == 0i32) {
                                        /* vertex 1 & 2 */
                                        xvertex0 = x0;

                                        xvertex1 = x1;

                                        yvertex0 = y0;

                                        yvertex1 = y1;

                                        zvertex0 = z0;

                                        zvertex1 = z1;

                                        tripi0 = (*trip).i1;

                                        tripi1 = (*trip).i2;
                                    } else {
                                        if (iedgemax == 1i32) {
                                            /* vertex 2 & 3 */
                                            xvertex0 = x1;

                                            xvertex1 = x2;

                                            yvertex0 = y1;

                                            yvertex1 = y2;

                                            zvertex0 = z1;

                                            zvertex1 = z2;

                                            tripi0 = (*trip).i2;

                                            tripi1 = (*trip).i3;
                                        } else {
                                            /* vertex 1 & 3 */
                                            xvertex0 = x0;

                                            xvertex1 = x2;

                                            yvertex0 = y0;

                                            yvertex1 = y2;

                                            zvertex0 = z0;

                                            zvertex1 = z2;

                                            tripi0 = (*trip).i1;

                                            tripi1 = (*trip).i3;
                                        }
                                    }

                                    if ((*(ray).offset((tripi0) as isize)).shootupward
                                        != (*(ray).offset((tripi1) as isize)).shootupward)
                                    {
                                        {
                                            trip = (*trip).next;
                                            trip
                                        };
                                        continue;
                                    }
                                    /* better not add one here */
                                    /***********************************************
                                    Add one more ray: nrays: Interpolation of ray
                                                data, whatever you need, could be traveltime;
                                                amplitude; phase, and even more.
                                          ***********************************************/
                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize)) =
                                        ((((xvertex0 + xvertex1) as f64) / 2f64) as f32);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((1i32) as isize)) =
                                        ((((yvertex0 + yvertex1) as f64) / 2f64) as f32);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize)) =
                                        ((((zvertex0 + zvertex1) as f64) / 2f64) as f32);
                                    /* Needs to determine which tetra it is in */
                                    fprintf!(&*jpfp, "Add one more ray\n");

                                    (*(ray).offset((nrays) as isize)).itetra = in_which_tetra(
                                        xmin,
                                        ymin,
                                        xmax,
                                        ymax,
                                        zmax,
                                        std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(nrays)) as isize))
                                                .x
                                        )
                                        .cast::<f32>(),
                                        point,
                                        tetra,
                                        ntetra,
                                    );

                                    if ((*(ray).offset((nrays) as isize)).itetra < 0i32) {
                                        (*(ray).offset((nrays) as isize)).icode = SPOILED;

                                        fprintf!(&*jpfp, "attempt to add this ray but failed\n");
                                    }

                                    (*(ray).offset((nrays) as isize)).shootupward =
                                        (*(ray).offset((tripi0) as isize)).shootupward;

                                    (*(ray).offset((nrays) as isize)).r =
                                        (((((*(ray).offset((tripi0) as isize)).r
                                            + (*(ray).offset((tripi1) as isize)).r)
                                            as f64)
                                            / 2f64)
                                            as f32);

                                    (*((*(ray).offset((nrays) as isize)).tseg)
                                        .offset((0i32) as isize)) = ((itwf as f32) * dtwf);

                                    (*(ray).offset((nrays) as isize)).ttotal =
                                        ((itwf as f32) * dtwf);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize)) = (((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(tripi0)) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize))
                                        + (*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(tripi1)) as isize))
                                                .p
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize)))
                                        as f64)
                                        / 2f64)
                                        as f32);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((1i32) as isize)) = (((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(tripi0)) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((1i32) as isize))
                                        + (*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(tripi1)) as isize))
                                                .p
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize)))
                                        as f64)
                                        / 2f64)
                                        as f32);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize)) = (((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(tripi0)) as isize)).p
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize))
                                        + (*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(tripi1)) as isize))
                                                .p
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize)))
                                        as f64)
                                        / 2f64)
                                        as f32);

                                    s = s_intp(
                                        point,
                                        std::ptr::addr_of_mut!(
                                            (*(tetra).offset((initetra) as isize))
                                        ),
                                        std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(nrays)) as isize))
                                                .x
                                        )
                                        .cast::<f32>(),
                                    );
                                    /* the current position */
                                    s_normalize(
                                        std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(nrays)) as isize))
                                                .p
                                        )
                                        .cast::<f32>(),
                                        s,
                                    );

                                    r = (((((*(ray).offset((tripi0) as isize)).r
                                        + (*(ray).offset((tripi1) as isize)).r)
                                        as f64)
                                        / 2f64) as f32);

                                    d2 = ((((xvertex0
                                        - (*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(nrays)) as isize))
                                                .x
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize)))
                                        * (xvertex0
                                            - (*(std::ptr::addr_of_mut!(
                                                (*(wf).offset(
                                                    ((ithis).wrapping_add(nrays)) as isize
                                                ))
                                                .x
                                            )
                                            .cast::<f32>())
                                            .offset((0i32) as isize))))
                                        + ((yvertex0
                                            - (*(std::ptr::addr_of_mut!(
                                                (*(wf).offset(
                                                    ((ithis).wrapping_add(nrays)) as isize
                                                ))
                                                .x
                                            )
                                            .cast::<f32>())
                                            .offset((1i32) as isize)))
                                            * (yvertex0
                                                - (*(std::ptr::addr_of_mut!(
                                                    (*(wf).offset(
                                                        ((ithis).wrapping_add(nrays)) as isize
                                                    ))
                                                    .x
                                                )
                                                .cast::<f32>())
                                                .offset((1i32) as isize)))))
                                        + ((zvertex0
                                            - (*(std::ptr::addr_of_mut!(
                                                (*(wf).offset(
                                                    ((ithis).wrapping_add(nrays)) as isize
                                                ))
                                                .x
                                            )
                                            .cast::<f32>())
                                            .offset((2i32) as isize)))
                                            * (zvertex0
                                                - (*(std::ptr::addr_of_mut!(
                                                    (*(wf).offset(
                                                        ((ithis).wrapping_add(nrays)) as isize
                                                    ))
                                                    .x
                                                )
                                                .cast::<f32>())
                                                .offset((2i32) as isize)))));

                                    d = ((d2 as f64).sqrt() as f32);

                                    sintheta = (((d as f64)
                                        / (if ((0.20000000000000001f64) > ((r) as f64)) {
                                            (0.20000000000000001f64)
                                        } else {
                                            ((r) as f64)
                                        })) as f32);

                                    sintheta = ((if ((1f64) < ((sintheta) as f64)) {
                                        (1f64)
                                    } else {
                                        ((sintheta) as f64)
                                    }) as f32);

                                    costheta =
                                        ((1f64 - ((sintheta * sintheta) as f64)).sqrt() as f32);

                                    adjust = (r - (r * costheta));

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize)) = ((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize))
                                        as f64)
                                        + (((adjust
                                            * (*(std::ptr::addr_of_mut!(
                                                (*(wf).offset(
                                                    ((ithis).wrapping_add(nrays)) as isize
                                                ))
                                                .p
                                            )
                                            .cast::<f32>())
                                            .offset((0i32) as isize)))
                                            as f64)
                                            / (s as f64).sqrt()))
                                        as f32);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((1i32) as isize)) = ((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((1i32) as isize))
                                        as f64)
                                        + (((adjust
                                            * (*(std::ptr::addr_of_mut!(
                                                (*(wf).offset(
                                                    ((ithis).wrapping_add(nrays)) as isize
                                                ))
                                                .p
                                            )
                                            .cast::<f32>())
                                            .offset((1i32) as isize)))
                                            as f64)
                                            / (s as f64).sqrt()))
                                        as f32);

                                    (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize)) = ((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize))
                                        as f64)
                                        + (((adjust
                                            * (*(std::ptr::addr_of_mut!(
                                                (*(wf).offset(
                                                    ((ithis).wrapping_add(nrays)) as isize
                                                ))
                                                .p
                                            )
                                            .cast::<f32>())
                                            .offset((2i32) as isize)))
                                            as f64)
                                            / (s as f64).sqrt()))
                                        as f32);

                                    (*(*((*(ray).offset((nrays) as isize)).xseg)
                                        .offset((0i32) as isize))
                                    .offset((0i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize));

                                    (*(*((*(ray).offset((nrays) as isize)).xseg)
                                        .offset((0i32) as isize))
                                    .offset((1i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((1i32) as isize));

                                    (*(*((*(ray).offset((nrays) as isize)).xseg)
                                        .offset((0i32) as isize))
                                    .offset((2i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add(nrays)) as isize)).x
                                    )
                                    .cast::<f32>())
                                    .offset((2i32) as isize));
                                    /* Needs to determine which tetra it is in */
                                    (*(ray).offset((nrays) as isize)).itetra = in_which_tetra(
                                        xmin,
                                        ymin,
                                        xmax,
                                        ymax,
                                        zmax,
                                        std::ptr::addr_of_mut!(
                                            (*(wf).offset(((ithis).wrapping_add(nrays)) as isize))
                                                .x
                                        )
                                        .cast::<f32>(),
                                        point,
                                        tetra,
                                        ntetra,
                                    );

                                    if ((*(ray).offset((nrays) as isize)).itetra < 0i32) {
                                        (*(ray).offset((nrays) as isize)).icode = SPOILED;

                                        fprintf!(&*jpfp, "tried to adjust this ray, failed\n");
                                    }
                                    /********************************************
                                    add 2 triangles.  Since any one edge is shared
                                                by 2 triangles, we need to add one triangle on
                                                one side, and the other one on the other side.
                                                i.e. add one triangle to trip, and one to
                                                trip0
                                                ********************************************/
                                    i1 = (*trip).i1;

                                    i2 = (*trip).i2;

                                    i3 = (*trip).i3;

                                    ii1 = {
                                        ii2 = {
                                            ii3 = 0i32;
                                            ii3
                                        };
                                        ii2
                                    };

                                    iedgemax0 = 999i32;
                                    /* find the next neighboring triangle */
                                    {
                                        trip0 = headp;
                                        loop {
                                            if !(trip0
                                                != ((std::ptr::null_mut::<std::ffi::c_void>())
                                                    as *mut TRI))
                                            {
                                                break;
                                            }
                                            {
                                                if (trip0 == trip) {
                                                    {
                                                        trip0 = (*trip0).next;
                                                        trip0
                                                    };
                                                    continue;
                                                }

                                                ii1 = (*trip0).i1;

                                                ii2 = (*trip0).i2;

                                                ii3 = (*trip0).i3;

                                                if ((((ii1 == i1) || (ii1 == i2)) || (ii1 == i3))
                                                    && (((ii2 == i1) || (ii2 == i2))
                                                        || (ii2 == i3)))
                                                {
                                                    iedgemax0 = 0i32;

                                                    break;
                                                } else {
                                                    if ((((ii2 == i1) || (ii2 == i2))
                                                        || (ii2 == i3))
                                                        && (((ii3 == i1) || (ii3 == i2))
                                                            || (ii3 == i3)))
                                                    {
                                                        iedgemax0 = 1i32;

                                                        break;
                                                    } else {
                                                        if ((((ii1 == i1) || (ii1 == i2))
                                                            || (ii1 == i3))
                                                            && (((ii3 == i1) || (ii3 == i2))
                                                                || (ii3 == i3)))
                                                        {
                                                            iedgemax0 = 2i32;

                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                            {
                                                trip0 = (*trip0).next;
                                                trip0
                                            };
                                        }
                                    }

                                    temtrip = (*trip).next;

                                    if (iedgemax == 0i32) {
                                        /* vertex 1 & 2 */
                                        (*trip).i2 = nrays;

                                        (*trip).next = (allocate(
                                            state,
                                            (std::mem::size_of::<TRI>() as u64) as usize,
                                        )
                                            as *mut TRI);

                                        trip = (*trip).next;

                                        (*trip).i1 = nrays;

                                        (*trip).i2 = i2;

                                        (*trip).i3 = i3;
                                    } else {
                                        if (iedgemax == 1i32) {
                                            /* vertex 2 & 3 */
                                            (*trip).i3 = nrays;

                                            (*trip).next = (allocate(
                                                state,
                                                (std::mem::size_of::<TRI>() as u64) as usize,
                                            )
                                                as *mut TRI);

                                            trip = (*trip).next;

                                            (*trip).i1 = i1;

                                            (*trip).i2 = nrays;

                                            (*trip).i3 = i3;
                                        } else {
                                            /* vertex 1 & 3 */
                                            (*trip).i3 = nrays;

                                            (*trip).next = (allocate(
                                                state,
                                                (std::mem::size_of::<TRI>() as u64) as usize,
                                            )
                                                as *mut TRI);

                                            trip = (*trip).next;

                                            (*trip).i1 = nrays;

                                            (*trip).i2 = i2;

                                            (*trip).i3 = i3;
                                        }
                                    }

                                    (*trip).next = temtrip;

                                    ntris = ntris.wrapping_add(1);

                                    if (((((trip0
                                        == ((std::ptr::null_mut::<std::ffi::c_void>())
                                            as *mut TRI))
                                        || (iedgemax0 == 999i32))
                                        || (ii1 == ii2))
                                        || (ii1 == ii3))
                                        || (ii2 == ii3))
                                    {
                                        nrays = nrays.wrapping_add(1);

                                        {
                                            trip = (*trip).next;
                                            trip
                                        };
                                        continue;
                                    }

                                    temtrip = (*trip0).next;

                                    if (iedgemax0 == 0i32) {
                                        /* vertex 1 & 2 */
                                        (*trip0).i2 = nrays;

                                        (*trip0).next = (allocate(
                                            state,
                                            (std::mem::size_of::<TRI>() as u64) as usize,
                                        )
                                            as *mut TRI);

                                        trip0 = (*trip0).next;

                                        (*trip0).i1 = nrays;

                                        (*trip0).i2 = ii2;

                                        (*trip0).i3 = ii3;
                                    } else {
                                        if (iedgemax0 == 1i32) {
                                            /* vertex 2 & 3 */
                                            (*trip0).i3 = nrays;

                                            (*trip0).next = (allocate(
                                                state,
                                                (std::mem::size_of::<TRI>() as u64) as usize,
                                            )
                                                as *mut TRI);

                                            trip0 = (*trip0).next;

                                            (*trip0).i1 = ii1;

                                            (*trip0).i2 = nrays;

                                            (*trip0).i3 = ii3;
                                        } else {
                                            /* vertex 1 & 3 */
                                            (*trip0).i3 = nrays;

                                            (*trip0).next = (allocate(
                                                state,
                                                (std::mem::size_of::<TRI>() as u64) as usize,
                                            )
                                                as *mut TRI);

                                            trip0 = (*trip0).next;

                                            (*trip0).i1 = nrays;

                                            (*trip0).i2 = ii2;

                                            (*trip0).i3 = ii3;
                                        }
                                    }

                                    (*trip0).next = temtrip;

                                    ntris = ntris.wrapping_add(1);

                                    nrays = nrays.wrapping_add(1);
                                }
                            }
                            {
                                trip = (*trip).next;
                                trip
                            };
                        }
                    }
                    /* end if check_tris */
                    /************************************************
                        Write the current wavefront to the wffile. This
                          can be visualized by viewer3
                    ************************************************/
                    if (((wffp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile))
                        && (itwf >= ifwf2dump))
                        && (itwf < (ifwf2dump).wrapping_add(nwf2dump)))
                    {
                        {
                            {
                                {
                                    trip = headp;
                                    trip
                                };
                                {
                                    ntr = 0i32;
                                    ntr
                                }
                            };
                            loop {
                                if !(trip
                                    != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI))
                                {
                                    break;
                                }
                                {
                                    if (((((((((((*(std::ptr::addr_of_mut!(
                                        (*(wf).offset(((ithis).wrapping_add((*trip).i1)) as isize))
                                            .x
                                    )
                                    .cast::<f32>())
                                    .offset((0i32) as isize))
                                        as f64)
                                        > 0f64)
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i1)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i1)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i2)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i2)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i2)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i3)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && (((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i3)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize))
                                            as f64)
                                            > 0f64))
                                        && ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i3)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            != 0.0))
                                    {
                                        ntr = ntr.wrapping_add(1);
                                    } else {
                                        ntr = 0i32;

                                        break;
                                    }
                                }
                                {
                                    trip = (*trip).next;
                                    trip
                                };
                            }
                        }

                        fprintf!(&*wffp, "%d = ntris\n", ntr);

                        if (ntr > 0i32) {
                            {
                                trip = headp;
                                loop {
                                    if !(trip
                                        != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI))
                                    {
                                        break;
                                    }
                                    fprintf!(
                                        &*wffp,
                                        "%f %f %f %f %f %f %f %f %f\n",
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i1)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i1)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i1)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i2)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i2)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i2)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i3)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((0i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i3)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((1i32) as isize))
                                            as f64),
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(wf).offset(
                                                ((ithis).wrapping_add((*trip).i3)) as isize
                                            ))
                                            .x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            as f64)
                                    );
                                    {
                                        trip = (*trip).next;
                                        trip
                                    };
                                }
                            }
                        }
                    }

                    if ((ntris < 2i32) || (nrays < 4i32)) {
                        break;
                    }
                }
                /* end if it!=0 */
                if (itwf >= (ntwf).wrapping_sub(1i32)) {
                    break;
                }

                done = ALREADY;

                {
                    iray = 0i32;
                    loop {
                        if !(iray < nrays) {
                            break;
                        }
                        {
                            if ((((*(ray).offset((iray) as isize)).icode as u32)
                                == (SPOILED as u32))
                                || (((*(ray).offset((iray) as isize)).icode as u32)
                                    == (HIT_SURFACE as u32)))
                            {
                                {
                                    let old = iray;
                                    iray = iray + 1;
                                    old
                                };
                                continue;
                            }

                            done = NOT_YET;

                            (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                (*(wf).offset(((ithis).wrapping_add(iray)) as isize)).x
                            )
                            .cast::<f32>())
                            .offset((0i32) as isize));

                            (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                (*(wf).offset(((ithis).wrapping_add(iray)) as isize)).x
                            )
                            .cast::<f32>())
                            .offset((1i32) as isize));

                            (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                (*(wf).offset(((ithis).wrapping_add(iray)) as isize)).x
                            )
                            .cast::<f32>())
                            .offset((2i32) as isize));

                            (*(ray).offset((iray) as isize)).t = (0f64 as f32);

                            (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).p)
                                .cast::<f32>())
                            .offset((0i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                (*(wf).offset(((ithis).wrapping_add(iray)) as isize)).p
                            )
                            .cast::<f32>())
                            .offset((0i32) as isize));

                            (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).p)
                                .cast::<f32>())
                            .offset((1i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                (*(wf).offset(((ithis).wrapping_add(iray)) as isize)).p
                            )
                            .cast::<f32>())
                            .offset((1i32) as isize));

                            (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).p)
                                .cast::<f32>())
                            .offset((2i32) as isize)) = (*(std::ptr::addr_of_mut!(
                                (*(wf).offset(((ithis).wrapping_add(iray)) as isize)).p
                            )
                            .cast::<f32>())
                            .offset((2i32) as isize));
                            /******************************************
                            Check if the current ray point source is out
                            of bounds
                            ******************************************/
                            while ((HIT_FACET_FIRST as u32)
                                == (({
                                    result = trace1tetra(
                                        maxnsegs,
                                        xmin,
                                        ymin,
                                        xmax,
                                        ymax,
                                        zmax,
                                        (ray).offset((iray) as isize),
                                        npoint,
                                        point,
                                        nfacet,
                                        facet,
                                        ntetra,
                                        tetra,
                                        irefseq,
                                        nrefseq,
                                        jpfp,
                                        dtwf,
                                    );
                                    result
                                }) as u32))
                            {}

                            if ((result as u32) == (OUT_OF_MODEL as u32)) {
                                (*(ray).offset((iray) as isize)).icode = SPOILED;
                            }

                            if (itwf < (ntwf).wrapping_sub(1i32)) {
                                (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((inext).wrapping_add(iray)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((0i32) as isize)) =
                                    (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                        .cast::<f32>())
                                    .offset((0i32) as isize));

                                (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((inext).wrapping_add(iray)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((1i32) as isize)) =
                                    (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                        .cast::<f32>())
                                    .offset((1i32) as isize));

                                (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((inext).wrapping_add(iray)) as isize)).x
                                )
                                .cast::<f32>())
                                .offset((2i32) as isize)) =
                                    (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                        .cast::<f32>())
                                    .offset((2i32) as isize));

                                (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((inext).wrapping_add(iray)) as isize)).p
                                )
                                .cast::<f32>())
                                .offset((0i32) as isize)) =
                                    (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).p)
                                        .cast::<f32>())
                                    .offset((0i32) as isize));

                                (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((inext).wrapping_add(iray)) as isize)).p
                                )
                                .cast::<f32>())
                                .offset((1i32) as isize)) =
                                    (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).p)
                                        .cast::<f32>())
                                    .offset((1i32) as isize));

                                (*(std::ptr::addr_of_mut!(
                                    (*(wf).offset(((inext).wrapping_add(iray)) as isize)).p
                                )
                                .cast::<f32>())
                                .offset((2i32) as isize)) =
                                    (*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).p)
                                        .cast::<f32>())
                                    .offset((2i32) as isize));
                            } else {
                                fprintf!(&*jpfp, "ntwf is not enough\n");

                                enough_ntwf = WF_FALSE;
                            }
                            /* plot the rays wherever they go */
                            if ((result as u32) == (HIT_SURFACE_FIRST as u32)) {
                                check_tris = PLEASE_DO_NOT;

                                (*(ray).offset((iray) as isize)).icode = HIT_SURFACE;

                                if (((*(std::ptr::addr_of_mut!((*(ray).offset((iray) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                                    as f64)
                                    .abs()
                                    > 1.0E-4f64)
                                {
                                    fprintf!(
                                        &*jpfp,
                                        "x[2]!=0.0, =%f\n",
                                        ((*(std::ptr::addr_of_mut!(
                                            (*(ray).offset((iray) as isize)).x
                                        )
                                        .cast::<f32>())
                                        .offset((2i32) as isize))
                                            as f64)
                                    );
                                }
                            }
                        }
                        {
                            let old = iray;
                            iray = iray + 1;
                            old
                        };
                    }
                }
            }
            {
                let old = itwf;
                itwf = itwf + 1;
                old
            };
        }
    }
    /* next it */
    /******************************************************************
    Now, let's output the ray information for viewer3 to check if all
    the rays are correct.
    ******************************************************************/
    if (rayfp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile)) {
        fprintf!(&*rayfp, "%d =Maximum number of segments\n", maxnsegs);

        ngoodrays = 0i32;

        {
            iray = 0i32;
            loop {
                if !(iray < nrays) {
                    break;
                }
                if (((*(ray).offset((iray) as isize)).icode as u32) != (SPOILED as u32)) {
                    ngoodrays = ngoodrays.wrapping_add(1);
                }
                {
                    let old = iray;
                    iray = iray + 1;
                    old
                };
            }
        }

        fprintf!(&*rayfp, "%d =Number of rays\n", ngoodrays);

        {
            iray = 0i32;
            loop {
                if !(iray < nrays) {
                    break;
                }
                {
                    if (((*(ray).offset((iray) as isize)).icode as u32) == (SPOILED as u32)) {
                        {
                            let old = iray;
                            iray = iray + 1;
                            old
                        };
                        continue;
                    }

                    fprintf!(
                        &*rayfp,
                        "%d=nseg %f=ttotal\n",
                        ((*(ray).offset((iray) as isize)).nseg).wrapping_add(1i32),
                        ((*(ray).offset((iray) as isize)).ttotal as f64)
                    );

                    tsum = (0f64 as f32);

                    {
                        iseg = 0i32;
                        loop {
                            if !(iseg <= (*(ray).offset((iray) as isize)).nseg) {
                                break;
                            }
                            {
                                tsum = ((tsum as f32)
                                    + (*((*(ray).offset((iray) as isize)).tseg)
                                        .offset((iseg) as isize)));

                                fprintf!(
                                    &*rayfp,
                                    "%f %f %f %f %f\n",
                                    ((*(*((*(ray).offset((iray) as isize)).xseg)
                                        .offset((iseg) as isize))
                                    .offset((0i32) as isize))
                                        as f64),
                                    ((*(*((*(ray).offset((iray) as isize)).xseg)
                                        .offset((iseg) as isize))
                                    .offset((1i32) as isize))
                                        as f64),
                                    ((*(*((*(ray).offset((iray) as isize)).xseg)
                                        .offset((iseg) as isize))
                                    .offset((2i32) as isize))
                                        as f64),
                                    (tsum as f64),
                                    ((*((*(ray).offset((iray) as isize)).tseg)
                                        .offset((iseg) as isize))
                                        as f64)
                                );
                            }
                            {
                                let old = iseg;
                                iseg = iseg + 1;
                                old
                            };
                        }
                    }
                }
                {
                    let old = iray;
                    iray = iray + 1;
                    old
                };
            }
        }
    }

    if ((enough_ntwf as u32) == (WF_FALSE as u32)) {
        return 1i32;
    }
    /* out the surface traveltime, only for visualization */
    if (sttfp != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile)) {
        {
            {
                {
                    trip = headp;
                    trip
                };
                {
                    ntr = 0i32;
                    ntr
                }
            };
            loop {
                if !(trip != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI)) {
                    break;
                }
                {
                    if ((((((((*(ray).offset(((*trip).i1) as isize)).icode as u32)
                        == (SPOILED as u32))
                        || (((*(ray).offset(((*trip).i2) as isize)).icode as u32)
                            == (SPOILED as u32)))
                        || (((*(ray).offset(((*trip).i3) as isize)).icode as u32)
                            == (SPOILED as u32)))
                        || (((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i1) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize)) as f64)
                            .abs()
                            > 0.01f64))
                        || (((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i2) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize)) as f64)
                            .abs()
                            > 0.01f64))
                        || (((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i3) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize)) as f64)
                            .abs()
                            > 0.01f64))
                    {
                        {
                            trip = (*trip).next;
                            trip
                        };
                        continue;
                    }

                    ntr = ntr.wrapping_add(1);
                }
                {
                    trip = (*trip).next;
                    trip
                };
            }
        }

        fprintf!(&*std::ptr::addr_of_mut!((*state).errors), "ntr=%d\n", ntr);

        fprintf!(&*sttfp, "%d = ntris\n", ntr);

        if (ntr > 0i32) {
            {
                trip = headp;
                loop {
                    if !(trip != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI)) {
                        break;
                    }
                    {
                        if ((((((((*(ray).offset(((*trip).i1) as isize)).icode as u32)
                            == (SPOILED as u32))
                            || (((*(ray).offset(((*trip).i2) as isize)).icode as u32)
                                == (SPOILED as u32)))
                            || (((*(ray).offset(((*trip).i3) as isize)).icode as u32)
                                == (SPOILED as u32)))
                            || (((*(std::ptr::addr_of_mut!(
                                (*(ray).offset(((*trip).i1) as isize)).x
                            )
                            .cast::<f32>())
                            .offset((2i32) as isize)) as f64)
                                .abs()
                                > 0.01f64))
                            || (((*(std::ptr::addr_of_mut!(
                                (*(ray).offset(((*trip).i2) as isize)).x
                            )
                            .cast::<f32>())
                            .offset((2i32) as isize)) as f64)
                                .abs()
                                > 0.01f64))
                            || (((*(std::ptr::addr_of_mut!(
                                (*(ray).offset(((*trip).i3) as isize)).x
                            )
                            .cast::<f32>())
                            .offset((2i32) as isize)) as f64)
                                .abs()
                                > 0.01f64))
                        {
                            {
                                trip = (*trip).next;
                                trip
                            };
                            continue;
                        }

                        fprintf!(
                            &*sttfp,
                            "%f %f %f %f %f %f %f %f %f %f %f %f\n",
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i1) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize)) as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i1) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize)) as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i1) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize)) as f64),
                            ((*(ray).offset(((*trip).i1) as isize)).ttotal as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize)) as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize)) as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize)) as f64),
                            ((*(ray).offset(((*trip).i2) as isize)).ttotal as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i3) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize)) as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i3) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize)) as f64),
                            ((*(std::ptr::addr_of_mut!((*(ray).offset(((*trip).i3) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize)) as f64),
                            ((*(ray).offset(((*trip).i3) as isize)).ttotal as f64)
                        );
                    }
                    {
                        trip = (*trip).next;
                        trip
                    };
                }
            }
        }

        fclose(&*sttfp);

        sttfp = ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile);
    }
    /**************************************************
    To do the evaluation on the surface, for each
    triangle, receivers inside this triangle will be
    evaluated. If a receiver is evaluated more than
    once, only the one with shortest traveltime
    is kept.
    *************************************************/
    {
        trip = headp;
        loop {
            if !(trip != ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut TRI)) {
                break;
            }
            {
                if (((((*(ray).offset(((*trip).i1) as isize)).icode as u32) == (SPOILED as u32))
                    || (((*(ray).offset(((*trip).i2) as isize)).icode as u32) == (SPOILED as u32)))
                    || (((*(ray).offset(((*trip).i3) as isize)).icode as u32) == (SPOILED as u32)))
                {
                    {
                        trip = (*trip).next;
                        trip
                    };
                    continue;
                }

                tri_multi_intp(
                    std::ptr::addr_of_mut!((*(ray).offset(((*trip).i1) as isize))),
                    std::ptr::addr_of_mut!((*(ray).offset(((*trip).i2) as isize))),
                    std::ptr::addr_of_mut!((*(ray).offset(((*trip).i3) as isize))),
                    fxgd,
                    fygd,
                    dxgd,
                    dygd,
                    nxgd,
                    nygd,
                    ttable,
                    ctable,
                    rtable,
                );
            }
            {
                trip = (*trip).next;
                trip
            };
        }
    }
    /******************************************************
    Now let us output the traveltime table to stdout, this
    will be used for sukdmig3d. Remember sx and sy are in
    dxgd and dygd already.
    ******************************************************/
    {
        iygd = iytf;
        loop {
            if !(iygd < (iytf).wrapping_add(nyt)) {
                break;
            }
            {
                {
                    ixgd = ixtf;
                    loop {
                        if !(ixgd < (ixtf).wrapping_add(nxt)) {
                            break;
                        }
                        (*((&mut (*state).tro).data.as_mut_ptr())
                            .offset(((ixgd).wrapping_sub(ixtf)) as isize)) = ((if ((0f64)
                            > (((*(*(ttable).offset((iygd) as isize)).offset((ixgd) as isize))
                                as f64)
                                * 1000f64))
                        {
                            (0f64)
                        } else {
                            (((*(*(ttable).offset((iygd) as isize)).offset((ixgd) as isize))
                                as f64)
                                * 1000f64)
                        })
                            as f32);
                        {
                            let old = ixgd;
                            ixgd = ixgd + 1;
                            old
                        };
                    }
                }
                /* source left/front range */
                (*state).tro.hdr.gx = ixtf;

                (*state).tro.hdr.gy = iytf;
                /* source positions */
                (*state).tro.hdr.sx =
                    ((((*(source).offset((0i32) as isize)) as f64) * 1000f64) as i32);

                (*state).tro.hdr.sy =
                    ((((*(source).offset((1i32) as isize)) as f64) * 1000f64) as i32);

                (*state).tro.hdr.unscale = (*(source).offset((2i32) as isize));
                /* this is a float */
                {
                    fputtr(
                        &*std::ptr::addr_of_mut!((*state).output),
                        &mut *std::ptr::addr_of_mut!((*state).tro),
                    );
                    0i32
                };

                if (crfp == ((std::ptr::null_mut::<std::ffi::c_void>()) as *mut CFile)) {
                    {
                        let old = iygd;
                        iygd = iygd + 1;
                        old
                    };
                    continue;
                }

                if ((fwrite(
                    std::slice::from_raw_parts(
                        (std::ptr::addr_of_mut!(
                            (*(*(ctable).offset((iygd) as isize)).offset((ixtf) as isize))
                        ) as *mut std::ffi::c_void)
                            .cast::<u8>(),
                        ((std::mem::size_of::<f32>() as u64) * (nxt as u64)) as usize,
                    ),
                    &*crfp,
                ) / (std::mem::size_of::<f32>() as u64) as usize) as u64
                    != (nxt as u64))
                {
                    err!("Can not write to crfile");
                }

                if ((fwrite(
                    std::slice::from_raw_parts(
                        (std::ptr::addr_of_mut!(
                            (*(*(rtable).offset((iygd) as isize)).offset((ixtf) as isize))
                        ) as *mut std::ffi::c_void)
                            .cast::<u8>(),
                        ((std::mem::size_of::<f32>() as u64) * (nxt as u64)) as usize,
                    ),
                    &*crfp,
                ) / (std::mem::size_of::<f32>() as u64) as usize) as u64
                    != (nxt as u64))
                {
                    err!("Can not write to crfile");
                }
            }
            {
                let old = iygd;
                iygd = iygd + 1;
                old
            };
        }
    }

    {
        fflush(&*std::ptr::addr_of_mut!((*state).output));
        0i32
    };

    fprintf!(&*jpfp, "traveltimes for this shot written.\n");

    return 0i32;
}

// zhorz (sutetraray.c:1687).

unsafe fn zhorz(
    mut xx: f32,
    mut yy: f32,
    mut n10: f32,
    mut n20: f32,
    mut n30: f32,
    mut x00: f32,
    mut y00: f32,
    mut z00: f32,
) -> f32 {
    if ((n30 as f64).abs() < 0.001f64) {
        return (9000i32 as f32);
    }

    return (z00 - ((((xx - x00) * n10) + ((yy - y00) * n20)) / n30));
}

// f (sutetraray.c:1725).

unsafe fn f(mut sigma: f32, mut s: f32, mut p: f32, mut a: f32) -> f32 {
    return ((((s + (p * sigma)) as f64)
        + (((0.25f64 * (a as f64)) * (sigma as f64)) * (sigma as f64))) as f32);
}

// cubicsolver (sutetraray.c:1750).

unsafe fn cubicsolver(mut d: f32, mut a: f32, mut b: f32, mut c: f32) -> f32 {
    let mut sigma: f32 = 0.0;
    let mut fderi: f32 = 0.0;
    let mut fsigma: f32 = 0.0;

    let mut MaxIterations: i32 = 20i32;

    let mut iter: i32 = 0;

    sigma = quadsolver(a, b, c);

    if ((d as f64) == 0f64) {
        return sigma;
    }

    if ((sigma as f64) > (0.5f64 * 1.0E+10f64)) {
        sigma = (1f64 as f32);

        MaxIterations = 30i32;
    }

    {
        iter = 0i32;
        loop {
            if !(iter < MaxIterations) {
                break;
            }
            {
                /* fderi=3dx^2+2ax+b */
                fsigma = ((((((d * sigma) + a) * sigma) + b) * sigma) + c);

                if ((fsigma as f64).abs() < 1.0000000000000001E-5f64) {
                    break;
                }

                fderi = ((((((3f64 * (d as f64)) * (sigma as f64)) + (2f64 * (a as f64)))
                    * (sigma as f64))
                    + (b as f64)) as f32);

                if ((fderi as f64).abs() < 1.0000000000000001E-5f64) {
                    break;
                }

                sigma = ((sigma as f32) - (fsigma / fderi));
            }
            {
                let old = iter;
                iter = iter + 1;
                old
            };
        }
    }

    if ((sigma as f64) < 0f64) {
        return (0f64 as f32);
    }

    return sigma;
}

// quadsolver (sutetraray.c:1796).

unsafe fn quadsolver(mut a: f32, mut b: f32, mut c: f32) -> f32 {
    let mut qdelta: f32 = 0.0;

    let mut sigma: f32 = 0.0;
    let mut sigma1: f32 = 0.0;
    let mut sigma2: f32 = 0.0;

    if ((c as f64) == 0f64) {
        if ((a as f64) == 0f64) {
            return (1.0E+10f64 as f32);
        }

        sigma = ((-b) / a);

        if ((sigma as f64) < 0f64) {
            return (1.0E+10f64 as f32);
        }

        return sigma;
    }

    if ((a as f64) == 0f64) {
        if ((b as f64) == 0f64) {
            return (1.0E+10f64 as f32);
        }

        sigma = ((-c) / b);

        if ((sigma as f64) < 0f64) {
            sigma = (1.0E+10f64 as f32);
        }

        return sigma;
    } else {
        qdelta = ((b * b) - (((4i32 as f32) * a) * c));

        if ((qdelta as f64) < (-1.0000000000000001E-5f64)) {
            return (1.0E+10f64 as f32);
        }

        qdelta = ((if ((0f64) > ((qdelta) as f64)) {
            (0f64)
        } else {
            ((qdelta) as f64)
        }) as f32);

        qdelta = ((qdelta as f64).sqrt() as f32);

        sigma1 = ((((qdelta - b) as f64) / (2f64 * (a as f64))) as f32);

        sigma2 = (((((-qdelta) - b) as f64) / (2f64 * (a as f64))) as f32);

        if ((sigma1 as f64) < (-1.0000000000000001E-5f64)) {
            sigma1 = (1.0E+10f64 as f32);
        }

        if ((sigma2 as f64) < (-1.0000000000000001E-5f64)) {
            sigma2 = (1.0E+10f64 as f32);
        }

        sigma = (if ((sigma1) < (sigma2)) {
            (sigma1)
        } else {
            (sigma2)
        });

        return ((if ((0f64) > ((sigma) as f64)) {
            (0f64)
        } else {
            ((sigma) as f64)
        }) as f32);
    }
}

// in_which_tetra (sutetraray.c:1840).

unsafe fn in_which_tetra(
    mut xmin: f32,
    mut ymin: f32,
    mut xmax: f32,
    mut ymax: f32,
    mut zmax: f32,
    mut x: *mut f32,
    mut point: *mut POINT,
    mut tetra: *mut TETRA,
    mut ntetra: i32,
) -> i32 {
    let mut it: i32 = 0;

    let mut ip0: i32 = 0;
    let mut ip1: i32 = 0;
    let mut ip2: i32 = 0;
    let mut ip3: i32 = 0;

    let mut vv: f32 = 0.0;
    let mut vv0: f32 = 0.0;
    let mut vv1: f32 = 0.0;
    let mut vv2: f32 = 0.0;
    let mut vv3: f32 = 0.0;

    let mut xxmin: [f32; 4] = [0.0; 4];
    let mut xxmax: [f32; 4] = [0.0; 4];

    let mut vdiff: f32 = 0.0;

    let mut itnearest: i32 = 0;

    if (((((((*(x).offset((0i32) as isize)) < xmin) || ((*(x).offset((0i32) as isize)) > xmax))
        || ((*(x).offset((1i32) as isize)) < ymin))
        || ((*(x).offset((1i32) as isize)) > ymax))
        || (((*(x).offset((2i32) as isize)) as f64) < 0f64))
        || ((*(x).offset((2i32) as isize)) > zmax))
    {
        return (-3i32);
    }

    itnearest = (-1i32);

    vdiff = (1000i32 as f32);

    {
        it = 0i32;
        loop {
            if !(it < ntetra) {
                break;
            }
            {
                vv = (*(tetra).offset((it) as isize)).v;

                if ((vv as f64) < 0f64) {
                    {
                        let old = it;
                        it = it + 1;
                        old
                    };
                    continue;
                }

                ip0 = (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                    .cast::<i32>())
                .offset((0i32) as isize));

                ip1 = (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                    .cast::<i32>())
                .offset((1i32) as isize));

                ip2 = (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                    .cast::<i32>())
                .offset((2i32) as isize));

                ip3 = (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                    .cast::<i32>())
                .offset((3i32) as isize));

                (*(std::ptr::addr_of_mut!(xxmin).cast::<f32>()).offset((0i32) as isize)) =
                    (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                        .cast::<f32>())
                    .offset((0i32) as isize))
                        < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((0i32) as isize))
                            < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            })
                        }))
                    {
                        (*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                            .cast::<f32>())
                        .offset((0i32) as isize))
                    } else {
                        (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((0i32) as isize))
                            < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            })
                        })
                    });

                (*(std::ptr::addr_of_mut!(xxmax).cast::<f32>()).offset((0i32) as isize)) =
                    (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                        .cast::<f32>())
                    .offset((0i32) as isize))
                        > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((0i32) as isize))
                            > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            })
                        }))
                    {
                        (*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                            .cast::<f32>())
                        .offset((0i32) as isize))
                    } else {
                        (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((0i32) as isize))
                            > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            })
                        })
                    });

                (*(std::ptr::addr_of_mut!(xxmin).cast::<f32>()).offset((1i32) as isize)) =
                    (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                        .cast::<f32>())
                    .offset((1i32) as isize))
                        < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((1i32) as isize))
                            < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            })
                        }))
                    {
                        (*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                            .cast::<f32>())
                        .offset((1i32) as isize))
                    } else {
                        (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((1i32) as isize))
                            < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            })
                        })
                    });

                (*(std::ptr::addr_of_mut!(xxmax).cast::<f32>()).offset((1i32) as isize)) =
                    (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                        .cast::<f32>())
                    .offset((1i32) as isize))
                        > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((1i32) as isize))
                            > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            })
                        }))
                    {
                        (*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                            .cast::<f32>())
                        .offset((1i32) as isize))
                    } else {
                        (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((1i32) as isize))
                            > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            })
                        })
                    });

                (*(std::ptr::addr_of_mut!(xxmin).cast::<f32>()).offset((2i32) as isize)) =
                    (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                        .cast::<f32>())
                    .offset((2i32) as isize))
                        < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize))
                            < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            })
                        }))
                    {
                        (*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize))
                    } else {
                        (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize))
                            < (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                < (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            })
                        })
                    });

                (*(std::ptr::addr_of_mut!(xxmax).cast::<f32>()).offset((2i32) as isize)) =
                    (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                        .cast::<f32>())
                    .offset((2i32) as isize))
                        > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize))
                            > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            })
                        }))
                    {
                        (*(std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize))
                    } else {
                        (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                            .cast::<f32>())
                        .offset((2i32) as isize))
                            > (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            }))
                        {
                            (*(std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                        } else {
                            (if ((*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize))
                                > (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize)))
                            {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            } else {
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            })
                        })
                    });

                if ((((((((*(x).offset((0i32) as isize)) as f64)
                    < (((*(std::ptr::addr_of_mut!(xxmin).cast::<f32>()).offset((0i32) as isize))
                        as f64)
                        - 0.10000000000000001f64))
                    || (((*(x).offset((0i32) as isize)) as f64)
                        > (((*(std::ptr::addr_of_mut!(xxmax).cast::<f32>()).offset((0i32) as isize))
                            as f64)
                            + 0.10000000000000001f64)))
                    || (((*(x).offset((1i32) as isize)) as f64)
                        < (((*(std::ptr::addr_of_mut!(xxmin).cast::<f32>()).offset((1i32) as isize))
                            as f64)
                            - 0.10000000000000001f64)))
                    || (((*(x).offset((1i32) as isize)) as f64)
                        > (((*(std::ptr::addr_of_mut!(xxmax).cast::<f32>()).offset((1i32) as isize))
                            as f64)
                            + 0.10000000000000001f64)))
                    || (((*(x).offset((2i32) as isize)) as f64)
                        < (((*(std::ptr::addr_of_mut!(xxmin).cast::<f32>()).offset((2i32) as isize))
                            as f64)
                            - 0.10000000000000001f64)))
                    || (((*(x).offset((2i32) as isize)) as f64)
                        > (((*(std::ptr::addr_of_mut!(xxmax).cast::<f32>()).offset((2i32) as isize))
                            as f64)
                            + 0.10000000000000001f64)))
                {
                    {
                        let old = it;
                        it = it + 1;
                        old
                    };
                    continue;
                }

                vv0 = tetra_volume(
                    std::slice::from_raw_parts(x, 3),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x).cast::<f32>(),
                        3,
                    ),
                );

                vv1 = tetra_volume(
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(x, 3),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x).cast::<f32>(),
                        3,
                    ),
                );

                vv2 = tetra_volume(
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(x, 3),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip3) as isize)).x).cast::<f32>(),
                        3,
                    ),
                );

                vv3 = tetra_volume(
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip0) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip1) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(
                        std::ptr::addr_of_mut!((*(point).offset((ip2) as isize)).x).cast::<f32>(),
                        3,
                    ),
                    std::slice::from_raw_parts(x, 3),
                );

                if (((((vv0 + vv1) + vv2) + vv3) - vv) < vdiff) {
                    vdiff = ((((vv0 + vv1) + vv2) + vv3) - vv);

                    itnearest = it;
                }
            }
            {
                let old = it;
                it = it + 1;
                old
            };
        }
    }

    if (itnearest == (-1i32)) {
        err!(
            "x=%f %f %f not in any tetra: vdiff=%f\n",
            ((*(x).offset((0i32) as isize)) as f64),
            ((*(x).offset((1i32) as isize)) as f64),
            ((*(x).offset((2i32) as isize)) as f64),
            (vdiff as f64)
        );
    }

    return itnearest;
}

// crossprod (sutetraray.c:1959).

unsafe fn crossprod(
    mut u1: f32,
    mut u2: f32,
    mut u3: f32,
    mut v1: f32,
    mut v2: f32,
    mut v3: f32,
    mut w1: *mut f32,
    mut w2: *mut f32,
    mut w3: *mut f32,
) {
    (*w1) = ((u2 * u3) - (v2 * v3));

    (*w2) = ((u3 * v1) - (u1 * v3));

    (*w3) = ((u1 * v2) - (v1 * u2));

    if (((((*w1) as f64) == 0f64) && (((*w2) as f64) == 0f64)) && (((*w3) as f64) == 0f64)) {
        (*w3) = (1f64 as f32);
    }
}

// tri_multi_intp (sutetraray.c:1985).

unsafe fn tri_multi_intp(
    mut ray0: *mut RAY,
    mut ray1: *mut RAY,
    mut ray2: *mut RAY,
    mut fxgd: f32,
    mut fygd: f32,
    mut dxgd: f32,
    mut dygd: f32,
    mut nxgd: i32,
    mut nygd: i32,
    mut ttable: *mut *mut f32,
    mut ctable: *mut *mut f32,
    mut rtable: *mut *mut f32,
) {
    let mut n1: f32 = 0.0;
    let mut n2: f32 = 0.0;
    let mut n3: f32 = 0.0;
    /* the normal to this triangle */
    let mut ixgd: i32 = 0;
    let mut iygd: i32 = 0;

    let mut ttemp: f32 = 0.0;
    let mut ctemp: f32 = 0.0;
    let mut rtemp: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut y: f32 = 0.0;

    let mut tmin: f32 = 0.0;
    let mut tmax: f32 = 0.0;

    let mut xmin: f32 = 0.0;
    let mut xmax: f32 = 0.0;

    let mut ymin: f32 = 0.0;
    let mut ymax: f32 = 0.0;

    let mut areat: f32 = 0.0;

    let mut sintheta: f32 = 0.0;
    let mut costheta: f32 = 0.0;

    let mut v0: f32 = 0.0;
    let mut v1: f32 = 0.0;
    let mut v2: f32 = 0.0;

    let mut r0: f32 = 0.0;
    let mut r1: f32 = 0.0;
    let mut r2: f32 = 0.0;

    let mut c0: f32 = 0.0;
    let mut c1: f32 = 0.0;
    let mut c2: f32 = 0.0;

    let mut x0: f32 = 0.0;
    let mut y0: f32 = 0.0;
    let mut t0: f32 = 0.0;

    let mut x1: f32 = 0.0;
    let mut y1: f32 = 0.0;
    let mut t1: f32 = 0.0;

    let mut x2: f32 = 0.0;
    let mut y2: f32 = 0.0;
    let mut t2: f32 = 0.0;

    let mut dist0: f32 = 0.0;
    let mut dist1: f32 = 0.0;
    let mut dist2: f32 = 0.0;

    v0 = (*ray0).v;

    v1 = (*ray1).v;

    v2 = (*ray2).v;

    r0 = (*ray0).r;

    r1 = (*ray1).r;

    r2 = (*ray2).r;

    c0 = (v0 * (*(std::ptr::addr_of_mut!((*ray0).p).cast::<f32>()).offset((2i32) as isize)));

    c1 = (v1 * (*(std::ptr::addr_of_mut!((*ray1).p).cast::<f32>()).offset((2i32) as isize)));

    c2 = (v2 * (*(std::ptr::addr_of_mut!((*ray2).p).cast::<f32>()).offset((2i32) as isize)));

    x0 = (*(std::ptr::addr_of_mut!((*ray0).x).cast::<f32>()).offset((0i32) as isize));

    y0 = (*(std::ptr::addr_of_mut!((*ray0).x).cast::<f32>()).offset((1i32) as isize));

    t0 = (*ray0).ttotal;

    x1 = (*(std::ptr::addr_of_mut!((*ray1).x).cast::<f32>()).offset((0i32) as isize));

    y1 = (*(std::ptr::addr_of_mut!((*ray1).x).cast::<f32>()).offset((1i32) as isize));

    t1 = (*ray1).ttotal;

    x2 = (*(std::ptr::addr_of_mut!((*ray2).x).cast::<f32>()).offset((0i32) as isize));

    y2 = (*(std::ptr::addr_of_mut!((*ray2).x).cast::<f32>()).offset((1i32) as isize));

    t2 = (*ray2).ttotal;

    areat = ((3f64 * (area(x0, y0, x1, y1, x2, y2) as f64)) as f32);

    if ((areat as f64) < 1.0000000000000001E-5f64) {
        return;
    }

    tmin = (if ((t0) < (if ((t1) < (t2)) { (t1) } else { (t2) })) {
        (t0)
    } else {
        (if ((t1) < (t2)) { (t1) } else { (t2) })
    });

    tmax = (if ((t0) > (if ((t1) > (t2)) { (t1) } else { (t2) })) {
        (t0)
    } else {
        (if ((t1) > (t2)) { (t1) } else { (t2) })
    });

    xmin = (if ((x0) < (if ((x1) < (x2)) { (x1) } else { (x2) })) {
        (x0)
    } else {
        (if ((x1) < (x2)) { (x1) } else { (x2) })
    });

    xmax = (if ((x0) > (if ((x1) > (x2)) { (x1) } else { (x2) })) {
        (x0)
    } else {
        (if ((x1) > (x2)) { (x1) } else { (x2) })
    });

    ymin = (if ((y0) < (if ((y1) < (y2)) { (y1) } else { (y2) })) {
        (y0)
    } else {
        (if ((y1) < (y2)) { (y1) } else { (y2) })
    });

    ymax = (if ((y0) > (if ((y1) > (y2)) { (y1) } else { (y2) })) {
        (y0)
    } else {
        (if ((y1) > (y2)) { (y1) } else { (y2) })
    });

    crossprod(
        (x1 - x0),
        (y1 - y0),
        (t1 - t0),
        (x2 - x0),
        (y2 - y0),
        (t2 - t0),
        std::ptr::addr_of_mut!(n1),
        std::ptr::addr_of_mut!(n2),
        std::ptr::addr_of_mut!(n3),
    );

    {
        ixgd = ((if (((0i32) as f32) > (((xmin - fxgd) / dxgd) - (1i32 as f32))) {
            ((0i32) as f32)
        } else {
            (((xmin - fxgd) / dxgd) - (1i32 as f32))
        }) as i32);
        loop {
            if !((ixgd as f32)
                <= (if ((((nxgd).wrapping_sub(1i32)) as f32)
                    < (((xmax - fxgd) / dxgd) + (1i32 as f32)))
                {
                    (((nxgd).wrapping_sub(1i32)) as f32)
                } else {
                    (((xmax - fxgd) / dxgd) + (1i32 as f32))
                }))
            {
                break;
            }
            {
                x = (fxgd + ((ixgd as f32) * dxgd));

                {
                    iygd = ((if (((0i32) as f32) > (((ymin - fygd) / dygd) - (1i32 as f32))) {
                        ((0i32) as f32)
                    } else {
                        (((ymin - fygd) / dygd) - (1i32 as f32))
                    }) as i32);
                    loop {
                        if !((iygd as f32)
                            <= (if ((((nygd).wrapping_sub(1i32)) as f32)
                                < (((ymax - fygd) / dygd) + (1i32 as f32)))
                            {
                                (((nygd).wrapping_sub(1i32)) as f32)
                            } else {
                                (((ymax - fygd) / dygd) + (1i32 as f32))
                            }))
                        {
                            break;
                        }
                        {
                            y = (fygd + ((iygd as f32) * dygd));

                            if ((WF_NO as u32)
                                == (in_tri(x0, y0, x1, y1, x2, y2, x, y, areat) as u32))
                            {
                                {
                                    let old = iygd;
                                    iygd = iygd + 1;
                                    old
                                };
                                continue;
                            }

                            dist0 = (((((x - x0) * (x - x0)) + ((y - y0) * (y - y0))) as f64).sqrt()
                                as f32);

                            dist1 = (((((x - x1) * (x - x1)) + ((y - y1) * (y - y1))) as f64).sqrt()
                                as f32);

                            dist2 = (((((x - x2) * (x - x2)) + ((y - y2) * (y - y2))) as f64).sqrt()
                                as f32);

                            if ((dist0 <= dist1) && (dist0 <= dist2)) {
                                ttemp = zhorz(x, y, n1, n2, n3, x0, y0, t0);

                                if (((ttemp as f64) < ((tmin as f64) - 1.0000000000000001E-5f64))
                                    || ((ttemp as f64)
                                        > ((tmax as f64) + 1.0000000000000001E-5f64)))
                                {
                                    {
                                        let old = iygd;
                                        iygd = iygd + 1;
                                        old
                                    };
                                    continue;
                                }

                                sintheta = (((dist0 as f64)
                                    / (if ((0.10000000000000001f64) > ((r0) as f64)) {
                                        (0.10000000000000001f64)
                                    } else {
                                        ((r0) as f64)
                                    })) as f32);

                                sintheta = ((if ((1f64) < ((sintheta) as f64)) {
                                    (1f64)
                                } else {
                                    ((sintheta) as f64)
                                }) as f32);

                                costheta = ((1f64 - ((sintheta * sintheta) as f64)).sqrt() as f32);

                                ttemp = (((ttemp as f64)
                                    - (((r0 as f64) * (1f64 - (costheta as f64))) / (v0 as f64)))
                                    as f32);

                                ctemp = zhorz(x, y, n1, n2, n3, x0, y0, c0);

                                rtemp = zhorz(x, y, n1, n2, n3, x0, y0, r0);
                            } else {
                                if ((dist1 <= dist0) && (dist1 <= dist2)) {
                                    ttemp = zhorz(x, y, n1, n2, n3, x1, y1, t1);

                                    if (((ttemp as f64)
                                        < ((tmin as f64) - 1.0000000000000001E-5f64))
                                        || ((ttemp as f64)
                                            > ((tmax as f64) + 1.0000000000000001E-5f64)))
                                    {
                                        {
                                            let old = iygd;
                                            iygd = iygd + 1;
                                            old
                                        };
                                        continue;
                                    }

                                    sintheta = (((dist1 as f64)
                                        / (if ((0.10000000000000001f64) > ((r1) as f64)) {
                                            (0.10000000000000001f64)
                                        } else {
                                            ((r1) as f64)
                                        })) as f32);

                                    sintheta = ((if ((1f64) < ((sintheta) as f64)) {
                                        (1f64)
                                    } else {
                                        ((sintheta) as f64)
                                    }) as f32);

                                    costheta =
                                        ((1f64 - ((sintheta * sintheta) as f64)).sqrt() as f32);

                                    ttemp = (((ttemp as f64)
                                        - (((r1 as f64) * (1f64 - (costheta as f64)))
                                            / (v1 as f64)))
                                        as f32);

                                    ctemp = zhorz(x, y, n1, n2, n3, x1, y1, c1);

                                    rtemp = zhorz(x, y, n1, n2, n3, x1, y1, r1);
                                } else {
                                    ttemp = zhorz(x, y, n1, n2, n3, x2, y2, t2);

                                    if (((ttemp as f64)
                                        < ((tmin as f64) - 1.0000000000000001E-5f64))
                                        || ((ttemp as f64)
                                            > ((tmax as f64) + 1.0000000000000001E-5f64)))
                                    {
                                        {
                                            let old = iygd;
                                            iygd = iygd + 1;
                                            old
                                        };
                                        continue;
                                    }

                                    sintheta = (((dist2 as f64)
                                        / (if ((0.10000000000000001f64) > ((r2) as f64)) {
                                            (0.10000000000000001f64)
                                        } else {
                                            ((r2) as f64)
                                        })) as f32);

                                    sintheta = ((if ((1f64) < ((sintheta) as f64)) {
                                        (1f64)
                                    } else {
                                        ((sintheta) as f64)
                                    }) as f32);

                                    costheta =
                                        ((1f64 - ((sintheta * sintheta) as f64)).sqrt() as f32);

                                    ttemp = (((ttemp as f64)
                                        - (((r2 as f64) * (1f64 - (costheta as f64)))
                                            / (v2 as f64)))
                                        as f32);

                                    ctemp = zhorz(x, y, n1, n2, n3, x2, y2, c2);

                                    rtemp = zhorz(x, y, n1, n2, n3, x2, y2, r2);
                                }
                            }

                            (*(*(ttable).offset((iygd) as isize)).offset((ixgd) as isize)) =
                                (if ((*(*(ttable).offset((iygd) as isize)).offset((ixgd) as isize))
                                    < (ttemp))
                                {
                                    (*(*(ttable).offset((iygd) as isize)).offset((ixgd) as isize))
                                } else {
                                    (ttemp)
                                });

                            (*(*(ctable).offset((iygd) as isize)).offset((ixgd) as isize)) =
                                ((if ((0f64)
                                    > (if ((1f64) < ((ctemp) as f64)) {
                                        (1f64)
                                    } else {
                                        ((ctemp) as f64)
                                    }))
                                {
                                    (0f64)
                                } else {
                                    (if ((1f64) < ((ctemp) as f64)) {
                                        (1f64)
                                    } else {
                                        ((ctemp) as f64)
                                    })
                                }) as f32);

                            (*(*(rtable).offset((iygd) as isize)).offset((ixgd) as isize)) =
                                ((if ((0.10000000000000001f64)
                                    > ((if (((1000i32) as f32) < (rtemp)) {
                                        ((1000i32) as f32)
                                    } else {
                                        (rtemp)
                                    }) as f64))
                                {
                                    (0.10000000000000001f64)
                                } else {
                                    ((if (((1000i32) as f32) < (rtemp)) {
                                        ((1000i32) as f32)
                                    } else {
                                        (rtemp)
                                    }) as f64)
                                }) as f32);
                        }
                        {
                            let old = iygd;
                            iygd = iygd + 1;
                            old
                        };
                    }
                }
            }
            {
                let old = ixgd;
                ixgd = ixgd + 1;
                old
            };
        }
    }
}

// in_tri (sutetraray.c:2138).

unsafe fn in_tri(
    mut x0: f32,
    mut y0: f32,
    mut x1: f32,
    mut y1: f32,
    mut x2: f32,
    mut y2: f32,
    mut x: f32,
    mut y: f32,
    mut areat: f32,
) -> i32 {
    let mut area0: f32 = 0.0;
    let mut area1: f32 = 0.0;
    let mut area2: f32 = 0.0;

    area0 = area(x, y, x1, y1, x2, y2);

    if (area0 > areat) {
        return WF_NO;
    }

    area1 = area(x0, y0, x, y, x2, y2);

    if ((area0 + area1) > areat) {
        return WF_NO;
    }

    area2 = area(x0, y0, x1, y1, x, y);

    if ((((area0 + area1) + area2) as f64) > ((areat as f64) + 0.10000000000000001f64)) {
        return WF_NO;
    }

    return WF_YES;
}

// tri_intp (sutetraray.c:2173).

unsafe fn tri_intp(
    mut x: f32,
    mut y: f32,
    mut x0: f32,
    mut y0: f32,
    mut t0: f32,
    mut x1: f32,
    mut y1: f32,
    mut t1: f32,
    mut x2: f32,
    mut y2: f32,
    mut t2: f32,
) -> f32 {
    let mut n1: f32 = 0.0;
    let mut n2: f32 = 0.0;
    let mut n3: f32 = 0.0;
    /* the normal to this triangle */
    crossprod(
        (x1 - x0),
        (y1 - y0),
        (t1 - t0),
        (x2 - x0),
        (y2 - y0),
        (t2 - t0),
        std::ptr::addr_of_mut!(n1),
        std::ptr::addr_of_mut!(n2),
        std::ptr::addr_of_mut!(n3),
    );

    return (t0 - ((((x - x0) * n1) + ((y - y0) * n2)) / n3));
}

// read_point (sutetraray.c:2214).

unsafe fn read_point(state: *mut State, mut point: *mut POINT, mut npoint: i32) {
    let mut ip: i32 = 0;
    let mut ipp: i32 = 0;
    /* index for npoints */
    {
        ip = 0i32;
        loop {
            if !(ip < npoint) {
                break;
            }
            {
                {
                    let mut scanned = 0i32;
                    'scan: {
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(ipp) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).x)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).x)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).x)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).s) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).n)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).n)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).n)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "\n") {
                            break 'scan;
                        }
                    }
                    scanned
                };
            }
            {
                let old = ip;
                ip = ip + 1;
                old
            };
        }
    }
}

// read_facet (sutetraray.c:2240).

unsafe fn read_facet(state: *mut State, mut facet: *mut FACET, mut nfacet: i32) {
    let mut ifs: i32 = 0;
    let mut ifss: i32 = 0;
    /* index for nfacets */
    {
        ifs = 0i32;
        loop {
            if !(ifs < nfacet) {
                break;
            }
            {
                {
                    let mut scanned = 0i32;
                    'scan: {
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(ifss) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).area) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((3i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((4i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((5i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((6i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((7i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((8i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ct)
                                    .cast::<f32>())
                                .offset((9i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ip)
                                    .cast::<i32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ip)
                                    .cast::<i32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).ip)
                                    .cast::<i32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).cn)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).cn)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(facet).offset((ifs) as isize)).cn)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!(
                                    (*(facet).offset((ifs) as isize)).itetra
                                )
                                .cast::<i32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!(
                                    (*(facet).offset((ifs) as isize)).itetra
                                )
                                .cast::<i32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "\n") {
                            break 'scan;
                        }
                    }
                    scanned
                };
            }
            {
                let old = ifs;
                ifs = ifs + 1;
                old
            };
        }
    }
}

// read_tetra (sutetraray.c:2277).

unsafe fn read_tetra(state: *mut State, mut tetra: *mut TETRA, mut ntetra: i32) {
    let mut it: i32 = 0;
    let mut itt: i32 = 0;
    /* index for ntetra */
    {
        it = 0i32;
        loop {
            if !(it < ntetra) {
                break;
            }
            {
                {
                    let mut scanned = 0i32;
                    'scan: {
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(itt) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).v) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                                    .cast::<i32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                                    .cast::<i32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                                    .cast::<i32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ip)
                                    .cast::<i32>())
                                .offset((3i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!(
                                    (*(tetra).offset((it) as isize)).ifacet
                                )
                                .cast::<i32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!(
                                    (*(tetra).offset((it) as isize)).ifacet
                                )
                                .cast::<i32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!(
                                    (*(tetra).offset((it) as isize)).ifacet
                                )
                                .cast::<i32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!(
                                    (*(tetra).offset((it) as isize)).ifacet
                                )
                                .cast::<i32>())
                                .offset((3i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_i32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).ireg) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " : ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).gs)
                                    .cast::<f32>())
                                .offset((0i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).gs)
                                    .cast::<f32>())
                                .offset((1i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), " ") {
                            break 'scan;
                        }
                        if let Some(v) = fscanf_f32(&*std::ptr::addr_of_mut!((*state).input)) {
                            *std::ptr::addr_of_mut!(
                                (*(std::ptr::addr_of_mut!((*(tetra).offset((it) as isize)).gs)
                                    .cast::<f32>())
                                .offset((2i32) as isize))
                            ) = v;
                            scanned += 1;
                        } else {
                            break 'scan;
                        }
                        if !scan_literal(&*std::ptr::addr_of_mut!((*state).input), "\n") {
                            break 'scan;
                        }
                    }
                    scanned
                };
            }
            {
                let old = it;
                it = it + 1;
                old
            };
        }
    }
}

// s_intp (sutetraray.c:2302).

unsafe fn s_intp(mut point: *mut POINT, mut tetra: *mut TETRA, mut x: *mut f32) -> f32 {
    let mut ss: f32 = 0.0;

    ss = ((((*(point).offset(
        (*(std::ptr::addr_of_mut!((*tetra).ip).cast::<i32>()).offset((0i32) as isize)) as isize,
    ))
    .s + ((*(std::ptr::addr_of_mut!((*tetra).gs).cast::<f32>())
        .offset((0i32) as isize))
        * ((*(x).offset((0i32) as isize))
            - (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*tetra).ip).cast::<i32>()).offset((0i32) as isize))
                        as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((0i32) as isize)))))
        + ((*(std::ptr::addr_of_mut!((*tetra).gs).cast::<f32>()).offset((1i32) as isize))
            * ((*(x).offset((1i32) as isize))
                - (*(std::ptr::addr_of_mut!(
                    (*(point).offset(
                        (*(std::ptr::addr_of_mut!((*tetra).ip).cast::<i32>())
                            .offset((0i32) as isize)) as isize
                    ))
                    .x
                )
                .cast::<f32>())
                .offset((1i32) as isize)))))
        + ((*(std::ptr::addr_of_mut!((*tetra).gs).cast::<f32>()).offset((2i32) as isize))
            * ((*(x).offset((2i32) as isize))
                - (*(std::ptr::addr_of_mut!(
                    (*(point).offset(
                        (*(std::ptr::addr_of_mut!((*tetra).ip).cast::<i32>())
                            .offset((0i32) as isize)) as isize
                    ))
                    .x
                )
                .cast::<f32>())
                .offset((2i32) as isize)))));

    return ss;
}

// trace1tetra (sutetraray.c:2329).

unsafe fn trace1tetra(
    mut maxnsegs: i32,
    mut xmin: f32,
    mut ymin: f32,
    mut xmax: f32,
    mut ymax: f32,
    mut zmax: f32,
    mut ray: *mut RAY,
    mut npoint: i32,
    mut point: *mut POINT,
    mut nfacet: i32,
    mut facet: *mut FACET,
    mut ntetra: i32,
    mut tetra: *mut TETRA,
    mut irefseq: *mut i32,
    mut nrefseq: i32,
    mut jpfp: *mut CFile,
    mut dtwf: f32,
) -> i32 {
    let mut sgdt: f32 = 0.0;
    /* sigma value when intersecting with it*dtwf face */
    let mut sgf4min: f32 = 0.0;
    /* minimum sigma in the 4 solutions */
    let mut ifhit: i32 = 0;
    /* the min facet */
    let mut sgf4: [f32; 4] = [0.0; 4];

    let mut xst: [f32; 3] = [0.0; 3];
    /* source position */
    let mut sx: f32 = 0.0;
    let mut sy: f32 = 0.0;
    let mut sz: f32 = 0.0;
    /* sloth gradient */
    let mut smst: f32 = 0.0;
    /* smallest sigma */
    let mut smst0: f32 = 0.0;
    /* almost the smallest sigma */
    let mut sursg: f32 = 0.0;
    /* sigma needed to hit the surface */
    let mut sincid: f32 = 0.0;

    let mut x0: f32 = (9000i32 as f32);

    let mut y0: f32 = (9000i32 as f32);

    let mut z0: f32 = (9000i32 as f32);

    let mut ip: i32 = 0;

    let mut cn0: f32 = 0.0;
    let mut cn1: f32 = 0.0;
    let mut cn2: f32 = 0.0;

    let mut dist: f32 = 0.0;
    let mut distmax: f32 = (0f64 as f32);

    let mut x00: f32 = 0.0;
    let mut y00: f32 = 0.0;
    let mut z00: f32 = 0.0;

    let mut ip0: i32 = 0;
    let mut never: i32 = 0;
    let mut never0: i32 = 0;

    let mut itetranew: i32 = 0;

    let mut ifacet: i32 = 0;
    let mut if0: i32 = 0;

    let mut tadjust: f32 = 0.0;

    let mut s: f32 = 0.0;
    let mut dlt: f32 = 0.0;

    let mut xcenter: f32 = 0.0;
    let mut ycenter: f32 = 0.0;
    let mut zcenter: f32 = 0.0;

    let mut alpha: f32 = 0.0;
    let mut alpha0: f32 = 0.0;

    if (((*ray).itetra > (ntetra).wrapping_sub(1i32)) || ((*ray).itetra < 0i32)) {
        err!("ray->itetra out of range");
    }

    s = s_intp(
        point,
        std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize))),
        std::ptr::addr_of_mut!((*ray).x).cast::<f32>(),
    );
    /* the current position */
    s_normalize(std::ptr::addr_of_mut!((*ray).p).cast::<f32>(), s);

    if (((*(tetra).offset(((*ray).itetra) as isize)).v as f64) < 0f64) {
        cross_slim_tetra(ray, point, facet, tetra, irefseq, nrefseq, s, jpfp);
        /* s needs to be reevaluated, because ray is changed */
        s = s_intp(
            point,
            std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize))),
            std::ptr::addr_of_mut!((*ray).x).cast::<f32>(),
        );
        /* the current position */
        s_normalize(std::ptr::addr_of_mut!((*ray).p).cast::<f32>(), s);
    }

    (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((0i32) as isize)) =
        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize));

    (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((1i32) as isize)) =
        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize));

    (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize)) =
        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((2i32) as isize));

    (*ray).v = ((1f64 / (s as f64).sqrt()) as f32);
    /*************************************************************
    We need to calculate sgdt. The traveltime needed to reach the
    next wavefront (from ray->t to dtwf).
    *************************************************************/
    /* sx,sy,sz is the sloth gradient */
    sx = (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).gs).cast::<f32>())
        .offset((0i32) as isize));

    sy = (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).gs).cast::<f32>())
        .offset((1i32) as isize));

    sz = (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).gs).cast::<f32>())
        .offset((2i32) as isize));

    sgdt = cubicsolver(
        ((((((sx * sx) + (sy * sy)) + (sz * sz)) as f64) / 12f64) as f32),
        ((0.5f64
            * ((((sx
                * (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)))
                + (sy
                    * (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize))))
                + (sz
                    * (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize))))
                as f64)) as f32),
        s,
        ((*ray).t - dtwf),
    );
    /*****************************************************
    The plane is:

         cn0*x+cn1*y+cn2*z=0,

    Suppose (x0,y0,z0) is a point on this plane.

    Then the linear transform from the old cartesian (x,y,z)
    to the now (x',y',z'), where cn0,cn1 are in plane, and
    cn2 is the normal to the plane, (x0,y0,z0) is the
    origin of the new cordinates, the transform can be written as:

       x'       x - x0
      (y') = T (y - y0),	   (*)
       z'       z - z0

    on the other hand, the ray can be written as

      x=xst + px*sigma + 0.25*sx*sigma^2,
      y=yst + py*sigma + 0.25*sy*sigma^2,
      z=zst + pz*sigma + 0.25*sz*sigma^2,

    Take the third equation of (*) specifically,

    z' = cn0*(x-x0) + cn1*(y-y0) + cn2*(z-z0)

    When the ray incides, z'=0: =>

      0.25*(cn0*sx+cn1*sy+cn2*sz)*sg^2
      +
      (cn0*px+cn1*py+cn2*pz)*sg
      +
      cn0*(xst-x0)+cn1*(ys-y0)+cn2*(zs-z0)=0

    Solve for sghorz => the intersection point.
    *****************************************************/
    {
        if0 = 0i32;
        loop {
            if !(if0 < 4i32) {
                break;
            }
            {
                /* for all facets of this tetra */
                ifacet =
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ifacet)
                        .cast::<i32>())
                    .offset((if0) as isize));

                if (ifacet > (nfacet).wrapping_sub(1i32)) {
                    err!("ifacet out of range");
                }

                if (ifacet == (*ray).lastfacet) {
                    (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize)) =
                        (1.0E+10f64 as f32);

                    {
                        let old = if0;
                        if0 = if0 + 1;
                        old
                    };
                    continue;
                }
                /* the normal to this facet */
                cn0 = (*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).cn)
                    .cast::<f32>())
                .offset((0i32) as isize));

                cn1 = (*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).cn)
                    .cast::<f32>())
                .offset((1i32) as isize));

                cn2 = (*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).cn)
                    .cast::<f32>())
                .offset((2i32) as isize));
                /* find the farthest vertex of this facet to the current point xst */
                distmax = (0f64 as f32);

                {
                    ip0 = 0i32;
                    loop {
                        if !(ip0 < 3i32) {
                            break;
                        }
                        {
                            ip = (*(std::ptr::addr_of_mut!(
                                (*(facet).offset((ifacet) as isize)).ip
                            )
                            .cast::<i32>())
                            .offset((ip0) as isize));

                            if (ip > (npoint).wrapping_sub(1i32)) {
                                err!("ip out of range");
                            }

                            x00 = (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).x)
                                .cast::<f32>())
                            .offset((0i32) as isize));

                            y00 = (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).x)
                                .cast::<f32>())
                            .offset((1i32) as isize));

                            z00 = (*(std::ptr::addr_of_mut!((*(point).offset((ip) as isize)).x)
                                .cast::<f32>())
                            .offset((2i32) as isize));

                            dist = ((((x00
                                - (*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                    .offset((0i32) as isize)))
                                * (x00
                                    - (*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                        .offset((0i32) as isize))))
                                + ((y00
                                    - (*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                        .offset((1i32) as isize)))
                                    * (y00
                                        - (*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                            .offset((1i32) as isize)))))
                                + ((z00
                                    - (*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                        .offset((2i32) as isize)))
                                    * (z00
                                        - (*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                            .offset((2i32) as isize)))));

                            if (dist > distmax) {
                                x0 = x00;

                                y0 = y00;

                                z0 = z00;

                                distmax = dist;
                            }
                        }
                        {
                            let old = ip0;
                            ip0 = ip0 + 1;
                            old
                        };
                    }
                }

                (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize)) = quadsolver(
                    ((0.25f64 * ((((cn0 * sx) + (cn1 * sy)) + (cn2 * sz)) as f64)) as f32),
                    (((cn0
                        * (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>())
                            .offset((0i32) as isize)))
                        + (cn1
                            * (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>())
                                .offset((1i32) as isize))))
                        + (cn2
                            * (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>())
                                .offset((2i32) as isize)))),
                    (((cn0
                        * ((*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                            .offset((0i32) as isize))
                            - x0))
                        + (cn1
                            * ((*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                .offset((1i32) as isize))
                                - y0)))
                        + (cn2
                            * ((*(std::ptr::addr_of_mut!(xst).cast::<f32>())
                                .offset((2i32) as isize))
                                - z0))),
                );

                if ((*(tetra).offset(
                    (*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).itetra)
                        .cast::<i32>())
                    .offset((0i32) as isize)) as isize,
                ))
                .ireg
                    == (*(tetra).offset(
                        (*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).itetra)
                            .cast::<i32>())
                        .offset((1i32) as isize)) as isize,
                    ))
                    .ireg)
                {
                    {
                        let old = if0;
                        if0 = if0 + 1;
                        old
                    };
                    continue;
                }

                if ((((((((*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                    .cast::<f32>())
                .offset((0i32) as isize)) as f64)
                    == 0f64)
                    && (((*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                        .cast::<f32>())
                    .offset((1i32) as isize)) as f64)
                        == 0f64))
                    && (((*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                        .cast::<f32>())
                    .offset((2i32) as isize)) as f64)
                        == 0f64))
                    && (((*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                        .cast::<f32>())
                    .offset((3i32) as isize)) as f64)
                        == 0f64))
                    && (((*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                        .cast::<f32>())
                    .offset((4i32) as isize)) as f64)
                        == 0f64))
                    && (((*(std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                        .cast::<f32>())
                    .offset((5i32) as isize)) as f64)
                        == 0f64))
                {
                    {
                        let old = if0;
                        if0 = if0 + 1;
                        old
                    };
                    continue;
                }

                (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize)) =
                    tetrasolver(
                        std::ptr::addr_of_mut!((*(facet).offset((ifacet) as isize)).ct)
                            .cast::<f32>(),
                        std::ptr::addr_of_mut!((*ray).p).cast::<f32>(),
                        std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).gs)
                            .cast::<f32>(),
                        std::ptr::addr_of_mut!(xst).cast::<f32>(),
                        (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize)),
                    );
            }
            {
                let old = if0;
                if0 = if0 + 1;
                old
            };
        }
    }

    sgf4min = (1.0E+10f64 as f32);

    never0 = 1i32;

    never = (-1i32);

    if ((*ray).shootupward == 1i32) {
        ifhit = 0i32;

        {
            if0 = 0i32;
            loop {
                if !(if0 < 4i32) {
                    break;
                }
                {
                    if (((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize))
                        as f64)
                        == 0f64)
                    {
                        never0 = 0i32;

                        never = if0;
                    }

                    if (((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize))
                        < sgf4min)
                        && (((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize))
                            as f64)
                            > 1.0000000000000001E-5f64))
                    {
                        sgf4min =
                            (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize));

                        ifhit = (*(std::ptr::addr_of_mut!(
                            (*(tetra).offset(((*ray).itetra) as isize)).ifacet
                        )
                        .cast::<i32>())
                        .offset((if0) as isize));
                    }
                }
                {
                    let old = if0;
                    if0 = if0 + 1;
                    old
                };
            }
        }
    } else {
        ifhit = 1i32;
        /* do not go over lower boundary */
        if ((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize)) < zmax) {
            (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((3i32) as isize)) =
                (if ((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((3i32) as isize))
                    < (quadsolver(
                        ((0.25f64 * (sz as f64)) as f32),
                        (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)),
                        ((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize))
                            - zmax),
                    )))
                {
                    (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((3i32) as isize))
                } else {
                    (quadsolver(
                        ((0.25f64 * (sz as f64)) as f32),
                        (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)),
                        ((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize))
                            - zmax),
                    ))
                });
        }

        {
            if0 = 1i32;
            loop {
                if !(if0 < 4i32) {
                    break;
                }
                {
                    if (((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize))
                        as f64)
                        == 0f64)
                    {
                        never0 = 0i32;

                        never = if0;
                    }

                    if (((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize))
                        < sgf4min)
                        && (((*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize))
                            as f64)
                            > 1.0000000000000001E-5f64))
                    {
                        sgf4min =
                            (*(std::ptr::addr_of_mut!(sgf4).cast::<f32>()).offset((if0) as isize));

                        ifhit = (*(std::ptr::addr_of_mut!(
                            (*(tetra).offset(((*ray).itetra) as isize)).ifacet
                        )
                        .cast::<i32>())
                        .offset((if0) as isize));
                    }
                }
                {
                    let old = if0;
                    if0 = if0 + 1;
                    old
                };
            }
        }
    }

    if (((never0 == 0i32) && ((sgf4min as f64) > (1.0E+10f64 / 2f64))) && (never >= 0i32)) {
        sgf4min = (0f64 as f32);

        ifhit = never;
    }

    smst0 = (if ((sgdt) < (sgf4min)) {
        (sgdt)
    } else {
        (sgf4min)
    });
    /* because some tetra intersect the boundary, so it is possible
    for a ray to hit the surface without hitting a facet on
    the surface. If we suppose the surface is z=0 plane, then
    the following trick solves this problem */
    if ((*ray).shootupward == 1i32) {
        sursg = quadsolver(
            ((0.25f64 * (sz as f64)) as f32),
            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)),
            (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize)),
        );
        /* sg needed to hit surface */
        smst = (if ((smst0) < (sursg)) {
            (smst0)
        } else {
            (sursg)
        });
    } else {
        smst = smst0;

        sursg = (1.0E+10f64 as f32);
    }

    if (((sgf4min as f64) > (1.0E+10f64 * 0.5f64)) && ((*ray).ionce >= 1000i32)) {
        fprintf!(&*jpfp, "Could not find any facet\n");
    }

    if (((((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) < xmin)
        || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) > xmax))
        || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) < ymin))
        || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) > ymax))
    {
        fprintf!(&*jpfp, "this ray out of model\n");

        (*ray).icode = SPOILED;

        return OUT_OF_MODEL;
    }

    if ((((smst as f64) > (1.0E+10f64 * 0.5f64)) || ((distmax as f64) < 1.0000000000000001E-5f64))
        && ((*ray).ionce < 1000i32))
    {
        (*ray).itetra = in_which_tetra(
            xmin,
            ymin,
            xmax,
            ymax,
            zmax,
            std::ptr::addr_of_mut!((*ray).x).cast::<f32>(),
            point,
            tetra,
            ntetra,
        );

        xcenter = (((((((*(std::ptr::addr_of_mut!(
            (*(point).offset(
                (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                    .cast::<i32>())
                .offset((0i32) as isize)) as isize
            ))
            .x
        )
        .cast::<f32>())
        .offset((0i32) as isize))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((1i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((0i32) as isize)))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((2i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((0i32) as isize)))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((3i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((0i32) as isize))) as f64)
            / 4f64) as f32);

        ycenter = (((((((*(std::ptr::addr_of_mut!(
            (*(point).offset(
                (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                    .cast::<i32>())
                .offset((0i32) as isize)) as isize
            ))
            .x
        )
        .cast::<f32>())
        .offset((1i32) as isize))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((1i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((1i32) as isize)))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((2i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((1i32) as isize)))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((3i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((1i32) as isize))) as f64)
            / 4f64) as f32);

        zcenter = (((((((*(std::ptr::addr_of_mut!(
            (*(point).offset(
                (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                    .cast::<i32>())
                .offset((0i32) as isize)) as isize
            ))
            .x
        )
        .cast::<f32>())
        .offset((2i32) as isize))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((1i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((2i32) as isize)))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((2i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((2i32) as isize)))
            + (*(std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize)).ip)
                        .cast::<i32>())
                    .offset((3i32) as isize)) as isize
                ))
                .x
            )
            .cast::<f32>())
            .offset((2i32) as isize))) as f64)
            / 4f64) as f32);

        (*ray).ionce = (*ray).ionce.wrapping_add(1);

        alpha = ((1f64 - ((((*ray).ionce as f64) * 0.5f64) / ((1000i32 as f32) as f64))) as f32);

        alpha0 = ((1f64 - (alpha as f64)) as f32);

        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) = ((alpha
            * (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((0i32) as isize)))
            + (alpha0 * xcenter));

        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) = ((alpha
            * (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((1i32) as isize)))
            + (alpha0 * ycenter));

        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((2i32) as isize)) = ((alpha
            * (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize)))
            + (alpha0 * zcenter));

        if (((((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize))
            < xmin)
            || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize))
                > xmax))
            || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize))
                < ymin))
            || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) > ymax))
        {
            fprintf!(&*jpfp, "this ray out of model\n");

            (*ray).icode = SPOILED;

            return OUT_OF_MODEL;
        }

        (*ray).lastfacet = (-1i32);

        return HIT_FACET_FIRST;
    }

    if ((*ray).ionce > 0i32) {
        (*ray).nseg = (*ray).nseg.wrapping_add(1);

        if ((*ray).nseg > maxnsegs) {
            err!("maxnsegs too small");
        }

        (*(*((*ray).xseg).offset(((*ray).nseg) as isize)).offset((0i32) as isize)) =
            (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((0i32) as isize));

        (*(*((*ray).xseg).offset(((*ray).nseg) as isize)).offset((1i32) as isize)) =
            (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((1i32) as isize));

        (*(*((*ray).xseg).offset(((*ray).nseg) as isize)).offset((2i32) as isize)) =
            (*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize));

        if ((*ray).shootupward == 1i32) {
            tadjust = ((-((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize))
                - (*(*((*ray).xseg).offset((((*ray).nseg).wrapping_sub(1i32)) as isize))
                    .offset((2i32) as isize))))
                / (*ray).v);
        } else {
            tadjust = (((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize))
                - (*(*((*ray).xseg).offset((((*ray).nseg).wrapping_sub(1i32)) as isize))
                    .offset((2i32) as isize)))
                / (*ray).v);
        }

        (*((*ray).tseg).offset(((*ray).nseg) as isize)) =
            (((*((*ray).tseg).offset(((*ray).nseg) as isize)) as f32) + tadjust);

        (*ray).ttotal = (((*ray).ttotal as f32) + tadjust);
    }

    (*ray).ionce = 0i32;

    (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) =
        ((((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((0i32) as isize)) as f64)
            + ((smst as f64)
                * (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize))
                    as f64)
                    + ((0.25f64 * (sx as f64)) * (smst as f64))))) as f32);

    (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) =
        ((((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((1i32) as isize)) as f64)
            + ((smst as f64)
                * (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize))
                    as f64)
                    + ((0.25f64 * (sy as f64)) * (smst as f64))))) as f32);

    (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((2i32) as isize)) =
        ((((*(std::ptr::addr_of_mut!(xst).cast::<f32>()).offset((2i32) as isize)) as f64)
            + ((smst as f64)
                * (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize))
                    as f64)
                    + ((0.25f64 * (sz as f64)) * (smst as f64))))) as f32);

    dlt = (((smst as f64)
        * ((s as f64)
            + ((smst as f64)
                * ((0.5f64
                    * (((((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>())
                        .offset((0i32) as isize))
                        * sx)
                        + ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>())
                            .offset((1i32) as isize))
                            * sy))
                        + ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>())
                            .offset((2i32) as isize))
                            * sz)) as f64))
                    + ((((((sx * sx) + (sy * sy)) + (sz * sz)) * smst) as f64) / 12f64)))))
        as f32);

    (*ray).t = (((*ray).t as f32) + dlt);

    (*ray).r = ((((*ray).r as f64) + ((dlt as f64) / (s as f64).sqrt())) as f32);

    if (((((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) < xmin)
        || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) > xmax))
        || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) < ymin))
        || ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)) > ymax))
    {
        fprintf!(&*jpfp, "this ray out of model\n");

        (*ray).icode = SPOILED;

        return OUT_OF_MODEL;
    }
    /*****************************************************
    Once the shortest sigma is found, we need to adjust
    the slowness vector (px,py,pz) using

    p_(sg) = p_i0 + 0.5*s_i*sg

    *****************************************************/
    (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) =
        ((((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) as f64)
            + ((0.5f64 * (sx as f64)) * (smst as f64))) as f32);

    (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) =
        ((((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) as f64)
            + ((0.5f64 * (sy as f64)) * (smst as f64))) as f32);

    (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) =
        ((((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) as f64)
            + ((0.5f64 * (sz as f64)) * (smst as f64))) as f32);

    sincid = s_intp(
        point,
        std::ptr::addr_of_mut!((*(tetra).offset(((*ray).itetra) as isize))),
        std::ptr::addr_of_mut!((*ray).x).cast::<f32>(),
    );
    /* the current position */
    s_normalize(std::ptr::addr_of_mut!((*ray).p).cast::<f32>(), sincid);

    (*ray).nseg = (*ray).nseg.wrapping_add(1);

    if ((*ray).nseg > maxnsegs) {
        err!("maxnsegs too small");
    }

    (*(*((*ray).xseg).offset(((*ray).nseg) as isize)).offset((0i32) as isize)) =
        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize));

    (*(*((*ray).xseg).offset(((*ray).nseg) as isize)).offset((1i32) as isize)) =
        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize));

    (*(*((*ray).xseg).offset(((*ray).nseg) as isize)).offset((2i32) as isize)) =
        (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((2i32) as isize));

    (*((*ray).tseg).offset(((*ray).nseg) as isize)) =
        (((*((*ray).tseg).offset(((*ray).nseg) as isize)) as f32) + dlt);

    (*ray).ttotal = (((*ray).ttotal as f32) + (*((*ray).tseg).offset(((*ray).nseg) as isize)));

    if ((sgdt < sgf4min) && (sursg >= smst0)) {
        (*ray).lastfacet = (-1i32);

        return HIT_WF_FIRST;
    }
    /* Now, ray hits surface first */
    if ((((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).itetra).cast::<i32>())
        .offset((0i32) as isize))
        == (-1i32))
        || ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).itetra).cast::<i32>())
            .offset((1i32) as isize))
            == (-1i32)))
        || (sursg < smst0))
    {
        (*ray).lastfacet = (-1i32);

        return HIT_SURFACE_FIRST;
    }

    itetranew = (*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).itetra)
        .cast::<i32>())
    .offset((0i32) as isize));

    if ((*ray).itetra == itetranew) {
        itetranew = (*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).itetra)
            .cast::<i32>())
        .offset((1i32) as isize));
    }

    if ((*ray).itetra == itetranew) {
        err!("ray->itetra==itetranew");
    }

    if (itetranew != (-2i32)) {
        if ((*(tetra).offset((itetranew) as isize)).ireg
            == (*(tetra).offset(((*ray).itetra) as isize)).ireg)
        {
            (*ray).itetra = itetranew;

            (*ray).lastfacet = ifhit;

            return HIT_FACET_FIRST;
        }
    }
    /* to apply Snell's law */
    return RTlaw(
        ray, point, facet, tetra, irefseq, nrefseq, ifhit, sincid, itetranew, jpfp,
    );
}

// area (sutetraray.c:2879).

unsafe fn area(
    mut x0: f32,
    mut y0: f32,
    mut x1: f32,
    mut y1: f32,
    mut x2: f32,
    mut y2: f32,
) -> f32 {
    let mut d1: f32 = 0.0;
    let mut d2: f32 = 0.0;
    let mut d3: f32 = 0.0;

    d1 = (((((x0 - x1) * (x0 - x1)) + ((y0 - y1) * (y0 - y1))) as f64).sqrt() as f32);

    d2 = (((((x0 - x2) * (x0 - x2)) + ((y0 - y2) * (y0 - y2))) as f64).sqrt() as f32);

    d3 = (((((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))) as f64).sqrt() as f32);

    return ((((((((d1 + d2) + d3) * (d2 + d3)) * (d1 + d3)) * (d1 + d2)) as f64).sqrt() / 12f64)
        as f32);
}

// cross_slim_tetra (sutetraray.c:2901).

unsafe fn cross_slim_tetra(
    mut ray: *mut RAY,
    mut point: *mut POINT,
    mut facet: *mut FACET,
    mut tetra: *mut TETRA,
    mut irefseq: *mut i32,
    mut nrefseq: i32,
    mut sincid: f32,
    mut jpfp: *mut CFile,
) {
    let mut ifhit: i32 = 0;
    /* last facet it hits */
    let mut ifs: i32 = 0;
    let mut newit: i32 = 0;
    let mut newifs: i32 = 0;

    let mut fit_err: f32 = 0.0;

    let mut ff: f32 = 0.0;

    let mut ireg0: i32 = 0;

    ireg0 = (*(tetra).offset(((*ray).itetra) as isize)).ireg;

    {
        loop {
            {
                fit_err = (9000i32 as f32);

                newifs = (-1i32);

                {
                    ifs = 0i32;
                    loop {
                        if !(ifs < 4i32) {
                            break;
                        }
                        {
                            ifhit = (*(std::ptr::addr_of_mut!(
                                (*(tetra).offset(((*ray).itetra) as isize)).ifacet
                            )
                            .cast::<i32>())
                            .offset((ifs) as isize));

                            if (ifhit == (*ray).lastfacet) {
                                {
                                    let old = ifs;
                                    ifs = ifs + 1;
                                    old
                                };
                                continue;
                            }

                            ff = facet_fit(
                                point,
                                std::ptr::addr_of_mut!((*ray).x).cast::<f32>(),
                                std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize))),
                            );

                            if (ff < fit_err) {
                                newifs = ifhit;

                                fit_err = ff;
                            }
                        }
                        {
                            let old = ifs;
                            ifs = ifs + 1;
                            old
                        };
                    }
                }

                if (newifs == (-1i32)) {
                    err!("newifs not found");
                }

                (*ray).lastfacet = newifs;

                newit = (*(std::ptr::addr_of_mut!((*(facet).offset((newifs) as isize)).itetra)
                    .cast::<i32>())
                .offset((0i32) as isize));

                if (newit == (*ray).itetra) {
                    newit =
                        (*(std::ptr::addr_of_mut!((*(facet).offset((newifs) as isize)).itetra)
                            .cast::<i32>())
                        .offset((1i32) as isize));
                }

                if (newit == (*ray).itetra) {
                    err!("newit wrong");
                }
                /* need to apply Snell's law here */
                if ((*(tetra).offset((newit) as isize)).ireg != ireg0) {
                    RTlaw(
                        ray, point, facet, tetra, irefseq, nrefseq, newifs, sincid, newit, jpfp,
                    );
                }

                (*ray).itetra = newit;

                if ((*ray).itetra == (-1i32)) {
                    break;
                }

                if (((*(tetra).offset(((*ray).itetra) as isize)).v as f64) > 0f64) {
                    break;
                }
            }
        }
    }
}

// facet_fit (sutetraray.c:3012).

unsafe fn facet_fit(mut point: *mut POINT, mut x: *mut f32, mut facet: *mut FACET) -> f32 {
    let mut areat: f32 = 0.0;
    let mut area0: f32 = 0.0;

    area0 = (*facet).area;

    areat = area3d(
        std::slice::from_raw_parts(x, 3),
        std::slice::from_raw_parts(
            std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*facet).ip).cast::<i32>()).offset((1i32) as isize))
                        as isize
                ))
                .x
            )
            .cast::<f32>(),
            3,
        ),
        std::slice::from_raw_parts(
            std::ptr::addr_of_mut!(
                (*(point).offset(
                    (*(std::ptr::addr_of_mut!((*facet).ip).cast::<i32>()).offset((2i32) as isize))
                        as isize
                ))
                .x
            )
            .cast::<f32>(),
            3,
        ),
    );

    areat = ((areat as f32)
        + area3d(
            std::slice::from_raw_parts(
                std::ptr::addr_of_mut!(
                    (*(point).offset(
                        (*(std::ptr::addr_of_mut!((*facet).ip).cast::<i32>())
                            .offset((0i32) as isize)) as isize
                    ))
                    .x
                )
                .cast::<f32>(),
                3,
            ),
            std::slice::from_raw_parts(x, 3),
            std::slice::from_raw_parts(
                std::ptr::addr_of_mut!(
                    (*(point).offset(
                        (*(std::ptr::addr_of_mut!((*facet).ip).cast::<i32>())
                            .offset((2i32) as isize)) as isize
                    ))
                    .x
                )
                .cast::<f32>(),
                3,
            ),
        ));

    areat = ((areat as f32)
        + area3d(
            std::slice::from_raw_parts(
                std::ptr::addr_of_mut!(
                    (*(point).offset(
                        (*(std::ptr::addr_of_mut!((*facet).ip).cast::<i32>())
                            .offset((0i32) as isize)) as isize
                    ))
                    .x
                )
                .cast::<f32>(),
                3,
            ),
            std::slice::from_raw_parts(
                std::ptr::addr_of_mut!(
                    (*(point).offset(
                        (*(std::ptr::addr_of_mut!((*facet).ip).cast::<i32>())
                            .offset((1i32) as isize)) as isize
                    ))
                    .x
                )
                .cast::<f32>(),
                3,
            ),
            std::slice::from_raw_parts(x, 3),
        ));

    if ((area0 as f64) < 1.0000000000000001E-5f64) {
        return (9000i32 as f32);
    }

    return (areat / area0);
}

// RTlaw (sutetraray.c:3052).

unsafe fn RTlaw(
    mut ray: *mut RAY,
    mut point: *mut POINT,
    mut facet: *mut FACET,
    mut tetra: *mut TETRA,
    mut irefseq: *mut i32,
    mut nrefseq: i32,
    mut ifhit: i32,
    mut sincid: f32,
    mut itetranew: i32,
    mut jpfp: *mut CFile,
) -> i32 {
    let mut b: f32 = 0.0;

    let mut strans: f32 = 0.0;

    let mut crossp: f32 = 0.0;

    let mut cn0: f32 = 0.0;
    let mut cn1: f32 = 0.0;
    let mut cn2: f32 = 0.0;

    let mut iregnew: i32 = 0;

    let mut cosa: f32 = 0.0;
    let mut cosb: f32 = 0.0;
    let mut sinb2: f32 = 0.0;

    let mut cosa2: f32 = 0.0;
    let mut sina2: f32 = 0.0;
    /******************************************************************
    Now adjust the direction to penetrate the horizon

    By the R/T law:

        P1  N
      ___________________
         /    \ |beta       /
        /      \|	  	 /
       /	      X		/
      /	      |\       /
     /	      | \alpha/
     --------------------
               P0

    As shown above, let P0 be the normalized initial slowness vector,
    while P0' be the normalized transmitted slowness vector,

    P0*N = cos(alpha)*sqrt(sincid),

    and

    P1*N = cos(beta)*sqrt(strans).

    By Snell's law,

    sin(alpha)/sin(beta) = vup/vdn;

    since P0' still lies in the same plane with P0 and N, suppose

    P1 = a*P0 + b*N,

    b = cos(beta)*sqrt(strans) - a*cos(alpha)*sqrt(sincid)

    Since b*N = P1 - a*P0, thus

    b = P1*N - a*P0*N

    using sloth, we have P0*P0 = sincid, P1*P1 = strans, thus we have

    P0*P0*a^2 + 2P0*N*a*b + b^2 N*N = strans

    => sincid*a^2 + 2cos(alpha)*sqrt(sincid)*a*b + b^2 = strans

    => a^2 = 1

    since a is positive for transmitted waves,

    a=sqrt(a^2)

    So we get P1.

    For reflection, just do:

    p0-p1=a N

    p1=p0 - a N

    where N points to upward here.

    sincid = sincid + a^2 - 2 a sqrt(sincid) cos(alpha)

    a^2 = 2 a sqrt(sincid) cos(alpha)

    a = 2 sqrt(sincid) cos(alpha)

    p1[0]=p0[0]-a*cn0;
    p1[1]=p0[1]-a*cn1;
    p1[2]=p0[2]-a*cn2;

    The derivative of the curved tile is:
    (2ax+by+cx+n0,bx+dz+2ey+n1,cx+dy+2fz+n2)
    *************************************************/
    cn0 = ((((((2f64
        * ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((0i32) as isize)) as f64))
        * ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)) as f64))
        + (((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((1i32) as isize))
            * (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize)))
            as f64))
        + (((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((2i32) as isize))
            * (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)))
            as f64))
        + ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((6i32) as isize)) as f64)) as f32);

    cn1 = (((((((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
        .offset((1i32) as isize))
        * (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)))
        + ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((3i32) as isize))
            * (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((2i32) as isize))))
        as f64)
        + ((2f64
            * ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
                .offset((4i32) as isize)) as f64))
            * ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize))
                as f64)))
        + ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((7i32) as isize)) as f64)) as f32);

    cn2 = (((((((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
        .offset((2i32) as isize))
        * (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((0i32) as isize)))
        + ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((3i32) as isize))
            * (*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((1i32) as isize))))
        as f64)
        + ((2f64
            * ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
                .offset((5i32) as isize)) as f64))
            * ((*(std::ptr::addr_of_mut!((*ray).x).cast::<f32>()).offset((2i32) as isize))
                as f64)))
        + ((*(std::ptr::addr_of_mut!((*(facet).offset((ifhit) as isize)).ct).cast::<f32>())
            .offset((8i32) as isize)) as f64)) as f32);
    /* make the normal facing downward */
    if ((cn2 as f64) > 0f64) {
        cn0 = (-cn0);

        cn1 = (-cn1);

        cn2 = (-cn2);
    }

    WF_normalize(
        std::ptr::addr_of_mut!(cn0),
        std::ptr::addr_of_mut!(cn1),
        std::ptr::addr_of_mut!(cn2),
    );

    crossp = ((((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize))
        * cn0)
        + ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) * cn1))
        + ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) * cn2));

    cosa2 = ((crossp * crossp) / sincid);

    cosa = ((cosa2 as f64).sqrt() as f32);

    sina2 = ((1f64 - (cosa2 as f64)) as f32);

    sina2 = ((if ((0f64) > ((sina2) as f64)) {
        (0f64)
    } else {
        ((sina2) as f64)
    }) as f32);

    if (itetranew == (-2i32)) {
        iregnew = (-2i32);
    } else {
        iregnew = (*(tetra).offset((itetranew) as isize)).ireg;
    }

    if ((((*ray).iref < nrefseq) && (iregnew == (*(irefseq).offset(((*ray).iref) as isize))))
        || (itetranew == (-2i32)))
    {
        fprintf!(&*jpfp, "To be reflected: iregnew=%d\n", iregnew);

        fprintf!(
            &*jpfp,
            "irefseq[%d]=%d\n",
            (*ray).iref,
            (*(irefseq).offset(((*ray).iref) as isize))
        );

        b = (((2f64 * (sincid as f64).sqrt()) * (cosa as f64)) as f32);

        fprintf!(
            &*jpfp,
            "b=%f p=%f %f %f\n",
            (b as f64),
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) as f64),
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) as f64),
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) as f64)
        );
        /* because n points downward */
        if ((*ray).shootupward == 1i32) {
            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) =
                (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize))
                    as f32)
                    - (b * cn0));

            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) =
                (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize))
                    as f32)
                    - (b * cn1));

            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) =
                (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize))
                    as f32)
                    - (b * cn2));
        } else {
            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) =
                (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize))
                    as f32)
                    + (b * cn0));

            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) =
                (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize))
                    as f32)
                    + (b * cn1));

            (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) =
                (((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize))
                    as f32)
                    + (b * cn2));
        }

        fprintf!(
            &*jpfp,
            "ref: p=%f %f %f\n",
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) as f64),
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) as f64),
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) as f64)
        );

        s_normalize(std::ptr::addr_of_mut!((*ray).p).cast::<f32>(), sincid);

        (*ray).iref = (*ray).iref.wrapping_add(1);

        (*ray).shootupward = (-(*ray).shootupward);
    } else {
        (*ray).itetra = itetranew;

        strans = s_intp(
            point,
            std::ptr::addr_of_mut!((*(tetra).offset((itetranew) as isize))),
            std::ptr::addr_of_mut!((*ray).x).cast::<f32>(),
        );
        /* the current position */
        (*ray).v = ((1f64 / (strans as f64).sqrt()) as f32);

        sinb2 = ((sina2 * sincid) / strans);
        /* that's sloth, you silly! */
        /******************************************************
        if sinb2 > 1, it is overcritical, return SPOILED
        *******************************************************/
        if ((sinb2 as f64) >= 1f64) {
            (*ray).icode = SPOILED;

            fprintf!(&*jpfp, "Overcritical\n");

            (*ray).lastfacet = (-1i32);

            (*ray).icode = SPOILED;

            return OUT_OF_MODEL;
        }

        cosb = ((1f64 - (sinb2 as f64)).sqrt() as f32);

        b = ((((cosb as f64) * (strans as f64).sqrt()) - ((cosa as f64) * (sincid as f64).sqrt()))
            as f32);

        (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize)) =
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((0i32) as isize))
                + (b * cn0));

        (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize)) =
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((1i32) as isize))
                + (b * cn1));

        (*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize)) =
            ((*(std::ptr::addr_of_mut!((*ray).p).cast::<f32>()).offset((2i32) as isize))
                + (b * cn2));

        s_normalize(std::ptr::addr_of_mut!((*ray).p).cast::<f32>(), strans);
    }

    (*ray).lastfacet = ifhit;

    return HIT_FACET_FIRST;
}

// s_normalize (sutetraray.c:3303).

unsafe fn s_normalize(mut p: *mut f32, mut s: f32) {
    let mut norm: f32 = 0.0;

    norm = (((((((*(p).offset((0i32) as isize)) * (*(p).offset((0i32) as isize)))
        + ((*(p).offset((1i32) as isize)) * (*(p).offset((1i32) as isize))))
        + ((*(p).offset((2i32) as isize)) * (*(p).offset((2i32) as isize))))
        / s) as f64)
        .sqrt() as f32);

    (*(p).offset((0i32) as isize)) = (((*(p).offset((0i32) as isize)) as f32) / norm);

    (*(p).offset((1i32) as isize)) = (((*(p).offset((1i32) as isize)) as f32) / norm);

    (*(p).offset((2i32) as isize)) = (((*(p).offset((2i32) as isize)) as f32) / norm);
}

// tetrasolver (sutetraray.c:3327).

unsafe fn tetrasolver(
    mut ct: *mut f32,
    mut p: *mut f32,
    mut gs: *mut f32,
    mut xst: *mut f32,
    mut sg0: f32,
) -> f32 {
    let mut sg: f32 = 0.0;
    let mut dfdsg: f32 = 0.0;
    let mut fsg: f32 = 0.0;
    let mut x0: f32 = 0.0;
    let mut y0: f32 = 0.0;
    let mut z0: f32 = 0.0;

    let mut iter: i32 = 0;

    let mut MaxIterations: i32 = 6i32;

    sg = sg0;

    {
        iter = 0i32;
        loop {
            if !(iter < MaxIterations) {
                break;
            }
            {
                x0 = ((((*(xst).offset((0i32) as isize)) as f64)
                    + ((((*(p).offset((0i32) as isize)) as f64)
                        + ((0.25f64 * ((*(gs).offset((0i32) as isize)) as f64)) * (sg as f64)))
                        * (sg as f64))) as f32);

                y0 = ((((*(xst).offset((1i32) as isize)) as f64)
                    + ((((*(p).offset((1i32) as isize)) as f64)
                        + ((0.25f64 * ((*(gs).offset((1i32) as isize)) as f64)) * (sg as f64)))
                        * (sg as f64))) as f32);

                z0 = ((((*(xst).offset((2i32) as isize)) as f64)
                    + ((((*(p).offset((2i32) as isize)) as f64)
                        + ((0.25f64 * ((*(gs).offset((2i32) as isize)) as f64)) * (sg as f64)))
                        * (sg as f64))) as f32);

                fsg = ((((((((((((*(ct).offset((0i32) as isize)) * x0) * x0)
                    + (((*(ct).offset((1i32) as isize)) * x0) * y0))
                    + (((*(ct).offset((2i32) as isize)) * x0) * z0))
                    + (((*(ct).offset((3i32) as isize)) * y0) * z0))
                    + (((*(ct).offset((4i32) as isize)) * y0) * y0))
                    + (((*(ct).offset((5i32) as isize)) * z0) * z0))
                    + ((*(ct).offset((6i32) as isize)) * x0))
                    + ((*(ct).offset((7i32) as isize)) * y0))
                    + ((*(ct).offset((8i32) as isize)) * z0))
                    + (*(ct).offset((9i32) as isize)));

                if ((fsg as f64).abs() < 1.0000000000000001E-5f64) {
                    break;
                }

                dfdsg = ((((((*(p).offset((0i32) as isize)) as f64)
                    * (((((2f64 * ((*(ct).offset((0i32) as isize)) as f64)) * (x0 as f64))
                        + (((*(ct).offset((1i32) as isize)) * y0) as f64))
                        + (((*(ct).offset((2i32) as isize)) * x0) as f64))
                        + ((*(ct).offset((6i32) as isize)) as f64)))
                    + (((*(p).offset((1i32) as isize)) as f64)
                        * ((((((*(ct).offset((1i32) as isize)) * y0)
                            + ((*(ct).offset((3i32) as isize)) * z0))
                            as f64)
                            + ((2f64 * ((*(ct).offset((4i32) as isize)) as f64)) * (y0 as f64)))
                            + ((*(ct).offset((7i32) as isize)) as f64))))
                    + (((*(p).offset((2i32) as isize)) as f64)
                        * ((((((*(ct).offset((2i32) as isize)) * x0)
                            + ((*(ct).offset((3i32) as isize)) * y0))
                            as f64)
                            + ((2f64 * ((*(ct).offset((4i32) as isize)) as f64)) * (z0 as f64)))
                            + ((*(ct).offset((8i32) as isize)) as f64))))
                    as f32);

                if ((dfdsg as f64).abs() < 1.0000000000000001E-5f64) {
                    break;
                }

                sg = ((sg as f32) - (fsg / dfdsg));
            }
            {
                let old = iter;
                iter = iter + 1;
                old
            };
        }
    }

    if ((sg as f64) <= 0f64) {
        return sg0;
    }

    return sg;
}

// Program storage keeps the original linked triangle and two dimensional
// array layout. All allocations are owned by the state in main.
struct Allocation {
    ptr: *mut u8,
    layout: std::alloc::Layout,
}
impl Drop for Allocation {
    fn drop(&mut self) {
        unsafe {
            std::alloc::dealloc(self.ptr, self.layout);
        }
    }
}
struct State {
    tro: Trace,
    input: CFile,
    output: CFile,
    errors: CFile,
    storage: Vec<Allocation>,
    files: Vec<Box<CFile>>,
    strings: Vec<std::ffi::CString>,
}
unsafe fn allocate(state: *mut State, bytes: usize) -> *mut std::ffi::c_void {
    // Boundary facets have tetra indices -1 and -2. The C tracer reads
    // their region fields before checking for a boundary. Reserve zeroed
    // guard storage so those otherwise undefined reads remain addressable.
    // Plain malloc callers in the C do not check failure. Return NULL for
    // impossible sizes too, rather than panicking on Rust's LayoutError.
    let Some(size) = bytes.checked_add(512) else {
        return std::ptr::null_mut();
    };
    let Ok(layout) = std::alloc::Layout::from_size_align(size, 16) else {
        return std::ptr::null_mut();
    };
    let ptr = std::alloc::alloc_zeroed(layout);
    if ptr.is_null() {
        return std::ptr::null_mut();
    }
    (*state).storage.push(Allocation { ptr, layout });
    ptr.add(256).cast()
}
unsafe fn array2(state: *mut State, n1: usize, n2: usize) -> *mut *mut f32 {
    let sizes = n2.checked_add(1).and_then(|nrows| {
        Some((
            nrows.checked_mul(std::mem::size_of::<*mut f32>())?,
            nrows.checked_mul(n1)?.checked_mul(std::mem::size_of::<f32>())?,
        ))
    });
    let Some((row_bytes, data_bytes)) = sizes else {
        err!("ealloc.c: malloc failed (%s)", oxwave_core::err::strerror(12));
    };
    let rows = allocate(state, row_bytes).cast::<*mut f32>();
    if rows.is_null() {
        err!("ealloc.c: malloc failed (%s)", oxwave_core::err::strerror(12));
    }
    let data = allocate(state, data_bytes).cast::<f32>();
    if data.is_null() {
        err!("ealloc.c: malloc failed (%s)", oxwave_core::err::strerror(12));
    }
    // The C test is nseg>maxnsegs, so nseg==maxnsegs accesses a row
    // outside the allocation. Reserve a usable extra row so the tracer
    // can continue to its own limit check. C alloc2 also sets p[0] when
    // n2==0; that row is needed to initialize the source position.
    for i in 0..=n2 {
        *rows.add(i) = data.add(i * n1);
    }
    rows
}
unsafe fn retire_triangle(_ptr: *mut std::ffi::c_void) {
    // The C loop frees itrip and then reads itrip->next. In the reference
    // macOS allocator this traversal faults when a consecutive spoiled
    // triangle is removed. Preserve the signal without a Rust use-after-free.
    extern "C" {
        fn signal(signal: i32, handler: usize) -> usize;
        fn raise(signal: i32) -> i32;
    }
    signal(11, 0);
    raise(11);
    std::process::abort();
}
unsafe fn cstr<'a>(p: *mut i8) -> &'a str {
    std::ffi::CStr::from_ptr(p).to_str().unwrap_or("")
}
unsafe fn open_file(state: *mut State, path: *mut i8, mode: *mut i8) -> *mut CFile {
    if let Some(f) = fopen(cstr(path), cstr(mode)) {
        let mut f = Box::new(f);
        let p = &mut *f as *mut CFile;
        (*state).files.push(f);
        p
    } else {
        std::ptr::null_mut()
    }
}
unsafe fn par_string(state: *mut State, name: &str, dest: *mut *mut i8) -> i32 {
    let mut s = String::new();
    let count = getparstring(name, &mut s);
    if count != 0 {
        let s = std::ffi::CString::new(s).unwrap();
        *dest = s.as_ptr() as *mut i8;
        (*state).strings.push(s);
    }
    count
}
unsafe fn par_int(name: &str, dest: *mut i32) -> i32 {
    if name == "irefseq" {
        let mut values = Vec::<i32>::new();
        let count = getparint(name, &mut values);
        for (i, v) in values.into_iter().enumerate() {
            *dest.add(i) = v;
        }
        count
    } else {
        getparint(name, &mut *dest)
    }
}
// Match only the literal portions of this program's model-file formats;
// numeric conversions use the frozen core fscanf functions.
fn scan_literal(fp: &CFile, text: &str) -> bool {
    for b in text.bytes() {
        if b.is_ascii_whitespace() {
            let mut c = fgetc(fp);
            while c != EOF && (c as u8).is_ascii_whitespace() {
                c = fgetc(fp);
            }
            if c != EOF {
                ungetc(c, fp);
            }
        } else {
            let c = fgetc(fp);
            if c != b as i32 {
                if c != EOF {
                    ungetc(c, fp);
                }
                return false;
            }
        }
    }
    true
}
fn main() {
    let argv = args();
    initargs(&argv);
    unsafe {
        let mut state = State {
            tro: Trace::new(),
            input: stdin(),
            output: stdout(),
            errors: stderr(),
            storage: Vec::new(),
            files: Vec::new(),
            strings: Vec::new(),
        };
        let status = run(&mut state);
        exit(status);
    }
}
