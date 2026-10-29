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

//! Dynamic ray tracing helpers, translated from tri/main/triray.c (2011/11/21).
//! Credits: Dave Hale (1991); modifications by Andreas Rueger (1993).

#![allow(
    non_snake_case,
    unused_parens,
    unused_assignments,
    unused_variables,
    unused_mut,
    unused_must_use,
    dead_code
)]
use super::*;

// shootRays (C line 238): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn shootRays(
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
    unsafe {
        /*****************************************************************************
        Shoot rays via dynamic ray tracing for a sloth model
        ******************************************************************************
        Input:
        m		trianglulated sloth model
        xs		horizontal coordinate of source
        zs		vertical coordinate (depth) of source
        nangle		number of ray takeoff angles
        dangle		increment in ray takeoff angle
        fangle		first ray takeoff angle
        rc		pointer to reflector structure
        kk		number of reflectors

        Output:
        rsp		array[nangle] of pointers to linked list of RaySteps
        re		array[nangle] of RayEnds
        ******************************************************************************
        Notes:
        Rays are traced from triangle to triangle.
        Whenever a ray crosses a triangle edge, two RaySteps are added
        to the linked list of RaySteps, one for each side of the edge.
        The RaySteps are useful if ray parameters along the entire raypath
        are required.  The RayEnds are useful if parameters at the ray
        endpoint are required.
        ******************************************************************************
        Author:  Dave Hale, Colorado School of Mines, 02/28/91
        Modified: Andreas Rueger, Colorado School of Mines, 08/12/93
        ******************************************************************************/

        let mut iangle: i32 = std::mem::zeroed();
        let mut refl: i32 = std::mem::zeroed();
        let mut stop: i32 = std::mem::zeroed();
        let mut i: i32 = std::mem::zeroed();
        let mut temp: i32 = std::mem::zeroed();
        let mut test: i32 = std::mem::zeroed();
        let mut angle: f32 = std::mem::zeroed();
        let mut sb: f32 = std::mem::zeroed();
        let mut zmax: f32 = std::mem::zeroed();
        let mut xmax: f32 = std::mem::zeroed();
        let mut eps: f32 = std::mem::zeroed();
        let mut tris: *mut Face = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
        let mut ea: *mut EdgeAttributes = std::mem::zeroed();
        let mut rs: *mut RayStep = std::mem::zeroed();
        let mut rslast: *mut RayStep = std::mem::zeroed();
        /* determine triangle in model containing source location */
        if m.is_null() {
            c_segv();
        }
        tris = insideTriInModel(m, std::ptr::null_mut::<Face>(), zs, xs);
        /* if source sits on edge, move it by eps */
        if ({
            test = ((checkIfSourceOnEdge(tris, zs, xs) != 0_i32) as i32);
            test
        }) != 0
        {
            zmax = (*(m)).xmax;
            eps = (*(m)).eps;
            if (test == 2_i32) {
                eps = (((5_i32) as f32) * eps);
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
            tris = insideTriInModel(m, std::ptr::null_mut::<Face>(), zs, xs);
        }
        /* determine sloth at source location (beginning of ray) */
        fa = (((*(tris)).fa) as *mut FaceAttributes);
        sb = (((*(fa)).s00 + ((*(fa)).dsdx * xs)) + ((*(fa)).dsdz * zs));
        /* loop over takeoff angles */
        {
            iangle = 0_i32;
            angle = fangle;
            while (iangle < nangle) {
                {
                    /* initialize linked list of ray steps */
                    rs = {
                        (*(rsp).offset((iangle) as isize)) =
                            ((Box::into_raw(Box::new(std::mem::zeroed::<RayStep>()))
                                .cast::<c_void>()) as *mut RayStep);
                        (*(rsp).offset((iangle) as isize))
                    };
                    (*(rs)).sigma = ((0.0_f64) as f32);
                    (*(rs)).x = xs;
                    (*(rs)).z = zs;
                    (*(rs)).px = ((((angle) as f64).sin() * ((sb) as f64).sqrt()) as f32);
                    (*(rs)).pz = ((((angle) as f64).cos() * ((sb) as f64).sqrt()) as f32);
                    (*(rs)).q1 = ((1.0_f64) as f32);
                    (*(rs)).p1 = ((0.0_f64) as f32);
                    (*(rs)).q2 = ((0.0_f64) as f32);
                    (*(rs)).p2 = ((1.0_f64) as f32);
                    (*(rs)).t = ((0.0_f64) as f32);
                    (*(rs)).kmah = 0_i32;
                    (*(rs)).nref = 0_i32;
                    (*(rs)).eu = std::ptr::null_mut::<EdgeUse>();
                    (*(rs)).f = tris;
                    (*(rs)).ampli = ((1_i32) as f32);
                    (*(rs)).ampliphase = ((0_i32) as f32);
                    (*(rs)).atten = ((0_i32) as f32);
                    (*(rs)).rsNext = std::ptr::null_mut::<RayStep>();
                    if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                        fprintf!(&*(infofp), "\n****************************************\nRAY WITH TAKEOFF ANGLE: %g (degrees)\n****************************************\n", ((((angle) as f64) * 180.0_f64) / (3.1415926535897931_f64)));
                    }
                    /* trace ray */
                    loop {
                        {
                            /* trace ray in triangle */
                            rs = traceRayInTri(rs);
                            /* remember last ray step */
                            rslast = rs;
                            /* edge attributes */
                            ea = (((*((*((*(rs)).eu)).e)).ea) as *mut EdgeAttributes);
                            if (ea == ((std::ptr::null_mut::<c_void>()) as *mut EdgeAttributes)) {
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
                                        {
                                            if (temp == (*(rc).offset((i) as isize)).k) {
                                                if ((*((*(rc).offset((i) as isize)).hitseq)
                                                    .offset((0_i32) as isize))
                                                    == 1_i32)
                                                {
                                                    refl = 1_i32;
                                                }
                                                if ((*((*(rc).offset((i) as isize)).hitseq)
                                                    .offset((0_i32) as isize))
                                                    == (-(1_i32)))
                                                {
                                                    stop = 1_i32;
                                                }
                                                ((*(rc).offset((i) as isize)).nhits) =
                                                    ((*(rc).offset((i) as isize)).nhits) + 1;
                                                ((*(rc).offset((i) as isize)).hitseq) =
                                                    ((*(rc).offset((i) as isize)).hitseq).offset(1);
                                                break;
                                            }
                                        }
                                        i = i + 1;
                                    }
                                }
                            }
                            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                                fprintf!(
                                    &*(infofp),
                                    "^^^\nInteraction with interface %d at (x=%g,z=%g).\n",
                                    temp,
                                    (((*(rslast)).x) as f64),
                                    (((*(rslast)).z) as f64)
                                );
                            }
                            /* if ray at stopping edge, stop */
                            if (stop) != 0 {
                                if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                                    fprintf!(&*(infofp), "Hit stopping edge.  Ray stopped.\n");
                                }
                                break;
                            } else {
                                if (refl) != 0 {
                                    rs = reflectRayFromEdge(rs, infofp);
                                } else {
                                    rs = traceRayAcrossEdge(rs, infofp);
                                }
                            }
                        }
                        if !(rs != ((std::ptr::null_mut::<c_void>()) as *mut RayStep)) {
                            break;
                        }
                    }
                    /* reinitialize RayChecks for new ray */
                    {
                        i = 0_i32;
                        while (i < kk) {
                            {
                                (*(rc).offset((i) as isize)).hitseq =
                                    ((*(rc).offset((i) as isize)).hitseq)
                                        .offset(((*(rc).offset((i) as isize)).nhits) as isize * -1);
                                (*(rc).offset((i) as isize)).nhits = 0_i32;
                            }
                            i = i + 1;
                        }
                    }
                    /* save ray end parameters */
                    fa = (((*((*(rslast)).f)).fa) as *mut FaceAttributes);
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
                    (*(re).offset((iangle) as isize)).se = (((*(fa)).s00
                        + ((*(fa)).dsdx * (*(rslast)).x))
                        + ((*(fa)).dsdz * (*(rslast)).z));
                    (*(re).offset((iangle) as isize)).dsdxe = (*(fa)).dsdx;
                    (*(re).offset((iangle) as isize)).dsdze = (*(fa)).dsdz;
                    (*(re).offset((iangle) as isize)).ab = angle;
                    (*(re).offset((iangle) as isize)).kend = temp;
                    (*(re).offset((iangle) as isize)).ampli = (*(rslast)).ampli;
                    (*(re).offset((iangle) as isize)).ampliphase = (*(rslast)).ampliphase;
                    (*(re).offset((iangle) as isize)).atten = (*(rslast)).atten;
                    (*(re).offset((iangle) as isize)).dangle = dangle;
                    /* optional output of rayend information */
                    if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                        fprintf!(&*(infofp), "---------------------------------------------------\n\nThe following information is stored at the rayend:\n\n");
                        fprintf!(
                            &*(infofp),
                            "Ray stops at (%g,%g) at interface %d\n",
                            (((*(re).offset((iangle) as isize)).x) as f64),
                            (((*(re).offset((iangle) as isize)).z) as f64),
                            (*(re).offset((iangle) as isize)).kend
                        );
                        fprintf!(
                            &*(infofp),
                            "Takeoff angle=%g   Emergence angle=%g\n",
                            (((((*(re).offset((iangle) as isize)).ab) as f64) * 180.0_f64)
                                / (3.1415926535897931_f64)),
                            (((-((((((*(re).offset((iangle) as isize)).px) as f64) * 1.0_f64)
                                / (((*(re).offset((iangle) as isize)).se) as f64).sqrt())
                            .asin()))
                                * 180.0_f64)
                                / (3.1415926535897931_f64))
                        );
                        fprintf!(
                            &*(infofp),
                            "Takeoff velocity=%g   Emergence velocity=%g\n",
                            (1.0_f64 / ((sb) as f64).sqrt()),
                            (1.0_f64 / (((*(re).offset((iangle) as isize)).se) as f64).sqrt())
                        );
                        fprintf!(
                            &*(infofp),
                            "Slowness components px=%g   pz=%g\n",
                            (((*(re).offset((iangle) as isize)).px) as f64),
                            (((*(re).offset((iangle) as isize)).pz) as f64)
                        );
                        fprintf!(
                            &*(infofp),
                            "Traveltime=%g   Sigma=%g\n",
                            (((*(re).offset((iangle) as isize)).t) as f64),
                            (((*(re).offset((iangle) as isize)).sigma) as f64)
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
                            (((*(re).offset((iangle) as isize)).q1) as f64),
                            (((*(re).offset((iangle) as isize)).q2) as f64)
                        );
                        fprintf!(
                            &*(infofp),
                            "Ray propagator p1=%g   p2=%g\n",
                            (((*(re).offset((iangle) as isize)).p1) as f64),
                            (((*(re).offset((iangle) as isize)).p2) as f64)
                        );
                        fprintf!(
                            &*(infofp),
                            "Attenuation factor=%g   Angle increment=%g\n",
                            (((*(re).offset((iangle) as isize)).atten) as f64),
                            ((((dangle) as f64) * 180.0_f64) / (3.1415926535897931_f64))
                        );
                        fprintf!(
                            &*(infofp),
                            "Refl/Transm Amplitude=%g   Phase=%g\n",
                            (((*(re).offset((iangle) as isize)).ampli) as f64),
                            (((*(re).offset((iangle) as isize)).ampliphase) as f64)
                        );
                        fprintf!(
                            &*(infofp),
                            "---------------------------------------------------\n\n"
                        );
                    }
                }
                iangle = iangle + 1;
                angle = ((((angle) as f32) + dangle) as f32);
            }
        }
    }
}

// writeRays (C line 460): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn writeRays(
    mut m: *mut Model,
    mut rayfp: *mut CFile,
    mut nxz: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
    mut re: *mut RayEnd,
    mut krecord: i32,
    mut prim: i32,
    mut infofp: *mut CFile,
    mut outparfp: *mut CFile,
) {
    unsafe {
        /* for each ray, write x,z pairs */

        let mut iray: i32 = std::mem::zeroed();
        let mut nrs: i32 = std::mem::zeroed();
        let mut irs: i32 = std::mem::zeroed();
        let mut nxzw: i32 = std::mem::zeroed();
        let mut ns: i32 = std::mem::zeroed();
        let mut is: i32 = std::mem::zeroed();
        let mut icount: i32 = std::mem::zeroed();
        let mut sigmamax: f32 = ((0.0_f64) as f32);
        let mut sigma1: f32 = std::mem::zeroed();
        let mut sigma2: f32 = std::mem::zeroed();
        let mut x: f32 = std::mem::zeroed();
        let mut z: f32 = std::mem::zeroed();
        let mut px: f32 = std::mem::zeroed();
        let mut pz: f32 = std::mem::zeroed();
        let mut dsdx: f32 = std::mem::zeroed();
        let mut dsdz: f32 = std::mem::zeroed();
        let mut ds: f32 = std::mem::zeroed();
        let mut xs: f32 = std::mem::zeroed();
        let mut zs: f32 = std::mem::zeroed();
        let mut dsigma: f32 = std::mem::zeroed();
        let mut rsi: *mut RayStep = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
        m = (m).offset((0_i32) as isize);
        /* dummy */
        /* initialize counter */
        icount = 0_i32;
        /* loop over rays */
        {
            iray = 0_i32;
            while (iray < nray) {
                {
                    if ((krecord == 2147483647_i32)
                        || ((krecord == (*(re).offset((iray) as isize)).kend)
                            && ((prim == 2147483647_i32)
                                || (prim == (*(re).offset((iray) as isize)).nref))))
                    {
                        /* count rays */
                        icount = icount + 1;
                        /* count ray steps while determining maximum sigma */
                        {
                            nrs = 0_i32;
                            rsi = (*(rs).offset((iray) as isize));
                            while (rsi != ((std::ptr::null_mut::<c_void>()) as *mut RayStep)) {
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
                                {
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
                                    fa = (((*((*(rsi)).f)).fa) as *mut FaceAttributes);
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
                                    ns = ((((1_i32) as f64)
                                        + ((((sigma2 - sigma1) * ((nxz - 1_i32) as f32)) as f64)
                                            / (2.0_f64 * ((sigmamax) as f64))))
                                        as i32);
                                    /* increment in sigma */
                                    ds = ((sigma2 - sigma1) / ((ns) as f32));
                                    /* loop over sigma */
                                    {
                                        is = 0_i32;
                                        dsigma = ((0.0_f64) as f32);
                                        while (is < ns) {
                                            {
                                                /* compute and write x and z */
                                                xs = ((((x) as f64)
                                                    + (((dsigma) as f64)
                                                        * (((px) as f64)
                                                            + ((((dsigma) as f64) * 0.25_f64)
                                                                * ((dsdx) as f64)))))
                                                    as f32);
                                                zs = ((((z) as f64)
                                                    + (((dsigma) as f64)
                                                        * (((pz) as f64)
                                                            + ((((dsigma) as f64) * 0.25_f64)
                                                                * ((dsdz) as f64)))))
                                                    as f32);
                                                if (efwrite(
                                                    std::slice::from_raw_parts(
                                                        ((&raw mut zs) as *mut c_void)
                                                            as *const f32,
                                                        ((1_i32) as usize),
                                                    ),
                                                    &*(rayfp),
                                                ) != ((1_i32) as usize))
                                                {
                                                    err!("error writing rays!\n");
                                                }
                                                if (efwrite(
                                                    std::slice::from_raw_parts(
                                                        ((&raw mut xs) as *mut c_void)
                                                            as *const f32,
                                                        ((1_i32) as usize),
                                                    ),
                                                    &*(rayfp),
                                                ) != ((1_i32) as usize))
                                                {
                                                    err!("error writing rays!\n");
                                                }
                                                /* if sufficient number written, */
                                                /* break */
                                                if ({
                                                    nxzw = nxzw + 1;
                                                    nxzw
                                                } >= nxz)
                                                {
                                                    break;
                                                }
                                            }
                                            is = is + 1;
                                            dsigma = ((((dsigma) as f32) + ds) as f32);
                                        }
                                    }
                                }
                                irs = irs + 1;
                                rsi = (*(rsi)).rsNext;
                            }
                        }
                        /* if necessary, repeat last (x,z) */
                        while (nxzw < nxz) {
                            if (efwrite(
                                std::slice::from_raw_parts(
                                    ((&raw mut zs) as *mut c_void) as *const f32,
                                    ((1_i32) as usize),
                                ),
                                &*(rayfp),
                            ) != ((1_i32) as usize))
                            {
                                err!("error writing rays!\n");
                            }
                            if (efwrite(
                                std::slice::from_raw_parts(
                                    ((&raw mut xs) as *mut c_void) as *const f32,
                                    ((1_i32) as usize),
                                ),
                                &*(rayfp),
                            ) != ((1_i32) as usize))
                            {
                                err!("error writing rays!\n");
                            }
                            nxzw = nxzw + 1;
                        }
                    } else {
                        /* manipulate rayend information */
                        (*(re).offset((iray) as isize)).kend = (-(1_i32));
                    }
                }
                iray = iray + 1;
            }
        }
        /* end loop over rays */
        if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
            fprintf!(&*(infofp), "\nWrote %d rays to ray file.\n", icount);
        }
        /* write number of rays to parameter file */
        fprintf!(&*(outparfp), "%d\n", icount);
    }
}

// writeWaves (C line 565): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn writeWaves(
    mut wavefp: *mut CFile,
    mut nt: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
) {
    unsafe {
        /* for each ray, write nt z(t),x(t) pairs uniformly sampled in time t */

        let mut iray: i32 = std::mem::zeroed();
        let mut it: i32 = std::mem::zeroed();
        let mut tmax: f32 = std::mem::zeroed();
        let mut dt: f32 = std::mem::zeroed();
        let mut t1: f32 = std::mem::zeroed();
        let mut t2: f32 = std::mem::zeroed();
        let mut x1: f32 = std::mem::zeroed();
        let mut z1: f32 = std::mem::zeroed();
        let mut px1: f32 = std::mem::zeroed();
        let mut pz1: f32 = std::mem::zeroed();
        let mut sigma1: f32 = std::mem::zeroed();
        let mut sigma2: f32 = std::mem::zeroed();
        let mut s00: f32 = std::mem::zeroed();
        let mut dsdx: f32 = std::mem::zeroed();
        let mut dsdz: f32 = std::mem::zeroed();
        let mut ax: f32 = std::mem::zeroed();
        let mut bx: f32 = std::mem::zeroed();
        let mut az: f32 = std::mem::zeroed();
        let mut bz: f32 = std::mem::zeroed();
        let mut at: f32 = std::mem::zeroed();
        let mut bt: f32 = std::mem::zeroed();
        let mut ct: f32 = std::mem::zeroed();
        let mut ti: f32 = std::mem::zeroed();
        let mut sigmah: f32 = std::mem::zeroed();
        let mut sigmam: f32 = std::mem::zeroed();
        let mut tm: f32 = std::mem::zeroed();
        let mut dsigma: f32 = std::mem::zeroed();
        let mut sigma: f32 = std::mem::zeroed();
        let mut t: f32 = std::mem::zeroed();
        let mut x: f32 = std::mem::zeroed();
        let mut z: f32 = std::mem::zeroed();
        let mut xstart: f32 = std::mem::zeroed();
        let mut zstart: f32 = std::mem::zeroed();
        let mut rsi: *mut RayStep = std::mem::zeroed();
        let mut rs1: *mut RayStep = std::mem::zeroed();
        let mut rs2: *mut RayStep = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
        /* determine maximum time for all rays */
        {
            iray = 0_i32;
            tmax = ((0.0_f64) as f32);
            while (iray < nray) {
                {
                    rsi = (*(rs).offset((iray) as isize));
                    while (rsi != ((std::ptr::null_mut::<c_void>()) as *mut RayStep)) {
                        tmax = (if ((tmax) > ((*(rsi)).t)) {
                            (tmax)
                        } else {
                            ((*(rsi)).t)
                        });
                        rsi = (*(rsi)).rsNext;
                    }
                }
                iray = iray + 1;
            }
        }
        /* determine time sampling interval */
        dt = (tmax / ((if (nt > 1_i32) { (nt - 1_i32) } else { 1_i32 }) as f32));
        fprintf!(
            &*(&mut stderr() as *mut CFile),
            "wavefront time increment = %g\n",
            ((dt) as f64)
        );
        /* loop over rays */
        {
            iray = 0_i32;
            while (iray < nray) {
                {
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
                    fa = (((*((*(rs1)).f)).fa) as *mut FaceAttributes);
                    s00 = (*(fa)).s00;
                    dsdx = (*(fa)).dsdx;
                    dsdz = (*(fa)).dsdz;
                    ax = px1;
                    bx = ((0.25_f64 * ((dsdx) as f64)) as f32);
                    az = pz1;
                    bz = ((0.25_f64 * ((dsdz) as f64)) as f32);
                    at = ((s00 + (dsdx * x1)) + (dsdz * z1));
                    bt = ((0.5_f64 * (((dsdx * px1) + (dsdz * pz1)) as f64)) as f32);
                    ct = ((0.083333299999999999_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64))
                        as f32);
                    sigma = sigma1;
                    t = t1;
                    /* remember start x,z */
                    xstart = x1;
                    zstart = z1;
                    /* loop over times */
                    {
                        it = 0_i32;
                        ti = ((0.0_f64) as f32);
                        while (it < nt) {
                            {
                                /* if necessary, go to next ray step */
                                if (t2 < ti) {
                                    loop {
                                        {
                                            rs1 = rs2;
                                            rs2 = (*(rs1)).rsNext;
                                            if (rs2
                                                == ((std::ptr::null_mut::<c_void>())
                                                    as *mut RayStep))
                                            {
                                                break;
                                            }
                                        }
                                        if !((*(rs2)).t < ti) {
                                            break;
                                        }
                                    }
                                    if (rs2 == ((std::ptr::null_mut::<c_void>()) as *mut RayStep)) {
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
                                    fa = (((*((*(rs1)).f)).fa) as *mut FaceAttributes);
                                    s00 = (*(fa)).s00;
                                    dsdx = (*(fa)).dsdx;
                                    dsdz = (*(fa)).dsdz;
                                    ax = px1;
                                    bx = ((0.25_f64 * ((dsdx) as f64)) as f32);
                                    az = pz1;
                                    bz = ((0.25_f64 * ((dsdz) as f64)) as f32);
                                    at = ((s00 + (dsdx * x1)) + (dsdz * z1));
                                    bt =
                                        ((0.5_f64 * (((dsdx * px1) + (dsdz * pz1)) as f64)) as f32);
                                    ct = ((0.083333299999999999_f64
                                        * (((dsdx * dsdx) + (dsdz * dsdz)) as f64))
                                        as f32);
                                    sigma = sigma1;
                                    t = t1;
                                }
                                /* determine dsigma for this time by bisection */
                                sigmah = sigma2;
                                while ((((t) as f64) + (0.01_f64 * ((dt) as f64))) < ((ti) as f64))
                                {
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
                                if (efwrite(
                                    std::slice::from_raw_parts(
                                        ((&raw mut z) as *mut c_void) as *const f32,
                                        ((1_i32) as usize),
                                    ),
                                    &*(wavefp),
                                ) != ((1_i32) as usize))
                                {
                                    err!("error writing waves!\n");
                                }
                                if (efwrite(
                                    std::slice::from_raw_parts(
                                        ((&raw mut x) as *mut c_void) as *const f32,
                                        ((1_i32) as usize),
                                    ),
                                    &*(wavefp),
                                ) != ((1_i32) as usize))
                                {
                                    err!("error writing waves!\n");
                                }
                            }
                            it = it + 1;
                            ti = ((((ti) as f32) + dt) as f32);
                        }
                    }
                    /* finish writing x and z */
                    while (it < nt) {
                        if (efwrite(
                            std::slice::from_raw_parts(
                                ((&raw mut zstart) as *mut c_void) as *const f32,
                                ((1_i32) as usize),
                            ),
                            &*(wavefp),
                        ) != ((1_i32) as usize))
                        {
                            err!("error writing waves!\n");
                        }
                        if (efwrite(
                            std::slice::from_raw_parts(
                                ((&raw mut xstart) as *mut c_void) as *const f32,
                                ((1_i32) as usize),
                            ),
                            &*(wavefp),
                        ) != ((1_i32) as usize))
                        {
                            err!("error writing waves!\n");
                        }
                        it = it + 1;
                    }
                }
                iray = iray + 1;
            }
        }
    }
}

// traceRayInTri (C line 687): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn traceRayInTri(mut rs: *mut RayStep) -> *mut RayStep {
    unsafe {
        /* Trace ray in triangle.  Return pointer to new RayStep. */

        let mut kmah: i32 = std::mem::zeroed();
        let mut nref: i32 = std::mem::zeroed();
        let mut x: f32 = std::mem::zeroed();
        let mut z: f32 = std::mem::zeroed();
        let mut px: f32 = std::mem::zeroed();
        let mut pz: f32 = std::mem::zeroed();
        let mut t: f32 = std::mem::zeroed();
        let mut q1: f32 = std::mem::zeroed();
        let mut p1: f32 = std::mem::zeroed();
        let mut q2: f32 = std::mem::zeroed();
        let mut p2: f32 = std::mem::zeroed();
        let mut ddt: f32 = std::mem::zeroed();
        let mut s00: f32 = std::mem::zeroed();
        let mut dsdx: f32 = std::mem::zeroed();
        let mut dsdz: f32 = std::mem::zeroed();
        let mut xa: f32 = std::mem::zeroed();
        let mut za: f32 = std::mem::zeroed();
        let mut xb: f32 = std::mem::zeroed();
        let mut zb: f32 = std::mem::zeroed();
        let mut dx: f32 = std::mem::zeroed();
        let mut dz: f32 = std::mem::zeroed();
        let mut a: f32 = std::mem::zeroed();
        let mut b: f32 = std::mem::zeroed();
        let mut c: f32 = std::mem::zeroed();
        let mut ds: f32 = std::mem::zeroed();
        let mut q: f32 = std::mem::zeroed();
        let mut sigma: f32 = std::mem::zeroed();
        let mut dsigma: f32 = std::mem::zeroed();
        let mut dsigma1: f32 = std::mem::zeroed();
        let mut dsigma2: f32 = std::mem::zeroed();
        let mut small: f32 = std::mem::zeroed();
        let mut frac: f32 = std::mem::zeroed();
        let mut ampli: f32 = std::mem::zeroed();
        let mut ampliphase: f32 = std::mem::zeroed();
        let mut atten: f32 = std::mem::zeroed();
        let mut qfac: f32 = std::mem::zeroed();
        let mut eu: *mut EdgeUse = std::mem::zeroed();
        let mut eut: *mut EdgeUse = std::mem::zeroed();
        let mut eusmall: *mut EdgeUse = std::mem::zeroed();
        let mut euzero: *mut EdgeUse = std::ptr::null_mut::<EdgeUse>();
        let mut f: *mut Face = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
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
        fa = (((*(f)).fa) as *mut FaceAttributes);
        s00 = (*(fa)).s00;
        dsdx = (*(fa)).dsdx;
        dsdz = (*(fa)).dsdz;
        /* determine Q-factor */
        qfac = (*(fa)).qfac;
        /* determine edge intersection with smallest positive dsigma */
        eusmall = std::ptr::null_mut::<EdgeUse>();
        small = 3.40282347E+38_f32;
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
                if ((((dsdx) as f64) == 0.0_f64) && (((dsdz) as f64) == 0.0_f64)) {
                    if (((b) as f64) != 0.0_f64) {
                        dsigma1 = ((-(c)) / b);
                    } else {
                        dsigma1 = 3.40282347E+38_f32;
                    }
                    dsigma2 = 3.40282347E+38_f32;
                } else {
                    a = ((0.25_f64 * (((dx * dsdz) - (dz * dsdx)) as f64)) as f32);
                    ds = ((((b * b) as f64) - ((4.0_f64 * ((a) as f64)) * ((c) as f64))) as f32);
                    if (((ds) as f64) < 0.0_f64) {
                        dsigma1 = {
                            dsigma2 = 3.40282347E+38_f32;
                            dsigma2
                        };
                    } else {
                        q = (((-(0.5_f64))
                            * (if (((b) as f64) < 0.0_f64) {
                                (((b) as f64) - ((ds) as f64).sqrt())
                            } else {
                                (((b) as f64) + ((ds) as f64).sqrt())
                            })) as f32);
                        if (((a) as f64) != 0.0_f64) {
                            dsigma1 = (q / a);
                        } else {
                            dsigma1 = 3.40282347E+38_f32;
                        }
                        if (((q) as f64) != 0.0_f64) {
                            dsigma2 = (c / q);
                        } else {
                            dsigma2 = 3.40282347E+38_f32;
                        }
                    }
                }
                /* remember edge with smallest positive dsigma */
                if ((0.0_f64 < ((dsigma1) as f64)) && (dsigma1 < small)) {
                    small = dsigma1;
                    eusmall = eut;
                }
                if ((0.0_f64 < ((dsigma2) as f64)) && (dsigma2 < small)) {
                    small = dsigma2;
                    eusmall = eut;
                }
                /* remember edge for which dsigma is zero */
                if (((dsigma1) as f64) == 0.0_f64) {
                    euzero = eut;
                }
                if (((dsigma2) as f64) == 0.0_f64) {
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
        if (eusmall != ((std::ptr::null_mut::<c_void>()) as *mut EdgeUse)) {
            dsigma = small;
            eu = eusmall;
        } else {
            dsigma = ((0.0_f64) as f32);
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
        sigma = ((((sigma) as f32) + dsigma) as f32);
        ddt = ((((dsigma) as f64)
            * ((((s00 + (dsdx * x)) + (dsdz * z)) as f64)
                + (((dsigma) as f64)
                    * ((0.5_f64 * (((dsdx * px) + (dsdz * pz)) as f64))
                        + (((dsigma) as f64)
                            * (0.083333299999999999_f64
                                * (((dsdx * dsdx) + (dsdz * dsdz)) as f64)))))))
            as f32);
        t = ((((t) as f32) + ddt) as f32);
        x = ((((x) as f64)
            + (((dsigma) as f64)
                * (((px) as f64) + ((0.25_f64 * ((dsdx) as f64)) * ((dsigma) as f64)))))
            as f32);
        z = ((((z) as f64)
            + (((dsigma) as f64)
                * (((pz) as f64) + ((0.25_f64 * ((dsdz) as f64)) * ((dsigma) as f64)))))
            as f32);
        px = ((((px) as f64) + ((0.5_f64 * ((dsdx) as f64)) * ((dsigma) as f64))) as f32);
        pz = ((((pz) as f64) + ((0.5_f64 * ((dsdz) as f64)) * ((dsigma) as f64))) as f32);
        /* don't let ray exit too close to a vertex */
        xa = (*((*((*(eu)).vu)).v)).y;
        dx = ((*((*((*((*(eu)).euCW)).vu)).v)).y - xa);
        za = (*((*((*(eu)).vu)).v)).x;
        dz = ((*((*((*((*(eu)).euCW)).vu)).v)).x - za);
        frac = (if ((if ((dx) < ((0_i32) as f32)) {
            (-(dx))
        } else {
            (dx)
        }) > (if ((dz) < ((0_i32) as f32)) {
            (-(dz))
        } else {
            (dz)
        })) {
            ((x - xa) / dx)
        } else {
            ((z - za) / dz)
        });
        if (((frac) as f64) < 1.0E-4_f64) {
            x = ((((xa) as f64) + (1.0E-4_f64 * ((dx) as f64))) as f32);
            z = ((((za) as f64) + (1.0E-4_f64 * ((dz) as f64))) as f32);
        } else {
            if (((frac) as f64) > 0.99990000000000001_f64) {
                x = ((((xa) as f64) + (0.99990000000000001_f64 * ((dx) as f64))) as f32);
                z = ((((za) as f64) + (0.99990000000000001_f64 * ((dz) as f64))) as f32);
            }
        }
        /* compute contribution due to attenuation */
        if (qfac != 3.40282347E+38_f32) {
            atten = ((((atten) as f32) + (ddt / qfac)) as f32);
        } else {
            atten = ((0_i32) as f32);
        }
        /* return new raystep */
        (*(rs)).rsNext = ((Box::into_raw(Box::new(std::mem::zeroed::<RayStep>())).cast::<c_void>())
            as *mut RayStep);
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
        (*(rs)).rsNext = std::ptr::null_mut::<RayStep>();
        (*(rs)).atten = atten;
        (*(rs)).ampli = ampli;
        (*(rs)).ampliphase = ampliphase;
        return rs;
    }
}

// traceRayAcrossEdge (C line 855): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn traceRayAcrossEdge(
    mut rs: *mut RayStep,
    mut infofp: *mut CFile,
) -> *mut RayStep {
    unsafe {
        /********************************************************************
         traceRayAcrossEdge - Trace ray across edge.
        *********************************************************************
        Input:
        rs		Pointer to a ray step
        infofp		info file pointer
        Returns:
        If successful, return pointer to new RayStep.
        If unsuccessful, return NULL.

        *********************************************************************
        Notes:
         Failure to trace across an edge is because:
        (1) the ray was incident with angle greater than the critical angle, or
        (2) the ray is incident at a boundary edge
        *********************************************************************
        Author: Andreas Rueger, 1993
        *********************************************************************/

        let mut kmah: i32 = std::mem::zeroed();
        let mut nref: i32 = std::mem::zeroed();
        let mut sigma: f32 = std::mem::zeroed();
        let mut x: f32 = std::mem::zeroed();
        let mut z: f32 = std::mem::zeroed();
        let mut px: f32 = std::mem::zeroed();
        let mut pz: f32 = std::mem::zeroed();
        let mut t: f32 = std::mem::zeroed();
        let mut q1: f32 = std::mem::zeroed();
        let mut p1: f32 = std::mem::zeroed();
        let mut q2: f32 = std::mem::zeroed();
        let mut p2: f32 = std::mem::zeroed();
        let mut temp1: f32 = std::mem::zeroed();
        let mut temp2: f32 = std::mem::zeroed();
        let mut s00: f32 = std::mem::zeroed();
        let mut s1: f32 = std::mem::zeroed();
        let mut ds1dx: f32 = std::mem::zeroed();
        let mut ds1dz: f32 = std::mem::zeroed();
        let mut s2: f32 = std::mem::zeroed();
        let mut ds2dx: f32 = std::mem::zeroed();
        let mut ds2dz: f32 = std::mem::zeroed();
        let mut px1: f32 = std::mem::zeroed();
        let mut pz1: f32 = std::mem::zeroed();
        let mut px1r: f32 = std::mem::zeroed();
        let mut pz1r: f32 = std::mem::zeroed();
        let mut px2r: f32 = std::mem::zeroed();
        let mut pz2r: f32 = std::mem::zeroed();
        let mut pz2rs: f32 = std::mem::zeroed();
        let mut c1ov1: f32 = std::mem::zeroed();
        let mut c2ov2: f32 = std::mem::zeroed();
        let mut oc2s: f32 = std::mem::zeroed();
        let mut g: f32 = std::mem::zeroed();
        let mut cterm: f32 = std::mem::zeroed();
        let mut iterm: f32 = std::mem::zeroed();
        let mut fac1: f32 = std::mem::zeroed();
        let mut fac2: f32 = std::mem::zeroed();
        let mut dv1dl: f32 = std::mem::zeroed();
        let mut dv2dl: f32 = std::mem::zeroed();
        let mut dv1dm: f32 = std::mem::zeroed();
        let mut dv2dm: f32 = std::mem::zeroed();
        let mut scale: f32 = std::mem::zeroed();
        let mut gx: f32 = std::mem::zeroed();
        let mut gz: f32 = std::mem::zeroed();
        let mut hx: f32 = std::mem::zeroed();
        let mut h_z: f32 = std::mem::zeroed();
        let mut frac: f32 = std::mem::zeroed();
        let mut dx: f32 = std::mem::zeroed();
        let mut dz: f32 = std::mem::zeroed();
        let mut taper: f32 = std::mem::zeroed();
        let mut ampli: f32 = std::mem::zeroed();
        let mut atten: f32 = std::mem::zeroed();
        let mut dens1: f32 = std::mem::zeroed();
        let mut dens2: f32 = std::mem::zeroed();
        let mut ampliphase: f32 = std::mem::zeroed();
        let mut coeff: f32 = std::mem::zeroed();
        let mut eu: *mut EdgeUse = std::mem::zeroed();
        let mut eum: *mut EdgeUse = std::mem::zeroed();
        let mut eua: *mut EdgeUseAttributes = std::mem::zeroed();
        let mut euma: *mut EdgeUseAttributes = std::mem::zeroed();
        let mut f: *mut Face = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
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
        if ((*((*(eu)).euMate)).f == ((std::ptr::null_mut::<c_void>()) as *mut Face)) {
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(&*(infofp), "Transmitted outside model.  Ray stopped.\n");
            }
            return std::ptr::null_mut::<RayStep>();
        }
        /* determine sloth on this side of edge */
        fa = (((*(f)).fa) as *mut FaceAttributes);
        s00 = (*(fa)).s00;
        ds1dx = (*(fa)).dsdx;
        ds1dz = (*(fa)).dsdz;
        s1 = ((s00 + (ds1dx * x)) + (ds1dz * z));
        /* determine density on this side */
        dens1 = (*(fa)).dens;
        /* determine sloth on other side of edge */
        eum = (*(eu)).euMate;
        f = (*(eum)).f;
        fa = (((*(f)).fa) as *mut FaceAttributes);
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
        frac = (if ((if ((dx) < ((0_i32) as f32)) {
            (-(dx))
        } else {
            (dx)
        }) > (if ((dz) < ((0_i32) as f32)) {
            (-(dz))
        } else {
            (dz)
        })) {
            ((x - (*((*((*(eu)).vu)).v)).y) / dx)
        } else {
            ((z - (*((*((*(eu)).vu)).v)).x) / dz)
        });
        /* here is speeding up potential */
        /* linearly interpolate unit vector g tangent to edge */
        eua = (((*(eu)).eua) as *mut EdgeUseAttributes);
        euma = (((*(eum)).eua) as *mut EdgeUseAttributes);
        if ((eua != ((std::ptr::null_mut::<c_void>()) as *mut EdgeUseAttributes))
            && (euma != ((std::ptr::null_mut::<c_void>()) as *mut EdgeUseAttributes)))
        {
            gx = ((((frac * (*(euma)).tx) as f64)
                - ((1.0_f64 - ((frac) as f64)) * (((*(eua)).tx) as f64))) as f32);
            gz = ((((frac * (*(euma)).tz) as f64)
                - ((1.0_f64 - ((frac) as f64)) * (((*(eua)).tz) as f64))) as f32);
        } else {
            gx = (-(dx));
            gz = (-(dz));
        }
        scale = ((1.0_f64 / (((gx * gx) + (gz * gz)) as f64).sqrt()) as f32);
        gx = ((((gx) as f32) * scale) as f32);
        gz = ((((gz) as f32) * scale) as f32);
        /* unit vector h normal to edge */
        hx = (-(gz));
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
        if (((pz2rs) as f64) <= 0.0_f64) {
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(&*(infofp), "Post-critical transmission.  Ray stopped.\n");
            }
            return std::ptr::null_mut::<RayStep>();
        }
        pz2r = ((((pz2rs) as f64).sqrt()) as f32);
        /* ray parameters on other side */
        px = ((px2r * h_z) + (pz2r * hx));
        pz = ((pz2r * h_z) - (px2r * hx));
        /* curvature term */
        c1ov1 = pz1r;
        c2ov2 = pz2r;
        oc2s = (s2 / pz2rs);
        if ((eua != ((std::ptr::null_mut::<c_void>()) as *mut EdgeUseAttributes))
            && (euma != ((std::ptr::null_mut::<c_void>()) as *mut EdgeUseAttributes)))
        {
            g = ((((frac * (*(euma)).c) as f64)
                - ((1.0_f64 - ((frac) as f64)) * (((*(eua)).c) as f64))) as f32);
            cterm = ((g * oc2s) * (c1ov1 - c2ov2));
        } else {
            cterm = ((0.0_f64) as f32);
        }
        /* update dynamic ray parameters */
        scale = (((((pz2r) as f64) * ((s1) as f64).sqrt())
            / (((pz1r) as f64) * ((s2) as f64).sqrt())) as f32);
        q1 = (q1 * scale);
        p1 = ((p1 / scale) + (cterm * q1));
        q2 = (q2 * scale);
        p2 = ((p2 / scale) + (cterm * q2));
        /* velocity derivatives tangent and normal to ray */
        fac1 = (((-(0.5_f64)) / ((s1 * s1) as f64)) as f32);
        fac2 = (((-(0.5_f64)) / ((s2 * s2) as f64)) as f32);
        dv1dl = (fac1 * ((px1 * ds1dx) + (pz1 * ds1dz)));
        dv2dl = (fac2 * ((px * ds2dx) + (pz * ds2dz)));
        dv1dm = (fac1 * ((pz1 * ds1dx) - (px1 * ds1dz)));
        dv2dm = (fac2 * ((pz * ds2dx) - (px * ds2dz)));
        /* inhomogeneity term */
        iterm = (((((-(px1r)) * oc2s) as f64)
            * ((2.0_f64 * (((dv1dm * c1ov1) - (dv2dm * c2ov2)) as f64))
                + ((px1r * (dv1dl - dv2dl)) as f64))) as f32);
        /* update dynamic ray parameters */
        p1 = ((((p1) as f32) + (iterm * q1)) as f32);
        p2 = ((((p2) as f32) + (iterm * q2)) as f32);
        /* transmission effects on amplitudes */
        if ((((if ((s1 - s2) < ((0_i32) as f32)) {
            (-(s1 - s2))
        } else {
            (s1 - s2)
        }) as f64)
            > 0.001_f64)
            || (dens1 != dens2))
        {
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(&*(infofp), "Transmitted WITH influence on amplitude.\n");
            }
            if (dens1 == 3.40282347E+38_f32) {
                dens1 = {
                    dens2 = ((1.0_f64) as f32);
                    dens2
                };
            }
            temp1 = (dens2 / dens1);
            temp2 = (pz2r / pz1r);
            coeff = ((1.0_f64 + (((temp1 - temp2) / (temp1 + temp2)) as f64)) as f32);
            /* check if too close to critical */
            taper = ((px1r * px1r) / s2);
            taper = ((if (((taper) as f64) > 0.90000000000000002_f64) {
                (((((taper) as f64) - 0.90000000000000002_f64) * 5.0_f64)
                    * (3.1415926535897931_f64))
                    .cos()
            } else {
                1.0_f64
            }) as f32);
            /* complete amplitude coeff */
            ampli = ((((ampli) as f64)
                * ((((coeff) as f64) * ((pz2r / pz1r) as f64).sqrt()) * ((taper) as f64)))
                as f32);
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(
                    &*(infofp),
                    "  --incident side: v1=%g  dens1=%g  incidence angle=%g\n",
                    (((1_i32) as f64) / ((s1) as f64).sqrt()),
                    ((dens1) as f64),
                    (((((pz1r) as f64) / ((s1) as f64).sqrt()).acos() * ((180_i32) as f64))
                        / (3.1415926535897931_f64))
                );
                fprintf!(
                    &*(infofp),
                    "  --opposite side: v2=%g  dens2=%g  refracted angle=%g\n",
                    (((1_i32) as f64) / ((s2) as f64).sqrt()),
                    ((dens2) as f64),
                    (((((pz2r) as f64) / ((s2) as f64).sqrt()).acos() * ((180_i32) as f64))
                        / (3.1415926535897931_f64))
                );
                fprintf!(&*(infofp), "  --transmission coeff: %g\n", ((coeff) as f64));
                fprintf!(
                    &*(infofp),
                    "  --critical transm. taper: %g\n",
                    ((taper) as f64)
                );
                fprintf!(&*(infofp), "  --total ampli. coeff: %g\n", ((ampli) as f64));
            }
            /* grazing incidence*/
            if (((pz1r * pz1r) as f64) <= (0.0080000000000000002_f64 * ((s1) as f64))) {
                if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                    fprintf!(&*(infofp), "Grazing incidence.  Ray stopped.\n");
                }
                return std::ptr::null_mut::<RayStep>();
            }
        } else {
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(&*(infofp), "Transmitted WITHOUT influencing amplitude.\n");
            }
        }
        /* return new raystep */
        (*(rs)).rsNext = ((Box::into_raw(Box::new(std::mem::zeroed::<RayStep>())).cast::<c_void>())
            as *mut RayStep);
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
        (*(rs)).rsNext = std::ptr::null_mut::<RayStep>();
        (*(rs)).ampli = ampli;
        (*(rs)).ampliphase = ampliphase;
        (*(rs)).atten = atten;
        return rs;
    }
}

// reflectRayFromEdge (C line 1096): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn reflectRayFromEdge(
    mut rs: *mut RayStep,
    mut infofp: *mut CFile,
) -> *mut RayStep {
    unsafe {
        /* Reflect ray from edge and  return pointer to new RayStep. */

        let mut kmah: i32 = std::mem::zeroed();
        let mut nref: i32 = std::mem::zeroed();
        let mut sigma: f32 = std::mem::zeroed();
        let mut x: f32 = std::mem::zeroed();
        let mut z: f32 = std::mem::zeroed();
        let mut px: f32 = std::mem::zeroed();
        let mut pz: f32 = std::mem::zeroed();
        let mut t: f32 = std::mem::zeroed();
        let mut q1: f32 = std::mem::zeroed();
        let mut p1: f32 = std::mem::zeroed();
        let mut q2: f32 = std::mem::zeroed();
        let mut p2: f32 = std::mem::zeroed();
        let mut s00: f32 = std::mem::zeroed();
        let mut dsdx: f32 = std::mem::zeroed();
        let mut dsdz: f32 = std::mem::zeroed();
        let mut s: f32 = std::mem::zeroed();
        let mut coeff: f32 = std::mem::zeroed();
        let mut px1: f32 = std::mem::zeroed();
        let mut pz1: f32 = std::mem::zeroed();
        let mut pxr: f32 = std::mem::zeroed();
        let mut pzr: f32 = std::mem::zeroed();
        let mut c1ov1: f32 = std::mem::zeroed();
        let mut c2ov2: f32 = std::mem::zeroed();
        let mut oc2s: f32 = std::mem::zeroed();
        let mut g: f32 = std::mem::zeroed();
        let mut cterm: f32 = std::mem::zeroed();
        let mut iterm: f32 = std::mem::zeroed();
        let mut dv1dl: f32 = std::mem::zeroed();
        let mut dv2dl: f32 = std::mem::zeroed();
        let mut dv1dm: f32 = std::mem::zeroed();
        let mut dv2dm: f32 = std::mem::zeroed();
        let mut scale: f32 = std::mem::zeroed();
        let mut gx: f32 = std::mem::zeroed();
        let mut gz: f32 = std::mem::zeroed();
        let mut hx: f32 = std::mem::zeroed();
        let mut h_z: f32 = std::mem::zeroed();
        let mut frac: f32 = std::mem::zeroed();
        let mut dx: f32 = std::mem::zeroed();
        let mut dz: f32 = std::mem::zeroed();
        let mut atten: f32 = std::mem::zeroed();
        let mut ampli: f32 = std::mem::zeroed();
        let mut ampliphase: f32 = std::mem::zeroed();
        let mut dens1: f32 = std::mem::zeroed();
        let mut dens2: f32 = std::mem::zeroed();
        let mut s2: f32 = std::mem::zeroed();
        let mut ds2dx: f32 = std::mem::zeroed();
        let mut ds2dz: f32 = std::mem::zeroed();
        let mut temp1: f32 = std::mem::zeroed();
        let mut temp2: f32 = std::mem::zeroed();
        let mut eu: *mut EdgeUse = std::mem::zeroed();
        let mut eum: *mut EdgeUse = std::mem::zeroed();
        let mut eua: *mut EdgeUseAttributes = std::mem::zeroed();
        let mut euma: *mut EdgeUseAttributes = std::mem::zeroed();
        let mut f: *mut Face = std::mem::zeroed();
        let mut r#fn: *mut Face = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
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
        fa = (((*(f)).fa) as *mut FaceAttributes);
        s00 = (*(fa)).s00;
        dsdx = (*(fa)).dsdx;
        dsdz = (*(fa)).dsdz;
        s = ((s00 + (dsdx * x)) + (dsdz * z));
        /* determine dens on incident side of edge */
        dens1 = (*(fa)).dens;
        /* edge vector */
        eum = (*(eu)).euMate;
        /* determine sloth on other side of edge */
        r#fn = (*(eum)).f;
        fa = (((*(r#fn)).fa) as *mut FaceAttributes);
        s00 = (*(fa)).s00;
        ds2dx = (*(fa)).dsdx;
        ds2dz = (*(fa)).dsdz;
        s2 = ((s00 + (ds2dx * x)) + (ds2dz * z));
        /* determine density on other side of the edge */
        dens2 = (*(fa)).dens;
        if (dens1 == 3.40282347E+38_f32) {
            dens1 = ((1.0_f64) as f32);
        }
        if (dens2 == 3.40282347E+38_f32) {
            dens2 = ((1.0_f64) as f32);
        }
        dx = ((*((*((*(eum)).vu)).v)).y - (*((*((*(eu)).vu)).v)).y);
        dz = ((*((*((*(eum)).vu)).v)).x - (*((*((*(eu)).vu)).v)).x);
        /* fractional distance along edge */
        frac = (if ((if ((dx) < ((0_i32) as f32)) {
            (-(dx))
        } else {
            (dx)
        }) > (if ((dz) < ((0_i32) as f32)) {
            (-(dz))
        } else {
            (dz)
        })) {
            ((x - (*((*((*(eu)).vu)).v)).y) / dx)
        } else {
            ((z - (*((*((*(eu)).vu)).v)).x) / dz)
        });
        /* linearly interpolate unit vector g tangent to edge */
        eua = (((*(eu)).eua) as *mut EdgeUseAttributes);
        euma = (((*(eum)).eua) as *mut EdgeUseAttributes);
        gx = ((((frac * (*(euma)).tx) as f64)
            - ((1.0_f64 - ((frac) as f64)) * (((*(eua)).tx) as f64))) as f32);
        gz = ((((frac * (*(euma)).tz) as f64)
            - ((1.0_f64 - ((frac) as f64)) * (((*(eua)).tz) as f64))) as f32);
        scale = ((1.0_f64 / (((gx * gx) + (gz * gz)) as f64).sqrt()) as f32);
        gx = ((((gx) as f32) * scale) as f32);
        gz = ((((gz) as f32) * scale) as f32);
        /* unit vector h normal to edge */
        hx = (-(gz));
        h_z = gx;
        /* remember incident ray parameters */
        px1 = px;
        pz1 = pz;
        /* rotated incident ray parameters */
        pxr = ((px * h_z) - (pz * hx));
        pzr = ((px * hx) + (pz * h_z));
        /* rotated reflected ray parameters */
        pxr = pxr;
        pzr = (-(pzr));
        /* reflected ray parameters */
        px = ((pxr * h_z) + (pzr * hx));
        pz = ((pzr * h_z) - (pxr * hx));
        /* curvature term */
        c1ov1 = {
            c2ov2 = (-(pzr));
            c2ov2
        };
        oc2s = (s / (pzr * pzr));
        g = ((((frac * (*(euma)).c) as f64) - ((1.0_f64 - ((frac) as f64)) * (((*(eua)).c) as f64)))
            as f32);
        cterm = ((g * oc2s) * (c1ov1 + c2ov2));
        /* update dynamic ray parameters */
        scale = ((-(c2ov2)) / c1ov1);
        q1 = (q1 * scale);
        p1 = ((p1 / scale) + (cterm * q1));
        q2 = (q2 * scale);
        p2 = ((p2 / scale) + (cterm * q2));
        /* if sloth not constant */
        if ((((dsdx) as f64) != 0.0_f64) || (((dsdz) as f64) != 0.0_f64)) {
            /* velocity derivatives tangent and normal to ray */
            scale = (((-(0.5_f64)) / ((s * s) as f64)) as f32);
            dv1dl = (scale * ((px1 * dsdx) + (pz1 * dsdz)));
            dv2dl = (scale * ((px * dsdx) + (pz * dsdz)));
            dv1dm = (scale * ((pz1 * dsdx) - (px1 * dsdz)));
            dv2dm = (scale * ((pz * dsdx) - (px * dsdz)));
            /* inhomogeneity term */
            iterm = (((((-(pxr)) * oc2s) as f64)
                * ((2.0_f64 * (((dv1dm * c1ov1) + (dv2dm * c2ov2)) as f64))
                    + ((pxr * (dv1dl - dv2dl)) as f64))) as f32);
            /* update dynamic ray parameters */
            p1 = ((((p1) as f32) + (iterm * q1)) as f32);
            p2 = ((((p2) as f32) + (iterm * q2)) as f32);
        }
        /* increment number of reflections */
        nref = nref + 1;
        if ((s2 - (pxr * pxr)) >= ((0_i32) as f32)) {
            /* pre-critical reflection */
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(
                    &*(infofp),
                    "Reflected pre-critically at %g degrees.\n",
                    (((((-(pzr)) as f64) / ((s) as f64).sqrt()).acos() * ((180_i32) as f64))
                        / (3.1415926535897931_f64))
                );
            }
            temp1 = ((((((-(pzr)) * dens2) / dens1) as f64) / ((s) as f64).sqrt()) as f32);
            temp2 = (((((s2 / s) - ((pxr * pxr) / s)) as f64).sqrt()) as f32);
            coeff = ((temp1 - temp2) / (temp1 + temp2));
            ampli = ((((ampli) as f32) * coeff) as f32);
        } else {
            /* post-critical reflection */
            coeff = ((1_i32) as f32);
            ampliphase = ((((ampliphase) as f64)
                - (((2_i32) as f64)
                    * (((-((((pxr * pxr) - s2) as f64).sqrt())) * ((dens1) as f64))
                        / ((dens2 * pzr) as f64))
                        .atan())) as f32);
            if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                fprintf!(
                    &*(infofp),
                    "Reflected post-critically at %g degrees.\n",
                    (((((-(pzr)) as f64) / ((s) as f64).sqrt()).acos() * ((180_i32) as f64))
                        / (3.1415926535897931_f64))
                );
            }
        }
        if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
            fprintf!(
                &*(infofp),
                "  --incident side: v1=%g  dens1=%g\n",
                (((1_i32) as f64) / ((s) as f64).sqrt()),
                ((dens1) as f64)
            );
            fprintf!(
                &*(infofp),
                "  --opposite side: v2=%g  dens2=%g\n",
                (((1_i32) as f64) / ((s2) as f64).sqrt()),
                ((dens2) as f64)
            );
            fprintf!(&*(infofp), "  --reflection coeff: %g\n", ((coeff) as f64));
            fprintf!(&*(infofp), "  --total ampl coeff: %g\n", ((ampli) as f64));
            fprintf!(
                &*(infofp),
                "  --total phase shift: %g\n",
                ((ampliphase) as f64)
            );
        }
        /* return new raystep */
        (*(rs)).rsNext = ((Box::into_raw(Box::new(std::mem::zeroed::<RayStep>())).cast::<c_void>())
            as *mut RayStep);
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
        (*(rs)).rsNext = std::ptr::null_mut::<RayStep>();
        (*(rs)).ampli = ampli;
        (*(rs)).ampliphase = ampliphase;
        (*(rs)).atten = atten;
        return rs;
    }
}

// evaluateDynamic (C line 1284): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn evaluateDynamic(
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
    unsafe {
        /* evaluate dynamic ray parameters */

        let mut s0: f64 = std::mem::zeroed();
        let mut s1: f64 = std::mem::zeroed();
        let mut s2: f64 = std::mem::zeroed();
        let mut ss: f64 = std::mem::zeroed();
        let mut scale: f64 = std::mem::zeroed();
        let mut a: f64 = std::mem::zeroed();
        let mut b: f64 = std::mem::zeroed();
        let mut c: f64 = std::mem::zeroed();
        let mut d: f64 = std::mem::zeroed();
        let mut e: f64 = std::mem::zeroed();
        let mut q1i: f64 = std::mem::zeroed();
        let mut p1i: f64 = std::mem::zeroed();
        let mut q2i: f64 = std::mem::zeroed();
        let mut p2i: f64 = std::mem::zeroed();
        let mut q1o: f64 = std::mem::zeroed();
        let mut p1o: f64 = std::mem::zeroed();
        let mut q2o: f64 = std::mem::zeroed();
        let mut p2o: f64 = std::mem::zeroed();
        /* get input dynamic ray parameters */
        q1i = ((*(q1)) as f64);
        p1i = ((*(p1)) as f64);
        q2i = ((*(q2)) as f64);
        p2i = ((*(p2)) as f64);
        /* trivial case:  constant sloth or ray parallel to sloth gradient */
        if ((pz * dsdx) == (px * dsdz)) {
            q1o = (q1i + (p1i * ((sigma) as f64)));
            q2o = (q2i + (p2i * ((sigma) as f64)));
            p2o = p2i;
            p1o = p1i;
        } else {
            /* constants */
            s0 = (((px * px) + (pz * pz)) as f64);
            s1 = (((px * dsdx) + (pz * dsdz)) as f64);
            s2 = (0.25_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64));
            ss = ((((s2 * ((sigma) as f64)) + s1) * ((sigma) as f64)) + s0);
            scale = (1.0_f64 / (s0 * ss).sqrt());
            a = (s0 + ((0.5_f64 * s1) * ((sigma) as f64)));
            b = (((((0.25_f64 * s1) * s1) / s0) - s2) * ((sigma) as f64));
            c = (s0 + (s1 * ((sigma) as f64)));
            d = ((0.5_f64 * s1) + (2.0_f64 * b));
            b = ((((b) as f64) * ((sigma) as f64)) as f64);
            e = (((0.5_f64 * s1) + (s2 * ((sigma) as f64))) / ss);
            /* update q1 and q2 */
            q1o = (scale * ((a * (q1i + (p1i * ((sigma) as f64)))) + (b * q1i)));
            q2o = (scale * ((a * (q2i + (p2i * ((sigma) as f64)))) + (b * q2i)));
            /* update p1 and p2 */
            p1o = ((scale * ((c * p1i) + (d * q1i))) - (e * q1o));
            p2o = ((scale * ((c * p2i) + (d * q2i))) - (e * q2o));
        }
        /* update kmah index */
        if ((((q2i * q2o) >= 0.0_f64) && ((p2i * p2o) < 0.0_f64)) && ((q2i * p2i) < 0.0_f64)) {
            (*(kmah)) = ((((*(kmah)) as i32) + 2_i32) as i32);
        } else {
            if (((q2o == 0.0_f64) || ((q2i * q2o) < 0.0_f64))
                || (((q2i == 0.0_f64) && ((p2i * p2o) < 0.0_f64)) && ((q2o * p2o) > 0.0_f64)))
            {
                (*(kmah)) = ((((*(kmah)) as i32) + 1_i32) as i32);
            }
        }
        /* return updated dynamic ray parameters */
        (*(q1)) = ((q1o) as f32);
        (*(p1)) = ((p1o) as f32);
        (*(q2)) = ((q2o) as f32);
        (*(p2)) = ((p2o) as f32);
    }
}

// evaluateDynamic2 (C line 1338): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn evaluateDynamic2(
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
    unsafe {
        /* evaluate dynamic ray parameters at each point of the ray */

        let mut ss: f64 = std::mem::zeroed();
        let mut scale: f64 = std::mem::zeroed();
        let mut a: f64 = std::mem::zeroed();
        let mut b: f64 = std::mem::zeroed();
        let mut c: f64 = std::mem::zeroed();
        let mut d: f64 = std::mem::zeroed();
        let mut e: f64 = std::mem::zeroed();
        let mut q1i: f64 = std::mem::zeroed();
        let mut p1i: f64 = std::mem::zeroed();
        let mut q2i: f64 = std::mem::zeroed();
        let mut p2i: f64 = std::mem::zeroed();
        let mut q1o: f64 = std::mem::zeroed();
        let mut p1o: f64 = std::mem::zeroed();
        let mut q2o: f64 = std::mem::zeroed();
        let mut p2o: f64 = std::mem::zeroed();
        /* get input dynamic ray parameters */
        q1i = ((*(q1)) as f64);
        p1i = ((*(p1)) as f64);
        q2i = ((*(q2)) as f64);
        p2i = ((*(p2)) as f64);
        /* trivial case:  constant sloth or ray parallel to sloth gradient */
        if ((pz * dsdx) == (px * dsdz)) {
            q1o = (q1i + (p1i * ((sigma) as f64)));
            q2o = (q2i + (p2i * ((sigma) as f64)));
            p2o = p2i;
            p1o = p1i;
        } else {
            /* constants */
            ss = (((((s2 * sigma) + s1) * sigma) + s0) as f64);
            scale = (1.0_f64 / (((s0) as f64) * ss).sqrt());
            a = (((s0) as f64) + ((0.5_f64 * ((s1) as f64)) * ((sigma) as f64)));
            b = (((((0.25_f64 * ((s1) as f64)) * ((s1) as f64)) / ((s0) as f64)) - ((s2) as f64))
                * ((sigma) as f64));
            c = ((s0 + (s1 * sigma)) as f64);
            d = ((0.5_f64 * ((s1) as f64)) + (2.0_f64 * b));
            b = ((((b) as f64) * ((sigma) as f64)) as f64);
            e = (((0.5_f64 * ((s1) as f64)) + ((s2 * sigma) as f64)) / ss);
            /* update q1 and q2 */
            q1o = (scale * ((a * (q1i + (p1i * ((sigma) as f64)))) + (b * q1i)));
            q2o = (scale * ((a * (q2i + (p2i * ((sigma) as f64)))) + (b * q2i)));
            /* update p1 and p2 */
            p1o = ((scale * ((c * p1i) + (d * q1i))) - (e * q1o));
            p2o = ((scale * ((c * p2i) + (d * q2i))) - (e * q2o));
        }
        /* return updated dynamic ray parameters */
        (*(q1)) = ((q1o) as f32);
        (*(p1)) = ((p1o) as f32);
        (*(q2)) = ((q2o) as f32);
        (*(p2)) = ((p2o) as f32);
    }
}

// checkIfSourceOnEdge (C line 1383): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn checkIfSourceOnEdge(mut tris: *mut Face, mut zs: f32, mut xs: f32) -> i32 {
    unsafe {
        let mut x1: f32 = std::mem::zeroed();
        let mut x2: f32 = std::mem::zeroed();
        let mut x3: f32 = std::mem::zeroed();
        let mut z1: f32 = std::mem::zeroed();
        let mut z2: f32 = std::mem::zeroed();
        let mut z3: f32 = std::mem::zeroed();
        let mut m12: f32 = std::mem::zeroed();
        let mut m13: f32 = std::mem::zeroed();
        let mut m23: f32 = std::mem::zeroed();
        let mut b12: f32 = std::mem::zeroed();
        let mut b13: f32 = std::mem::zeroed();
        let mut b23: f32 = std::mem::zeroed();
        let mut eps: f32 = std::mem::zeroed();
        let mut eu: *mut EdgeUse = std::mem::zeroed();
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
        if (((if ((x1 - x2) < ((0_i32) as f32)) {
            (-(x1 - x2))
        } else {
            (x1 - x2)
        }) as f64)
            < (0.10000000000000001_f64 * ((eps) as f64)))
        {
            if (((if ((xs - x1) < ((0_i32) as f32)) {
                (-(xs - x1))
            } else {
                (xs - x1)
            }) as f64)
                < (0.10000000000000001_f64 * ((eps) as f64)))
            {
                return 1_i32;
            } else {
                m12 = ((0.01_f64 * ((3.40282347E+38_f32) as f64)) as f32);
            }
        } else {
            m12 = ((z1 - z2) / (x1 - x2));
        }
        if (((if ((x1 - x3) < ((0_i32) as f32)) {
            (-(x1 - x3))
        } else {
            (x1 - x3)
        }) as f64)
            < (0.10000000000000001_f64 * ((eps) as f64)))
        {
            if (((if ((xs - x1) < ((0_i32) as f32)) {
                (-(xs - x1))
            } else {
                (xs - x1)
            }) as f64)
                < (0.10000000000000001_f64 * ((eps) as f64)))
            {
                return 1_i32;
            } else {
                m13 = ((0.01_f64 * ((3.40282347E+38_f32) as f64)) as f32);
            }
        } else {
            m13 = ((z1 - z3) / (x1 - x3));
        }
        if (((if ((x2 - x3) < ((0_i32) as f32)) {
            (-(x2 - x3))
        } else {
            (x2 - x3)
        }) as f64)
            < (0.10000000000000001_f64 * ((eps) as f64)))
        {
            if (((if ((xs - x2) < ((0_i32) as f32)) {
                (-(xs - x2))
            } else {
                (xs - x2)
            }) as f64)
                < (0.10000000000000001_f64 * ((eps) as f64)))
            {
                return 1_i32;
            } else {
                m23 = ((0.01_f64 * ((3.40282347E+38_f32) as f64)) as f32);
            }
        } else {
            m23 = ((z2 - z3) / (x2 - x3));
        }
        b12 = (z1 - (m12 * x1));
        b13 = (z1 - (m13 * x1));
        b23 = (z2 - (m23 * x2));
        /* source on edge? */
        if (((((if (((zs - (m12 * xs)) - b12) < ((0_i32) as f32)) {
            (-((zs - (m12 * xs)) - b12))
        } else {
            ((zs - (m12 * xs)) - b12)
        }) as f64)
            < (0.10000000000000001_f64 * ((eps) as f64)))
            || (((if (((zs - (m13 * xs)) - b13) < ((0_i32) as f32)) {
                (-((zs - (m13 * xs)) - b13))
            } else {
                ((zs - (m13 * xs)) - b13)
            }) as f64)
                < (0.10000000000000001_f64 * ((eps) as f64))))
            || (((if (((zs - (m23 * xs)) - b23) < ((0_i32) as f32)) {
                (-((zs - (m23 * xs)) - b23))
            } else {
                ((zs - (m23 * xs)) - b23)
            }) as f64)
                < (0.10000000000000001_f64 * ((eps) as f64))))
        {
            return 1_i32;
        }
        return 0_i32;
    }
}

// writeFresnel (C line 1444): preserve statement order and C arithmetic conversions.
pub(super) unsafe fn writeFresnel(
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
    mut outparfp: *mut CFile,
) {
    unsafe {
        /* for each ray, write x,z pairs and FresnelVolume */

        let mut iray: i32 = std::mem::zeroed();
        let mut nrs: i32 = std::mem::zeroed();
        let mut irs: i32 = std::mem::zeroed();
        let mut nxzw: i32 = std::mem::zeroed();
        let mut ns: i32 = std::mem::zeroed();
        let mut is: i32 = std::mem::zeroed();
        let mut icount: i32 = std::mem::zeroed();
        let mut sigmamax: f32 = ((0.0_f64) as f32);
        let mut sigma1: f32 = std::mem::zeroed();
        let mut sigma2: f32 = std::mem::zeroed();
        let mut x: f32 = std::mem::zeroed();
        let mut z: f32 = std::mem::zeroed();
        let mut px: f32 = std::mem::zeroed();
        let mut pz: f32 = std::mem::zeroed();
        let mut dsdx: f32 = std::mem::zeroed();
        let mut dsdz: f32 = std::mem::zeroed();
        let mut dummy: f32 = std::mem::zeroed();
        let mut ds: f32 = std::mem::zeroed();
        let mut xs: f32 = std::mem::zeroed();
        let mut zs: f32 = std::mem::zeroed();
        let mut dsigma: f32 = std::mem::zeroed();
        let mut mfa: f32 = std::mem::zeroed();
        let mut mfb: f32 = std::mem::zeroed();
        let mut scale: f32 = std::mem::zeroed();
        let mut q1end: f32 = ((0.0_f64) as f32);
        let mut q2end: f32 = ((0.0_f64) as f32);
        let mut q1: f32 = std::mem::zeroed();
        let mut q2: f32 = std::mem::zeroed();
        let mut p1: f32 = std::mem::zeroed();
        let mut p2: f32 = std::mem::zeroed();
        let mut s0: f32 = std::mem::zeroed();
        let mut s1: f32 = std::mem::zeroed();
        let mut s2: f32 = std::mem::zeroed();
        let mut q1i: f32 = std::mem::zeroed();
        let mut q2i: f32 = std::mem::zeroed();
        let mut p1i: f32 = std::mem::zeroed();
        let mut p2i: f32 = std::mem::zeroed();
        let mut xfresnel1: f32 = std::mem::zeroed();
        let mut xfresnel2: f32 = std::mem::zeroed();
        let mut zfresnel1: f32 = std::mem::zeroed();
        let mut zfresnel2: f32 = std::mem::zeroed();
        let mut xa: *mut f32 = std::mem::zeroed();
        let mut za: *mut f32 = std::mem::zeroed();
        let mut r: *mut f32 = std::mem::zeroed();
        let mut rsi: *mut RayStep = std::mem::zeroed();
        let mut fa: *mut FaceAttributes = std::mem::zeroed();
        m = (m).offset((0_i32) as isize);
        /* dummy */
        /* allocate space for auxiliary arrays */
        xa = alloc_float(((nxz) as usize));
        za = alloc_float(((nxz) as usize));
        r = alloc_float(((nxz) as usize));
        /* initialize counter */
        icount = 0_i32;
        /* loop over rays */
        {
            iray = 0_i32;
            while (iray < nray) {
                {
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
                            while (rsi != ((std::ptr::null_mut::<c_void>()) as *mut RayStep)) {
                                {
                                    sigmamax = (*(rsi)).sigma;
                                    q1end = (*(rsi)).q1;
                                    q2end = (*(rsi)).q2;
                                }
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
                                {
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
                                    fa = (((*((*(rsi)).f)).fa) as *mut FaceAttributes);
                                    dsdx = (*(fa)).dsdx;
                                    dsdz = (*(fa)).dsdz;
                                    /* constants for nontrivial case */
                                    s0 = ((px * px) + (pz * pz));
                                    s1 = ((px * dsdx) + (pz * dsdz));
                                    s2 = ((0.25_f64 * (((dsdx * dsdx) + (dsdz * dsdz)) as f64))
                                        as f32);
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
                                    ns = ((((1_i32) as f64)
                                        + ((((sigma2 - sigma1) * ((nxz - 1_i32) as f32)) as f64)
                                            / (2.0_f64 * ((sigmamax) as f64))))
                                        as i32);
                                    /* increment in sigma */
                                    ds = ((sigma2 - sigma1) / ((ns) as f32));
                                    /* loop over sigma */
                                    {
                                        is = 0_i32;
                                        dsigma = ((0.0_f64) as f32);
                                        while (is < ns) {
                                            {
                                                /* compute and write x and z */
                                                xs = ((((x) as f64)
                                                    + (((dsigma) as f64)
                                                        * (((px) as f64)
                                                            + ((((dsigma) as f64) * 0.25_f64)
                                                                * ((dsdx) as f64)))))
                                                    as f32);
                                                zs = ((((z) as f64)
                                                    + (((dsigma) as f64)
                                                        * (((pz) as f64)
                                                            + ((((dsigma) as f64) * 0.25_f64)
                                                                * ((dsdz) as f64)))))
                                                    as f32);
                                                (*(xa).offset((nxzw) as isize)) = xs;
                                                (*(za).offset((nxzw) as isize)) = zs;
                                                if (efwrite(
                                                    std::slice::from_raw_parts(
                                                        ((&raw mut zs) as *mut c_void)
                                                            as *const f32,
                                                        ((1_i32) as usize),
                                                    ),
                                                    &*(rayfp),
                                                ) != ((1_i32) as usize))
                                                {
                                                    err!("error writing ray file!\n");
                                                }
                                                if (efwrite(
                                                    std::slice::from_raw_parts(
                                                        ((&raw mut xs) as *mut c_void)
                                                            as *const f32,
                                                        ((1_i32) as usize),
                                                    ),
                                                    &*(rayfp),
                                                ) != ((1_i32) as usize))
                                                {
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
                                                if (q2 == ((0_i32) as f32)) {
                                                    mfa = 3.40282347E+38_f32;
                                                } else {
                                                    mfa = (p2 / q2);
                                                }
                                                /* compute m(Of,B) */
                                                mfb = (((-(q1)) * q2end) + (q2 * q1end));
                                                if (mfb == ((0_i32) as f32)) {
                                                    mfb = 3.40282347E+38_f32;
                                                } else {
                                                    mfb = (((p2 * q1end) - (p1 * q2end)) / mfb);
                                                }
                                                /* compute fresnel radius */
                                                dummy = ((if ((1.0_f64
                                                    / (((mfa - mfb) * freq) as f64))
                                                    < ((0_i32) as f64))
                                                {
                                                    (-(1.0_f64 / (((mfa - mfb) * freq) as f64)))
                                                } else {
                                                    (1.0_f64 / (((mfa - mfb) * freq) as f64))
                                                })
                                                    as f32);
                                                (*(r).offset((nxzw) as isize)) =
                                                    ((((dummy) as f64).sqrt()) as f32);
                                                /* if sufficient number */
                                                /* written, break */
                                                if ({
                                                    nxzw = nxzw + 1;
                                                    nxzw
                                                } >= nxz)
                                                {
                                                    break;
                                                }
                                            }
                                            is = is + 1;
                                            dsigma = ((((dsigma) as f32) + ds) as f32);
                                        }
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
                            (*(r).offset((nxzw) as isize)) = ((0_i32) as f32);
                            if (efwrite(
                                std::slice::from_raw_parts(
                                    ((&raw mut zs) as *mut c_void) as *const f32,
                                    ((1_i32) as usize),
                                ),
                                &*(rayfp),
                            ) != ((1_i32) as usize))
                            {
                                err!("error writing ray file!\n");
                            }
                            if (efwrite(
                                std::slice::from_raw_parts(
                                    ((&raw mut xs) as *mut c_void) as *const f32,
                                    ((1_i32) as usize),
                                ),
                                &*(rayfp),
                            ) != ((1_i32) as usize))
                            {
                                err!("error writing ray file!\n");
                            }
                            nxzw = nxzw + 1;
                        }
                        /* compute fresnel points */
                        if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
                            fprintf!(&*(infofp), "Fresnel points:\n");
                        }
                        {
                            is = 0_i32;
                            while (is < nxz) {
                                {
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
                                        scale = (((((((*(za).offset((is - 1_i32) as isize))
                                            - (*(za).offset((is + 1_i32) as isize)))
                                            * ((*(za).offset((is - 1_i32) as isize))
                                                - (*(za).offset((is + 1_i32) as isize))))
                                            + (((*(xa).offset((is + 1_i32) as isize))
                                                - (*(xa).offset((is - 1_i32) as isize)))
                                                * ((*(xa).offset((is + 1_i32) as isize))
                                                    - (*(xa).offset((is - 1_i32) as isize)))))
                                            as f64)
                                            .sqrt())
                                            as f32);
                                        if (((scale == ((0_i32) as f32)) || (is == 0_i32))
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
                                            dummy = ((((1.0_f64 / ((scale) as f64))
                                                * (((*(za).offset((is - 1_i32) as isize))
                                                    - (*(za).offset((is + 1_i32) as isize)))
                                                    as f64))
                                                * ((*(r).offset((is) as isize)) as f64))
                                                as f32);
                                            xfresnel1 = ((*(xa).offset((is) as isize)) + dummy);
                                            xfresnel2 = ((*(xa).offset((is) as isize)) - dummy);
                                            dummy = ((((1.0_f64 / ((scale) as f64))
                                                * (((*(xa).offset((is + 1_i32) as isize))
                                                    - (*(xa).offset((is - 1_i32) as isize)))
                                                    as f64))
                                                * ((*(r).offset((is) as isize)) as f64))
                                                as f32);
                                            zfresnel1 = ((*(za).offset((is) as isize)) + dummy);
                                            zfresnel2 = ((*(za).offset((is) as isize)) - dummy);
                                        }
                                    }
                                    if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile))
                                    {
                                        fprintf!(
                                            &*(infofp),
                                            "x1=%g  z1=%g  x2=%g  z2=%g  r=%g\n",
                                            ((xfresnel1) as f64),
                                            ((zfresnel1) as f64),
                                            ((xfresnel2) as f64),
                                            ((zfresnel2) as f64),
                                            ((*(r).offset((is) as isize)) as f64)
                                        );
                                    }
                                    if (efwrite(
                                        std::slice::from_raw_parts(
                                            ((&raw mut zfresnel1) as *mut c_void) as *const f32,
                                            ((1_i32) as usize),
                                        ),
                                        &*(fresnelfp),
                                    ) != ((1_i32) as usize))
                                    {
                                        err!("error writing Fresnel file!\n");
                                    }
                                    if (efwrite(
                                        std::slice::from_raw_parts(
                                            ((&raw mut xfresnel1) as *mut c_void) as *const f32,
                                            ((1_i32) as usize),
                                        ),
                                        &*(fresnelfp),
                                    ) != ((1_i32) as usize))
                                    {
                                        err!("error writing Fresnel file!\n");
                                    }
                                    if (efwrite(
                                        std::slice::from_raw_parts(
                                            ((&raw mut zfresnel2) as *mut c_void) as *const f32,
                                            ((1_i32) as usize),
                                        ),
                                        &*(fresnelfp),
                                    ) != ((1_i32) as usize))
                                    {
                                        err!("error writing Fresnel file!\n");
                                    }
                                    if (efwrite(
                                        std::slice::from_raw_parts(
                                            ((&raw mut xfresnel2) as *mut c_void) as *const f32,
                                            ((1_i32) as usize),
                                        ),
                                        &*(fresnelfp),
                                    ) != ((1_i32) as usize))
                                    {
                                        err!("error writing Fresnel file!\n");
                                    }
                                }
                                is = is + 1;
                            }
                        }
                    } else {
                        /* do not write rays and fresnel volumes; */
                        /* manipulate rayend information */
                        (*(re).offset((iray) as isize)).kend = (-(1_i32));
                    }
                }
                iray = iray + 1;
            }
        }
        /* close loop over rays */
        if (infofp != ((std::ptr::null_mut::<c_void>()) as *mut CFile)) {
            fprintf!(&*(infofp), "\nWrote %d rays to ray file.\n", icount);
        }
        /* write number of rays to parameter file */
        fprintf!(&*(outparfp), "%d\n", icount);
        /* free space */
        drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
            xa,
            nxz as usize,
        )));
        drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
            za,
            nxz as usize,
        )));
        drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
            r,
            nxz as usize,
        )));
    }
}
