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

//! Private helpers from Trielas/main/elaray.c; Andreas Rueger, 1994–1995.
//! Based on gbray.c (Andreas Rueger) and sdray.c (Dave Hale).
//! Original additional notices: Copyright (c) Colorado School of Mines,
//! 1999 and 2011. All rights reserved.

use super::*;
/// Shoot rays from triangle to triangle, adding a step on each side of an edge.
pub(super) unsafe fn shootRaysAni(
    mut m: *mut Model,
    mut xs: f32,
    mut zs: f32,
    mut nangle: i32,
    mut reftrans: i32,
    mut dangle: f32,
    mut fangle: f32,
    mut rc: *mut RefCheck,
    mut kk: i32,
    mut rsp: *mut *mut RayStep,
    mut re: *mut RayEnd,
    mut nrefseq: *mut i32,
    mut ifp: *mut CFile,
    mut mode: i32,
    mut fdip: f32,
    mut f0: f32,
    mut fazi: f32,
) {
    // The C reflection diagnostic reuses a stack slot across calls and rays.
    let mut coeff_slot = 0.0f32;
    let mut iangle: i32 = 0;
    let mut refl: i32 = 0;
    let mut stop: i32 = 0;
    let mut i: i32 = 0;
    let mut temp: i32 = 0;
    let mut nameref: i32 = 0;
    let mut free_surface: i32 = 0;
    let mut refmod: i32 = 0;
    let mut transmod: i32 = 0;
    let mut ireset: *mut i32 = std::ptr::null_mut();
    let mut mindex: i32 = 0;
    let mut a1111: f32 = 0.0;
    let mut a3333: f32 = 0.0;
    let mut a1133: f32 = 0.0;
    let mut a1313: f32 = 0.0;
    let mut a1113: f32 = 0.0;
    let mut a3313: f32 = 0.0;
    let mut a1212: f32 = 0.0;
    let mut a1223: f32 = 0.0;
    let mut a2323: f32 = 0.0;
    let mut polx: f32 = 0.0;
    let mut polz: f32 = 0.0;
    let mut rhob: f32 = 0.0;
    let mut s: f32 = 0.0;
    let mut ss: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut cc: f32 = 0.0;
    let mut sc: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut sa: f32 = 0.0;
    let mut soa: f32 = 0.0;
    let mut angle: f32 = 0.0;
    let mut gamm11: f32 = 0.0;
    let mut gamm33: f32 = 0.0;
    let mut gamm13: f32 = 0.0;
    let mut sqr: f32 = 0.0;
    let mut vp2: f32 = 0.0;
    let mut vp: f32 = 0.0;
    let mut ovp: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g33: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut fx: f32 = 0.0;
    let mut fz: f32 = 0.0;
    let mut fy: f32 = 0.0;
    let mut junkfp: *mut CFile = std::ptr::null_mut();
    let mut tris: *mut Tri = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    let mut ea: *mut EdgeAttributes = std::ptr::null_mut();
    let mut rs: *mut RayStep = std::ptr::null_mut();
    let mut rslast: *mut RayStep = std::ptr::null_mut();
    ireset = Vec::leak(vec![0i32; ((kk) as usize) as usize]).as_mut_ptr();
    junkfp = std::ptr::null_mut();
    pz = ((0f64) as f32);
    px = pz;
    nameref = 0;
    tris = insideTriInModel(m, std::ptr::null_mut(), zs, xs);
    temp = 1;
    i = 0;
    while (temp == 1) {
        temp = testIfSourceOnEdge(tris, (&raw mut zs), (&raw mut xs));
        if (temp) != 0 {
            tris = insideTriInModel(m, std::ptr::null_mut(), zs, xs);
            i = i + 1;
            if (i == 5) {
                err!("\n ERROR in testIfSourceOnEdge\n");
            }
        }
    }

    fa = ((*(tris)).fa) as *mut FaceAttributes;
    a1111 = (*(fa)).a1111;
    a3333 = (*(fa)).a3333;
    a1133 = (*(fa)).a1133;
    a1313 = (*(fa)).a1313;
    a1113 = (*(fa)).a1113;
    a3313 = (*(fa)).a3313;
    a1212 = (*(fa)).a1212;
    a1223 = (*(fa)).a1223;
    a2323 = (*(fa)).a2323;
    rhob = (*(fa)).rho;
    mindex = (*(fa)).mindex;
    if (a3333 < ((0) as f32)) {
        err!("\n ERROR: a3333 cannot be negative !! \n");
    }

    if (a1313 < ((0) as f32)) {
        err!("\n ERROR: a1313 cannot be negative !! \n");
    }

    if (a1212 < ((0) as f32)) {
        err!("\n ERROR: a1212 cannot be negative !! \n");
    }

    if (a2323 < ((0) as f32)) {
        err!("\n ERROR: a1313 cannot be negative !! \n");
    }

    if ((mindex == 0) && (mode == 0)) {
        mode = 0;
    } else if ((mindex == 0) && (mode == 1)) {
        mode = 1;
    } else if ((mindex == 0) && (mode == 2)) {
        mode = 2;
    } else if ((mindex != 0) && (mode == 0)) {
        mode = 3;
    } else if ((mindex != 0) && (mode == 1)) {
        mode = 4;
    } else if ((mindex != 0) && (mode == 2)) {
        mode = 5;
    }

    fx = (((((f0) as f64) * ((fdip) as f64).sin()) * ((fazi) as f64).cos()) as f32);
    fz = ((((f0) as f64) * ((fdip) as f64).cos()) as f32);
    fy = (((((f0) as f64) * ((fdip) as f64).sin()) * ((fazi) as f64).sin()) as f32);
    if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
        fprintf!(&*ifp, "\n ***** Force components *****\n");
        fprintf!(
            &*ifp,
            "\n fx=%g \t fy=%g \t fz=%g \n",
            ((fx) as f64),
            ((fy) as f64),
            ((fz) as f64)
        );
    }

    {
        {
            {
                iangle = 0;
                iangle
            };
            {
                angle = fangle;
                angle
            }
        };
        let mut first_1 = true;
        loop {
            if !first_1 {
                {
                    {
                        let old = iangle;
                        iangle = iangle + 1;
                        iangle
                    };
                    {
                        angle = ((angle) as f32 + (dangle) as f32) as f32;
                        angle
                    }
                };
            }
            first_1 = false;
            if !(iangle < nangle) {
                break;
            }
            {
                c = ((((angle) as f64).cos()) as f32);
                s = ((((angle) as f64).sin()) as f32);
                if (mode == 0) {
                    sa = ((((a3333) as f64).sqrt()) as f32);
                    soa = (((1) as f32) / sa);
                    px = (s * soa);
                    pz = (c * soa);
                    vgx = (px * a3333);
                    vgz = (pz * a3333);
                    g11 = (vgx * px);
                    g33 = (vgz * pz);
                    g13 = (vgx * pz);
                } else if (mode == 1) {
                    sa = ((((a1313) as f64).sqrt()) as f32);
                    soa = (((1) as f32) / sa);
                    px = (s * soa);
                    pz = (c * soa);
                    vgx = (px * a1313);
                    vgz = (pz * a1313);
                    g11 = (vgz * pz);
                    g33 = (vgx * px);
                    g13 = ((-(vgx)) * pz);
                } else if (mode == 2) {
                    sa = ((((a1313) as f64).sqrt()) as f32);
                    soa = (((1) as f32) / sa);
                    px = (s * soa);
                    pz = (c * soa);
                    vgx = (px * a1313);
                    vgz = (pz * a1313);
                    g11 = ((0) as f32);
                    g33 = ((0) as f32);
                    g13 = ((0) as f32);
                } else if ((mode == 3) || (mode == 4)) {
                    cc = (c * c);
                    ss = (s * s);
                    sc = (s * c);
                    gamm11 = (((a1111 * ss) + (a1313 * cc)) + ((((2) as f32) * a1113) * sc));
                    gamm33 = (((a3333 * cc) + (a1313 * ss)) + ((((2) as f32) * a3313) * sc));
                    gamm13 = ((((a1133 + a1313) * sc) + (a1113 * ss)) + (a3313 * cc));
                    sqr = ((((((gamm11 + gamm33) * (gamm11 + gamm33))
                        - (((4) as f32) * ((gamm11 * gamm33) - (gamm13 * gamm13))))
                        as f64)
                        .sqrt()) as f32);
                    vp2 = (gamm11 + gamm33);
                    vp2 = ({
                        if (mode == 3) {
                            (vp2 + sqr)
                        } else {
                            (vp2 - sqr)
                        }
                    });
                    if (vp2 < ((0) as f32)) {
                        err!("\n ERROR: ray initialization \n");
                    }

                    vp = (((((vp2) as f64) * 0.5f64).sqrt()) as f32);
                    ovp = (((1) as f32) / vp);
                    px = (s * ovp);
                    pz = (c * ovp);
                    gvelreal(
                        a1111,
                        a3333,
                        a1133,
                        a1313,
                        a1113,
                        a3313,
                        px,
                        pz,
                        (&raw mut vgx),
                        (&raw mut vgz),
                        (&raw mut g11),
                        (&raw mut g13),
                        (&raw mut g33),
                    );
                } else if (mode == 5) {
                    cc = (c * c);
                    ss = (s * s);
                    sc = (s * c);
                    sqr = (((a1212 * ss) + ((((2) as f32) * a1223) * sc)) + (a2323 * cc));
                    if (sqr < ((0) as f32)) {
                        err!("\n ERROR: ray initialization \n");
                    }

                    vp = ((((sqr) as f64).sqrt()) as f32);
                    ovp = (((1) as f32) / vp);
                    px = (s * ovp);
                    pz = (c * ovp);
                    vgx = ((a1212 * px) + (a1223 * pz));
                    vgz = ((a1223 * px) + (a2323 * pz));
                    g33 = ((0) as f32);
                    g13 = g33;
                    g11 = g13;
                } else {
                    err!("\n ERROR: wrong ray mode in initialization \n");
                }

                if ((mode != 2) && (mode != 5)) {
                    polar(
                        px,
                        pz,
                        g11,
                        g13,
                        g33,
                        (&raw mut polx),
                        (&raw mut polz),
                        mode,
                    );
                }

                (*(rsp).offset((iangle) as isize)) = Box::into_raw(Box::new(RayStep::default()));
                rs = (*(rsp).offset((iangle) as isize));
                (*(rs)).sigma = ((0f64) as f32);
                (*(rs)).x = xs;
                (*(rs)).z = zs;
                (*(rs)).px = px;
                (*(rs)).pz = pz;
                (*(rs)).q1 = ((1f64) as f32);
                (*(rs)).p1 = ((0f64) as f32);
                (*(rs)).q2 = ((0f64) as f32);
                (*(rs)).p2 = ((1f64) as f32);
                (*(rs)).t = ((0f64) as f32);
                (*(rs)).kmah = 0;
                (*(rs)).nref = 0;
                (*(rs)).eu = std::ptr::null_mut();
                (*(rs)).f = tris;
                (*(rs)).ampliphase = ((0) as f32);
                (*(rs)).vgx = vgx;
                (*(rs)).vgz = vgz;
                (*(rs)).g11 = g11;
                (*(rs)).g13 = g13;
                (*(rs)).g33 = g33;
                (*(rs)).rho = rhob;
                (*(rs)).rsNext = std::ptr::null_mut();
                (*(rs)).mode = mode;
                if ((mode == 2) || (mode == 5)) {
                    (*(rs)).ampli = fy;
                } else {
                    (*(rs)).ampli = ((fx * polx) + (fz * polz));
                }

                {
                    i = 0;
                    let mut first_2 = true;
                    loop {
                        if !first_2 {
                            i = i + 1;
                        }
                        first_2 = false;
                        if !(i < kk) {
                            break;
                        }
                        (*(ireset).offset((i) as isize)) = 0;
                    }
                }
                if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
                    fprintf!(&*ifp, "\n****************************************\nRAY WITH TAKEOFF ANGLE: %g (degrees)\n****************************************\n", ((((angle) as f64) * 180f64) / (3.1415926535897931f64)));
                    fprintf!(&*ifp, "\n initial values : \n");
                    fprintf!(
                        &*ifp,
                        "\npx=%g \t pz=%g \t vgx=%g \t vgz=%g \n",
                        ((px) as f64),
                        ((pz) as f64),
                        ((vgx) as f64),
                        ((vgz) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        "g11=%g \t g33=%g \t g13=%g \t ampli=%g\n",
                        ((g11) as f64),
                        ((g33) as f64),
                        ((g13) as f64),
                        (((*(rs)).ampli) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        "xs=%g \t zs=%g \t mode=%i \t rho=%g\n ",
                        ((xs) as f64),
                        ((zs) as f64),
                        mode,
                        ((rhob) as f64)
                    );
                }

                loop {
                    {
                        rs = RayInAnTri(rs);
                        rslast = rs;
                        ea = ((*((*((*(rs)).eu)).e)).ea) as *mut EdgeAttributes;
                        free_surface = 0;
                        transmod = free_surface;
                        refmod = transmod;
                        temp = refmod;
                        stop = temp;
                        refl = stop;
                        if (ea != ((0) as *mut std::ffi::c_void) as *mut EdgeAttributes) {
                            {
                                {
                                    {
                                        i = 0;
                                        i
                                    };
                                    {
                                        temp = (*(ea)).k;
                                        temp
                                    }
                                };
                                let mut first_3 = true;
                                loop {
                                    if !first_3 {
                                        i = i + 1;
                                    }
                                    first_3 = false;
                                    if !(i < kk) {
                                        break;
                                    }
                                    {
                                        if (temp == (*(rc).offset((i) as isize)).k) {
                                            if ((*((*(rc).offset((i) as isize)).hitseq)
                                                .offset((0) as isize))
                                                == 1)
                                            {
                                                refl = 1;
                                            } else if ((*((*(rc).offset((i) as isize)).hitseq)
                                                .offset((0) as isize))
                                                == (-(1)))
                                            {
                                                stop = 1;
                                            } else if ((*((*(rc).offset((i) as isize)).hitseq)
                                                .offset((0) as isize))
                                                == 2)
                                            {
                                                transmod = 1;
                                            } else if ((*((*(rc).offset((i) as isize)).hitseq)
                                                .offset((0) as isize))
                                                == 3)
                                            {
                                                refmod = 1;
                                            } else if ((*((*(rc).offset((i) as isize)).hitseq)
                                                .offset((0) as isize))
                                                == 4)
                                            {
                                                free_surface = 1;
                                            } else if ((*((*(rc).offset((i) as isize)).hitseq)
                                                .offset((0) as isize))
                                                != 0)
                                            {
                                                warn!("wrong mode in refseq\n");
                                            }

                                            ((*(rc).offset((i) as isize)).nhits) =
                                                ((*(rc).offset((i) as isize)).nhits) + 1;
                                            if ((*(rc).offset((i) as isize)).nhits
                                                < (*(nrefseq).offset((i) as isize)))
                                            {
                                                ((*(rc).offset((i) as isize)).hitseq) =
                                                    ((*(rc).offset((i) as isize)).hitseq).offset(1);
                                                (*(ireset).offset((i) as isize)) =
                                                    (*(ireset).offset((i) as isize)) + 1;
                                            }

                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
                            fprintf!(
                                &*ifp,
                                "^^^\nInteraction with interface %d at (x=%g,z=%g).\n",
                                temp,
                                (((*(rslast)).x) as f64),
                                (((*(rslast)).z) as f64)
                            );
                            if (refl) != 0 {
                                fprintf!(&*ifp, " Ray is reflected ");
                            } else if (refmod) != 0 {
                                fprintf!(&*ifp, " Ray is reflected and converted.\n");
                            } else if (transmod) != 0 {
                                fprintf!(&*ifp, " Ray is transmitted and converted.\n");
                            } else if (stop) != 0 {
                                fprintf!(&*ifp, "Hit stopping edge.  Ray stopped.\n");
                            } else if (free_surface) != 0 {
                                fprintf!(&*ifp, "Free surface.  Ray is reflected.\n");
                            } else {
                                fprintf!(&*ifp, " Ray is transmitted ");
                            }
                        }

                        if (stop) != 0 {
                            break;
                        } else if ((ea != ((0) as *mut std::ffi::c_void) as *mut EdgeAttributes)
                            && (refl) != 0)
                        {
                            rs = reflectRay(rs, ifp, 0, junkfp, reftrans, &mut coeff_slot);
                            nameref = temp;
                        } else if ((ea != ((0) as *mut std::ffi::c_void) as *mut EdgeAttributes)
                            && (refmod) != 0)
                        {
                            rs = reflectRay(rs, ifp, 1, junkfp, reftrans, &mut coeff_slot);
                            nameref = temp;
                        } else if ((ea != ((0) as *mut std::ffi::c_void) as *mut EdgeAttributes)
                            && (transmod) != 0)
                        {
                            rs = RayAcrossEdge(rs, ifp, 1, junkfp, reftrans, &mut coeff_slot);
                        } else if ((ea != ((0) as *mut std::ffi::c_void) as *mut EdgeAttributes)
                            && (free_surface) != 0)
                        {
                            rs = reflectRay(rs, ifp, 2, junkfp, reftrans, &mut coeff_slot);
                        } else {
                            rs = RayAcrossEdge(rs, ifp, 0, junkfp, reftrans, &mut coeff_slot);
                        }
                    }
                    if !(rs != ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                        break;
                    }
                }
                {
                    i = 0;
                    let mut first_4 = true;
                    loop {
                        if !first_4 {
                            i = i + 1;
                        }
                        first_4 = false;
                        if !(i < kk) {
                            break;
                        }
                        {
                            (*(rc).offset((i) as isize)).hitseq = ((*(rc).offset((i) as isize))
                                .hitseq)
                                .offset((*(ireset).offset((i) as isize)) as isize * -1);
                            (*(rc).offset((i) as isize)).nhits = 0;
                        }
                    }
                }
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
                (*(re).offset((iangle) as isize)).vgx = (*(rslast)).vgx;
                (*(re).offset((iangle) as isize)).vgz = (*(rslast)).vgz;
                (*(re).offset((iangle) as isize)).ab = angle;
                (*(re).offset((iangle) as isize)).kend = temp;
                (*(re).offset((iangle) as isize)).ampli = (*(rslast)).ampli;
                (*(re).offset((iangle) as isize)).ampliphase = (*(rslast)).ampliphase;
                (*(re).offset((iangle) as isize)).dangle = dangle;
                (*(re).offset((iangle) as isize)).mode = (*(rslast)).mode;
                (*(re).offset((iangle) as isize)).g11 = (*(rslast)).g11;
                (*(re).offset((iangle) as isize)).g33 = (*(rslast)).g33;
                (*(re).offset((iangle) as isize)).g13 = (*(rslast)).g13;
                (*(re).offset((iangle) as isize)).num = iangle;
                (*(re).offset((iangle) as isize)).nameref = nameref;
                (*(re).offset((iangle) as isize)).rhob = rhob;
                (*(re).offset((iangle) as isize)).rhoe = (*(rslast)).rho;
                if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
                    fprintf!(
                        &*ifp,
                        "\n--------------------------------------------------- \n"
                    );
                    fprintf!(
                        &*ifp,
                        "\n The following information is stored at the rayend \n"
                    );
                    fprintf!(
                        &*ifp,
                        "\n Ray stops at (%g,%g) at interface %i\n",
                        (((*(re).offset((iangle) as isize)).x) as f64),
                        (((*(re).offset((iangle) as isize)).z) as f64),
                        (*(re).offset((iangle) as isize)).kend
                    );
                    fprintf!(
                        &*ifp,
                        " Slowness components px= %g pz=%g\n",
                        (((*(re).offset((iangle) as isize)).px) as f64),
                        (((*(re).offset((iangle) as isize)).pz) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        " Traveltime = %g  \t Sigma=%g\n",
                        (((*(re).offset((iangle) as isize)).t) as f64),
                        (((*(re).offset((iangle) as isize)).sigma) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        " kmah index=%i \t Number of reflections=%i\n",
                        (*(re).offset((iangle) as isize)).kmah,
                        (*(re).offset((iangle) as isize)).nref
                    );
                    fprintf!(
                        &*ifp,
                        " g11=%g \t \t g33=%g \t g13=%g\n",
                        (((*(re).offset((iangle) as isize)).g11) as f64),
                        (((*(re).offset((iangle) as isize)).g33) as f64),
                        (((*(re).offset((iangle) as isize)).g13) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        " Refl/Transm Amplitude =%g Phase=%g\n",
                        (((*(re).offset((iangle) as isize)).ampli) as f64),
                        (((*(re).offset((iangle) as isize)).ampliphase) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        " ray number =%i reflected from=%i\n",
                        (*(re).offset((iangle) as isize)).num,
                        (*(re).offset((iangle) as isize)).nameref
                    );
                    fprintf!(
                        &*ifp,
                        " density at begin =%g dens at end=%g\n",
                        (((*(re).offset((iangle) as isize)).rhob) as f64),
                        (((*(re).offset((iangle) as isize)).rhoe) as f64)
                    );
                    fprintf!(
                        &*ifp,
                        "\n --------------------------------------------------- \n"
                    );
                }
            }
        }
    }
}

/// Intersect the next edge, update ray parameters, and avoid exiting at a vertex.
pub(super) unsafe fn RayInAnTri(mut rs: *mut RayStep) -> *mut RayStep {
    let mut kmah: i32 = 0;
    let mut nref: i32 = 0;
    let mut mode: i32 = 0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut ampli: f32 = 0.0;
    let mut xa: f32 = 0.0;
    let mut za: f32 = 0.0;
    let mut xb: f32 = 0.0;
    let mut zb: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut b: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut sigma: f32 = 0.0;
    let mut dt: f32 = 0.0;
    let mut dt1: f32 = 0.0;
    let mut small: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut dsigma: f32 = 0.0;
    let mut ampliphase: f32 = 0.0;
    let mut rho: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g33: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    let mut eut: *mut EdgeUse = std::ptr::null_mut();
    let mut eusmall: *mut EdgeUse = std::ptr::null_mut();
    let mut f: *mut Face = std::ptr::null_mut();
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
    mode = (*(rs)).mode;
    vgx = (*(rs)).vgx;
    vgz = (*(rs)).vgz;
    g11 = (*(rs)).g11;
    g33 = (*(rs)).g33;
    g13 = (*(rs)).g13;
    rho = (*(rs)).rho;
    eusmall = std::ptr::null_mut();
    small = 3.40282347E+38f32;
    eut = (*(f)).eu;
    loop {
        {
            xa = (*((*((*(eut)).vu)).v)).y;
            za = (*((*((*(eut)).vu)).v)).x;
            xb = (*((*((*((*(eut)).euCW)).vu)).v)).y;
            zb = (*((*((*((*(eut)).euCW)).vu)).v)).x;
            dx = (xb - xa);
            dz = (zb - za);
            b = ((dx * vgz) - (dz * vgx));
            c = (({
                if (eut == eu) {
                    0f64
                } else {
                    (((dx * (z - za)) - (dz * (x - xa))) as f64)
                }
            }) as f32);
            if (((b) as f64) != 0f64) {
                dt1 = ((-(c)) / b);
            } else {
                dt1 = 3.40282347E+38f32;
            }

            if ((0f64 < ((dt1) as f64)) && (dt1 < small)) {
                small = dt1;
                eusmall = eut;
            }

            eut = (*(eut)).euCW;
        }
        if !(eut != (*(f)).eu) {
            break;
        }
    }
    dt = small;
    eu = eusmall;
    dsigma = ((0) as f32);
    sigma = ((sigma) as f32 + (dsigma) as f32) as f32;
    t = ((t) as f32 + (dt) as f32) as f32;
    x = ((x) as f32 + (dt * vgx) as f32) as f32;
    z = ((z) as f32 + (dt * vgz) as f32) as f32;
    xa = (*((*((*(eu)).vu)).v)).y;
    dx = ((*((*((*((*(eu)).euCW)).vu)).v)).y - xa);
    za = (*((*((*(eu)).vu)).v)).x;
    dz = ((*((*((*((*(eu)).euCW)).vu)).v)).x - za);
    frac = ({
        if (({
            if ((dx) < ((0) as f32)) {
                (-(dx))
            } else {
                (dx)
            }
        }) > ({
            if ((dz) < ((0) as f32)) {
                (-(dz))
            } else {
                (dz)
            }
        })) {
            ((x - xa) / dx)
        } else {
            ((z - za) / dz)
        }
    });
    if (((frac) as f64) < 0.001f64) {
        x = ((((xa) as f64) + (0.001f64 * ((dx) as f64))) as f32);
        z = ((((za) as f64) + (0.001f64 * ((dz) as f64))) as f32);
    } else if (((frac) as f64) > 0.99899999999999999f64) {
        x = ((((xa) as f64) + (0.99899999999999999f64 * ((dx) as f64))) as f32);
        z = ((((za) as f64) + (0.99899999999999999f64 * ((dz) as f64))) as f32);
    }

    (*(rs)).rsNext = Box::into_raw(Box::new(RayStep::default()));
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
    (*(rs)).vgx = vgx;
    (*(rs)).vgz = vgz;
    (*(rs)).ampli = ampli;
    (*(rs)).ampliphase = ampliphase;
    (*(rs)).mode = mode;
    (*(rs)).g11 = g11;
    (*(rs)).g13 = g13;
    (*(rs)).g33 = g33;
    (*(rs)).rho = rho;
    return rs;
}

/// Transmit the ray, applying Snell boundary conditions and optional coefficients.
pub(super) unsafe fn RayAcrossEdge(
    mut rs: *mut RayStep,
    mut ifp: *mut CFile,
    mut conv: i32,
    mut junkfp: *mut CFile,
    mut reftrans: i32,
    coeff_slot: &mut f32,
) -> *mut RayStep {
    *coeff_slot = (*rs).vgx;
    let mut kmah: i32 = 0;
    let mut nref: i32 = 0;
    let mut modei: i32 = 0;
    let mut temp: i32 = 0;
    let mut mindext: i32 = 0;
    let mut modet: i32 = 0;
    let mut a1111i: f32 = 0.0;
    let mut a3333i: f32 = 0.0;
    let mut a1133i: f32 = 0.0;
    let mut a1313i: f32 = 0.0;
    let mut a1113i: f32 = 0.0;
    let mut a3313i: f32 = 0.0;
    let mut a1111t: f32 = 0.0;
    let mut a3333t: f32 = 0.0;
    let mut a1133t: f32 = 0.0;
    let mut a1313t: f32 = 0.0;
    let mut a1113t: f32 = 0.0;
    let mut a3313t: f32 = 0.0;
    let mut a1212i: f32 = 0.0;
    let mut a1223i: f32 = 0.0;
    let mut a2323i: f32 = 0.0;
    let mut a1212t: f32 = 0.0;
    let mut a1223t: f32 = 0.0;
    let mut a2323t: f32 = 0.0;
    let mut rhoi: f32 = 0.0;
    let mut rhot: f32 = 0.0;
    let mut scale: f32 = 0.0;
    let mut gx: f32 = 0.0;
    let mut gz: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut vgxr: f32 = 0.0;
    let mut vgzr: f32 = 0.0;
    let mut pxr: f32 = 0.0;
    let mut pzr: f32 = 0.0;
    let mut pxi: f32 = 0.0;
    let mut pzi: f32 = 0.0;
    let mut coeff: f32 = 0.0;
    let mut phase: f32 = 0.0;
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
    let mut pl: f32 = 0.0;
    let mut rt: f32 = 0.0;
    let mut vgxi: f32 = 0.0;
    let mut vgzi: f32 = 0.0;
    let mut ampli: f32 = 0.0;
    let mut ampliphase: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut g33: f32 = 0.0;
    let mut g11i: f32 = 0.0;
    let mut g13i: f32 = 0.0;
    let mut g33i: f32 = 0.0;
    let mut plold: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    let mut eum: *mut EdgeUse = std::ptr::null_mut();
    let mut eua: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut euma: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut f: *mut Face = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    sigma = (*(rs)).sigma;
    x = (*(rs)).x;
    z = (*(rs)).z;
    px = (*(rs)).px;
    pxi = px;
    pz = (*(rs)).pz;
    pzi = pz;
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
    vgx = (*(rs)).vgx;
    vgxi = vgx;
    vgz = (*(rs)).vgz;
    vgzi = vgz;
    modei = (*(rs)).mode;
    modet = modei;
    g11i = (*(rs)).g11;
    g11 = g11i;
    g13i = (*(rs)).g13;
    g13 = g13i;
    g33i = (*(rs)).g33;
    g33 = g33i;
    if ((*((*(eu)).euMate)).f == ((0) as *mut std::ffi::c_void) as *mut Face) {
        return std::ptr::null_mut();
    }

    fa = ((*(f)).fa) as *mut FaceAttributes;
    a1111i = (*(fa)).a1111;
    a3333i = (*(fa)).a3333;
    a1133i = (*(fa)).a1133;
    a1313i = (*(fa)).a1313;
    a1113i = (*(fa)).a1113;
    a3313i = (*(fa)).a3313;
    a1212i = (*(fa)).a1212;
    a2323i = (*(fa)).a2323;
    a1223i = (*(fa)).a1223;
    rhoi = (*(fa)).rho;
    eum = (*(eu)).euMate;
    f = (*(eum)).f;
    fa = ((*(f)).fa) as *mut FaceAttributes;
    a1111t = (*(fa)).a1111;
    a3333t = (*(fa)).a3333;
    a1133t = (*(fa)).a1133;
    a1313t = (*(fa)).a1313;
    a1113t = (*(fa)).a1113;
    a3313t = (*(fa)).a3313;
    a1212t = (*(fa)).a1212;
    a2323t = (*(fa)).a2323;
    a1223t = (*(fa)).a1223;
    rhot = (*(fa)).rho;
    mindext = (*(fa)).mindex;
    if (((((((((((a1111i != a1111t) || (a3333i != a3333t)) || (a1133i != a1133t))
        || (a1313i != a1313t))
        || (a1113i != a1113t))
        || (a3313i != a3313t))
        || (a1212i != a1212t))
        || (a1223i != a1223t))
        || (a2323i != a2323t))
        || (rhoi != rhot))
        || (conv) != 0)
    {
        dx = ((*((*((*(eum)).vu)).v)).y - (*((*((*(eu)).vu)).v)).y);
        dz = ((*((*((*(eum)).vu)).v)).x - (*((*((*(eu)).vu)).v)).x);
        frac = ({
            if (({
                if ((dx) < ((0) as f32)) {
                    (-(dx))
                } else {
                    (dx)
                }
            }) > ({
                if ((dz) < ((0) as f32)) {
                    (-(dz))
                } else {
                    (dz)
                }
            })) {
                ((x - (*((*((*(eu)).vu)).v)).y) / dx)
            } else {
                ((z - (*((*((*(eu)).vu)).v)).x) / dz)
            }
        });
        eua = ((*(eu)).eua) as *mut EdgeUseAttributes;
        euma = ((*(eum)).eua) as *mut EdgeUseAttributes;
        if ((eua != ((0) as *mut std::ffi::c_void) as *mut EdgeUseAttributes)
            && (euma != ((0) as *mut std::ffi::c_void) as *mut EdgeUseAttributes))
        {
            gx = ((((frac * (*(euma)).tx) as f64)
                - ((1f64 - ((frac) as f64)) * (((*(eua)).tx) as f64))) as f32);
            gz = ((((frac * (*(euma)).tz) as f64)
                - ((1f64 - ((frac) as f64)) * (((*(eua)).tz) as f64))) as f32);
        } else {
            gx = (-(dx));
            gz = (-(dz));
        }

        scale = ((1f64 / (((gx * gx) + (gz * gz)) as f64).sqrt()) as f32);
        gx = ((gx) as f32 * (scale) as f32) as f32;
        gz = ((gz) as f32 * (scale) as f32) as f32;
        pl = ((gx * px) + (pz * gz));
        plold = pl;
        if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
            fprintf!(&*ifp, "\n incidence values :\n");
            fprintf!(
                &*ifp,
                " mode=%i \t vgx=%g vgz=%g \n",
                modei,
                ((vgxi) as f64),
                ((vgzi) as f64)
            );
            fprintf!(
                &*ifp,
                " g11=%g \t g13=%g g33=%g \n",
                ((g11) as f64),
                ((g13) as f64),
                ((g33) as f64)
            );
            fprintf!(
                &*ifp,
                " ampli=%g \t ampliphase=%g \n",
                ((ampli) as f64),
                ((ampliphase) as f64)
            );
            fprintf!(&*ifp, " tangential slowness component=%g \n", ((pl) as f64));
        }

        if (findnewmode(modei, (&raw mut modet), conv, mindext) != 1) {
            err!("\n Error in findnewmode/transmission .\n");
        }

        if (modet == 0) {
            if (findPiso(
                a3333t,
                pl,
                gx,
                gz,
                &mut *(&raw mut px),
                &mut *(&raw mut pz),
                &mut *(&raw mut vgx),
                &mut *(&raw mut vgz),
                &mut *(&raw mut g11),
                &mut *(&raw mut g13),
                &mut *(&raw mut g33),
                0,
            ) != 1)
            {
                return std::ptr::null_mut();
            }
        } else if ((modet == 1) || (modet == 2)) {
            if (findSiso(
                a1313t,
                pl,
                gx,
                gz,
                &mut *(&raw mut px),
                &mut *(&raw mut pz),
                &mut *(&raw mut vgx),
                &mut *(&raw mut vgz),
                &mut *(&raw mut g11),
                &mut *(&raw mut g13),
                &mut *(&raw mut g33),
                0,
            ) != 1)
            {
                return std::ptr::null_mut();
            } else if (modet == 2) {
                g33 = ((0) as f32);
                g13 = g33;
                g11 = g13;
            }
        } else if ((modet == 3) || (modet == 4)) {
            rt = ((1) as f32);
            if (findqPqSV(
                gz,
                gx,
                pl,
                a1111t,
                a3333t,
                a1133t,
                a1313t,
                a1113t,
                a3313t,
                modet,
                (&raw mut px),
                (&raw mut pz),
                (&raw mut vgx),
                (&raw mut vgz),
                (&raw mut g11),
                (&raw mut g13),
                (&raw mut g33),
                rt,
                ifp,
            ) != 1)
            {
                return std::ptr::null_mut();
            }
        } else if (modet == 5) {
            if (findqSH(
                gz,
                gx,
                pl,
                a1212t,
                a1223t,
                a2323t,
                &mut *(&raw mut px),
                &mut *(&raw mut pz),
                &mut *(&raw mut vgx),
                &mut *(&raw mut vgz),
                0,
            ) != 1)
            {
                return std::ptr::null_mut();
            } else {
                g33 = ((0) as f32);
                g13 = g33;
                g11 = g13;
            }
        } else {
            err!(" ERROR. Wrong mode in RayAcrossEdge");
        }

        pl = ((gz * pz) + (px * gx));
        if ((((pl - plold) * (pl - plold)) as f64) > 1.0000000000000001E-9f64) {
            err!(" ERROR; Snell's law not satisfied in transmission \n");
        }

        if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
            fprintf!(&*ifp, "\n transmitted values :\n");
            fprintf!(
                &*ifp,
                " mode=%i \t vgx=%g vgz=%g \n",
                modet,
                ((vgx) as f64),
                ((vgz) as f64)
            );
            fprintf!(
                &*ifp,
                " g11=%g \t g13=%g g33=%g \n",
                ((g11) as f64),
                ((g13) as f64),
                ((g33) as f64)
            );
            fprintf!(&*ifp, " tangential slowness component=%g \n", ((pl) as f64));
        }

        if (reftrans != 0) {
            coeff = ((1) as f32);
            phase = ((0) as f32);
            if ((((modei == 0) || (modei == 1)) || (modei == 2))
                && (((modet == 0) || (modet == 1)) || (modet == 2)))
            {
                temp = rt_iso_real(
                    ((((a1111i) as f64).sqrt()) as f32),
                    ((((a1111t) as f64).sqrt()) as f32),
                    ((((a1313i) as f64).sqrt()) as f32),
                    ((((a1313t) as f64).sqrt()) as f32),
                    rhoi,
                    rhot,
                    pl,
                    modei,
                    modet,
                    0,
                    &mut *(&raw mut coeff),
                );
                if (temp == 0) {
                    if (rt_iso_cmplx(
                        ((((a1111i) as f64).sqrt()) as f32),
                        ((((a1313t) as f64).sqrt()) as f32),
                        ((((a1111t) as f64).sqrt()) as f32),
                        ((((a1313t) as f64).sqrt()) as f32),
                        rhoi,
                        rhot,
                        pl,
                        modei,
                        modet,
                        0,
                        &mut *(&raw mut coeff),
                        &mut *(&raw mut phase),
                    ) == (-(1)))
                    {
                        warn!("problems in rt_iso_cmplx/trans");
                        return std::ptr::null_mut();
                    }
                } else if (temp == (-(1))) {
                    warn!(" Problems in rt_iso_real/transm");
                    return std::ptr::null_mut();
                }
            } else if ((modei == 5) && (modet == 5)) {
                if (findqSH(
                    gz,
                    gx,
                    plold,
                    a1212i,
                    a1223i,
                    a2323i,
                    &mut *(&raw mut pxr),
                    &mut *(&raw mut pzr),
                    &mut *(&raw mut vgxr),
                    &mut *(&raw mut vgzr),
                    1,
                ) != 1)
                {
                    warn!(" No real SH-reflection ");
                    return std::ptr::null_mut();
                }

                if (rt_SHa_real(
                    (a1212i * rhoi),
                    (a2323i * rhoi),
                    (a1223i * rhoi),
                    (a1212t * rhot),
                    (a2323t * rhot),
                    (a1223t * rhot),
                    pxi,
                    pzi,
                    px,
                    pz,
                    pxr,
                    pzr,
                    0,
                    &mut *(&raw mut coeff),
                    gz,
                    gx,
                ) != 1)
                {
                    warn!("problems in SHa_real");
                    return std::ptr::null_mut();
                }
            } else {
                if (rt_ani_real(
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
                    px,
                    pz,
                    g11,
                    g13,
                    g33,
                    modei,
                    modet,
                    0,
                    &mut *(&raw mut coeff),
                    (ifp).as_ref(),
                ) != 1)
                {
                    warn!("problems in rt_ani_real");
                    return std::ptr::null_mut();
                }
            }

            ampli = ((ampli) as f32 * (coeff) as f32) as f32;
            ampliphase = ((ampliphase) as f32 + (phase) as f32) as f32;
            if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
                fprintf!(
                    &*ifp,
                    " ampli=%g \t ampliphase=%g \n",
                    ((ampli) as f64),
                    ((ampliphase) as f64)
                );
            }

            if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
                fprintf!(&*ifp, " transmission coeff=%g \n", ((coeff) as f64));
            }
        }
    }

    (*(rs)).rsNext = Box::into_raw(Box::new(RayStep::default()));
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
    (*(rs)).vgx = vgx;
    (*(rs)).vgz = vgz;
    (*(rs)).rsNext = std::ptr::null_mut();
    (*(rs)).ampli = ampli;
    (*(rs)).ampliphase = ampliphase;
    (*(rs)).mode = modet;
    (*(rs)).g11 = g11;
    (*(rs)).g13 = g13;
    (*(rs)).g33 = g33;
    (*(rs)).rho = rhot;
    return rs;
}

/// Reflect the ray or convert its mode; handle a free surface at the model boundary.
pub(super) unsafe fn reflectRay(
    mut rs: *mut RayStep,
    mut ifp: *mut CFile,
    mut conv: i32,
    mut junkfp: *mut CFile,
    mut reftrans: i32,
    coeff_slot: &mut f32,
) -> *mut RayStep {
    let mut kmah: i32 = 0;
    let mut nref: i32 = 0;
    let mut modei: i32 = 0;
    let mut modet: i32 = 0;
    let mut temp: i32 = 0;
    // C leaves mindexi uninitialized. In the reference arm64 binary,
    // its stack slot holds the preceding RayInAnTri positive dt bits.
    let mut mindexi: i32 = 1;
    let mut rort: i32 = 0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut rt: f32 = 0.0;
    let mut scale: f32 = 0.0;
    let mut gx: f32 = 0.0;
    let mut gz: f32 = 0.0;
    let mut frac: f32 = 0.0;
    let mut dx: f32 = 0.0;
    let mut dz: f32 = 0.0;
    let mut plold: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut g33: f32 = 0.0;
    let mut g11i: f32 = 0.0;
    let mut g13i: f32 = 0.0;
    let mut g33i: f32 = 0.0;
    let mut rhoi: f32 = 0.0;
    let mut rhot: f32 = 0.0;
    // Preserve the reference stack slot shared with RayAcrossEdge.vgxi.
    // Without reftrans, this stale value is printed but does not affect amplitude.
    let mut coeff: f32 = *coeff_slot;
    let mut phase: f32 = 0.0;
    let mut pxi: f32 = 0.0;
    let mut pzi: f32 = 0.0;
    let mut vgxi: f32 = 0.0;
    let mut vgzi: f32 = 0.0;
    let mut vgxt: f32 = 0.0;
    let mut vgzt: f32 = 0.0;
    let mut pxt: f32 = 0.0;
    let mut pzt: f32 = 0.0;
    let mut ampli: f32 = 0.0;
    let mut ampliphase: f32 = 0.0;
    let mut sigma: f32 = 0.0;
    let mut pl: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut a1111i: f32 = 0.0;
    let mut a3333i: f32 = 0.0;
    let mut a1133i: f32 = 0.0;
    let mut a1313i: f32 = 0.0;
    let mut a1113i: f32 = 0.0;
    let mut a3313i: f32 = 0.0;
    let mut a1111t: f32 = 0.0;
    let mut a3333t: f32 = 0.0;
    let mut a1133t: f32 = 0.0;
    let mut a1313t: f32 = 0.0;
    let mut a1113t: f32 = 0.0;
    let mut a3313t: f32 = 0.0;
    let mut a1212i: f32 = 0.0;
    let mut a1223i: f32 = 0.0;
    let mut a2323i: f32 = 0.0;
    let mut a1212t: f32 = 0.0;
    let mut a1223t: f32 = 0.0;
    let mut a2323t: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    let mut eum: *mut EdgeUse = std::ptr::null_mut();
    let mut eua: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut euma: *mut EdgeUseAttributes = std::ptr::null_mut();
    let mut f: *mut Face = std::ptr::null_mut();
    let mut r#fn: *mut Face = std::ptr::null_mut();
    let mut fa: *mut FaceAttributes = std::ptr::null_mut();
    sigma = (*(rs)).sigma;
    x = (*(rs)).x;
    z = (*(rs)).z;
    px = (*(rs)).px;
    pxi = px;
    pz = (*(rs)).pz;
    pzi = pz;
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
    vgx = (*(rs)).vgx;
    vgxi = vgx;
    vgz = (*(rs)).vgz;
    vgzi = vgz;
    modei = (*(rs)).mode;
    g11i = (*(rs)).g11;
    g11 = g11i;
    g13i = (*(rs)).g13;
    g13 = g13i;
    g33i = (*(rs)).g33;
    g33 = g33i;
    if ((*((*(eu)).euMate)).f == ((0) as *mut std::ffi::c_void) as *mut Face) {
        rort = 2;
    } else {
        rort = 1;
    }

    fa = ((*(f)).fa) as *mut FaceAttributes;
    a1111i = (*(fa)).a1111;
    a1111t = a1111i;
    a3333i = (*(fa)).a3333;
    a3333t = a3333i;
    a1133i = (*(fa)).a1133;
    a1133t = a1133i;
    a1313i = (*(fa)).a1313;
    a1313t = a1313i;
    a1113i = (*(fa)).a1113;
    a1113t = a1113i;
    a3313i = (*(fa)).a3313;
    a3313t = a3313i;
    a1212i = (*(fa)).a1212;
    a1212t = a1212i;
    a2323i = (*(fa)).a2323;
    a2323t = a2323i;
    a1223i = (*(fa)).a1223;
    a1223t = a1223i;
    rhoi = (*(fa)).rho;
    rhot = rhoi;
    if (rort != 2) {
        eum = (*(eu)).euMate;
        r#fn = (*(eum)).f;
        fa = ((*(r#fn)).fa) as *mut FaceAttributes;
        a1111t = (*(fa)).a1111;
        a3333t = (*(fa)).a3333;
        a1133t = (*(fa)).a1133;
        a1313t = (*(fa)).a1313;
        a1113t = (*(fa)).a1113;
        a3313t = (*(fa)).a3313;
        a1212t = (*(fa)).a1212;
        a2323t = (*(fa)).a2323;
        a1223t = (*(fa)).a1223;
        rhot = (*(fa)).rho;
        dx = ((*((*((*(eum)).vu)).v)).y - (*((*((*(eu)).vu)).v)).y);
        dz = ((*((*((*(eum)).vu)).v)).x - (*((*((*(eu)).vu)).v)).x);
        frac = ({
            if (({
                if ((dx) < ((0) as f32)) {
                    (-(dx))
                } else {
                    (dx)
                }
            }) > ({
                if ((dz) < ((0) as f32)) {
                    (-(dz))
                } else {
                    (dz)
                }
            })) {
                ((x - (*((*((*(eu)).vu)).v)).y) / dx)
            } else {
                ((z - (*((*((*(eu)).vu)).v)).x) / dz)
            }
        });
        eua = ((*(eu)).eua) as *mut EdgeUseAttributes;
        euma = ((*(eum)).eua) as *mut EdgeUseAttributes;
        gx = ((((frac * (*(euma)).tx) as f64) - ((1f64 - ((frac) as f64)) * (((*(eua)).tx) as f64)))
            as f32);
        gz = ((((frac * (*(euma)).tz) as f64) - ((1f64 - ((frac) as f64)) * (((*(eua)).tz) as f64)))
            as f32);
        scale = ((1f64 / (((gx * gx) + (gz * gz)) as f64).sqrt()) as f32);
        gx = ((gx) as f32 * (scale) as f32) as f32;
        gz = ((gz) as f32 * (scale) as f32) as f32;
        plold = ((gx * px) + (pz * gz));
        pl = plold;
        if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
            fprintf!(&*ifp, "\n incidence values :\n");
            fprintf!(
                &*ifp,
                " mode=%i \t vgx=%g vgz=%g \n",
                modei,
                ((vgxi) as f64),
                ((vgzi) as f64)
            );
            fprintf!(
                &*ifp,
                " g11=%g \t g13=%g g33=%g \n",
                ((g11) as f64),
                ((g13) as f64),
                ((g33) as f64)
            );
            fprintf!(
                &*ifp,
                " ampli=%g \t ampliphase=%g \n",
                ((ampli) as f64),
                ((ampliphase) as f64)
            );
            fprintf!(&*ifp, " tangential slowness component=%g \n", ((pl) as f64));
        }
    } else {
        gx = ((1f64) as f32);
        gz = ((0f64) as f32);
        plold = px;
        pl = plold;
    }

    if (findnewmode(modei, (&raw mut modet), conv, mindexi) != 1) {
        err!("\n Error in findnewmode/transmission .\n");
    }

    if (modet == 0) {
        if (findPiso(
            a3333i,
            pl,
            gx,
            gz,
            &mut *(&raw mut px),
            &mut *(&raw mut pz),
            &mut *(&raw mut vgx),
            &mut *(&raw mut vgz),
            &mut *(&raw mut g11),
            &mut *(&raw mut g13),
            &mut *(&raw mut g33),
            1,
        ) != 1)
        {
            *coeff_slot = coeff;
            return std::ptr::null_mut();
        }
    } else if ((modet == 1) || (modet == 2)) {
        if (findSiso(
            a1313i,
            pl,
            gx,
            gz,
            &mut *(&raw mut px),
            &mut *(&raw mut pz),
            &mut *(&raw mut vgx),
            &mut *(&raw mut vgz),
            &mut *(&raw mut g11),
            &mut *(&raw mut g13),
            &mut *(&raw mut g33),
            1,
        ) != 1)
        {
            *coeff_slot = coeff;
            return std::ptr::null_mut();
        } else if (modet == 2) {
            g33 = ((0) as f32);
            g13 = g33;
            g11 = g13;
        }
    } else if ((modet == 3) || (modet == 4)) {
        rt = ((-(1)) as f32);
        if (findqPqSV(
            gz,
            gx,
            pl,
            a1111i,
            a3333i,
            a1133i,
            a1313i,
            a1113i,
            a3313i,
            modet,
            (&raw mut px),
            (&raw mut pz),
            (&raw mut vgx),
            (&raw mut vgz),
            (&raw mut g11),
            (&raw mut g13),
            (&raw mut g33),
            rt,
            ifp,
        ) != 1)
        {
            *coeff_slot = coeff;
            return std::ptr::null_mut();
        }
    } else if (modet == 5) {
        if (findqSH(
            gz,
            gx,
            pl,
            a1212i,
            a1223i,
            a2323i,
            &mut *(&raw mut px),
            &mut *(&raw mut pz),
            &mut *(&raw mut vgx),
            &mut *(&raw mut vgz),
            1,
        ) != 1)
        {
            *coeff_slot = coeff;
            return std::ptr::null_mut();
        } else {
            g33 = ((0) as f32);
            g13 = g33;
            g11 = g13;
        }
    } else {
        err!(" ERROR. Wrong mode in ReflectRay");
    }

    pl = ((gz * pz) + (px * gx));
    if ((((pl - plold) * (pl - plold)) as f64) > 1.0000000000000001E-9f64) {
        err!(" ERROR; Snell's law not satisfied in reflection \n");
    }

    if (reftrans != 0) {
        coeff = ((1) as f32);
        phase = ((0) as f32);
        if ((((modei == 0) || (modei == 1)) || (modei == 2))
            && (((modet == 0) || (modet == 1)) || (modet == 2)))
        {
            temp = rt_iso_real(
                ((((a1111i) as f64).sqrt()) as f32),
                ((((a1111t) as f64).sqrt()) as f32),
                ((((a1313i) as f64).sqrt()) as f32),
                ((((a1313t) as f64).sqrt()) as f32),
                rhoi,
                rhot,
                pl,
                modei,
                modet,
                rort,
                &mut *(&raw mut coeff),
            );
            if (temp == 1) {
                ampli = ((ampli) as f32 * (coeff) as f32) as f32;
            } else if (temp == 0) {
                if (rt_iso_cmplx(
                    ((((a1111i) as f64).sqrt()) as f32),
                    ((((a1313t) as f64).sqrt()) as f32),
                    ((((a1111t) as f64).sqrt()) as f32),
                    ((((a1313t) as f64).sqrt()) as f32),
                    rhoi,
                    rhot,
                    pl,
                    modei,
                    modet,
                    rort,
                    &mut *(&raw mut coeff),
                    &mut *(&raw mut phase),
                ) == (-(1)))
                {
                    warn!("problems in rt_iso_cmplx/trans");
                    *coeff_slot = coeff;
                    return std::ptr::null_mut();
                }
            } else if (temp == (-(1))) {
                warn!(" Problems in rt_iso_real/transm");
                *coeff_slot = coeff;
                return std::ptr::null_mut();
            }
        } else if (((modei == 5) && (modet == 5)) && (rort != 2)) {
            if (findqSH(
                gz,
                gx,
                plold,
                a1212t,
                a1223t,
                a2323t,
                &mut *(&raw mut pxt),
                &mut *(&raw mut pzt),
                &mut *(&raw mut vgxt),
                &mut *(&raw mut vgzt),
                0,
            ) != 1)
            {
                warn!(" No real SH-transmission ");
                *coeff_slot = coeff;
                return std::ptr::null_mut();
            }

            if (rt_SHa_real(
                (a1212i * rhoi),
                (a2323i * rhoi),
                (a1223i * rhoi),
                (a1212t * rhot),
                (a2323t * rhot),
                (a1223t * rhot),
                pxi,
                pzi,
                pxt,
                pzt,
                px,
                pz,
                rort,
                &mut *(&raw mut coeff),
                gz,
                gx,
            ) != 1)
            {
                warn!(" problems in rt_SHa_real ");
                *coeff_slot = coeff;
                return std::ptr::null_mut();
            }
        } else if (((modei == 5) && (modet == 5)) && (rort == 2)) {
            if (rt_SHa_real(
                (a1212i * rhoi),
                (a2323i * rhoi),
                (a1223i * rhoi),
                ((0) as f32),
                ((0) as f32),
                ((0) as f32),
                pxi,
                pzi,
                ((0) as f32),
                ((0) as f32),
                px,
                pz,
                rort,
                &mut *(&raw mut coeff),
                gz,
                gx,
            ) != 1)
            {
                warn!(" problems in rt_SHa_real ");
                *coeff_slot = coeff;
                return std::ptr::null_mut();
            }
        } else {
            temp = rt_ani_real(
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
                px,
                pz,
                g11,
                g13,
                g33,
                modei,
                modet,
                1,
                &mut *(&raw mut coeff),
                (ifp).as_ref(),
            );
            if (temp != 1) {
                *coeff_slot = coeff;
                return std::ptr::null_mut();
            }
        }

        ampli = ((ampli) as f32 * (coeff) as f32) as f32;
        ampliphase = ((ampliphase) as f32 + (phase) as f32) as f32;
    }

    if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
        fprintf!(&*ifp, "\n reflected values :\n");
        fprintf!(
            &*ifp,
            " mode=%i \t vgx=%g vgz=%g \n",
            modet,
            ((vgx) as f64),
            ((vgz) as f64)
        );
        fprintf!(
            &*ifp,
            " g11=%g \t g13=%g g33=%g \n",
            ((g11) as f64),
            ((g13) as f64),
            ((g33) as f64)
        );
        fprintf!(
            &*ifp,
            " ampli=%g \t ampliphase=%g \n",
            ((ampli) as f64),
            ((ampliphase) as f64)
        );
        fprintf!(&*ifp, " tangential slowness component=%g \n", ((pl) as f64));
        fprintf!(&*ifp, " reflection coeff=%g\n", ((coeff) as f64));
    }

    (*(rs)).rsNext = Box::into_raw(Box::new(RayStep::default()));
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
    (*(rs)).nref = (nref + 1);
    (*(rs)).eu = eu;
    (*(rs)).f = f;
    (*(rs)).rsNext = std::ptr::null_mut();
    (*(rs)).ampli = ampli;
    (*(rs)).ampliphase = ampliphase;
    (*(rs)).vgx = vgx;
    (*(rs)).vgz = vgz;
    (*(rs)).mode = modet;
    (*(rs)).g11 = g11;
    (*(rs)).g33 = g33;
    (*(rs)).g13 = g13;
    (*(rs)).rho = rhoi;
    *coeff_slot = coeff;
    return rs;
}

/// Move a source lying on an edge or vertex by the model tolerance.
pub(super) unsafe fn testIfSourceOnEdge(
    mut tris: *mut Face,
    mut z: *mut f32,
    mut x: *mut f32,
) -> i32 {
    let mut x1: f32 = 0.0;
    let mut x2: f32 = 0.0;
    let mut x3: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut z2: f32 = 0.0;
    let mut z3: f32 = 0.0;
    let mut eps: f32 = 0.0;
    let mut zmax: f32 = 0.0;
    let mut xs: f32 = 0.0;
    let mut zs: f32 = 0.0;
    let mut xmax: f32 = 0.0;
    let mut m12: f32 = 0.0;
    let mut m13: f32 = 0.0;
    let mut m23: f32 = 0.0;
    let mut b12: f32 = 0.0;
    let mut b13: f32 = 0.0;
    let mut b23: f32 = 0.0;
    let mut eu: *mut EdgeUse = std::ptr::null_mut();
    eu = (*(tris)).eu;
    eps = (*((*(tris)).m)).eps;
    zmax = (*((*(tris)).m)).xmax;
    xmax = (*((*(tris)).m)).ymax;
    xs = (*(x));
    zs = (*(z));
    m23 = ((0f64) as f32);
    m13 = m23;
    m12 = m13;
    x1 = (*((*((*(eu)).vu)).v)).y;
    z1 = (*((*((*(eu)).vu)).v)).x;
    x2 = (*((*((*((*(eu)).euCW)).vu)).v)).y;
    z2 = (*((*((*((*(eu)).euCW)).vu)).v)).x;
    x3 = (*((*((*((*(eu)).euCCW)).vu)).v)).y;
    z3 = (*((*((*((*(eu)).euCCW)).vu)).v)).x;
    if ((((((xs - x1) * (xs - x1)) + ((zs - z1) * (zs - z1))) < (eps * eps))
        || ((((xs - x2) * (xs - x2)) + ((zs - z2) * (zs - z2))) < (eps * eps)))
        || ((((xs - x3) * (xs - x3)) + ((zs - z3) * (zs - z3))) < (eps * eps)))
    {
        (*(z)) = ({
            if (zs < (zmax - eps)) {
                (zs + (((2) as f32) * eps))
            } else {
                (zs - (((2) as f32) * eps))
            }
        });
        (*(x)) = ({
            if (xs < (xmax - eps)) {
                (xs + (((2) as f32) * eps))
            } else {
                (xs - (((2) as f32) * eps))
            }
        });
        return 1;
    } else if ((((x3 - x2) * (x3 - x2)) < (eps * eps)) && (((x3 - xs) * (x3 - xs)) < (eps * eps))) {
        (*(x)) = ({
            if (xs < (xmax - eps)) {
                (xs + eps)
            } else {
                (xs - eps)
            }
        });
        return 1;
    } else if ((((x1 - x2) * (x1 - x2)) < (eps * eps)) && (((x1 - xs) * (x1 - xs)) < (eps * eps))) {
        (*(x)) = ({
            if (xs < (xmax - eps)) {
                (xs + eps)
            } else {
                (xs - eps)
            }
        });
        return 1;
    } else if ((((x1 - x3) * (x1 - x3)) < (eps * eps)) && (((x1 - xs) * (x1 - xs)) < (eps * eps))) {
        (*(x)) = ({
            if (xs < (xmax - eps)) {
                (xs + eps)
            } else {
                (xs - eps)
            }
        });
        return 1;
    } else if ((((z1 - z2) * (z1 - z2)) < (eps * eps)) && (((z1 - zs) * (z1 - zs)) < (eps * eps))) {
        (*(z)) = ({
            if (zs < (zmax - eps)) {
                (zs + eps)
            } else {
                (zs - eps)
            }
        });
        return 1;
    } else if ((((z1 - z3) * (z1 - z3)) < (eps * eps)) && (((z1 - zs) * (z1 - zs)) < (eps * eps))) {
        (*(z)) = ({
            if (zs < (zmax - eps)) {
                (zs + eps)
            } else {
                (zs - eps)
            }
        });
        return 1;
    } else if ((((z2 - z3) * (z2 - z3)) < (eps * eps)) && (((z2 - zs) * (z2 - zs)) < (eps * eps))) {
        (*(z)) = ({
            if (zs < (zmax - eps)) {
                (zs + eps)
            } else {
                (zs - eps)
            }
        });
        return 1;
    } else if (((x1 - x2) * (x1 - x2)) < (eps * eps)) {
        m12 = 3.40282347E+38f32;
    } else if (((x1 - x3) * (x1 - x3)) < (eps * eps)) {
        m13 = 3.40282347E+38f32;
    } else if (((x2 - x3) * (x2 - x3)) < (eps * eps)) {
        m23 = 3.40282347E+38f32;
    }

    if (m13 != 3.40282347E+38f32) {
        m13 = ((z1 - z3) / (x1 - x3));
    }

    if (m23 != 3.40282347E+38f32) {
        m23 = ((z2 - z3) / (x2 - x3));
    }

    if (m12 != 3.40282347E+38f32) {
        m12 = ((z1 - z2) / (x1 - x2));
    }

    b12 = (z1 - (m12 * x1));
    b13 = (z1 - (m13 * x1));
    b23 = (z2 - (m23 * x2));
    if ((((((zs - (m12 * xs)) - b12) * ((zs - (m12 * xs)) - b12)) < (eps * eps))
        || ((((zs - (m13 * xs)) - b13) * ((zs - (m13 * xs)) - b13)) < (eps * eps)))
        || ((((zs - (m23 * xs)) - b23) * ((zs - (m23 * xs)) - b23)) < (eps * eps)))
    {
        (*(z)) = ({
            if (zs < (zmax - eps)) {
                (zs + eps)
            } else {
                (zs - eps)
            }
        });
        return 1;
    } else {
        return 0;
    }
}

/// For each ray, write nt (z,x) pairs uniformly sampled in time.
pub(super) unsafe fn writeWaves(
    mut wavefp: *mut CFile,
    mut nt: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
) {
    let mut iray: i32 = 0;
    let mut it: i32 = 0;
    let mut tmax: f32 = 0.0;
    let mut dt: f32 = 0.0;
    let mut t1: f32 = 0.0;
    let mut t2: f32 = 0.0;
    let mut x1: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut ti: f32 = 0.0;
    let mut ddt: f32 = 0.0;
    let mut rsi: *mut RayStep = std::ptr::null_mut();
    let mut rs1: *mut RayStep = std::ptr::null_mut();
    let mut rs2: *mut RayStep = std::ptr::null_mut();
    {
        {
            {
                iray = 0;
                iray
            };
            {
                tmax = ((0f64) as f32);
                tmax
            }
        };
        let mut first_5 = true;
        loop {
            if !first_5 {
                iray = iray + 1;
            }
            first_5 = false;
            if !(iray < nray) {
                break;
            }
            {
                rsi = (*(rs).offset((iray) as isize));
                let mut first_6 = true;
                loop {
                    if !first_6 {
                        rsi = (*(rsi)).rsNext;
                    }
                    first_6 = false;
                    if !(rsi != ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                        break;
                    }
                    tmax = ({
                        if ((tmax) > ((*(rsi)).t)) {
                            (tmax)
                        } else {
                            ((*(rsi)).t)
                        }
                    });
                }
            }
        }
    }
    dt = (tmax
        / (({
            if (nt > 1) {
                (nt - 1)
            } else {
                1
            }
        }) as f32));
    fprintf!(stderr(), "wavefront time increment = %g\n", ((dt) as f64));
    {
        iray = 0;
        let mut first_7 = true;
        loop {
            if !first_7 {
                iray = iray + 1;
            }
            first_7 = false;
            if !(iray < nray) {
                break;
            }
            {
                rs1 = (*(rs).offset((iray) as isize));
                rs2 = (*(rs1)).rsNext;
                t1 = (*(rs1)).t;
                t2 = (*(rs2)).t;
                x1 = (*(rs1)).x;
                z1 = (*(rs1)).z;
                vgx = (*(rs1)).vgx;
                vgz = (*(rs1)).vgz;
                {
                    {
                        {
                            it = 0;
                            it
                        };
                        {
                            ti = ((0f64) as f32);
                            ti
                        }
                    };
                    let mut first_8 = true;
                    loop {
                        if !first_8 {
                            {
                                {
                                    let old = it;
                                    it = it + 1;
                                    it
                                };
                                {
                                    ti = ((ti) as f32 + (dt) as f32) as f32;
                                    ti
                                }
                            };
                        }
                        first_8 = false;
                        if !(it < nt) {
                            break;
                        }
                        {
                            if (t2 < ti) {
                                loop {
                                    {
                                        rs1 = rs2;
                                        rs2 = (*(rs1)).rsNext;
                                        if (rs2 == ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                                            break;
                                        }
                                    }
                                    if !((*(rs2)).t < ti) {
                                        break;
                                    }
                                }
                                if (rs2 == ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                                    break;
                                }

                                t1 = (*(rs1)).t;
                                t2 = (*(rs2)).t;
                                x1 = (*(rs1)).x;
                                z1 = (*(rs1)).z;
                                vgx = (*(rs1)).vgx;
                                vgz = (*(rs1)).vgz;
                            }

                            ddt = (ti - t1);
                            x = (x1 + (ddt * vgx));
                            z = (z1 + (ddt * vgz));
                            fwrite(
                                std::slice::from_raw_parts(
                                    ((&raw mut z) as *mut std::ffi::c_void) as *const f32,
                                    ((1) as usize) as usize,
                                ),
                                &*(wavefp),
                            );
                            fwrite(
                                std::slice::from_raw_parts(
                                    ((&raw mut x) as *mut std::ffi::c_void) as *const f32,
                                    ((1) as usize) as usize,
                                ),
                                &*(wavefp),
                            );
                        }
                    }
                }
                while (it < nt) {
                    fwrite(
                        std::slice::from_raw_parts(
                            ((&raw mut z) as *mut std::ffi::c_void) as *const f32,
                            ((1) as usize) as usize,
                        ),
                        &*(wavefp),
                    );
                    fwrite(
                        std::slice::from_raw_parts(
                            ((&raw mut x) as *mut std::ffi::c_void) as *const f32,
                            ((1) as usize) as usize,
                        ),
                        &*(wavefp),
                    );
                    it = it + 1;
                }
            }
        }
    }
}

/// For each ray, write its (z,x) at tw, or its start if tw is past the end.
pub(super) unsafe fn writeWaves2(
    mut wavefp: *mut CFile,
    mut tw: f32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
) {
    let mut iray: i32 = 0;
    let mut tmax: f32 = 0.0;
    let mut t1: f32 = 0.0;
    let mut t2: f32 = 0.0;
    let mut x1: f32 = 0.0;
    let mut z1: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut xstart: f32 = 0.0;
    let mut zstart: f32 = 0.0;
    let mut vgx: f32 = 0.0;
    let mut vgz: f32 = 0.0;
    let mut ddt: f32 = 0.0;
    let mut rsi: *mut RayStep = std::ptr::null_mut();
    let mut rs1: *mut RayStep = std::ptr::null_mut();
    let mut rs2: *mut RayStep = std::ptr::null_mut();
    fprintf!(stderr(), "\n wavefront plotted for tw=%g\n", ((tw) as f64));
    {
        {
            {
                iray = 0;
                iray
            };
            {
                tmax = ((0f64) as f32);
                tmax
            }
        };
        let mut first_9 = true;
        loop {
            if !first_9 {
                iray = iray + 1;
            }
            first_9 = false;
            if !(iray < nray) {
                break;
            }
            {
                rsi = (*(rs).offset((iray) as isize));
                let mut first_10 = true;
                loop {
                    if !first_10 {
                        rsi = (*(rsi)).rsNext;
                    }
                    first_10 = false;
                    if !(rsi != ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                        break;
                    }
                    tmax = ({
                        if ((tmax) > ((*(rsi)).t)) {
                            (tmax)
                        } else {
                            ((*(rsi)).t)
                        }
                    });
                }
            }
        }
    }
    {
        iray = 0;
        let mut first_11 = true;
        loop {
            if !first_11 {
                iray = iray + 1;
            }
            first_11 = false;
            if !(iray < nray) {
                break;
            }
            {
                rs1 = (*(rs).offset((iray) as isize));
                rs2 = (*(rs1)).rsNext;
                t1 = (*(rs1)).t;
                t2 = (*(rs2)).t;
                x1 = (*(rs1)).x;
                z1 = (*(rs1)).z;
                vgx = (*(rs1)).vgx;
                vgz = (*(rs1)).vgz;
                xstart = x1;
                zstart = z1;
                if (t2 < tw) {
                    loop {
                        {
                            rs1 = rs2;
                            rs2 = (*(rs1)).rsNext;
                            if (rs2 == ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                                break;
                            }
                        }
                        if !((*(rs2)).t < tw) {
                            break;
                        }
                    }
                    if (rs2 == ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                        fwrite(
                            std::slice::from_raw_parts(
                                ((&raw mut zstart) as *mut std::ffi::c_void) as *const f32,
                                ((1) as usize) as usize,
                            ),
                            &*(wavefp),
                        );
                        fwrite(
                            std::slice::from_raw_parts(
                                ((&raw mut xstart) as *mut std::ffi::c_void) as *const f32,
                                ((1) as usize) as usize,
                            ),
                            &*(wavefp),
                        );
                        continue;
                    }

                    t1 = (*(rs1)).t;
                    t2 = (*(rs2)).t;
                    x1 = (*(rs1)).x;
                    z1 = (*(rs1)).z;
                    vgx = (*(rs1)).vgx;
                    vgz = (*(rs1)).vgz;
                }

                ddt = (tw - t1);
                x = (x1 + (ddt * vgx));
                z = (z1 + (ddt * vgz));
                fwrite(
                    std::slice::from_raw_parts(
                        ((&raw mut z) as *mut std::ffi::c_void) as *const f32,
                        ((1) as usize) as usize,
                    ),
                    &*(wavefp),
                );
                fwrite(
                    std::slice::from_raw_parts(
                        ((&raw mut x) as *mut std::ffi::c_void) as *const f32,
                        ((1) as usize) as usize,
                    ),
                    &*(wavefp),
                );
            }
        }
    }
}

/// Write selected paths, padding with their last point, and update ray-end filters.
pub(super) unsafe fn writeRays(
    mut m: *mut Model,
    mut rayfp: *mut CFile,
    mut nxz: i32,
    mut nray: i32,
    mut rs: *mut *mut RayStep,
    mut re: *mut RayEnd,
    mut krecord: i32,
    mut prim: i32,
    mut ifp: *mut CFile,
    mut outparfp: *mut CFile,
) {
    let mut iray: i32 = 0;
    let mut nxzw: i32 = 0;
    let mut icount: i32 = 0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut rsi: *mut RayStep = std::ptr::null_mut();
    icount = 0;
    {
        iray = 0;
        let mut first_12 = true;
        loop {
            if !first_12 {
                iray = iray + 1;
            }
            first_12 = false;
            if !(iray < nray) {
                break;
            }
            {
                if ((krecord == 2147483647)
                    || ((krecord == (*(re).offset((iray) as isize)).kend)
                        && ((prim == 2147483647)
                            || (prim == (*(re).offset((iray) as isize)).nref))))
                {
                    icount = ((icount) as i32 + (1) as i32) as i32;
                    {
                        {
                            {
                                rsi = (*(rs).offset((iray) as isize));
                                rsi
                            };
                            {
                                nxzw = 0;
                                nxzw
                            }
                        };
                        let mut first_13 = true;
                        loop {
                            if !first_13 {
                                {
                                    {
                                        rsi = (*(rsi)).rsNext;
                                        rsi
                                    };
                                    {
                                        let old = nxzw;
                                        nxzw = nxzw + 1;
                                        nxzw
                                    }
                                };
                            }
                            first_13 = false;
                            if !(rsi != ((0) as *mut std::ffi::c_void) as *mut RayStep) {
                                break;
                            }
                            {
                                if (nxzw >= nxz) {
                                    break;
                                }

                                x = (*(rsi)).x;
                                z = (*(rsi)).z;
                                fwrite(
                                    std::slice::from_raw_parts(
                                        ((&raw mut z) as *mut std::ffi::c_void) as *const f32,
                                        ((1) as usize) as usize,
                                    ),
                                    &*(rayfp),
                                );
                                fwrite(
                                    std::slice::from_raw_parts(
                                        ((&raw mut x) as *mut std::ffi::c_void) as *const f32,
                                        ((1) as usize) as usize,
                                    ),
                                    &*(rayfp),
                                );
                            }
                        }
                    }
                    while (nxzw < nxz) {
                        fwrite(
                            std::slice::from_raw_parts(
                                ((&raw mut z) as *mut std::ffi::c_void) as *const f32,
                                ((1) as usize) as usize,
                            ),
                            &*(rayfp),
                        );
                        fwrite(
                            std::slice::from_raw_parts(
                                ((&raw mut x) as *mut std::ffi::c_void) as *const f32,
                                ((1) as usize) as usize,
                            ),
                            &*(rayfp),
                        );
                        nxzw = nxzw + 1;
                    }
                } else {
                    (*(re).offset((iray) as isize)).kend = (-(1));
                }
            }
        }
    }
    if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
        fprintf!(
            &*ifp,
            "\n Number of rays written into output : %i \n",
            icount
        );
    }

    fprintf!(&*outparfp, "\n %i\n", icount);
}

/// Compute group velocity and polarization, retaining the C operation order.
pub(super) unsafe fn gvelreal(
    mut a1111: f32,
    mut a3333: f32,
    mut a1133: f32,
    mut a1313: f32,
    mut a1113: f32,
    mut a3313: f32,
    mut px: f32,
    mut pz: f32,
    mut vgx: *mut f32,
    mut vgz: *mut f32,
    mut g11n: *mut f32,
    mut g13n: *mut f32,
    mut g33n: *mut f32,
) {
    let mut gamm11: f32 = 0.0;
    let mut gamm13: f32 = 0.0;
    let mut gamm33: f32 = 0.0;
    let mut px2: f32 = 0.0;
    let mut pz2: f32 = 0.0;
    let mut pxz: f32 = 0.0;
    let mut den: f32 = 0.0;
    let mut g11: f32 = 0.0;
    let mut g13: f32 = 0.0;
    let mut g33: f32 = 0.0;
    px2 = (px * px);
    pz2 = (pz * pz);
    pxz = (px * pz);
    gamm11 = (((a1111 * px2) + (a1313 * pz2)) + ((((2) as f32) * a1113) * pxz));
    gamm33 = (((a3333 * pz2) + (a1313 * px2)) + ((((2) as f32) * a3313) * pxz));
    gamm13 = ((((a1133 + a1313) * pxz) + (a1113 * px2)) + (a3313 * pz2));
    den = (((1) as f32) / ((gamm11 + gamm33) - ((2) as f32)));
    g11 = ((gamm33 - ((1) as f32)) * den);
    g33 = ((gamm11 - ((1) as f32)) * den);
    g13 = ((-(gamm13)) * den);
    (*(vgx)) = ((((((a1111 * px) * g11) + (((a1133 + a1313) * pz) * g13)) + ((a3313 * pz) * g33))
        + (a1113 * ((pz * g11) + ((((2) as f32) * px) * g13))))
        + ((a1313 * g33) * px));
    (*(vgz)) = ((((((a3333 * pz) * g33) + (((a1133 + a1313) * px) * g13)) + ((a1113 * px) * g11))
        + ((a3313) * ((px * g33) + ((((2) as f32) * pz) * g13))))
        + ((a1313 * g11) * pz));
    if ((g11 < (((-(10)) as f32) * 1.1920929E-7f32))
        || (g33 < (((-(10)) as f32) * 1.1920929E-7f32)))
    {
        err!("\n ERROR: g11 or g33 <0 in initialization \n");
    } else if (((g11) as f64) < 0f64) {
        g11 = ((0) as f32);
    } else if (((g33) as f64) < 0f64) {
        g33 = ((0) as f32);
    }

    (*(g11n)) = g11;
    (*(g13n)) = g13;
    (*(g33n)) = g33;
}

/// Orient polarization along slowness (Aki and Richards convention).
pub(super) unsafe fn polar(
    mut px: f32,
    mut pz: f32,
    mut g11: f32,
    mut g13: f32,
    mut g33: f32,
    mut polx: *mut f32,
    mut polz: *mut f32,
    mut mode: i32,
) {
    let mut polarx: f32 = 0.0;
    let mut polarz: f32 = 0.0;
    if (g11 < 1.1920929E-7f32) {
        polarx = ((0) as f32);
        polarz = ((((g33) as f64).sqrt()) as f32);
    } else if (g33 < 1.1920929E-7f32) {
        polarz = ((0) as f32);
        polarx = ((((g11) as f64).sqrt()) as f32);
    } else {
        polarx = ((((g11) as f64).sqrt()
            * ({
                if ((g13) < ((0) as f32)) {
                    (-(1f64))
                } else {
                    1f64
                }
            })) as f32);
        polarz = ((((g33) as f64).sqrt()) as f32);
    }

    if (((mode == 0) || (mode == 3)) && (((px * polarx) + (pz * polarz)) < ((0) as f32))) {
        polarx = (-(polarx));
        polarz = (-(polarz));
    } else if (((mode == 4) || (mode == 1)) && ((px * polarx) < ((0) as f32))) {
        polarx = (-(polarx));
        polarz = (-(polarz));
    }

    (*(polx)) = polarx;
    (*(polz)) = polarz;
}

/// Select the isotropic or anisotropic mode after optional conversion.
pub(super) unsafe fn findnewmode(
    mut mode: i32,
    mut newmode: *mut i32,
    mut conv: i32,
    mut mindex: i32,
) -> i32 {
    if ((((((mode == 0) && (mindex == 0)) && (conv == 0))
        || (((mode == 3) && (mindex == 0)) && (conv == 0)))
        || (((mode == 4) && (mindex == 0)) && (conv == 1)))
        || (((mode == 1) && (mindex == 0)) && (conv == 1)))
    {
        (*(newmode)) = 0;
    } else if ((((((mode == 1) && (mindex == 0)) && (conv == 0))
        || (((mode == 4) && (mindex == 0)) && (conv == 0)))
        || (((mode == 0) && (mindex == 0)) && (conv == 1)))
        || (((mode == 3) && (mindex == 0)) && (conv == 1)))
    {
        (*(newmode)) = 1;
    } else if ((((((mode == 2) && (mindex == 0)) && (conv == 0))
        || (((mode == 2) && (mindex == 0)) && (conv == 1)))
        || (((mode == 5) && (mindex == 0)) && (conv == 0)))
        || (((mode == 5) && (mindex == 0)) && (conv == 1)))
    {
        (*(newmode)) = 2;
    } else if ((((((mode == 0) && (mindex > 0)) && (conv == 0))
        || (((mode == 1) && (mindex > 0)) && (conv == 1)))
        || (((mode == 3) && (mindex > 0)) && (conv == 0)))
        || (((mode == 4) && (mindex > 0)) && (conv == 1)))
    {
        (*(newmode)) = 3;
    } else if ((((((mode == 0) && (mindex > 0)) && (conv == 1))
        || (((mode == 1) && (mindex > 0)) && (conv == 0)))
        || (((mode == 3) && (mindex > 0)) && (conv == 1)))
        || (((mode == 4) && (mindex > 0)) && (conv == 0)))
    {
        (*(newmode)) = 4;
    } else if ((((((mode == 2) && (mindex > 0)) && (conv == 1))
        || (((mode == 2) && (mindex > 0)) && (conv == 0)))
        || (((mode == 5) && (mindex > 0)) && (conv == 1)))
        || (((mode == 5) && (mindex > 0)) && (conv == 0)))
    {
        (*(newmode)) = 5;
    } else {
        return (-(1));
    }

    return 1;
}

/// Continue slowness across an interface; retry the other root when necessary.
pub(super) unsafe fn findqPqSV(
    mut s: f32,
    mut c: f32,
    mut pl: f32,
    mut a1111: f32,
    mut a3333: f32,
    mut a1133: f32,
    mut a1313: f32,
    mut a1113: f32,
    mut a3313: f32,
    mut mode: i32,
    mut pxnew: *mut f32,
    mut pznew: *mut f32,
    mut vgx: *mut f32,
    mut vgz: *mut f32,
    mut g11: *mut f32,
    mut g13: *mut f32,
    mut g33: *mut f32,
    mut rt: f32,
    mut ifp: *mut CFile,
) -> i32 {
    let mut check1: i32 = 0;
    let mut check2: i32 = 0;
    let mut vgm: f32 = 0.0;
    let mut vgxd: f32 = 0.0;
    let mut vgzd: f32 = 0.0;
    let mut pxnewd: f32 = 0.0;
    let mut pznewd: f32 = 0.0;
    let mut plold: f32 = 0.0;
    check2 = 0;
    check1 = check2;
    plold = pl;
    if (mode == 3) {
        check1 = solveForSlowqP(
            s,
            c,
            pl,
            a1111,
            a3333,
            a1133,
            a1313,
            a1113,
            a3313,
            &mut *(pxnew),
            &mut *(pznew),
            rt,
        );
    } else if (mode == 4) {
        check1 = solveForSlowqSV(
            s,
            c,
            pl,
            a1111,
            a3333,
            a1133,
            a1313,
            a1113,
            a3313,
            &mut *(pxnew),
            &mut *(pznew),
            rt,
        );
    } else {
        err!(" wrong mode in findqP \n ");
    }

    if (check1 == 1) {
        gvelreal(
            a1111,
            a3333,
            a1133,
            a1313,
            a1113,
            a3313,
            (*(pxnew)),
            (*(pznew)),
            vgx,
            vgz,
            g11,
            g13,
            g33,
        );
        vgxd = (*(vgx));
        vgzd = (*(vgz));
        vgm = (((-(s)) * vgxd) + (c * vgzd));
        if ((rt == ((1) as f32)) && (vgm > ((0) as f32))) {
            check2 = 1;
        } else if ((rt == ((-(1)) as f32)) && (vgm < ((0) as f32))) {
            check2 = 1;
        } else {
            check2 = 0;
        }

        if ((ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) && (check2 == 1)) {
            fprintf!(&*ifp, "\n first try in rootfinder was successfull \n");
        }
    }

    if ((check1 != 1) || (check2 != 1)) {
        if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
            fprintf!(&*ifp, "\n need second try in rootfinder \n");
        }

        if (mode == 3) {
            check1 = solveForSlowqP(
                s,
                c,
                pl,
                a1111,
                a3333,
                a1133,
                a1313,
                a1113,
                a3313,
                &mut *(pxnew),
                &mut *(pznew),
                (((-(1)) as f32) * rt),
            );
        } else if (mode == 4) {
            check1 = solveForSlowqSV(
                s,
                c,
                pl,
                a1111,
                a3333,
                a1133,
                a1313,
                a1113,
                a3313,
                &mut *(pxnew),
                &mut *(pznew),
                (((-(1)) as f32) * rt),
            );
        }

        if (check1 != 1) {
            if ifp.is_null() {
                segv();
            }
            fprintf!(&*ifp, " rootfinder not sucessful !!! \n");
            return 0;
        } else {
            gvelreal(
                a1111,
                a3333,
                a1133,
                a1313,
                a1113,
                a3313,
                (*(pxnew)),
                (*(pznew)),
                vgx,
                vgz,
                g11,
                g13,
                g33,
            );
            vgxd = (*(vgx));
            vgzd = (*(vgz));
            vgm = (((-(s)) * vgxd) + (c * vgzd));
            if ((rt == ((1) as f32)) && (vgm < ((0) as f32))) {
                return 0;
            } else if ((rt == ((-(1)) as f32)) && (vgm > ((0) as f32))) {
                return 0;
            } else {
                check2 = 1;
                if (ifp != ((0) as *mut std::ffi::c_void) as *mut CFile) {
                    fprintf!(&*ifp, "\n second try in rootfinder was successfull \n");
                }
            }
        }
    }

    pxnewd = (*(pxnew));
    pznewd = (*(pznew));
    pl = ((c * pxnewd) + (pznewd * s));
    if ((((((pl - plold) * (pl - plold)) as f64) < 1.0000000000000001E-9f64) && (check2 == 1))
        && (check1 == 1))
    {
        return 1;
    } else {
        return 0;
    }
}
