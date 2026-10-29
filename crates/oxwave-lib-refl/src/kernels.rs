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

//! Port of Refl/lib/auxil.c, compute.c, convolve.c, psvreflect.c, shreflect.c,
//! and synthetic.c. Credits: Copyright (c) Colorado School of Mines, 2011.
//! The original Credits sections contain no individual author names.
//!
//! Private pointer-based kernels preserve the C scratch-array aliases and
//! expression types. Public slice interfaces live in lib.rs.
#![allow(unused, non_snake_case, unused_parens, clippy::all)]
use super::bounds::{checked_at, BoundScope};
use super::exact::{cinv, cwp_csqrt, rcabs};
use oxwave_core::prelude::CFile;
use oxwave_cwp::complex::*;
use oxwave_f_complexf::{cipow, crpow, cwp_cexp1};
use std::any::Any;
#[derive(Default)]
struct Arena {
    storage: Vec<Box<dyn Any>>,
    bounds: BoundScope,
}
impl Arena {
    fn alloc<T: Default + Clone + 'static>(&mut self, n: usize) -> *mut T {
        // The C singularity branch starts its zeroing loop at ibl-1,
        // which is -1 when the first slowness is singular. A guard element
        // makes that otherwise undefined scratch write explicit and isolated.
        let mut data = vec![T::default(); n.checked_add(1).unwrap()].into_boxed_slice();
        let p = unsafe { data.as_mut_ptr().add(1) };
        self.bounds.track_workspace(p, n);
        self.storage.push(Box::new(data));
        p
    }
    unsafe fn alloc3(&mut self, n1: i32, n2: i32, n3: i32) -> *mut *mut *mut Complex {
        let a = usize::try_from(n1).unwrap();
        let b = usize::try_from(n2).unwrap();
        let c = usize::try_from(n3).unwrap();
        let data = self.alloc::<Complex>(a.checked_mul(b).unwrap().checked_mul(c).unwrap());
        let rows = self.alloc::<*mut Complex>(b.checked_mul(c).unwrap());
        let planes = self.alloc::<*mut *mut Complex>(c);
        for k in 0..n3 {
            *planes.add(k as usize) = rows.add(k as usize * b);
            self.bounds
                .track_slice(std::slice::from_raw_parts(rows.add(k as usize * b), b));
            for j in 0..n2 {
                *rows.add(k as usize * b + j as usize) =
                    data.add((k as usize * b + j as usize) * a);
                self.bounds.track_slice(std::slice::from_raw_parts(
                    data.add((k as usize * b + j as usize) * a),
                    a,
                ));
            }
        }
        planes
    }
}
// Refl/lib/auxil.c: compute_w_aux_arrays
pub(super) unsafe fn raw_compute_w_aux_arrays(
    mut wtype: i32,
    mut layern: i32,
    mut nlayers: i32,
    mut nor: i32,
    mut lsource: i32,
    mut np: *mut i32,
    mut block_size: *mut i32,
    mut nblock: *mut i32,
    mut left: *mut i32,
    mut fref: f32,
    mut wrefp: f32,
    mut wrefs: f32,
    mut lobs: *mut i32,
    mut tsec: f32,
    mut p2w: f32,
    mut fs: f32,
    mut xmax: f32,
    mut w: f32,
    mut decay: f32,
    mut epsp: f32,
    mut epss: f32,
    mut sigp: f32,
    mut sigs: f32,
    mut pw1: *mut f32,
    mut pw2: *mut f32,
    mut pw3: *mut f32,
    mut pw4: *mut f32,
    mut dp: *mut f32,
    mut fp: *mut f32,
    mut wpie: Complex,
    mut cdp: *mut Complex,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut ql: *mut f32,
    mut qt: *mut f32,
    mut rho: *mut f32,
    mut t: *mut f32,
    mut al: *mut Complex,
    mut at: *mut Complex,
    mut prs: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        /* float wref;	*/
        /* reference frequency in rads/sec */
        let mut eps: f32 = 0.0;
        let mut sigma: f32 = 0.0;
        let mut e4: f32 = 0.0;
        let mut dk: f32 = 0.0;
        /* auxiliary variables */
        /* initialize varaibles */
        eps = (0.001f64 as f32);
        sigma = (0.10000000000000001f64 as f32);
        dk = (0.80000000000000004f64 as f32);
        e4 = (0.80000000000000004f64 as f32);
        /* wref=2*PI*fref; */
        /* compute al and at */
        raw_compute_al_at(
            wtype,
            nlayers,
            layern,
            eps,
            epsp,
            epss,
            sigma,
            sigp,
            sigs,
            (&raw mut wrefp),
            (&raw mut wrefs),
            fref,
            cl,
            ct,
            ql,
            qt,
            rho,
            wpie,
            al,
            at,
        );
        /* compute prs for every receiver */
        raw_compute_prs(nor, lobs, al, at, rho, prs);
        /* compute number of slowness, slowness increment and maximum slowness */
        raw_compute_slowness(
            wtype, lsource, w, p2w, np, fp, fs, e4, dk, decay, xmax, dp, pw1, pw2, pw3, pw4, cdp,
            cl, ct, t,
        );
        /* compute number of blocks and block size for matrix computations */
        raw_compute_block_size((*np), block_size, nblock, left);
    }
}

// Refl/lib/auxil.c: compute_p_aux_arrays
pub(super) unsafe fn raw_compute_p_aux_arrays(
    mut wtype: i32,
    mut nlayers: i32,
    mut lsource: i32,
    mut block_size: i32,
    mut bp: f32,
    mut dp: f32,
    mut w: f32,
    mut decay: f32,
    mut pw1: f32,
    mut pw2: f32,
    mut pw3: f32,
    mut pw4: f32,
    mut pwin: *mut Complex,
    mut m1: f32,
    mut m2: f32,
    mut m3: f32,
    mut h1: f32,
    mut h2: f32,
    mut wpie: Complex,
    mut divfac: *mut Complex,
    mut p: *mut Complex,
    mut pp: *mut Complex,
    mut al: *mut Complex,
    mut at: *mut Complex,
    mut rho: *mut f32,
    mut gl: *mut Complex,
    mut gt: *mut Complex,
    mut gam: *mut Complex,
    mut alpha: *mut Complex,
    mut betha: *mut Complex,
    mut sigmad1: *mut Complex,
    mut sigmad2: *mut Complex,
    mut sigmau1: *mut Complex,
    mut sigmau2: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut rhos: f32 = 0.0;
        let mut als: Complex = cmplx(0.0, 0.0);
        let mut ats: Complex = cmplx(0.0, 0.0);
        let mut alphasq: Complex = cmplx(0.0, 0.0);
        let mut bethasq: Complex = cmplx(0.0, 0.0);
        let mut meu: Complex = cmplx(0.0, 0.0);
        let mut lambda: Complex = cmplx(0.0, 0.0);
        let mut s1: Complex = cmplx(0.0, 0.0);
        let mut s2: Complex = cmplx(0.0, 0.0);
        let mut s4: Complex = cmplx(0.0, 0.0);
        let mut a1: Complex = cmplx(0.0, 0.0);
        let mut a2: Complex = cmplx(0.0, 0.0);
        let mut a3: Complex = cmplx(0.0, 0.0);
        /* save al, at and rho values for the source layer */
        als = (*checked_at(al, (lsource - 1i32) as isize));
        ats = (*checked_at(at, (lsource - 1i32) as isize));
        rhos = (*checked_at(rho, (lsource - 1i32) as isize));
        /* compute complex auxiliary variables */
        alphasq = cdiv(cmplx((1f64 as f32), (0f64 as f32)), cmul(als, als));
        if ((((ats).r as f64) == 0f64) && (((ats).i as f64) == 0f64)) {
            bethasq = cmplx((0f64 as f32), (0f64 as f32));
        } else {
            bethasq = cdiv(cmplx((1f64 as f32), (0f64 as f32)), cmul(ats, ats));
        }

        meu = crmul(bethasq, rhos);
        if ((((meu).r as f64) != 0f64) || (((meu).i as f64) != 0f64)) {
            a1 = cdiv(cmplx((1f64 as f32), (0f64 as f32)), meu);
        } else {
            a1 = cmplx((0f64 as f32), (0f64 as f32));
        }

        lambda = crmul(
            cmul(csub(alphasq, cmplx((2i32 as f32), (0i32 as f32))), bethasq),
            rhos,
        );
        (*divfac) = cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie);
        s1 = crmul(a1, m2);
        s2 = crmul(crmul(cmul(als, als), ((1f64 / (rhos as f64)) as f32)), m3);
        s4 = cdiv(cmplx(h2, (0f64 as f32)), (*divfac));
        a2 = cdiv(cmplx(h1, (0f64 as f32)), (*divfac));
        a3 = csub(cmul(lambda, s2), cmplx(m1, (0f64 as f32)));
        /* compute pwin */
        raw_compute_pwin(
            wtype, block_size, bp, dp, w, decay, pw1, pw2, pw3, pw4, pwin, p, pp,
        );
        /* compute gl and gt and gam */
        raw_compute_gl_gt_gam(wtype, nlayers, block_size, al, at, rho, pp, gl, gt, gam);
        /* compute alpha and betha */
        raw_compute_alpha_betha(lsource, block_size, ats, rhos, gl, gt, alpha, betha);
        /* compute sigmau1,sigmau2,sigmad1,sigmad2 */
        raw_compute_sigmas(
            wtype, block_size, lsource, m1, m2, meu, alpha, betha, p, gl, gt, gam, a2, a3, pwin,
            s1, s2, s4, sigmau1, sigmau2, sigmad1, sigmad2,
        );
    }
}

// Refl/lib/auxil.c: source_receiver_type
pub(super) unsafe fn raw_source_receiver_type(
    mut nor: i32,
    mut nlayers: i32,
    mut lsource: i32,
    mut lobs: *mut i32,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut acoustic: *mut i32,
    mut flag: *mut i32,
) {
    let mut arena = Arena::default();
    {
        let mut k: i32 = 0i32;
        let mut iz: i32 = 0;
        /* loop counters */
        let mut jkl: i32 = 0;
        let mut lmn: i32 = 0;
        /* auxiliary indices */
        /* loop over the receivers */
        {
            iz = 0i32;
            while (iz < nor) {
                {
                    jkl = (lsource - 1i32);
                    lmn = ((*checked_at(lobs, (iz) as isize)) - 1i32);
                    (*checked_at(acoustic, (iz) as isize)) = 0i32;
                    (*checked_at(flag, (iz) as isize)) = 0i32;
                    /* test condition for flag */
                    if (((*checked_at(cl, (nlayers - 1i32) as isize))
                        == (*checked_at(cl, (lmn) as isize)))
                        && ((*checked_at(ct, (nlayers - 1i32) as isize))
                            == (*checked_at(ct, (lmn) as isize))))
                    {
                        (*checked_at(flag, (iz) as isize)) = 1i32;
                    }

                    /* test condition for source and receiver type */
                    if (((*checked_at(ct, (jkl) as isize)) as f64) == 0f64) {
                        /* source acoustic */
                        if (((*checked_at(ct, (lmn) as isize)) as f64) == 0f64) {
                            /* receiver acoustic */
                            (*checked_at(acoustic, (iz) as isize)) = 1i32;
                        } else {
                            if (lmn > (lsource - 1i32)) {
                                k += 1;
                            }

                            if (k == 1i32) {
                                /* 1st elastic rec below sorce*/
                                (*checked_at(acoustic, (iz) as isize)) = 2i32;
                            } else {
                                /* elastic receiver below surface */
                                (*checked_at(acoustic, (iz) as isize)) = 3i32;
                            }
                        }
                    } else {
                        /* source elastic */
                        if (((*checked_at(ct, (lmn) as isize)) as f64) == 0f64) {
                            /* receiver acoustic */
                            (*checked_at(acoustic, (iz) as isize)) = 1i32;
                        }
                    }
                }
                iz += 1;
            }
        }
        /* only for Po/So */
        if (((nor <= 1i32) && (lsource > (*checked_at(lobs, (0i32) as isize))))
            && (((*checked_at(ct, (0i32) as isize)) as f64) == 0f64))
        {
            (*checked_at(acoustic, (0i32) as isize)) = 4i32;
        }
    }
}

// Refl/lib/auxil.c: compute_slowness
pub(super) unsafe fn raw_compute_slowness(
    mut wtype: i32,
    mut lsource: i32,
    mut w: f32,
    mut p2w: f32,
    mut np: *mut i32,
    mut fp: *mut f32,
    mut fs: f32,
    mut e4: f32,
    mut dk: f32,
    mut decay: f32,
    mut xmax: f32,
    mut dp: *mut f32,
    mut pw1: *mut f32,
    mut pw2: *mut f32,
    mut pw3: *mut f32,
    mut pw4: *mut f32,
    mut cdp: *mut Complex,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut t: *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut numw: i32 = 0;
        /* number of frequencies */
        let mut index: i32 = 0i32;
        let mut fq: f32 = 0.0;
        /* frequency in hertz */
        let mut fp1: f32 = 0.0;
        let mut fp2: f32 = 0.0;
        let mut dkt: f32 = 0.0;
        /* auxiliary variables */
        let mut bp: f32 = 0.0;
        /* beginning and final ray parameters */
        /* initialize variables */
        fq = (((w as f64) / ((3.1415926535897931f64) * 2f64)) as f32);
        bp = (0f64 as f32);
        if (wtype == 1i32) {
            index = (lsource - 1i32);
        } else {
            if (wtype == 2i32) {
                index = 0i32;
            }
        }

        /* compute fp1 and fp2 */
        if (((*checked_at(ct, (index) as isize)) as f64) <= 5.0000000000000001E-4f64) {
            fp1 = (((1f64 / ((*checked_at(cl, (lsource - 1i32) as isize)) as f64))
                .powf((2i32 as f64))
                + (50f64
                    / ((w as f64)
                        * (((*checked_at(t, (lsource - 1i32) as isize)) as f64)
                            + 1.0000000000000001E-5f64)))
                    .powf((2i32 as f64)))
            .sqrt() as f32);
            fp2 = (((1f64 / ((*checked_at(cl, (lsource - 2i32) as isize)) as f64))
                .powf((2i32 as f64))
                + (50f64
                    / ((w as f64)
                        * (((*checked_at(t, (lsource - 2i32) as isize)) as f64)
                            + 1.0000000000000001E-5f64)))
                    .powf((2i32 as f64)))
            .sqrt() as f32);
        } else {
            fp1 = (((1f64 / ((*checked_at(ct, (lsource - 1i32) as isize)) as f64))
                .powf((2i32 as f64))
                + (50f64
                    / ((w as f64)
                        * (((*checked_at(t, (lsource - 1i32) as isize)) as f64)
                            + 1.0000000000000001E-5f64)))
                    .powf((2i32 as f64)))
            .sqrt() as f32);
            fp2 = (((1f64 / ((*checked_at(ct, (lsource - 2i32) as isize)) as f64))
                .powf((2i32 as f64))
                + (50f64
                    / ((w as f64)
                        * (((*checked_at(t, (lsource - 2i32) as isize)) as f64)
                            + 1.0000000000000001E-5f64)))
                    .powf((2i32 as f64)))
            .sqrt() as f32);
        }

        /* compute the slowness parameters */
        (*fp) = (if ((fp1) > (fp2)) { (fp1) } else { (fp2) });
        dkt = ((((dk as f64) * (3.1415926535897931f64)) / (xmax as f64)) as f32);
        (*dp) =
            ((((dkt / w) as f64) * (1f64 + ((fq / fs) as f64)).powf(((e4 as f64) / 2f64))) as f32);
        (*dp) = ((if ((0.002f64) < ((*dp) as f64)) {
            (0.002f64)
        } else {
            ((*dp) as f64)
        }) as f32);
        if (((p2w as f64) > 0f64) && ((*fp) >= p2w)) {
            (*fp) = p2w;
        }

        (*np) = (((*fp) / (*dp)) as i32);
        (*fp) = (((*np) as f32) * (*dp));
        (*cdp) = cmplx((*dp), ((-(*dp)) * (decay / w)));
        /* if required, use default values for window slowness parameters */
        if ((((*pw2) as f64) <= 0f64) && (((*pw3) as f64) <= 0f64)) {
            numw = (((((*np) as f64) * 0.14999999999999999f64) - (1i32 as f64)) as i32);
            (*pw1) = bp;
            (*pw2) = (bp + ((*dp) * (numw as f32)));
            (*pw3) = ((*fp) - ((*dp) * (numw as f32)));
            (*pw4) = (*fp);
        }
    }
}

// Refl/lib/auxil.c: compute_al_at
pub(super) unsafe fn raw_compute_al_at(
    mut wtype: i32,
    mut nlayers: i32,
    mut layern: i32,
    mut eps: f32,
    mut epsp: f32,
    mut epss: f32,
    mut sigma: f32,
    mut sigp: f32,
    mut sigs: f32,
    mut wrefp: *mut f32,
    mut wrefs: *mut f32,
    mut fw: f32,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut ql: *mut f32,
    mut qt: *mut f32,
    mut rho: *mut f32,
    mut wpie: Complex,
    mut al: *mut Complex,
    mut at: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut i: i32 = 0;
        /* loop counter */
        let mut qwif: f32 = 0.0;
        let mut qwifs: f32 = 0.0;
        let mut qwifp: f32 = 0.0;
        /* auxilairy variables */
        let mut wref: f32 = 0.0;
        /* reference frequency in rad/sec */
        let mut wrefpp: f32 = 0.0;
        let mut wrefss: f32 = 0.0;
        /* reference p ands frequencies in rad/sec */
        let mut eiw: Complex = cmplx(0.0, 0.0);
        let mut eiwp: Complex = cmplx(0.0, 0.0);
        let mut eiws: Complex = cmplx(0.0, 0.0);
        /* auxiliary variables */
        /* compute reference frequencies for p and s waves in rad/s */
        wrefpp = (*wrefp);
        wrefss = (*wrefs);
        wref = (((2f64 * (3.1415926535897931f64)) * (fw as f64)) as f32);
        wrefpp = ((wrefpp as f64) * ((2f64 * (3.1415926535897931f64)) as f64)) as f32;
        wrefss = ((wrefss as f64) * ((2f64 * (3.1415926535897931f64)) as f64)) as f32;
        /* compute the new p and s q values */
        qwif = (((((eps * eps) + (wref * wref)) as f64).powf(((sigma as f64) / 2f64))
            / ((((sigma as f64) * (3.1415926535897931f64)) / 2f64).sin() * 2f64))
            as f32);
        qwifp = (((((epsp * epsp) + (wrefpp * wrefpp)) as f64).powf(((sigp as f64) / 2f64))
            / ((((sigp as f64) * (3.1415926535897931f64)) / 2f64).sin() * 2f64))
            as f32);
        qwifs = (((((epss * epss) + (wrefss * wrefss)) as f64).powf(((sigs as f64) / 2f64))
            / ((((sigs as f64) * (3.1415926535897931f64)) / 2f64).sin() * 2f64))
            as f32);
        /* compute wie */
        eiw = crpow(
            csub(
                cmplx(eps, (0f64 as f32)),
                cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie),
            ),
            sigma,
        );
        /* check to see if layern was provided as an input */
        if (layern != 0i32) {
            eiwp = crpow(
                csub(
                    cmplx(epsp, (0f64 as f32)),
                    cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie),
                ),
                sigp,
            );
            eiws = crpow(
                csub(
                    cmplx(epss, (0f64 as f32)),
                    cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie),
                ),
                sigs,
            );
            /* compute al and at until just before layern */
            {
                i = 0i32;
                while (i < (layern - 1i32)) {
                    {
                        (*checked_at(al, (i) as isize)) = cmul(
                            cmplx(
                                ((1f64 / ((*checked_at(cl, (i) as isize)) as f64)) as f32),
                                (0f64 as f32),
                            ),
                            cadd(
                                cdiv(
                                    cmplx((qwif / (*checked_at(ql, (i) as isize))), (0f64 as f32)),
                                    eiw,
                                ),
                                cmplx((1f64 as f32), (0f64 as f32)),
                            ),
                        );
                        if (((*checked_at(ct, (i) as isize)) as f64) == 0f64) {
                            (*checked_at(at, (i) as isize)) = cmplx((0f64 as f32), (0f64 as f32));
                        } else {
                            (*checked_at(at, (i) as isize)) = cmul(
                                cmplx(
                                    ((1f64 / ((*checked_at(ct, (i) as isize)) as f64)) as f32),
                                    (0f64 as f32),
                                ),
                                cadd(
                                    cdiv(
                                        cmplx(
                                            (qwif / (*checked_at(qt, (i) as isize))),
                                            (0f64 as f32),
                                        ),
                                        eiw,
                                    ),
                                    cmplx((1f64 as f32), (0f64 as f32)),
                                ),
                            );
                        }
                    }
                    i += 1;
                }
            }
            /* compute al and at for the remaining layers */
            {
                i = (layern - 1i32);
                while (i < nlayers) {
                    {
                        (*checked_at(al, (i) as isize)) = cmul(
                            cmplx(
                                ((1f64 / ((*checked_at(cl, (i) as isize)) as f64)) as f32),
                                (0f64 as f32),
                            ),
                            cadd(
                                cdiv(
                                    cmplx((qwifp / (*checked_at(ql, (i) as isize))), (0f64 as f32)),
                                    eiwp,
                                ),
                                cmplx((1f64 as f32), (0f64 as f32)),
                            ),
                        );
                        if (((*checked_at(ct, (i) as isize)) as f64) == 0f64) {
                            (*checked_at(at, (i) as isize)) = cmplx((0f64 as f32), (0f64 as f32));
                        } else {
                            (*checked_at(at, (i) as isize)) = cmul(
                                cmplx(
                                    ((1f64 / ((*checked_at(ct, (i) as isize)) as f64)) as f32),
                                    (0f64 as f32),
                                ),
                                cadd(
                                    cdiv(
                                        cmplx(
                                            (qwifs / (*checked_at(qt, (i) as isize))),
                                            (0f64 as f32),
                                        ),
                                        eiws,
                                    ),
                                    cmplx((1f64 as f32), (0f64 as f32)),
                                ),
                            );
                        }
                    }
                    i += 1;
                }
            }
        } else {
            /* compute al and at for all layers */
            {
                i = 0i32;
                while (i < nlayers) {
                    {
                        (*checked_at(al, (i) as isize)) = cmul(
                            cmplx(
                                ((1f64 / ((*checked_at(cl, (i) as isize)) as f64)) as f32),
                                (0f64 as f32),
                            ),
                            cadd(
                                cdiv(
                                    cmplx((qwif / (*checked_at(ql, (i) as isize))), (0f64 as f32)),
                                    eiw,
                                ),
                                cmplx((1f64 as f32), (0f64 as f32)),
                            ),
                        );
                        if (((*checked_at(ct, (i) as isize)) as f64) == 0f64) {
                            (*checked_at(at, (i) as isize)) = cmplx((0f64 as f32), (0f64 as f32));
                        } else {
                            (*checked_at(at, (i) as isize)) = cmul(
                                cmplx(
                                    ((1f64 / ((*checked_at(ct, (i) as isize)) as f64)) as f32),
                                    (0f64 as f32),
                                ),
                                cadd(
                                    cdiv(
                                        cmplx(
                                            (qwif / (*checked_at(qt, (i) as isize))),
                                            (0f64 as f32),
                                        ),
                                        eiw,
                                    ),
                                    cmplx((1f64 as f32), (0f64 as f32)),
                                ),
                            );
                        }
                    }
                    i += 1;
                }
            }
        }

        /* update pointers */
        (*wrefp) = wrefpp;
        (*wrefs) = wrefss;
    }
}

// Refl/lib/auxil.c: compute_prs
pub(super) unsafe fn raw_compute_prs(
    mut nor: i32,
    mut lobs: *mut i32,
    mut al: *mut Complex,
    mut at: *mut Complex,
    mut rho: *mut f32,
    mut prs: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut ijk: i32 = 0;
        let mut ijk1: i32 = 0;
        let mut ctemp: Complex = cmplx(0.0, 0.0);
        let mut rp: f32 = 0.0;
        /* loop over the receivers */
        {
            ijk = 0i32;
            while (ijk < nor) {
                {
                    ijk1 = ((*checked_at(lobs, (ijk) as isize)) - 1i32);
                    rp = (*checked_at(at, (ijk1) as isize)).r;
                    if ((rp as f64) == 0f64) {
                        (*checked_at(prs, (ijk) as isize)) = crmul(
                            cdiv(
                                cmplx((1f64 as f32), (0f64 as f32)),
                                cmul(
                                    (*checked_at(al, (ijk1) as isize)),
                                    (*checked_at(al, (ijk1) as isize)),
                                ),
                            ),
                            (*checked_at(rho, (ijk1) as isize)),
                        );
                    } else {
                        ctemp = cmul(
                            cdiv(
                                (*checked_at(al, (ijk1) as isize)),
                                (*checked_at(at, (ijk1) as isize)),
                            ),
                            cdiv(
                                (*checked_at(al, (ijk1) as isize)),
                                (*checked_at(at, (ijk1) as isize)),
                            ),
                        );
                        (*checked_at(prs, (ijk) as isize)) = crmul(
                            csub(
                                cmplx((1f64 as f32), (0f64 as f32)),
                                crmul(ctemp, ((4f64 / 3f64) as f32)),
                            ),
                            (*checked_at(rho, (ijk1) as isize)),
                        );
                    }
                }
                ijk += 1;
            }
        }
    }
}

// Refl/lib/auxil.c: compute_block_size
pub(super) unsafe fn raw_compute_block_size(
    mut np: i32,
    mut block_size: *mut i32,
    mut nblock: *mut i32,
    mut left: *mut i32,
) {
    let mut arena = Arena::default();
    {
        if (np <= 100i32) {
            /* one block is enough and nothing is left */
            (*block_size) = np;
            (*nblock) = 1i32;
            (*left) = 0i32;
        } else {
            /* define block size and compute number of required blocks */
            (*block_size) = 100i32;
            (*left) = (np % 100i32);
            (*nblock) = ((np / 100i32) + 1i32);
        }
    }
}

// Refl/lib/auxil.c: compute_pwin
pub(super) unsafe fn raw_compute_pwin(
    mut wtype: i32,
    mut block_size: i32,
    mut bp: f32,
    mut dp: f32,
    mut w: f32,
    mut decay: f32,
    mut pw1: f32,
    mut pw2: f32,
    mut pw3: f32,
    mut pw4: f32,
    mut pwin: *mut Complex,
    mut p: *mut Complex,
    mut pp: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut ip: i32 = 0;
        /* loop counter			*/
        let mut rp: f32 = 0.0;
        /* current ray parameter	*/
        let mut win: f32 = (0f64 as f32);
        let mut win1: f32 = 0.0;
        let mut win2: f32 = 0.0;
        /* auxiliary variables		*/
        /* loop over blocks */
        {
            ip = 0i32;
            while (ip < block_size) {
                {
                    /* compute array of complex slownesses */
                    rp = (bp + ((ip as f32) * dp));
                    if (wtype == 1i32) {
                        (*checked_at(p, (ip) as isize)) = cmplx(rp, (0f64 as f32));
                    } else {
                        if (wtype == 2i32) {
                            (*checked_at(p, (ip) as isize)) = cmplx(rp, (((-decay) * rp) / w));
                        }
                    }

                    (*checked_at(pp, (ip) as isize)) = cmul(
                        (*checked_at(p, (ip) as isize)),
                        (*checked_at(p, (ip) as isize)),
                    );
                    /* compute Hanning tapered window */
                    if ((rp >= pw1) && (rp < pw2)) {
                        win1 = ((0.5f64
                            * (1f64
                                + ((3.1415926535897931f64)
                                    * ((((rp - pw1) / (pw2 - pw1)) as f64) - 1f64))
                                    .cos())) as f32);
                        win = win1;
                    }

                    if ((rp >= pw2) && (rp <= pw3)) {
                        win = (1f64 as f32);
                    }

                    if ((rp > pw3) && (rp < pw4)) {
                        win2 = ((0.5f64
                            * (1f64
                                + ((3.1415926535897931f64) * (((rp - pw3) / (pw4 - pw3)) as f64))
                                    .cos())) as f32);
                        win = win2;
                    }

                    if ((rp < pw1) || (rp >= pw4)) {
                        win = (0f64 as f32);
                    }

                    /* multiply complex ray parameter by window to compute pwin */
                    (*checked_at(pwin, (ip) as isize)) =
                        crmul(cwp_csqrt((*checked_at(p, (ip) as isize))), win);
                }
                ip += 1;
            }
        }
    }
}

// Refl/lib/auxil.c: compute_gl_gt_gam
pub(super) unsafe fn raw_compute_gl_gt_gam(
    mut wtype: i32,
    mut nlayers: i32,
    mut block_size: i32,
    mut al: *mut Complex,
    mut at: *mut Complex,
    mut rho: *mut f32,
    mut pp: *mut Complex,
    mut gl: *mut Complex,
    mut gt: *mut Complex,
    mut gam: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut il: i32 = 0;
        let mut ip: i32 = 0;
        let mut ijk: i32 = 0i32;
        /* loop counters */
        let mut ik1: i32 = 0;
        let mut ijk1: i32 = 0;
        /* auxiliary indices */
        let mut rho1: f32 = 0.0;
        /* current density */
        let mut p1: Complex = cmplx(0.0, 0.0);
        let mut p2: Complex = cmplx(0.0, 0.0);
        /* loop over layers */
        {
            il = 0i32;
            while (il < nlayers) {
                {
                    p1 = cmul(
                        (*checked_at(al, (il) as isize)),
                        (*checked_at(al, (il) as isize)),
                    );
                    p2 = cmul(
                        (*checked_at(at, (il) as isize)),
                        (*checked_at(at, (il) as isize)),
                    );
                    rho1 = (*checked_at(rho, (il) as isize));
                    ik1 = (il * block_size);
                    /* if p2 is not complex number zero */
                    if ((((p2).r as f64) != 0f64) || (((p2).i as f64) != 0f64)) {
                        {
                            ip = 0i32;
                            while (ip < block_size) {
                                {
                                    ijk1 = (ik1 + ip);
                                    /* compute gl and gt */
                                    (*checked_at(gl, (ijk1) as isize)) =
                                        cwp_csqrt(csub(p1, (*checked_at(pp, (ip) as isize))));
                                    (*checked_at(gt, (ijk1) as isize)) =
                                        cwp_csqrt(csub(p2, (*checked_at(pp, (ip) as isize))));
                                    /* make imaginary parts of gl and gt positive */
                                    if (((*checked_at(gl, (ijk1) as isize)).i as f64) < 0f64) {
                                        (*checked_at(gl, (ijk1) as isize)).i =
                                            (((*checked_at(gl, (ijk1) as isize)).i as f64)
                                                * ((-1f64) as f64))
                                                as f32;
                                    }

                                    if (((*checked_at(gt, (ijk1) as isize)).i as f64) < 0f64) {
                                        (*checked_at(gt, (ijk1) as isize)).i =
                                            (((*checked_at(gt, (ijk1) as isize)).i as f64)
                                                * ((-1f64) as f64))
                                                as f32;
                                    }

                                    if (wtype == 2i32) {
                                        /* make imaginary parts of gl and gt positive */
                                        if (((*checked_at(gl, (ijk1) as isize)).i as f64) < 0f64) {
                                            (*checked_at(gl, (ijk1) as isize)).r =
                                                (((*checked_at(gl, (ijk1) as isize)).r as f64)
                                                    * ((-1f64) as f64))
                                                    as f32;
                                        }

                                        if (((*checked_at(gt, (ijk1) as isize)).i as f64) < 0f64) {
                                            (*checked_at(gt, (ijk1) as isize)).r =
                                                (((*checked_at(gt, (ijk1) as isize)).r as f64)
                                                    * ((-1f64) as f64))
                                                    as f32;
                                        }
                                    }

                                    /* compute gam */
                                    (*checked_at(gam, (ijk1) as isize)) = crmul(
                                        csub(
                                            cmplx((1f64 as f32), (0f64 as f32)),
                                            crmul(
                                                cdiv((*checked_at(pp, (ip) as isize)), p2),
                                                (2f64 as f32),
                                            ),
                                        ),
                                        rho1,
                                    );
                                }
                                ip += 1;
                            }
                        }
                    } else {
                        /* p2 is complex number zero */
                        {
                            ip = 0i32;
                            while (ip < block_size) {
                                {
                                    ijk1 = (ik1 + ip);
                                    /* compute gl and set gt to zero */
                                    (*checked_at(gl, (ijk1) as isize)) =
                                        cwp_csqrt(csub(p1, (*checked_at(pp, (ip) as isize))));
                                    (*checked_at(gt, (ijk1) as isize)) =
                                        cmplx((0f64 as f32), (0f64 as f32));
                                    /* if necessary, make imaginary part of gl positive */
                                    if (((*checked_at(gl, (ijk) as isize)).i as f64) < 0f64) {
                                        (*checked_at(gl, (ijk1) as isize)) =
                                            cneg((*checked_at(gl, (ijk1) as isize)));
                                    }

                                    (*checked_at(gam, (ijk1) as isize)) =
                                        cmplx(rho1, (0f64 as f32));
                                }
                                ip += 1;
                            }
                        }
                    }
                }
                il += 1;
            }
        }
    }
}

// Refl/lib/auxil.c: compute_alpha_betha
pub(super) unsafe fn raw_compute_alpha_betha(
    mut lsource: i32,
    mut block_size: i32,
    mut ats: Complex,
    mut rhos: f32,
    mut gl: *mut Complex,
    mut gt: *mut Complex,
    mut alpha: *mut Complex,
    mut betha: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut ip: i32 = 0;
        /* loop counter */
        let mut ik1: i32 = 0;
        let mut ijk1: i32 = 0;
        /* auxiliary indices */
        ik1 = ((lsource - 1i32) * block_size);
        if ((((ats).r as f64) == 0f64) && (((ats).i as f64) == 0f64)) {
            {
                ip = 0i32;
                while (ip < block_size) {
                    {
                        ijk1 = (ik1 + ip);
                        (*checked_at(alpha, (ip) as isize)) = cdiv(
                            cmplx((1f64 as f32), (0f64 as f32)),
                            crmul(
                                crmul((*checked_at(gl, (ijk1) as isize)), rhos),
                                (2f64 as f32),
                            ),
                        );
                        (*checked_at(betha, (ip) as isize)) = cmplx((0f64 as f32), (0f64 as f32));
                    }
                    ip += 1;
                }
            }
        } else {
            {
                ip = 0i32;
                while (ip < block_size) {
                    {
                        ijk1 = (ik1 + ip);
                        (*checked_at(alpha, (ip) as isize)) = cdiv(
                            cmplx((1f64 as f32), (0f64 as f32)),
                            crmul(
                                crmul((*checked_at(gl, (ijk1) as isize)), rhos),
                                (2f64 as f32),
                            ),
                        );
                        (*checked_at(betha, (ip) as isize)) = cdiv(
                            cmplx((1f64 as f32), (0f64 as f32)),
                            crmul(
                                crmul((*checked_at(gt, (ijk1) as isize)), rhos),
                                (2f64 as f32),
                            ),
                        );
                    }
                    ip += 1;
                }
            }
        }
    }
}

// Refl/lib/auxil.c: compute_sigmas
pub(super) unsafe fn raw_compute_sigmas(
    mut wtype: i32,
    mut block_size: i32,
    mut lsource: i32,
    mut m1: f32,
    mut m2: f32,
    mut meu: Complex,
    mut alpha: *mut Complex,
    mut betha: *mut Complex,
    mut p: *mut Complex,
    mut gl: *mut Complex,
    mut gt: *mut Complex,
    mut gam: *mut Complex,
    mut a2: Complex,
    mut a3: Complex,
    mut pwin: *mut Complex,
    mut s1: Complex,
    mut s2: Complex,
    mut s4: Complex,
    mut sigmau1: *mut Complex,
    mut sigmau2: *mut Complex,
    mut sigmad1: *mut Complex,
    mut sigmad2: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut ip: i32 = 0;
        /* loop counter */
        let mut ijk1: i32 = 0;
        let mut ik1: i32 = 0;
        /* auxiliary indices */
        let mut d1: Complex = cmplx(0.0, 0.0);
        let mut d2: Complex = cmplx(0.0, 0.0);
        let mut d3: Complex = cmplx(0.0, 0.0);
        let mut d4: Complex = cmplx(0.0, 0.0);
        let mut d5: Complex = cmplx(0.0, 0.0);
        let mut d6: Complex = cmplx(0.0, 0.0);
        let mut d7: Complex = cmplx(0.0, 0.0);
        let mut d8: Complex = cmplx(0.0, 0.0);
        let mut s3: Complex = cmplx(0.0, 0.0);
        /* auxiliary variables */
        /* compute auxiliary index */
        ik1 = ((lsource - 1i32) * block_size);
        /* loop over layers in a block */
        {
            ip = 0i32;
            while (ip < block_size) {
                {
                    ijk1 = (ik1 + ip);
                    if (wtype == 1i32) {
                        /* compute auxiliary variables */
                        d1 = crmul(
                            cmul(
                                cmul(meu, (*checked_at(gl, (ijk1) as isize))),
                                cmul(
                                    (*checked_at(alpha, (ip) as isize)),
                                    (*checked_at(p, (ip) as isize)),
                                ),
                            ),
                            (2f64 as f32),
                        );
                        d2 = cneg(cmul(
                            (*checked_at(gam, (ijk1) as isize)),
                            (*checked_at(alpha, (ip) as isize)),
                        ));
                        d3 = cmul(
                            (*checked_at(p, (ip) as isize)),
                            (*checked_at(alpha, (ip) as isize)),
                        );
                        d4 = cneg(cmul(
                            (*checked_at(gl, (ijk1) as isize)),
                            (*checked_at(alpha, (ip) as isize)),
                        ));
                        d5 = cmul(
                            (*checked_at(gam, (ijk1) as isize)),
                            (*checked_at(betha, (ip) as isize)),
                        );
                        d6 = crmul(
                            cmul(
                                cmul(meu, (*checked_at(gt, (ijk1) as isize))),
                                cmul(
                                    (*checked_at(betha, (ip) as isize)),
                                    (*checked_at(p, (ip) as isize)),
                                ),
                            ),
                            (2f64 as f32),
                        );
                        d7 = cmul(
                            (*checked_at(gt, (ijk1) as isize)),
                            (*checked_at(betha, (ip) as isize)),
                        );
                        d8 = cmul(
                            (*checked_at(p, (ip) as isize)),
                            (*checked_at(betha, (ip) as isize)),
                        );
                        s3 = cadd(a2, cmul((*checked_at(p, (ip) as isize)), a3));
                        /* compute sigmas */
                        (*checked_at(sigmau1, (ip) as isize)) = cneg(cmul(
                            cadd(
                                cadd(cmul(d1, s1), cmul(d2, s2)),
                                cadd(cmul(d3, s3), cmul(d4, s4)),
                            ),
                            (*checked_at(pwin, (ip) as isize)),
                        ));
                        (*checked_at(sigmau2, (ip) as isize)) = cneg(cmul(
                            cadd(
                                cadd(cmul(d5, s1), cmul(d6, s2)),
                                cadd(cmul(d7, s3), cmul(d8, s4)),
                            ),
                            (*checked_at(pwin, (ip) as isize)),
                        ));
                        (*checked_at(sigmad1, (ip) as isize)) = cmul(
                            cadd(
                                csub(cmul(d1, s1), cmul(d2, s2)),
                                csub(cmul(d4, s4), cmul(d3, s3)),
                            ),
                            (*checked_at(pwin, (ip) as isize)),
                        );
                        (*checked_at(sigmad2, (ip) as isize)) = cmul(
                            cadd(
                                csub(cmul(d6, s2), cmul(d5, s1)),
                                csub(cmul(d7, s3), cmul(d8, s4)),
                            ),
                            (*checked_at(pwin, (ip) as isize)),
                        );
                    } else {
                        if (wtype == 2i32) {
                            /* avoid p and gt being zero */
                            if ((((*checked_at(p, (ip) as isize)).r as f64) == 0f64)
                                && (((*checked_at(p, (ip) as isize)).i as f64) == 0f64))
                            {
                                (*checked_at(p, (ip) as isize)) =
                                    cmplx((1.0E-4f64 as f32), (0f64 as f32));
                            }

                            if ((((*checked_at(gt, (ijk1) as isize)).r as f64) == 0f64)
                                && (((*checked_at(gt, (ijk1) as isize)).i as f64) == 0f64))
                            {
                                (*checked_at(gt, (ijk1) as isize)) =
                                    cmplx((1.0E-4f64 as f32), (0f64 as f32));
                            }

                            /* compute sigmas */
                            (*checked_at(sigmau2, (ip) as isize)) =
                                cmplx((0f64 as f32), (0f64 as f32));
                            (*checked_at(sigmad2, (ip) as isize)) =
                                cmplx((0f64 as f32), (0f64 as f32));
                            (*checked_at(sigmau1, (ip) as isize)) = cneg(cmul(
                                cadd(
                                    crmul(
                                        cinv(cmul(meu, (*checked_at(p, (ip) as isize)))),
                                        (((m2 as f64) / 2f64) as f32),
                                    ),
                                    crmul(
                                        cinv(cmul(meu, (*checked_at(gt, (ijk1) as isize)))),
                                        (((m1 as f64) / 2f64) as f32),
                                    ),
                                ),
                                (*checked_at(pwin, (ip) as isize)),
                            ));
                            (*checked_at(sigmad1, (ip) as isize)) =
                                (*checked_at(sigmau1, (ip) as isize));
                        }
                    }
                }
                ip += 1;
            }
        }
    }
}

// Refl/lib/auxil.c: compute_moment_tensor
pub(super) unsafe fn raw_compute_moment_tensor(
    mut wtype: i32,
    mut phi: f32,
    mut lambda: f32,
    mut delta: f32,
    mut phis: f32,
    mut m0: f32,
    mut m1: *mut f32,
    mut m2: *mut f32,
    mut m3: *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut a: f32 = 0.0;
        let mut b: f32 = 0.0;
        let mut c: f32 = 0.0;
        let mut d: f32 = 0.0;
        let mut e: f32 = 0.0;
        let mut f: f32 = 0.0;
        /* auxiliary variables */
        /* update fault plane azimuth */
        phis = ((phis as f32) - (phi as f32)) as f32;
        /* convert angles to radians */
        delta = ((delta as f64) * (((3.1415926535897931f64) / 180f64) as f64)) as f32;
        lambda = ((lambda as f64) * (((3.1415926535897931f64) / 180f64) as f64)) as f32;
        phis = ((phis as f64) * (((3.1415926535897931f64) / 180f64) as f64)) as f32;
        if (wtype == 1i32) {
            /* compute auxiliary variables */
            a = ((phis as f64).sin() as f32);
            b = (a * a);
            c = ((((delta as f64).sin() * (lambda as f64).cos()) * (2f64 * (phis as f64)).sin())
                as f32);
            d = (((2f64 * (delta as f64)).sin() * (lambda as f64).sin()) as f32);
            e = ((((delta as f64).cos() * (lambda as f64).cos()) * (phis as f64).cos()) as f32);
            f = ((((2f64 * (delta as f64)).cos() * (lambda as f64).sin()) * (a as f64)) as f32);
            /* compute the moment tensor components */
            (*m1) = ((-m0) * (c + (d * b)));
            (*m2) = ((-m0) * (e + f));
            (*m3) = (m0 * d);
        } else {
            if (wtype == 2i32) {
                (*m1) = (((m0 as f64)
                    * ((((delta as f64).sin() * (lambda as f64).cos())
                        * (2f64 * (phis as f64)).cos())
                        + (((0.5f64 * (2f64 * (delta as f64)).sin()) * (lambda as f64).sin())
                            * (2f64 * (phis as f64)).sin()))) as f32);
                (*m2) = ((((-m0) as f64)
                    * ((((delta as f64).cos() * (lambda as f64).cos()) * (phis as f64).sin())
                        - (((2f64 * (delta as f64)).cos() * (lambda as f64).sin())
                            * (phis as f64).cos()))) as f32);
                (*m3) = (0f64 as f32);
            }
        }
    }
}

// Refl/lib/auxil.c: parameter_interpolation
pub(super) unsafe fn raw_parameter_interpolation(
    mut nlayers: i32,
    mut intlayers: *mut i32,
    mut nintlayers: *mut i32,
    mut intlayth: *mut f32,
    mut cl: *mut f32,
    mut ql: *mut f32,
    mut ct: *mut f32,
    mut qt: *mut f32,
    mut rho: *mut f32,
    mut t: *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut i: i32 = 0i32;
        let mut ii: i32 = 0i32;
        let mut j: i32 = 0;
        /* loop counters */
        let mut nintlay: i32 = 0;
        /* number of layer to interpolate */
        let mut tintlay: f32 = 0.0;
        /* thickness at which to interpolate */
        let mut incrcl: f32 = 0.0;
        /* interpolation step for cl */
        let mut incrql: f32 = 0.0;
        /* interpolation step for ql */
        let mut incrct: f32 = 0.0;
        /* interpolation step for ct */
        let mut incrqt: f32 = 0.0;
        /* interpolation step for qt */
        let mut incrrho: f32 = 0.0;
        /* interpolation step for rho */
        let mut crust: f32 = (0f64 as f32);
        /* auxiliary variable */
        while (i < nlayers) {
            if ((*checked_at(intlayers, (ii) as isize)) == i) {
                /* interpolation requested for this layer */
                /* get number of layers to interpolate */
                nintlay = (*checked_at(nintlayers, (ii) as isize));
                tintlay = (*checked_at(intlayth, (ii) as isize));
                /* check if thickness is negative */
                if ((tintlay as f64) < 0f64) {
                    tintlay = ((-tintlay) - crust);
                }

                crust = ((crust as f32) + (tintlay as f32)) as f32;
                /* if required, adjust input parameters */
                if (((*checked_at(ct, (i + nintlay) as isize)) as f64) == (-1f64)) {
                    (*checked_at(ct, (i + nintlay) as isize)) =
                        ((((*checked_at(cl, (i + nintlay) as isize)) as f64) / (3f64).sqrt())
                            as f32);
                }

                if (((*checked_at(qt, (i + nintlay) as isize)) as f64) == (-1f64)) {
                    (*checked_at(qt, (i + nintlay) as isize)) =
                        ((((*checked_at(ql, (i + nintlay) as isize)) as f64) / 2f64) as f32);
                }

                if (((*checked_at(rho, ((i + nintlay) - 1i32) as isize)) as f64) == (-1f64)) {
                    (*checked_at(rho, (i + nintlay) as isize)) =
                        (((((*checked_at(cl, (i + nintlay) as isize)) as f64) + 1.5f64) / 3f64)
                            as f32);
                }

                /* compute the interpolation steps */
                incrcl = (((*checked_at(cl, (i + nintlay) as isize))
                    - (*checked_at(cl, (i - 1i32) as isize)))
                    / ((nintlay + 1i32) as f32));
                incrql = (((*checked_at(ql, (i + nintlay) as isize))
                    - (*checked_at(ql, (i - 1i32) as isize)))
                    / ((nintlay + 1i32) as f32));
                incrct = (((*checked_at(ct, (i + nintlay) as isize))
                    - (*checked_at(ct, (i - 1i32) as isize)))
                    / ((nintlay + 1i32) as f32));
                incrqt = (((*checked_at(qt, (i + nintlay) as isize))
                    - (*checked_at(qt, (i - 1i32) as isize)))
                    / ((nintlay + 1i32) as f32));
                incrrho = (((*checked_at(rho, (i + nintlay) as isize))
                    - (*checked_at(rho, (i - 1i32) as isize)))
                    / ((nintlay + 1i32) as f32));
                /* do the actual interpolation */
                {
                    j = 0i32;
                    while (j < nintlay) {
                        {
                            (*checked_at(cl, ((j + i) - 1i32) as isize)) =
                                ((*checked_at(cl, (i - 1i32) as isize))
                                    + (((j + 1i32) as f32) * incrcl));
                            (*checked_at(ql, ((j + i) - 1i32) as isize)) =
                                ((*checked_at(ql, (i - 1i32) as isize))
                                    + (((j + 1i32) as f32) * incrql));
                            (*checked_at(ct, ((j + i) - 1i32) as isize)) =
                                ((*checked_at(ct, (i - 1i32) as isize))
                                    + (((j + 1i32) as f32) * incrct));
                            (*checked_at(qt, ((j + i) - 1i32) as isize)) =
                                ((*checked_at(qt, (i - 1i32) as isize))
                                    + (((j + 1i32) as f32) * incrqt));
                            (*checked_at(rho, ((j + i) - 1i32) as isize)) =
                                ((*checked_at(rho, (i - 1i32) as isize))
                                    + (((j + 1i32) as f32) * incrrho));
                            (*checked_at(t, ((j + i) - 1i32) as isize)) =
                                (tintlay / (nintlay as f32));
                        }
                        j += 1;
                    }
                }
                /* update counter by number of interopolated layers */
                i = ((i as i32) + (nintlay as i32)) as i32;
                ii += 1;
            } else {
                /* no layer interpolation requested */
                /* if required, adjust input parameters */
                if (((*checked_at(ct, (i) as isize)) as f64) == (-1f64)) {
                    (*checked_at(ct, (i) as isize)) =
                        ((((*checked_at(cl, (i) as isize)) as f64) / (3f64).sqrt()) as f32);
                }

                if (((*checked_at(qt, (i) as isize)) as f64) == (-1f64)) {
                    (*checked_at(qt, (i) as isize)) =
                        ((((*checked_at(ql, (i) as isize)) as f64) / 2f64) as f32);
                }

                if (((*checked_at(rho, (i) as isize)) as f64) == (-1f64)) {
                    (*checked_at(rho, (i) as isize)) =
                        (((((*checked_at(cl, (i) as isize)) as f64) + 1.5f64) / 3f64) as f32);
                }

                if (((*checked_at(t, (i) as isize)) as f64) < 0f64) {
                    (*checked_at(t, (i) as isize)) = ((-(*checked_at(t, (i) as isize))) - crust);
                }

                if (i < (nlayers - 1i32)) {
                    crust = ((crust as f32) + ((*checked_at(t, (i) as isize)) as f32)) as f32;
                } else {
                    (*checked_at(t, (i) as isize)) = (1.0E+4f64 as f32);
                }

                /* make last layer a subspace */
                i += 1;
            }
        }
    }
}

// Refl/lib/auxil.c: random_velocity_layers
pub(super) unsafe fn raw_random_velocity_layers(
    mut nlayers: *mut i32,
    mut lsource: *mut i32,
    mut nrand_layers: i32,
    mut sdcl: f32,
    mut sdct: f32,
    mut layer: f32,
    mut zlayer: f32,
    mut cl: *mut f32,
    mut ql: *mut f32,
    mut ct: *mut f32,
    mut qt: *mut f32,
    mut rho: *mut f32,
    mut t: *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut i: i32 = 0i32;
        let mut il: i32 = 0;
        let mut j: i32 = 0i32;
        /* loop counters */
        let mut icnt: i32 = 0;
        /* count number of layers to insert */
        let mut ijk: i32 = 0;
        let mut ksource: i32 = 0i32;
        /* index to update source layer number*/
        let mut nl: i32 = (*nlayers);
        /* current number of layers */
        let mut nlmax: i32 = 0;
        let mut d1: f32 = 0.0;
        let mut d2: f32 = 0.0;
        let mut tlast: f32 = 0.0;
        let mut x: f32 = 0.0;
        let mut zl: f32 = 0.0;
        /* auxiliary varibles */
        let mut cls: f32 = 0.0;
        let mut cts: f32 = 0.0;
        let mut qls: f32 = 0.0;
        let mut qts: f32 = 0.0;
        let mut rhos: f32 = 0.0;
        /* auxiliary variables */
        let mut cll: *mut f32 = std::ptr::null_mut();
        let mut ctt: *mut f32 = std::ptr::null_mut();
        let mut qll: *mut f32 = std::ptr::null_mut();
        /* scratch arrays */
        let mut qtt: *mut f32 = std::ptr::null_mut();
        let mut rhoo: *mut f32 = std::ptr::null_mut();
        let mut tt: *mut f32 = std::ptr::null_mut();
        /* more scratch arrays */
        oxwave_core::fprintf!(
            &oxwave_core::prelude::stderr(),
            "before rand: nlayers=%d ",
            nl
        );
        /* allocate working space */
        cll = arena.alloc::<f32>((nl as usize) as usize);
        ctt = arena.alloc::<f32>((nl as usize) as usize);
        qll = arena.alloc::<f32>((nl as usize) as usize);
        qtt = arena.alloc::<f32>((nl as usize) as usize);
        rhoo = arena.alloc::<f32>((nl as usize) as usize);
        tt = arena.alloc::<f32>((nl as usize) as usize);
        /* initialize variables */
        d1 = (((2f64 * (3f64).sqrt()) * ((sdcl as f64) / 100f64)) as f32);
        d2 = (((2f64 * (3f64).sqrt()) * ((sdct as f64) / 100f64)) as f32);
        /* save information from last layer */
        cls = (*checked_at(cl, (nl - 1i32) as isize));
        cts = (*checked_at(ct, (nl - 1i32) as isize));
        qls = (*checked_at(ql, (nl - 1i32) as isize));
        qts = (*checked_at(qt, (nl - 1i32) as isize));
        rhos = (*checked_at(rho, (nl - 1i32) as isize));
        nlmax = (((((nrand_layers + nl) as f32) - layer) + (1i32 as f32)) as i32);
        /* main loop to compute the random layers */
        {
            il = ((layer - (1i32 as f32)) as i32);
            while (il < (nl - 1i32)) {
                {
                    if ((*checked_at(t, (il) as isize)) <= zlayer) {
                        icnt = 1i32;
                        /* compute just one layer */
                        zl = (*checked_at(t, (il) as isize));
                    } else
                    /* with thickness t[il] */
                    {
                        tlast =
                            ((((*checked_at(t, (il) as isize)) as f64) % (zlayer as f64)) as f32);
                        icnt = (((*checked_at(t, (il) as isize)) / zlayer) as i32);
                        /* make icnt layers with */
                        zl = (zlayer + (tlast / (icnt as f32)));
                    }

                    /* thickness t[il]/zlayer */
                    if (il == ((*lsource) - 1i32)) {
                        ksource = j;
                    }

                    {
                        ijk = 0i32;
                        while (ijk < icnt) {
                            {
                                (*checked_at(tt, (i) as isize)) = zl;
                                /* all layers, same thickness */
                                x = (((oxwave_f_franuni::franuni() as f64) - 0.5f64) as f32);
                                (*checked_at(cll, (i) as isize)) =
                                    ((((*checked_at(cl, (il) as isize)) as f64)
                                        * (1f64 + ((d1 * x) as f64)))
                                        as f32);
                                /* random velocities */
                                (*checked_at(ctt, (i) as isize)) =
                                    ((((*checked_at(ct, (il) as isize)) as f64)
                                        * (1f64 + ((d2 * x) as f64)))
                                        as f32);
                                (*checked_at(qll, (i) as isize)) = (*checked_at(ql, (il) as isize));
                                /* same q's and rho */
                                (*checked_at(qtt, (i) as isize)) = (*checked_at(qt, (il) as isize));
                                (*checked_at(rhoo, (i) as isize)) =
                                    (*checked_at(rho, (il) as isize));
                                /* if maximum number of layers reached, exit inner loop */
                                if (i == nlmax) {
                                    break;
                                }

                                /* update counters */
                                i += 1;
                                if (il < ((*lsource) - 1i32)) {
                                    j += 1;
                                }
                            }
                            ijk += 1;
                        }
                    }
                    /* if maximum number of layers reached, exit outer loop */
                    if (i == nlmax) {
                        break;
                    }
                }
                il += 1;
            }
        }
        /* update number of layers and source layer index */
        nl = ((((i as f32) + layer) - (1i32 as f32)) as i32);
        (*lsource) = (((*lsource) as i32) + ((ksource - 1i32) as i32)) as i32;
        /* restore information for last layer */
        (*checked_at(cl, (nl - 1i32) as isize)) = cls;
        (*checked_at(ql, (nl - 1i32) as isize)) = qls;
        (*checked_at(ct, (nl - 1i32) as isize)) = cts;
        (*checked_at(qt, (nl - 1i32) as isize)) = qts;
        (*checked_at(rho, (nl - 1i32) as isize)) = rhos;
        /* main loop to insert the computed information for the random layers */
        {
            {
                il = ((layer - (1i32 as f32)) as i32);
                i = 0i32;
            };
            while (il < (nl - 1i32)) {
                {
                    (*checked_at(cl, (il) as isize)) = (*checked_at(cll, (i) as isize));
                    (*checked_at(ct, (il) as isize)) = (*checked_at(ctt, (i) as isize));
                    (*checked_at(rho, (il) as isize)) = (*checked_at(rhoo, (i) as isize));
                    (*checked_at(ql, (il) as isize)) = (*checked_at(qll, (i) as isize));
                    (*checked_at(qt, (il) as isize)) = (*checked_at(qtt, (i) as isize));
                    (*checked_at(t, (il) as isize)) = (*checked_at(tt, (i) as isize));
                    if (il == ((*lsource) - 1i32)) {
                        (*checked_at(cl, (il) as isize)) = (*checked_at(cl, (il - 1i32) as isize));
                        (*checked_at(ct, (il) as isize)) = (*checked_at(ct, (il - 1i32) as isize));
                        (*checked_at(ql, (il) as isize)) = (*checked_at(ql, (il - 1i32) as isize));
                        (*checked_at(qt, (il) as isize)) = (*checked_at(qt, (il - 1i32) as isize));
                        (*checked_at(rho, (il) as isize)) =
                            (*checked_at(rho, (il - 1i32) as isize));
                    }
                }
                {
                    il += 1;
                    i += 1;
                };
            }
        }
        (*checked_at(t, (nl - 1i32) as isize)) = (1.0E+4f64 as f32);
        /* make last layer a subspace */
        /* update pointer for number of layers */
        (*nlayers) = nl;
        oxwave_core::fprintf!(
            &oxwave_core::prelude::stderr(),
            "after rand: nlayers=%d\n",
            nl
        );
        /* free allocated space */
        ();
        ();
        ();
        ();
        ();
        ();
    }
}

// Refl/lib/auxil.c: apply_earth_flattening
pub(super) unsafe fn raw_apply_earth_flattening(
    mut nlayers: i32,
    mut z0: f32,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut rho: *mut f32,
    mut t: *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut il: i32 = 0;
        /* loop counter */
        let mut y0: f32 = 0.0;
        let mut w0: f32 = 0.0;
        /* auxiliary variables */
        let mut z: f32 = z0;
        /* layer depth */
        let mut zz: *mut f32 = std::ptr::null_mut();
        /* scratch array */
        /* allocate working space */
        zz = arena.alloc::<f32>((nlayers as usize) as usize);
        /* create and array of absolute depths from the thicknesses */
        (*checked_at(zz, (0i32) as isize)) = z;
        {
            il = 0i32;
            while (il < (nlayers - 1i32)) {
                {
                    (*checked_at(zz, (il + 1i32) as isize)) = (z + (*checked_at(t, (il) as isize)));
                    z = (*checked_at(zz, (il + 1i32) as isize));
                }
                il += 1;
            }
        }
        /* apply the earth flattening correction */
        {
            il = 0i32;
            while (il < nlayers) {
                {
                    /* compute the correction */
                    y0 = ((6371f64 / (6371f64 - ((*checked_at(zz, (il) as isize)) as f64))) as f32);
                    (*checked_at(zz, (il) as isize)) = ((6371f64 * (y0 as f64).ln()) as f32);
                    w0 = ((6371f64 / (6371f64 - ((*checked_at(zz, (il) as isize)) as f64))) as f32);
                    /* apply the correction */
                    (*checked_at(cl, (il) as isize)) =
                        (((*checked_at(cl, (il) as isize)) as f32) * (w0 as f32)) as f32;
                    (*checked_at(ct, (il) as isize)) =
                        (((*checked_at(ct, (il) as isize)) as f32) * (w0 as f32)) as f32;
                    (*checked_at(rho, (il) as isize)) =
                        (((*checked_at(rho, (il) as isize)) as f32) / (w0 as f32)) as f32;
                }
                il += 1;
            }
        }
        /* update the thicknesses */
        {
            il = 1i32;
            while (il < nlayers) {
                (*checked_at(t, (il - 1i32) as isize)) =
                    ((*checked_at(zz, (il) as isize)) - (*checked_at(zz, (il - 1i32) as isize)));
                il += 1;
            }
        }
        /* make last layer a half space*/
        (*checked_at(t, (nlayers - 1i32) as isize)) = (1.0E+4f64 as f32);
        /* free allocated space */
        ();
    }
}

// Refl/lib/compute.c: compute_reflectivities
pub(super) unsafe fn raw_compute_reflectivities(
    mut int_type: i32,
    mut verbose: i32,
    mut wtype: i32,
    mut wfield: i32,
    mut vsp: i32,
    mut flt: i32,
    mut win: i32,
    mut nx: i32,
    mut nt: i32,
    mut ntc: i32,
    mut nor: i32,
    mut nf: i32,
    mut nlayers: i32,
    mut lsource: i32,
    mut layern: i32,
    mut nfilters: i32,
    mut filters_phase: *mut i32,
    mut nw: i32,
    mut np: i32,
    mut bp: f32,
    mut tlag: f32,
    mut red_vel: f32,
    mut w1: f32,
    mut w2: f32,
    mut fx: f32,
    mut dx: f32,
    mut bx: f32,
    mut fs: f32,
    mut decay: f32,
    mut p2w: f32,
    mut tsec: f32,
    mut fw: f32,
    mut wrefp: f32,
    mut wrefs: f32,
    mut epsp: f32,
    mut epss: f32,
    mut sigp: f32,
    mut sigs: f32,
    mut pw1: f32,
    mut pw2: f32,
    mut pw3: f32,
    mut pw4: f32,
    mut h1: f32,
    mut h2: f32,
    mut m1: f32,
    mut m2: f32,
    mut m3: f32,
    mut fref: f32,
    mut lobs: *mut i32,
    mut filters_type: *mut i32,
    mut dbpo: *mut f32,
    mut f1: *mut f32,
    mut f2: *mut f32,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut ql: *mut f32,
    mut qt: *mut f32,
    mut rho: *mut f32,
    mut t: *mut f32,
    mut wavefield1: *mut *mut f32,
    mut wavefield2: *mut *mut f32,
    mut wavefield3: *mut *mut f32,
    mut outfp: *mut CFile,
) {
    let mut arena = Arena::default();
    {
        let mut il: i32 = 0;
        /* loop counters */
        let mut xmax: f32 = 0.0;
        /* maximum range */
        let mut acoustic: *mut i32 = std::ptr::null_mut();
        /* flags for receiver type */
        let mut flag: *mut i32 = std::ptr::null_mut();
        /* ??? */
        let mut response1: *mut *mut *mut Complex = std::ptr::null_mut();
        /* scratch array for pressure */
        let mut response2: *mut *mut *mut Complex = std::ptr::null_mut();
        /* scratch array for x-wavefield */
        let mut response3: *mut *mut *mut Complex = std::ptr::null_mut();
        /* scratch array for z-wavefield */
        /* initialize variables */
        oxwave_core::fprintf!(
            &oxwave_core::prelude::stderr(),
            "original decay factor: %g nx=%d\n",
            (decay as f64),
            nx
        );
        fx = (bx + (dx * ((nx - 1i32) as f32)));
        xmax = ((6i32 as f32) * fx);
        if ((decay as f64) <= 0f64) {
            decay = (50f64 as f32);
        }

        decay = (((decay as f64).ln() / (tsec as f64)) as f32);
        /* allocate working space */
        acoustic = arena.alloc::<i32>((nor as usize) as usize);
        flag = arena.alloc::<i32>((nor as usize) as usize);
        response1 = arena.alloc3(
            (nor as usize) as i32,
            (nx as usize) as i32,
            (nw as usize) as i32,
        );
        if (wtype == 1i32) {
            response2 = arena.alloc3(
                (nor as usize) as i32,
                (nx as usize) as i32,
                (nw as usize) as i32,
            );
            response3 = arena.alloc3(
                (nor as usize) as i32,
                (nx as usize) as i32,
                (nw as usize) as i32,
            );
        }

        /* determine type of sources and receivers */
        raw_source_receiver_type(nor, nlayers, lsource, lobs, cl, ct, acoustic, flag);
        if ((verbose == 1i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(&oxwave_core::prelude::stderr(), "Number of traces %d\n", nx);
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Offset of first trace %g\n",
                (bx as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Offset of last trace %g\n",
                (fx as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Trace to trace offset %g\n",
                (dx as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Number of samples in output traces %d\n",
                nt
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Time sampling interval %g\n",
                ((tsec / (nt as f32)) as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Trace length in seconds %g\n",
                (tsec as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Time lag applied to each trace (sec)%g\n",
                (tlag as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Total number of layers processed %d\n",
                nlayers
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Total number of frequencies processed %d\n",
                nw
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Decay factor applied=%g\n",
                (decay as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Number of frequencies in out traces %d\n",
                nf
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "First frequency in out traces %g\n",
                ((fref / tsec) as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "Last frequency in output traces %g\n",
                (((nw as f32) / tsec) as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "\nPROCESSING INFORMATION FOR EACH FREQUENCY\n"
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "  iw   w\tbp    fp    np    dp    pw1\tpw2\tpw3\tpw4   \n"
            );
        }

        if ((verbose == 2i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(&*outfp, "Number of traces %d\n", nx);
            oxwave_core::fprintf!(&*outfp, "Offset of first trace %g\n", (bx as f64));
            oxwave_core::fprintf!(&*outfp, "Offset of last trace %g\n", (fx as f64));
            oxwave_core::fprintf!(&*outfp, "Trace to trace offset %g\n", (dx as f64));
            oxwave_core::fprintf!(&*outfp, "Number of samples in output traces %d\n", nt);
            oxwave_core::fprintf!(
                &*outfp,
                "Time sampling interval %g\n",
                ((tsec / (nt as f32)) as f64)
            );
            oxwave_core::fprintf!(&*outfp, "Trace length in seconds %g\n", (tsec as f64));
            oxwave_core::fprintf!(
                &*outfp,
                "Time lag applied to each trace in sec%g\n",
                (tlag as f64)
            );
            oxwave_core::fprintf!(&*outfp, "Total number of layers processed %d\n", nlayers);
            oxwave_core::fprintf!(&*outfp, "Total number of frequencies processed %d\n", nw);
            oxwave_core::fprintf!(&*outfp, "Decay factor applied=%g\n", (decay as f64));
            oxwave_core::fprintf!(&*outfp, "Number of frequencies in out traces %d\n", nf);
            oxwave_core::fprintf!(
                &*outfp,
                "First frequency in out traces %g\n",
                ((fref / tsec) as f64)
            );
            oxwave_core::fprintf!(
                &*outfp,
                "Last frequency in output traces %g\n",
                (((nw as f32) / tsec) as f64)
            );
            oxwave_core::fprintf!(&*outfp, "\nPROCESSING INFORMATION FOR EACH FREQUENCY\n");
            oxwave_core::fprintf!(
                &*outfp,
                "  iw   w\tbp    fp   np    dp    pw1\tpw2    pw3\tpw4\n"
            );
        }

        /* compute reflectivities */
        if (wtype == 1i32) {
            /* PSV */
            raw_psv_reflectivities(
                int_type, verbose, wtype, nw, nlayers, nx, layern, nor, np, bp, m1, m2, m3, h1, h2,
                lsource, bx, dx, xmax, decay, fref, wrefp, wrefs, tsec, p2w, fs, epsp, epss, sigp,
                sigs, pw1, pw2, pw3, pw4, acoustic, flag, lobs, rho, t, cl, ct, ql, qt, response1,
                response2, response3, outfp,
            );
        } else {
            if (wtype == 2i32) {
                raw_sh_reflectivities(
                    int_type, verbose, wtype, nw, nlayers, nx, layern, nor, np, bp, m1, m2, m3, h1,
                    h2, lsource, bx, dx, xmax, decay, fref, wrefp, wrefs, tsec, p2w, fs, epsp,
                    epss, sigp, sigs, pw1, pw2, pw3, pw4, flag, lobs, rho, t, cl, ct, ql, qt,
                    response1, outfp,
                );
            }
        }

        /* if requested, output processing information */
        if ((verbose == 1i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "\nMODIFIED INPUT PARAMETERS\n"
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "il   cl    ct\tql\tqt\trho\tt\n"
            );
            {
                il = 0i32;
                while (il < nlayers) {
                    {
                        oxwave_core::fprintf!(
                            &oxwave_core::prelude::stderr(),
                            "%3d%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                            il,
                            ((*checked_at(cl, (il) as isize)) as f64),
                            ((*checked_at(ct, (il) as isize)) as f64),
                            ((*checked_at(ql, (il) as isize)) as f64),
                            ((*checked_at(qt, (il) as isize)) as f64),
                            ((*checked_at(rho, (il) as isize)) as f64),
                            ((*checked_at(t, (il) as isize)) as f64)
                        );
                    }
                    il += 1;
                }
            }
            oxwave_core::fprintf!(&oxwave_core::prelude::stderr(), "\n");
        }

        if ((verbose == 2i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(&*outfp, "\nMODIFIED INPUT PARAMETERS\n");
            oxwave_core::fprintf!(&*outfp, "il    cl\tct\tql\tqt\trho    t  \n");
            {
                il = 0i32;
                while (il < nlayers) {
                    {
                        oxwave_core::fprintf!(
                            &*outfp,
                            "%3d%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                            il,
                            ((*checked_at(cl, (il) as isize)) as f64),
                            ((*checked_at(ct, (il) as isize)) as f64),
                            ((*checked_at(ql, (il) as isize)) as f64),
                            ((*checked_at(qt, (il) as isize)) as f64),
                            ((*checked_at(rho, (il) as isize)) as f64),
                            ((*checked_at(t, (il) as isize)) as f64)
                        );
                    }
                    il += 1;
                }
            }
            oxwave_core::fprintf!(&*outfp, "\n");
        }

        /* compute actual synthetic in t-x domain */
        if (wtype == 1i32) {
            /* pressure */
            raw_compute_synthetics(
                verbose,
                nt,
                ntc,
                nx,
                nor,
                nw,
                nlayers,
                lsource,
                nf,
                flt,
                vsp,
                win,
                wfield,
                tlag,
                red_vel,
                decay,
                tsec,
                bx,
                dx,
                w1,
                w2,
                filters_phase,
                nfilters,
                filters_type,
                dbpo,
                f1,
                f2,
                lobs,
                cl,
                t,
                response1,
                wavefield1,
                outfp,
            );
            /* radial component */
            raw_compute_synthetics(
                verbose,
                nt,
                ntc,
                nx,
                nor,
                nw,
                nlayers,
                lsource,
                nf,
                flt,
                vsp,
                win,
                wfield,
                tlag,
                red_vel,
                decay,
                tsec,
                bx,
                dx,
                w1,
                w2,
                filters_phase,
                nfilters,
                filters_type,
                dbpo,
                f1,
                f2,
                lobs,
                cl,
                t,
                response2,
                wavefield2,
                outfp,
            );
            /* vertical component */
            raw_compute_synthetics(
                verbose,
                nt,
                ntc,
                nx,
                nor,
                nw,
                nlayers,
                lsource,
                nf,
                flt,
                vsp,
                win,
                wfield,
                tlag,
                red_vel,
                decay,
                tsec,
                bx,
                dx,
                w1,
                w2,
                filters_phase,
                nfilters,
                filters_type,
                dbpo,
                f1,
                f2,
                lobs,
                cl,
                t,
                response3,
                wavefield3,
                outfp,
            );
            /* free allocated space */
            ();
            ();
            ();
        } else {
            if (wtype == 2i32) {
                /* tangential component */
                raw_compute_synthetics(
                    verbose,
                    nt,
                    nx,
                    ntc,
                    nor,
                    nw,
                    nlayers,
                    lsource,
                    nf,
                    flt,
                    vsp,
                    win,
                    wfield,
                    tlag,
                    red_vel,
                    decay,
                    tsec,
                    bx,
                    dx,
                    w1,
                    w2,
                    filters_phase,
                    nfilters,
                    filters_type,
                    dbpo,
                    f1,
                    f2,
                    lobs,
                    cl,
                    t,
                    response1,
                    wavefield1,
                    outfp,
                );
                /* free allocated space */
                ();
            } else {
                oxwave_core::err!("wtype flag has to be one ro two");
            }
        }

        /* free allocated space */
        ();
        ();
    }
}

// Refl/lib/convolve.c: convolve_wavelet
pub(super) unsafe fn raw_convolve_wavelet(
    mut wavelet_type: i32,
    mut nx: i32,
    mut nt: i32,
    mut dt: f32,
    mut fpeak: f32,
    mut wfield: *mut *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut ix: i32 = 0;
        let mut it: i32 = 0;
        /* loop counters */
        let mut wavelet: *mut f32 = std::ptr::null_mut();
        /* wavelet array */
        let mut temp: *mut f32 = std::ptr::null_mut();
        /* scratch array */
        /* allocate working space */
        wavelet = arena.alloc::<f32>((nt as usize) as usize);
        temp = arena.alloc::<f32>((nt as usize) as usize);
        /* compute wavelet */
        if (wavelet_type == 1i32) {
            oxwave_waveforms::spike_wavelet(
                (1i32) as usize,
                std::slice::from_raw_parts_mut(wavelet, (nt) as usize),
            );
        } else {
            if (wavelet_type == 2i32) {
                oxwave_waveforms::ricker1_wavelet(
                    dt,
                    fpeak,
                    std::slice::from_raw_parts_mut(wavelet, (nt) as usize),
                );
            } else {
                if (wavelet_type == 3i32) {
                    oxwave_waveforms::ricker2_wavelet(
                        (nt / 2i32) as usize,
                        dt,
                        ((1f64 / (fpeak as f64)) as f32),
                        (1f64 as f32),
                        (0f64 as f32),
                        std::slice::from_raw_parts_mut(wavelet, (2 * (nt / 2i32) - 1) as usize),
                    );
                } else {
                    if (wavelet_type == 4i32) {
                        oxwave_waveforms::akb_wavelet(
                            dt,
                            fpeak,
                            std::slice::from_raw_parts_mut(wavelet, (nt) as usize),
                        );
                    } else {
                        oxwave_core::err!("wavelet_type has to be 1,2,3 or 4");
                    }
                }
            }
        }

        /* loop over traces to convolve with wavelet */
        {
            ix = 0i32;
            while (ix < nx) {
                {
                    /* save input reflectivity trace */
                    {
                        it = 0i32;
                        while (it < nt) {
                            (*checked_at(temp, (it) as isize)) =
                                (*checked_at((*checked_at(wfield, (ix) as isize)), (it) as isize));
                            it += 1;
                        }
                    }
                    /* convolve with wavelet */
                    oxwave_f_convolution::convolve_cwp(
                        0i32,
                        std::slice::from_raw_parts(temp, (nt) as usize),
                        0i32,
                        std::slice::from_raw_parts(wavelet, (nt) as usize),
                        0i32,
                        std::slice::from_raw_parts_mut(
                            (*checked_at(wfield, (ix) as isize)),
                            (nt) as usize,
                        ),
                    );
                }
                ix += 1;
            }
        }
        /* free allocated space */
        ();
        ();
    }
}

// Refl/lib/psvreflect.c: psv_reflectivities
pub(super) unsafe fn raw_psv_reflectivities(
    mut int_type: i32,
    mut verbose: i32,
    mut wtype: i32,
    mut nw: i32,
    mut nlayers: i32,
    mut nx: i32,
    mut layern: i32,
    mut nor: i32,
    mut np: i32,
    mut bp1: f32,
    mut m1: f32,
    mut m2: f32,
    mut m3: f32,
    mut h1: f32,
    mut h2: f32,
    mut lsource: i32,
    mut bx: f32,
    mut dx: f32,
    mut xmax: f32,
    mut decay: f32,
    mut fref: f32,
    mut wrefp: f32,
    mut wrefs: f32,
    mut tsec: f32,
    mut p2w: f32,
    mut fs: f32,
    mut epsp: f32,
    mut epss: f32,
    mut sigp: f32,
    mut sigs: f32,
    mut pw1: f32,
    mut pw2: f32,
    mut pw3: f32,
    mut pw4: f32,
    mut acoustic: *mut i32,
    mut flag: *mut i32,
    mut lobs: *mut i32,
    mut rho: *mut f32,
    mut t: *mut f32,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut ql: *mut f32,
    mut qt: *mut f32,
    mut response1: *mut *mut *mut Complex,
    mut response2: *mut *mut *mut Complex,
    mut response3: *mut *mut *mut Complex,
    mut outfp: *mut CFile,
) {
    let mut arena = Arena::default();
    {
        let mut iw: i32 = 0;
        let mut iz: i32 = 0;
        let mut ip: i32 = 0;
        let mut ix: i32 = 0;
        /* loop counters */
        let mut ijk1: i32 = 0i32;
        let mut ik1: i32 = 0;
        let mut ik2: i32 = 0;
        let mut il: i32 = 0;
        let mut jl: i32 = 0;
        /* auxiliary indices */
        let mut ijk: i32 = 0i32;
        let mut ijk2: i32 = 0;
        let mut jj: i32 = 0;
        /* more auxiliary indices */
        /* int flg1,lrec; */
        let mut lrec: i32 = 0;
        let mut bp: f32 = (0f64 as f32);
        /* smallest ray parameter */
        let mut fp: f32 = 0.0;
        /* final ray parameter */
        let mut dp: f32 = 0.0;
        /* ray parameter increment */
        let mut nblock: i32 = 0;
        /* number of blocks to process */
        let mut block_size: i32 = 0;
        /* size of current block */
        let mut left: i32 = 0;
        /* remainder */
        let mut prv: *mut i32 = std::ptr::null_mut();
        let mut ibl: *mut i32 = std::ptr::null_mut();
        /* scratch arrays */
        let mut x1: f32 = 0.0;
        let mut w: f32 = 0.0;
        let mut wbeg: f32 = (0f64 as f32);
        let mut wfin: f32 = (0f64 as f32);
        /* aux variables for output */
        let mut cdp: Complex = cmplx(0.0, 0.0);
        let mut wpie: Complex = cmplx(0.0, 0.0);
        /* complex cdp and frequency */
        /* complex auxiliary variables */
        let mut divfac: Complex = cmplx(0.0, 0.0);
        let mut g: Complex = cmplx(0.0, 0.0);
        let mut h: Complex = cmplx(0.0, 0.0);
        let mut e: Complex = cmplx(0.0, 0.0);
        let mut f: Complex = cmplx(0.0, 0.0);
        let mut mu: Complex = cmplx(0.0, 0.0);
        let mut x: Complex = cmplx(0.0, 0.0);
        let mut y: Complex = cmplx(0.0, 0.0);
        let mut z: Complex = cmplx(0.0, 0.0);
        let mut a: Complex = cmplx(0.0, 0.0);
        let mut b: Complex = cmplx(0.0, 0.0);
        let mut c: Complex = cmplx(0.0, 0.0);
        let mut d: Complex = cmplx(0.0, 0.0);
        let mut t1: Complex = cmplx(0.0, 0.0);
        let mut tem: Complex = cmplx(0.0, 0.0);
        let mut temr: Complex = cmplx(0.0, 0.0);
        let mut sig: Complex = cmplx(0.0, 0.0);
        let mut sigh: Complex = cmplx(0.0, 0.0);
        let mut rho1: Complex = cmplx(0.0, 0.0);
        let mut rho2: Complex = cmplx(0.0, 0.0);
        let mut det: Complex = cmplx(0.0, 0.0);
        let mut psq: Complex = cmplx(0.0, 0.0);
        let mut gl1: Complex = cmplx(0.0, 0.0);
        let mut gt1: Complex = cmplx(0.0, 0.0);
        let mut gl2: Complex = cmplx(0.0, 0.0);
        let mut gt2: Complex = cmplx(0.0, 0.0);
        let mut gam1: Complex = cmplx(0.0, 0.0);
        let mut gam2: Complex = cmplx(0.0, 0.0);
        let mut at1: Complex = cmplx(0.0, 0.0);
        let mut at2: Complex = cmplx(0.0, 0.0);
        let mut al1: Complex = cmplx(0.0, 0.0);
        let mut al2: Complex = cmplx(0.0, 0.0);
        let mut rvrb111: Complex = cmplx(0.0, 0.0);
        let mut rvrb112: Complex = cmplx(0.0, 0.0);
        let mut rvrb121: Complex = cmplx(0.0, 0.0);
        let mut rvrb122: Complex = cmplx(0.0, 0.0);
        let mut rvrb211: Complex = cmplx(0.0, 0.0);
        let mut rvrb212: Complex = cmplx(0.0, 0.0);
        let mut rvrb221: Complex = cmplx(0.0, 0.0);
        let mut ewigh11: Complex = cmplx(0.0, 0.0);
        let mut ewigh22: Complex = cmplx(0.0, 0.0);
        let mut rvrb222: Complex = cmplx(0.0, 0.0);
        let mut wrn011: Complex = cmplx(0.0, 0.0);
        let mut wrn012: Complex = cmplx(0.0, 0.0);
        let mut wrn021: Complex = cmplx(0.0, 0.0);
        let mut wrn022: Complex = cmplx(0.0, 0.0);
        let mut ewrd11: Complex = cmplx(0.0, 0.0);
        let mut ewrd12: Complex = cmplx(0.0, 0.0);
        let mut ewrd21: Complex = cmplx(0.0, 0.0);
        let mut ewrd22: Complex = cmplx(0.0, 0.0);
        let mut ewd0211: Complex = cmplx(0.0, 0.0);
        let mut ewd0212: Complex = cmplx(0.0, 0.0);
        let mut ewd0221: Complex = cmplx(0.0, 0.0);
        let mut ewd0222: Complex = cmplx(0.0, 0.0);
        let mut ewtu211: Complex = cmplx(0.0, 0.0);
        let mut ewtu212: Complex = cmplx(0.0, 0.0);
        let mut ewtu221: Complex = cmplx(0.0, 0.0);
        let mut ewtu222: Complex = cmplx(0.0, 0.0);
        let mut rtb111: Complex = cmplx(0.0, 0.0);
        let mut rtb112: Complex = cmplx(0.0, 0.0);
        let mut rtb121: Complex = cmplx(0.0, 0.0);
        let mut rtb122: Complex = cmplx(0.0, 0.0);
        let mut rtb211: Complex = cmplx(0.0, 0.0);
        let mut rtb212: Complex = cmplx(0.0, 0.0);
        let mut rtb221: Complex = cmplx(0.0, 0.0);
        let mut rtb222: Complex = cmplx(0.0, 0.0);
        let mut rtb311: Complex = cmplx(0.0, 0.0);
        let mut rtb312: Complex = cmplx(0.0, 0.0);
        let mut rtb321: Complex = cmplx(0.0, 0.0);
        let mut rtb322: Complex = cmplx(0.0, 0.0);
        let mut func1: Complex = cmplx(0.0, 0.0);
        let mut func2: Complex = cmplx(0.0, 0.0);
        let mut func3: Complex = cmplx(0.0, 0.0);
        let mut func4: Complex = cmplx(0.0, 0.0);
        let mut func5: Complex = cmplx(0.0, 0.0);
        let mut func6: Complex = cmplx(0.0, 0.0);
        let mut tun0p11: Complex = cmplx(0.0, 0.0);
        let mut tun0p12: Complex = cmplx(0.0, 0.0);
        let mut tun0p21: Complex = cmplx(0.0, 0.0);
        let mut tun0p22: Complex = cmplx(0.0, 0.0);
        let mut rd0np11: Complex = cmplx(0.0, 0.0);
        let mut rd0np12: Complex = cmplx(0.0, 0.0);
        let mut rd0np21: Complex = cmplx(0.0, 0.0);
        let mut rd0np22: Complex = cmplx(0.0, 0.0);
        let mut td0np11: Complex = cmplx(0.0, 0.0);
        let mut td0np12: Complex = cmplx(0.0, 0.0);
        let mut td0np21: Complex = cmplx(0.0, 0.0);
        let mut td0np22: Complex = cmplx(0.0, 0.0);
        let mut rd11: Complex = cmplx(0.0, 0.0);
        let mut rd12: Complex = cmplx(0.0, 0.0);
        let mut rd21: Complex = cmplx(0.0, 0.0);
        let mut rd22: Complex = cmplx(0.0, 0.0);
        let mut td11: Complex = cmplx(0.0, 0.0);
        let mut td12: Complex = cmplx(0.0, 0.0);
        let mut td21: Complex = cmplx(0.0, 0.0);
        let mut td22: Complex = cmplx(0.0, 0.0);
        let mut tu11: Complex = cmplx(0.0, 0.0);
        let mut tu12: Complex = cmplx(0.0, 0.0);
        let mut tu21: Complex = cmplx(0.0, 0.0);
        let mut tu22: Complex = cmplx(0.0, 0.0);
        let mut ru11: Complex = cmplx(0.0, 0.0);
        let mut ru12: Complex = cmplx(0.0, 0.0);
        let mut ru21: Complex = cmplx(0.0, 0.0);
        let mut ru22: Complex = cmplx(0.0, 0.0);
        let mut ruf11: Complex = cmplx(0.0, 0.0);
        let mut ruf12: Complex = cmplx(0.0, 0.0);
        let mut ruf21: Complex = cmplx(0.0, 0.0);
        let mut ruf22: Complex = cmplx(0.0, 0.0);
        let mut fact11: Complex = cmplx(0.0, 0.0);
        let mut fact12: Complex = cmplx(0.0, 0.0);
        let mut vuzs1: Complex = cmplx(0.0, 0.0);
        let mut vuzs2: Complex = cmplx(0.0, 0.0);
        let mut vdzs1: Complex = cmplx(0.0, 0.0);
        let mut vdzs2: Complex = cmplx(0.0, 0.0);
        /* complex scratch arrays */
        let mut rd0n11: *mut Complex = std::ptr::null_mut();
        let mut rd0n12: *mut Complex = std::ptr::null_mut();
        let mut rd0n21: *mut Complex = std::ptr::null_mut();
        let mut rd0n22: *mut Complex = std::ptr::null_mut();
        let mut td0n11: *mut Complex = std::ptr::null_mut();
        let mut td0n12: *mut Complex = std::ptr::null_mut();
        let mut td0n21: *mut Complex = std::ptr::null_mut();
        let mut td0n22: *mut Complex = std::ptr::null_mut();
        let mut tun011: *mut Complex = std::ptr::null_mut();
        let mut tun012: *mut Complex = std::ptr::null_mut();
        let mut tun021: *mut Complex = std::ptr::null_mut();
        let mut tun022: *mut Complex = std::ptr::null_mut();
        let mut run011: *mut Complex = std::ptr::null_mut();
        let mut run012: *mut Complex = std::ptr::null_mut();
        let mut run021: *mut Complex = std::ptr::null_mut();
        let mut run022: *mut Complex = std::ptr::null_mut();
        let mut tusf11: *mut Complex = std::ptr::null_mut();
        let mut tusf12: *mut Complex = std::ptr::null_mut();
        let mut tusf21: *mut Complex = std::ptr::null_mut();
        let mut tusf22: *mut Complex = std::ptr::null_mut();
        let mut rdfs11: *mut Complex = std::ptr::null_mut();
        let mut rdfs12: *mut Complex = std::ptr::null_mut();
        let mut rdfs21: *mut Complex = std::ptr::null_mut();
        let mut rdfs22: *mut Complex = std::ptr::null_mut();
        let mut tdfs11: *mut Complex = std::ptr::null_mut();
        let mut tdfs12: *mut Complex = std::ptr::null_mut();
        let mut tdfs21: *mut Complex = std::ptr::null_mut();
        let mut tdfs22: *mut Complex = std::ptr::null_mut();
        let mut rurf11: *mut Complex = std::ptr::null_mut();
        let mut rurf12: *mut Complex = std::ptr::null_mut();
        let mut rurf21: *mut Complex = std::ptr::null_mut();
        let mut rurf22: *mut Complex = std::ptr::null_mut();
        let mut turf11: *mut Complex = std::ptr::null_mut();
        let mut turf12: *mut Complex = std::ptr::null_mut();
        let mut turf21: *mut Complex = std::ptr::null_mut();
        let mut turf22: *mut Complex = std::ptr::null_mut();
        let mut rdfr11: *mut Complex = std::ptr::null_mut();
        let mut rdfr12: *mut Complex = std::ptr::null_mut();
        let mut rdfr21: *mut Complex = std::ptr::null_mut();
        let mut rdfr22: *mut Complex = std::ptr::null_mut();
        let mut tdfr11: *mut Complex = std::ptr::null_mut();
        let mut tdfr12: *mut Complex = std::ptr::null_mut();
        let mut tdfr21: *mut Complex = std::ptr::null_mut();
        let mut tdfr22: *mut Complex = std::ptr::null_mut();
        let mut rnr11: *mut Complex = std::ptr::null_mut();
        let mut rnr12: *mut Complex = std::ptr::null_mut();
        let mut rnr21: *mut Complex = std::ptr::null_mut();
        let mut rnr22: *mut Complex = std::ptr::null_mut();
        let mut rdsr11: *mut Complex = std::ptr::null_mut();
        let mut rdsr12: *mut Complex = std::ptr::null_mut();
        let mut rdsr21: *mut Complex = std::ptr::null_mut();
        let mut rdsr22: *mut Complex = std::ptr::null_mut();
        /* complex *r11=NULL,*r12=NULL,*r21=NULL,*r22=NULL; */
        let mut texp: *mut Complex = std::ptr::null_mut();
        let mut vos1: *mut Complex = std::ptr::null_mut();
        let mut vos2: *mut Complex = std::ptr::null_mut();
        let mut vos3: *mut Complex = std::ptr::null_mut();
        let mut ux: *mut Complex = std::ptr::null_mut();
        let mut uz: *mut Complex = std::ptr::null_mut();
        let mut up: *mut Complex = std::ptr::null_mut();
        let mut alpha: *mut Complex = std::ptr::null_mut();
        let mut betha: *mut Complex = std::ptr::null_mut();
        let mut sigmad1: *mut Complex = std::ptr::null_mut();
        let mut sigmad2: *mut Complex = std::ptr::null_mut();
        let mut sigmau1: *mut Complex = std::ptr::null_mut();
        let mut sigmau2: *mut Complex = std::ptr::null_mut();
        let mut al: *mut Complex = std::ptr::null_mut();
        let mut at: *mut Complex = std::ptr::null_mut();
        let mut gl: *mut Complex = std::ptr::null_mut();
        let mut gt: *mut Complex = std::ptr::null_mut();
        let mut gam: *mut Complex = std::ptr::null_mut();
        let mut p: *mut Complex = std::ptr::null_mut();
        let mut pp: *mut Complex = std::ptr::null_mut();
        let mut prs: *mut Complex = std::ptr::null_mut();
        let mut pwin: *mut Complex = std::ptr::null_mut();
        /* complex pointers */
        let mut rusf11: *mut Complex = std::ptr::null_mut();
        let mut rusf12: *mut Complex = std::ptr::null_mut();
        let mut rusf21: *mut Complex = std::ptr::null_mut();
        let mut rusf22: *mut Complex = std::ptr::null_mut();
        let mut tusr11: *mut Complex = std::ptr::null_mut();
        let mut tusr12: *mut Complex = std::ptr::null_mut();
        let mut tusr21: *mut Complex = std::ptr::null_mut();
        let mut tusr22: *mut Complex = std::ptr::null_mut();
        let mut rusr11: *mut Complex = std::ptr::null_mut();
        let mut rusr12: *mut Complex = std::ptr::null_mut();
        let mut rusr21: *mut Complex = std::ptr::null_mut();
        let mut rusr22: *mut Complex = std::ptr::null_mut();
        let mut tdsr11: *mut Complex = std::ptr::null_mut();
        let mut tdsr12: *mut Complex = std::ptr::null_mut();
        let mut tdsr21: *mut Complex = std::ptr::null_mut();
        let mut tdsr22: *mut Complex = std::ptr::null_mut();
        let mut rdnr11: *mut Complex = std::ptr::null_mut();
        let mut rdnr12: *mut Complex = std::ptr::null_mut();
        let mut rdnr21: *mut Complex = std::ptr::null_mut();
        let mut rdnr22: *mut Complex = std::ptr::null_mut();
        /* allocate working space */
        ux = arena.alloc::<Complex>((100i32 as usize) as usize);
        uz = arena.alloc::<Complex>((100i32 as usize) as usize);
        up = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun011 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun012 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun021 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun022 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run011 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run012 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run021 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run022 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tdfs11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tdfs12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tdfs21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tdfs22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rurf11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rurf12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rurf21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rurf22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rnr11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rnr12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rnr21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rnr22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdsr11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdsr12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdsr21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdsr22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        vos1 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rusf21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rusf22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        vos2 = arena.alloc::<Complex>((100i32 as usize) as usize);
        vos3 = arena.alloc::<Complex>((100i32 as usize) as usize);
        /*	r11=alloc1complex(BLK);r12=alloc1complex(BLK);
            r21=alloc1complex(BLK);r22=alloc1complex(BLK);
        */
        texp = arena.alloc::<Complex>((100i32 as usize) as usize);
        /* allocate other working space */
        rdnr11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdnr12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdnr21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdnr22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        prv = arena.alloc::<i32>((nor as usize) as usize);
        prs = arena.alloc::<Complex>((nor as usize) as usize);
        ibl = arena.alloc::<i32>((100i32 as usize) as usize);
        p = arena.alloc::<Complex>((100i32 as usize) as usize);
        pp = arena.alloc::<Complex>((100i32 as usize) as usize);
        pwin = arena.alloc::<Complex>((100i32 as usize) as usize);
        al = arena.alloc::<Complex>((nlayers as usize) as usize);
        at = arena.alloc::<Complex>((nlayers as usize) as usize);
        gl = arena.alloc::<Complex>(((100i32 * nlayers) as usize) as usize);
        gt = arena.alloc::<Complex>(((100i32 * nlayers) as usize) as usize);
        gam = arena.alloc::<Complex>(((100i32 * nlayers) as usize) as usize);
        alpha = arena.alloc::<Complex>((100i32 as usize) as usize);
        betha = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmau1 = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmau2 = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmad1 = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmad2 = arena.alloc::<Complex>((100i32 as usize) as usize);
        /* set pointers */
        rusf11 = alpha;
        rusf12 = betha;
        tusr11 = tun011;
        tusr12 = tun012;
        tusr21 = tun021;
        tusr22 = tun022;
        rusr11 = run011;
        rusr12 = run012;
        rusr21 = run021;
        rusr22 = run022;
        tdsr11 = td0n11;
        tdsr12 = td0n12;
        tdsr21 = td0n21;
        tdsr22 = td0n22;
        /*
            rdnr11=r11;
            rdnr12=r12;
            rdnr21=r21;
            rdnr22=r22;
        */
        /* assign initial value to working variables */
        tem = cdiv(
            cmplx((1f64 as f32), (0f64 as f32)),
            cwp_cexp(crmul(
                cmplx((0f64 as f32), (1f64 as f32)),
                (((3.1415926535897931f64) / (4i32 as f64)) as f32),
            )),
        );
        ewigh11 = cmplx((0f64 as f32), (0f64 as f32));
        ewigh22 = cmplx((0f64 as f32), (0f64 as f32));
        /* initialize complex output arrays */
        {
            iw = 0i32;
            while (iw < nw) {
                {
                    ix = 0i32;
                    while (ix < nx) {
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    (*checked_at(
                                        (*checked_at(
                                            (*checked_at(response1, (iw) as isize)),
                                            (ix) as isize,
                                        )),
                                        (iz) as isize,
                                    )) = cmplx((0f64 as f32), (0f64 as f32));
                                    (*checked_at(
                                        (*checked_at(
                                            (*checked_at(response2, (iw) as isize)),
                                            (ix) as isize,
                                        )),
                                        (iz) as isize,
                                    )) = cmplx((0f64 as f32), (0f64 as f32));
                                    (*checked_at(
                                        (*checked_at(
                                            (*checked_at(response3, (iw) as isize)),
                                            (ix) as isize,
                                        )),
                                        (iz) as isize,
                                    )) = cmplx((0f64 as f32), (0f64 as f32));
                                }
                                iz += 1;
                            }
                        }
                        ix += 1;
                    }
                }
                iw += 1;
            }
        }
        /**************************************************************************
        	*						Main loop over frequencies			*
        	**************************************************************************/
        {
            iw = 0i32;
            while (iw < nw) {
                {
                    w = ((((((iw + 1i32) as f64) * (3.1415926535897931f64)) * (2i32 as f64))
                        / (tsec as f64)) as f32);
                    /* w in radians/sec */
                    if (iw == 0i32) {
                        wbeg = w;
                    }

                    if (iw == (nw - 1i32)) {
                        wfin = w;
                    }

                    wpie = cmplx(w, decay);
                    /* complex w */
                    /* compute frequency-dependent auxiliary variables for reflectivity */
                    raw_compute_w_aux_arrays(
                        wtype,
                        layern,
                        nlayers,
                        nor,
                        lsource,
                        (&raw mut np),
                        (&raw mut block_size),
                        (&raw mut nblock),
                        (&raw mut left),
                        fref,
                        wrefp,
                        wrefs,
                        lobs,
                        tsec,
                        p2w,
                        fs,
                        xmax,
                        w,
                        decay,
                        epsp,
                        epss,
                        sigp,
                        sigs,
                        (&raw mut pw1),
                        (&raw mut pw2),
                        (&raw mut pw3),
                        (&raw mut pw4),
                        (&raw mut dp),
                        (&raw mut fp),
                        wpie,
                        (&raw mut cdp),
                        cl,
                        ct,
                        ql,
                        qt,
                        rho,
                        t,
                        al,
                        at,
                        prs,
                    );
                    /* if requested, output processing information */
                    if ((verbose == 1i32) || (verbose == 3i32)) {
                        oxwave_core::fprintf!(
                            &oxwave_core::prelude::stderr(),
                            "%3d%6.1f%7.3f%7.3f%6d%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                            iw,
                            (w as f64),
                            (bp as f64),
                            (fp as f64),
                            np,
                            (dp as f64),
                            (pw1 as f64),
                            (pw2 as f64),
                            (pw3 as f64),
                            (pw4 as f64)
                        );
                    }

                    if ((verbose == 2i32) || (verbose == 3i32)) {
                        oxwave_core::fprintf!(
                            &*outfp,
                            "%3d%6.1f%7.3f%7.3f%6d%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                            iw,
                            (w as f64),
                            (bp as f64),
                            (fp as f64),
                            np,
                            (dp as f64),
                            (pw1 as f64),
                            (pw2 as f64),
                            (pw3 as f64),
                            (pw4 as f64)
                        );
                    }

                    /**********************************************************************
                    	* 		loop over slopes to compute the reflectivities (by blocks)    *
                     **********************************************************************/
                    while (1i32 != 0) {
                        if ((nblock == 1i32) && (left != 0i32)) {
                            block_size = left;
                        }

                        /* compute p_aux_arrays */
                        raw_compute_p_aux_arrays(
                            wtype,
                            nlayers,
                            lsource,
                            block_size,
                            bp,
                            dp,
                            w,
                            decay,
                            pw1,
                            pw2,
                            pw3,
                            pw4,
                            pwin,
                            m1,
                            m2,
                            m3,
                            h1,
                            h2,
                            wpie,
                            (&raw mut divfac),
                            p,
                            pp,
                            al,
                            at,
                            rho,
                            gl,
                            gt,
                            gam,
                            alpha,
                            betha,
                            sigmad1,
                            sigmad2,
                            sigmau1,
                            sigmau2,
                        );
                        /******************************************************************
                        	*			start computation of reflectivity series			  *
                        	******************************************************************/
                        /* loop over reflecting layers */
                        ijk = 0i32;
                        {
                            il = 0i32;
                            while (il < (nlayers - 1i32)) {
                                {
                                    jl = (il + 1i32);
                                    ik1 = (il * block_size);
                                    ik2 = ((il + 1i32) * block_size);
                                    al1 = (*checked_at(al, (il) as isize));
                                    al2 = (*checked_at(al, (jl) as isize));
                                    at1 = (*checked_at(at, (il) as isize));
                                    at2 = (*checked_at(at, (jl) as isize));
                                    rho1 = cmplx((*checked_at(rho, (il) as isize)), (0f64 as f32));
                                    rho2 = cmplx((*checked_at(rho, (jl) as isize)), (0f64 as f32));
                                    t1 = cmplx((*checked_at(t, (il) as isize)), (0f64 as f32));
                                    if (il == 0i32) {
                                        /* initialize, first layer only */
                                        if ((((at1).r as f64) == 0f64)
                                            && (((at1).i as f64) == 0f64))
                                        {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        gl1 = (*checked_at(gl, (ijk1) as isize));
                                                        ewigh11 = cwp_cexp(cmul(
                                                            cmplx((0f64 as f32), (1f64 as f32)),
                                                            cmul(wpie, cmul(gl1, t1)),
                                                        ));
                                                        (*checked_at(rd0n11, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n11, (ip) as isize)) =
                                                            ewigh11;
                                                        (*checked_at(td0n12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun011, (ip) as isize)) =
                                                            ewigh11;
                                                        (*checked_at(tun012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run011, (ip) as isize)) =
                                                            cneg(cmul(ewigh11, ewigh11));
                                                        (*checked_at(run012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        } else {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        psq = (*checked_at(pp, (ip) as isize));
                                                        gl1 = (*checked_at(gl, (ijk1) as isize));
                                                        gt1 = (*checked_at(gt, (ijk1) as isize));
                                                        gam1 = (*checked_at(gam, (ijk1) as isize));
                                                        mu = cmul(
                                                            cdiv(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                cmul(at1, at1),
                                                            ),
                                                            rho1,
                                                        );
                                                        rtb211 = crmul(
                                                            cmul(
                                                                (*checked_at(p, (ip) as isize)),
                                                                cmul(mu, gl1),
                                                            ),
                                                            (2f64 as f32),
                                                        );
                                                        rtb212 = gam1;
                                                        rtb221 = cneg(gam1);
                                                        rtb222 = crmul(
                                                            cmul(
                                                                (*checked_at(p, (ip) as isize)),
                                                                cmul(mu, gt1),
                                                            ),
                                                            (2f64 as f32),
                                                        );
                                                        det = csub(
                                                            cmul(rtb211, rtb222),
                                                            cmul(rtb212, rtb221),
                                                        );
                                                        rtb111 = cdiv(rtb222, det);
                                                        rtb112 = cdiv(rtb212, det);
                                                        rtb121 = cneg(cdiv(rtb221, det));
                                                        rtb122 = cdiv(rtb211, det);
                                                        ruf11 = cadd(
                                                            cmul(rtb211, rtb222),
                                                            cmul(rtb212, rtb221),
                                                        );
                                                        ruf12 = cmul(
                                                            rtb222,
                                                            crmul(rtb212, (2f64 as f32)),
                                                        );
                                                        ruf21 = cmul(
                                                            rtb211,
                                                            crmul(rtb221, (2f64 as f32)),
                                                        );
                                                        ruf22 = ruf11;
                                                        ruf11 = cdiv(ruf11, det);
                                                        ruf12 = cdiv(ruf12, det);
                                                        ruf21 = cdiv(ruf21, det);
                                                        ruf22 = cdiv(ruf22, det);
                                                        ewigh11 = cwp_cexp(cmul(
                                                            gl1,
                                                            cmul(
                                                                wpie,
                                                                cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    t1,
                                                                ),
                                                            ),
                                                        ));
                                                        ewigh22 = cwp_cexp(cmul(
                                                            gt1,
                                                            cmul(
                                                                wpie,
                                                                cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    t1,
                                                                ),
                                                            ),
                                                        ));
                                                        (*checked_at(rd0n11, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n11, (ip) as isize)) =
                                                            ewigh11;
                                                        (*checked_at(td0n12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n22, (ip) as isize)) =
                                                            ewigh22;
                                                        (*checked_at(tun011, (ip) as isize)) =
                                                            ewigh11;
                                                        (*checked_at(tun012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun022, (ip) as isize)) =
                                                            ewigh22;
                                                        rtb111 = cmul(ewigh11, ruf11);
                                                        rtb112 = cmul(ewigh11, ruf12);
                                                        rtb121 = cmul(ewigh22, ruf21);
                                                        rtb122 = cmul(ewigh22, ruf22);
                                                        (*checked_at(run011, (ip) as isize)) =
                                                            cmul(ewigh11, rtb111);
                                                        (*checked_at(run012, (ip) as isize)) =
                                                            cmul(ewigh22, rtb112);
                                                        (*checked_at(run021, (ip) as isize)) =
                                                            cmul(ewigh11, rtb121);
                                                        (*checked_at(run022, (ip) as isize)) =
                                                            cmul(ewigh22, rtb122);
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        }
                                    } else {
                                        /* layers after the first one */
                                        /* case 1: liquid-liquid */
                                        if ((((((at1).r as f64) == 0f64)
                                            && (((at1).i as f64) == 0f64))
                                            && (((at2).r as f64) == 0f64))
                                            && (((at2).i as f64) == 0f64))
                                        {
                                            if (((al1).r == (al2).r) && ((al1).i == (al2).i)) {
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            gl1 =
                                                                (*checked_at(gl, (ijk1) as isize));
                                                            ewigh11 = cwp_cexp(cmul(
                                                                wpie,
                                                                cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    cmul(t1, gl1),
                                                                ),
                                                            ));
                                                            tun0p11 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    tun011,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            tun0p21 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    tun021,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb111 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    run011,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            (*checked_at(run011, (ip) as isize)) =
                                                                cmul(ewigh11, rtb111);
                                                            (*checked_at(run012, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(run021, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(run022, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            td0np11 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    td0n11,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            td0np12 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    td0n12,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            (*checked_at(tun011, (ip) as isize)) =
                                                                tun0p11;
                                                            (*checked_at(tun012, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(tun021, (ip) as isize)) =
                                                                tun0p21;
                                                            (*checked_at(tun022, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(td0n11, (ip) as isize)) =
                                                                td0np11;
                                                            (*checked_at(td0n12, (ip) as isize)) =
                                                                td0np12;
                                                            (*checked_at(td0n21, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(td0n22, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            } else {
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            ijk2 = (ik2 + ip);
                                                            psq = (*checked_at(pp, (ip) as isize));
                                                            gl1 =
                                                                (*checked_at(gl, (ijk1) as isize));
                                                            gt1 =
                                                                (*checked_at(gt, (ijk1) as isize));
                                                            gam1 =
                                                                (*checked_at(gam, (ijk1) as isize));
                                                            gam2 =
                                                                (*checked_at(gam, (ijk2) as isize));
                                                            gl2 =
                                                                (*checked_at(gl, (ijk2) as isize));
                                                            gt2 =
                                                                (*checked_at(gt, (ijk2) as isize));
                                                            ewigh11 = cwp_cexp(cmul(
                                                                wpie,
                                                                cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    cmul(t1, gl1),
                                                                ),
                                                            ));
                                                            a = csub(gam2, gam1);
                                                            b = cadd(rho1, a);
                                                            c = csub(rho2, a);
                                                            e = cadd(cmul(b, gl1), cmul(c, gl2));
                                                            y = cadd(
                                                                cmul(gl2, rho1),
                                                                cmul(gl1, rho2),
                                                            );
                                                            td11 = cmul(
                                                                cdiv(gl1, y),
                                                                crmul(rho1, (2f64 as f32)),
                                                            );
                                                            rd11 = cdiv(
                                                                csub(
                                                                    cmul(gl1, rho2),
                                                                    cmul(gl2, rho1),
                                                                ),
                                                                y,
                                                            );
                                                            tu11 = cmul(
                                                                cdiv(gl2, y),
                                                                crmul(rho2, (2f64 as f32)),
                                                            );
                                                            ru11 = cneg(rd11);
                                                            ewrd11 = cmul(ewigh11, rd11);
                                                            wrn011 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    run011,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rvrb111 = cdiv(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                csub(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    cmul(wrn011, ewrd11),
                                                                ),
                                                            );
                                                            rvrb211 = cadd(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                cmul(rvrb111, cmul(ewrd11, wrn011)),
                                                            );
                                                            ewd0211 = cmul(
                                                                (*checked_at(
                                                                    td0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(rvrb111, ewigh11),
                                                            );
                                                            ewd0212 = cmul(
                                                                (*checked_at(
                                                                    td0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(rvrb111, ewigh11),
                                                            );
                                                            (*checked_at(td0n11, (ip) as isize)) =
                                                                cmul(td11, ewd0211);
                                                            (*checked_at(td0n12, (ip) as isize)) =
                                                                cmul(td11, ewd0212);
                                                            (*checked_at(td0n21, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(td0n22, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            ewtu211 =
                                                                cmul(rvrb211, cmul(ewigh11, tu11));
                                                            tun0p11 = cmul(
                                                                ewtu211,
                                                                (*checked_at(
                                                                    tun011,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            tun0p21 = cmul(
                                                                ewtu211,
                                                                (*checked_at(
                                                                    tun021,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb111 = cmul(ewrd11, ewd0211);
                                                            rtb112 = cmul(ewrd11, ewd0212);
                                                            rd0np11 = cadd(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun011,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb111,
                                                                ),
                                                            );
                                                            rd0np12 = cadd(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun011,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb112,
                                                                ),
                                                            );
                                                            rd0np21 = cadd(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun021,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb111,
                                                                ),
                                                            );
                                                            rd0np22 = cadd(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun022,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb112,
                                                                ),
                                                            );
                                                            (*checked_at(run011, (ip) as isize)) =
                                                                cadd(
                                                                    ru11,
                                                                    cmul(
                                                                        ewtu211,
                                                                        cmul(wrn011, td11),
                                                                    ),
                                                                );
                                                            (*checked_at(run012, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(run021, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(run022, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(tun011, (ip) as isize)) =
                                                                tun0p11;
                                                            (*checked_at(tun012, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(tun021, (ip) as isize)) =
                                                                tun0p21;
                                                            (*checked_at(tun022, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rd0n11, (ip) as isize)) =
                                                                rd0np11;
                                                            (*checked_at(rd0n12, (ip) as isize)) =
                                                                rd0np12;
                                                            (*checked_at(rd0n21, (ip) as isize)) =
                                                                rd0np21;
                                                            (*checked_at(rd0n22, (ip) as isize)) =
                                                                rd0np22;
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            }
                                        } else {
                                            if (((((at1).r as f64) == 0f64)
                                                && (((at1).i as f64) == 0f64))
                                                && ((((at2).r as f64) != 0f64)
                                                    || (((at2).i as f64) != 0f64)))
                                            {
                                                /* case 2: liquid-solid */
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            ijk2 = (ik2 + ip);
                                                            psq = (*checked_at(pp, (ip) as isize));
                                                            gl1 =
                                                                (*checked_at(gl, (ijk1) as isize));
                                                            gt1 =
                                                                (*checked_at(gt, (ijk1) as isize));
                                                            gam1 =
                                                                (*checked_at(gam, (ijk1) as isize));
                                                            gl2 =
                                                                (*checked_at(gl, (ijk2) as isize));
                                                            gt2 =
                                                                (*checked_at(gt, (ijk2) as isize));
                                                            gam2 =
                                                                (*checked_at(gam, (ijk2) as isize));
                                                            ewigh11 = cwp_cexp(cmul(
                                                                wpie,
                                                                cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    cmul(t1, gl1),
                                                                ),
                                                            ));
                                                            ewigh22 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            a = csub(gam2, gam1);
                                                            b = cadd(rho1, a);
                                                            c = csub(rho2, a);
                                                            e = cadd(cmul(gl1, b), cmul(gl2, c));
                                                            d = cdiv(
                                                                crmul(rho2, (2f64 as f32)),
                                                                cmul(at2, at2),
                                                            );
                                                            f = b;
                                                            g = csub(a, cmul(d, cmul(gl1, gt2)));
                                                            h = cneg(cmul(d, gl2));
                                                            z = cdiv(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                cadd(
                                                                    cmul(e, f),
                                                                    cmul(g, cmul(h, psq)),
                                                                ),
                                                            );
                                                            y = cmul(
                                                                cmul(z, gl1),
                                                                crmul(rho1, (2f64 as f32)),
                                                            );
                                                            td11 = cmul(y, f);
                                                            td12 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            td21 = cneg(cmul(
                                                                y,
                                                                cmul(
                                                                    h,
                                                                    (*checked_at(p, (ip) as isize)),
                                                                ),
                                                            ));
                                                            td22 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            rd11 = cmul(
                                                                z,
                                                                csub(
                                                                    cmul(
                                                                        f,
                                                                        csub(
                                                                            cmul(gl1, b),
                                                                            cmul(gl2, c),
                                                                        ),
                                                                    ),
                                                                    cmul(
                                                                        cmul(psq, h),
                                                                        cadd(
                                                                            a,
                                                                            cmul(gl1, cmul(gt2, d)),
                                                                        ),
                                                                    ),
                                                                ),
                                                            );
                                                            rd12 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            rd21 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            rd22 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            y = cmul(z, crmul(rho2, (2f64 as f32)));
                                                            tu11 = cmul(y, cmul(f, gl2));
                                                            tu12 = cmul(
                                                                (*checked_at(p, (ip) as isize)),
                                                                cmul(y, cmul(gt2, h)),
                                                            );
                                                            tu21 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            tu22 =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            y = crmul(
                                                                cmul(
                                                                    (*checked_at(p, (ip) as isize)),
                                                                    cmul(gl1, cmul(b, cmul(d, z))),
                                                                ),
                                                                (2f64 as f32),
                                                            );
                                                            ru11 = cmul(
                                                                z,
                                                                csub(
                                                                    cmul(
                                                                        f,
                                                                        csub(
                                                                            cmul(gl2, c),
                                                                            cmul(gl1, b),
                                                                        ),
                                                                    ),
                                                                    cmul(
                                                                        cmul(gl2, d),
                                                                        cmul(g, psq),
                                                                    ),
                                                                ),
                                                            );
                                                            ru12 = cmul(gt2, y);
                                                            ru21 = cneg(cmul(gl2, y));
                                                            ru22 = cneg(cmul(
                                                                z,
                                                                cadd(
                                                                    cmul(b, e),
                                                                    cmul(
                                                                        cmul(psq, h),
                                                                        cadd(
                                                                            a,
                                                                            cmul(gl1, cmul(gt2, d)),
                                                                        ),
                                                                    ),
                                                                ),
                                                            ));
                                                            ewrd11 = cmul(ewigh11, rd11);
                                                            ewrd12 = cmul(ewigh11, rd12);
                                                            ewrd21 = cmul(ewigh22, rd21);
                                                            ewrd22 = cmul(ewigh22, rd22);
                                                            wrn011 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    run011,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            wrn012 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    run012,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            wrn021 = cmul(
                                                                ewigh22,
                                                                (*checked_at(
                                                                    run021,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            wrn022 = cmul(
                                                                ewigh22,
                                                                (*checked_at(
                                                                    run022,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb111 = cadd(
                                                                cmul(wrn011, ewrd11),
                                                                cmul(wrn012, ewrd21),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(wrn011, ewrd12),
                                                                cmul(wrn012, ewrd22),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(wrn021, ewrd11),
                                                                cmul(wrn022, ewrd21),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(wrn021, ewrd12),
                                                                cmul(wrn022, ewrd22),
                                                            );
                                                            x = csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                rtb111,
                                                            );
                                                            y = csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                rtb122,
                                                            );
                                                            det = csub(
                                                                cmul(x, y),
                                                                cmul(rtb112, rtb121),
                                                            );
                                                            rvrb111 = cdiv(y, det);
                                                            rvrb112 = cdiv(rtb112, det);
                                                            rvrb121 = cdiv(rtb121, det);
                                                            rvrb122 = cdiv(x, det);
                                                            rtb111 = cadd(
                                                                cmul(rvrb111, wrn011),
                                                                cmul(rvrb112, wrn021),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(rvrb111, wrn012),
                                                                cmul(rvrb112, wrn022),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(rvrb121, wrn011),
                                                                cmul(rvrb122, wrn021),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(rvrb121, wrn012),
                                                                cmul(rvrb122, wrn022),
                                                            );
                                                            rtb211 = cadd(
                                                                cmul(ewrd11, rtb111),
                                                                cmul(ewrd12, rtb121),
                                                            );
                                                            rtb212 = cadd(
                                                                cmul(ewrd11, rtb112),
                                                                cmul(ewrd12, rtb122),
                                                            );
                                                            rtb221 = cadd(
                                                                cmul(ewrd21, rtb111),
                                                                cmul(ewrd22, rtb121),
                                                            );
                                                            rtb222 = cadd(
                                                                cmul(ewrd21, rtb112),
                                                                cmul(ewrd22, rtb122),
                                                            );
                                                            rvrb211 = cadd(
                                                                rtb211,
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                            );
                                                            rvrb212 = rtb212;
                                                            rvrb221 = rtb221;
                                                            rvrb222 = cadd(
                                                                rtb222,
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                            );
                                                            rtb111 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    td0n11,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb112 = cmul(
                                                                ewigh11,
                                                                (*checked_at(
                                                                    td0n12,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb121 = cmul(
                                                                ewigh22,
                                                                (*checked_at(
                                                                    td0n21,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb122 = cmul(
                                                                ewigh22,
                                                                (*checked_at(
                                                                    td0n22,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            ewd0211 = cadd(
                                                                cmul(rtb111, rvrb111),
                                                                cmul(rvrb112, rtb121),
                                                            );
                                                            ewd0212 = cadd(
                                                                cmul(rtb112, rvrb111),
                                                                cmul(rvrb112, rtb122),
                                                            );
                                                            ewd0221 = cadd(
                                                                cmul(rtb111, rvrb121),
                                                                cmul(rvrb122, rtb121),
                                                            );
                                                            ewd0222 = cadd(
                                                                cmul(rtb112, rvrb121),
                                                                cmul(rvrb122, rtb122),
                                                            );
                                                            (*checked_at(td0n11, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(td11, ewd0211),
                                                                    cmul(td12, ewd0221),
                                                                );
                                                            (*checked_at(td0n12, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(td11, ewd0212),
                                                                    cmul(td12, ewd0222),
                                                                );
                                                            (*checked_at(td0n21, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(td21, ewd0211),
                                                                    cmul(td22, ewd0221),
                                                                );
                                                            (*checked_at(td0n22, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(td21, ewd0212),
                                                                    cmul(td22, ewd0222),
                                                                );
                                                            rtb111 = cadd(ewigh11, tu11);
                                                            rtb112 = cadd(ewigh11, tu12);
                                                            rtb121 = cadd(ewigh22, tu21);
                                                            rtb122 = cadd(ewigh22, tu22);
                                                            ewtu211 = cadd(
                                                                cmul(rvrb211, rtb111),
                                                                cmul(rvrb212, rtb121),
                                                            );
                                                            ewtu212 = cadd(
                                                                cmul(rvrb211, rtb112),
                                                                cmul(rvrb212, rtb122),
                                                            );
                                                            ewtu221 = cadd(
                                                                cmul(rvrb221, rtb111),
                                                                cmul(rvrb222, rtb121),
                                                            );
                                                            ewtu222 = cadd(
                                                                cmul(rvrb221, rtb112),
                                                                cmul(rvrb222, rtb122),
                                                            );
                                                            tun0p11 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun011,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu211,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun012,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu221,
                                                                ),
                                                            );
                                                            tun0p12 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun011,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu212,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun012,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu222,
                                                                ),
                                                            );
                                                            tun0p21 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun021,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu211,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun022,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu221,
                                                                ),
                                                            );
                                                            tun0p22 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun021,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu212,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tun022,
                                                                        (ip) as isize,
                                                                    )),
                                                                    ewtu222,
                                                                ),
                                                            );
                                                            rtb111 = cadd(
                                                                cmul(ewrd11, ewd0211),
                                                                cmul(ewrd12, ewd0221),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(ewrd11, ewd0212),
                                                                cmul(ewrd12, ewd0222),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(ewrd21, ewd0211),
                                                                cmul(ewrd22, ewd0221),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(ewrd21, ewd0212),
                                                                cmul(ewrd22, ewd0222),
                                                            );
                                                            rtb211 = cadd(
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        tun011,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb121,
                                                                    (*checked_at(
                                                                        tun012,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb212 = cadd(
                                                                cmul(
                                                                    rtb112,
                                                                    (*checked_at(
                                                                        tun011,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb122,
                                                                    (*checked_at(
                                                                        tun012,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb221 = cadd(
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        tun021,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb121,
                                                                    (*checked_at(
                                                                        tun022,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb222 = cadd(
                                                                cmul(
                                                                    rtb112,
                                                                    (*checked_at(
                                                                        tun021,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb122,
                                                                    (*checked_at(
                                                                        tun022,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rd0np11 = cadd(
                                                                rtb211,
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rd0np12 = cadd(
                                                                rtb212,
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rd0np21 = cadd(
                                                                rtb221,
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rd0np22 = cadd(
                                                                rtb222,
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                            );
                                                            rtb111 = cadd(
                                                                cmul(wrn011, ewtu211),
                                                                cmul(ewtu221, wrn012),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(wrn011, ewtu212),
                                                                cmul(ewtu222, wrn012),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(wrn021, ewtu211),
                                                                cmul(ewtu221, wrn022),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(wrn021, ewtu212),
                                                                cmul(ewtu222, wrn022),
                                                            );
                                                            rtb211 = cadd(
                                                                cmul(td11, rtb111),
                                                                cmul(td12, rtb121),
                                                            );
                                                            rtb212 = cadd(
                                                                cmul(td11, rtb112),
                                                                cmul(td12, rtb122),
                                                            );
                                                            rtb221 = cadd(
                                                                cmul(td21, rtb111),
                                                                cmul(td22, rtb121),
                                                            );
                                                            rtb222 = cadd(
                                                                cmul(td21, rtb112),
                                                                cmul(td22, rtb122),
                                                            );
                                                            (*checked_at(run011, (ip) as isize)) =
                                                                cadd(ru11, rtb211);
                                                            (*checked_at(run012, (ip) as isize)) =
                                                                cadd(ru12, rtb212);
                                                            (*checked_at(run021, (ip) as isize)) =
                                                                cadd(ru21, rtb221);
                                                            (*checked_at(run022, (ip) as isize)) =
                                                                cadd(ru22, rtb222);
                                                            (*checked_at(rd0n11, (ip) as isize)) =
                                                                rd0np11;
                                                            (*checked_at(rd0n12, (ip) as isize)) =
                                                                rd0np12;
                                                            (*checked_at(rd0n21, (ip) as isize)) =
                                                                rd0np21;
                                                            (*checked_at(rd0n22, (ip) as isize)) =
                                                                rd0np22;
                                                            (*checked_at(tun011, (ip) as isize)) =
                                                                tun0p11;
                                                            (*checked_at(tun012, (ip) as isize)) =
                                                                tun0p12;
                                                            (*checked_at(tun021, (ip) as isize)) =
                                                                tun0p21;
                                                            (*checked_at(tun022, (ip) as isize)) =
                                                                tun0p22;
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            } else {
                                                if ((((((at1).r as f64) != 0f64)
                                                    && (((at1).i as f64) != 0f64))
                                                    && (((at2).r as f64) == 0f64))
                                                    && (((at2).i as f64) == 0f64))
                                                {
                                                    /* case 3: solid-liquid */
                                                    {
                                                        ip = 0i32;
                                                        while (ip < block_size) {
                                                            {
                                                                ijk1 = (ik1 + ip);
                                                                ijk2 = (ik2 + ip);
                                                                psq = (*checked_at(
                                                                    pp,
                                                                    (ip) as isize,
                                                                ));
                                                                gl1 = (*checked_at(
                                                                    gl,
                                                                    (ijk1) as isize,
                                                                ));
                                                                gt1 = (*checked_at(
                                                                    gt,
                                                                    (ijk1) as isize,
                                                                ));
                                                                gam1 = (*checked_at(
                                                                    gam,
                                                                    (ijk1) as isize,
                                                                ));
                                                                gl2 = (*checked_at(
                                                                    gl,
                                                                    (ik2) as isize,
                                                                ));
                                                                gt2 = (*checked_at(
                                                                    gt,
                                                                    (ik2) as isize,
                                                                ));
                                                                gam2 = (*checked_at(
                                                                    gam,
                                                                    (ijk2) as isize,
                                                                ));
                                                                ewigh11 = cwp_cexp(cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    cmul(wpie, cmul(gl1, t1)),
                                                                ));
                                                                ewigh22 = cwp_cexp(cmul(
                                                                    cmplx(
                                                                        (0f64 as f32),
                                                                        (1f64 as f32),
                                                                    ),
                                                                    cmul(wpie, cmul(gt1, t1)),
                                                                ));
                                                                a = csub(gam2, gam1);
                                                                b = cadd(rho1, a);
                                                                c = csub(rho2, a);
                                                                e = cadd(
                                                                    cmul(b, gl1),
                                                                    cmul(c, gl2),
                                                                );
                                                                d = cdiv(
                                                                    crmul(rho1, ((-2f64) as f32)),
                                                                    cmul(at1, at1),
                                                                );
                                                                f = c;
                                                                g = cneg(cmul(gl1, d));
                                                                h = csub(
                                                                    a,
                                                                    cmul(gl2, cmul(d, gt1)),
                                                                );
                                                                z = cdiv(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    cadd(
                                                                        cmul(e, f),
                                                                        cmul(g, cmul(h, psq)),
                                                                    ),
                                                                );
                                                                y = cmul(
                                                                    z,
                                                                    crmul(rho1, (2i32 as f32)),
                                                                );
                                                                td11 = cmul(y, cmul(gl1, f));
                                                                td12 = cmul(
                                                                    y,
                                                                    cmul(
                                                                        (*checked_at(
                                                                            p,
                                                                            (ip) as isize,
                                                                        )),
                                                                        cmul(gt1, g),
                                                                    ),
                                                                );
                                                                td21 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                td22 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                y = crmul(cmul(z,cmul(c,cmul(d,cmul((*checked_at(p,(ip) as isize)),gl2)))),(2f64 as f32));
                                                                rd11 = cmul(
                                                                    z,
                                                                    csub(
                                                                        cmul(
                                                                            f,
                                                                            csub(
                                                                                cmul(gl1, b),
                                                                                cmul(gl2, c),
                                                                            ),
                                                                        ),
                                                                        cmul(
                                                                            cmul(psq, h),
                                                                            cmul(gl1, d),
                                                                        ),
                                                                    ),
                                                                );
                                                                rd12 = cmul(y, gt1);
                                                                rd21 = cneg(cmul(gl1, y));
                                                                rd22 = cneg(cmul(
                                                                    z,
                                                                    cadd(
                                                                        cmul(c, e),
                                                                        cmul(
                                                                            cadd(
                                                                                a,
                                                                                cmul(
                                                                                    gt1,
                                                                                    cmul(gl2, d),
                                                                                ),
                                                                            ),
                                                                            cmul(psq, g),
                                                                        ),
                                                                    ),
                                                                ));
                                                                y = cmul(
                                                                    cmul(z, gl2),
                                                                    crmul(rho2, (2f64 as f32)),
                                                                );
                                                                tu11 = cmul(y, f);
                                                                tu12 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                tu21 = cneg(cmul(
                                                                    y,
                                                                    cmul(
                                                                        (*checked_at(
                                                                            p,
                                                                            (ip) as isize,
                                                                        )),
                                                                        y,
                                                                    ),
                                                                ));
                                                                tu22 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                ru11 = cmul(
                                                                    z,
                                                                    csub(
                                                                        cmul(
                                                                            f,
                                                                            csub(
                                                                                cmul(gl2, c),
                                                                                cmul(gl1, b),
                                                                            ),
                                                                        ),
                                                                        cmul(
                                                                            cadd(
                                                                                a,
                                                                                cmul(
                                                                                    gl2,
                                                                                    cmul(d, gt1),
                                                                                ),
                                                                            ),
                                                                            cmul(psq, g),
                                                                        ),
                                                                    ),
                                                                );
                                                                ru12 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                ru21 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                ru22 = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                ewrd11 = cmul(ewigh11, rd11);
                                                                ewrd12 = cmul(ewigh11, rd12);
                                                                ewrd21 = cmul(ewigh22, rd21);
                                                                ewrd22 = cmul(ewigh22, rd22);
                                                                wrn011 = cmul(
                                                                    ewigh11,
                                                                    (*checked_at(
                                                                        run011,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                wrn012 = cmul(
                                                                    ewigh11,
                                                                    (*checked_at(
                                                                        run012,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                wrn021 = cmul(
                                                                    ewigh22,
                                                                    (*checked_at(
                                                                        run021,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                wrn022 = cmul(
                                                                    ewigh22,
                                                                    (*checked_at(
                                                                        run022,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                rtb111 = cadd(
                                                                    cmul(wrn011, ewrd11),
                                                                    cmul(ewrd21, wrn012),
                                                                );
                                                                rtb112 = cadd(
                                                                    cmul(wrn011, ewrd12),
                                                                    cmul(ewrd22, wrn012),
                                                                );
                                                                rtb121 = cadd(
                                                                    cmul(wrn021, ewrd11),
                                                                    cmul(ewrd21, wrn022),
                                                                );
                                                                rtb122 = cadd(
                                                                    cmul(wrn021, ewrd12),
                                                                    cmul(ewrd22, wrn022),
                                                                );
                                                                x = csub(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    rtb111,
                                                                );
                                                                y = csub(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    rtb122,
                                                                );
                                                                det = csub(
                                                                    cmul(x, y),
                                                                    cmul(rtb112, rtb121),
                                                                );
                                                                rvrb111 = cdiv(y, det);
                                                                rvrb112 = cdiv(rtb112, det);
                                                                rvrb121 = cdiv(rtb121, det);
                                                                rvrb122 = cdiv(x, det);
                                                                rtb111 = cadd(
                                                                    cmul(rvrb111, wrn011),
                                                                    cmul(rvrb112, wrn021),
                                                                );
                                                                rtb112 = cadd(
                                                                    cmul(rvrb111, wrn012),
                                                                    cmul(rvrb112, wrn022),
                                                                );
                                                                rtb121 = cadd(
                                                                    cmul(rvrb121, wrn011),
                                                                    cmul(rvrb122, wrn021),
                                                                );
                                                                rtb122 = cadd(
                                                                    cmul(rvrb121, wrn012),
                                                                    cmul(rvrb122, wrn022),
                                                                );
                                                                rtb211 = cadd(
                                                                    cmul(ewrd11, rtb111),
                                                                    cmul(ewrd12, rtb121),
                                                                );
                                                                rtb212 = cadd(
                                                                    cmul(ewrd11, rtb112),
                                                                    cmul(ewrd12, rtb122),
                                                                );
                                                                rtb221 = cadd(
                                                                    cmul(ewrd21, rtb111),
                                                                    cmul(ewrd22, rtb121),
                                                                );
                                                                rtb222 = cadd(
                                                                    cmul(ewrd21, rtb112),
                                                                    cmul(ewrd22, rtb122),
                                                                );
                                                                rvrb211 = cadd(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    rtb211,
                                                                );
                                                                rvrb212 = rtb212;
                                                                rvrb221 = rtb221;
                                                                rvrb222 = cadd(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    rtb222,
                                                                );
                                                                rtb111 = cadd(
                                                                    ewigh11,
                                                                    (*checked_at(
                                                                        td0n11,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                rtb112 = cadd(
                                                                    ewigh11,
                                                                    (*checked_at(
                                                                        td0n12,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                rtb121 = cadd(
                                                                    ewigh22,
                                                                    (*checked_at(
                                                                        td0n21,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                rtb122 = cadd(
                                                                    ewigh22,
                                                                    (*checked_at(
                                                                        td0n22,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                                ewd0211 = cadd(
                                                                    cmul(rvrb111, rtb111),
                                                                    cmul(rvrb112, rtb121),
                                                                );
                                                                ewd0212 = cadd(
                                                                    cmul(rvrb111, rtb112),
                                                                    cmul(rvrb112, rtb122),
                                                                );
                                                                ewd0221 = cadd(
                                                                    cmul(rvrb121, rtb111),
                                                                    cmul(rvrb122, rtb121),
                                                                );
                                                                ewd0222 = cadd(
                                                                    cmul(rvrb121, rtb112),
                                                                    cmul(rvrb122, rtb122),
                                                                );
                                                                (*checked_at(
                                                                    td0n11,
                                                                    (ip) as isize,
                                                                )) = cadd(
                                                                    cmul(ewd0211, td11),
                                                                    cmul(ewd0221, td12),
                                                                );
                                                                (*checked_at(
                                                                    td0n12,
                                                                    (ip) as isize,
                                                                )) = cadd(
                                                                    cmul(ewd0212, td11),
                                                                    cmul(ewd0222, td12),
                                                                );
                                                                (*checked_at(
                                                                    td0n21,
                                                                    (ip) as isize,
                                                                )) = cadd(
                                                                    cmul(ewd0211, td21),
                                                                    cmul(ewd0221, td22),
                                                                );
                                                                (*checked_at(
                                                                    td0n22,
                                                                    (ip) as isize,
                                                                )) = cadd(
                                                                    cmul(ewd0212, td21),
                                                                    cmul(ewd0222, td22),
                                                                );
                                                                rtb111 = cmul(ewigh11, tu11);
                                                                rtb112 = cmul(ewigh11, tu12);
                                                                rtb121 = cmul(ewigh22, tu21);
                                                                rtb122 = cmul(ewigh22, tu22);
                                                                ewtu211 = cadd(
                                                                    cmul(rvrb211, rtb111),
                                                                    cmul(rvrb212, rtb121),
                                                                );
                                                                ewtu212 = cadd(
                                                                    cmul(rvrb211, rtb112),
                                                                    cmul(rvrb212, rtb122),
                                                                );
                                                                ewtu221 = cadd(
                                                                    cmul(rvrb221, rtb111),
                                                                    cmul(rvrb222, rtb121),
                                                                );
                                                                ewtu222 = cadd(
                                                                    cmul(rvrb221, rtb112),
                                                                    cmul(rvrb222, rtb122),
                                                                );
                                                                tun0p11 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun011,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu211,
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun012,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu221,
                                                                    ),
                                                                );
                                                                tun0p12 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun011,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu212,
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun012,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu222,
                                                                    ),
                                                                );
                                                                tun0p21 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun021,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu211,
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun022,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu221,
                                                                    ),
                                                                );
                                                                tun0p22 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun021,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu212,
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tun022,
                                                                            (ip) as isize,
                                                                        )),
                                                                        ewtu222,
                                                                    ),
                                                                );
                                                                rtb111 = cadd(
                                                                    cmul(ewrd11, ewd0211),
                                                                    cmul(ewd0221, ewrd12),
                                                                );
                                                                rtb112 = cadd(
                                                                    cmul(ewrd11, ewd0212),
                                                                    cmul(ewd0222, ewrd12),
                                                                );
                                                                rtb121 = cadd(
                                                                    cmul(ewrd21, ewd0211),
                                                                    cmul(ewd0221, ewrd22),
                                                                );
                                                                rtb122 = cadd(
                                                                    cmul(ewrd21, ewd0212),
                                                                    cmul(ewd0222, ewrd22),
                                                                );
                                                                rtb211 = cadd(
                                                                    cmul(
                                                                        rtb111,
                                                                        (*checked_at(
                                                                            tun011,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb121,
                                                                        (*checked_at(
                                                                            tun012,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb212 = cadd(
                                                                    cmul(
                                                                        rtb112,
                                                                        (*checked_at(
                                                                            tun011,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb122,
                                                                        (*checked_at(
                                                                            tun012,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb221 = cadd(
                                                                    cmul(
                                                                        rtb111,
                                                                        (*checked_at(
                                                                            tun021,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb121,
                                                                        (*checked_at(
                                                                            tun022,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb222 = cadd(
                                                                    cmul(
                                                                        rtb112,
                                                                        (*checked_at(
                                                                            tun021,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb122,
                                                                        (*checked_at(
                                                                            tun022,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rd0np11 = cadd(
                                                                    (*checked_at(
                                                                        rd0n11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb211,
                                                                );
                                                                rd0np12 = cadd(
                                                                    (*checked_at(
                                                                        rd0n12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb212,
                                                                );
                                                                rd0np21 = cadd(
                                                                    (*checked_at(
                                                                        rd0n21,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb221,
                                                                );
                                                                rd0np22 = cadd(
                                                                    (*checked_at(
                                                                        rd0n22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb222,
                                                                );
                                                                rtb111 = cadd(
                                                                    cmul(wrn011, ewtu211),
                                                                    cmul(wrn012, ewtu221),
                                                                );
                                                                rtb112 = cadd(
                                                                    cmul(wrn011, ewtu212),
                                                                    cmul(wrn012, ewtu222),
                                                                );
                                                                rtb121 = cadd(
                                                                    cmul(wrn021, ewtu211),
                                                                    cmul(wrn022, ewtu221),
                                                                );
                                                                rtb122 = cadd(
                                                                    cmul(wrn021, ewtu212),
                                                                    cmul(wrn022, ewtu222),
                                                                );
                                                                rtb211 = cadd(
                                                                    cmul(td11, rtb111),
                                                                    cmul(td12, rtb121),
                                                                );
                                                                rtb212 = cadd(
                                                                    cmul(td11, rtb112),
                                                                    cmul(td12, rtb122),
                                                                );
                                                                rtb221 = cadd(
                                                                    cmul(td21, rtb111),
                                                                    cmul(td22, rtb121),
                                                                );
                                                                rtb222 = cadd(
                                                                    cmul(td21, rtb112),
                                                                    cmul(td22, rtb122),
                                                                );
                                                                (*checked_at(
                                                                    run011,
                                                                    (ip) as isize,
                                                                )) = cadd(ru11, rtb211);
                                                                (*checked_at(
                                                                    run012,
                                                                    (ip) as isize,
                                                                )) = cadd(ru12, rtb212);
                                                                (*checked_at(
                                                                    run021,
                                                                    (ip) as isize,
                                                                )) = cadd(ru21, rtb221);
                                                                (*checked_at(
                                                                    run022,
                                                                    (ip) as isize,
                                                                )) = cadd(ru22, rtb222);
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )) = rd0np11;
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )) = rd0np12;
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )) = rd0np21;
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )) = rd0np22;
                                                                (*checked_at(
                                                                    tun011,
                                                                    (ip) as isize,
                                                                )) = tun0p11;
                                                                (*checked_at(
                                                                    tun012,
                                                                    (ip) as isize,
                                                                )) = tun0p12;
                                                                (*checked_at(
                                                                    tun021,
                                                                    (ip) as isize,
                                                                )) = tun0p21;
                                                                (*checked_at(
                                                                    tun022,
                                                                    (ip) as isize,
                                                                )) = tun0p22;
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                } else {
                                                    if (((((at1).r as f64) != 0f64)
                                                        || (((at1).i as f64) != 0f64))
                                                        && ((((at2).r as f64) != 0f64)
                                                            || (((at2).i as f64) != 0f64)))
                                                    {
                                                        /* case 4: solid-solid */
                                                        if (((((al1).r == (al2).r)
                                                            && ((al1).i == (al2).i))
                                                            && ((at1).r == (at2).r))
                                                            && ((at1).i == (at2).i))
                                                        {
                                                            {
                                                                ip = 0i32;
                                                                while (ip < block_size) {
                                                                    {
                                                                        ijk1 = (ik1 + ip);
                                                                        gl1 = (*checked_at(
                                                                            gl,
                                                                            (ijk1) as isize,
                                                                        ));
                                                                        gt1 = (*checked_at(
                                                                            gt,
                                                                            (ijk1) as isize,
                                                                        ));
                                                                        ewigh11 = cwp_cexp(cmul(
                                                                            cmplx(
                                                                                (0f64 as f32),
                                                                                (1f64 as f32),
                                                                            ),
                                                                            cmul(
                                                                                t1,
                                                                                cmul(wpie, gl1),
                                                                            ),
                                                                        ));
                                                                        ewigh22 = cwp_cexp(cmul(
                                                                            cmplx(
                                                                                (0f64 as f32),
                                                                                (1f64 as f32),
                                                                            ),
                                                                            cmul(
                                                                                t1,
                                                                                cmul(wpie, gt1),
                                                                            ),
                                                                        ));
                                                                        tun0p11 = cmul(
                                                                            (*checked_at(
                                                                                tun011,
                                                                                (ip) as isize,
                                                                            )),
                                                                            ewigh11,
                                                                        );
                                                                        tun0p12 = cmul(
                                                                            (*checked_at(
                                                                                tun012,
                                                                                (ip) as isize,
                                                                            )),
                                                                            ewigh22,
                                                                        );
                                                                        tun0p21 = cmul(
                                                                            (*checked_at(
                                                                                tun021,
                                                                                (ip) as isize,
                                                                            )),
                                                                            ewigh11,
                                                                        );
                                                                        tun0p22 = cmul(
                                                                            (*checked_at(
                                                                                tun022,
                                                                                (ip) as isize,
                                                                            )),
                                                                            ewigh22,
                                                                        );
                                                                        rtb111 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                run011,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb112 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                run012,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb121 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                run021,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb122 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                run022,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        (*checked_at(
                                                                            run011,
                                                                            (ip) as isize,
                                                                        )) = cmul(ewigh11, rtb111);
                                                                        (*checked_at(
                                                                            run012,
                                                                            (ip) as isize,
                                                                        )) = cmul(ewigh22, rtb112);
                                                                        (*checked_at(
                                                                            run021,
                                                                            (ip) as isize,
                                                                        )) = cmul(ewigh11, rtb121);
                                                                        (*checked_at(
                                                                            run022,
                                                                            (ip) as isize,
                                                                        )) = cmul(ewigh22, rtb122);
                                                                        td0np11 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                td0n11,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        td0np12 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                td0n12,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        td0np21 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                td0n21,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        td0np22 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                td0n22,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        (*checked_at(
                                                                            tun011,
                                                                            (ip) as isize,
                                                                        )) = tun0p11;
                                                                        (*checked_at(
                                                                            tun012,
                                                                            (ip) as isize,
                                                                        )) = tun0p12;
                                                                        (*checked_at(
                                                                            tun021,
                                                                            (ip) as isize,
                                                                        )) = tun0p21;
                                                                        (*checked_at(
                                                                            tun022,
                                                                            (ip) as isize,
                                                                        )) = tun0p22;
                                                                        (*checked_at(
                                                                            td0n11,
                                                                            (ip) as isize,
                                                                        )) = td0np11;
                                                                        (*checked_at(
                                                                            td0n12,
                                                                            (ip) as isize,
                                                                        )) = td0np12;
                                                                        (*checked_at(
                                                                            td0n21,
                                                                            (ip) as isize,
                                                                        )) = td0np21;
                                                                        (*checked_at(
                                                                            td0n22,
                                                                            (ip) as isize,
                                                                        )) = td0np22;
                                                                    }
                                                                    ip += 1;
                                                                }
                                                            }
                                                        } else {
                                                            d = crmul(
                                                                csub(
                                                                    cdiv(rho2, cmul(at2, at2)),
                                                                    cdiv(rho1, cmul(at1, at1)),
                                                                ),
                                                                (2f64 as f32),
                                                            );
                                                            {
                                                                ip = 0i32;
                                                                while (ip < block_size) {
                                                                    {
                                                                        ijk1 = (ik1 + ip);
                                                                        ijk2 = (ik2 + ip);
                                                                        psq = (*checked_at(
                                                                            pp,
                                                                            (ip) as isize,
                                                                        ));
                                                                        gl1 = (*checked_at(
                                                                            gl,
                                                                            (ijk1) as isize,
                                                                        ));
                                                                        gt1 = (*checked_at(
                                                                            gt,
                                                                            (ijk1) as isize,
                                                                        ));
                                                                        gam1 = (*checked_at(
                                                                            gam,
                                                                            (ijk1) as isize,
                                                                        ));
                                                                        gl2 = (*checked_at(
                                                                            gl,
                                                                            (ijk2) as isize,
                                                                        ));
                                                                        gt2 = (*checked_at(
                                                                            gt,
                                                                            (ijk2) as isize,
                                                                        ));
                                                                        gam2 = (*checked_at(
                                                                            gam,
                                                                            (ijk2) as isize,
                                                                        ));
                                                                        ewigh11 = cwp_cexp(cmul(
                                                                            cmplx(
                                                                                (0f64 as f32),
                                                                                (1f64 as f32),
                                                                            ),
                                                                            cmul(
                                                                                wpie,
                                                                                cmul(gl1, t1),
                                                                            ),
                                                                        ));
                                                                        ewigh22 = cwp_cexp(cmul(
                                                                            cmplx(
                                                                                (0f64 as f32),
                                                                                (1f64 as f32),
                                                                            ),
                                                                            cmul(
                                                                                wpie,
                                                                                cmul(gt1, t1),
                                                                            ),
                                                                        ));
                                                                        a = csub(gam2, gam1);
                                                                        b = cadd(rho1, a);
                                                                        c = csub(rho2, a);
                                                                        e = cadd(
                                                                            cmul(gl1, b),
                                                                            cmul(gl2, c),
                                                                        );
                                                                        f = cadd(
                                                                            cmul(gt1, b),
                                                                            cmul(gt2, c),
                                                                        );
                                                                        g = csub(
                                                                            a,
                                                                            cmul(gl1, cmul(gt2, d)),
                                                                        );
                                                                        h = csub(
                                                                            a,
                                                                            cmul(gl2, cmul(gt1, d)),
                                                                        );
                                                                        z = cdiv(
                                                                            cmplx(
                                                                                (1f64 as f32),
                                                                                (0f64 as f32),
                                                                            ),
                                                                            cadd(
                                                                                cmul(e, f),
                                                                                cmul(
                                                                                    g,
                                                                                    cmul(psq, h),
                                                                                ),
                                                                            ),
                                                                        );
                                                                        y = cmul(
                                                                            z,
                                                                            crmul(
                                                                                rho1,
                                                                                (2f64 as f32),
                                                                            ),
                                                                        );
                                                                        td11 =
                                                                            cmul(y, cmul(gl1, f));
                                                                        td12 = cmul(y,cmul(gt1,cmul((*checked_at(p,(ip) as isize)),g)));
                                                                        td21 = cneg(cmul(y,cmul(gl1,cmul((*checked_at(p,(ip) as isize)),h))));
                                                                        td22 =
                                                                            cmul(y, cmul(gt1, e));
                                                                        y = crmul(cmul(z,cmul((*checked_at(p,(ip) as isize)),cadd(cmul(a,b),cmul(gl2,cmul(gt2,cmul(c,d)))))),(2f64 as f32));
                                                                        rd11 = cmul(
                                                                            z,
                                                                            csub(
                                                                                cmul(
                                                                                    f,
                                                                                    csub(
                                                                                        cmul(
                                                                                            gl1, b,
                                                                                        ),
                                                                                        cmul(
                                                                                            gl2, c,
                                                                                        ),
                                                                                    ),
                                                                                ),
                                                                                cmul(
                                                                                    cadd(
                                                                                        a,
                                                                                        cmul(
                                                                                            gl1,
                                                                                            cmul(
                                                                                                gt2,
                                                                                                d,
                                                                                            ),
                                                                                        ),
                                                                                    ),
                                                                                    cmul(psq, h),
                                                                                ),
                                                                            ),
                                                                        );
                                                                        rd12 = cmul(y, gt1);
                                                                        rd21 = cneg(cmul(y, gl1));
                                                                        rd22 = cmul(
                                                                            z,
                                                                            csub(
                                                                                cmul(
                                                                                    e,
                                                                                    csub(
                                                                                        cmul(
                                                                                            gt1, b,
                                                                                        ),
                                                                                        cmul(
                                                                                            gt2, c,
                                                                                        ),
                                                                                    ),
                                                                                ),
                                                                                cmul(
                                                                                    cadd(
                                                                                        a,
                                                                                        cmul(
                                                                                            gt1,
                                                                                            cmul(
                                                                                                gl2,
                                                                                                d,
                                                                                            ),
                                                                                        ),
                                                                                    ),
                                                                                    cmul(psq, g),
                                                                                ),
                                                                            ),
                                                                        );
                                                                        y = cmul(
                                                                            z,
                                                                            crmul(
                                                                                rho2,
                                                                                (2f64 as f32),
                                                                            ),
                                                                        );
                                                                        tu11 =
                                                                            cmul(y, cmul(gl2, f));
                                                                        tu12 = cmul(y,cmul(gt2,cmul((*checked_at(p,(ip) as isize)),h)));
                                                                        tu21 = cneg(cmul(y,cmul(gl2,cmul((*checked_at(p,(ip) as isize)),g))));
                                                                        tu22 =
                                                                            cmul(y, cmul(gt2, e));
                                                                        y = crmul(cmul(z,cmul((*checked_at(p,(ip) as isize)),cadd(cmul(a,c),cmul(gl1,cmul(gt1,cmul(b,d)))))),(2f64 as f32));
                                                                        ru11 = cmul(
                                                                            z,
                                                                            csub(
                                                                                cmul(
                                                                                    f,
                                                                                    csub(
                                                                                        cmul(
                                                                                            gl2, c,
                                                                                        ),
                                                                                        cmul(
                                                                                            gl1, b,
                                                                                        ),
                                                                                    ),
                                                                                ),
                                                                                cmul(
                                                                                    cadd(
                                                                                        a,
                                                                                        cmul(
                                                                                            gl2,
                                                                                            cmul(
                                                                                                gt1,
                                                                                                d,
                                                                                            ),
                                                                                        ),
                                                                                    ),
                                                                                    cmul(psq, g),
                                                                                ),
                                                                            ),
                                                                        );
                                                                        ru12 = cmul(y, gt2);
                                                                        ru21 = cneg(cmul(y, gl2));
                                                                        ru22 = cmul(
                                                                            z,
                                                                            csub(
                                                                                cmul(
                                                                                    e,
                                                                                    csub(
                                                                                        cmul(
                                                                                            gt2, c,
                                                                                        ),
                                                                                        cmul(
                                                                                            gt1, b,
                                                                                        ),
                                                                                    ),
                                                                                ),
                                                                                cmul(
                                                                                    cadd(
                                                                                        a,
                                                                                        cmul(
                                                                                            gl1,
                                                                                            cmul(
                                                                                                gt2,
                                                                                                d,
                                                                                            ),
                                                                                        ),
                                                                                    ),
                                                                                    cmul(psq, h),
                                                                                ),
                                                                            ),
                                                                        );
                                                                        ewrd11 =
                                                                            cmul(ewigh11, rd11);
                                                                        ewrd12 =
                                                                            cmul(ewigh11, rd12);
                                                                        ewrd21 =
                                                                            cmul(ewigh22, rd21);
                                                                        ewrd22 =
                                                                            cmul(ewigh22, rd22);
                                                                        wrn011 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                run011,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        wrn012 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                run012,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        wrn021 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                run021,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        wrn022 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                run022,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb111 = cadd(
                                                                            cmul(wrn011, ewrd11),
                                                                            cmul(wrn012, ewrd21),
                                                                        );
                                                                        rtb112 = cadd(
                                                                            cmul(wrn011, ewrd12),
                                                                            cmul(wrn012, ewrd22),
                                                                        );
                                                                        rtb121 = cadd(
                                                                            cmul(wrn021, ewrd11),
                                                                            cmul(wrn022, ewrd21),
                                                                        );
                                                                        rtb122 = cadd(
                                                                            cmul(wrn021, ewrd12),
                                                                            cmul(wrn022, ewrd22),
                                                                        );
                                                                        x = csub(
                                                                            cmplx(
                                                                                (1f64 as f32),
                                                                                (0f64 as f32),
                                                                            ),
                                                                            rtb111,
                                                                        );
                                                                        y = csub(
                                                                            cmplx(
                                                                                (1f64 as f32),
                                                                                (0f64 as f32),
                                                                            ),
                                                                            rtb122,
                                                                        );
                                                                        det = csub(
                                                                            cmul(x, y),
                                                                            cmul(rtb112, rtb121),
                                                                        );
                                                                        rvrb111 = cdiv(y, det);
                                                                        rvrb112 = cdiv(rtb112, det);
                                                                        rvrb121 = cdiv(rtb121, det);
                                                                        rvrb122 = cdiv(x, det);
                                                                        rtb111 = cadd(
                                                                            cmul(wrn011, rvrb111),
                                                                            cmul(rvrb112, wrn021),
                                                                        );
                                                                        rtb112 = cadd(
                                                                            cmul(wrn012, rvrb111),
                                                                            cmul(rvrb112, wrn022),
                                                                        );
                                                                        rtb121 = cadd(
                                                                            cmul(wrn011, rvrb121),
                                                                            cmul(rvrb122, wrn021),
                                                                        );
                                                                        rtb122 = cadd(
                                                                            cmul(wrn012, rvrb121),
                                                                            cmul(rvrb122, wrn022),
                                                                        );
                                                                        rtb211 = cadd(
                                                                            cmul(ewrd11, rtb111),
                                                                            cmul(ewrd12, rtb121),
                                                                        );
                                                                        rtb212 = cadd(
                                                                            cmul(ewrd11, rtb112),
                                                                            cmul(ewrd12, rtb122),
                                                                        );
                                                                        rtb221 = cadd(
                                                                            cmul(ewrd21, rtb111),
                                                                            cmul(ewrd22, rtb121),
                                                                        );
                                                                        rtb222 = cadd(
                                                                            cmul(ewrd21, rtb112),
                                                                            cmul(ewrd22, rtb122),
                                                                        );
                                                                        rvrb211 = cadd(
                                                                            cmplx(
                                                                                (1f64 as f32),
                                                                                (0f64 as f32),
                                                                            ),
                                                                            rtb211,
                                                                        );
                                                                        rvrb212 = rtb212;
                                                                        rvrb221 = rtb221;
                                                                        rvrb222 = cadd(
                                                                            cmplx(
                                                                                (1f64 as f32),
                                                                                (0f64 as f32),
                                                                            ),
                                                                            rtb222,
                                                                        );
                                                                        rtb111 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                td0n11,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb112 = cmul(
                                                                            ewigh11,
                                                                            (*checked_at(
                                                                                td0n12,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb121 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                td0n21,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        rtb122 = cmul(
                                                                            ewigh22,
                                                                            (*checked_at(
                                                                                td0n22,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        ewd0211 = cadd(
                                                                            cmul(rvrb111, rtb111),
                                                                            cmul(rvrb112, rtb121),
                                                                        );
                                                                        ewd0212 = cadd(
                                                                            cmul(rvrb111, rtb112),
                                                                            cmul(rvrb112, rtb122),
                                                                        );
                                                                        ewd0221 = cadd(
                                                                            cmul(rvrb121, rtb111),
                                                                            cmul(rvrb122, rtb121),
                                                                        );
                                                                        ewd0222 = cadd(
                                                                            cmul(rvrb121, rtb112),
                                                                            cmul(rvrb122, rtb122),
                                                                        );
                                                                        (*checked_at(
                                                                            td0n11,
                                                                            (ip) as isize,
                                                                        )) = cadd(
                                                                            cmul(td11, ewd0211),
                                                                            cmul(ewd0221, td12),
                                                                        );
                                                                        (*checked_at(
                                                                            td0n12,
                                                                            (ip) as isize,
                                                                        )) = cadd(
                                                                            cmul(td11, ewd0212),
                                                                            cmul(ewd0222, td12),
                                                                        );
                                                                        (*checked_at(
                                                                            td0n21,
                                                                            (ip) as isize,
                                                                        )) = cadd(
                                                                            cmul(td21, ewd0211),
                                                                            cmul(ewd0221, td22),
                                                                        );
                                                                        (*checked_at(
                                                                            td0n22,
                                                                            (ip) as isize,
                                                                        )) = cadd(
                                                                            cmul(td21, ewd0212),
                                                                            cmul(ewd0222, td22),
                                                                        );
                                                                        rtb111 =
                                                                            cmul(ewigh11, tu11);
                                                                        rtb112 =
                                                                            cmul(ewigh11, tu12);
                                                                        rtb121 =
                                                                            cmul(ewigh22, tu21);
                                                                        rtb122 =
                                                                            cmul(ewigh22, tu22);
                                                                        ewtu211 = cadd(
                                                                            cmul(rvrb211, rtb111),
                                                                            cmul(rvrb212, rtb121),
                                                                        );
                                                                        ewtu212 = cadd(
                                                                            cmul(rvrb211, rtb112),
                                                                            cmul(rvrb212, rtb122),
                                                                        );
                                                                        ewtu221 = cadd(
                                                                            cmul(rvrb221, rtb111),
                                                                            cmul(rvrb222, rtb121),
                                                                        );
                                                                        ewtu222 = cadd(
                                                                            cmul(rvrb221, rtb112),
                                                                            cmul(rvrb222, rtb122),
                                                                        );
                                                                        tun0p11 = cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun011,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu211,
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun012,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu221,
                                                                            ),
                                                                        );
                                                                        tun0p12 = cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun011,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu212,
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun012,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu222,
                                                                            ),
                                                                        );
                                                                        tun0p21 = cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun021,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu211,
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun022,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu221,
                                                                            ),
                                                                        );
                                                                        tun0p22 = cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun021,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu212,
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    tun022,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                ewtu222,
                                                                            ),
                                                                        );
                                                                        rtb111 = cadd(
                                                                            cmul(ewrd11, ewd0211),
                                                                            cmul(ewrd12, ewd0221),
                                                                        );
                                                                        rtb112 = cadd(
                                                                            cmul(ewrd11, ewd0212),
                                                                            cmul(ewrd12, ewd0222),
                                                                        );
                                                                        rtb121 = cadd(
                                                                            cmul(ewrd21, ewd0211),
                                                                            cmul(ewrd22, ewd0221),
                                                                        );
                                                                        rtb122 = cadd(
                                                                            cmul(ewrd21, ewd0212),
                                                                            cmul(ewrd22, ewd0222),
                                                                        );
                                                                        rtb211 = cadd(
                                                                            cmul(
                                                                                rtb111,
                                                                                (*checked_at(
                                                                                    tun011,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb121,
                                                                                (*checked_at(
                                                                                    tun012,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb212 = cadd(
                                                                            cmul(
                                                                                rtb112,
                                                                                (*checked_at(
                                                                                    tun011,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb122,
                                                                                (*checked_at(
                                                                                    tun012,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb221 = cadd(
                                                                            cmul(
                                                                                rtb111,
                                                                                (*checked_at(
                                                                                    tun021,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb121,
                                                                                (*checked_at(
                                                                                    tun022,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb222 = cadd(
                                                                            cmul(
                                                                                rtb112,
                                                                                (*checked_at(
                                                                                    tun021,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb122,
                                                                                (*checked_at(
                                                                                    tun022,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rd0np11 = cadd(
                                                                            (*checked_at(
                                                                                rd0n11,
                                                                                (ip) as isize,
                                                                            )),
                                                                            rtb211,
                                                                        );
                                                                        rd0np12 = cadd(
                                                                            (*checked_at(
                                                                                rd0n12,
                                                                                (ip) as isize,
                                                                            )),
                                                                            rtb212,
                                                                        );
                                                                        rd0np21 = cadd(
                                                                            (*checked_at(
                                                                                rd0n21,
                                                                                (ip) as isize,
                                                                            )),
                                                                            rtb221,
                                                                        );
                                                                        rd0np22 = cadd(
                                                                            (*checked_at(
                                                                                rd0n22,
                                                                                (ip) as isize,
                                                                            )),
                                                                            rtb222,
                                                                        );
                                                                        rtb111 = cadd(
                                                                            cmul(wrn011, ewtu211),
                                                                            cmul(ewtu221, wrn012),
                                                                        );
                                                                        rtb112 = cadd(
                                                                            cmul(wrn011, ewtu212),
                                                                            cmul(ewtu222, wrn012),
                                                                        );
                                                                        rtb121 = cadd(
                                                                            cmul(wrn021, ewtu211),
                                                                            cmul(ewtu221, wrn022),
                                                                        );
                                                                        rtb122 = cadd(
                                                                            cmul(wrn021, ewtu212),
                                                                            cmul(ewtu222, wrn022),
                                                                        );
                                                                        rtb211 = cadd(
                                                                            cmul(td11, rtb111),
                                                                            cmul(td12, rtb121),
                                                                        );
                                                                        rtb212 = cadd(
                                                                            cmul(td11, rtb112),
                                                                            cmul(td12, rtb122),
                                                                        );
                                                                        rtb221 = cadd(
                                                                            cmul(td21, rtb111),
                                                                            cmul(td22, rtb121),
                                                                        );
                                                                        rtb222 = cadd(
                                                                            cmul(td21, rtb112),
                                                                            cmul(td22, rtb122),
                                                                        );
                                                                        (*checked_at(
                                                                            run011,
                                                                            (ip) as isize,
                                                                        )) = cadd(ru11, rtb211);
                                                                        (*checked_at(
                                                                            run012,
                                                                            (ip) as isize,
                                                                        )) = cadd(ru12, rtb212);
                                                                        (*checked_at(
                                                                            run021,
                                                                            (ip) as isize,
                                                                        )) = cadd(ru21, rtb221);
                                                                        (*checked_at(
                                                                            run022,
                                                                            (ip) as isize,
                                                                        )) = cadd(ru22, rtb222);
                                                                        (*checked_at(
                                                                            rd0n11,
                                                                            (ip) as isize,
                                                                        )) = rd0np11;
                                                                        (*checked_at(
                                                                            rd0n12,
                                                                            (ip) as isize,
                                                                        )) = rd0np12;
                                                                        (*checked_at(
                                                                            rd0n21,
                                                                            (ip) as isize,
                                                                        )) = rd0np21;
                                                                        (*checked_at(
                                                                            rd0n22,
                                                                            (ip) as isize,
                                                                        )) = rd0np22;
                                                                        (*checked_at(
                                                                            tun011,
                                                                            (ip) as isize,
                                                                        )) = tun0p11;
                                                                        (*checked_at(
                                                                            tun012,
                                                                            (ip) as isize,
                                                                        )) = tun0p12;
                                                                        (*checked_at(
                                                                            tun021,
                                                                            (ip) as isize,
                                                                        )) = tun0p21;
                                                                        (*checked_at(
                                                                            tun022,
                                                                            (ip) as isize,
                                                                        )) = tun0p22;
                                                                    }
                                                                    ip += 1;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    /* closes "if" of layers after first one*/
                                    if (il == ((*checked_at(lobs, (ijk) as isize)) - 2i32)) {
                                        ik1 = (ijk * block_size);
                                        if ((((*checked_at(acoustic, (ijk) as isize)) as f64)
                                            == 2f64)
                                            || (((*checked_at(acoustic, (ijk) as isize)) as f64)
                                                == 4f64))
                                        {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        iz = (ik1 + ip);
                                                        (*checked_at(rurf11, (iz) as isize)) =
                                                            (*checked_at(run011, (ip) as isize));
                                                        (*checked_at(rurf12, (iz) as isize)) =
                                                            (*checked_at(run012, (ip) as isize));
                                                        (*checked_at(rurf21, (iz) as isize)) =
                                                            (*checked_at(run021, (ip) as isize));
                                                        (*checked_at(rurf22, (iz) as isize)) =
                                                            (*checked_at(run022, (ip) as isize));
                                                        (*checked_at(rdfr11, (iz) as isize)) =
                                                            (*checked_at(rd0n11, (ip) as isize));
                                                        (*checked_at(rdfr12, (iz) as isize)) =
                                                            (*checked_at(rd0n12, (ip) as isize));
                                                        (*checked_at(rdfr21, (iz) as isize)) =
                                                            (*checked_at(rd0n21, (ip) as isize));
                                                        (*checked_at(rdfr22, (iz) as isize)) =
                                                            (*checked_at(rd0n22, (ip) as isize));
                                                        (*checked_at(turf11, (iz) as isize)) =
                                                            (*checked_at(tun011, (ip) as isize));
                                                        (*checked_at(turf12, (iz) as isize)) =
                                                            (*checked_at(tun012, (ip) as isize));
                                                        (*checked_at(turf21, (iz) as isize)) =
                                                            (*checked_at(tun021, (ip) as isize));
                                                        (*checked_at(turf22, (iz) as isize)) =
                                                            (*checked_at(tun022, (ip) as isize));
                                                        (*checked_at(tdfr11, (iz) as isize)) =
                                                            (*checked_at(td0n11, (ip) as isize));
                                                        (*checked_at(tdfr12, (iz) as isize)) =
                                                            (*checked_at(td0n12, (ip) as isize));
                                                        (*checked_at(tdfr21, (iz) as isize)) =
                                                            (*checked_at(td0n21, (ip) as isize));
                                                        (*checked_at(tdfr22, (iz) as isize)) =
                                                            (*checked_at(td0n22, (ip) as isize));
                                                        (*checked_at(tun011, (ip) as isize)) =
                                                            cmplx((1f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun022, (ip) as isize)) =
                                                            cmplx((1f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n11, (ip) as isize)) =
                                                            cmplx((1f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n22, (ip) as isize)) =
                                                            cmplx((1f64 as f32), (0f64 as f32));
                                                        (*checked_at(run011, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n11, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        } else {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        iz = (ik1 + ip);
                                                        (*checked_at(rurf11, (iz) as isize)) =
                                                            (*checked_at(run011, (ip) as isize));
                                                        (*checked_at(rurf12, (iz) as isize)) =
                                                            (*checked_at(run012, (ip) as isize));
                                                        (*checked_at(rurf21, (iz) as isize)) =
                                                            (*checked_at(run021, (ip) as isize));
                                                        (*checked_at(rurf22, (iz) as isize)) =
                                                            (*checked_at(run022, (ip) as isize));
                                                        (*checked_at(rdfr11, (ip) as isize)) =
                                                            (*checked_at(rd0n11, (ip) as isize));
                                                        (*checked_at(rdfr12, (ip) as isize)) =
                                                            (*checked_at(rd0n12, (ip) as isize));
                                                        (*checked_at(rdfr21, (ip) as isize)) =
                                                            (*checked_at(rd0n21, (ip) as isize));
                                                        (*checked_at(rdfr22, (ip) as isize)) =
                                                            (*checked_at(rd0n22, (ip) as isize));
                                                        (*checked_at(turf11, (iz) as isize)) =
                                                            (*checked_at(tun011, (ip) as isize));
                                                        (*checked_at(turf12, (iz) as isize)) =
                                                            (*checked_at(tun012, (ip) as isize));
                                                        (*checked_at(turf21, (iz) as isize)) =
                                                            (*checked_at(tun021, (ip) as isize));
                                                        (*checked_at(turf22, (iz) as isize)) =
                                                            (*checked_at(tun022, (ip) as isize));
                                                        (*checked_at(tdfr11, (iz) as isize)) =
                                                            (*checked_at(td0n11, (ip) as isize));
                                                        (*checked_at(tdfr12, (iz) as isize)) =
                                                            (*checked_at(td0n12, (ip) as isize));
                                                        (*checked_at(tdfr21, (iz) as isize)) =
                                                            (*checked_at(td0n21, (ip) as isize));
                                                        (*checked_at(tdfr22, (iz) as isize)) =
                                                            (*checked_at(td0n22, (ip) as isize));
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        }

                                        ijk += 1;
                                    }

                                    if (il == (lsource - 2i32)) {
                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    (*checked_at(rusf11, (ip) as isize)) =
                                                        (*checked_at(run011, (ip) as isize));
                                                    (*checked_at(rusf12, (ip) as isize)) =
                                                        (*checked_at(run012, (ip) as isize));
                                                    (*checked_at(rusf21, (ip) as isize)) =
                                                        (*checked_at(run021, (ip) as isize));
                                                    (*checked_at(rusf22, (ip) as isize)) =
                                                        (*checked_at(run022, (ip) as isize));
                                                    (*checked_at(rdfs11, (ip) as isize)) =
                                                        (*checked_at(rd0n11, (ip) as isize));
                                                    (*checked_at(rdfs12, (ip) as isize)) =
                                                        (*checked_at(rd0n12, (ip) as isize));
                                                    (*checked_at(rdfs21, (ip) as isize)) =
                                                        (*checked_at(rd0n21, (ip) as isize));
                                                    (*checked_at(rdfs22, (ip) as isize)) =
                                                        (*checked_at(rd0n22, (ip) as isize));
                                                    (*checked_at(run011, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run012, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run021, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run022, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n11, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n12, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n22, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tusf11, (ip) as isize)) =
                                                        (*checked_at(tun011, (ip) as isize));
                                                    (*checked_at(tusf12, (ip) as isize)) =
                                                        (*checked_at(tun012, (ip) as isize));
                                                    (*checked_at(tusf21, (ip) as isize)) =
                                                        (*checked_at(tun021, (ip) as isize));
                                                    (*checked_at(tusf22, (ip) as isize)) =
                                                        (*checked_at(tun022, (ip) as isize));
                                                    (*checked_at(tdfs11, (ip) as isize)) =
                                                        (*checked_at(td0n11, (ip) as isize));
                                                    (*checked_at(tdfs12, (ip) as isize)) =
                                                        (*checked_at(td0n12, (ip) as isize));
                                                    (*checked_at(tdfs21, (ip) as isize)) =
                                                        (*checked_at(td0n21, (ip) as isize));
                                                    (*checked_at(tdfs22, (ip) as isize)) =
                                                        (*checked_at(td0n22, (ip) as isize));
                                                    (*checked_at(tun011, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun012, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun021, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun022, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n11, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n12, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n22, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                }
                                                ip += 1;
                                            }
                                        }
                                    }
                                }
                                il += 1;
                            }
                        }
                        /* loop over receivers */
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    ik1 = (iz * block_size);
                                    {
                                        ip = 0i32;
                                        while (ip < block_size) {
                                            {
                                                ijk1 = (ik1 + ip);
                                                if ((*checked_at(acoustic, (iz) as isize)) == 1i32)
                                                {
                                                    det = csub(
                                                        (*checked_at(rd0n11, (ip) as isize)),
                                                        (*checked_at(rdfr11, (ijk1) as isize)),
                                                    );
                                                } else {
                                                    det = csub(
                                                        cmul(
                                                            csub(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rdfr11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            csub(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rdfr22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        ),
                                                        cmul(
                                                            csub(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rdfr12,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            csub(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rdfr21,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        ),
                                                    );
                                                }

                                                if ((rcabs(det) as f64)
                                                    <= 9.9999999999999998E-13f64)
                                                {
                                                    (*checked_at(prv, (iz) as isize)) =
                                                        (*checked_at(flag, (iz) as isize));
                                                    (*checked_at(flag, (iz) as isize)) = 2i32;
                                                    (*checked_at(ibl, (iz) as isize)) = ip;
                                                    break;
                                                }
                                            }
                                            ip += 1;
                                        }
                                    }
                                    if ((*checked_at(acoustic, (iz) as isize)) == 2i32) {
                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    ijk1 = (ik1 + ip);
                                                    (*checked_at(rnr11, (ip) as isize)) =
                                                        (*checked_at(rd0n11, (ip) as isize));
                                                    (*checked_at(rnr12, (ip) as isize)) =
                                                        (*checked_at(rd0n12, (ip) as isize));
                                                    (*checked_at(rnr21, (ip) as isize)) =
                                                        (*checked_at(rd0n21, (ip) as isize));
                                                    (*checked_at(rnr22, (ip) as isize)) =
                                                        (*checked_at(rd0n22, (ip) as isize));
                                                    (*checked_at(tdsr11, (ip) as isize)) =
                                                        (*checked_at(tdfr11, (ijk1) as isize));
                                                    (*checked_at(tdsr12, (ip) as isize)) =
                                                        (*checked_at(tdfr12, (ijk1) as isize));
                                                    (*checked_at(tdsr21, (ip) as isize)) =
                                                        (*checked_at(tdfr21, (ijk1) as isize));
                                                    (*checked_at(tdsr22, (ip) as isize)) =
                                                        (*checked_at(tdfr22, (ijk1) as isize));
                                                    (*checked_at(tusr11, (ip) as isize)) =
                                                        (*checked_at(turf11, (ijk1) as isize));
                                                    (*checked_at(tusr12, (ip) as isize)) =
                                                        (*checked_at(turf12, (ijk1) as isize));
                                                    (*checked_at(tusr21, (ip) as isize)) =
                                                        (*checked_at(turf21, (ijk1) as isize));
                                                    (*checked_at(tusr22, (ip) as isize)) =
                                                        (*checked_at(turf22, (ijk1) as isize));
                                                    (*checked_at(rdsr11, (ip) as isize)) =
                                                        (*checked_at(rdfr11, (ijk1) as isize));
                                                    (*checked_at(rdsr12, (ip) as isize)) =
                                                        (*checked_at(rdfr12, (ijk1) as isize));
                                                    (*checked_at(rdsr21, (ip) as isize)) =
                                                        (*checked_at(rdfr21, (ijk1) as isize));
                                                    (*checked_at(rdsr22, (ip) as isize)) =
                                                        (*checked_at(rdfr22, (ijk1) as isize));
                                                    (*checked_at(rusr11, (ip) as isize)) =
                                                        (*checked_at(rurf11, (ijk1) as isize));
                                                    (*checked_at(rusr12, (ip) as isize)) =
                                                        (*checked_at(rurf12, (ijk1) as isize));
                                                    (*checked_at(rusr21, (ip) as isize)) =
                                                        (*checked_at(rurf21, (ijk1) as isize));
                                                    (*checked_at(rusr22, (ip) as isize)) =
                                                        (*checked_at(rurf22, (ijk1) as isize));
                                                    rtb111 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            (*checked_at(rnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            (*checked_at(rnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb112 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            (*checked_at(rnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            (*checked_at(rnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb121 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf21, (ijk1) as isize)),
                                                            (*checked_at(rnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf22, (ijk1) as isize)),
                                                            (*checked_at(rnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb122 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf21, (ijk1) as isize)),
                                                            (*checked_at(rnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf22, (ijk1) as isize)),
                                                            (*checked_at(rnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    x = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb111,
                                                    );
                                                    y = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb122,
                                                    );
                                                    det = csub(cmul(x, y), cmul(rtb112, rtb121));
                                                    rvrb111 = cdiv(y, det);
                                                    rvrb112 = cdiv(rtb112, det);
                                                    rvrb121 = cdiv(rtb121, det);
                                                    rvrb122 = cdiv(x, det);
                                                    rtb111 = cadd(
                                                        cmul(
                                                            rvrb111,
                                                            (*checked_at(tdfr11, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb112,
                                                            (*checked_at(tdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb112 = cadd(
                                                        cmul(
                                                            rvrb111,
                                                            (*checked_at(tdfr12, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb112,
                                                            (*checked_at(tdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb121 = cadd(
                                                        cmul(
                                                            rvrb121,
                                                            (*checked_at(tdfr11, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb122,
                                                            (*checked_at(tdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb122 = cadd(
                                                        cmul(
                                                            rvrb121,
                                                            (*checked_at(tdfr12, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb122,
                                                            (*checked_at(tdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb211 = cadd(
                                                        cmul(
                                                            (*checked_at(rnr11, (ip) as isize)),
                                                            rtb111,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rnr12, (ip) as isize)),
                                                            rtb121,
                                                        ),
                                                    );
                                                    rtb212 = cadd(
                                                        cmul(
                                                            (*checked_at(rnr11, (ip) as isize)),
                                                            rtb112,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rnr12, (ip) as isize)),
                                                            rtb122,
                                                        ),
                                                    );
                                                    rtb221 = cadd(
                                                        cmul(
                                                            (*checked_at(rnr21, (ip) as isize)),
                                                            rtb111,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rnr22, (ip) as isize)),
                                                            rtb121,
                                                        ),
                                                    );
                                                    rtb222 = cadd(
                                                        cmul(
                                                            (*checked_at(rnr21, (ip) as isize)),
                                                            rtb112,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rnr22, (ip) as isize)),
                                                            rtb122,
                                                        ),
                                                    );
                                                    (*checked_at(rd0n11, (ip) as isize)) = cadd(
                                                        cmul(
                                                            rtb211,
                                                            (*checked_at(turf11, (ijk1) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                rtb221,
                                                                (*checked_at(
                                                                    turf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(rdfr11, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    (*checked_at(rd0n12, (ip) as isize)) = cadd(
                                                        cmul(
                                                            rtb212,
                                                            (*checked_at(turf11, (ijk1) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                rtb222,
                                                                (*checked_at(
                                                                    turf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(rdfr12, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    (*checked_at(rd0n21, (ip) as isize)) = cadd(
                                                        cmul(
                                                            rtb211,
                                                            (*checked_at(turf21, (ijk1) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                rtb221,
                                                                (*checked_at(
                                                                    turf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(rdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    (*checked_at(rd0n22, (ip) as isize)) = cadd(
                                                        cmul(
                                                            rtb212,
                                                            (*checked_at(turf21, (ijk1) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                rtb222,
                                                                (*checked_at(
                                                                    turf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(rdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                }
                                                ip += 1;
                                            }
                                        }
                                    }
                                }
                                iz += 1;
                            }
                        }
                        /* loop over receivers */
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    ik1 = (iz * block_size);
                                    lrec = (*checked_at(lobs, (iz) as isize));
                                    ik2 = ((lrec - 1i32) * block_size);
                                    /* receiver above the source */
                                    if (lsource > lrec) {
                                        if ((*checked_at(acoustic, (iz) as isize)) == 4i32) {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        ijk2 = (ik2 + ip);
                                                        rtb111 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb112 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb121 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb122 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rdfs22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        x = csub(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            rtb111,
                                                        );
                                                        y = csub(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            rtb122,
                                                        );
                                                        det =
                                                            csub(cmul(x, y), cmul(rtb112, rtb121));
                                                        rvrb111 = cdiv(y, det);
                                                        rvrb112 = cdiv(rtb112, det);
                                                        rvrb121 = cdiv(rtb121, det);
                                                        rvrb122 = cdiv(x, det);
                                                        rtb111 = cadd(
                                                            cmul(
                                                                rvrb111,
                                                                (*checked_at(
                                                                    tusf11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                rvrb112,
                                                                (*checked_at(
                                                                    tusf21,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb112 = cadd(
                                                            cmul(
                                                                rvrb111,
                                                                (*checked_at(
                                                                    tusf12,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                rvrb112,
                                                                (*checked_at(
                                                                    tusf22,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb121 = cadd(
                                                            cmul(
                                                                rvrb121,
                                                                (*checked_at(
                                                                    tusf11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                rvrb122,
                                                                (*checked_at(
                                                                    tusf21,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb122 = cadd(
                                                            cmul(
                                                                rvrb121,
                                                                (*checked_at(
                                                                    tusf12,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                rvrb122,
                                                                (*checked_at(
                                                                    tusf22,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb211 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb111,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb121,
                                                            ),
                                                        );
                                                        rtb212 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb112,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb122,
                                                            ),
                                                        );
                                                        rtb221 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb111,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb121,
                                                            ),
                                                        );
                                                        rtb222 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb112,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb122,
                                                            ),
                                                        );
                                                        rtb311 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs11,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb211,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs12,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb221,
                                                            ),
                                                        );
                                                        rtb312 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs11,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb212,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs12,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb222,
                                                            ),
                                                        );
                                                        rtb321 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs21,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb211,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs22,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb221,
                                                            ),
                                                        );
                                                        rtb322 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs21,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb212,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    tdfs22,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb222,
                                                            ),
                                                        );
                                                        (*checked_at(rusf11, (ip) as isize)) = cadd(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            rtb311,
                                                        );
                                                        (*checked_at(rusf12, (ip) as isize)) = cadd(
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                            rtb312,
                                                        );
                                                        (*checked_at(rusf21, (ip) as isize)) = cadd(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            rtb321,
                                                        );
                                                        (*checked_at(rusf22, (ip) as isize)) = cadd(
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                            rtb322,
                                                        );
                                                        rtb311 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf21,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb312 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf12,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf22,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb321 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf21,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb322 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf12,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf22,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        x = csub(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            rtb311,
                                                        );
                                                        y = csub(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            rtb322,
                                                        );
                                                        det =
                                                            csub(cmul(x, y), cmul(rtb312, rtb321));
                                                        rvrb111 = cdiv(y, det);
                                                        rvrb112 = cdiv(rtb312, det);
                                                        rvrb121 = cdiv(rtb321, det);
                                                        rvrb122 = cdiv(x, det);
                                                        fact11 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    sigmad1,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cadd(
                                                                (*checked_at(
                                                                    sigmau1,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rd0n12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        sigmad2,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            ),
                                                        );
                                                        fact12 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmad1,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cadd(
                                                                (*checked_at(
                                                                    sigmau2,
                                                                    (ip) as isize,
                                                                )),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rd0n22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        sigmad2,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            ),
                                                        );
                                                        vuzs1 = cadd(
                                                            cmul(rvrb111, fact11),
                                                            cmul(rvrb112, fact12),
                                                        );
                                                        vuzs2 = cadd(
                                                            cmul(rvrb121, fact11),
                                                            cmul(rvrb122, fact12),
                                                        );
                                                        (*checked_at(up, (ip) as isize)) = cmul(
                                                            (*checked_at(prs, (iz) as isize)),
                                                            cadd(
                                                                cmul(vuzs1, cadd(rtb111, rtb211)),
                                                                cmul(vuzs2, cadd(rtb112, rtb212)),
                                                            ),
                                                        );
                                                        rtb311 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                            csub(
                                                                cmul(
                                                                    (*checked_at(p, (ip) as isize)),
                                                                    (*checked_at(
                                                                        gt,
                                                                        (ijk2) as isize,
                                                                    )),
                                                                ),
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb312 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                            csub(
                                                                cmul(
                                                                    (*checked_at(p, (ip) as isize)),
                                                                    (*checked_at(
                                                                        gt,
                                                                        (ijk2) as isize,
                                                                    )),
                                                                ),
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb321 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(gl, (ijk2) as isize)),
                                                            ),
                                                            csub(
                                                                cmul(
                                                                    (*checked_at(p, (ip) as isize)),
                                                                    (*checked_at(
                                                                        rurf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                (*checked_at(gl, (ijk2) as isize)),
                                                            ),
                                                        );
                                                        rtb322 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(gl, (ijk2) as isize)),
                                                            ),
                                                            cadd(
                                                                cmul(
                                                                    (*checked_at(p, (ip) as isize)),
                                                                    (*checked_at(
                                                                        rurf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                        );
                                                        rtb211 = cadd(
                                                            cmul(rtb311, rtb111),
                                                            cmul(rtb312, rtb121),
                                                        );
                                                        rtb212 = cadd(
                                                            cmul(rtb311, rtb112),
                                                            cmul(rtb312, rtb122),
                                                        );
                                                        rtb221 = cadd(
                                                            cmul(rtb321, rtb111),
                                                            cmul(rtb322, rtb121),
                                                        );
                                                        rtb222 = cadd(
                                                            cmul(rtb321, rtb112),
                                                            cmul(rtb322, rtb122),
                                                        );
                                                        (*checked_at(ux, (ip) as isize)) = cadd(
                                                            cmul(rtb211, vuzs1),
                                                            cmul(rtb212, vuzs2),
                                                        );
                                                        (*checked_at(uz, (ip) as isize)) = cadd(
                                                            cmul(rtb221, vuzs1),
                                                            cmul(rtb222, vuzs2),
                                                        );
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        } else {
                                            /* source elastic,receiver elastic or*/
                                            /* source acoustic, receiver elastic */
                                            if (((*checked_at(acoustic, (iz) as isize)) == 0i32)
                                                || ((*checked_at(acoustic, (iz) as isize)) == 3i32))
                                            {
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            det = csub(
                                                                cmul(
                                                                    (*checked_at(
                                                                        turf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        turf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        turf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        turf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb211 = cdiv(
                                                                (*checked_at(
                                                                    turf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                                det,
                                                            );
                                                            rtb212 = cneg(cdiv(
                                                                (*checked_at(
                                                                    turf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                det,
                                                            ));
                                                            rtb221 = cneg(cdiv(
                                                                (*checked_at(
                                                                    turf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                                det,
                                                            ));
                                                            rtb222 = cdiv(
                                                                (*checked_at(
                                                                    turf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                det,
                                                            );
                                                            (*checked_at(rdnr11, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(
                                                                        rtb211,
                                                                        (*checked_at(
                                                                            tusf11,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb212,
                                                                        (*checked_at(
                                                                            tusf21,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                            (*checked_at(rdnr12, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(
                                                                        rtb211,
                                                                        (*checked_at(
                                                                            tusf12,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb212,
                                                                        (*checked_at(
                                                                            tusf22,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                            (*checked_at(rdnr21, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(
                                                                        rtb221,
                                                                        (*checked_at(
                                                                            tusf11,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb222,
                                                                        (*checked_at(
                                                                            tusf21,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                            (*checked_at(rdnr22, (ip) as isize)) =
                                                                cadd(
                                                                    cmul(
                                                                        rtb221,
                                                                        (*checked_at(
                                                                            tusf12,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb222,
                                                                        (*checked_at(
                                                                            tusf22,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            } else {
                                                /* source elastic,receiver acoustic */
                                                /* or source and receiver acoustic */
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            rtb111 = cdiv(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                (*checked_at(
                                                                    turf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            );
                                                            (*checked_at(rdnr11, (ip) as isize)) =
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        tusf11,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                            (*checked_at(rdnr12, (ip) as isize)) =
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        tusf12,
                                                                        (ip) as isize,
                                                                    )),
                                                                );
                                                            (*checked_at(rdnr21, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr22, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            }

                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        ijk2 = (ik2 + ip);
                                                        rtb111 =
                                                            (*checked_at(rdnr11, (ip) as isize));
                                                        rtb112 =
                                                            (*checked_at(rdnr12, (ip) as isize));
                                                        rtb121 =
                                                            (*checked_at(rdnr21, (ip) as isize));
                                                        rtb122 =
                                                            (*checked_at(rdnr22, (ip) as isize));
                                                        rtb211 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb111,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb121,
                                                            ),
                                                        );
                                                        rtb212 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb112,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb122,
                                                            ),
                                                        );
                                                        rtb221 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb111,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb121,
                                                            ),
                                                        );
                                                        rtb222 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb112,
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                                rtb122,
                                                            ),
                                                        );
                                                        rtb311 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf21,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb312 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf12,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf22,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb321 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf11,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf21,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb322 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf12,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rusf22,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        x = csub(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            rtb311,
                                                        );
                                                        y = csub(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            rtb322,
                                                        );
                                                        det =
                                                            csub(cmul(x, y), cmul(rtb312, rtb321));
                                                        rvrb111 = cdiv(y, det);
                                                        rvrb112 = cdiv(rtb312, det);
                                                        rvrb121 = cdiv(rtb321, det);
                                                        rvrb122 = cdiv(x, det);
                                                        fact11 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmad1,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rd0n12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        sigmad2,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                (*checked_at(
                                                                    sigmau1,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        fact12 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n21,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmad1,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rd0n22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        sigmad2,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                (*checked_at(
                                                                    sigmau2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        vuzs1 = cadd(
                                                            cmul(rvrb111, fact11),
                                                            cmul(rvrb112, fact12),
                                                        );
                                                        vuzs2 = cadd(
                                                            cmul(rvrb121, fact11),
                                                            cmul(rvrb122, fact12),
                                                        );
                                                        (*checked_at(up, (ip) as isize)) = cmul(
                                                            (*checked_at(prs, (iz) as isize)),
                                                            cadd(
                                                                cmul(cadd(rtb111, rtb211), vuzs1),
                                                                cmul(cadd(rtb112, rtb212), vuzs2),
                                                            ),
                                                        );
                                                        rtb211 = csub(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                            cmul(
                                                                (*checked_at(gt, (ijk2) as isize)),
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb212 = csub(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                            cmul(
                                                                (*checked_at(gt, (ijk2) as isize)),
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        rtb221 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(gl, (ijk2) as isize)),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                        );
                                                        rtb222 = cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(gl, (ijk2) as isize)),
                                                            ),
                                                            cmul(
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )),
                                                                (*checked_at(p, (ip) as isize)),
                                                            ),
                                                        );
                                                        rtb311 = cadd(
                                                            (*checked_at(p, (ip) as isize)),
                                                            rtb211,
                                                        );
                                                        rtb312 = cadd(
                                                            (*checked_at(gt, (ijk2) as isize)),
                                                            rtb212,
                                                        );
                                                        rtb321 = csub(
                                                            rtb221,
                                                            (*checked_at(gl, (ijk2) as isize)),
                                                        );
                                                        rtb322 = cadd(
                                                            (*checked_at(p, (ip) as isize)),
                                                            rtb222,
                                                        );
                                                        rtb211 = cadd(
                                                            cmul(rtb311, rtb111),
                                                            cmul(rtb312, rtb121),
                                                        );
                                                        rtb212 = cadd(
                                                            cmul(rtb311, rtb112),
                                                            cmul(rtb312, rtb122),
                                                        );
                                                        rtb221 = cadd(
                                                            cmul(rtb321, rtb111),
                                                            cmul(rtb322, rtb121),
                                                        );
                                                        rtb222 = cadd(
                                                            cmul(rtb321, rtb112),
                                                            cmul(rtb322, rtb122),
                                                        );
                                                        (*checked_at(ux, (ip) as isize)) = cadd(
                                                            cmul(rtb211, vuzs1),
                                                            cmul(rtb212, vuzs2),
                                                        );
                                                        (*checked_at(uz, (ip) as isize)) = cadd(
                                                            cmul(rtb221, vuzs1),
                                                            cmul(rtb222, vuzs2),
                                                        );
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        }
                                    } else {
                                        /* receiver below the source */
                                        if ((*checked_at(flag, (iz) as isize)) == 2i32) {
                                            if ((*checked_at(acoustic, (iz) as isize)) == 3i32) {
                                                {
                                                    ip = 0i32;
                                                    while (ip
                                                        < ((*checked_at(ibl, (iz) as isize))
                                                            - 1i32))
                                                    {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            rtb111 = csub(
                                                                (*checked_at(rnr11, (ip) as isize)),
                                                                (*checked_at(
                                                                    rdfr11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            );
                                                            rtb112 = csub(
                                                                (*checked_at(rnr12, (ip) as isize)),
                                                                (*checked_at(
                                                                    rdfr12,
                                                                    (ijk1) as isize,
                                                                )),
                                                            );
                                                            rtb121 = csub(
                                                                (*checked_at(rnr21, (ip) as isize)),
                                                                (*checked_at(
                                                                    rdfr21,
                                                                    (ijk1) as isize,
                                                                )),
                                                            );
                                                            rtb122 = csub(
                                                                (*checked_at(rnr22, (ip) as isize)),
                                                                (*checked_at(
                                                                    rdfr22,
                                                                    (ijk1) as isize,
                                                                )),
                                                            );
                                                            det = csub(
                                                                cmul(rtb111, rtb122),
                                                                cmul(rtb112, rtb121),
                                                            );
                                                            rtb211 = cdiv(rtb122, det);
                                                            rtb212 = cneg(cdiv(rtb112, det));
                                                            rtb221 = cneg(cdiv(rtb121, det));
                                                            rtb222 = cdiv(rtb111, det);
                                                            rtb311 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        turf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb211,
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        turf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb312 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        turf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb211,
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        turf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb321 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        turf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb221,
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        turf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb322 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        turf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb221,
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        turf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb111 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb311,
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb321,
                                                                        (*checked_at(
                                                                            tdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb312,
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb322,
                                                                        (*checked_at(
                                                                            tdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb311,
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb321,
                                                                        (*checked_at(
                                                                            tdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb312,
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb322,
                                                                        (*checked_at(
                                                                            tdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            det = csub(
                                                                cmul(rtb111, rtb122),
                                                                cmul(rtb112, rtb121),
                                                            );
                                                            (*checked_at(rdnr11, (ip) as isize)) =
                                                                cdiv(rtb122, det);
                                                            (*checked_at(rdnr12, (ip) as isize)) =
                                                                cneg(cdiv(rtb112, det));
                                                            (*checked_at(rdnr21, (ip) as isize)) =
                                                                cneg(cdiv(rtb121, det));
                                                            (*checked_at(rdnr22, (ip) as isize)) =
                                                                cdiv(rtb111, det);
                                                            rtb111 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            x = csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                rtb111,
                                                            );
                                                            y = csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                rtb122,
                                                            );
                                                            det = csub(
                                                                cmul(x, y),
                                                                cmul(rtb112, rtb121),
                                                            );
                                                            rvrb111 = cdiv(y, det);
                                                            rvrb112 = cdiv(rtb112, det);
                                                            rvrb121 = cdiv(rtb121, det);
                                                            rvrb122 = cdiv(x, det);
                                                            rtb211 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb111,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb121,
                                                                ),
                                                            );
                                                            rtb212 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb112,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb122,
                                                                ),
                                                            );
                                                            rtb221 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb111,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb121,
                                                                ),
                                                            );
                                                            rtb222 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb112,
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb122,
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr11,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb211,
                                                                    (*checked_at(
                                                                        tdsr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        tdsr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr12,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb211,
                                                                    (*checked_at(
                                                                        tdsr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        tdsr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr21,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb221,
                                                                    (*checked_at(
                                                                        tdsr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        tdsr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr22,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb221,
                                                                    (*checked_at(
                                                                        tdsr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        tdsr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb111 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb211,
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb211,
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb221,
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    rtb221,
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf11,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        rurf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb112,
                                                                        (*checked_at(
                                                                            turf21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf12,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        rurf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb112,
                                                                        (*checked_at(
                                                                            turf22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf21,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb121,
                                                                    (*checked_at(
                                                                        rurf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb122,
                                                                        (*checked_at(
                                                                            turf21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf22,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb121,
                                                                    (*checked_at(
                                                                        rurf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb122,
                                                                        (*checked_at(
                                                                            turf22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                                {
                                                    ip = ((*checked_at(ibl, (iz) as isize)) - 1i32);
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            (*checked_at(rdnr11, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr12, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr21, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr22, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            rtb111 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            x = csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                rtb111,
                                                            );
                                                            y = csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                rtb122,
                                                            );
                                                            det = csub(
                                                                cmul(x, y),
                                                                cmul(rtb112, rtb121),
                                                            );
                                                            rvrb111 = cdiv(y, det);
                                                            rvrb112 = cdiv(rtb112, det);
                                                            rvrb121 = cdiv(rtb121, det);
                                                            rvrb122 = cdiv(x, det);
                                                            rtb211 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb111,
                                                                ),
                                                                cmul(
                                                                    rvrb121,
                                                                    (*checked_at(
                                                                        tdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb212 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb112,
                                                                ),
                                                                cmul(
                                                                    rvrb122,
                                                                    (*checked_at(
                                                                        tdfr12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb221 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb111,
                                                                ),
                                                                cmul(
                                                                    rvrb121,
                                                                    (*checked_at(
                                                                        tdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb222 = cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rvrb112,
                                                                ),
                                                                cmul(
                                                                    rvrb122,
                                                                    (*checked_at(
                                                                        tdfr22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr11,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb211,
                                                                    (*checked_at(
                                                                        tdsr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        tdsr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr12,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb211,
                                                                    (*checked_at(
                                                                        tdsr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        tdsr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr21,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb221,
                                                                    (*checked_at(
                                                                        tdsr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        tdsr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                tdfr22,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb221,
                                                                    (*checked_at(
                                                                        tdsr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        tdsr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb111 = cadd(
                                                                cmul(
                                                                    rtb211,
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb112 = cadd(
                                                                cmul(
                                                                    rtb211,
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb212,
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb121 = cadd(
                                                                cmul(
                                                                    rtb221,
                                                                    (*checked_at(
                                                                        rusr11,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        rusr21,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            rtb122 = cadd(
                                                                cmul(
                                                                    rtb221,
                                                                    (*checked_at(
                                                                        rusr12,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                                cmul(
                                                                    rtb222,
                                                                    (*checked_at(
                                                                        rusr22,
                                                                        (ip) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf11,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        turf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb112,
                                                                        (*checked_at(
                                                                            turf21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf12,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb111,
                                                                    (*checked_at(
                                                                        turf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb112,
                                                                        (*checked_at(
                                                                            turf22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf21,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb121,
                                                                    (*checked_at(
                                                                        turf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb122,
                                                                        (*checked_at(
                                                                            turf21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf21,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(
                                                                rurf22,
                                                                (ijk1) as isize,
                                                            )) = cadd(
                                                                cmul(
                                                                    rtb121,
                                                                    (*checked_at(
                                                                        turf12,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                                cadd(
                                                                    cmul(
                                                                        rtb122,
                                                                        (*checked_at(
                                                                            turf22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    (*checked_at(
                                                                        rurf22,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            } else {
                                                if ((*checked_at(acoustic, (iz) as isize)) == 1i32)
                                                {
                                                    {
                                                        ip = 0i32;
                                                        while (ip
                                                            < ((*checked_at(ibl, (iz) as isize))
                                                                - 1i32))
                                                        {
                                                            {
                                                                ijk1 = (ik1 + ip);
                                                                rtb211 = cdiv(
                                                                    (*checked_at(
                                                                        turf11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    csub(
                                                                        (*checked_at(
                                                                            rd0n11,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    rdnr11,
                                                                    (ip) as isize,
                                                                )) = cdiv(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb211,
                                                                        ),
                                                                        (*checked_at(
                                                                            rurf11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    rdnr12,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr21,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr22,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                    {
                                                        ip = ((*checked_at(ibl, (iz) as isize))
                                                            - 1i32);
                                                        while (ip < block_size) {
                                                            {
                                                                (*checked_at(
                                                                    rdnr11,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr12,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr21,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr22,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                } else {
                                                    if ((*checked_at(acoustic, (iz) as isize))
                                                        == 2i32)
                                                    {
                                                        {
                                                            ip = 0i32;
                                                            while (ip < block_size) {
                                                                {
                                                                    (*checked_at(
                                                                        rdnr11,
                                                                        (ip) as isize,
                                                                    )) = (*checked_at(
                                                                        rnr11,
                                                                        (ip) as isize,
                                                                    ));
                                                                    (*checked_at(
                                                                        rdnr12,
                                                                        (ip) as isize,
                                                                    )) = (*checked_at(
                                                                        rnr12,
                                                                        (ip) as isize,
                                                                    ));
                                                                    (*checked_at(
                                                                        rdnr21,
                                                                        (ip) as isize,
                                                                    )) = (*checked_at(
                                                                        rnr21,
                                                                        (ip) as isize,
                                                                    ));
                                                                    (*checked_at(
                                                                        rdnr22,
                                                                        (ip) as isize,
                                                                    )) = (*checked_at(
                                                                        rnr22,
                                                                        (ip) as isize,
                                                                    ));
                                                                }
                                                                ip += 1;
                                                            }
                                                        }
                                                    } else {
                                                        {
                                                            ip = 0i32;
                                                            while (ip
                                                                < ((*checked_at(
                                                                    ibl,
                                                                    (iz) as isize,
                                                                )) - 1i32))
                                                            {
                                                                {
                                                                    ijk1 = (ik1 + ip);
                                                                    rtb111 = csub(
                                                                        (*checked_at(
                                                                            rd0n11,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    rtb112 = csub(
                                                                        (*checked_at(
                                                                            rd0n12,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    rtb121 = csub(
                                                                        (*checked_at(
                                                                            rd0n21,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    rtb122 = csub(
                                                                        (*checked_at(
                                                                            rd0n22,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    det = csub(
                                                                        cmul(rtb111, rtb122),
                                                                        cmul(rtb112, rtb121),
                                                                    );
                                                                    rtb211 = cdiv(rtb122, det);
                                                                    rtb212 =
                                                                        cneg(cdiv(rtb112, det));
                                                                    rtb221 =
                                                                        cneg(cdiv(rtb121, det));
                                                                    rtb222 = cdiv(rtb111, det);
                                                                    rtb311 = cadd(
                                                                        cmul(
                                                                            rtb211,
                                                                            (*checked_at(
                                                                                turf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                turf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb312 = cadd(
                                                                        cmul(
                                                                            rtb211,
                                                                            (*checked_at(
                                                                                turf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                turf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb321 = cadd(
                                                                        cmul(
                                                                            rtb221,
                                                                            (*checked_at(
                                                                                turf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                turf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb322 = cadd(
                                                                        cmul(
                                                                            rtb221,
                                                                            (*checked_at(
                                                                                turf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                turf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb111 = cadd(
                                                                        cmul(
                                                                            rtb311,
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb321,
                                                                                (*checked_at(
                                                                                    tdfr12,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb112 = cadd(
                                                                        cmul(
                                                                            rtb312,
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb322,
                                                                                (*checked_at(
                                                                                    tdfr12,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb121 = cadd(
                                                                        cmul(
                                                                            rtb311,
                                                                            (*checked_at(
                                                                                tdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb321,
                                                                                (*checked_at(
                                                                                    tdfr22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb122 = cadd(
                                                                        cmul(
                                                                            rtb312,
                                                                            (*checked_at(
                                                                                tdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb322,
                                                                                (*checked_at(
                                                                                    tdfr22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    det = csub(
                                                                        cmul(rtb111, rtb122),
                                                                        cmul(rtb112, rtb121),
                                                                    );
                                                                    (*checked_at(
                                                                        rdnr11,
                                                                        (ip) as isize,
                                                                    )) = cdiv(rtb122, det);
                                                                    (*checked_at(
                                                                        rdnr12,
                                                                        (ip) as isize,
                                                                    )) = cneg(cdiv(rtb112, det));
                                                                    (*checked_at(
                                                                        rdnr21,
                                                                        (ip) as isize,
                                                                    )) = cneg(cdiv(rtb121, det));
                                                                    (*checked_at(
                                                                        rdnr22,
                                                                        (ip) as isize,
                                                                    )) = cdiv(rtb111, det);
                                                                }
                                                                ip += 1;
                                                            }
                                                        }
                                                        {
                                                            ip =
                                                                ((*checked_at(ibl, (iz) as isize))
                                                                    - 1i32);
                                                            while (ip < block_size) {
                                                                {
                                                                    (*checked_at(
                                                                        rdnr11,
                                                                        (ip) as isize,
                                                                    )) = cmplx(
                                                                        (0f64 as f32),
                                                                        (0f64 as f32),
                                                                    );
                                                                    (*checked_at(
                                                                        rdnr12,
                                                                        (ip) as isize,
                                                                    )) = cmplx(
                                                                        (0f64 as f32),
                                                                        (0f64 as f32),
                                                                    );
                                                                    (*checked_at(
                                                                        rdnr21,
                                                                        (ip) as isize,
                                                                    )) = cmplx(
                                                                        (0f64 as f32),
                                                                        (0f64 as f32),
                                                                    );
                                                                    (*checked_at(
                                                                        rdnr22,
                                                                        (ip) as isize,
                                                                    )) = cmplx(
                                                                        (0f64 as f32),
                                                                        (0f64 as f32),
                                                                    );
                                                                }
                                                                ip += 1;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        } else {
                                            if ((*checked_at(flag, (iz) as isize)) == 0i32) {
                                                /* source acoustic, 1st elastic rec */
                                                if ((*checked_at(acoustic, (iz) as isize)) == 2i32)
                                                {
                                                    {
                                                        ip = 0i32;
                                                        while (ip < block_size) {
                                                            {
                                                                (*checked_at(
                                                                    rdnr11,
                                                                    (ip) as isize,
                                                                )) = (*checked_at(
                                                                    rnr11,
                                                                    (ip) as isize,
                                                                ));
                                                                (*checked_at(
                                                                    rdnr12,
                                                                    (ip) as isize,
                                                                )) = (*checked_at(
                                                                    rnr12,
                                                                    (ip) as isize,
                                                                ));
                                                                (*checked_at(
                                                                    rdnr21,
                                                                    (ip) as isize,
                                                                )) = (*checked_at(
                                                                    rnr21,
                                                                    (ip) as isize,
                                                                ));
                                                                (*checked_at(
                                                                    rdnr22,
                                                                    (ip) as isize,
                                                                )) = (*checked_at(
                                                                    rnr22,
                                                                    (ip) as isize,
                                                                ));
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                } else {
                                                    if ((*checked_at(acoustic, (iz) as isize))
                                                        == 3i32)
                                                    {
                                                        {
                                                            ip = 0i32;
                                                            while (ip < block_size) {
                                                                {
                                                                    ijk1 = (ik1 + ip);
                                                                    rtb111 = csub(
                                                                        (*checked_at(
                                                                            rnr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    rtb112 = csub(
                                                                        (*checked_at(
                                                                            rnr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    rtb121 = csub(
                                                                        (*checked_at(
                                                                            rnr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    rtb122 = csub(
                                                                        (*checked_at(
                                                                            rnr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    );
                                                                    det = csub(
                                                                        cmul(rtb111, rtb122),
                                                                        cmul(rtb112, rtb121),
                                                                    );
                                                                    rtb211 = cdiv(rtb122, det);
                                                                    rtb212 =
                                                                        cneg(cdiv(rtb112, det));
                                                                    rtb221 =
                                                                        cneg(cdiv(rtb121, det));
                                                                    rtb222 = cdiv(rtb111, det);
                                                                    rtb311 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                turf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb211,
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                turf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb312 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                turf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb211,
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                turf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb321 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                turf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb221,
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                turf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb322 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                turf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb221,
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                turf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb111 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb311,
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb321,
                                                                                (*checked_at(
                                                                                    tdfr12,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb112 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb312,
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb322,
                                                                                (*checked_at(
                                                                                    tdfr12,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb121 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb311,
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb321,
                                                                                (*checked_at(
                                                                                    tdfr22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb122 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb312,
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb322,
                                                                                (*checked_at(
                                                                                    tdfr22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    det = csub(
                                                                        cmul(rtb111, rtb122),
                                                                        cmul(rtb112, rtb121),
                                                                    );
                                                                    (*checked_at(
                                                                        rdnr11,
                                                                        (ip) as isize,
                                                                    )) = cdiv(rtb122, det);
                                                                    (*checked_at(
                                                                        rdnr12,
                                                                        (ip) as isize,
                                                                    )) = cneg(cdiv(rtb112, det));
                                                                    (*checked_at(
                                                                        rdnr21,
                                                                        (ip) as isize,
                                                                    )) = cneg(cdiv(rtb121, det));
                                                                    (*checked_at(
                                                                        rdnr22,
                                                                        (ip) as isize,
                                                                    )) = cdiv(rtb111, det);
                                                                    rtb111 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr11,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr12,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb112 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr11,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr12,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb121 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr21,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr22,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb122 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr21,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            (*checked_at(
                                                                                rusr22,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    x = csub(
                                                                        cmplx(
                                                                            (1f64 as f32),
                                                                            (0f64 as f32),
                                                                        ),
                                                                        rtb111,
                                                                    );
                                                                    y = csub(
                                                                        cmplx(
                                                                            (1f64 as f32),
                                                                            (0f64 as f32),
                                                                        ),
                                                                        rtb122,
                                                                    );
                                                                    det = csub(
                                                                        cmul(x, y),
                                                                        cmul(rtb112, rtb121),
                                                                    );
                                                                    rvrb111 = cdiv(y, det);
                                                                    rvrb112 = cdiv(rtb112, det);
                                                                    rvrb121 = cdiv(rtb121, det);
                                                                    rvrb122 = cdiv(x, det);
                                                                    rtb211 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rvrb111,
                                                                        ),
                                                                        cmul(
                                                                            rvrb121,
                                                                            (*checked_at(
                                                                                tdfr12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb212 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rvrb112,
                                                                        ),
                                                                        cmul(
                                                                            rvrb122,
                                                                            (*checked_at(
                                                                                tdfr12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb221 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rvrb111,
                                                                        ),
                                                                        cmul(
                                                                            rvrb121,
                                                                            (*checked_at(
                                                                                tdfr22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb222 = cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rvrb112,
                                                                        ),
                                                                        cmul(
                                                                            rvrb122,
                                                                            (*checked_at(
                                                                                tdfr22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb211,
                                                                            (*checked_at(
                                                                                tdsr11,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                tdsr21,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        tdfr12,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb211,
                                                                            (*checked_at(
                                                                                tdsr12,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                tdsr22,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        tdfr21,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb221,
                                                                            (*checked_at(
                                                                                tdsr11,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                tdsr21,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        tdfr22,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb221,
                                                                            (*checked_at(
                                                                                tdsr12,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                tdsr22,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb111 = cadd(
                                                                        cmul(
                                                                            rtb211,
                                                                            (*checked_at(
                                                                                rusr11,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                rusr21,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb112 = cadd(
                                                                        cmul(
                                                                            rtb211,
                                                                            (*checked_at(
                                                                                rusr12,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb212,
                                                                            (*checked_at(
                                                                                rusr22,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb121 = cadd(
                                                                        cmul(
                                                                            rtb221,
                                                                            (*checked_at(
                                                                                rusr11,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                rusr21,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    rtb122 = cadd(
                                                                        cmul(
                                                                            rtb221,
                                                                            (*checked_at(
                                                                                rusr12,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                        cmul(
                                                                            rtb222,
                                                                            (*checked_at(
                                                                                rusr22,
                                                                                (ip) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        rurf11,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb111,
                                                                            (*checked_at(
                                                                                turf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb112,
                                                                                (*checked_at(
                                                                                    turf21,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        rurf12,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb111,
                                                                            (*checked_at(
                                                                                turf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb112,
                                                                                (*checked_at(
                                                                                    turf22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        rurf21,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb121,
                                                                            (*checked_at(
                                                                                turf11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb122,
                                                                                (*checked_at(
                                                                                    turf21,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        rurf22,
                                                                        (ijk1) as isize,
                                                                    )) = cadd(
                                                                        cmul(
                                                                            rtb121,
                                                                            (*checked_at(
                                                                                turf12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        cadd(
                                                                            cmul(
                                                                                rtb122,
                                                                                (*checked_at(
                                                                                    turf22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            (*checked_at(
                                                                                rurf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                    );
                                                                }
                                                                ip += 1;
                                                            }
                                                        }
                                                    } else {
                                                        if ((*checked_at(acoustic, (iz) as isize))
                                                            == 0i32)
                                                        {
                                                            {
                                                                ip = 0i32;
                                                                while (ip < block_size) {
                                                                    {
                                                                        ijk1 = (ik1 + ip);
                                                                        rtb111 = csub(
                                                                            (*checked_at(
                                                                                rd0n11,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        );
                                                                        rtb112 = csub(
                                                                            (*checked_at(
                                                                                rd0n12,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr12,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        );
                                                                        rtb121 = csub(
                                                                            (*checked_at(
                                                                                rd0n21,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        );
                                                                        rtb122 = csub(
                                                                            (*checked_at(
                                                                                rd0n22,
                                                                                (ip) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                rdfr22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        );
                                                                        det = csub(
                                                                            cmul(rtb111, rtb122),
                                                                            cmul(rtb112, rtb121),
                                                                        );
                                                                        rtb211 = cdiv(rtb122, det);
                                                                        rtb212 =
                                                                            cneg(cdiv(rtb112, det));
                                                                        rtb221 =
                                                                            cneg(cdiv(rtb121, det));
                                                                        rtb222 = cdiv(rtb111, det);
                                                                        rtb311 = cadd(
                                                                            cmul(
                                                                                rtb211,
                                                                                (*checked_at(
                                                                                    turf11,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb212,
                                                                                (*checked_at(
                                                                                    turf21,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb312 = cadd(
                                                                            cmul(
                                                                                rtb211,
                                                                                (*checked_at(
                                                                                    turf12,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb212,
                                                                                (*checked_at(
                                                                                    turf22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb321 = cadd(
                                                                            cmul(
                                                                                rtb221,
                                                                                (*checked_at(
                                                                                    turf11,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb222,
                                                                                (*checked_at(
                                                                                    turf21,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb322 = cadd(
                                                                            cmul(
                                                                                rtb221,
                                                                                (*checked_at(
                                                                                    turf12,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                rtb222,
                                                                                (*checked_at(
                                                                                    turf22,
                                                                                    (ijk1) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        rtb111 = cadd(cmul((*checked_at(tdfr11,(ijk1) as isize)),rtb311),cadd(cmul((*checked_at(tdfr12,(ijk1) as isize)),rtb321),(*checked_at(rurf11,(ijk1) as isize))));
                                                                        rtb112 = cadd(cmul((*checked_at(tdfr11,(ijk1) as isize)),rtb312),cadd(cmul((*checked_at(tdfr12,(ijk1) as isize)),rtb322),(*checked_at(rurf12,(ijk1) as isize))));
                                                                        rtb121 = cadd(cmul((*checked_at(tdfr21,(ijk1) as isize)),rtb311),cadd(cmul((*checked_at(tdfr22,(ijk1) as isize)),rtb321),(*checked_at(rurf21,(ijk1) as isize))));
                                                                        rtb122 = cadd(cmul((*checked_at(tdfr21,(ijk1) as isize)),rtb312),cadd(cmul((*checked_at(tdfr22,(ijk1) as isize)),rtb322),(*checked_at(rurf22,(ijk1) as isize))));
                                                                        det = csub(
                                                                            cmul(rtb111, rtb122),
                                                                            cmul(rtb112, rtb121),
                                                                        );
                                                                        (*checked_at(
                                                                            rdnr11,
                                                                            (ip) as isize,
                                                                        )) = cdiv(rtb122, det);
                                                                        (*checked_at(
                                                                            rdnr12,
                                                                            (ip) as isize,
                                                                        )) =
                                                                            cneg(cdiv(rtb112, det));
                                                                        (*checked_at(
                                                                            rdnr21,
                                                                            (ip) as isize,
                                                                        )) =
                                                                            cneg(cdiv(rtb121, det));
                                                                        (*checked_at(
                                                                            rdnr22,
                                                                            (ip) as isize,
                                                                        )) = cdiv(rtb111, det);
                                                                    }
                                                                    ip += 1;
                                                                }
                                                            }
                                                        } else {
                                                            if ((*checked_at(
                                                                acoustic,
                                                                (iz) as isize,
                                                            )) == 1i32)
                                                            {
                                                                {
                                                                    ip = 0i32;
                                                                    while (ip < block_size) {
                                                                        {
                                                                            ijk1 = (ik1 + ip);
                                                                            rtb211 = cdiv((*checked_at(turf11,(ijk1) as isize)),csub((*checked_at(rd0n11,(ip) as isize)),(*checked_at(rdfr11,(ijk1) as isize))));
                                                                            (*checked_at(rdnr11,(ip) as isize)) = cdiv(cmplx((1f64 as f32),(0f64 as f32)),cadd(cmul((*checked_at(tdfr11,(ijk1) as isize)),rtb211),(*checked_at(rurf11,(ijk1) as isize))));
                                                                            (*checked_at(
                                                                                rdnr12,
                                                                                (ip) as isize,
                                                                            )) = cmplx(
                                                                                (0f64 as f32),
                                                                                (0f64 as f32),
                                                                            );
                                                                            (*checked_at(
                                                                                rdnr21,
                                                                                (ip) as isize,
                                                                            )) = cmplx(
                                                                                (0f64 as f32),
                                                                                (0f64 as f32),
                                                                            );
                                                                            (*checked_at(
                                                                                rdnr22,
                                                                                (ip) as isize,
                                                                            )) = cmplx(
                                                                                (0f64 as f32),
                                                                                (0f64 as f32),
                                                                            );
                                                                        }
                                                                        ip += 1;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            } else {
                                                if ((((*checked_at(acoustic, (iz) as isize))
                                                    == 0i32)
                                                    || ((*checked_at(acoustic, (iz) as isize))
                                                        == 1i32))
                                                    || ((*checked_at(acoustic, (iz) as isize))
                                                        == 2i32))
                                                {
                                                    {
                                                        ip = 0i32;
                                                        while (ip < block_size) {
                                                            {
                                                                (*checked_at(
                                                                    rdnr11,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr12,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr21,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr22,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                } else {
                                                    {
                                                        ip = 0i32;
                                                        while (ip < block_size) {
                                                            {
                                                                ijk1 = (ik1 + ip);
                                                                (*checked_at(
                                                                    rdnr11,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr12,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr21,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                (*checked_at(
                                                                    rdnr22,
                                                                    (ip) as isize,
                                                                )) = cmplx(
                                                                    (0f64 as f32),
                                                                    (0f64 as f32),
                                                                );
                                                                rtb111 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb112 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb121 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb122 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        (*checked_at(
                                                                            rusr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                        (*checked_at(
                                                                            rdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                x = csub(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    rtb111,
                                                                );
                                                                y = csub(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    rtb122,
                                                                );
                                                                det = csub(
                                                                    cmul(x, y),
                                                                    cmul(rtb112, rtb121),
                                                                );
                                                                rvrb111 = cdiv(y, det);
                                                                rvrb112 = cdiv(rtb112, det);
                                                                rvrb121 = cdiv(rtb121, det);
                                                                rvrb122 = cdiv(x, det);
                                                                rtb211 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                        rvrb111,
                                                                    ),
                                                                    cmul(
                                                                        rvrb121,
                                                                        (*checked_at(
                                                                            tdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb212 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tdfr11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                        rvrb112,
                                                                    ),
                                                                    cmul(
                                                                        rvrb122,
                                                                        (*checked_at(
                                                                            tdfr12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb221 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tdfr21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                        rvrb111,
                                                                    ),
                                                                    cmul(
                                                                        rvrb121,
                                                                        (*checked_at(
                                                                            tdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb222 = cadd(
                                                                    cmul(
                                                                        (*checked_at(
                                                                            tdfr21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                        rvrb112,
                                                                    ),
                                                                    cmul(
                                                                        rvrb122,
                                                                        (*checked_at(
                                                                            tdfr22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    tdfr11,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb211,
                                                                        (*checked_at(
                                                                            tdsr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb212,
                                                                        (*checked_at(
                                                                            tdsr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    tdfr12,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb211,
                                                                        (*checked_at(
                                                                            tdsr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb212,
                                                                        (*checked_at(
                                                                            tdsr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    tdfr21,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb221,
                                                                        (*checked_at(
                                                                            tdsr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb222,
                                                                        (*checked_at(
                                                                            tdsr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    tdfr22,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb221,
                                                                        (*checked_at(
                                                                            tdsr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb222,
                                                                        (*checked_at(
                                                                            tdsr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb111 = cadd(
                                                                    cmul(
                                                                        rtb211,
                                                                        (*checked_at(
                                                                            rusr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb212,
                                                                        (*checked_at(
                                                                            rusr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb112 = cadd(
                                                                    cmul(
                                                                        rtb211,
                                                                        (*checked_at(
                                                                            rusr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb212,
                                                                        (*checked_at(
                                                                            rusr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb121 = cadd(
                                                                    cmul(
                                                                        rtb221,
                                                                        (*checked_at(
                                                                            rusr11,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb222,
                                                                        (*checked_at(
                                                                            rusr21,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                rtb122 = cadd(
                                                                    cmul(
                                                                        rtb221,
                                                                        (*checked_at(
                                                                            rusr12,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                    cmul(
                                                                        rtb222,
                                                                        (*checked_at(
                                                                            rusr22,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb111,
                                                                        (*checked_at(
                                                                            turf11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cadd(
                                                                        cmul(
                                                                            rtb112,
                                                                            (*checked_at(
                                                                                turf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        (*checked_at(
                                                                            rurf11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    rurf12,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb111,
                                                                        (*checked_at(
                                                                            turf12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cadd(
                                                                        cmul(
                                                                            rtb112,
                                                                            (*checked_at(
                                                                                turf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        (*checked_at(
                                                                            rurf12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    rurf21,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb121,
                                                                        (*checked_at(
                                                                            turf11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cadd(
                                                                        cmul(
                                                                            rtb122,
                                                                            (*checked_at(
                                                                                turf21,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        (*checked_at(
                                                                            rurf21,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    rurf22,
                                                                    (ijk1) as isize,
                                                                )) = cadd(
                                                                    cmul(
                                                                        rtb121,
                                                                        (*checked_at(
                                                                            turf12,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                    cadd(
                                                                        cmul(
                                                                            rtb122,
                                                                            (*checked_at(
                                                                                turf22,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                        ),
                                                                        (*checked_at(
                                                                            rurf22,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        if ((*checked_at(flag, (iz) as isize)) == 2i32) {
                                            (*checked_at(flag, (iz) as isize)) =
                                                (*checked_at(prv, (iz) as isize));
                                        }

                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    ijk1 = (ik1 + ip);
                                                    ijk2 = (ik2 + ip);
                                                    rtb111 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb112 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb121 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf21, (ijk1) as isize)),
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf22, (ijk1) as isize)),
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb122 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf21, (ijk1) as isize)),
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf22, (ijk1) as isize)),
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    x = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb111,
                                                    );
                                                    y = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb122,
                                                    );
                                                    det = csub(cmul(x, y), cmul(rtb112, rtb121));
                                                    rvrb111 = cdiv(y, det);
                                                    rvrb112 = cdiv(rtb112, det);
                                                    rvrb121 = cdiv(rtb121, det);
                                                    rvrb122 = cdiv(x, det);
                                                    rtb111 = cadd(
                                                        cmul(
                                                            rvrb111,
                                                            (*checked_at(tdfr11, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb112,
                                                            (*checked_at(tdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb112 = cadd(
                                                        cmul(
                                                            rvrb111,
                                                            (*checked_at(tdfr12, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb112,
                                                            (*checked_at(tdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb121 = cadd(
                                                        cmul(
                                                            rvrb121,
                                                            (*checked_at(tdfr11, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb122,
                                                            (*checked_at(tdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb122 = cadd(
                                                        cmul(
                                                            rvrb121,
                                                            (*checked_at(tdfr12, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb122,
                                                            (*checked_at(tdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb211 = cadd(
                                                        cmul(
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                            rtb111,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                            rtb121,
                                                        ),
                                                    );
                                                    rtb212 = cadd(
                                                        cmul(
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                            rtb112,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                            rtb122,
                                                        ),
                                                    );
                                                    rtb221 = cadd(
                                                        cmul(
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                            rtb111,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                            rtb121,
                                                        ),
                                                    );
                                                    rtb222 = cadd(
                                                        cmul(
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                            rtb112,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                            rtb122,
                                                        ),
                                                    );
                                                    rtb311 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb312 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb321 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb322 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                        ),
                                                    );
                                                    x = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb311,
                                                    );
                                                    y = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb322,
                                                    );
                                                    det = csub(cmul(x, y), cmul(rtb312, rtb321));
                                                    rvrb111 = cdiv(y, det);
                                                    rvrb112 = cdiv(rtb312, det);
                                                    rvrb121 = cdiv(rtb321, det);
                                                    rvrb122 = cdiv(x, det);
                                                    fact11 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            (*checked_at(sigmau1, (ip) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rusf12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmau2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(sigmad1, (ip) as isize)),
                                                        ),
                                                    );
                                                    fact12 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            (*checked_at(sigmau1, (ip) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rusf22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmau2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(sigmad2, (ip) as isize)),
                                                        ),
                                                    );
                                                    vdzs1 = cadd(
                                                        cmul(rvrb111, fact11),
                                                        cmul(rvrb112, fact12),
                                                    );
                                                    vdzs2 = cadd(
                                                        cmul(rvrb121, fact11),
                                                        cmul(rvrb122, fact12),
                                                    );
                                                    (*checked_at(up, (ip) as isize)) = cmul(
                                                        (*checked_at(prs, (iz) as isize)),
                                                        cadd(
                                                            cmul(cadd(rtb111, rtb211), vdzs1),
                                                            cmul(cadd(rtb112, rtb212), vdzs2),
                                                        ),
                                                    );
                                                    rtb211 = cadd(
                                                        cmul(
                                                            (*checked_at(p, (ip) as isize)),
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(gt, (ijk2) as isize)),
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb212 = cadd(
                                                        cmul(
                                                            (*checked_at(p, (ip) as isize)),
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(gt, (ijk2) as isize)),
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb221 = csub(
                                                        cmul(
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                            (*checked_at(p, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(gl, (ijk2) as isize)),
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb222 = csub(
                                                        cmul(
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                            (*checked_at(p, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(gl, (ijk2) as isize)),
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb311 = cadd(
                                                        (*checked_at(p, (ip) as isize)),
                                                        rtb211,
                                                    );
                                                    rtb312 = csub(
                                                        rtb212,
                                                        (*checked_at(gt, (ijk2) as isize)),
                                                    );
                                                    rtb321 = cadd(
                                                        (*checked_at(gl, (ijk2) as isize)),
                                                        rtb221,
                                                    );
                                                    rtb322 = cadd(
                                                        (*checked_at(p, (ip) as isize)),
                                                        rtb222,
                                                    );
                                                    rtb211 = cadd(
                                                        cmul(rtb311, rtb111),
                                                        cmul(rtb312, rtb121),
                                                    );
                                                    rtb212 = cadd(
                                                        cmul(rtb311, rtb112),
                                                        cmul(rtb312, rtb122),
                                                    );
                                                    rtb221 = cadd(
                                                        cmul(rtb321, rtb111),
                                                        cmul(rtb322, rtb121),
                                                    );
                                                    rtb222 = cadd(
                                                        cmul(rtb321, rtb112),
                                                        cmul(rtb322, rtb122),
                                                    );
                                                    (*checked_at(ux, (ip) as isize)) = cadd(
                                                        cmul(rtb211, vdzs1),
                                                        cmul(rtb212, vdzs2),
                                                    );
                                                    (*checked_at(uz, (ip) as isize)) = cadd(
                                                        cmul(rtb221, vdzs1),
                                                        cmul(rtb222, vdzs2),
                                                    );
                                                }
                                                ip += 1;
                                            }
                                        }
                                    }

                                    /* loop over ranges (traces) */
                                    {
                                        ix = 0i32;
                                        while (ix < nx) {
                                            {
                                                x1 = (bx + (dx * (ix as f32)));
                                                if ((x1 as f64) == 0f64) {
                                                    /* zero offset traces are a special case */
                                                    {
                                                        ip = 1i32;
                                                        while (ip < (block_size - 1i32)) {
                                                            {
                                                                t1 = cdiv(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    cwp_csqrt(
                                                                        (*checked_at(
                                                                            p,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    vos1,
                                                                    (ip) as isize,
                                                                )) = cmul(
                                                                    (*checked_at(
                                                                        up,
                                                                        (ip) as isize,
                                                                    )),
                                                                    t1,
                                                                );
                                                                (*checked_at(
                                                                    vos2,
                                                                    (ip) as isize,
                                                                )) = cmul(
                                                                    (*checked_at(
                                                                        ux,
                                                                        (ip) as isize,
                                                                    )),
                                                                    t1,
                                                                );
                                                                (*checked_at(
                                                                    vos3,
                                                                    (ip) as isize,
                                                                )) = cmul(
                                                                    (*checked_at(
                                                                        uz,
                                                                        (ip) as isize,
                                                                    )),
                                                                    t1,
                                                                );
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                    (*checked_at(vos1, (0i32) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(vos1, (1i32) as isize)) = crmul(
                                                        (*checked_at(vos1, (1i32) as isize)),
                                                        (0.5f64 as f32),
                                                    );
                                                    (*checked_at(
                                                        vos1,
                                                        (block_size - 2i32) as isize,
                                                    )) = crmul(
                                                        (*checked_at(
                                                            vos1,
                                                            (block_size - 2i32) as isize,
                                                        )),
                                                        (0.5f64 as f32),
                                                    );
                                                } else {
                                                    temr = crmul(
                                                        wpie,
                                                        ((((2i32 as f64) * (3.1415926535897931f64))
                                                            * (x1 as f64))
                                                            as f32),
                                                    );
                                                    temr = cdiv(tem, cwp_csqrt(temr));
                                                    sig = crmul(divfac, x1);
                                                    sigh = cmul(sig, cdp);
                                                    {
                                                        ip = 0i32;
                                                        while (ip < block_size) {
                                                            {
                                                                (*checked_at(
                                                                    texp,
                                                                    (ip) as isize,
                                                                )) = cwp_cexp(cmul(
                                                                    sig,
                                                                    (*checked_at(p, (ip) as isize)),
                                                                ));
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                    /******************************************************
                                                    	*		compute the slowness inytegral				*
                                                    	******************************************************/
                                                    if (int_type == 1i32) {
                                                        /* use trapezoidal rule */
                                                        temr =
                                                            crmul(cmul(temr, cdp), (0.5f64 as f32));
                                                        {
                                                            ip = 0i32;
                                                            while (ip < (block_size - 1i32)) {
                                                                {
                                                                    jj = (ip + 1i32);
                                                                    (*checked_at(
                                                                        vos1,
                                                                        (ip) as isize,
                                                                    )) = cmul(
                                                                        temr,
                                                                        cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        vos2,
                                                                        (ip) as isize,
                                                                    )) = cmul(
                                                                        temr,
                                                                        cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    ux,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    ux,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                        ),
                                                                    );
                                                                    (*checked_at(
                                                                        vos3,
                                                                        (ip) as isize,
                                                                    )) = cmul(
                                                                        temr,
                                                                        cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    uz,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    uz,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                        ),
                                                                    );
                                                                }
                                                                ip += 1;
                                                            }
                                                        }
                                                    } else {
                                                        if (int_type == 2i32) {
                                                            /* use first order Filon */
                                                            {
                                                                ip = 0i32;
                                                                while (ip < (block_size - 1i32)) {
                                                                    {
                                                                        jj = (ip + 1i32);
                                                                        t1 = csub(
                                                                            (*checked_at(
                                                                                texp,
                                                                                (jj) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                texp,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        func1 = csub(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        func2 = cmul(
                                                                            csub(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            t1,
                                                                        );
                                                                        func3 = csub(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    ux,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    ux,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        func4 = cmul(
                                                                            csub(
                                                                                (*checked_at(
                                                                                    ux,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    ux,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            t1,
                                                                        );
                                                                        func5 = csub(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    uz,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    uz,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        func6 = cmul(
                                                                            csub(
                                                                                (*checked_at(
                                                                                    uz,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    uz,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            t1,
                                                                        );
                                                                        (*checked_at(
                                                                            vos1,
                                                                            (ip) as isize,
                                                                        )) = cmul(
                                                                            cdiv(cdp, sigh),
                                                                            cmul(
                                                                                csub(
                                                                                    func1,
                                                                                    cdiv(
                                                                                        func2, sigh,
                                                                                    ),
                                                                                ),
                                                                                temr,
                                                                            ),
                                                                        );
                                                                        (*checked_at(
                                                                            vos2,
                                                                            (ip) as isize,
                                                                        )) = cmul(
                                                                            cdiv(cdp, sigh),
                                                                            cmul(
                                                                                csub(
                                                                                    func3,
                                                                                    cdiv(
                                                                                        func4, sigh,
                                                                                    ),
                                                                                ),
                                                                                temr,
                                                                            ),
                                                                        );
                                                                        (*checked_at(
                                                                            vos3,
                                                                            (ip) as isize,
                                                                        )) = cmul(
                                                                            cdiv(cdp, sigh),
                                                                            cmul(
                                                                                csub(
                                                                                    func5,
                                                                                    cdiv(
                                                                                        func6, sigh,
                                                                                    ),
                                                                                ),
                                                                                temr,
                                                                            ),
                                                                        );
                                                                    }
                                                                    ip += 1;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }

                                                /**********************************************************
                                                	* 				update output arrays
                                                	**********************************************************/
                                                {
                                                    ip = 0i32;
                                                    while (ip < (block_size - 1i32)) {
                                                        {
                                                            (*checked_at(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        response1,
                                                                        (iw) as isize,
                                                                    )),
                                                                    (ix) as isize,
                                                                )),
                                                                (iz) as isize,
                                                            )) = cadd(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        (*checked_at(
                                                                            response1,
                                                                            (iw) as isize,
                                                                        )),
                                                                        (ix) as isize,
                                                                    )),
                                                                    (iz) as isize,
                                                                )),
                                                                (*checked_at(vos1, (ip) as isize)),
                                                            );
                                                            (*checked_at(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        response2,
                                                                        (iw) as isize,
                                                                    )),
                                                                    (ix) as isize,
                                                                )),
                                                                (iz) as isize,
                                                            )) = cadd(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        (*checked_at(
                                                                            response2,
                                                                            (iw) as isize,
                                                                        )),
                                                                        (ix) as isize,
                                                                    )),
                                                                    (iz) as isize,
                                                                )),
                                                                (*checked_at(vos2, (ip) as isize)),
                                                            );
                                                            (*checked_at(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        response3,
                                                                        (iw) as isize,
                                                                    )),
                                                                    (ix) as isize,
                                                                )),
                                                                (iz) as isize,
                                                            )) = cadd(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        (*checked_at(
                                                                            response3,
                                                                            (iw) as isize,
                                                                        )),
                                                                        (ix) as isize,
                                                                    )),
                                                                    (iz) as isize,
                                                                )),
                                                                (*checked_at(vos3, (ip) as isize)),
                                                            );
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            }
                                            ix += 1;
                                        }
                                    }
                                }
                                iz += 1;
                            }
                        }
                        /******************************************************************
                        	* 						update loop variables					*
                        	*****************************************************************/
                        nblock -= 1;
                        if (nblock != 0i32) {
                            bp = (*checked_at(p, (block_size - 1i32) as isize)).r;
                            left += 1;
                            if (left > block_size) {
                                left = ((left as i32) - (block_size as i32)) as i32;
                                nblock += 1;
                            }
                        } else {
                            if (int_type == 1i32) {
                                bp = bp1;
                            } else {
                                if (int_type == 2i32) {
                                    bp = (0f64 as f32);
                                }
                            }

                            break;
                        }
                    }
                }
                iw += 1;
            }
        }
        /* last block, exit while loop */
        /* if requested, output some important processing information */
        if ((verbose == 1i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(&oxwave_core::prelude::stderr(),"\nA total of %d frequencies, from %g to %g (Hz)\nand a total of %d ray parameters, from %g to %g (s/km) were processed\n",nw,(wbeg as f64),(wfin as f64),np,(bp as f64),(fp as f64));
        }

        if ((verbose == 2i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(&*outfp,"\nA total of %d frequencies, from %g to %g (Hz)\nand a total of %d ray parameters, from %g to %g (s/km) were processed\n",nw,(wbeg as f64),(wfin as f64),np,(bp as f64),(fp as f64));
        }

        /* free allocated space */
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        /*	free1complex(r11);free1complex(r12);free1complex(r21);free1complex(r22);
         */
        ();
        ();
        ();
        /* free more working space */
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
    }
}

// Refl/lib/shreflect.c: sh_reflectivities
pub(super) unsafe fn raw_sh_reflectivities(
    mut int_type: i32,
    mut verbose: i32,
    mut wtype: i32,
    mut nw: i32,
    mut nlayers: i32,
    mut nx: i32,
    mut layern: i32,
    mut nor: i32,
    mut np: i32,
    mut bp1: f32,
    mut m1: f32,
    mut m2: f32,
    mut m3: f32,
    mut h1: f32,
    mut h2: f32,
    mut lsource: i32,
    mut bx: f32,
    mut dx: f32,
    mut xmax: f32,
    mut decay: f32,
    mut fref: f32,
    mut wrefp: f32,
    mut wrefs: f32,
    mut tsec: f32,
    mut p2w: f32,
    mut fs: f32,
    mut epsp: f32,
    mut epss: f32,
    mut sigp: f32,
    mut sigs: f32,
    mut pw1: f32,
    mut pw2: f32,
    mut pw3: f32,
    mut pw4: f32,
    mut flag: *mut i32,
    mut lobs: *mut i32,
    mut rho: *mut f32,
    mut t: *mut f32,
    mut cl: *mut f32,
    mut ct: *mut f32,
    mut ql: *mut f32,
    mut qt: *mut f32,
    mut response1: *mut *mut *mut Complex,
    mut outfp: *mut CFile,
) {
    let mut arena = Arena::default();
    {
        let mut iw: i32 = 0;
        let mut iz: i32 = 0;
        let mut ip: i32 = 0;
        let mut ix: i32 = 0;
        /* loop counters */
        let mut ijk1: i32 = 0i32;
        let mut ik1: i32 = 0;
        let mut ik2: i32 = 0;
        let mut il: i32 = 0;
        let mut jl: i32 = 0;
        /* auxiliary indices */
        let mut ijk: i32 = 0i32;
        let mut ijk2: i32 = 0;
        let mut jj: i32 = 0i32;
        /* more auxiliary indices */
        let mut lrec: i32 = 0;
        let mut bp: f32 = (0f64 as f32);
        /* beginning ray parameter */
        let mut fp: f32 = 0.0;
        /* final ray parameter */
        let mut dp: f32 = 0.0;
        /* ray parameter increment */
        let mut nblock: i32 = 0;
        /* number of blocks to process */
        let mut block_size: i32 = 0;
        /* size of processing blocks */
        let mut left: i32 = 0;
        /* remainder */
        let mut prv: *mut i32 = std::ptr::null_mut();
        let mut ibl: *mut i32 = std::ptr::null_mut();
        /* scratch arrays */
        let mut x1: f32 = 0.0;
        let mut w: f32 = 0.0;
        let mut wpie: Complex = cmplx(0.0, 0.0);
        let mut cdp: Complex = cmplx(0.0, 0.0);
        /* complex frecuency and ray parameter inc */
        /* complex auxiliary variables */
        let mut x: Complex = cmplx(0.0, 0.0);
        let mut y: Complex = cmplx(0.0, 0.0);
        let mut t1: Complex = cmplx(0.0, 0.0);
        let mut divfac: Complex = cmplx(0.0, 0.0);
        /* complex det,psq,gl1,gt1,gl2,gt2,gam1,gam2,at1,at2,al1,al2; */
        let mut det: Complex = cmplx(0.0, 0.0);
        let mut gt1: Complex = cmplx(0.0, 0.0);
        let mut gt2: Complex = cmplx(0.0, 0.0);
        let mut at1: Complex = cmplx(0.0, 0.0);
        let mut at2: Complex = cmplx(0.0, 0.0);
        let mut tem: Complex = cmplx(0.0, 0.0);
        let mut temr: Complex = cmplx(0.0, 0.0);
        let mut sig: Complex = cmplx(0.0, 0.0);
        let mut sigh: Complex = cmplx(0.0, 0.0);
        let mut rho1: Complex = cmplx(0.0, 0.0);
        let mut rho2: Complex = cmplx(0.0, 0.0);
        let mut rvrb111: Complex = cmplx(0.0, 0.0);
        let mut rvrb112: Complex = cmplx(0.0, 0.0);
        let mut rvrb121: Complex = cmplx(0.0, 0.0);
        let mut rvrb122: Complex = cmplx(0.0, 0.0);
        let mut rvrb211: Complex = cmplx(0.0, 0.0);
        /* complex ewigh11,ewigh22,wrn011,ewrd11; */
        let mut ewigh11: Complex = cmplx(0.0, 0.0);
        let mut wrn011: Complex = cmplx(0.0, 0.0);
        let mut ewrd11: Complex = cmplx(0.0, 0.0);
        let mut ewd0211: Complex = cmplx(0.0, 0.0);
        let mut ewd0212: Complex = cmplx(0.0, 0.0);
        let mut ewtu211: Complex = cmplx(0.0, 0.0);
        let mut rtb111: Complex = cmplx(0.0, 0.0);
        let mut rtb112: Complex = cmplx(0.0, 0.0);
        let mut rtb121: Complex = cmplx(0.0, 0.0);
        let mut rtb122: Complex = cmplx(0.0, 0.0);
        let mut rtb211: Complex = cmplx(0.0, 0.0);
        /* complex rtb212,rtb221,rtb222,rtb311,rtb312,rtb321,rtb322; */
        let mut rtb212: Complex = cmplx(0.0, 0.0);
        let mut rtb311: Complex = cmplx(0.0, 0.0);
        let mut rtb312: Complex = cmplx(0.0, 0.0);
        let mut rtb321: Complex = cmplx(0.0, 0.0);
        let mut rtb322: Complex = cmplx(0.0, 0.0);
        let mut func1: Complex = cmplx(0.0, 0.0);
        let mut func2: Complex = cmplx(0.0, 0.0);
        let mut tun0p11: Complex = cmplx(0.0, 0.0);
        let mut tun0p21: Complex = cmplx(0.0, 0.0);
        let mut rd0np11: Complex = cmplx(0.0, 0.0);
        let mut rd0np12: Complex = cmplx(0.0, 0.0);
        let mut rd0np21: Complex = cmplx(0.0, 0.0);
        let mut rd0np22: Complex = cmplx(0.0, 0.0);
        let mut td0np11: Complex = cmplx(0.0, 0.0);
        let mut td0np12: Complex = cmplx(0.0, 0.0);
        let mut rd11: Complex = cmplx(0.0, 0.0);
        let mut td11: Complex = cmplx(0.0, 0.0);
        let mut tu11: Complex = cmplx(0.0, 0.0);
        let mut ru11: Complex = cmplx(0.0, 0.0);
        let mut fact11: Complex = cmplx(0.0, 0.0);
        let mut fact12: Complex = cmplx(0.0, 0.0);
        let mut vuzs1: Complex = cmplx(0.0, 0.0);
        let mut vuzs2: Complex = cmplx(0.0, 0.0);
        let mut vdzs1: Complex = cmplx(0.0, 0.0);
        let mut vdzs2: Complex = cmplx(0.0, 0.0);
        /* complex scratch arrays */
        let mut rd0n11: *mut Complex = std::ptr::null_mut();
        let mut rd0n12: *mut Complex = std::ptr::null_mut();
        let mut rd0n21: *mut Complex = std::ptr::null_mut();
        let mut rd0n22: *mut Complex = std::ptr::null_mut();
        let mut td0n11: *mut Complex = std::ptr::null_mut();
        let mut td0n12: *mut Complex = std::ptr::null_mut();
        let mut td0n21: *mut Complex = std::ptr::null_mut();
        let mut td0n22: *mut Complex = std::ptr::null_mut();
        let mut tun011: *mut Complex = std::ptr::null_mut();
        let mut tun012: *mut Complex = std::ptr::null_mut();
        let mut tun021: *mut Complex = std::ptr::null_mut();
        let mut tun022: *mut Complex = std::ptr::null_mut();
        let mut run011: *mut Complex = std::ptr::null_mut();
        let mut run012: *mut Complex = std::ptr::null_mut();
        let mut run021: *mut Complex = std::ptr::null_mut();
        let mut run022: *mut Complex = std::ptr::null_mut();
        let mut rusf21: *mut Complex = std::ptr::null_mut();
        let mut rusf22: *mut Complex = std::ptr::null_mut();
        let mut tusf11: *mut Complex = std::ptr::null_mut();
        let mut tusf12: *mut Complex = std::ptr::null_mut();
        let mut tusf21: *mut Complex = std::ptr::null_mut();
        let mut tusf22: *mut Complex = std::ptr::null_mut();
        let mut rdfs11: *mut Complex = std::ptr::null_mut();
        let mut rdfs12: *mut Complex = std::ptr::null_mut();
        let mut rdfs21: *mut Complex = std::ptr::null_mut();
        let mut rdfs22: *mut Complex = std::ptr::null_mut();
        let mut rurf11: *mut Complex = std::ptr::null_mut();
        let mut rurf12: *mut Complex = std::ptr::null_mut();
        let mut rurf21: *mut Complex = std::ptr::null_mut();
        let mut rurf22: *mut Complex = std::ptr::null_mut();
        let mut turf11: *mut Complex = std::ptr::null_mut();
        let mut turf12: *mut Complex = std::ptr::null_mut();
        let mut turf21: *mut Complex = std::ptr::null_mut();
        let mut turf22: *mut Complex = std::ptr::null_mut();
        let mut rdfr11: *mut Complex = std::ptr::null_mut();
        let mut rdfr12: *mut Complex = std::ptr::null_mut();
        let mut rdfr21: *mut Complex = std::ptr::null_mut();
        let mut rdfr22: *mut Complex = std::ptr::null_mut();
        let mut tdfr11: *mut Complex = std::ptr::null_mut();
        let mut tdfr12: *mut Complex = std::ptr::null_mut();
        let mut tdfr21: *mut Complex = std::ptr::null_mut();
        let mut tdfr22: *mut Complex = std::ptr::null_mut();
        let mut r11: *mut Complex = std::ptr::null_mut();
        let mut r12: *mut Complex = std::ptr::null_mut();
        let mut r21: *mut Complex = std::ptr::null_mut();
        let mut r22: *mut Complex = std::ptr::null_mut();
        let mut texp: *mut Complex = std::ptr::null_mut();
        let mut vos1: *mut Complex = std::ptr::null_mut();
        let mut vos2: *mut Complex = std::ptr::null_mut();
        let mut vos3: *mut Complex = std::ptr::null_mut();
        let mut up: *mut Complex = std::ptr::null_mut();
        let mut sigmau1: *mut Complex = std::ptr::null_mut();
        let mut sigmau2: *mut Complex = std::ptr::null_mut();
        let mut sigmad1: *mut Complex = std::ptr::null_mut();
        let mut sigmad2: *mut Complex = std::ptr::null_mut();
        let mut alpha: *mut Complex = std::ptr::null_mut();
        let mut betha: *mut Complex = std::ptr::null_mut();
        let mut al: *mut Complex = std::ptr::null_mut();
        let mut at: *mut Complex = std::ptr::null_mut();
        let mut gl: *mut Complex = std::ptr::null_mut();
        let mut gt: *mut Complex = std::ptr::null_mut();
        let mut gam: *mut Complex = std::ptr::null_mut();
        let mut p: *mut Complex = std::ptr::null_mut();
        let mut pp: *mut Complex = std::ptr::null_mut();
        let mut prs: *mut Complex = std::ptr::null_mut();
        let mut pwin: *mut Complex = std::ptr::null_mut();
        /* complex pointers */
        let mut rusf11: *mut Complex = std::ptr::null_mut();
        let mut rusf12: *mut Complex = std::ptr::null_mut();
        /*
            complex *tusr11,*tusr12,*tusr21,*tusr22;
            complex *rusr11,*rusr12,*rusr21,*rusr22;
            complex *tdsr11,*tdsr12,*tdsr21,*tdsr22;
        */
        let mut rdnr11: *mut Complex = std::ptr::null_mut();
        let mut rdnr12: *mut Complex = std::ptr::null_mut();
        let mut rdnr21: *mut Complex = std::ptr::null_mut();
        let mut rdnr22: *mut Complex = std::ptr::null_mut();
        /* allocate working space */
        /* allocate working space */
        up = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rd0n22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        td0n22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun011 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun012 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun021 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tun022 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run011 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run012 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run021 = arena.alloc::<Complex>((100i32 as usize) as usize);
        run022 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        tusf22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rdfs22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rurf11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rurf12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rurf21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rurf22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        turf22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rdfr22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr11 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr12 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr21 = arena.alloc::<Complex>((300i32 as usize) as usize);
        tdfr22 = arena.alloc::<Complex>((300i32 as usize) as usize);
        rusf21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        rusf22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        vos1 = arena.alloc::<Complex>((100i32 as usize) as usize);
        vos2 = arena.alloc::<Complex>((100i32 as usize) as usize);
        vos3 = arena.alloc::<Complex>((100i32 as usize) as usize);
        r11 = arena.alloc::<Complex>((100i32 as usize) as usize);
        r12 = arena.alloc::<Complex>((100i32 as usize) as usize);
        r21 = arena.alloc::<Complex>((100i32 as usize) as usize);
        r22 = arena.alloc::<Complex>((100i32 as usize) as usize);
        texp = arena.alloc::<Complex>((100i32 as usize) as usize);
        /* allocate other working space */
        prv = arena.alloc::<i32>((nor as usize) as usize);
        ibl = arena.alloc::<i32>((nor as usize) as usize);
        /* ibl=alloc1int(BLK); */
        p = arena.alloc::<Complex>((100i32 as usize) as usize);
        pp = arena.alloc::<Complex>((100i32 as usize) as usize);
        pwin = arena.alloc::<Complex>((100i32 as usize) as usize);
        /* gl=alloc1complex(BLK); */
        gl = arena.alloc::<Complex>(((100i32 * nlayers) as usize) as usize);
        gt = arena.alloc::<Complex>(((100i32 * nlayers) as usize) as usize);
        gam = arena.alloc::<Complex>(((100i32 * nlayers) as usize) as usize);
        alpha = arena.alloc::<Complex>((100i32 as usize) as usize);
        betha = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmau1 = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmau2 = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmad1 = arena.alloc::<Complex>((100i32 as usize) as usize);
        sigmad2 = arena.alloc::<Complex>((100i32 as usize) as usize);
        at = arena.alloc::<Complex>((nlayers as usize) as usize);
        al = arena.alloc::<Complex>((nlayers as usize) as usize);
        /* prs=ealloc1complex(nlayers); */
        prs = arena.alloc::<Complex>((nor as usize) as usize);
        /* initialize pointers */
        rusf11 = alpha;
        rusf12 = betha;
        /*
            tusr11=tun011;
            tusr12=tun012;
            tusr21=tun021;
            tusr22=tun022;
            rusr11=run011;
            rusr12=run012;
            rusr21=run021;
            rusr22=run022;
            tdsr11=td0n11;
            tdsr12=td0n12;
            tdsr21=td0n21;
            tdsr22=td0n22;
        */
        rdnr11 = r11;
        rdnr12 = r12;
        rdnr21 = r21;
        rdnr22 = r22;
        /* assign initial value to working variables */
        tem = cdiv(
            cmplx((1f64 as f32), (0f64 as f32)),
            cwp_cexp(crmul(
                cmplx((0f64 as f32), (1f64 as f32)),
                (((3.1415926535897931f64) / (4i32 as f64)) as f32),
            )),
        );
        ewigh11 = cmplx((0f64 as f32), (0f64 as f32));
        /* ewigh22=cmplx(0.0,0.0); */
        /* initialize complex output arrays */
        {
            iw = 0i32;
            while (iw < nw) {
                {
                    ix = 0i32;
                    while (ix < nx) {
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    (*checked_at(
                                        (*checked_at(
                                            (*checked_at(response1, (iw) as isize)),
                                            (ix) as isize,
                                        )),
                                        (iz) as isize,
                                    )) = cmplx((0f64 as f32), (0f64 as f32));
                                }
                                iz += 1;
                            }
                        }
                        ix += 1;
                    }
                }
                iw += 1;
            }
        }
        /**************************************************************************
        	*					Main loop over frequencies			    *
        	**************************************************************************/
        {
            iw = 0i32;
            while (iw < nw) {
                {
                    w = ((((((iw + 1i32) as f64) * (3.1415926535897931f64)) * (2i32 as f64))
                        / (tsec as f64)) as f32);
                    /* w in radians/sec */
                    wpie = cmplx(w, decay);
                    /* complex w */
                    /* compute frequency-dependent auxiliary variables for reflectivity */
                    raw_compute_w_aux_arrays(
                        wtype,
                        layern,
                        nlayers,
                        nor,
                        lsource,
                        (&raw mut np),
                        (&raw mut block_size),
                        (&raw mut nblock),
                        (&raw mut left),
                        fref,
                        wrefp,
                        wrefs,
                        lobs,
                        tsec,
                        p2w,
                        fs,
                        xmax,
                        w,
                        decay,
                        epsp,
                        epss,
                        sigp,
                        sigs,
                        (&raw mut pw1),
                        (&raw mut pw2),
                        (&raw mut pw3),
                        (&raw mut pw4),
                        (&raw mut dp),
                        (&raw mut fp),
                        wpie,
                        (&raw mut cdp),
                        cl,
                        ct,
                        ql,
                        qt,
                        rho,
                        t,
                        al,
                        at,
                        prs,
                    );
                    /**********************************************************************
                    	* 	loop over slopes to compute the reflectivities (by blocks)		  *
                    	**********************************************************************/
                    while (1i32 != 0) {
                        if ((nblock == 1i32) && (left != 0i32)) {
                            block_size = left;
                        }

                        /* compute p_aux_arrays */
                        raw_compute_p_aux_arrays(
                            wtype,
                            nlayers,
                            lsource,
                            block_size,
                            bp,
                            dp,
                            w,
                            decay,
                            pw1,
                            pw2,
                            pw3,
                            pw4,
                            pwin,
                            m1,
                            m2,
                            m3,
                            h1,
                            h2,
                            wpie,
                            (&raw mut divfac),
                            p,
                            pp,
                            al,
                            at,
                            rho,
                            gl,
                            gt,
                            gam,
                            alpha,
                            betha,
                            sigmad1,
                            sigmad2,
                            sigmau1,
                            sigmau2,
                        );
                        /* if requested, output processing information */
                        if ((verbose == 1i32) || (verbose == 3i32)) {
                            oxwave_core::fprintf!(
                                &oxwave_core::prelude::stderr(),
                                "%3d%6.1f%7.3f%7.3f%6d%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                                iw,
                                (w as f64),
                                (bp as f64),
                                (fp as f64),
                                np,
                                (dp as f64),
                                (pw1 as f64),
                                (pw2 as f64),
                                (pw3 as f64),
                                (pw4 as f64)
                            );
                        }

                        if ((verbose == 2i32) || (verbose == 3i32)) {
                            oxwave_core::fprintf!(
                                &*outfp,
                                "%3d%6.1f%7.3f%7.3f%6d%7.3f%7.3f%7.3f%7.3f%7.3f\n",
                                iw,
                                (w as f64),
                                (bp as f64),
                                (fp as f64),
                                np,
                                (dp as f64),
                                (pw1 as f64),
                                (pw2 as f64),
                                (pw3 as f64),
                                (pw4 as f64)
                            );
                        }

                        /******************************************************************
                        	* 			start the computation of the reflectivities 		*
                        	******************************************************************/
                        /* loop over reflecting layers */
                        /* for (il=0; il<nlayers; il++) { */
                        ijk = 0i32;
                        {
                            il = 0i32;
                            while (il < (nlayers - 1i32)) {
                                {
                                    jl = (il + 1i32);
                                    ik1 = (il * block_size);
                                    ik2 = ((il + 1i32) * block_size);
                                    /*
                                                    al1=al[il];
                                                    al2=al[jl];
                                    */
                                    at1 = (*checked_at(at, (il) as isize));
                                    at2 = (*checked_at(at, (jl) as isize));
                                    rho1 = cmplx((*checked_at(rho, (il) as isize)), (0f64 as f32));
                                    rho2 = cmplx((*checked_at(rho, (jl) as isize)), (0f64 as f32));
                                    t1 = cmplx((*checked_at(t, (il) as isize)), (0f64 as f32));
                                    if (il == 0i32) {
                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    ijk1 = (ik1 + ip);
                                                    gt1 = (*checked_at(gt, (ijk1) as isize));
                                                    ewigh11 = cwp_cexp(cmul(
                                                        cmplx((0f64 as f32), (1f64 as f32)),
                                                        cmul(wpie, cmul(gt1, t1)),
                                                    ));
                                                    (*checked_at(rd0n11, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n12, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n22, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    /* free surface */
                                                    (*checked_at(td0n11, (ip) as isize)) = ewigh11;
                                                    (*checked_at(td0n12, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n22, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun011, (ip) as isize)) = ewigh11;
                                                    (*checked_at(tun012, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun021, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun022, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run011, (ip) as isize)) =
                                                        cmul(ewigh11, ewigh11);
                                                    (*checked_at(run012, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run021, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run022, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                }
                                                ip += 1;
                                            }
                                        }
                                    } else {
                                        if (((at1).r == (at2).r) && ((at1).i == (at2).i)) {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        gt1 = (*checked_at(gt, (ijk1) as isize));
                                                        ewigh11 = cwp_cexp(cmul(
                                                            cmplx((0f64 as f32), (1f64 as f32)),
                                                            cmul(wpie, cmul(gt1, t1)),
                                                        ));
                                                        tun0p11 = cmul(
                                                            (*checked_at(tun011, (ip) as isize)),
                                                            ewigh11,
                                                        );
                                                        tun0p21 = cmul(
                                                            (*checked_at(tun021, (ip) as isize)),
                                                            ewigh11,
                                                        );
                                                        /* -R/T coefficients betha1=betha2 */
                                                        rtb111 = cmul(
                                                            ewigh11,
                                                            (*checked_at(run011, (ip) as isize)),
                                                        );
                                                        (*checked_at(run011, (ip) as isize)) =
                                                            cmul(rtb111, ewigh11);
                                                        (*checked_at(run012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        td0np11 = cmul(
                                                            ewigh11,
                                                            (*checked_at(td0n11, (ip) as isize)),
                                                        );
                                                        td0np12 = cmul(
                                                            ewigh11,
                                                            (*checked_at(td0n12, (ip) as isize)),
                                                        );
                                                        (*checked_at(tun011, (ip) as isize)) =
                                                            tun0p11;
                                                        (*checked_at(tun012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun021, (ip) as isize)) =
                                                            tun0p21;
                                                        (*checked_at(tun022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n11, (ip) as isize)) =
                                                            td0np11;
                                                        (*checked_at(td0n12, (ip) as isize)) =
                                                            td0np12;
                                                        (*checked_at(td0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        } else {
                                            {
                                                ip = 0i32;
                                                while (ip < block_size) {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        ijk2 = (ik2 + ip);
                                                        /*							psq=pp[ip];
                                                                                gl1=gl[ijk1];
                                                        */
                                                        gt1 = (*checked_at(gt, (ijk1) as isize));
                                                        /*
                                                                                    gam1=gam[ijk1];
                                                                                    gl2=gl[ijk2];
                                                        */
                                                        gt2 = (*checked_at(gt, (ijk2) as isize));
                                                        /*
                                                                    gam2=gam[ijk2];
                                                        */
                                                        ewigh11 = cwp_cexp(cmul(
                                                            wpie,
                                                            cmul(
                                                                cmplx((0f64 as f32), (1f64 as f32)),
                                                                cmul(t1, gt1),
                                                            ),
                                                        ));
                                                        y = cadd(
                                                            cmul(rho1, cdiv(gt1, cmul(at1, at1))),
                                                            cmul(rho2, cdiv(gt2, cmul(at2, at2))),
                                                        );
                                                        td11 = crmul(
                                                            cmul(
                                                                rho1,
                                                                cdiv(gt1, cmul(y, cmul(at1, at1))),
                                                            ),
                                                            (2f64 as f32),
                                                        );
                                                        rd11 = cdiv(
                                                            csub(
                                                                cmul(rho1, cdiv(at1, at1)),
                                                                cmul(
                                                                    rho1,
                                                                    cdiv(gt2, cmul(at2, at2)),
                                                                ),
                                                            ),
                                                            y,
                                                        );
                                                        tu11 = crmul(
                                                            cmul(
                                                                rho2,
                                                                cdiv(gt2, cmul(y, cmul(at2, at2))),
                                                            ),
                                                            (2f64 as f32),
                                                        );
                                                        ru11 = cneg(rd11);
                                                        ewrd11 = cmul(ewigh11, rd11);
                                                        wrn011 = cmul(
                                                            ewigh11,
                                                            (*checked_at(run011, (ip) as isize)),
                                                        );
                                                        rvrb111 = cdiv(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            csub(
                                                                cmplx((1f64 as f32), (0f64 as f32)),
                                                                cmul(wrn011, ewrd11),
                                                            ),
                                                        );
                                                        rvrb211 = cadd(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            cmul(rvrb111, cmul(ewrd11, wrn011)),
                                                        );
                                                        ewd0211 = cmul(
                                                            (*checked_at(td0n11, (ip) as isize)),
                                                            cmul(rvrb111, ewigh11),
                                                        );
                                                        ewd0212 = cmul(
                                                            (*checked_at(td0n12, (ip) as isize)),
                                                            cmul(rvrb111, ewigh11),
                                                        );
                                                        (*checked_at(td0n11, (ip) as isize)) =
                                                            cmul(td11, ewd0211);
                                                        (*checked_at(td0n12, (ip) as isize)) =
                                                            cmul(td11, ewd0212);
                                                        (*checked_at(td0n21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(td0n22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        ewtu211 =
                                                            cmul(rvrb211, cmul(ewigh11, tu11));
                                                        tun0p11 = cmul(
                                                            ewtu211,
                                                            (*checked_at(tun011, (ip) as isize)),
                                                        );
                                                        tun0p21 = cmul(
                                                            ewtu211,
                                                            (*checked_at(tun021, (ip) as isize)),
                                                        );
                                                        rtb111 = cmul(ewrd11, ewd0211);
                                                        rtb112 = cmul(ewrd11, ewd0212);
                                                        rd0np11 = cadd(
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                            cmul(
                                                                (*checked_at(
                                                                    tun011,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb111,
                                                            ),
                                                        );
                                                        rd0np12 = cadd(
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                            cmul(
                                                                (*checked_at(
                                                                    tun011,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb112,
                                                            ),
                                                        );
                                                        rd0np21 = cadd(
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                            cmul(
                                                                (*checked_at(
                                                                    tun021,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb111,
                                                            ),
                                                        );
                                                        rd0np22 = cadd(
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                            cmul(
                                                                (*checked_at(
                                                                    tun022,
                                                                    (ip) as isize,
                                                                )),
                                                                rtb112,
                                                            ),
                                                        );
                                                        (*checked_at(run011, (ip) as isize)) = cadd(
                                                            ru11,
                                                            cmul(ewtu211, cmul(wrn011, td11)),
                                                        );
                                                        (*checked_at(run012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run021, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(run022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun011, (ip) as isize)) =
                                                            tun0p11;
                                                        (*checked_at(tun012, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(tun021, (ip) as isize)) =
                                                            tun0p21;
                                                        (*checked_at(tun022, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rd0n11, (ip) as isize)) =
                                                            rd0np11;
                                                        (*checked_at(rd0n12, (ip) as isize)) =
                                                            rd0np12;
                                                        (*checked_at(rd0n21, (ip) as isize)) =
                                                            rd0np21;
                                                        (*checked_at(rd0n22, (ip) as isize)) =
                                                            rd0np22;
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        }
                                    }

                                    /* if (il==lobs[ijk]-1) { */
                                    if (il == ((*checked_at(lobs, (ijk) as isize)) - 2i32)) {
                                        ik1 = (ijk * block_size);
                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    iz = (ik1 + ip);
                                                    (*checked_at(rurf11, (iz) as isize)) =
                                                        (*checked_at(run011, (ip) as isize));
                                                    (*checked_at(rurf12, (iz) as isize)) =
                                                        (*checked_at(run012, (ip) as isize));
                                                    (*checked_at(rurf21, (iz) as isize)) =
                                                        (*checked_at(run021, (ip) as isize));
                                                    (*checked_at(rurf22, (iz) as isize)) =
                                                        (*checked_at(run022, (ip) as isize));
                                                    (*checked_at(rdfr11, (iz) as isize)) =
                                                        (*checked_at(rd0n11, (ip) as isize));
                                                    (*checked_at(rdfr12, (iz) as isize)) =
                                                        (*checked_at(rd0n12, (ip) as isize));
                                                    (*checked_at(rdfr21, (iz) as isize)) =
                                                        (*checked_at(rd0n21, (ip) as isize));
                                                    (*checked_at(rdfr22, (iz) as isize)) =
                                                        (*checked_at(rd0n22, (ip) as isize));
                                                    (*checked_at(turf11, (iz) as isize)) =
                                                        (*checked_at(tun011, (ip) as isize));
                                                    (*checked_at(turf12, (iz) as isize)) =
                                                        (*checked_at(tun012, (ip) as isize));
                                                    (*checked_at(turf21, (iz) as isize)) =
                                                        (*checked_at(tun021, (ip) as isize));
                                                    (*checked_at(turf22, (iz) as isize)) =
                                                        (*checked_at(tun022, (ip) as isize));
                                                    (*checked_at(tdfr11, (iz) as isize)) =
                                                        (*checked_at(td0n11, (ip) as isize));
                                                    (*checked_at(tdfr12, (iz) as isize)) =
                                                        (*checked_at(td0n12, (ip) as isize));
                                                    (*checked_at(tdfr21, (iz) as isize)) =
                                                        (*checked_at(td0n21, (ip) as isize));
                                                    (*checked_at(tdfr22, (iz) as isize)) =
                                                        (*checked_at(td0n22, (ip) as isize));
                                                }
                                                ip += 1;
                                            }
                                        }
                                        ijk += 1;
                                    }

                                    /* if (il==lsource-1) { */
                                    if (il == (lsource - 2i32)) {
                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    (*checked_at(rusf11, (ip) as isize)) =
                                                        (*checked_at(run011, (ip) as isize));
                                                    (*checked_at(rusf12, (ip) as isize)) =
                                                        (*checked_at(run012, (ip) as isize));
                                                    (*checked_at(rusf21, (ip) as isize)) =
                                                        (*checked_at(run021, (ip) as isize));
                                                    (*checked_at(rusf22, (ip) as isize)) =
                                                        (*checked_at(run022, (ip) as isize));
                                                    (*checked_at(rdfs11, (ip) as isize)) =
                                                        (*checked_at(rd0n11, (ip) as isize));
                                                    (*checked_at(rdfs12, (ip) as isize)) =
                                                        (*checked_at(rd0n12, (ip) as isize));
                                                    (*checked_at(rdfs21, (ip) as isize)) =
                                                        (*checked_at(rd0n21, (ip) as isize));
                                                    (*checked_at(rdfs22, (ip) as isize)) =
                                                        (*checked_at(rd0n22, (ip) as isize));
                                                    (*checked_at(run011, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run012, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run021, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(run022, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n11, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n12, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rd0n22, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tusf11, (ip) as isize)) =
                                                        (*checked_at(tun011, (ip) as isize));
                                                    (*checked_at(tusf12, (ip) as isize)) =
                                                        (*checked_at(tun012, (ip) as isize));
                                                    (*checked_at(tusf21, (ip) as isize)) =
                                                        (*checked_at(tun021, (ip) as isize));
                                                    (*checked_at(tusf22, (ip) as isize)) =
                                                        (*checked_at(tun022, (ip) as isize));
                                                    (*checked_at(tun011, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun012, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun021, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(tun022, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n11, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n12, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(td0n22, (ip) as isize)) =
                                                        cmplx((1f64 as f32), (0f64 as f32));
                                                }
                                                ip += 1;
                                            }
                                        }
                                    }
                                }
                                il += 1;
                            }
                        }
                        /* loop over receiver depths */
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    ik1 = (iz * block_size);
                                    {
                                        ip = 0i32;
                                        while (ip < block_size) {
                                            {
                                                ijk1 = (ik1 + ip);
                                                det = csub(
                                                    (*checked_at(rd0n11, (ip) as isize)),
                                                    (*checked_at(rdfr11, (ijk1) as isize)),
                                                );
                                                if ((rcabs(det) as f64)
                                                    <= 9.9999999999999998E-13f64)
                                                {
                                                    (*checked_at(prv, (iz) as isize)) =
                                                        (*checked_at(flag, (iz) as isize));
                                                    (*checked_at(flag, (iz) as isize)) = 2i32;
                                                    (*checked_at(ibl, (iz) as isize)) = ip;
                                                    break;
                                                }
                                            }
                                            ip += 1;
                                        }
                                    }
                                }
                                iz += 1;
                            }
                        }
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    ik1 = (iz * block_size);
                                    lrec = (*checked_at(lobs, (iz) as isize));
                                    ik2 = ((lrec - 1i32) * block_size);
                                    /* receiver above the source */
                                    if (lsource > lrec) {
                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    ijk1 = (ik1 + ip);
                                                    ijk2 = (ik2 + ip);
                                                    rtb111 = cdiv(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        (*checked_at(turf11, (ijk1) as isize)),
                                                    );
                                                    (*checked_at(rdnr11, (ip) as isize)) = cmul(
                                                        rtb111,
                                                        (*checked_at(tusf11, (ip) as isize)),
                                                    );
                                                    (*checked_at(rdnr12, (ip) as isize)) = cmul(
                                                        rtb111,
                                                        (*checked_at(tusf12, (ip) as isize)),
                                                    );
                                                    (*checked_at(rdnr21, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(rdnr22, (ip) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    rtb111 = (*checked_at(rdnr11, (ip) as isize));
                                                    rtb112 = (*checked_at(rdnr12, (ip) as isize));
                                                    rtb121 = (*checked_at(rdnr21, (ip) as isize));
                                                    rtb122 = (*checked_at(rdnr22, (ip) as isize));
                                                    rtb211 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            rtb111,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            rtb121,
                                                        ),
                                                    );
                                                    rtb212 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            rtb112,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            rtb122,
                                                        ),
                                                    );
                                                    /*
                                                                            rtb221=cadd(cmul(rurf21[ijk1],rtb111),cmul(rurf22[ijk1],
                                                                                rtb121));
                                                                            rtb222=cadd(cmul(rurf21[ijk1],rtb112),cmul(rurf22[ijk1],
                                                                                rtb122));
                                                    */
                                                    rtb311 = cadd(
                                                        cmul(
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb312 = cadd(
                                                        cmul(
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb321 = cadd(
                                                        cmul(
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb322 = cadd(
                                                        cmul(
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                        ),
                                                    );
                                                    x = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb311,
                                                    );
                                                    y = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb322,
                                                    );
                                                    det = csub(cmul(x, y), cmul(rtb312, rtb321));
                                                    rvrb111 = cdiv(y, det);
                                                    rvrb112 = cdiv(rtb312, det);
                                                    rvrb121 = cdiv(rtb321, det);
                                                    rvrb122 = cdiv(x, det);
                                                    fact11 = cadd(
                                                        cmul(
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                            (*checked_at(sigmad1, (ip) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmad2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(sigmau1, (ip) as isize)),
                                                        ),
                                                    );
                                                    fact12 = cadd(
                                                        cmul(
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                            (*checked_at(sigmad1, (ip) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rd0n22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmad2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(sigmau2, (ip) as isize)),
                                                        ),
                                                    );
                                                    vuzs1 = cadd(
                                                        cmul(rvrb111, fact11),
                                                        cmul(rvrb112, fact12),
                                                    );
                                                    vuzs2 = cadd(
                                                        cmul(rvrb121, fact11),
                                                        cmul(rvrb122, fact12),
                                                    );
                                                    (*checked_at(up, (ip) as isize)) = cmul(
                                                        (*checked_at(p, (ip) as isize)),
                                                        cadd(
                                                            cmul(cadd(rtb111, rtb211), vuzs1),
                                                            cmul(cadd(rtb112, rtb212), vuzs2),
                                                        ),
                                                    );
                                                }
                                                ip += 1;
                                            }
                                        }
                                    } else {
                                        /* receiver below the source */
                                        if ((*checked_at(flag, (iz) as isize)) == 2i32) {
                                            {
                                                ip = 0i32;
                                                while (ip
                                                    < ((*checked_at(ibl, (iz) as isize)) - 1i32))
                                                {
                                                    {
                                                        ijk1 = (ik1 + ip);
                                                        rtb211 = cdiv(
                                                            (*checked_at(turf11, (ijk1) as isize)),
                                                            csub(
                                                                (*checked_at(
                                                                    rd0n11,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    rdfr11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        (*checked_at(rdnr11, (ip) as isize)) = cdiv(
                                                            cmplx((1f64 as f32), (0f64 as f32)),
                                                            cadd(
                                                                cmul(
                                                                    (*checked_at(
                                                                        tdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                    rtb211,
                                                                ),
                                                                (*checked_at(
                                                                    rurf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                            ),
                                                        );
                                                        (*checked_at(rdnr12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rdnr21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rdnr22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                            {
                                                ip = ((*checked_at(ibl, (iz) as isize)) - 1i32);
                                                while (ip < block_size) {
                                                    {
                                                        (*checked_at(rdnr11, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rdnr12, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rdnr21, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                        (*checked_at(rdnr22, (ip) as isize)) =
                                                            cmplx((0f64 as f32), (0f64 as f32));
                                                    }
                                                    ip += 1;
                                                }
                                            }
                                        } else {
                                            if ((*checked_at(flag, (iz) as isize)) == 0i32) {
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            ijk1 = (ik1 + ip);
                                                            rtb211 = cdiv(
                                                                (*checked_at(
                                                                    turf11,
                                                                    (ijk1) as isize,
                                                                )),
                                                                csub(
                                                                    (*checked_at(
                                                                        rd0n11,
                                                                        (ip) as isize,
                                                                    )),
                                                                    (*checked_at(
                                                                        rdfr11,
                                                                        (ijk1) as isize,
                                                                    )),
                                                                ),
                                                            );
                                                            (*checked_at(rdnr11, (ip) as isize)) =
                                                                cdiv(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    cadd(
                                                                        cmul(
                                                                            (*checked_at(
                                                                                tdfr11,
                                                                                (ijk1) as isize,
                                                                            )),
                                                                            rtb211,
                                                                        ),
                                                                        (*checked_at(
                                                                            rurf11,
                                                                            (ijk1) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                            (*checked_at(rdnr12, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr21, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr22, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            } else {
                                                {
                                                    ip = 0i32;
                                                    while (ip < block_size) {
                                                        {
                                                            (*checked_at(rdnr11, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr12, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr21, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                            (*checked_at(rdnr22, (ip) as isize)) =
                                                                cmplx((0f64 as f32), (0f64 as f32));
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            }
                                        }

                                        if ((*checked_at(flag, (iz) as isize)) == 2i32) {
                                            (*checked_at(flag, (iz) as isize)) =
                                                (*checked_at(prv, (iz) as isize));
                                        }

                                        {
                                            ip = 0i32;
                                            while (ip < block_size) {
                                                {
                                                    ijk1 = (ik1 + ip);
                                                    ijk2 = (ik2 + ip);
                                                    rtb111 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb112 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf11, (ijk1) as isize)),
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf12, (ijk1) as isize)),
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb121 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf21, (ijk1) as isize)),
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf22, (ijk1) as isize)),
                                                            (*checked_at(rdnr21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb122 = cadd(
                                                        cmul(
                                                            (*checked_at(rurf21, (ijk1) as isize)),
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rurf22, (ijk1) as isize)),
                                                            (*checked_at(rdnr22, (ip) as isize)),
                                                        ),
                                                    );
                                                    x = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb111,
                                                    );
                                                    y = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb122,
                                                    );
                                                    det = csub(cmul(x, y), cmul(rtb112, rtb121));
                                                    rvrb111 = cdiv(y, det);
                                                    rvrb112 = cdiv(rtb112, det);
                                                    rvrb121 = cdiv(rtb121, det);
                                                    rvrb122 = cdiv(x, det);
                                                    rtb111 = cadd(
                                                        cmul(
                                                            rvrb111,
                                                            (*checked_at(tdfr11, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb112,
                                                            (*checked_at(tdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb112 = cadd(
                                                        cmul(
                                                            rvrb111,
                                                            (*checked_at(tdfr12, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb112,
                                                            (*checked_at(tdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb121 = cadd(
                                                        cmul(
                                                            rvrb121,
                                                            (*checked_at(tdfr11, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb122,
                                                            (*checked_at(tdfr21, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb122 = cadd(
                                                        cmul(
                                                            rvrb121,
                                                            (*checked_at(tdfr12, (ijk1) as isize)),
                                                        ),
                                                        cmul(
                                                            rvrb122,
                                                            (*checked_at(tdfr22, (ijk1) as isize)),
                                                        ),
                                                    );
                                                    rtb211 = cadd(
                                                        cmul(
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                            rtb111,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                            rtb121,
                                                        ),
                                                    );
                                                    rtb212 = cadd(
                                                        cmul(
                                                            (*checked_at(rdnr11, (ip) as isize)),
                                                            rtb112,
                                                        ),
                                                        cmul(
                                                            (*checked_at(rdnr12, (ip) as isize)),
                                                            rtb122,
                                                        ),
                                                    );
                                                    /*
                                                                            rtb221=cadd(cmul(rdnr21[ip],rtb111),cmul(rdnr22[ip],
                                                                                rtb121));
                                                                            rtb222=cadd(cmul(rdnr21[ip],rtb112),cmul(rdnr22[ip],
                                                                                rtb122));
                                                    */
                                                    rtb311 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb312 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf12, (ip) as isize)),
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb321 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            (*checked_at(rd0n11, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                            (*checked_at(rd0n21, (ip) as isize)),
                                                        ),
                                                    );
                                                    rtb322 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            (*checked_at(rd0n12, (ip) as isize)),
                                                        ),
                                                        cmul(
                                                            (*checked_at(rusf22, (ip) as isize)),
                                                            (*checked_at(rd0n22, (ip) as isize)),
                                                        ),
                                                    );
                                                    x = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb311,
                                                    );
                                                    y = csub(
                                                        cmplx((1f64 as f32), (0f64 as f32)),
                                                        rtb322,
                                                    );
                                                    det = csub(cmul(x, y), cmul(rtb312, rtb321));
                                                    rvrb111 = cdiv(y, det);
                                                    rvrb112 = cdiv(rtb312, det);
                                                    rvrb121 = cdiv(rtb321, det);
                                                    rvrb122 = cdiv(x, det);
                                                    fact11 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf11, (ip) as isize)),
                                                            (*checked_at(sigmau1, (ip) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rusf12,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmau2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(sigmad1, (ip) as isize)),
                                                        ),
                                                    );
                                                    fact12 = cadd(
                                                        cmul(
                                                            (*checked_at(rusf21, (ip) as isize)),
                                                            (*checked_at(sigmau1, (ip) as isize)),
                                                        ),
                                                        cadd(
                                                            cmul(
                                                                (*checked_at(
                                                                    rusf22,
                                                                    (ip) as isize,
                                                                )),
                                                                (*checked_at(
                                                                    sigmau2,
                                                                    (ip) as isize,
                                                                )),
                                                            ),
                                                            (*checked_at(sigmad2, (ip) as isize)),
                                                        ),
                                                    );
                                                    vdzs1 = cadd(
                                                        cmul(rvrb111, fact11),
                                                        cmul(rvrb112, fact12),
                                                    );
                                                    vdzs2 = cadd(
                                                        cmul(rvrb121, fact11),
                                                        cmul(rvrb122, fact12),
                                                    );
                                                    (*checked_at(up, (ip) as isize)) = cmul(
                                                        (*checked_at(p, (ip) as isize)),
                                                        cadd(
                                                            cmul(cadd(rtb111, rtb211), vdzs1),
                                                            cmul(cadd(rtb112, rtb212), vdzs2),
                                                        ),
                                                    );
                                                }
                                                ip += 1;
                                            }
                                        }
                                    }

                                    {
                                        ix = 0i32;
                                        while (ix < nx) {
                                            {
                                                x1 = (bx + (dx * (ix as f32)));
                                                if ((x1 as f64) == 0f64) {
                                                    {
                                                        ip = 1i32;
                                                        while (ip < (block_size - 1i32)) {
                                                            {
                                                                t1 = cdiv(
                                                                    cmplx(
                                                                        (1f64 as f32),
                                                                        (0f64 as f32),
                                                                    ),
                                                                    cwp_csqrt(
                                                                        (*checked_at(
                                                                            p,
                                                                            (ip) as isize,
                                                                        )),
                                                                    ),
                                                                );
                                                                (*checked_at(
                                                                    vos1,
                                                                    (ip) as isize,
                                                                )) = cmul(
                                                                    (*checked_at(
                                                                        up,
                                                                        (ip) as isize,
                                                                    )),
                                                                    t1,
                                                                );
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                    (*checked_at(vos1, (0i32) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                    (*checked_at(vos1, (1i32) as isize)) = crmul(
                                                        (*checked_at(vos1, (1i32) as isize)),
                                                        (0.5f64 as f32),
                                                    );
                                                    (*checked_at(
                                                        vos1,
                                                        (block_size - 2i32) as isize,
                                                    )) = crmul(
                                                        (*checked_at(
                                                            vos1,
                                                            (block_size - 2i32) as isize,
                                                        )),
                                                        (0.5f64 as f32),
                                                    );
                                                } else {
                                                    temr = crmul(
                                                        wpie,
                                                        ((((2i32 as f64) * (3.1415926535897931f64))
                                                            * (x1 as f64))
                                                            as f32),
                                                    );
                                                    temr = cdiv(tem, cwp_csqrt(temr));
                                                    sig = crmul(divfac, x1);
                                                    sigh = cmul(sig, cdp);
                                                    {
                                                        ip = 0i32;
                                                        while (ip < (block_size - 1i32)) {
                                                            {
                                                                (*checked_at(
                                                                    texp,
                                                                    (ip) as isize,
                                                                )) = cwp_cexp(cmul(
                                                                    sig,
                                                                    (*checked_at(p, (ip) as isize)),
                                                                ));
                                                            }
                                                            ip += 1;
                                                        }
                                                    }
                                                    /******************************************************
                                                    	*			Compute the slowness integral			  *
                                                    	******************************************************/
                                                    if (int_type == 1i32) {
                                                        /* use trapezoidal rule */
                                                        temr =
                                                            crmul(cmul(temr, cdp), (0.5f64 as f32));
                                                        {
                                                            ip = 0i32;
                                                            while (ip < (block_size - 1i32)) {
                                                                {
                                                                    jj = (ip + 1i32);
                                                                    (*checked_at(
                                                                        vos1,
                                                                        (ip) as isize,
                                                                    )) = cmul(
                                                                        temr,
                                                                        cadd(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                        ),
                                                                    );
                                                                }
                                                                ip += 1;
                                                            }
                                                        }
                                                    } else {
                                                        if (int_type == 2i32) {
                                                            {
                                                                ip = 0i32;
                                                                while (ip < (block_size - 1i32)) {
                                                                    {
                                                                        jj = (ip + 1i32);
                                                                        t1 = csub(
                                                                            (*checked_at(
                                                                                texp,
                                                                                (jj) as isize,
                                                                            )),
                                                                            (*checked_at(
                                                                                texp,
                                                                                (ip) as isize,
                                                                            )),
                                                                        );
                                                                        func1 = csub(
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (jj) as isize,
                                                                                )),
                                                                            ),
                                                                            cmul(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (ip) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    texp,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                        );
                                                                        func2 = cmul(
                                                                            csub(
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (jj) as isize,
                                                                                )),
                                                                                (*checked_at(
                                                                                    up,
                                                                                    (ip) as isize,
                                                                                )),
                                                                            ),
                                                                            t1,
                                                                        );
                                                                        (*checked_at(
                                                                            vos1,
                                                                            (ip) as isize,
                                                                        )) = cmul(
                                                                            cdiv(cdp, sigh),
                                                                            cmul(
                                                                                csub(
                                                                                    func1,
                                                                                    cdiv(
                                                                                        func2, sigh,
                                                                                    ),
                                                                                ),
                                                                                temr,
                                                                            ),
                                                                        );
                                                                    }
                                                                    ip += 1;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }

                                                /**********************************************************
                                                 *		update output array			*
                                                 **********************************************************/
                                                {
                                                    ip = 0i32;
                                                    while (ip < (block_size - 1i32)) {
                                                        {
                                                            (*checked_at(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        response1,
                                                                        (iw) as isize,
                                                                    )),
                                                                    (ix) as isize,
                                                                )),
                                                                (iz) as isize,
                                                            )) = cadd(
                                                                (*checked_at(
                                                                    (*checked_at(
                                                                        (*checked_at(
                                                                            response1,
                                                                            (iw) as isize,
                                                                        )),
                                                                        (ix) as isize,
                                                                    )),
                                                                    (iz) as isize,
                                                                )),
                                                                (*checked_at(vos1, (ip) as isize)),
                                                            );
                                                        }
                                                        ip += 1;
                                                    }
                                                }
                                            }
                                            ix += 1;
                                        }
                                    }
                                }
                                iz += 1;
                            }
                        }
                        /* update loop variables */
                        nblock -= 1;
                        if (nblock != 0i32) {
                            /*	bp=p[block_size].r; */
                            bp = (*checked_at(p, (block_size - 1i32) as isize)).r;
                            left += 1;
                            if (left > block_size) {
                                left = ((left as i32) - (block_size as i32)) as i32;
                                nblock += 1;
                            }
                        } else {
                            if (int_type == 1i32) {
                                bp = bp1;
                            } else {
                                if (int_type == 2i32) {
                                    bp = (0f64 as f32);
                                }
                            }

                            break;
                        }
                    }
                }
                iw += 1;
            }
        }
        /* last block, exit while loop */
        /* free allocated space */
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        /* free working space */
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
        ();
    }
}

// Refl/lib/synthetic.c: compute_synthetics
pub(super) unsafe fn raw_compute_synthetics(
    mut verbose: i32,
    mut nt: i32,
    mut ntc: i32,
    mut nx: i32,
    mut nor: i32,
    mut nw: i32,
    mut nlayers: i32,
    mut lsource: i32,
    mut nf: i32,
    mut flt: i32,
    mut vsp: i32,
    mut win: i32,
    mut wtype: i32,
    mut tlag: f32,
    mut red_vel: f32,
    mut decay: f32,
    mut tsec: f32,
    mut bx: f32,
    mut dx: f32,
    mut w1: f32,
    mut w2: f32,
    mut fil_phase: *mut i32,
    mut nfilters: i32,
    mut fil_type: *mut i32,
    mut dbpo: *mut f32,
    mut f1: *mut f32,
    mut f2: *mut f32,
    mut lobs: *mut i32,
    mut cl: *mut f32,
    mut t: *mut f32,
    mut response: *mut *mut *mut Complex,
    mut reflectivity: *mut *mut f32,
    mut outfp: *mut CFile,
) {
    let mut arena = Arena::default();
    {
        let mut il: i32 = 0;
        let mut iw: i32 = 0;
        let mut ix: i32 = 0;
        let mut iz: i32 = 0i32;
        /* loop counters */
        let mut kz: i32 = 0;
        /* auxiliary index */
        let mut ntfft: i32 = 0;
        /* number of FFT samples */
        let mut nw1: i32 = 0;
        /* number of frequencies for FFT */
        let mut sum: f32 = (0f64 as f32);
        /* auxiliary variable */
        let mut zr: *mut f32 = std::ptr::null_mut();
        let mut zs: f32 = (0f64 as f32);
        /* depths of source and receivers */
        let mut sphrd: f32 = (1f64 as f32);
        /* factor for flattening earth correction */
        let mut if1: i32 = 0i32;
        let mut if2: i32 = 0i32;
        /* freq indices to start and end Hanning window */
        let mut unexp: f32 = 0.0;
        /* parameter to undo decay factor */
        /* float x,z;*/
        /* range and depth */
        let mut x: f32 = 0.0;
        /* range  */
        let mut w: f32 = 0.0;
        /* frequency in rad/sec */
        let mut wpie: Complex = cmplx(0.0, 0.0);
        /* complex frequency */
        let mut wsq: Complex = cmplx(0.0, 0.0);
        /* complex frequency squared */
        let mut rvfac: Complex = cmplx(0.0, 0.0);
        /* reducing velocity factor */
        let mut refw: *mut Complex = std::ptr::null_mut();
        /* scratch array for complex wavefield */
        let mut window: Complex = cmplx(0.0, 0.0);
        /* frequency domain Hanning window */
        /* if requested, check input parameters */
        if ((verbose == 1i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "in compute reflectivity series, checking input ...\n"
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "nt=%d nx=%d nor=%d nw=%d",
                nt,
                nx,
                nor,
                nw
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                " nlayers=%d lsource=%d nf=%d\n",
                nlayers,
                lsource,
                nf
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "flt=%d vsp=%d win=%d wtype=%d\n",
                flt,
                vsp,
                win,
                wtype
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "tlag=%g red_vel=%g decay=%g ",
                (tlag as f64),
                (red_vel as f64),
                (decay as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "tsec=%g bx=%g dx=%g w1=%g w2=%g\n",
                (tsec as f64),
                (bx as f64),
                (dx as f64),
                (w1 as f64),
                (w2 as f64)
            );
            oxwave_core::fprintf!(&oxwave_core::prelude::stderr(), "nfilters=%d\n", nfilters);
        }

        if ((verbose == 2i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(
                &*outfp,
                "in compute reflectivity series, checking input ...\n"
            );
            oxwave_core::fprintf!(&*outfp, "nt=%d nx=%d nor=%d nw=%d", nt, nx, nor, nw);
            oxwave_core::fprintf!(
                &*outfp,
                " nlayers=%d lsource=%d nf=%d\n",
                nlayers,
                lsource,
                nf
            );
            oxwave_core::fprintf!(
                &*outfp,
                "flt=%d vsp=%d win=%d wtype=%d\n",
                flt,
                vsp,
                win,
                wtype
            );
            oxwave_core::fprintf!(
                &*outfp,
                "tlag=%g red_vel=%g decay=%g ",
                (tlag as f64),
                (red_vel as f64),
                (decay as f64)
            );
            oxwave_core::fprintf!(
                &*outfp,
                "tsec=%g bx=%g dx=%g w1=%g w2=%g\n",
                (tsec as f64),
                (bx as f64),
                (dx as f64),
                (w1 as f64),
                (w2 as f64)
            );
            oxwave_core::fprintf!(&*outfp, "nfilters=%d\n", nfilters);
        }

        /* compute number of fft samples and frequencies */
        ntfft = oxwave_pfa_fft::npfa(ntc);
        nw1 = (ntfft / 2i32);
        if ((nf > 0i32) && (nf <= nw)) {
            nw = nf;
        }

        oxwave_core::fprintf!(
            &oxwave_core::prelude::stderr(),
            "\n\n nt=%d ntfft=%d nw=%d nw1=%d\n",
            nt,
            ntfft,
            nw,
            nw1
        );
        /* allocate working space */
        zr = arena.alloc::<f32>((nor as usize) as usize);
        refw = arena.alloc::<Complex>((ntfft as usize) as usize);
        /* if not provided, set reducing velocity to highest pwave velocity */
        if ((red_vel as f64) < 0f64) {
            red_vel = (1.51f64 as f32);
            {
                il = 0i32;
                while (il < nlayers) {
                    {
                        if ((*checked_at(cl, (il) as isize)) > red_vel) {
                            red_vel = (*checked_at(cl, (il) as isize));
                        }
                    }
                    il += 1;
                }
            }
        }

        /* compute depths of source and receivers */
        iz = 0i32;
        {
            il = 0i32;
            while (il < (nlayers - 1i32)) {
                {
                    kz = (il + 1i32);
                    sum = ((sum as f32) + ((*checked_at(t, (il) as isize)) as f32)) as f32;
                    if (kz == (lsource - 1i32)) {
                        zs = sum;
                    }

                    /* depth of the source */
                    if ((iz < nor) && (kz == ((*checked_at(lobs, (iz) as isize)) - 1i32))) {
                        (*checked_at(zr, (iz) as isize)) = sum;
                        /* depth of the receivers */
                        iz += 1;
                    }
                }
                il += 1;
            }
        }
        /* if requested, output processing information */
        if ((verbose == 1i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "\nSOURCE AND RECEIVER INFORMATION\n"
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "\nReducing velocity=%g source depth=%g\n",
                (red_vel as f64),
                (zs as f64)
            );
            oxwave_core::fprintf!(
                &oxwave_core::prelude::stderr(),
                "\nReceiver Number  Depth   \n"
            );
            {
                iz = 0i32;
                while (iz < nor) {
                    {
                        oxwave_core::fprintf!(
                            &oxwave_core::prelude::stderr(),
                            "%8d%15.3f\n",
                            iz,
                            ((*checked_at(zr, (iz) as isize)) as f64)
                        );
                    }
                    iz += 1;
                }
            }
        }

        if ((verbose == 2i32) || (verbose == 3i32)) {
            oxwave_core::fprintf!(&*outfp, "\nSOURCE AND RECEIVER INFORMATION\n");
            oxwave_core::fprintf!(
                &*outfp,
                "\nReducing velocity=%g source depth=%g\n",
                (red_vel as f64),
                (zs as f64)
            );
            oxwave_core::fprintf!(&*outfp, "\nReceiver Number  Depth   \n");
            {
                iz = 0i32;
                while (iz < nor) {
                    {
                        oxwave_core::fprintf!(
                            &*outfp,
                            "%8d%15.3f\n",
                            iz,
                            ((*checked_at(zr, (iz) as isize)) as f64)
                        );
                    }
                    iz += 1;
                }
            }
        }

        /* compute windowed frequencies */
        if (win == 1i32) {
            if (w1 == (0i32 as f32)) {
                if1 = ((0.14999999999999999f64 * (nf as f64)) as i32);
            } else {
                if1 = ((w1 * tsec) as i32);
            }

            if (w2 == (0i32 as f32)) {
                if2 = ((0.84999999999999998f64 * (nf as f64)) as i32);
            } else {
                if2 = ((w2 * tsec) as i32);
            }
        }

        unexp = (decay * tsec);
        /* to undo exponential decay */
        /* main loop to compute a vsp or a seismogram */
        if (vsp == 1i32) {
            {
                ix = 0i32;
                while (ix < nx) {
                    {
                        /* compute range and flattening correction factor */
                        x = (bx + (dx * (ix as f32)));
                        if (flt == 1i32) {
                            sphrd = ((((x as f64) / 6371f64)
                                / (if ((((x as f64) / 6371f64).sin()) < (0i32 as f64)) {
                                    (-(((x as f64) / 6371f64).sin()))
                                } else {
                                    (((x as f64) / 6371f64).sin())
                                }))
                            .sqrt() as f32);
                        }

                        /* loop over the receivers */
                        {
                            iz = 0i32;
                            while (iz < nor) {
                                {
                                    /* z=zr[iz]; */
                                    /* loop over frequencies */
                                    {
                                        iw = 0i32;
                                        while (iw < nw1) {
                                            {
                                                /* initialize output array */
                                                if (iw < nw) {
                                                    w = ((((((iw + 1i32) as f64)
                                                        * (3.1415926535897931f64))
                                                        * (2i32 as f64))
                                                        / (tsec as f64))
                                                        as f32);
                                                    (*checked_at(refw, (iw + 1i32) as isize)) =
                                                        (*checked_at(
                                                            (*checked_at(
                                                                (*checked_at(
                                                                    response,
                                                                    (iw) as isize,
                                                                )),
                                                                (ix) as isize,
                                                            )),
                                                            (iz) as isize,
                                                        ));
                                                } else {
                                                    w = (0f64 as f32);
                                                    (*checked_at(refw, (iw + 1i32) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                }

                                                /* compute complex frequencies */
                                                wpie = cmplx(w, decay);
                                                wsq = cneg(cmul(
                                                    cmplx((0f64 as f32), (1f64 as f32)),
                                                    wpie,
                                                ));
                                                /* compute reducing velocity factor */
                                                raw_red_vel_factor(
                                                    x,
                                                    red_vel,
                                                    tlag,
                                                    wpie,
                                                    wsq,
                                                    (&raw mut rvfac),
                                                );
                                                /* compute Hanning window */
                                                raw_compute_Hanning_window(
                                                    win,
                                                    iw,
                                                    if1,
                                                    if2,
                                                    nw,
                                                    (&raw mut window),
                                                );
                                                /* apply Hanning window */
                                                (*checked_at(refw, (iw + 1i32) as isize)) = cmul(
                                                    (*checked_at(refw, (iw + 1i32) as isize)),
                                                    cmul(rvfac, window),
                                                );
                                            }
                                            iw += 1;
                                        }
                                    }
                                    /* construct time image */
                                    raw_construct_tx_trace(
                                        nw1,
                                        nt,
                                        ntc,
                                        ntfft,
                                        nfilters,
                                        fil_phase,
                                        tsec,
                                        unexp,
                                        sphrd,
                                        refw,
                                        fil_type,
                                        f1,
                                        f2,
                                        dbpo,
                                        (*checked_at(reflectivity, (ix) as isize)),
                                    );
                                }
                                iz += 1;
                            }
                        }
                    }
                    ix += 1;
                }
            }
        } else {
            /* compute "normal" seismogram */
            /* loop over recievers */
            {
                iz = 0i32;
                while (iz < nor) {
                    {
                        /* z=zr[iz]; */
                        /* loop over ranges (seismogram traces) */
                        {
                            ix = 0i32;
                            while (ix < nx) {
                                {
                                    x = (bx + (dx * (ix as f32)));
                                    /* if earth flattening correction requested */
                                    if (flt == 1i32) {
                                        sphrd = ((((x as f64) / 6371f64)
                                            / (if ((((x as f64) / 6371f64).sin()) < (0i32 as f64)) {
                                                (-(((x as f64) / 6371f64).sin()))
                                            } else {
                                                (((x as f64) / 6371f64).sin())
                                            }))
                                        .sqrt()
                                            as f32);
                                    }

                                    /* loop over frequencies */
                                    {
                                        iw = 0i32;
                                        while (iw < nw1) {
                                            {
                                                /* initialize output array */
                                                if (iw < nw) {
                                                    w = ((((((iw + 1i32) as f64)
                                                        * (3.1415926535897931f64))
                                                        * (2i32 as f64))
                                                        / (tsec as f64))
                                                        as f32);
                                                    (*checked_at(refw, (iw + 1i32) as isize)) =
                                                        (*checked_at(
                                                            (*checked_at(
                                                                (*checked_at(
                                                                    response,
                                                                    (iw) as isize,
                                                                )),
                                                                (ix) as isize,
                                                            )),
                                                            (iz) as isize,
                                                        ));
                                                } else {
                                                    w = (0f64 as f32);
                                                    (*checked_at(refw, (iw + 1i32) as isize)) =
                                                        cmplx((0f64 as f32), (0f64 as f32));
                                                }

                                                /* compute complex frequencies */
                                                wpie = cmplx(w, decay);
                                                wsq =
                                                    cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie);
                                                wsq = cipow(wsq, wtype);
                                                /* compute reducing velocity factor */
                                                raw_red_vel_factor(
                                                    x,
                                                    red_vel,
                                                    tlag,
                                                    wpie,
                                                    wsq,
                                                    (&raw mut rvfac),
                                                );
                                                /* compute Hanning window */
                                                raw_compute_Hanning_window(
                                                    win,
                                                    iw,
                                                    if1,
                                                    if2,
                                                    nw,
                                                    (&raw mut window),
                                                );
                                                /* apply Hanning window and rvfac */
                                                (*checked_at(refw, (iw + 1i32) as isize)) = cmul(
                                                    (*checked_at(refw, (iw + 1i32) as isize)),
                                                    cmul(rvfac, window),
                                                );
                                            }
                                            iw += 1;
                                        }
                                    }
                                    /* construct time image */
                                    raw_construct_tx_trace(
                                        nw1,
                                        nt,
                                        ntc,
                                        ntfft,
                                        nfilters,
                                        fil_phase,
                                        tsec,
                                        unexp,
                                        sphrd,
                                        refw,
                                        fil_type,
                                        f1,
                                        f2,
                                        dbpo,
                                        (*checked_at(reflectivity, (ix) as isize)),
                                    );
                                }
                                ix += 1;
                            }
                        }
                    }
                    iz += 1;
                }
            }
        }

        /* free working space */
        ();
        ();
    }
}

// Refl/lib/synthetic.c: red_vel_factor
pub(super) unsafe fn raw_red_vel_factor(
    mut x: f32,
    mut red_vel: f32,
    mut tlag: f32,
    mut wpie: Complex,
    mut wsq: Complex,
    mut rvfac: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        /* compute reducing velocity factor */
        if ((red_vel as f64) == 0f64) {
            (*rvfac) = crmul(cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie), tlag);
        } else {
            (*rvfac) = crmul(
                cmul(cmplx((0f64 as f32), (1f64 as f32)), wpie),
                (tlag - (x / red_vel)),
            );
        }

        /* multiply by complex factor */
        (*rvfac) = cmul(cwp_cexp1((*rvfac)), wsq);
    }
}

// Refl/lib/synthetic.c: construct_tx_trace
pub(super) unsafe fn raw_construct_tx_trace(
    mut nw: i32,
    mut nt: i32,
    mut ntc: i32,
    mut ntfft: i32,
    mut nfilters: i32,
    mut filters_phase: *mut i32,
    mut tsec: f32,
    mut unexp: f32,
    mut sphrd: f32,
    mut refw: *mut Complex,
    mut fil_type: *mut i32,
    mut f1: *mut f32,
    mut f2: *mut f32,
    mut dbpo: *mut f32,
    mut reft: *mut f32,
) {
    let mut arena = Arena::default();
    {
        let mut it: i32 = 0;
        /* loop counter */
        let mut a: f32 = 0.0;
        let mut b: f32 = 0.0;
        let mut c: f32 = 0.0;
        /* auxiliary variables */
        /* if requested, apply filters */
        if (nfilters != 0i32) {
            raw_apply_filters(
                filters_phase,
                nfilters,
                nw,
                tsec,
                f1,
                f2,
                dbpo,
                fil_type,
                refw,
            );
        }

        /* inverse fft to transfer data to t-x domain */
        oxwave_pfa_fft::pfacc(
            (-1i32),
            std::slice::from_raw_parts_mut(refw, (ntfft) as usize),
        );
        /* construct the t-x image */
        {
            it = 0i32;
            while (it < nt) {
                {
                    a = ((it as f32) / (ntc as f32));
                    b = (unexp * a);
                    c = ((b as f64).exp() as f32);
                    (*checked_at(reft, (it) as isize)) =
                        (((*checked_at(refw, (it) as isize)).r * c) * sphrd);
                }
                it += 1;
            }
        }
        /* reinitialize working arrays */
        {
            it = 0i32;
            while (it < ntfft) {
                {
                    (*checked_at(refw, (it) as isize)) = cmplx((0f64 as f32), (0f64 as f32));
                }
                it += 1;
            }
        }
    }
}

// Refl/lib/synthetic.c: compute_Hanning_window
pub(super) unsafe fn raw_compute_Hanning_window(
    mut iwin: i32,
    mut iw: i32,
    mut if1: i32,
    mut if2: i32,
    mut nf: i32,
    mut win: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        /* compute Hanning window */
        if (iwin == 1i32) {
            if (iw < (if1 - 1i32)) {
                (*win) = crmul(
                    cmplx((1f64 as f32), (0f64 as f32)),
                    ((0.5f64
                        * (1f64
                            + (((3.1415926535897931f64) * (((if1 - iw) + 1i32) as f64))
                                / (if1 as f64))
                                .cos())) as f32),
                );
            } else {
                if ((iw >= (if1 - 1i32)) && (iw <= (if2 - 1i32))) {
                    (*win) = cmplx((1f64 as f32), (0f64 as f32));
                } else {
                    if ((iw > (if2 - 1i32)) && (iw < (nf - 1i32))) {
                        (*win) = crmul(
                            cmplx((1f64 as f32), (0f64 as f32)),
                            ((0.5f64
                                * (1f64
                                    + (((3.1415926535897931f64) * (((iw - if2) + 1i32) as f64))
                                        / ((nf - if2) as f64))
                                        .cos())) as f32),
                        );
                    } else {
                        (*win) = cmplx((0f64 as f32), (0f64 as f32));
                    }
                }
            }
        } else {
            if (iw < (nf - 1i32)) {
                (*win) = crmul(
                    cmplx((1f64 as f32), (0f64 as f32)),
                    ((0.5f64
                        * (1f64
                            + (((3.1415926535897931f64) * ((iw + 1i32) as f64)) / (nf as f64))
                                .cos())) as f32),
                );
            } else {
                (*win) = cmplx((0f64 as f32), (0f64 as f32));
            }
        }
    }
}

// Refl/lib/synthetic.c: apply_filters
pub(super) unsafe fn raw_apply_filters(
    mut min_phase: *mut i32,
    mut nfilters: i32,
    mut nw: i32,
    mut tsec: f32,
    mut f1: *mut f32,
    mut f2: *mut f32,
    mut dbpo: *mut f32,
    mut filter_type: *mut i32,
    mut refw: *mut Complex,
) {
    let mut arena = Arena::default();
    {
        let mut i: i32 = 0;
        let mut j: i32 = 0;
        let mut iw: i32 = 0;
        /* loop counters */
        let mut nw1: i32 = (nw + 1i32);
        /* number of frequencies (positive+zero) */
        let mut nwt: i32 = (2i32 * nw);
        /* total number of frequencies */
        let mut index1: i32 = 0;
        let mut index2: i32 = 0;
        /* auxiliary indices for notch filter */
        let mut delfq: f32 = (0f64 as f32);
        let mut poles: f32 = 0.0;
        /* filter poles */
        let mut r: f32 = 0.0;
        let mut factor: f32 = 0.0;
        /* auxiliary variables */
        let mut maxabs: f32 = 0.0;
        let mut filter: *mut Complex = std::ptr::null_mut();
        /* scratch array for filter coefficients */
        let mut filter1: *mut Complex = std::ptr::null_mut();
        /* scratch array for filter coefficients */
        /* defensive programming */
        if (nfilters == 0i32) {
            return;
        }

        /* allocate working space */
        filter = arena.alloc::<Complex>(((2i32 * nwt) as usize) as usize);
        filter1 = arena.alloc::<Complex>(((2i32 * nwt) as usize) as usize);
        /* initialize filter array */
        {
            iw = 0i32;
            while (iw < nw1) {
                (*checked_at(filter, (iw) as isize)) = cmplx((1f64 as f32), (0f64 as f32));
                iw += 1;
            }
        }
        /* main loop over number of requested filters */
        {
            i = 0i32;
            while (i < nfilters) {
                {
                    /* smearing range */
                    if ((*checked_at(min_phase, (i) as isize)) == 1i32) {
                        delfq =
                            ((((2i32 as f64) * (3.1415926535897931f64)) / (tsec as f64)) as f32);
                    }

                    if ((*checked_at(filter_type, (i) as isize)) == 1i32) {
                        /* if low-pass Butterworth */
                        poles = (((((*checked_at(dbpo, (i) as isize)) as f64) / 6f64)
                            + (1i32 as f64)) as f32);
                        r = ((1f64 / ((((*checked_at(f1, (i) as isize)) - delfq) * tsec) as f64))
                            as f32);
                        /* compute filter coefficients */
                        {
                            iw = 0i32;
                            while (iw < nw1) {
                                {
                                    factor = ((1f64
                                        / (1f64 + (((iw as f32) * r) as f64).powf((poles as f64)))
                                            .sqrt())
                                        as f32);
                                    (*checked_at(filter, (iw) as isize)) =
                                        crmul((*checked_at(filter, (iw) as isize)), factor);
                                }
                                iw += 1;
                            }
                        }
                    } else {
                        if ((*checked_at(filter_type, (i) as isize)) == 2i32) {
                            /* hi-pass Butterworth */
                            poles = (((((*checked_at(dbpo, (i) as isize)) as f64) / 6f64)
                                + (1i32 as f64)) as f32);
                            r = ((1f64
                                / ((((*checked_at(f1, (i) as isize)) + delfq) * tsec) as f64))
                                as f32);
                            /* compute filter coefficients */
                            {
                                iw = 0i32;
                                while (iw < nw1) {
                                    {
                                        factor = ((((iw as f32) * r) as f64).powf((poles as f64))
                                            as f32);
                                        (*checked_at(filter, (iw) as isize)) = crmul(
                                            (*checked_at(filter, (iw) as isize)),
                                            (((factor / ((1i32 as f32) + factor)) as f64).sqrt()
                                                as f32),
                                        );
                                    }
                                    iw += 1;
                                }
                            }
                        } else {
                            if ((*checked_at(filter_type, (i) as isize)) == 3i32) {
                                /* notch filter */
                                index1 = (((((*checked_at(f1, (i) as isize)) - delfq) * tsec)
                                    + (1i32 as f32))
                                    as i32);
                                index2 = (((((*checked_at(f2, (i) as isize)) + delfq) * tsec)
                                    + (2i32 as f32))
                                    as i32);
                                if ((index2 - index1) < 16i32) {
                                    oxwave_core::warn!("in notch filter %d f1 too close to f2, this filter skiped\n",i);
                                    i += 1;
                                    continue;
                                }

                                /* db reduction Neils design */
                                r = ((1f64
                                    - (10f64).powf(
                                        (((-(*checked_at(dbpo, (i) as isize))) as f64) / 20f64),
                                    )) as f32);
                                {
                                    j = index1;
                                    while (j <= index2) {
                                        {
                                            factor = ((((((2i32 as f64)
                                                * (3.1415926535897931f64))
                                                * ((j - index1) as f64))
                                                / ((index2 - index1) as f64))
                                                .cos()
                                                - (1i32 as f64))
                                                as f32);
                                            (*checked_at(filter, (j - 1i32) as isize)) = crmul(
                                                (*checked_at(filter, (j - 1i32) as isize)),
                                                (((((1i32 as f32) + (r * factor)) as f64) / 2f64)
                                                    as f32),
                                            );
                                        }
                                        j += 1;
                                    }
                                }
                            } else {
                                oxwave_core::warn!("unknown filter type for filter number %d, this filter skiped\n",i);
                            }
                        }
                    }

                    /* apply minimum phase correction */
                    if ((*checked_at(min_phase, (i) as isize)) == 1i32) {
                        {
                            iw = 0i32;
                            while (iw < nw1) {
                                {
                                    if (((*checked_at(filter, (iw) as isize)).r as f64) == 0f64) {
                                        (*checked_at(filter, (iw) as isize)) =
                                            cmplx(((-30f64) as f32), (0f64 as f32));
                                    } else {
                                        /* set A[w]=ln(f[w]) */
                                        (*checked_at(filter, (iw) as isize)) = cmplx(
                                            (((*checked_at(filter, (iw) as isize)).r as f64).ln()
                                                as f32),
                                            (0f64 as f32),
                                        );
                                    }
                                }
                                iw += 1;
                            }
                        }
                        /* take care of negative frequencies needed for FFT */
                        {
                            iw = nw1;
                            while (iw < nwt) {
                                {
                                    (*checked_at(filter, (iw) as isize)) =
                                        conjg((*checked_at(filter, ((nwt - iw) + 1i32) as isize)));
                                }
                                iw += 1;
                            }
                        }
                        /* take inverse FFT to go to time domain */
                        oxwave_pfa_fft::pfacc(
                            (-1i32),
                            std::slice::from_raw_parts_mut(filter, (nwt) as usize),
                        );
                        {
                            iw = 1i32;
                            while (iw < nw1) {
                                {
                                    r = (((1f64
                                        + (((3.1415926535897931f64) * (iw as f64)) / (nw as f64))
                                            .cos())
                                        / 2f64) as f32);
                                    (*checked_at(filter, (iw) as isize)) =
                                        crmul((*checked_at(filter, (iw) as isize)), r);
                                    (*checked_at(filter, ((nwt - iw) + 1i32) as isize)) = crmul(
                                        (*checked_at(filter, ((nwt - iw) + 1i32) as isize)),
                                        r,
                                    );
                                }
                                iw += 1;
                            }
                        }
                        /* construct filter coefficients for all frequencies */
                        {
                            iw = 0i32;
                            while (iw < nw1) {
                                {
                                    (*checked_at(filter1, (iw) as isize)) =
                                        (*checked_at(filter, (iw) as isize));
                                }
                                iw += 1;
                            }
                        }
                        {
                            iw = nw1;
                            while (iw < nwt) {
                                {
                                    (*checked_at(filter1, (iw) as isize)) =
                                        cneg(conjg((*checked_at(filter, (iw) as isize))));
                                }
                                iw += 1;
                            }
                        }
                        /* take forward FFT to get back to frequency domain */
                        oxwave_pfa_fft::pfacc(
                            1i32,
                            std::slice::from_raw_parts_mut(filter, (nwt) as usize),
                        );
                        oxwave_pfa_fft::pfacc(
                            1i32,
                            std::slice::from_raw_parts_mut(filter1, (nwt) as usize),
                        );
                        maxabs = (0f64 as f32);
                        {
                            iw = 0i32;
                            while (iw < nw1) {
                                {
                                    (*checked_at(filter, (iw) as isize)) = cwp_cexp1(cadd(
                                        (*checked_at(filter, (iw) as isize)),
                                        (*checked_at(filter1, (iw) as isize)),
                                    ));
                                    maxabs = (if ((rcabs((*checked_at(filter, (iw) as isize))))
                                        > (maxabs))
                                    {
                                        (rcabs((*checked_at(filter, (iw) as isize))))
                                    } else {
                                        (maxabs)
                                    });
                                }
                                iw += 1;
                            }
                        }
                        /* normalize */
                        maxabs = ((1f64 / (maxabs as f64)) as f32);
                        {
                            iw = 0i32;
                            while (iw < nw1) {
                                {
                                    (*checked_at(filter, (iw) as isize)) =
                                        crmul((*checked_at(filter, (iw) as isize)), maxabs);
                                }
                                iw += 1;
                            }
                        }
                    }
                }
                i += 1;
            }
        }
        /* apply the filter */
        {
            iw = 0i32;
            while (iw < nw1) {
                {
                    (*checked_at(refw, (iw) as isize)) = cmul(
                        (*checked_at(refw, (iw) as isize)),
                        (*checked_at(filter, (iw) as isize)),
                    );
                }
                iw += 1;
            }
        }
        {
            iw = nw1;
            while (iw < nwt) {
                {
                    (*checked_at(refw, (iw) as isize)) =
                        conjg((*checked_at(refw, (nwt - iw) as isize)));
                }
                iw += 1;
            }
        }
        /* free allocated space */
        ();
        ();
    }
}
