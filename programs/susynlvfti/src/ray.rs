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

//! Anisotropic ray tracer from su/main/synthetics_waveforms_testpatterns/susynlvfti.c.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Dave Hale, 09/17/91, Colorado School of Mines.
//! Upgrade to Transversely Isotropic medium: T. Alkhallfah, 4/15/93.

// Keep the C intermediate types, including float operations inside double
// expressions and float rounding when rotating coordinates. Uninitialized C
// scratch variables are initialized to zero so Rust has defined behavior.
#![allow(unused_parens, unused_assignments, unused_mut, unused_variables)]

/// Trace a ray for v(x,z) = v00 + dvdx*x + dvdz*z in a factorized TI medium.
/// Alkhalifah, T., 1993, Efficient synthetic seismograms for transversely
/// isotropic media with constant velocity gradient, CWP project review.
#[allow(clippy::too_many_arguments)]
pub(super) fn raylvt(
    mut v00: f32,
    mut dvdx: f32,
    mut dvdz: f32,
    mut x0: f32,
    mut z0: f32,
    mut x: f32,
    mut z: f32,
    mut a: f32,
    mut f: f32,
    mut l: f32,
    mut ntries: i32,
    mut nitmax: f32,
    mut epst: f32,
    mut epsx: f32,
    mut angxs: f32,
    mut c: &mut f32,
    mut s: &mut f32,
    mut t: &mut f32,
    mut q: &mut f32,
)
/*****************************************************************************
Trace ray between two points, for linear velocity v(x,z) = v00+dvdx*x+dvdz*z
in a transversely isotropic medium
Reference:
Alkhalifah, T., 1993, Efficient synthetic seismograms for transversely isotropic media with constant velocity gradient:CWP project review.
******************************************************************************
Input:
v00	     velocity at (x=0,z=0)
dvdx	    derivative dv/dx of velocity with respect to x
dvdz	    derivative dv/dz of velocity with respect to z
x0		x coordinate at beginning of ray
z0		z coordinate at beginning of ray
x		x coordinate at end of ray
z		z coordinate at end of ray
ntries	  number of iterations in Snell's law search
epsx	    lateral offset tolerance
a,f,l	   ratios of velocities for the FAI transversely isotropic medium
angxs		angle of symmetry axis of anisotriopy

Output:
c		cosine of propagation angle at end of ray
s		sine of propagation angle at end of ray
t		time along raypath
q		integral of velocity along raypath
*****************************************************************************/
{
    let mut px0: f64 = 0.0;
    let mut v0: f64 = 0.0;
    let mut aa: f64 = 0.0;
    let mut alfa: f64 = 0.0;
    let mut beta: f64 = 0.0f64;
    let mut vcp: f64 = 0.0f64;
    let mut oa: f64 = 0.0;
    let mut v: f64 = 0.0;
    let mut c1: f64 = 0.0;
    let mut c2: f64 = 0.0;
    let mut f1: f64 = 0.0;
    let mut fx: f64 = 0.0;
    let mut cr: f64 = 0.0;
    let mut sr: f64 = 0.0;
    let mut r: f64 = 0.0;
    let mut or: f64 = 0.0;
    let mut vcpi: f64 = 0.0f64;
    let mut vcps: f64 = 0.0f64;
    let mut p11: f64 = 0.0;
    let mut so: f64 = 0.0;
    let mut co: f64 = 0.0f64;
    let mut ci: f64 = 0.0f64;
    let mut cz2: f64 = 0.0;
    let mut dso: f64 = 0.0;
    let mut so1: f64 = 0.0f64;
    let mut ss: f64 = 0.0f64;
    let mut cs: f64 = 0.0f64;
    let mut si2: f64 = 0.0;
    let mut co2: f64 = 0.0;
    let mut sc: f64 = 0.0;
    let mut sx: f64 = 0.0f64;
    let mut cx: f64 = 0.0f64;
    let mut so3: f64 = 0.0;
    let mut diff: f64 = 0.0;
    let mut so4: f64 = 0.0;
    let mut diff1: f64 = 0.0f64;
    let mut sz2: f64 = 0.0;
    let mut w1: f64 = 0.0;
    let mut w2: f64 = 0.0;
    let mut w3: f64 = 0.0;
    let mut px2: f64 = 0.0;
    let mut pz2: f64 = 0.0;
    let mut betaa: f64 = 0.0;
    let mut betab: f64 = 0.0;
    let mut alpha: f64 = 0.0;
    let mut pz0: f64 = 0.0;
    let mut betat: f64 = 0.0;
    let mut xx: f64 = 0.0;
    let mut dvcp: f64 = 0.0;
    let mut bbb: f64 = 0.0;
    let mut soo: f64 = 0.0;
    let mut umax: f64 = 0.0;
    let mut time: f64 = 0.0f64;
    let mut fctn1: f64 = 0.0;
    let mut fctn2: f64 = 0.0;
    let mut fctnu: f64 = 0.0;
    let mut temp: f64 = 0.0;
    let mut mx: f64 = 0.0;
    let mut sm: f64 = 0.0f64;
    let mut vg: f64 = 0.0;
    let mut gang: f64 = 0.0f64;
    let mut rt: f64 = 0.0;
    let mut timenew: f64 = 0.0f64;
    let mut tnm: f64 = 0.0;
    let mut u: f64 = 0.0;
    let mut sum: f64 = 0.0;
    let mut del: f64 = 0.0;
    let mut si: f64 = 0.0f64;
    let mut vrat: f64 = 0.0;
    let mut zz: f64 = 0.0;
    let mut uv: f64 = 0.0;
    let mut cw: f64 = 0.0f64;
    let mut sw: f64 = 0.0f64;
    let mut uv0: f64 = 0.0;
    let mut px: f64 = 0.0;
    let mut pz: f64 = 0.0;
    let mut j: i32 = 0;
    let mut count: i32 = 0i32;
    let mut it: i32 = 0i32;
    let mut nit: i32 = 0;
    let mut sig: i32 = 0;
    let mut zneg: i32 = 0;
    let mut err: i32 = 0;
    let mut tfturn: i32 = 0;
    let mut sign: i32 = 0;
    /* if (x,z) same as (x0,z0) */
    if ((x == x0) && (z == z0)) {
        (*c) = (1.0f64 as f32);
        (*s) = (0.0f64 as f32);
        (*t) = (0.0f64 as f32);
        (*q) = (0.0f64 as f32);
        return;
    }
    /*defining variables variables*/
    c1 = ((1.0f32 + l) as f64);
    c2 = ((a + l) as f64);
    w1 = ((a - l) as f64);
    w2 = ((1.0f32 - l) as f64);
    w3 = (((f + l) * (f + l)) as f64);
    /*if necessary, rotate coordinates until the symmetry axis is in z direction*/
    if ((angxs as f64) != 0.0f64) {
        sx = (angxs as f64).sin();
        cx = (angxs as f64).cos();
        co2 = (cx * cx);
        si2 = (sx * sx);
        alfa = ((c2 * si2) + (c1 * co2));
        beta = (((w1 * si2) - (w2 * co2)).powf(2.0f64) + (((4.0f64 * si2) * co2) * w3)).sqrt();
        vrat = (0.5f64 * (alfa + beta)).sqrt();
        temp = ((cx * (x0 as f64)) - (sx * (z0 as f64)));
        z0 = (((cx * (z0 as f64)) + (sx * (x0 as f64))) as f32);
        x0 = (temp as f32);
        temp = ((cx * (x as f64)) - (sx * (z as f64)));
        z = (((cx * (z as f64)) + (sx * (x as f64))) as f32);
        x = (temp as f32);
        temp = ((cx * (dvdx as f64)) - (sx * (dvdz as f64)));
        dvdz = ((((cx * (dvdz as f64)) + (sx * (dvdx as f64))) / vrat) as f32);
        dvdx = ((temp / vrat) as f32);
        v00 = (((v00 as f64) / vrat) as f32);
    }
    /* if constant velocity */
    if (((dvdx as f64) == 0.0f64) && ((dvdz as f64) == 0.0f64)) {
        x -= x0;
        z -= z0;
        r = (((x * x) + (z * z)) as f64).sqrt();
        or = (1.0f64 / r);
        (*s) = (((x as f64) * or) as f32);
        so = ((*s) as f64);
        dso = 0.1f64;
        so4 = 0.0f64;
        /*secant search for takeoff ray angle*/
        j = 0i32;
        while (j < ntries) {
            if (so > 1.0f64) {
                so = (0.5f64 * (1.0f64 + so1));
            }
            if (so < ((-1i32) as f64)) {
                so = (0.5f64 * (((-1i32) as f64) + so1));
            }
            so4 = so;
            co = (1.0f64 - (so * so)).sqrt();
            co2 = (co * co);
            si2 = (so * so);
            sc = (so * co);
            alfa = ((c2 * si2) + (c1 * co2));
            beta = (((w1 * si2) - (w2 * co2)).powf(2.0f64) + (((4.0f64 * si2) * co2) * w3)).sqrt();
            vcp = (0.5f64 * (alfa + beta)).sqrt();
            dvcp = (((((0.5f64 * ((a + l) as f64)) * sc)
                - ((0.5f64 * ((1.0f32 + l) as f64)) * sc))
                + ((0.25f64
                    * ((((2.0f64
                        * ((((a - l) as f64) * si2) - (((1.0f32 - l) as f64) * co2)))
                        * ((((a - l) as f64) * sc) + (((1.0f32 - l) as f64) * sc)))
                        + ((((4.0f64 * sc) * co2) * ((l + f) as f64)) * ((l + f) as f64)))
                        - ((((4.0f64 * si2) * sc) * ((l + f) as f64)) * ((l + f) as f64))))
                    / beta))
                / vcp);
            gang = (dvcp / vcp).atan();
            diff = (((z as f64) * (gang + (so).asin()).tan()) - (x as f64));
            if (j == 0i32) {
                if ((if (diff < 0.0f64) { (-diff) } else { diff }) < (epsx as f64)) {
                    err = 0i32;
                    break;
                }
                so1 = so;
                so = (so + dso);
                diff1 = diff;
            } else {
                count += 1i32;
                so3 = so;
                if ((if (diff < 0.0f64) { (-diff) } else { diff }) < (epsx as f64)) {
                    err = 0i32;
                    break;
                }
                so = (((diff1 * so3) - (diff * so1)) / (diff1 - diff));
                diff1 = diff;
                so1 = so3;
            }
            j += 1;
        }
        (*s) = ((((*s) as f64).asin() + gang).sin() as f32);
        (*c) = (((1.0f32 - ((*s) * (*s))) as f64).sqrt() as f32);
        vg = (vcp / (gang).cos());
        (*t) = ((r / ((v00 as f64) * vg)) as f32);
        (*q) = (((r * (v00 as f64)) * vg) as f32);
        /* if coordinates were rotated, unrotate cosine and sine */
        if ((angxs as f64) != 0.0f64) {
            temp = ((cx * ((*s) as f64)) + (sx * ((*c) as f64)));
            (*c) = (((cx * ((*c) as f64)) - (sx * ((*s) as f64))) as f32);
            (*s) = (temp as f32);
        }
        return;
    }
    /*define distances before axes rotation*/
    xx = ((x - x0) as f64);
    zz = ((z - z0) as f64);
    /*if not rotated*/
    cr = 1.0f64;
    sr = 0.0f64;
    /* if necessary, rotate coordinates so that v(x,z) = v0+a*(z-z0) */
    if (((dvdx as f64) != 0.0f64) || (dvdz < 0.0f32)) {
        aa = (((dvdx * dvdx) + (dvdz * dvdz)) as f64).sqrt();
        oa = (1.0f64 / aa);
        cr = ((dvdz as f64) * oa);
        sr = ((dvdx as f64) * oa);
        temp = ((cr * (x0 as f64)) - (sr * (z0 as f64)));
        z0 = (((cr * (z0 as f64)) + (sr * (x0 as f64))) as f32);
        x0 = (temp as f32);
        temp = ((cr * (x as f64)) - (sr * (z as f64)));
        z = (((cr * (z as f64)) + (sr * (x as f64))) as f32);
        x = (temp as f32);
    } else {
        aa = (dvdz as f64);
    }
    /*velocities and distances after possible rotation of coordinates*/
    v0 = ((v00 as f64) + (aa * (z0 as f64)));
    v = ((v00 as f64) + (aa * (z as f64)));
    x -= x0;
    z -= z0;
    z0 = ((v0 / aa) as f32);
    /* defining regions depending on wether we needed the equation reversed */
    if ((((x as f64) < 0.0f64) && (aa > 0.0f64)) || (((x as f64) > 0.0f64) && (aa < 0.0f64))) {
        sig = (-1i32);
    } else {
        sig = 1i32;
    }
    zneg = 1i32;
    /* if raypath is parallel to velocity gradient */
    if (((x * x) as f64) < ((1e-05f64 * (z as f64)) * (z as f64))) {
        so = sr;
        dso = 0.1f64;
        so4 = 0.0f64;
        /*secant search for takeoff ray angle*/
        j = 0i32;
        while (j < ntries) {
            if ((if (so < 0.0f64) { (-so) } else { so }) > 1.0f64) {
                so = (0.5f64 * (((sig as f64) * sm) + so4));
            }
            so4 = so;
            co = (1.0f64 - (so * so)).sqrt();
            co2 = (co * co);
            si2 = (so * so);
            sc = (so * co);
            alfa = ((c2 * si2) + (c1 * co2));
            beta = (((w1 * si2) - (w2 * co2)).powf(2.0f64) + (((4.0f64 * si2) * co2) * w3)).sqrt();
            vcp = (0.5f64 * (alfa + beta)).sqrt();
            dvcp = (((((0.5f64 * ((a + l) as f64)) * sc)
                - ((0.5f64 * ((1.0f32 + l) as f64)) * sc))
                + ((0.25f64
                    * ((((2.0f64
                        * ((((a - l) as f64) * si2) - (((1.0f32 - l) as f64) * co2)))
                        * ((((a - l) as f64) * sc) + (((1.0f32 - l) as f64) * sc)))
                        + ((((4.0f64 * sc) * co2) * ((l + f) as f64)) * ((l + f) as f64)))
                        - ((((4.0f64 * si2) * sc) * ((l + f) as f64)) * ((l + f) as f64))))
                    / beta))
                / vcp);
            gang = (dvcp / vcp).atan();
            diff = ((zz * (gang + (so).asin()).tan()) - xx);
            if (j == 0i32) {
                if ((if (diff < 0.0f64) { (-diff) } else { diff }) < (epsx as f64)) {
                    err = 0i32;
                    break;
                }
                so1 = so;
                so = (so + dso);
                diff1 = diff;
            } else {
                count += 1i32;
                so3 = so;
                if ((if (diff < 0.0f64) { (-diff) } else { diff }) < (epsx as f64)) {
                    err = 0i32;
                    break;
                }
                so = (((diff1 * so3) - (diff * so1)) / (diff1 - diff));
                diff1 = diff;
                so1 = so3;
            }
            j += 1;
        }
        /* if ray is going down */
        if ((z as f64) > 0.0f64) {
            (*s) = (((sr).asin() - gang).sin() as f32);
            (*c) = (((1.0f32 - ((*s) * (*s))) as f64).sqrt() as f32);
            vg = (vcp / (gang).cos());
            (*t) = (((v / v0).ln() / (aa * vg)) as f32);
            (*q) = ((((0.5f64 * (z as f64)) * (v + v0)) * vg) as f32);
        } else {
            (*s) = (((sr).asin() - gang).sin() as f32);
            (*c) = ((-((1.0f32 - ((*s) * (*s))) as f64).sqrt()) as f32);
            vg = (vcp / (gang).cos());
            (*t) = (((-(v / v0).ln()) / (aa * vg)) as f32);
            (*q) = (((((-0.5f64) * (z as f64)) * (v + v0)) * vg) as f32);
        }
    } else {
        /*calculating the center for a circular raypath for isotropic*/
        z0 = ((v0 / aa) as f32);
        z += z0;
        x0 = ((((-(((x * x) + (z * z)) - (z0 * z0))) as f64) / (2.0f64 * (x as f64))) as f32);
        rt = (((x0 * x0) + (z0 * z0)) as f64).sqrt();
        f1 = (1.0f64 / ((((aa * aa) * (z0 as f64)) * (z0 as f64)) / (((rt * rt) * v0) * v0)));
        so = ((((sig as f32) * z0) as f64) / rt);
        soo = so;
        /*intial trial parameters*/
        count = 1i32;
        err = 1i32;
        dso = 0.1f64;
        so4 = 0.0f64;
        sm = 1.0f64;
        /*for search effeciency*/
        if (((dvdz as f64) * cx) < 0.0f64) {
            zneg = (-1i32);
        }
        /* looping over ntries to find the spatial root through Secant method*/
        j = 0i32;
        while (j < ntries) {
            if ((if (so < 0.0f64) { (-so) } else { so }) > sm) {
                so = (0.5f64 * (((sig as f64) * sm) + so4));
            }
            if (((sig as f64) * so) < 0.0f64) {
                so = (0.5f64 * (so4 + 0.0f64));
            }
            so4 = so;
            co = ((zneg as f64) * (1.0f64 - (so * so)).sqrt());
            sw = ((cr * so) + (sr * co));
            cw = ((cr * co) - (sr * so));
            sz2 = (sw * sw);
            cz2 = (cw * cw);
            alfa = ((c2 * sz2) + (c1 * cz2));
            bbb = ((w1 * sz2) - (w2 * cz2));
            beta = ((bbb * bbb) + (((4.0f64 * sz2) * cz2) * w3)).sqrt();
            vcpi = (0.5f64 * (alfa + beta)).sqrt();
            if (((z0 as f64) - (((x as f64) * so) / co)) == 0.0f64) {
                ss = (sig as f64);
            } else {
                p11 = (((z as f64) * so) / (((z0 as f64) * co) - ((x as f64) * so)));
                ss = ((sig as f64) * ((p11 * p11) / (1.0f64 + (p11 * p11))).sqrt());
            }
            /*if beyond the ray turning point*/
            if (((z0 as f64) - (((x as f64) * so) / co)) < 0.0f64) {
                tfturn = (-1i32);
            } else {
                tfturn = 1i32;
            }
            cs = (((zneg * tfturn) as f64) * (1.0f64 - (ss * ss)).sqrt());
            si = ((cr * ss) + (sr * cs));
            ci = ((cr * cs) - (sr * ss));
            sz2 = (si * si);
            cz2 = (ci * ci);
            alfa = ((c2 * sz2) + (c1 * cz2));
            bbb = ((w1 * sz2) - (w2 * cz2));
            beta = ((bbb * bbb) + (((4.0f64 * sz2) * cz2) * w3)).sqrt();
            vcps = (0.5f64 * (alfa + beta)).sqrt();
            diff = (((z as f64) / ss) - (((z0 as f64) * (vcpi / vcps)) / so));
            if (j == 0i32) {
                if ((if (diff < 0.0f64) { (-diff) } else { diff }) < (epsx as f64)) {
                    err = 0i32;
                    break;
                }
                so1 = so;
                so = (so + ((dso * (tfturn as f64)) * (sig as f64)));
                diff1 = diff;
            } else {
                count += 1i32;
                so3 = so;
                if ((if (diff < 0.0f64) { (-diff) } else { diff }) < (epsx as f64)) {
                    err = 0i32;
                    break;
                }
                if (j == (ntries - 1i32)) {
                    zneg *= (-1i32);
                    j = (-1i32);
                    so = soo;
                    if (err == 3i32) {
                        err = 1i32;
                        break;
                    }
                    err = 3i32;
                } else {
                    so = (((diff1 * so3) - (diff * so1)) / (diff1 - diff));
                    diff1 = diff;
                    so1 = so3;
                }
            }
            j += 1;
        }
        /*ignoring unfound rays*/
        if (err != 0) {
            (*t) = 0.0f32;
            (*q) = 10000000.0f32;
            (*s) = 0.0f32;
            (*c) = 1.0f32;
        } else {
            fx = ((((-z) as f64) * cs) / ss);
            f1 = ((vcps * vcps) * ((fx * fx) + ((z * z) as f64)));
            /* signed radius of raypath; if ray going clockwise, r > 0 */
            if ((x as f64) < 0.0f64) {
                r = ((fx * fx) + ((z * z) as f64)).sqrt();
            } else {
                r = (-((fx * fx) + ((z * z) as f64)).sqrt());
            }
            /*initial ray parameters for traveltime calculations*/
            uv0 = (1.0f64 / (vcpi * v0));
            uv = (1.0f64 / (vcps * v));
            px0 = (sw * uv0);
            pz0 = (cw * uv0);
            px = (si * uv);
            pz = (ci * uv);
            /*calculating the time*/
            if (dvdx == 0.0f32) {
                umax = ((pz0 - pz) / (dvdz as f64));
            } else {
                umax = ((px0 - px) / (dvdx as f64));
            }
            nit = 0i32;
            while ((nit as f32) < nitmax) {
                if (nit == 0i32) {
                    u = 0.0f64;
                    px2 = (px0 * px0);
                    pz2 = (pz0 * pz0);
                    alpha = ((c2 * px2) + (c1 * pz2));
                    betaa = (((w1 * px2) - (w2 * pz2)) * ((w1 * px2) - (w2 * pz2)));
                    betab = (((4.0f64 * w3) * px2) * pz2);
                    betat = (betaa + betab).sqrt();
                    fctn1 = (1.0f64 / (0.5f64 * (alpha + betat)).sqrt());
                    /*fctn1 = v0;*/
                    u = umax;
                    px2 = ((px0 - ((dvdx as f64) * u)) * (px0 - ((dvdx as f64) * u)));
                    pz2 = ((pz0 - ((dvdz as f64) * u)) * (pz0 - ((dvdz as f64) * u)));
                    alpha = ((c2 * px2) + (c1 * pz2));
                    betaa = (((w1 * px2) - (w2 * pz2)) * ((w1 * px2) - (w2 * pz2)));
                    betab = (((4.0f64 * w3) * px2) * pz2);
                    betat = (betaa + betab).sqrt();
                    fctn2 = (1.0f64 / (0.5f64 * (alpha + betat)).sqrt());
                    /*fctn2 = v;*/
                    timenew = ((0.5f64 * umax) * (fctn2 + fctn1));
                    it = 1i32;
                } else {
                    tnm = (1.0f64 / (it as f64));
                    del = (umax * tnm);
                    u = (0.5f64 * del);
                    sum = 0.0f64;
                    j = 0i32;
                    while (j < it) {
                        px2 = ((px0 - ((dvdx as f64) * u)) * (px0 - ((dvdx as f64) * u)));
                        pz2 = ((pz0 - ((dvdz as f64) * u)) * (pz0 - ((dvdz as f64) * u)));
                        alpha = ((c2 * px2) + (c1 * pz2));
                        betaa = (((w1 * px2) - (w2 * pz2)) * ((w1 * px2) - (w2 * pz2)));
                        betab = (((4.0f64 * w3) * px2) * pz2);
                        betat = (betaa + betab).sqrt();
                        fctnu = (1.0f64 / (0.5f64 * (alpha + betat)).sqrt());
                        sum += fctnu;
                        u += del;
                        j += 1;
                    }
                    timenew = (0.5f64 * (timenew + ((umax * sum) * tnm)));
                    it = (2i32 * it);
                }
                if ((if ((timenew - time) < 0.0f64) {
                    (-(timenew - time))
                } else {
                    (timenew - time)
                }) < (epst as f64))
                {
                    time = (if (timenew < 0.0f64) {
                        (-timenew)
                    } else {
                        timenew
                    });
                    break;
                }
                time = timenew;
                nit += 1;
            }
            (*t) = ((if (time < 0.0f64) { (-time) } else { time }) as f32);
            /*calculating the geometrical spreading*/
            sign = ((umax / (if (umax < 0.0f64) { (-umax) } else { umax })) as i32);
            (*s) = ((si * (sign as f64)) as f32);
            (*c) = ((ci * (sign as f64)) as f32);
            mx = ((f1).sqrt() * (((cw - ci) * cr) + ((sw - si) * sr)));
            co2 = (ci * ci);
            si2 = (si * si);
            sc = (si * ci);
            dvcp = (((((0.5f64 * c2) * sc) - ((0.5f64 * c1) * sc))
                + ((0.25f64
                    * ((((2.0f64 * ((w1 * si2) - (w2 * co2))) * ((w1 * sc) + (w2 * sc)))
                        + (((4.0f64 * sc) * co2) * w3))
                        - (((4.0f64 * si2) * sc) * w3)))
                    / beta))
                / vcps);
            (*q) = ((if ((((aa * ((dvcp / vcps).atan()).cos()) * r) * mx) < 0.0f64) {
                (-(((aa * ((dvcp / vcps).atan()).cos()) * r) * mx))
            } else {
                (((aa * ((dvcp / vcps).atan()).cos()) * r) * mx)
            }) as f32);
        }
    }
    /* if coordinates were rotated, unrotate cosine and sine */
    if ((angxs as f64) != 0.0f64) {
        temp = ((cx * ((*s) as f64)) + (sx * ((*c) as f64)));
        (*c) = (((cx * ((*c) as f64)) - (sx * ((*s) as f64))) as f32);
        (*s) = (temp as f32);
    }
}
