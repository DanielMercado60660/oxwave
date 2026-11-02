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

//! Port of cwp/lib/pfafft.c.
//! Credits: Dave Hale, Colorado School of Mines, 04/27/89 (complex FFT),
//! 04/28/89 and 06/13/89 (complex size selection), 06/16/89 (real size
//! selection), 06/13/89 (real and 2-D real FFT), and 06/15/89 (multiple
//! and 2-D complex FFT); modified 08/05/91 and 11/26/89.
//! FFT cost table: Dave Hale, 08/04/91; rebuilt by Jack Cohen, 05/03/94,
//! on advice of David and Gregory Chudnovsky.
//! References: C. Temperton, Journal of Computational Physics 58 (1985),
//! 283–299, and 75 (1988), 190–198; Press et al., Numerical Recipes in C
//! (1988), p. 417.

// These kernels retain the C's float-slot indexing and radix statement order.
// Address arithmetic uses i64 to avoid signed C int overflow.
// Unsuffixed C floating literals promote only their enclosing subexpressions
// to f64; wide()/narrow() reproduce those promotions and assignment rounding.
// Single/Double also retain operand order for native NaN sign/payload selection.
#![allow(unused_parens)]
use super::float::{slots, Double, Single};

const P120: Single = Single(0.120536680);
const P142: Single = Single(0.142314838);
const P173: Single = Single(0.173648178);
const P222: Single = Single(0.222520934);
const P239: Single = Single(0.239315664);
const P281: Single = Single(0.281732557);
const P342: Single = Single(0.342020143);
const P354: Single = Single(0.354604887);
const P382: Single = Single(0.382683432);
const P415: Single = Single(0.415415013);
const P433: Single = Single(0.433883739);
const P464: Single = Single(0.464723172);
const P540: Single = Single(0.540640817);
const P559: Single = Single(0.559016994);
const P568: Single = Single(0.568064747);
const P587: Single = Single(0.587785252);
const P623: Single = Single(0.623489802);
const P642: Single = Single(0.642787610);
const P654: Single = Single(0.654860734);
const P663: Single = Single(0.663122658);
const P707: Single = Single(0.707106781);
const P748: Single = Single(0.748510748);
const P755: Single = Single(0.755749574);
const P766: Single = Single(0.766044443);
const P781: Single = Single(0.781831482);
const P822: Single = Single(0.822983866);
const P841: Single = Single(0.841253533);
const P866: Single = Single(0.866025404);
const P885: Single = Single(0.885456026);
const P900: Single = Single(0.900968868);
const P909: Single = Single(0.909631995);
const P923: Single = Single(0.923879533);
const P935: Single = Single(0.935016243);
const P939: Single = Single(0.939692621);
const P951: Single = Single(0.951056516);
const P959: Single = Single(0.959492974);
const P970: Single = Single(0.970941817);
const P974: Single = Single(0.974927912);
const P984: Single = Single(0.984807753);
const P989: Single = Single(0.989821442);
const P992: Single = Single(0.992708874);
const NFAX: i32 = 10;

pub(super) fn pfacc(isign: i32, n: i32, z: &mut [f32]) {
    let z = slots(z);
    let kfax: [i32; 10] = [16, 13, 11, 9, 8, 7, 5, 4, 3, 2];
    let mut j00: i64;
    let mut j01: i64;
    let mut j2: i64;
    let mut j3: i64;
    let mut j4: i64;
    let mut j5: i64;
    let mut j6: i64;
    let mut j7: i64;
    let mut j8: i64;
    let mut j9: i64;
    let mut j10: i64;
    let mut j11: i64;
    let mut j12: i64;
    let mut j13: i64;
    let mut j14: i64;
    let mut j15: i64;
    let mut jt: i64;
    let mut nleft: i32;
    let mut ifac: i32;
    let mut jinc: i64;
    let mut jmax: i64;
    let mut ndiv: i32;
    let mut m: i32;
    let mut mm: i32 = 0;
    let mut mu: i32 = 0;
    let mut t1r: Single;
    let mut t1i: Single;
    let mut t2r: Single;
    let mut t2i: Single;
    let mut t3r: Single;
    let mut t3i: Single;
    let mut t4r: Single;
    let mut t4i: Single;
    let mut t5r: Single;
    let mut t5i: Single;
    let mut t6r: Single;
    let mut t6i: Single;
    let mut t7r: Single;
    let mut t7i: Single;
    let mut t8r: Single;
    let mut t8i: Single;
    let mut t9r: Single;
    let mut t9i: Single;
    let mut t10r: Single;
    let mut t10i: Single;
    let mut t11r: Single;
    let mut t11i: Single;
    let mut t12r: Single;
    let mut t12i: Single;
    let mut t13r: Single;
    let mut t13i: Single;
    let mut t14r: Single;
    let mut t14i: Single;
    let mut t15r: Single;
    let mut t15i: Single;
    let mut t16r: Single;
    let mut t16i: Single;
    let mut t17r: Single;
    let mut t17i: Single;
    let mut t18r: Single;
    let mut t18i: Single;
    let mut t19r: Single;
    let mut t19i: Single;
    let mut t20r: Single;
    let mut t20i: Single;
    let mut t21r: Single;
    let mut t21i: Single;
    let mut t22r: Single;
    let mut t22i: Single;
    let mut t23r: Single;
    let mut t23i: Single;
    let mut t24r: Single;
    let mut t24i: Single;
    let mut t25r: Single;
    let mut t25i: Single;
    let mut t26r: Single;
    let mut t26i: Single;
    let mut t27r: Single;
    let mut t27i: Single;
    let mut t28r: Single;
    let mut t28i: Single;
    let mut t29r: Single;
    let mut t29i: Single;
    let mut t30r: Single;
    let mut t30i: Single;
    let mut t31r: Single;
    let mut t31i: Single;
    let mut t32r: Single;
    let mut t32i: Single;
    let mut t33r: Single;
    let mut t33i: Single;
    let mut t34r: Single;
    let mut t34i: Single;
    let mut t35r: Single;
    let mut t35i: Single;
    let mut t36r: Single;
    let mut t36i: Single;
    let mut t37r: Single;
    let mut t37i: Single;
    let mut t38r: Single;
    let mut t38i: Single;
    let mut t39r: Single;
    let mut t39i: Single;
    let mut t40r: Single;
    let mut t40i: Single;
    let mut t41r: Single;
    let mut t41i: Single;
    let mut t42r: Single;
    let mut t42i: Single;
    let mut y1r: Single;
    let mut y1i: Single;
    let mut y2r: Single;
    let mut y2i: Single;
    let mut y3r: Single;
    let mut y3i: Single;
    let mut y4r: Single;
    let mut y4i: Single;
    let mut y5r: Single;
    let mut y5i: Single;
    let mut y6r: Single;
    let mut y6i: Single;
    let mut y7r: Single;
    let mut y7i: Single;
    let mut y8r: Single;
    let mut y8i: Single;
    let mut y9r: Single;
    let mut y9i: Single;
    let mut y10r: Single;
    let mut y10i: Single;
    let mut y11r: Single;
    let mut y11i: Single;
    let mut y12r: Single;
    let mut y12i: Single;
    let mut y13r: Single;
    let mut y13i: Single;
    let mut y14r: Single;
    let mut y14i: Single;
    let mut y15r: Single;
    let mut y15i: Single;
    let mut c1: Single;
    let mut c2: Single;
    let mut c3: Single;
    let mut c4: Single;
    let mut c5: Single;
    let mut c6: Single;
    let mut c7: Single;
    let mut c8: Single;
    let mut c9: Single;
    let mut c10: Single;
    let mut c11: Single;
    let mut c12: Single;
    /* keep track of n left after dividing by factors */
    nleft = n;
    /* begin loop over possible factors (from biggest to smallest) */
    for jfax in 0..NFAX {
        /* skip if not a mutually prime factor of n */
        ifac = kfax[(jfax) as usize];
        ndiv = (nleft / ifac);
        if ((ndiv * ifac) != nleft) {
            continue;
        }
        /* update n left and determine n divided by factor */
        nleft = ndiv;
        m = (n / ifac);
        /* determine rotation factor mu and stride mm */
        for jfac in 1..=ifac {
            mu = jfac;
            mm = (jfac * m);
            if ((mm % ifac) == 1) {
                break;
            }
        }
        // Without a modular inverse C accesses beyond the transform.
        // Keep length-zero internal kernels (used by unit real transforms).
        if m != 0 && mm % ifac != 1 {
            oxwave_core::err!("pfafft.c: pfacc: invalid transform length");
        }
        /* adjust rotation factor for sign of transform */
        if (isign < 0) {
            mu = (ifac - mu);
        }
        /* compute stride, limit, and pointers */
        jinc = (2 * i64::from(mm));
        jmax = (2 * i64::from(n));
        j00 = 0;
        j01 = (j00 + jinc);
        /* if factor is 2 */
        if (ifac == 2) {
            for _l in 0..m {
                t1r = (z[(j00) as usize] - z[(j01) as usize]);
                t1i = (z[(j00 + 1) as usize] - z[(j01 + 1) as usize]);
                z[(j00) as usize] = (z[(j00) as usize] + z[(j01) as usize]);
                z[(j00 + 1) as usize] = (z[(j00 + 1) as usize] + z[(j01 + 1) as usize]);
                z[(j01) as usize] = t1r;
                z[(j01 + 1) as usize] = t1i;
                jt = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j2 = (j01 + jinc);
        if (j2 >= jmax) {
            j2 = (j2 - jmax);
        }
        /* if factor is 3 */
        if (ifac == 3) {
            if (mu == 1) {
                c1 = P866;
            } else {
                c1 = (-P866);
            }
            for _l in 0..m {
                t1r = (z[(j01) as usize] + z[(j2) as usize]);
                t1i = (z[(j01 + 1) as usize] + z[(j2 + 1) as usize]);
                y1r = (z[(j00) as usize].wide() - Double(0.5) * t1r.wide()).narrow();
                y1i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t1i.wide()).narrow();
                y2r = (c1 * (z[(j01) as usize] - z[(j2) as usize]));
                y2i = (c1 * (z[(j01 + 1) as usize] - z[(j2 + 1) as usize]));
                z[(j00) as usize] = (z[(j00) as usize] + t1r);
                z[(j00 + 1) as usize] = (z[(j00 + 1) as usize] + t1i);
                z[(j01) as usize] = (y1r - y2i);
                z[(j01 + 1) as usize] = (y1i + y2r);
                z[(j2) as usize] = (y1r + y2i);
                z[(j2 + 1) as usize] = (y1i - y2r);
                jt = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j3 = (j2 + jinc);
        if (j3 >= jmax) {
            j3 = (j3 - jmax);
        }
        /* if factor is 4 */
        if (ifac == 4) {
            if (mu == 1) {
                c1 = Double(1.0).narrow();
            } else {
                c1 = (-Double(1.0)).narrow();
            }
            for _l in 0..m {
                t1r = (z[(j00) as usize] + z[(j2) as usize]);
                t1i = (z[(j00 + 1) as usize] + z[(j2 + 1) as usize]);
                t2r = (z[(j01) as usize] + z[(j3) as usize]);
                t2i = (z[(j01 + 1) as usize] + z[(j3 + 1) as usize]);
                y1r = (z[(j00) as usize] - z[(j2) as usize]);
                y1i = (z[(j00 + 1) as usize] - z[(j2 + 1) as usize]);
                y3r = (c1 * (z[(j01) as usize] - z[(j3) as usize]));
                y3i = (c1 * (z[(j01 + 1) as usize] - z[(j3 + 1) as usize]));
                z[(j00) as usize] = (t1r + t2r);
                z[(j00 + 1) as usize] = (t1i + t2i);
                z[(j01) as usize] = (y1r - y3i);
                z[(j01 + 1) as usize] = (y1i + y3r);
                z[(j2) as usize] = (t1r - t2r);
                z[(j2 + 1) as usize] = (t1i - t2i);
                z[(j3) as usize] = (y1r + y3i);
                z[(j3 + 1) as usize] = (y1i - y3r);
                jt = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j4 = (j3 + jinc);
        if (j4 >= jmax) {
            j4 = (j4 - jmax);
        }
        /* if factor is 5 */
        if (ifac == 5) {
            if (mu == 1) {
                c1 = P559;
                c2 = P951;
                c3 = P587;
            } else if (mu == 2) {
                c1 = (-P559);
                c2 = P587;
                c3 = (-P951);
            } else if (mu == 3) {
                c1 = (-P559);
                c2 = (-P587);
                c3 = P951;
            } else {
                c1 = P559;
                c2 = (-P951);
                c3 = (-P587);
            }
            for _l in 0..m {
                t1r = (z[(j01) as usize] + z[(j4) as usize]);
                t1i = (z[(j01 + 1) as usize] + z[(j4 + 1) as usize]);
                t2r = (z[(j2) as usize] + z[(j3) as usize]);
                t2i = (z[(j2 + 1) as usize] + z[(j3 + 1) as usize]);
                t3r = (z[(j01) as usize] - z[(j4) as usize]);
                t3i = (z[(j01 + 1) as usize] - z[(j4 + 1) as usize]);
                t4r = (z[(j2) as usize] - z[(j3) as usize]);
                t4i = (z[(j2 + 1) as usize] - z[(j3 + 1) as usize]);
                t5r = (t1r + t2r);
                t5i = (t1i + t2i);
                t6r = (c1 * (t1r - t2r));
                t6i = (c1 * (t1i - t2i));
                t7r = (z[(j00) as usize].wide() - Double(0.25) * t5r.wide()).narrow();
                t7i = (z[(j00 + 1) as usize].wide() - Double(0.25) * t5i.wide()).narrow();
                y1r = (t7r + t6r);
                y1i = (t7i + t6i);
                y2r = (t7r - t6r);
                y2i = (t7i - t6i);
                y3r = ((c3 * t3r) - (c2 * t4r));
                y3i = ((c3 * t3i) - (c2 * t4i));
                y4r = ((c2 * t3r) + (c3 * t4r));
                y4i = ((c2 * t3i) + (c3 * t4i));
                z[(j00) as usize] = (z[(j00) as usize] + t5r);
                z[(j00 + 1) as usize] = (z[(j00 + 1) as usize] + t5i);
                z[(j01) as usize] = (y1r - y4i);
                z[(j01 + 1) as usize] = (y1i + y4r);
                z[(j2) as usize] = (y2r - y3i);
                z[(j2 + 1) as usize] = (y2i + y3r);
                z[(j3) as usize] = (y2r + y3i);
                z[(j3 + 1) as usize] = (y2i - y3r);
                z[(j4) as usize] = (y1r + y4i);
                z[(j4 + 1) as usize] = (y1i - y4r);
                jt = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j5 = (j4 + jinc);
        if (j5 >= jmax) {
            j5 = (j5 - jmax);
        }
        j6 = (j5 + jinc);
        if (j6 >= jmax) {
            j6 = (j6 - jmax);
        }
        /* if factor is 7 */
        if (ifac == 7) {
            if (mu == 1) {
                c1 = P623;
                c2 = (-P222);
                c3 = (-P900);
                c4 = P781;
                c5 = P974;
                c6 = P433;
            } else if (mu == 2) {
                c1 = (-P222);
                c2 = (-P900);
                c3 = P623;
                c4 = P974;
                c5 = (-P433);
                c6 = (-P781);
            } else if (mu == 3) {
                c1 = (-P900);
                c2 = P623;
                c3 = (-P222);
                c4 = P433;
                c5 = (-P781);
                c6 = P974;
            } else if (mu == 4) {
                c1 = (-P900);
                c2 = P623;
                c3 = (-P222);
                c4 = (-P433);
                c5 = P781;
                c6 = (-P974);
            } else if (mu == 5) {
                c1 = (-P222);
                c2 = (-P900);
                c3 = P623;
                c4 = (-P974);
                c5 = P433;
                c6 = P781;
            } else {
                c1 = P623;
                c2 = (-P222);
                c3 = (-P900);
                c4 = (-P781);
                c5 = (-P974);
                c6 = (-P433);
            }
            for _l in 0..m {
                t1r = (z[(j01) as usize] + z[(j6) as usize]);
                t1i = (z[(j01 + 1) as usize] + z[(j6 + 1) as usize]);
                t2r = (z[(j2) as usize] + z[(j5) as usize]);
                t2i = (z[(j2 + 1) as usize] + z[(j5 + 1) as usize]);
                t3r = (z[(j3) as usize] + z[(j4) as usize]);
                t3i = (z[(j3 + 1) as usize] + z[(j4 + 1) as usize]);
                t4r = (z[(j01) as usize] - z[(j6) as usize]);
                t4i = (z[(j01 + 1) as usize] - z[(j6 + 1) as usize]);
                t5r = (z[(j2) as usize] - z[(j5) as usize]);
                t5i = (z[(j2 + 1) as usize] - z[(j5 + 1) as usize]);
                t6r = (z[(j3) as usize] - z[(j4) as usize]);
                t6i = (z[(j3 + 1) as usize] - z[(j4 + 1) as usize]);
                t7r = (z[(j00) as usize].wide() - Double(0.5) * t3r.wide()).narrow();
                t7i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t3i.wide()).narrow();
                t8r = (t1r - t3r);
                t8i = (t1i - t3i);
                t9r = (t2r - t3r);
                t9i = (t2i - t3i);
                y1r = ((t7r + (c1 * t8r)) + (c2 * t9r));
                y1i = ((t7i + (c1 * t8i)) + (c2 * t9i));
                y2r = ((t7r + (c2 * t8r)) + (c3 * t9r));
                y2i = ((t7i + (c2 * t8i)) + (c3 * t9i));
                y3r = ((t7r + (c3 * t8r)) + (c1 * t9r));
                y3i = ((t7i + (c3 * t8i)) + (c1 * t9i));
                y4r = (((c6 * t4r) - (c4 * t5r)) + (c5 * t6r));
                y4i = (((c6 * t4i) - (c4 * t5i)) + (c5 * t6i));
                y5r = (((c5 * t4r) - (c6 * t5r)) - (c4 * t6r));
                y5i = (((c5 * t4i) - (c6 * t5i)) - (c4 * t6i));
                y6r = (((c4 * t4r) + (c5 * t5r)) + (c6 * t6r));
                y6i = (((c4 * t4i) + (c5 * t5i)) + (c6 * t6i));
                z[(j00) as usize] = (((z[(j00) as usize] + t1r) + t2r) + t3r);
                z[(j00 + 1) as usize] = (((z[(j00 + 1) as usize] + t1i) + t2i) + t3i);
                z[(j01) as usize] = (y1r - y6i);
                z[(j01 + 1) as usize] = (y1i + y6r);
                z[(j2) as usize] = (y2r - y5i);
                z[(j2 + 1) as usize] = (y2i + y5r);
                z[(j3) as usize] = (y3r - y4i);
                z[(j3 + 1) as usize] = (y3i + y4r);
                z[(j4) as usize] = (y3r + y4i);
                z[(j4 + 1) as usize] = (y3i - y4r);
                z[(j5) as usize] = (y2r + y5i);
                z[(j5 + 1) as usize] = (y2i - y5r);
                z[(j6) as usize] = (y1r + y6i);
                z[(j6 + 1) as usize] = (y1i - y6r);
                jt = (j6 + 2);
                j6 = (j5 + 2);
                j5 = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j7 = (j6 + jinc);
        if (j7 >= jmax) {
            j7 = (j7 - jmax);
        }
        /* if factor is 8 */
        if (ifac == 8) {
            if (mu == 1) {
                c1 = Double(1.0).narrow();
                c2 = P707;
            } else if (mu == 3) {
                c1 = (-Double(1.0)).narrow();
                c2 = (-P707);
            } else if (mu == 5) {
                c1 = Double(1.0).narrow();
                c2 = (-P707);
            } else {
                c1 = (-Double(1.0)).narrow();
                c2 = P707;
            }
            c3 = (c1 * c2);
            for _l in 0..m {
                t1r = (z[(j00) as usize] + z[(j4) as usize]);
                t1i = (z[(j00 + 1) as usize] + z[(j4 + 1) as usize]);
                t2r = (z[(j00) as usize] - z[(j4) as usize]);
                t2i = (z[(j00 + 1) as usize] - z[(j4 + 1) as usize]);
                t3r = (z[(j01) as usize] + z[(j5) as usize]);
                t3i = (z[(j01 + 1) as usize] + z[(j5 + 1) as usize]);
                t4r = (z[(j01) as usize] - z[(j5) as usize]);
                t4i = (z[(j01 + 1) as usize] - z[(j5 + 1) as usize]);
                t5r = (z[(j2) as usize] + z[(j6) as usize]);
                t5i = (z[(j2 + 1) as usize] + z[(j6 + 1) as usize]);
                t6r = (c1 * (z[(j2) as usize] - z[(j6) as usize]));
                t6i = (c1 * (z[(j2 + 1) as usize] - z[(j6 + 1) as usize]));
                t7r = (z[(j3) as usize] + z[(j7) as usize]);
                t7i = (z[(j3 + 1) as usize] + z[(j7 + 1) as usize]);
                t8r = (z[(j3) as usize] - z[(j7) as usize]);
                t8i = (z[(j3 + 1) as usize] - z[(j7 + 1) as usize]);
                t9r = (t1r + t5r);
                t9i = (t1i + t5i);
                t10r = (t3r + t7r);
                t10i = (t3i + t7i);
                t11r = (c2 * (t4r - t8r));
                t11i = (c2 * (t4i - t8i));
                t12r = (c3 * (t4r + t8r));
                t12i = (c3 * (t4i + t8i));
                y1r = (t2r + t11r);
                y1i = (t2i + t11i);
                y2r = (t1r - t5r);
                y2i = (t1i - t5i);
                y3r = (t2r - t11r);
                y3i = (t2i - t11i);
                y5r = (t12r - t6r);
                y5i = (t12i - t6i);
                y6r = (c1 * (t3r - t7r));
                y6i = (c1 * (t3i - t7i));
                y7r = (t12r + t6r);
                y7i = (t12i + t6i);
                z[(j00) as usize] = (t9r + t10r);
                z[(j00 + 1) as usize] = (t9i + t10i);
                z[(j01) as usize] = (y1r - y7i);
                z[(j01 + 1) as usize] = (y1i + y7r);
                z[(j2) as usize] = (y2r - y6i);
                z[(j2 + 1) as usize] = (y2i + y6r);
                z[(j3) as usize] = (y3r - y5i);
                z[(j3 + 1) as usize] = (y3i + y5r);
                z[(j4) as usize] = (t9r - t10r);
                z[(j4 + 1) as usize] = (t9i - t10i);
                z[(j5) as usize] = (y3r + y5i);
                z[(j5 + 1) as usize] = (y3i - y5r);
                z[(j6) as usize] = (y2r + y6i);
                z[(j6 + 1) as usize] = (y2i - y6r);
                z[(j7) as usize] = (y1r + y7i);
                z[(j7 + 1) as usize] = (y1i - y7r);
                jt = (j7 + 2);
                j7 = (j6 + 2);
                j6 = (j5 + 2);
                j5 = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j8 = (j7 + jinc);
        if (j8 >= jmax) {
            j8 = (j8 - jmax);
        }
        /* if factor is 9 */
        if (ifac == 9) {
            if (mu == 1) {
                c1 = P866;
                c2 = P766;
                c3 = P642;
                c4 = P173;
                c5 = P984;
            } else if (mu == 2) {
                c1 = (-P866);
                c2 = P173;
                c3 = P984;
                c4 = (-P939);
                c5 = P342;
            } else if (mu == 4) {
                c1 = P866;
                c2 = (-P939);
                c3 = P342;
                c4 = P766;
                c5 = (-P642);
            } else if (mu == 5) {
                c1 = (-P866);
                c2 = (-P939);
                c3 = (-P342);
                c4 = P766;
                c5 = P642;
            } else if (mu == 7) {
                c1 = P866;
                c2 = P173;
                c3 = (-P984);
                c4 = (-P939);
                c5 = (-P342);
            } else {
                c1 = (-P866);
                c2 = P766;
                c3 = (-P642);
                c4 = P173;
                c5 = (-P984);
            }
            c6 = (c1 * c2);
            c7 = (c1 * c3);
            c8 = (c1 * c4);
            c9 = (c1 * c5);
            for _l in 0..m {
                t1r = (z[(j3) as usize] + z[(j6) as usize]);
                t1i = (z[(j3 + 1) as usize] + z[(j6 + 1) as usize]);
                t2r = (z[(j00) as usize].wide() - Double(0.5) * t1r.wide()).narrow();
                t2i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t1i.wide()).narrow();
                t3r = (c1 * (z[(j3) as usize] - z[(j6) as usize]));
                t3i = (c1 * (z[(j3 + 1) as usize] - z[(j6 + 1) as usize]));
                t4r = (z[(j00) as usize] + t1r);
                t4i = (z[(j00 + 1) as usize] + t1i);
                t5r = (z[(j4) as usize] + z[(j7) as usize]);
                t5i = (z[(j4 + 1) as usize] + z[(j7 + 1) as usize]);
                t6r = (z[(j01) as usize].wide() - Double(0.5) * t5r.wide()).narrow();
                t6i = (z[(j01 + 1) as usize].wide() - Double(0.5) * t5i.wide()).narrow();
                t7r = (z[(j4) as usize] - z[(j7) as usize]);
                t7i = (z[(j4 + 1) as usize] - z[(j7 + 1) as usize]);
                t8r = (z[(j01) as usize] + t5r);
                t8i = (z[(j01 + 1) as usize] + t5i);
                t9r = (z[(j2) as usize] + z[(j5) as usize]);
                t9i = (z[(j2 + 1) as usize] + z[(j5 + 1) as usize]);
                t10r = (z[(j8) as usize].wide() - Double(0.5) * t9r.wide()).narrow();
                t10i = (z[(j8 + 1) as usize].wide() - Double(0.5) * t9i.wide()).narrow();
                t11r = (z[(j2) as usize] - z[(j5) as usize]);
                t11i = (z[(j2 + 1) as usize] - z[(j5 + 1) as usize]);
                t12r = (z[(j8) as usize] + t9r);
                t12i = (z[(j8 + 1) as usize] + t9i);
                t13r = (t8r + t12r);
                t13i = (t8i + t12i);
                t14r = (t6r + t10r);
                t14i = (t6i + t10i);
                t15r = (t6r - t10r);
                t15i = (t6i - t10i);
                t16r = (t7r + t11r);
                t16i = (t7i + t11i);
                t17r = (t7r - t11r);
                t17i = (t7i - t11i);
                t18r = ((c2 * t14r) - (c7 * t17r));
                t18i = ((c2 * t14i) - (c7 * t17i));
                t19r = ((c4 * t14r) + (c9 * t17r));
                t19i = ((c4 * t14i) + (c9 * t17i));
                t20r = ((c3 * t15r) + (c6 * t16r));
                t20i = ((c3 * t15i) + (c6 * t16i));
                t21r = ((c5 * t15r) - (c8 * t16r));
                t21i = ((c5 * t15i) - (c8 * t16i));
                t22r = (t18r + t19r);
                t22i = (t18i + t19i);
                t23r = (t20r - t21r);
                t23i = (t20i - t21i);
                y1r = (t2r + t18r);
                y1i = (t2i + t18i);
                y2r = (t2r + t19r);
                y2i = (t2i + t19i);
                y3r = (t4r.wide() - Double(0.5) * t13r.wide()).narrow();
                y3i = (t4i.wide() - Double(0.5) * t13i.wide()).narrow();
                y4r = (t2r - t22r);
                y4i = (t2i - t22i);
                y5r = (t3r - t23r);
                y5i = (t3i - t23i);
                y6r = (c1 * (t8r - t12r));
                y6i = (c1 * (t8i - t12i));
                y7r = (t21r - t3r);
                y7i = (t21i - t3i);
                y8r = (t3r + t20r);
                y8i = (t3i + t20i);
                z[(j00) as usize] = (t4r + t13r);
                z[(j00 + 1) as usize] = (t4i + t13i);
                z[(j01) as usize] = (y1r - y8i);
                z[(j01 + 1) as usize] = (y1i + y8r);
                z[(j2) as usize] = (y2r - y7i);
                z[(j2 + 1) as usize] = (y2i + y7r);
                z[(j3) as usize] = (y3r - y6i);
                z[(j3 + 1) as usize] = (y3i + y6r);
                z[(j4) as usize] = (y4r - y5i);
                z[(j4 + 1) as usize] = (y4i + y5r);
                z[(j5) as usize] = (y4r + y5i);
                z[(j5 + 1) as usize] = (y4i - y5r);
                z[(j6) as usize] = (y3r + y6i);
                z[(j6 + 1) as usize] = (y3i - y6r);
                z[(j7) as usize] = (y2r + y7i);
                z[(j7 + 1) as usize] = (y2i - y7r);
                z[(j8) as usize] = (y1r + y8i);
                z[(j8 + 1) as usize] = (y1i - y8r);
                jt = (j8 + 2);
                j8 = (j7 + 2);
                j7 = (j6 + 2);
                j6 = (j5 + 2);
                j5 = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j9 = (j8 + jinc);
        if (j9 >= jmax) {
            j9 = (j9 - jmax);
        }
        j10 = (j9 + jinc);
        if (j10 >= jmax) {
            j10 = (j10 - jmax);
        }
        /* if factor is 11 */
        if (ifac == 11) {
            if (mu == 1) {
                c1 = P841;
                c2 = P415;
                c3 = (-P142);
                c4 = (-P654);
                c5 = (-P959);
                c6 = P540;
                c7 = P909;
                c8 = P989;
                c9 = P755;
                c10 = P281;
            } else if (mu == 2) {
                c1 = P415;
                c2 = (-P654);
                c3 = (-P959);
                c4 = (-P142);
                c5 = P841;
                c6 = P909;
                c7 = P755;
                c8 = (-P281);
                c9 = (-P989);
                c10 = (-P540);
            } else if (mu == 3) {
                c1 = (-P142);
                c2 = (-P959);
                c3 = P415;
                c4 = P841;
                c5 = (-P654);
                c6 = P989;
                c7 = (-P281);
                c8 = (-P909);
                c9 = P540;
                c10 = P755;
            } else if (mu == 4) {
                c1 = (-P654);
                c2 = (-P142);
                c3 = P841;
                c4 = (-P959);
                c5 = P415;
                c6 = P755;
                c7 = (-P989);
                c8 = P540;
                c9 = P281;
                c10 = (-P909);
            } else if (mu == 5) {
                c1 = (-P959);
                c2 = P841;
                c3 = (-P654);
                c4 = P415;
                c5 = (-P142);
                c6 = P281;
                c7 = (-P540);
                c8 = P755;
                c9 = (-P909);
                c10 = P989;
            } else if (mu == 6) {
                c1 = (-P959);
                c2 = P841;
                c3 = (-P654);
                c4 = P415;
                c5 = (-P142);
                c6 = (-P281);
                c7 = P540;
                c8 = (-P755);
                c9 = P909;
                c10 = (-P989);
            } else if (mu == 7) {
                c1 = (-P654);
                c2 = (-P142);
                c3 = P841;
                c4 = (-P959);
                c5 = P415;
                c6 = (-P755);
                c7 = P989;
                c8 = (-P540);
                c9 = (-P281);
                c10 = P909;
            } else if (mu == 8) {
                c1 = (-P142);
                c2 = (-P959);
                c3 = P415;
                c4 = P841;
                c5 = (-P654);
                c6 = (-P989);
                c7 = P281;
                c8 = P909;
                c9 = (-P540);
                c10 = (-P755);
            } else if (mu == 9) {
                c1 = P415;
                c2 = (-P654);
                c3 = (-P959);
                c4 = (-P142);
                c5 = P841;
                c6 = (-P909);
                c7 = (-P755);
                c8 = P281;
                c9 = P989;
                c10 = P540;
            } else {
                c1 = P841;
                c2 = P415;
                c3 = (-P142);
                c4 = (-P654);
                c5 = (-P959);
                c6 = (-P540);
                c7 = (-P909);
                c8 = (-P989);
                c9 = (-P755);
                c10 = (-P281);
            }
            for _l in 0..m {
                t1r = (z[(j01) as usize] + z[(j10) as usize]);
                t1i = (z[(j01 + 1) as usize] + z[(j10 + 1) as usize]);
                t2r = (z[(j2) as usize] + z[(j9) as usize]);
                t2i = (z[(j2 + 1) as usize] + z[(j9 + 1) as usize]);
                t3r = (z[(j3) as usize] + z[(j8) as usize]);
                t3i = (z[(j3 + 1) as usize] + z[(j8 + 1) as usize]);
                t4r = (z[(j4) as usize] + z[(j7) as usize]);
                t4i = (z[(j4 + 1) as usize] + z[(j7 + 1) as usize]);
                t5r = (z[(j5) as usize] + z[(j6) as usize]);
                t5i = (z[(j5 + 1) as usize] + z[(j6 + 1) as usize]);
                t6r = (z[(j01) as usize] - z[(j10) as usize]);
                t6i = (z[(j01 + 1) as usize] - z[(j10 + 1) as usize]);
                t7r = (z[(j2) as usize] - z[(j9) as usize]);
                t7i = (z[(j2 + 1) as usize] - z[(j9 + 1) as usize]);
                t8r = (z[(j3) as usize] - z[(j8) as usize]);
                t8i = (z[(j3 + 1) as usize] - z[(j8 + 1) as usize]);
                t9r = (z[(j4) as usize] - z[(j7) as usize]);
                t9i = (z[(j4 + 1) as usize] - z[(j7 + 1) as usize]);
                t10r = (z[(j5) as usize] - z[(j6) as usize]);
                t10i = (z[(j5 + 1) as usize] - z[(j6 + 1) as usize]);
                t11r = (z[(j00) as usize].wide() - Double(0.5) * t5r.wide()).narrow();
                t11i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t5i.wide()).narrow();
                t12r = (t1r - t5r);
                t12i = (t1i - t5i);
                t13r = (t2r - t5r);
                t13i = (t2i - t5i);
                t14r = (t3r - t5r);
                t14i = (t3i - t5i);
                t15r = (t4r - t5r);
                t15i = (t4i - t5i);
                y1r = ((((t11r + (c1 * t12r)) + (c2 * t13r)) + (c3 * t14r)) + (c4 * t15r));
                y1i = ((((t11i + (c1 * t12i)) + (c2 * t13i)) + (c3 * t14i)) + (c4 * t15i));
                y2r = ((((t11r + (c2 * t12r)) + (c4 * t13r)) + (c5 * t14r)) + (c3 * t15r));
                y2i = ((((t11i + (c2 * t12i)) + (c4 * t13i)) + (c5 * t14i)) + (c3 * t15i));
                y3r = ((((t11r + (c3 * t12r)) + (c5 * t13r)) + (c2 * t14r)) + (c1 * t15r));
                y3i = ((((t11i + (c3 * t12i)) + (c5 * t13i)) + (c2 * t14i)) + (c1 * t15i));
                y4r = ((((t11r + (c4 * t12r)) + (c3 * t13r)) + (c1 * t14r)) + (c5 * t15r));
                y4i = ((((t11i + (c4 * t12i)) + (c3 * t13i)) + (c1 * t14i)) + (c5 * t15i));
                y5r = ((((t11r + (c5 * t12r)) + (c1 * t13r)) + (c4 * t14r)) + (c2 * t15r));
                y5i = ((((t11i + (c5 * t12i)) + (c1 * t13i)) + (c4 * t14i)) + (c2 * t15i));
                y6r = (((((c10 * t6r) - (c6 * t7r)) + (c9 * t8r)) - (c7 * t9r)) + (c8 * t10r));
                y6i = (((((c10 * t6i) - (c6 * t7i)) + (c9 * t8i)) - (c7 * t9i)) + (c8 * t10i));
                y7r = (((((c9 * t6r) - (c8 * t7r)) + (c6 * t8r)) + (c10 * t9r)) - (c7 * t10r));
                y7i = (((((c9 * t6i) - (c8 * t7i)) + (c6 * t8i)) + (c10 * t9i)) - (c7 * t10i));
                y8r = (((((c8 * t6r) - (c10 * t7r)) - (c7 * t8r)) + (c6 * t9r)) + (c9 * t10r));
                y8i = (((((c8 * t6i) - (c10 * t7i)) - (c7 * t8i)) + (c6 * t9i)) + (c9 * t10i));
                y9r = (((((c7 * t6r) + (c9 * t7r)) - (c10 * t8r)) - (c8 * t9r)) - (c6 * t10r));
                y9i = (((((c7 * t6i) + (c9 * t7i)) - (c10 * t8i)) - (c8 * t9i)) - (c6 * t10i));
                y10r = (((((c6 * t6r) + (c7 * t7r)) + (c8 * t8r)) + (c9 * t9r)) + (c10 * t10r));
                y10i = (((((c6 * t6i) + (c7 * t7i)) + (c8 * t8i)) + (c9 * t9i)) + (c10 * t10i));
                z[(j00) as usize] = (((((z[(j00) as usize] + t1r) + t2r) + t3r) + t4r) + t5r);
                z[(j00 + 1) as usize] =
                    (((((z[(j00 + 1) as usize] + t1i) + t2i) + t3i) + t4i) + t5i);
                z[(j01) as usize] = (y1r - y10i);
                z[(j01 + 1) as usize] = (y1i + y10r);
                z[(j2) as usize] = (y2r - y9i);
                z[(j2 + 1) as usize] = (y2i + y9r);
                z[(j3) as usize] = (y3r - y8i);
                z[(j3 + 1) as usize] = (y3i + y8r);
                z[(j4) as usize] = (y4r - y7i);
                z[(j4 + 1) as usize] = (y4i + y7r);
                z[(j5) as usize] = (y5r - y6i);
                z[(j5 + 1) as usize] = (y5i + y6r);
                z[(j6) as usize] = (y5r + y6i);
                z[(j6 + 1) as usize] = (y5i - y6r);
                z[(j7) as usize] = (y4r + y7i);
                z[(j7 + 1) as usize] = (y4i - y7r);
                z[(j8) as usize] = (y3r + y8i);
                z[(j8 + 1) as usize] = (y3i - y8r);
                z[(j9) as usize] = (y2r + y9i);
                z[(j9 + 1) as usize] = (y2i - y9r);
                z[(j10) as usize] = (y1r + y10i);
                z[(j10 + 1) as usize] = (y1i - y10r);
                jt = (j10 + 2);
                j10 = (j9 + 2);
                j9 = (j8 + 2);
                j8 = (j7 + 2);
                j7 = (j6 + 2);
                j6 = (j5 + 2);
                j5 = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j11 = (j10 + jinc);
        if (j11 >= jmax) {
            j11 = (j11 - jmax);
        }
        j12 = (j11 + jinc);
        if (j12 >= jmax) {
            j12 = (j12 - jmax);
        }
        /* if factor is 13 */
        if (ifac == 13) {
            if (mu == 1) {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = (-P354);
                c5 = (-P748);
                c6 = (-P970);
                c7 = P464;
                c8 = P822;
                c9 = P992;
                c10 = P935;
                c11 = P663;
                c12 = P239;
            } else if (mu == 2) {
                c1 = P568;
                c2 = (-P354);
                c3 = (-P970);
                c4 = (-P748);
                c5 = P120;
                c6 = P885;
                c7 = P822;
                c8 = P935;
                c9 = P239;
                c10 = (-P663);
                c11 = (-P992);
                c12 = (-P464);
            } else if (mu == 3) {
                c1 = P120;
                c2 = (-P970);
                c3 = (-P354);
                c4 = P885;
                c5 = P568;
                c6 = (-P748);
                c7 = P992;
                c8 = P239;
                c9 = (-P935);
                c10 = (-P464);
                c11 = P822;
                c12 = P663;
            } else if (mu == 4) {
                c1 = (-P354);
                c2 = (-P748);
                c3 = P885;
                c4 = P120;
                c5 = (-P970);
                c6 = P568;
                c7 = P935;
                c8 = (-P663);
                c9 = (-P464);
                c10 = P992;
                c11 = (-P239);
                c12 = (-P822);
            } else if (mu == 5) {
                c1 = (-P748);
                c2 = P120;
                c3 = P568;
                c4 = (-P970);
                c5 = P885;
                c6 = (-P354);
                c7 = P663;
                c8 = (-P992);
                c9 = P822;
                c10 = (-P239);
                c11 = (-P464);
                c12 = P935;
            } else if (mu == 6) {
                c1 = (-P970);
                c2 = P885;
                c3 = (-P748);
                c4 = P568;
                c5 = (-P354);
                c6 = P120;
                c7 = P239;
                c8 = (-P464);
                c9 = P663;
                c10 = (-P822);
                c11 = P935;
                c12 = (-P992);
            } else if (mu == 7) {
                c1 = (-P970);
                c2 = P885;
                c3 = (-P748);
                c4 = P568;
                c5 = (-P354);
                c6 = P120;
                c7 = (-P239);
                c8 = P464;
                c9 = (-P663);
                c10 = P822;
                c11 = (-P935);
                c12 = P992;
            } else if (mu == 8) {
                c1 = (-P748);
                c2 = P120;
                c3 = P568;
                c4 = (-P970);
                c5 = P885;
                c6 = (-P354);
                c7 = (-P663);
                c8 = P992;
                c9 = (-P822);
                c10 = P239;
                c11 = P464;
                c12 = (-P935);
            } else if (mu == 9) {
                c1 = (-P354);
                c2 = (-P748);
                c3 = P885;
                c4 = P120;
                c5 = (-P970);
                c6 = P568;
                c7 = (-P935);
                c8 = P663;
                c9 = P464;
                c10 = (-P992);
                c11 = P239;
                c12 = P822;
            } else if (mu == 10) {
                c1 = P120;
                c2 = (-P970);
                c3 = (-P354);
                c4 = P885;
                c5 = P568;
                c6 = (-P748);
                c7 = (-P992);
                c8 = (-P239);
                c9 = P935;
                c10 = P464;
                c11 = (-P822);
                c12 = (-P663);
            } else if (mu == 11) {
                c1 = P568;
                c2 = (-P354);
                c3 = (-P970);
                c4 = (-P748);
                c5 = P120;
                c6 = P885;
                c7 = (-P822);
                c8 = (-P935);
                c9 = (-P239);
                c10 = P663;
                c11 = P992;
                c12 = P464;
            } else {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = (-P354);
                c5 = (-P748);
                c6 = (-P970);
                c7 = (-P464);
                c8 = (-P822);
                c9 = (-P992);
                c10 = (-P935);
                c11 = (-P663);
                c12 = (-P239);
            }
            for _l in 0..m {
                t1r = (z[(j01) as usize] + z[(j12) as usize]);
                t1i = (z[(j01 + 1) as usize] + z[(j12 + 1) as usize]);
                t2r = (z[(j2) as usize] + z[(j11) as usize]);
                t2i = (z[(j2 + 1) as usize] + z[(j11 + 1) as usize]);
                t3r = (z[(j3) as usize] + z[(j10) as usize]);
                t3i = (z[(j3 + 1) as usize] + z[(j10 + 1) as usize]);
                t4r = (z[(j4) as usize] + z[(j9) as usize]);
                t4i = (z[(j4 + 1) as usize] + z[(j9 + 1) as usize]);
                t5r = (z[(j5) as usize] + z[(j8) as usize]);
                t5i = (z[(j5 + 1) as usize] + z[(j8 + 1) as usize]);
                t6r = (z[(j6) as usize] + z[(j7) as usize]);
                t6i = (z[(j6 + 1) as usize] + z[(j7 + 1) as usize]);
                t7r = (z[(j01) as usize] - z[(j12) as usize]);
                t7i = (z[(j01 + 1) as usize] - z[(j12 + 1) as usize]);
                t8r = (z[(j2) as usize] - z[(j11) as usize]);
                t8i = (z[(j2 + 1) as usize] - z[(j11 + 1) as usize]);
                t9r = (z[(j3) as usize] - z[(j10) as usize]);
                t9i = (z[(j3 + 1) as usize] - z[(j10 + 1) as usize]);
                t10r = (z[(j4) as usize] - z[(j9) as usize]);
                t10i = (z[(j4 + 1) as usize] - z[(j9 + 1) as usize]);
                t11r = (z[(j5) as usize] - z[(j8) as usize]);
                t11i = (z[(j5 + 1) as usize] - z[(j8 + 1) as usize]);
                t12r = (z[(j6) as usize] - z[(j7) as usize]);
                t12i = (z[(j6 + 1) as usize] - z[(j7 + 1) as usize]);
                t13r = (z[(j00) as usize].wide() - Double(0.5) * t6r.wide()).narrow();
                t13i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t6i.wide()).narrow();
                t14r = (t1r - t6r);
                t14i = (t1i - t6i);
                t15r = (t2r - t6r);
                t15i = (t2i - t6i);
                t16r = (t3r - t6r);
                t16i = (t3i - t6i);
                t17r = (t4r - t6r);
                t17i = (t4i - t6i);
                t18r = (t5r - t6r);
                t18i = (t5i - t6i);
                y1r = (((((t13r + (c1 * t14r)) + (c2 * t15r)) + (c3 * t16r)) + (c4 * t17r))
                    + (c5 * t18r));
                y1i = (((((t13i + (c1 * t14i)) + (c2 * t15i)) + (c3 * t16i)) + (c4 * t17i))
                    + (c5 * t18i));
                y2r = (((((t13r + (c2 * t14r)) + (c4 * t15r)) + (c6 * t16r)) + (c5 * t17r))
                    + (c3 * t18r));
                y2i = (((((t13i + (c2 * t14i)) + (c4 * t15i)) + (c6 * t16i)) + (c5 * t17i))
                    + (c3 * t18i));
                y3r = (((((t13r + (c3 * t14r)) + (c6 * t15r)) + (c4 * t16r)) + (c1 * t17r))
                    + (c2 * t18r));
                y3i = (((((t13i + (c3 * t14i)) + (c6 * t15i)) + (c4 * t16i)) + (c1 * t17i))
                    + (c2 * t18i));
                y4r = (((((t13r + (c4 * t14r)) + (c5 * t15r)) + (c1 * t16r)) + (c3 * t17r))
                    + (c6 * t18r));
                y4i = (((((t13i + (c4 * t14i)) + (c5 * t15i)) + (c1 * t16i)) + (c3 * t17i))
                    + (c6 * t18i));
                y5r = (((((t13r + (c5 * t14r)) + (c3 * t15r)) + (c2 * t16r)) + (c6 * t17r))
                    + (c1 * t18r));
                y5i = (((((t13i + (c5 * t14i)) + (c3 * t15i)) + (c2 * t16i)) + (c6 * t17i))
                    + (c1 * t18i));
                y6r = (((((t13r + (c6 * t14r)) + (c1 * t15r)) + (c5 * t16r)) + (c2 * t17r))
                    + (c4 * t18r));
                y6i = (((((t13i + (c6 * t14i)) + (c1 * t15i)) + (c5 * t16i)) + (c2 * t17i))
                    + (c4 * t18i));
                y7r = ((((((c12 * t7r) - (c7 * t8r)) + (c11 * t9r)) - (c8 * t10r)) + (c10 * t11r))
                    - (c9 * t12r));
                y7i = ((((((c12 * t7i) - (c7 * t8i)) + (c11 * t9i)) - (c8 * t10i)) + (c10 * t11i))
                    - (c9 * t12i));
                y8r = ((((((c11 * t7r) - (c9 * t8r)) + (c8 * t9r)) - (c12 * t10r)) - (c7 * t11r))
                    + (c10 * t12r));
                y8i = ((((((c11 * t7i) - (c9 * t8i)) + (c8 * t9i)) - (c12 * t10i)) - (c7 * t11i))
                    + (c10 * t12i));
                y9r = ((((((c10 * t7r) - (c11 * t8r)) - (c7 * t9r)) + (c9 * t10r)) - (c12 * t11r))
                    - (c8 * t12r));
                y9i = ((((((c10 * t7i) - (c11 * t8i)) - (c7 * t9i)) + (c9 * t10i)) - (c12 * t11i))
                    - (c8 * t12i));
                y10r = ((((((c9 * t7r) + (c12 * t8r)) - (c10 * t9r)) - (c7 * t10r)) + (c8 * t11r))
                    + (c11 * t12r));
                y10i = ((((((c9 * t7i) + (c12 * t8i)) - (c10 * t9i)) - (c7 * t10i)) + (c8 * t11i))
                    + (c11 * t12i));
                y11r = ((((((c8 * t7r) + (c10 * t8r)) + (c12 * t9r)) - (c11 * t10r))
                    - (c9 * t11r))
                    - (c7 * t12r));
                y11i = ((((((c8 * t7i) + (c10 * t8i)) + (c12 * t9i)) - (c11 * t10i))
                    - (c9 * t11i))
                    - (c7 * t12i));
                y12r = ((((((c7 * t7r) + (c8 * t8r)) + (c9 * t9r)) + (c10 * t10r)) + (c11 * t11r))
                    + (c12 * t12r));
                y12i = ((((((c7 * t7i) + (c8 * t8i)) + (c9 * t9i)) + (c10 * t10i)) + (c11 * t11i))
                    + (c12 * t12i));
                z[(j00) as usize] =
                    ((((((z[(j00) as usize] + t1r) + t2r) + t3r) + t4r) + t5r) + t6r);
                z[(j00 + 1) as usize] =
                    ((((((z[(j00 + 1) as usize] + t1i) + t2i) + t3i) + t4i) + t5i) + t6i);
                z[(j01) as usize] = (y1r - y12i);
                z[(j01 + 1) as usize] = (y1i + y12r);
                z[(j2) as usize] = (y2r - y11i);
                z[(j2 + 1) as usize] = (y2i + y11r);
                z[(j3) as usize] = (y3r - y10i);
                z[(j3 + 1) as usize] = (y3i + y10r);
                z[(j4) as usize] = (y4r - y9i);
                z[(j4 + 1) as usize] = (y4i + y9r);
                z[(j5) as usize] = (y5r - y8i);
                z[(j5 + 1) as usize] = (y5i + y8r);
                z[(j6) as usize] = (y6r - y7i);
                z[(j6 + 1) as usize] = (y6i + y7r);
                z[(j7) as usize] = (y6r + y7i);
                z[(j7 + 1) as usize] = (y6i - y7r);
                z[(j8) as usize] = (y5r + y8i);
                z[(j8 + 1) as usize] = (y5i - y8r);
                z[(j9) as usize] = (y4r + y9i);
                z[(j9 + 1) as usize] = (y4i - y9r);
                z[(j10) as usize] = (y3r + y10i);
                z[(j10 + 1) as usize] = (y3i - y10r);
                z[(j11) as usize] = (y2r + y11i);
                z[(j11 + 1) as usize] = (y2i - y11r);
                z[(j12) as usize] = (y1r + y12i);
                z[(j12 + 1) as usize] = (y1i - y12r);
                jt = (j12 + 2);
                j12 = (j11 + 2);
                j11 = (j10 + 2);
                j10 = (j9 + 2);
                j9 = (j8 + 2);
                j8 = (j7 + 2);
                j7 = (j6 + 2);
                j6 = (j5 + 2);
                j5 = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
        j13 = (j12 + jinc);
        if (j13 >= jmax) {
            j13 = (j13 - jmax);
        }
        j14 = (j13 + jinc);
        if (j14 >= jmax) {
            j14 = (j14 - jmax);
        }
        j15 = (j14 + jinc);
        if (j15 >= jmax) {
            j15 = (j15 - jmax);
        }
        /* if factor is 16 */
        if (ifac == 16) {
            if (mu == 1) {
                c1 = Double(1.0).narrow();
                c2 = P923;
                c3 = P382;
                c4 = P707;
            } else if (mu == 3) {
                c1 = (-Double(1.0)).narrow();
                c2 = P382;
                c3 = P923;
                c4 = (-P707);
            } else if (mu == 5) {
                c1 = Double(1.0).narrow();
                c2 = (-P382);
                c3 = P923;
                c4 = (-P707);
            } else if (mu == 7) {
                c1 = (-Double(1.0)).narrow();
                c2 = (-P923);
                c3 = P382;
                c4 = P707;
            } else if (mu == 9) {
                c1 = Double(1.0).narrow();
                c2 = (-P923);
                c3 = (-P382);
                c4 = P707;
            } else if (mu == 11) {
                c1 = (-Double(1.0)).narrow();
                c2 = (-P382);
                c3 = (-P923);
                c4 = (-P707);
            } else if (mu == 13) {
                c1 = Double(1.0).narrow();
                c2 = P382;
                c3 = (-P923);
                c4 = (-P707);
            } else {
                c1 = (-Double(1.0)).narrow();
                c2 = P923;
                c3 = (-P382);
                c4 = P707;
            }
            c5 = (c1 * c4);
            c6 = (c1 * c3);
            c7 = (c1 * c2);
            for _l in 0..m {
                t1r = (z[(j00) as usize] + z[(j8) as usize]);
                t1i = (z[(j00 + 1) as usize] + z[(j8 + 1) as usize]);
                t2r = (z[(j4) as usize] + z[(j12) as usize]);
                t2i = (z[(j4 + 1) as usize] + z[(j12 + 1) as usize]);
                t3r = (z[(j00) as usize] - z[(j8) as usize]);
                t3i = (z[(j00 + 1) as usize] - z[(j8 + 1) as usize]);
                t4r = (c1 * (z[(j4) as usize] - z[(j12) as usize]));
                t4i = (c1 * (z[(j4 + 1) as usize] - z[(j12 + 1) as usize]));
                t5r = (t1r + t2r);
                t5i = (t1i + t2i);
                t6r = (t1r - t2r);
                t6i = (t1i - t2i);
                t7r = (z[(j01) as usize] + z[(j9) as usize]);
                t7i = (z[(j01 + 1) as usize] + z[(j9 + 1) as usize]);
                t8r = (z[(j5) as usize] + z[(j13) as usize]);
                t8i = (z[(j5 + 1) as usize] + z[(j13 + 1) as usize]);
                t9r = (z[(j01) as usize] - z[(j9) as usize]);
                t9i = (z[(j01 + 1) as usize] - z[(j9 + 1) as usize]);
                t10r = (z[(j5) as usize] - z[(j13) as usize]);
                t10i = (z[(j5 + 1) as usize] - z[(j13 + 1) as usize]);
                t11r = (t7r + t8r);
                t11i = (t7i + t8i);
                t12r = (t7r - t8r);
                t12i = (t7i - t8i);
                t13r = (z[(j2) as usize] + z[(j10) as usize]);
                t13i = (z[(j2 + 1) as usize] + z[(j10 + 1) as usize]);
                t14r = (z[(j6) as usize] + z[(j14) as usize]);
                t14i = (z[(j6 + 1) as usize] + z[(j14 + 1) as usize]);
                t15r = (z[(j2) as usize] - z[(j10) as usize]);
                t15i = (z[(j2 + 1) as usize] - z[(j10 + 1) as usize]);
                t16r = (z[(j6) as usize] - z[(j14) as usize]);
                t16i = (z[(j6 + 1) as usize] - z[(j14 + 1) as usize]);
                t17r = (t13r + t14r);
                t17i = (t13i + t14i);
                t18r = (c4 * (t15r - t16r));
                t18i = (c4 * (t15i - t16i));
                t19r = (c5 * (t15r + t16r));
                t19i = (c5 * (t15i + t16i));
                t20r = (c1 * (t13r - t14r));
                t20i = (c1 * (t13i - t14i));
                t21r = (z[(j3) as usize] + z[(j11) as usize]);
                t21i = (z[(j3 + 1) as usize] + z[(j11 + 1) as usize]);
                t22r = (z[(j7) as usize] + z[(j15) as usize]);
                t22i = (z[(j7 + 1) as usize] + z[(j15 + 1) as usize]);
                t23r = (z[(j3) as usize] - z[(j11) as usize]);
                t23i = (z[(j3 + 1) as usize] - z[(j11 + 1) as usize]);
                t24r = (z[(j7) as usize] - z[(j15) as usize]);
                t24i = (z[(j7 + 1) as usize] - z[(j15 + 1) as usize]);
                t25r = (t21r + t22r);
                t25i = (t21i + t22i);
                t26r = (t21r - t22r);
                t26i = (t21i - t22i);
                t27r = (t9r + t24r);
                t27i = (t9i + t24i);
                t28r = (t10r + t23r);
                t28i = (t10i + t23i);
                t29r = (t9r - t24r);
                t29i = (t9i - t24i);
                t30r = (t10r - t23r);
                t30i = (t10i - t23i);
                t31r = (t5r + t17r);
                t31i = (t5i + t17i);
                t32r = (t11r + t25r);
                t32i = (t11i + t25i);
                t33r = (t3r + t18r);
                t33i = (t3i + t18i);
                t34r = ((c2 * t29r) - (c6 * t30r));
                t34i = ((c2 * t29i) - (c6 * t30i));
                t35r = (t3r - t18r);
                t35i = (t3i - t18i);
                t36r = ((c7 * t27r) - (c3 * t28r));
                t36i = ((c7 * t27i) - (c3 * t28i));
                t37r = (t4r + t19r);
                t37i = (t4i + t19i);
                t38r = ((c3 * t27r) + (c7 * t28r));
                t38i = ((c3 * t27i) + (c7 * t28i));
                t39r = (t4r - t19r);
                t39i = (t4i - t19i);
                t40r = ((c6 * t29r) + (c2 * t30r));
                t40i = ((c6 * t29i) + (c2 * t30i));
                t41r = (c4 * (t12r - t26r));
                t41i = (c4 * (t12i - t26i));
                t42r = (c5 * (t12r + t26r));
                t42i = (c5 * (t12i + t26i));
                y1r = (t33r + t34r);
                y1i = (t33i + t34i);
                y2r = (t6r + t41r);
                y2i = (t6i + t41i);
                y3r = (t35r + t40r);
                y3i = (t35i + t40i);
                y4r = (t5r - t17r);
                y4i = (t5i - t17i);
                y5r = (t35r - t40r);
                y5i = (t35i - t40i);
                y6r = (t6r - t41r);
                y6i = (t6i - t41i);
                y7r = (t33r - t34r);
                y7i = (t33i - t34i);
                y9r = (t38r - t37r);
                y9i = (t38i - t37i);
                y10r = (t42r - t20r);
                y10i = (t42i - t20i);
                y11r = (t36r + t39r);
                y11i = (t36i + t39i);
                y12r = (c1 * (t11r - t25r));
                y12i = (c1 * (t11i - t25i));
                y13r = (t36r - t39r);
                y13i = (t36i - t39i);
                y14r = (t42r + t20r);
                y14i = (t42i + t20i);
                y15r = (t38r + t37r);
                y15i = (t38i + t37i);
                z[(j00) as usize] = (t31r + t32r);
                z[(j00 + 1) as usize] = (t31i + t32i);
                z[(j01) as usize] = (y1r - y15i);
                z[(j01 + 1) as usize] = (y1i + y15r);
                z[(j2) as usize] = (y2r - y14i);
                z[(j2 + 1) as usize] = (y2i + y14r);
                z[(j3) as usize] = (y3r - y13i);
                z[(j3 + 1) as usize] = (y3i + y13r);
                z[(j4) as usize] = (y4r - y12i);
                z[(j4 + 1) as usize] = (y4i + y12r);
                z[(j5) as usize] = (y5r - y11i);
                z[(j5 + 1) as usize] = (y5i + y11r);
                z[(j6) as usize] = (y6r - y10i);
                z[(j6 + 1) as usize] = (y6i + y10r);
                z[(j7) as usize] = (y7r - y9i);
                z[(j7 + 1) as usize] = (y7i + y9r);
                z[(j8) as usize] = (t31r - t32r);
                z[(j8 + 1) as usize] = (t31i - t32i);
                z[(j9) as usize] = (y7r + y9i);
                z[(j9 + 1) as usize] = (y7i - y9r);
                z[(j10) as usize] = (y6r + y10i);
                z[(j10 + 1) as usize] = (y6i - y10r);
                z[(j11) as usize] = (y5r + y11i);
                z[(j11 + 1) as usize] = (y5i - y11r);
                z[(j12) as usize] = (y4r + y12i);
                z[(j12 + 1) as usize] = (y4i - y12r);
                z[(j13) as usize] = (y3r + y13i);
                z[(j13 + 1) as usize] = (y3i - y13r);
                z[(j14) as usize] = (y2r + y14i);
                z[(j14 + 1) as usize] = (y2i - y14r);
                z[(j15) as usize] = (y1r + y15i);
                z[(j15 + 1) as usize] = (y1i - y15r);
                jt = (j15 + 2);
                j15 = (j14 + 2);
                j14 = (j13 + 2);
                j13 = (j12 + 2);
                j12 = (j11 + 2);
                j11 = (j10 + 2);
                j10 = (j9 + 2);
                j9 = (j8 + 2);
                j8 = (j7 + 2);
                j7 = (j6 + 2);
                j6 = (j5 + 2);
                j5 = (j4 + 2);
                j4 = (j3 + 2);
                j3 = (j2 + 2);
                j2 = (j01 + 2);
                j01 = (j00 + 2);
                j00 = jt;
            }
            continue;
        }
    }
}

pub(super) fn pfamcc(isign: i32, n: i32, nt: i32, k: i32, kt: i32, z: &mut [f32]) {
    let z = slots(z);
    let kfax: [i32; 10] = [16, 13, 11, 9, 8, 7, 5, 4, 3, 2];
    let mut j00: i64;
    let mut j01: i64;
    let mut j2: i64;
    let mut j3: i64;
    let mut j4: i64;
    let mut j5: i64;
    let mut j6: i64;
    let mut j7: i64;
    let mut j8: i64;
    let mut j9: i64;
    let mut j10: i64;
    let mut j11: i64;
    let mut j12: i64;
    let mut j13: i64;
    let mut j14: i64;
    let mut j15: i64;
    let mut nleft: i32;
    let mut ifac: i32;
    let mut iinc: i64;
    let mut imax: i64;
    let mut ndiv: i32;
    let mut m: i32;
    let mut mm: i32 = 0;
    let mut mu: i32 = 0;
    let istep: i64;
    let jstep: i64;
    let mut i0: i64;
    let mut i1: i64;
    let mut i2: i64;
    let mut i3: i64;
    let mut i4: i64;
    let mut i5: i64;
    let mut i6: i64;
    let mut i7: i64;
    let mut i8: i64;
    let mut i9: i64;
    let mut i10: i64;
    let mut i11: i64;
    let mut i12: i64;
    let mut i13: i64;
    let mut i14: i64;
    let mut i15: i64;
    let mut it: i64;
    let mut t1r: Single;
    let mut t1i: Single;
    let mut t2r: Single;
    let mut t2i: Single;
    let mut t3r: Single;
    let mut t3i: Single;
    let mut t4r: Single;
    let mut t4i: Single;
    let mut t5r: Single;
    let mut t5i: Single;
    let mut t6r: Single;
    let mut t6i: Single;
    let mut t7r: Single;
    let mut t7i: Single;
    let mut t8r: Single;
    let mut t8i: Single;
    let mut t9r: Single;
    let mut t9i: Single;
    let mut t10r: Single;
    let mut t10i: Single;
    let mut t11r: Single;
    let mut t11i: Single;
    let mut t12r: Single;
    let mut t12i: Single;
    let mut t13r: Single;
    let mut t13i: Single;
    let mut t14r: Single;
    let mut t14i: Single;
    let mut t15r: Single;
    let mut t15i: Single;
    let mut t16r: Single;
    let mut t16i: Single;
    let mut t17r: Single;
    let mut t17i: Single;
    let mut t18r: Single;
    let mut t18i: Single;
    let mut t19r: Single;
    let mut t19i: Single;
    let mut t20r: Single;
    let mut t20i: Single;
    let mut t21r: Single;
    let mut t21i: Single;
    let mut t22r: Single;
    let mut t22i: Single;
    let mut t23r: Single;
    let mut t23i: Single;
    let mut t24r: Single;
    let mut t24i: Single;
    let mut t25r: Single;
    let mut t25i: Single;
    let mut t26r: Single;
    let mut t26i: Single;
    let mut t27r: Single;
    let mut t27i: Single;
    let mut t28r: Single;
    let mut t28i: Single;
    let mut t29r: Single;
    let mut t29i: Single;
    let mut t30r: Single;
    let mut t30i: Single;
    let mut t31r: Single;
    let mut t31i: Single;
    let mut t32r: Single;
    let mut t32i: Single;
    let mut t33r: Single;
    let mut t33i: Single;
    let mut t34r: Single;
    let mut t34i: Single;
    let mut t35r: Single;
    let mut t35i: Single;
    let mut t36r: Single;
    let mut t36i: Single;
    let mut t37r: Single;
    let mut t37i: Single;
    let mut t38r: Single;
    let mut t38i: Single;
    let mut t39r: Single;
    let mut t39i: Single;
    let mut t40r: Single;
    let mut t40i: Single;
    let mut t41r: Single;
    let mut t41i: Single;
    let mut t42r: Single;
    let mut t42i: Single;
    let mut y1r: Single;
    let mut y1i: Single;
    let mut y2r: Single;
    let mut y2i: Single;
    let mut y3r: Single;
    let mut y3i: Single;
    let mut y4r: Single;
    let mut y4i: Single;
    let mut y5r: Single;
    let mut y5i: Single;
    let mut y6r: Single;
    let mut y6i: Single;
    let mut y7r: Single;
    let mut y7i: Single;
    let mut y8r: Single;
    let mut y8i: Single;
    let mut y9r: Single;
    let mut y9i: Single;
    let mut y10r: Single;
    let mut y10i: Single;
    let mut y11r: Single;
    let mut y11i: Single;
    let mut y12r: Single;
    let mut y12i: Single;
    let mut y13r: Single;
    let mut y13i: Single;
    let mut y14r: Single;
    let mut y14i: Single;
    let mut y15r: Single;
    let mut y15i: Single;
    let mut c1: Single;
    let mut c2: Single;
    let mut c3: Single;
    let mut c4: Single;
    let mut c5: Single;
    let mut c6: Single;
    let mut c7: Single;
    let mut c8: Single;
    let mut c9: Single;
    let mut c10: Single;
    let mut c11: Single;
    let mut c12: Single;
    // Float-slot strides and indices are widened to avoid C int overflow,
    // including unused strides and the final increment after a transform.
    // Floating-point operations and their order remain identical to C.
    /* determine step within and between transforms */
    istep = (2 * i64::from(k));
    jstep = (2 * i64::from(kt));
    /* keep track of n left after dividing by factors */
    nleft = n;
    /* begin loop over possible factors (from biggest to smallest) */
    for jfax in 0..NFAX {
        /* skip if not a mutually prime factor of n */
        ifac = kfax[(jfax) as usize];
        ndiv = (nleft / ifac);
        if ((ndiv * ifac) != nleft) {
            continue;
        }
        /* update n left and determine n divided by factor */
        nleft = ndiv;
        m = (n / ifac);
        /* determine rotation factor mu and stride mm */
        for jfac in 1..=ifac {
            mu = jfac;
            mm = (jfac * m);
            if ((mm % ifac) == 1) {
                break;
            }
        }
        // With zero sample stride all accesses alias slot zero, as in C.
        if m != 0 && istep != 0 && mm % ifac != 1 {
            oxwave_core::err!("pfafft.c: pfamcc: invalid transform length");
        }
        /* adjust rotation factor for sign of transform */
        if (isign < 0) {
            mu = (ifac - mu);
        }
        /* compute stride, limit, and pointers */
        iinc = (istep * i64::from(mm));
        imax = (istep * i64::from(n));
        i0 = 0;
        i1 = (i0 + iinc);
        /* if factor is 2 */
        if (ifac == 2) {
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                for _jt in 0..nt {
                    t1r = (z[(j00) as usize] - z[(j01) as usize]);
                    t1i = (z[(j00 + 1) as usize] - z[(j01 + 1) as usize]);
                    z[(j00) as usize] = (z[(j00) as usize] + z[(j01) as usize]);
                    z[(j00 + 1) as usize] = (z[(j00 + 1) as usize] + z[(j01 + 1) as usize]);
                    z[(j01) as usize] = t1r;
                    z[(j01 + 1) as usize] = t1i;
                    j00 += jstep;
                    j01 += jstep;
                }
                it = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i2 = (i1 + iinc);
        if (i2 >= imax) {
            i2 = (i2 - imax);
        }
        /* if factor is 3 */
        if (ifac == 3) {
            if (mu == 1) {
                c1 = P866;
            } else {
                c1 = (-P866);
            }
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                for _jt in 0..nt {
                    t1r = (z[(j01) as usize] + z[(j2) as usize]);
                    t1i = (z[(j01 + 1) as usize] + z[(j2 + 1) as usize]);
                    y1r = (z[(j00) as usize].wide() - Double(0.5) * t1r.wide()).narrow();
                    y1i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t1i.wide()).narrow();
                    y2r = (c1 * (z[(j01) as usize] - z[(j2) as usize]));
                    y2i = (c1 * (z[(j01 + 1) as usize] - z[(j2 + 1) as usize]));
                    z[(j00) as usize] = (z[(j00) as usize] + t1r);
                    z[(j00 + 1) as usize] = (z[(j00 + 1) as usize] + t1i);
                    z[(j01) as usize] = (y1r - y2i);
                    z[(j01 + 1) as usize] = (y1i + y2r);
                    z[(j2) as usize] = (y1r + y2i);
                    z[(j2 + 1) as usize] = (y1i - y2r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                }
                it = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i3 = (i2 + iinc);
        if (i3 >= imax) {
            i3 = (i3 - imax);
        }
        /* if factor is 4 */
        if (ifac == 4) {
            if (mu == 1) {
                c1 = Double(1.0).narrow();
            } else {
                c1 = (-Double(1.0)).narrow();
            }
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                for _jt in 0..nt {
                    t1r = (z[(j00) as usize] + z[(j2) as usize]);
                    t1i = (z[(j00 + 1) as usize] + z[(j2 + 1) as usize]);
                    t2r = (z[(j01) as usize] + z[(j3) as usize]);
                    t2i = (z[(j01 + 1) as usize] + z[(j3 + 1) as usize]);
                    y1r = (z[(j00) as usize] - z[(j2) as usize]);
                    y1i = (z[(j00 + 1) as usize] - z[(j2 + 1) as usize]);
                    y3r = (c1 * (z[(j01) as usize] - z[(j3) as usize]));
                    y3i = (c1 * (z[(j01 + 1) as usize] - z[(j3 + 1) as usize]));
                    z[(j00) as usize] = (t1r + t2r);
                    z[(j00 + 1) as usize] = (t1i + t2i);
                    z[(j01) as usize] = (y1r - y3i);
                    z[(j01 + 1) as usize] = (y1i + y3r);
                    z[(j2) as usize] = (t1r - t2r);
                    z[(j2 + 1) as usize] = (t1i - t2i);
                    z[(j3) as usize] = (y1r + y3i);
                    z[(j3 + 1) as usize] = (y1i - y3r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                }
                it = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i4 = (i3 + iinc);
        if (i4 >= imax) {
            i4 = (i4 - imax);
        }
        /* if factor is 5 */
        if (ifac == 5) {
            if (mu == 1) {
                c1 = P559;
                c2 = P951;
                c3 = P587;
            } else if (mu == 2) {
                c1 = (-P559);
                c2 = P587;
                c3 = (-P951);
            } else if (mu == 3) {
                c1 = (-P559);
                c2 = (-P587);
                c3 = P951;
            } else {
                c1 = P559;
                c2 = (-P951);
                c3 = (-P587);
            }
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                for _jt in 0..nt {
                    t1r = (z[(j01) as usize] + z[(j4) as usize]);
                    t1i = (z[(j01 + 1) as usize] + z[(j4 + 1) as usize]);
                    t2r = (z[(j2) as usize] + z[(j3) as usize]);
                    t2i = (z[(j2 + 1) as usize] + z[(j3 + 1) as usize]);
                    t3r = (z[(j01) as usize] - z[(j4) as usize]);
                    t3i = (z[(j01 + 1) as usize] - z[(j4 + 1) as usize]);
                    t4r = (z[(j2) as usize] - z[(j3) as usize]);
                    t4i = (z[(j2 + 1) as usize] - z[(j3 + 1) as usize]);
                    t5r = (t1r + t2r);
                    t5i = (t1i + t2i);
                    t6r = (c1 * (t1r - t2r));
                    t6i = (c1 * (t1i - t2i));
                    t7r = (z[(j00) as usize].wide() - Double(0.25) * t5r.wide()).narrow();
                    t7i = (z[(j00 + 1) as usize].wide() - Double(0.25) * t5i.wide()).narrow();
                    y1r = (t7r + t6r);
                    y1i = (t7i + t6i);
                    y2r = (t7r - t6r);
                    y2i = (t7i - t6i);
                    y3r = ((c3 * t3r) - (c2 * t4r));
                    y3i = ((c3 * t3i) - (c2 * t4i));
                    y4r = ((c2 * t3r) + (c3 * t4r));
                    y4i = ((c2 * t3i) + (c3 * t4i));
                    z[(j00) as usize] = (z[(j00) as usize] + t5r);
                    z[(j00 + 1) as usize] = (z[(j00 + 1) as usize] + t5i);
                    z[(j01) as usize] = (y1r - y4i);
                    z[(j01 + 1) as usize] = (y1i + y4r);
                    z[(j2) as usize] = (y2r - y3i);
                    z[(j2 + 1) as usize] = (y2i + y3r);
                    z[(j3) as usize] = (y2r + y3i);
                    z[(j3 + 1) as usize] = (y2i - y3r);
                    z[(j4) as usize] = (y1r + y4i);
                    z[(j4 + 1) as usize] = (y1i - y4r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                }
                it = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i5 = (i4 + iinc);
        if (i5 >= imax) {
            i5 = (i5 - imax);
        }
        i6 = (i5 + iinc);
        if (i6 >= imax) {
            i6 = (i6 - imax);
        }
        /* if factor is 7 */
        if (ifac == 7) {
            if (mu == 1) {
                c1 = P623;
                c2 = (-P222);
                c3 = (-P900);
                c4 = P781;
                c5 = P974;
                c6 = P433;
            } else if (mu == 2) {
                c1 = (-P222);
                c2 = (-P900);
                c3 = P623;
                c4 = P974;
                c5 = (-P433);
                c6 = (-P781);
            } else if (mu == 3) {
                c1 = (-P900);
                c2 = P623;
                c3 = (-P222);
                c4 = P433;
                c5 = (-P781);
                c6 = P974;
            } else if (mu == 4) {
                c1 = (-P900);
                c2 = P623;
                c3 = (-P222);
                c4 = (-P433);
                c5 = P781;
                c6 = (-P974);
            } else if (mu == 5) {
                c1 = (-P222);
                c2 = (-P900);
                c3 = P623;
                c4 = (-P974);
                c5 = P433;
                c6 = P781;
            } else {
                c1 = P623;
                c2 = (-P222);
                c3 = (-P900);
                c4 = (-P781);
                c5 = (-P974);
                c6 = (-P433);
            }
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                for _jt in 0..nt {
                    t1r = (z[(j01) as usize] + z[(j6) as usize]);
                    t1i = (z[(j01 + 1) as usize] + z[(j6 + 1) as usize]);
                    t2r = (z[(j2) as usize] + z[(j5) as usize]);
                    t2i = (z[(j2 + 1) as usize] + z[(j5 + 1) as usize]);
                    t3r = (z[(j3) as usize] + z[(j4) as usize]);
                    t3i = (z[(j3 + 1) as usize] + z[(j4 + 1) as usize]);
                    t4r = (z[(j01) as usize] - z[(j6) as usize]);
                    t4i = (z[(j01 + 1) as usize] - z[(j6 + 1) as usize]);
                    t5r = (z[(j2) as usize] - z[(j5) as usize]);
                    t5i = (z[(j2 + 1) as usize] - z[(j5 + 1) as usize]);
                    t6r = (z[(j3) as usize] - z[(j4) as usize]);
                    t6i = (z[(j3 + 1) as usize] - z[(j4 + 1) as usize]);
                    t7r = (z[(j00) as usize].wide() - Double(0.5) * t3r.wide()).narrow();
                    t7i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t3i.wide()).narrow();
                    t8r = (t1r - t3r);
                    t8i = (t1i - t3i);
                    t9r = (t2r - t3r);
                    t9i = (t2i - t3i);
                    y1r = ((t7r + (c1 * t8r)) + (c2 * t9r));
                    y1i = ((t7i + (c1 * t8i)) + (c2 * t9i));
                    y2r = ((t7r + (c2 * t8r)) + (c3 * t9r));
                    y2i = ((t7i + (c2 * t8i)) + (c3 * t9i));
                    y3r = ((t7r + (c3 * t8r)) + (c1 * t9r));
                    y3i = ((t7i + (c3 * t8i)) + (c1 * t9i));
                    y4r = (((c6 * t4r) - (c4 * t5r)) + (c5 * t6r));
                    y4i = (((c6 * t4i) - (c4 * t5i)) + (c5 * t6i));
                    y5r = (((c5 * t4r) - (c6 * t5r)) - (c4 * t6r));
                    y5i = (((c5 * t4i) - (c6 * t5i)) - (c4 * t6i));
                    y6r = (((c4 * t4r) + (c5 * t5r)) + (c6 * t6r));
                    y6i = (((c4 * t4i) + (c5 * t5i)) + (c6 * t6i));
                    z[(j00) as usize] = (((z[(j00) as usize] + t1r) + t2r) + t3r);
                    z[(j00 + 1) as usize] = (((z[(j00 + 1) as usize] + t1i) + t2i) + t3i);
                    z[(j01) as usize] = (y1r - y6i);
                    z[(j01 + 1) as usize] = (y1i + y6r);
                    z[(j2) as usize] = (y2r - y5i);
                    z[(j2 + 1) as usize] = (y2i + y5r);
                    z[(j3) as usize] = (y3r - y4i);
                    z[(j3 + 1) as usize] = (y3i + y4r);
                    z[(j4) as usize] = (y3r + y4i);
                    z[(j4 + 1) as usize] = (y3i - y4r);
                    z[(j5) as usize] = (y2r + y5i);
                    z[(j5 + 1) as usize] = (y2i - y5r);
                    z[(j6) as usize] = (y1r + y6i);
                    z[(j6 + 1) as usize] = (y1i - y6r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                }
                it = (i6 + istep);
                i6 = (i5 + istep);
                i5 = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i7 = (i6 + iinc);
        if (i7 >= imax) {
            i7 = (i7 - imax);
        }
        /* if factor is 8 */
        if (ifac == 8) {
            if (mu == 1) {
                c1 = Double(1.0).narrow();
                c2 = P707;
            } else if (mu == 3) {
                c1 = (-Double(1.0)).narrow();
                c2 = (-P707);
            } else if (mu == 5) {
                c1 = Double(1.0).narrow();
                c2 = (-P707);
            } else {
                c1 = (-Double(1.0)).narrow();
                c2 = P707;
            }
            c3 = (c1 * c2);
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                for _jt in 0..nt {
                    t1r = (z[(j00) as usize] + z[(j4) as usize]);
                    t1i = (z[(j00 + 1) as usize] + z[(j4 + 1) as usize]);
                    t2r = (z[(j00) as usize] - z[(j4) as usize]);
                    t2i = (z[(j00 + 1) as usize] - z[(j4 + 1) as usize]);
                    t3r = (z[(j01) as usize] + z[(j5) as usize]);
                    t3i = (z[(j01 + 1) as usize] + z[(j5 + 1) as usize]);
                    t4r = (z[(j01) as usize] - z[(j5) as usize]);
                    t4i = (z[(j01 + 1) as usize] - z[(j5 + 1) as usize]);
                    t5r = (z[(j2) as usize] + z[(j6) as usize]);
                    t5i = (z[(j2 + 1) as usize] + z[(j6 + 1) as usize]);
                    t6r = (c1 * (z[(j2) as usize] - z[(j6) as usize]));
                    t6i = (c1 * (z[(j2 + 1) as usize] - z[(j6 + 1) as usize]));
                    t7r = (z[(j3) as usize] + z[(j7) as usize]);
                    t7i = (z[(j3 + 1) as usize] + z[(j7 + 1) as usize]);
                    t8r = (z[(j3) as usize] - z[(j7) as usize]);
                    t8i = (z[(j3 + 1) as usize] - z[(j7 + 1) as usize]);
                    t9r = (t1r + t5r);
                    t9i = (t1i + t5i);
                    t10r = (t3r + t7r);
                    t10i = (t3i + t7i);
                    t11r = (c2 * (t4r - t8r));
                    t11i = (c2 * (t4i - t8i));
                    t12r = (c3 * (t4r + t8r));
                    t12i = (c3 * (t4i + t8i));
                    y1r = (t2r + t11r);
                    y1i = (t2i + t11i);
                    y2r = (t1r - t5r);
                    y2i = (t1i - t5i);
                    y3r = (t2r - t11r);
                    y3i = (t2i - t11i);
                    y5r = (t12r - t6r);
                    y5i = (t12i - t6i);
                    y6r = (c1 * (t3r - t7r));
                    y6i = (c1 * (t3i - t7i));
                    y7r = (t12r + t6r);
                    y7i = (t12i + t6i);
                    z[(j00) as usize] = (t9r + t10r);
                    z[(j00 + 1) as usize] = (t9i + t10i);
                    z[(j01) as usize] = (y1r - y7i);
                    z[(j01 + 1) as usize] = (y1i + y7r);
                    z[(j2) as usize] = (y2r - y6i);
                    z[(j2 + 1) as usize] = (y2i + y6r);
                    z[(j3) as usize] = (y3r - y5i);
                    z[(j3 + 1) as usize] = (y3i + y5r);
                    z[(j4) as usize] = (t9r - t10r);
                    z[(j4 + 1) as usize] = (t9i - t10i);
                    z[(j5) as usize] = (y3r + y5i);
                    z[(j5 + 1) as usize] = (y3i - y5r);
                    z[(j6) as usize] = (y2r + y6i);
                    z[(j6 + 1) as usize] = (y2i - y6r);
                    z[(j7) as usize] = (y1r + y7i);
                    z[(j7 + 1) as usize] = (y1i - y7r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                }
                it = (i7 + istep);
                i7 = (i6 + istep);
                i6 = (i5 + istep);
                i5 = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i8 = (i7 + iinc);
        if (i8 >= imax) {
            i8 = (i8 - imax);
        }
        /* if factor is 9 */
        if (ifac == 9) {
            if (mu == 1) {
                c1 = P866;
                c2 = P766;
                c3 = P642;
                c4 = P173;
                c5 = P984;
            } else if (mu == 2) {
                c1 = (-P866);
                c2 = P173;
                c3 = P984;
                c4 = (-P939);
                c5 = P342;
            } else if (mu == 4) {
                c1 = P866;
                c2 = (-P939);
                c3 = P342;
                c4 = P766;
                c5 = (-P642);
            } else if (mu == 5) {
                c1 = (-P866);
                c2 = (-P939);
                c3 = (-P342);
                c4 = P766;
                c5 = P642;
            } else if (mu == 7) {
                c1 = P866;
                c2 = P173;
                c3 = (-P984);
                c4 = (-P939);
                c5 = (-P342);
            } else {
                c1 = (-P866);
                c2 = P766;
                c3 = (-P642);
                c4 = P173;
                c5 = (-P984);
            }
            c6 = (c1 * c2);
            c7 = (c1 * c3);
            c8 = (c1 * c4);
            c9 = (c1 * c5);
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                for _jt in 0..nt {
                    t1r = (z[(j3) as usize] + z[(j6) as usize]);
                    t1i = (z[(j3 + 1) as usize] + z[(j6 + 1) as usize]);
                    t2r = (z[(j00) as usize].wide() - Double(0.5) * t1r.wide()).narrow();
                    t2i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t1i.wide()).narrow();
                    t3r = (c1 * (z[(j3) as usize] - z[(j6) as usize]));
                    t3i = (c1 * (z[(j3 + 1) as usize] - z[(j6 + 1) as usize]));
                    t4r = (z[(j00) as usize] + t1r);
                    t4i = (z[(j00 + 1) as usize] + t1i);
                    t5r = (z[(j4) as usize] + z[(j7) as usize]);
                    t5i = (z[(j4 + 1) as usize] + z[(j7 + 1) as usize]);
                    t6r = (z[(j01) as usize].wide() - Double(0.5) * t5r.wide()).narrow();
                    t6i = (z[(j01 + 1) as usize].wide() - Double(0.5) * t5i.wide()).narrow();
                    t7r = (z[(j4) as usize] - z[(j7) as usize]);
                    t7i = (z[(j4 + 1) as usize] - z[(j7 + 1) as usize]);
                    t8r = (z[(j01) as usize] + t5r);
                    t8i = (z[(j01 + 1) as usize] + t5i);
                    t9r = (z[(j2) as usize] + z[(j5) as usize]);
                    t9i = (z[(j2 + 1) as usize] + z[(j5 + 1) as usize]);
                    t10r = (z[(j8) as usize].wide() - Double(0.5) * t9r.wide()).narrow();
                    t10i = (z[(j8 + 1) as usize].wide() - Double(0.5) * t9i.wide()).narrow();
                    t11r = (z[(j2) as usize] - z[(j5) as usize]);
                    t11i = (z[(j2 + 1) as usize] - z[(j5 + 1) as usize]);
                    t12r = (z[(j8) as usize] + t9r);
                    t12i = (z[(j8 + 1) as usize] + t9i);
                    t13r = (t8r + t12r);
                    t13i = (t8i + t12i);
                    t14r = (t6r + t10r);
                    t14i = (t6i + t10i);
                    t15r = (t6r - t10r);
                    t15i = (t6i - t10i);
                    t16r = (t7r + t11r);
                    t16i = (t7i + t11i);
                    t17r = (t7r - t11r);
                    t17i = (t7i - t11i);
                    t18r = ((c2 * t14r) - (c7 * t17r));
                    t18i = ((c2 * t14i) - (c7 * t17i));
                    t19r = ((c4 * t14r) + (c9 * t17r));
                    t19i = ((c4 * t14i) + (c9 * t17i));
                    t20r = ((c3 * t15r) + (c6 * t16r));
                    t20i = ((c3 * t15i) + (c6 * t16i));
                    t21r = ((c5 * t15r) - (c8 * t16r));
                    t21i = ((c5 * t15i) - (c8 * t16i));
                    t22r = (t18r + t19r);
                    t22i = (t18i + t19i);
                    t23r = (t20r - t21r);
                    t23i = (t20i - t21i);
                    y1r = (t2r + t18r);
                    y1i = (t2i + t18i);
                    y2r = (t2r + t19r);
                    y2i = (t2i + t19i);
                    y3r = (t4r.wide() - Double(0.5) * t13r.wide()).narrow();
                    y3i = (t4i.wide() - Double(0.5) * t13i.wide()).narrow();
                    y4r = (t2r - t22r);
                    y4i = (t2i - t22i);
                    y5r = (t3r - t23r);
                    y5i = (t3i - t23i);
                    y6r = (c1 * (t8r - t12r));
                    y6i = (c1 * (t8i - t12i));
                    y7r = (t21r - t3r);
                    y7i = (t21i - t3i);
                    y8r = (t3r + t20r);
                    y8i = (t3i + t20i);
                    z[(j00) as usize] = (t4r + t13r);
                    z[(j00 + 1) as usize] = (t4i + t13i);
                    z[(j01) as usize] = (y1r - y8i);
                    z[(j01 + 1) as usize] = (y1i + y8r);
                    z[(j2) as usize] = (y2r - y7i);
                    z[(j2 + 1) as usize] = (y2i + y7r);
                    z[(j3) as usize] = (y3r - y6i);
                    z[(j3 + 1) as usize] = (y3i + y6r);
                    z[(j4) as usize] = (y4r - y5i);
                    z[(j4 + 1) as usize] = (y4i + y5r);
                    z[(j5) as usize] = (y4r + y5i);
                    z[(j5 + 1) as usize] = (y4i - y5r);
                    z[(j6) as usize] = (y3r + y6i);
                    z[(j6 + 1) as usize] = (y3i - y6r);
                    z[(j7) as usize] = (y2r + y7i);
                    z[(j7 + 1) as usize] = (y2i - y7r);
                    z[(j8) as usize] = (y1r + y8i);
                    z[(j8 + 1) as usize] = (y1i - y8r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                }
                it = (i8 + istep);
                i8 = (i7 + istep);
                i7 = (i6 + istep);
                i6 = (i5 + istep);
                i5 = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i9 = (i8 + iinc);
        if (i9 >= imax) {
            i9 = (i9 - imax);
        }
        i10 = (i9 + iinc);
        if (i10 >= imax) {
            i10 = (i10 - imax);
        }
        /* if factor is 11 */
        if (ifac == 11) {
            if (mu == 1) {
                c1 = P841;
                c2 = P415;
                c3 = (-P142);
                c4 = (-P654);
                c5 = (-P959);
                c6 = P540;
                c7 = P909;
                c8 = P989;
                c9 = P755;
                c10 = P281;
            } else if (mu == 2) {
                c1 = P415;
                c2 = (-P654);
                c3 = (-P959);
                c4 = (-P142);
                c5 = P841;
                c6 = P909;
                c7 = P755;
                c8 = (-P281);
                c9 = (-P989);
                c10 = (-P540);
            } else if (mu == 3) {
                c1 = (-P142);
                c2 = (-P959);
                c3 = P415;
                c4 = P841;
                c5 = (-P654);
                c6 = P989;
                c7 = (-P281);
                c8 = (-P909);
                c9 = P540;
                c10 = P755;
            } else if (mu == 4) {
                c1 = (-P654);
                c2 = (-P142);
                c3 = P841;
                c4 = (-P959);
                c5 = P415;
                c6 = P755;
                c7 = (-P989);
                c8 = P540;
                c9 = P281;
                c10 = (-P909);
            } else if (mu == 5) {
                c1 = (-P959);
                c2 = P841;
                c3 = (-P654);
                c4 = P415;
                c5 = (-P142);
                c6 = P281;
                c7 = (-P540);
                c8 = P755;
                c9 = (-P909);
                c10 = P989;
            } else if (mu == 6) {
                c1 = (-P959);
                c2 = P841;
                c3 = (-P654);
                c4 = P415;
                c5 = (-P142);
                c6 = (-P281);
                c7 = P540;
                c8 = (-P755);
                c9 = P909;
                c10 = (-P989);
            } else if (mu == 7) {
                c1 = (-P654);
                c2 = (-P142);
                c3 = P841;
                c4 = (-P959);
                c5 = P415;
                c6 = (-P755);
                c7 = P989;
                c8 = (-P540);
                c9 = (-P281);
                c10 = P909;
            } else if (mu == 8) {
                c1 = (-P142);
                c2 = (-P959);
                c3 = P415;
                c4 = P841;
                c5 = (-P654);
                c6 = (-P989);
                c7 = P281;
                c8 = P909;
                c9 = (-P540);
                c10 = (-P755);
            } else if (mu == 9) {
                c1 = P415;
                c2 = (-P654);
                c3 = (-P959);
                c4 = (-P142);
                c5 = P841;
                c6 = (-P909);
                c7 = (-P755);
                c8 = P281;
                c9 = P989;
                c10 = P540;
            } else {
                c1 = P841;
                c2 = P415;
                c3 = (-P142);
                c4 = (-P654);
                c5 = (-P959);
                c6 = (-P540);
                c7 = (-P909);
                c8 = (-P989);
                c9 = (-P755);
                c10 = (-P281);
            }
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                for _jt in 0..nt {
                    t1r = (z[(j01) as usize] + z[(j10) as usize]);
                    t1i = (z[(j01 + 1) as usize] + z[(j10 + 1) as usize]);
                    t2r = (z[(j2) as usize] + z[(j9) as usize]);
                    t2i = (z[(j2 + 1) as usize] + z[(j9 + 1) as usize]);
                    t3r = (z[(j3) as usize] + z[(j8) as usize]);
                    t3i = (z[(j3 + 1) as usize] + z[(j8 + 1) as usize]);
                    t4r = (z[(j4) as usize] + z[(j7) as usize]);
                    t4i = (z[(j4 + 1) as usize] + z[(j7 + 1) as usize]);
                    t5r = (z[(j5) as usize] + z[(j6) as usize]);
                    t5i = (z[(j5 + 1) as usize] + z[(j6 + 1) as usize]);
                    t6r = (z[(j01) as usize] - z[(j10) as usize]);
                    t6i = (z[(j01 + 1) as usize] - z[(j10 + 1) as usize]);
                    t7r = (z[(j2) as usize] - z[(j9) as usize]);
                    t7i = (z[(j2 + 1) as usize] - z[(j9 + 1) as usize]);
                    t8r = (z[(j3) as usize] - z[(j8) as usize]);
                    t8i = (z[(j3 + 1) as usize] - z[(j8 + 1) as usize]);
                    t9r = (z[(j4) as usize] - z[(j7) as usize]);
                    t9i = (z[(j4 + 1) as usize] - z[(j7 + 1) as usize]);
                    t10r = (z[(j5) as usize] - z[(j6) as usize]);
                    t10i = (z[(j5 + 1) as usize] - z[(j6 + 1) as usize]);
                    t11r = (z[(j00) as usize].wide() - Double(0.5) * t5r.wide()).narrow();
                    t11i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t5i.wide()).narrow();
                    t12r = (t1r - t5r);
                    t12i = (t1i - t5i);
                    t13r = (t2r - t5r);
                    t13i = (t2i - t5i);
                    t14r = (t3r - t5r);
                    t14i = (t3i - t5i);
                    t15r = (t4r - t5r);
                    t15i = (t4i - t5i);
                    y1r = ((((t11r + (c1 * t12r)) + (c2 * t13r)) + (c3 * t14r)) + (c4 * t15r));
                    y1i = ((((t11i + (c1 * t12i)) + (c2 * t13i)) + (c3 * t14i)) + (c4 * t15i));
                    y2r = ((((t11r + (c2 * t12r)) + (c4 * t13r)) + (c5 * t14r)) + (c3 * t15r));
                    y2i = ((((t11i + (c2 * t12i)) + (c4 * t13i)) + (c5 * t14i)) + (c3 * t15i));
                    y3r = ((((t11r + (c3 * t12r)) + (c5 * t13r)) + (c2 * t14r)) + (c1 * t15r));
                    y3i = ((((t11i + (c3 * t12i)) + (c5 * t13i)) + (c2 * t14i)) + (c1 * t15i));
                    y4r = ((((t11r + (c4 * t12r)) + (c3 * t13r)) + (c1 * t14r)) + (c5 * t15r));
                    y4i = ((((t11i + (c4 * t12i)) + (c3 * t13i)) + (c1 * t14i)) + (c5 * t15i));
                    y5r = ((((t11r + (c5 * t12r)) + (c1 * t13r)) + (c4 * t14r)) + (c2 * t15r));
                    y5i = ((((t11i + (c5 * t12i)) + (c1 * t13i)) + (c4 * t14i)) + (c2 * t15i));
                    y6r = (((((c10 * t6r) - (c6 * t7r)) + (c9 * t8r)) - (c7 * t9r)) + (c8 * t10r));
                    y6i = (((((c10 * t6i) - (c6 * t7i)) + (c9 * t8i)) - (c7 * t9i)) + (c8 * t10i));
                    y7r = (((((c9 * t6r) - (c8 * t7r)) + (c6 * t8r)) + (c10 * t9r)) - (c7 * t10r));
                    y7i = (((((c9 * t6i) - (c8 * t7i)) + (c6 * t8i)) + (c10 * t9i)) - (c7 * t10i));
                    y8r = (((((c8 * t6r) - (c10 * t7r)) - (c7 * t8r)) + (c6 * t9r)) + (c9 * t10r));
                    y8i = (((((c8 * t6i) - (c10 * t7i)) - (c7 * t8i)) + (c6 * t9i)) + (c9 * t10i));
                    y9r = (((((c7 * t6r) + (c9 * t7r)) - (c10 * t8r)) - (c8 * t9r)) - (c6 * t10r));
                    y9i = (((((c7 * t6i) + (c9 * t7i)) - (c10 * t8i)) - (c8 * t9i)) - (c6 * t10i));
                    y10r = (((((c6 * t6r) + (c7 * t7r)) + (c8 * t8r)) + (c9 * t9r)) + (c10 * t10r));
                    y10i = (((((c6 * t6i) + (c7 * t7i)) + (c8 * t8i)) + (c9 * t9i)) + (c10 * t10i));
                    z[(j00) as usize] = (((((z[(j00) as usize] + t1r) + t2r) + t3r) + t4r) + t5r);
                    z[(j00 + 1) as usize] =
                        (((((z[(j00 + 1) as usize] + t1i) + t2i) + t3i) + t4i) + t5i);
                    z[(j01) as usize] = (y1r - y10i);
                    z[(j01 + 1) as usize] = (y1i + y10r);
                    z[(j2) as usize] = (y2r - y9i);
                    z[(j2 + 1) as usize] = (y2i + y9r);
                    z[(j3) as usize] = (y3r - y8i);
                    z[(j3 + 1) as usize] = (y3i + y8r);
                    z[(j4) as usize] = (y4r - y7i);
                    z[(j4 + 1) as usize] = (y4i + y7r);
                    z[(j5) as usize] = (y5r - y6i);
                    z[(j5 + 1) as usize] = (y5i + y6r);
                    z[(j6) as usize] = (y5r + y6i);
                    z[(j6 + 1) as usize] = (y5i - y6r);
                    z[(j7) as usize] = (y4r + y7i);
                    z[(j7 + 1) as usize] = (y4i - y7r);
                    z[(j8) as usize] = (y3r + y8i);
                    z[(j8 + 1) as usize] = (y3i - y8r);
                    z[(j9) as usize] = (y2r + y9i);
                    z[(j9 + 1) as usize] = (y2i - y9r);
                    z[(j10) as usize] = (y1r + y10i);
                    z[(j10 + 1) as usize] = (y1i - y10r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                }
                it = (i10 + istep);
                i10 = (i9 + istep);
                i9 = (i8 + istep);
                i8 = (i7 + istep);
                i7 = (i6 + istep);
                i6 = (i5 + istep);
                i5 = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i11 = (i10 + iinc);
        if (i11 >= imax) {
            i11 = (i11 - imax);
        }
        i12 = (i11 + iinc);
        if (i12 >= imax) {
            i12 = (i12 - imax);
        }
        /* if factor is 13 */
        if (ifac == 13) {
            if (mu == 1) {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = (-P354);
                c5 = (-P748);
                c6 = (-P970);
                c7 = P464;
                c8 = P822;
                c9 = P992;
                c10 = P935;
                c11 = P663;
                c12 = P239;
            } else if (mu == 2) {
                c1 = P568;
                c2 = (-P354);
                c3 = (-P970);
                c4 = (-P748);
                c5 = P120;
                c6 = P885;
                c7 = P822;
                c8 = P935;
                c9 = P239;
                c10 = (-P663);
                c11 = (-P992);
                c12 = (-P464);
            } else if (mu == 3) {
                c1 = P120;
                c2 = (-P970);
                c3 = (-P354);
                c4 = P885;
                c5 = P568;
                c6 = (-P748);
                c7 = P992;
                c8 = P239;
                c9 = (-P935);
                c10 = (-P464);
                c11 = P822;
                c12 = P663;
            } else if (mu == 4) {
                c1 = (-P354);
                c2 = (-P748);
                c3 = P885;
                c4 = P120;
                c5 = (-P970);
                c6 = P568;
                c7 = P935;
                c8 = (-P663);
                c9 = (-P464);
                c10 = P992;
                c11 = (-P239);
                c12 = (-P822);
            } else if (mu == 5) {
                c1 = (-P748);
                c2 = P120;
                c3 = P568;
                c4 = (-P970);
                c5 = P885;
                c6 = (-P354);
                c7 = P663;
                c8 = (-P992);
                c9 = P822;
                c10 = (-P239);
                c11 = (-P464);
                c12 = P935;
            } else if (mu == 6) {
                c1 = (-P970);
                c2 = P885;
                c3 = (-P748);
                c4 = P568;
                c5 = (-P354);
                c6 = P120;
                c7 = P239;
                c8 = (-P464);
                c9 = P663;
                c10 = (-P822);
                c11 = P935;
                c12 = (-P992);
            } else if (mu == 7) {
                c1 = (-P970);
                c2 = P885;
                c3 = (-P748);
                c4 = P568;
                c5 = (-P354);
                c6 = P120;
                c7 = (-P239);
                c8 = P464;
                c9 = (-P663);
                c10 = P822;
                c11 = (-P935);
                c12 = P992;
            } else if (mu == 8) {
                c1 = (-P748);
                c2 = P120;
                c3 = P568;
                c4 = (-P970);
                c5 = P885;
                c6 = (-P354);
                c7 = (-P663);
                c8 = P992;
                c9 = (-P822);
                c10 = P239;
                c11 = P464;
                c12 = (-P935);
            } else if (mu == 9) {
                c1 = (-P354);
                c2 = (-P748);
                c3 = P885;
                c4 = P120;
                c5 = (-P970);
                c6 = P568;
                c7 = (-P935);
                c8 = P663;
                c9 = P464;
                c10 = (-P992);
                c11 = P239;
                c12 = P822;
            } else if (mu == 10) {
                c1 = P120;
                c2 = (-P970);
                c3 = (-P354);
                c4 = P885;
                c5 = P568;
                c6 = (-P748);
                c7 = (-P992);
                c8 = (-P239);
                c9 = P935;
                c10 = P464;
                c11 = (-P822);
                c12 = (-P663);
            } else if (mu == 11) {
                c1 = P568;
                c2 = (-P354);
                c3 = (-P970);
                c4 = (-P748);
                c5 = P120;
                c6 = P885;
                c7 = (-P822);
                c8 = (-P935);
                c9 = (-P239);
                c10 = P663;
                c11 = P992;
                c12 = P464;
            } else {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = (-P354);
                c5 = (-P748);
                c6 = (-P970);
                c7 = (-P464);
                c8 = (-P822);
                c9 = (-P992);
                c10 = (-P935);
                c11 = (-P663);
                c12 = (-P239);
            }
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                j11 = i11;
                j12 = i12;
                for _jt in 0..nt {
                    t1r = (z[(j01) as usize] + z[(j12) as usize]);
                    t1i = (z[(j01 + 1) as usize] + z[(j12 + 1) as usize]);
                    t2r = (z[(j2) as usize] + z[(j11) as usize]);
                    t2i = (z[(j2 + 1) as usize] + z[(j11 + 1) as usize]);
                    t3r = (z[(j3) as usize] + z[(j10) as usize]);
                    t3i = (z[(j3 + 1) as usize] + z[(j10 + 1) as usize]);
                    t4r = (z[(j4) as usize] + z[(j9) as usize]);
                    t4i = (z[(j4 + 1) as usize] + z[(j9 + 1) as usize]);
                    t5r = (z[(j5) as usize] + z[(j8) as usize]);
                    t5i = (z[(j5 + 1) as usize] + z[(j8 + 1) as usize]);
                    t6r = (z[(j6) as usize] + z[(j7) as usize]);
                    t6i = (z[(j6 + 1) as usize] + z[(j7 + 1) as usize]);
                    t7r = (z[(j01) as usize] - z[(j12) as usize]);
                    t7i = (z[(j01 + 1) as usize] - z[(j12 + 1) as usize]);
                    t8r = (z[(j2) as usize] - z[(j11) as usize]);
                    t8i = (z[(j2 + 1) as usize] - z[(j11 + 1) as usize]);
                    t9r = (z[(j3) as usize] - z[(j10) as usize]);
                    t9i = (z[(j3 + 1) as usize] - z[(j10 + 1) as usize]);
                    t10r = (z[(j4) as usize] - z[(j9) as usize]);
                    t10i = (z[(j4 + 1) as usize] - z[(j9 + 1) as usize]);
                    t11r = (z[(j5) as usize] - z[(j8) as usize]);
                    t11i = (z[(j5 + 1) as usize] - z[(j8 + 1) as usize]);
                    t12r = (z[(j6) as usize] - z[(j7) as usize]);
                    t12i = (z[(j6 + 1) as usize] - z[(j7 + 1) as usize]);
                    t13r = (z[(j00) as usize].wide() - Double(0.5) * t6r.wide()).narrow();
                    t13i = (z[(j00 + 1) as usize].wide() - Double(0.5) * t6i.wide()).narrow();
                    t14r = (t1r - t6r);
                    t14i = (t1i - t6i);
                    t15r = (t2r - t6r);
                    t15i = (t2i - t6i);
                    t16r = (t3r - t6r);
                    t16i = (t3i - t6i);
                    t17r = (t4r - t6r);
                    t17i = (t4i - t6i);
                    t18r = (t5r - t6r);
                    t18i = (t5i - t6i);
                    y1r = (((((t13r + (c1 * t14r)) + (c2 * t15r)) + (c3 * t16r)) + (c4 * t17r))
                        + (c5 * t18r));
                    y1i = (((((t13i + (c1 * t14i)) + (c2 * t15i)) + (c3 * t16i)) + (c4 * t17i))
                        + (c5 * t18i));
                    y2r = (((((t13r + (c2 * t14r)) + (c4 * t15r)) + (c6 * t16r)) + (c5 * t17r))
                        + (c3 * t18r));
                    y2i = (((((t13i + (c2 * t14i)) + (c4 * t15i)) + (c6 * t16i)) + (c5 * t17i))
                        + (c3 * t18i));
                    y3r = (((((t13r + (c3 * t14r)) + (c6 * t15r)) + (c4 * t16r)) + (c1 * t17r))
                        + (c2 * t18r));
                    y3i = (((((t13i + (c3 * t14i)) + (c6 * t15i)) + (c4 * t16i)) + (c1 * t17i))
                        + (c2 * t18i));
                    y4r = (((((t13r + (c4 * t14r)) + (c5 * t15r)) + (c1 * t16r)) + (c3 * t17r))
                        + (c6 * t18r));
                    y4i = (((((t13i + (c4 * t14i)) + (c5 * t15i)) + (c1 * t16i)) + (c3 * t17i))
                        + (c6 * t18i));
                    y5r = (((((t13r + (c5 * t14r)) + (c3 * t15r)) + (c2 * t16r)) + (c6 * t17r))
                        + (c1 * t18r));
                    y5i = (((((t13i + (c5 * t14i)) + (c3 * t15i)) + (c2 * t16i)) + (c6 * t17i))
                        + (c1 * t18i));
                    y6r = (((((t13r + (c6 * t14r)) + (c1 * t15r)) + (c5 * t16r)) + (c2 * t17r))
                        + (c4 * t18r));
                    y6i = (((((t13i + (c6 * t14i)) + (c1 * t15i)) + (c5 * t16i)) + (c2 * t17i))
                        + (c4 * t18i));
                    y7r = ((((((c12 * t7r) - (c7 * t8r)) + (c11 * t9r)) - (c8 * t10r))
                        + (c10 * t11r))
                        - (c9 * t12r));
                    y7i = ((((((c12 * t7i) - (c7 * t8i)) + (c11 * t9i)) - (c8 * t10i))
                        + (c10 * t11i))
                        - (c9 * t12i));
                    y8r = ((((((c11 * t7r) - (c9 * t8r)) + (c8 * t9r)) - (c12 * t10r))
                        - (c7 * t11r))
                        + (c10 * t12r));
                    y8i = ((((((c11 * t7i) - (c9 * t8i)) + (c8 * t9i)) - (c12 * t10i))
                        - (c7 * t11i))
                        + (c10 * t12i));
                    y9r = ((((((c10 * t7r) - (c11 * t8r)) - (c7 * t9r)) + (c9 * t10r))
                        - (c12 * t11r))
                        - (c8 * t12r));
                    y9i = ((((((c10 * t7i) - (c11 * t8i)) - (c7 * t9i)) + (c9 * t10i))
                        - (c12 * t11i))
                        - (c8 * t12i));
                    y10r = ((((((c9 * t7r) + (c12 * t8r)) - (c10 * t9r)) - (c7 * t10r))
                        + (c8 * t11r))
                        + (c11 * t12r));
                    y10i = ((((((c9 * t7i) + (c12 * t8i)) - (c10 * t9i)) - (c7 * t10i))
                        + (c8 * t11i))
                        + (c11 * t12i));
                    y11r = ((((((c8 * t7r) + (c10 * t8r)) + (c12 * t9r)) - (c11 * t10r))
                        - (c9 * t11r))
                        - (c7 * t12r));
                    y11i = ((((((c8 * t7i) + (c10 * t8i)) + (c12 * t9i)) - (c11 * t10i))
                        - (c9 * t11i))
                        - (c7 * t12i));
                    y12r = ((((((c7 * t7r) + (c8 * t8r)) + (c9 * t9r)) + (c10 * t10r))
                        + (c11 * t11r))
                        + (c12 * t12r));
                    y12i = ((((((c7 * t7i) + (c8 * t8i)) + (c9 * t9i)) + (c10 * t10i))
                        + (c11 * t11i))
                        + (c12 * t12i));
                    z[(j00) as usize] =
                        ((((((z[(j00) as usize] + t1r) + t2r) + t3r) + t4r) + t5r) + t6r);
                    z[(j00 + 1) as usize] =
                        ((((((z[(j00 + 1) as usize] + t1i) + t2i) + t3i) + t4i) + t5i) + t6i);
                    z[(j01) as usize] = (y1r - y12i);
                    z[(j01 + 1) as usize] = (y1i + y12r);
                    z[(j2) as usize] = (y2r - y11i);
                    z[(j2 + 1) as usize] = (y2i + y11r);
                    z[(j3) as usize] = (y3r - y10i);
                    z[(j3 + 1) as usize] = (y3i + y10r);
                    z[(j4) as usize] = (y4r - y9i);
                    z[(j4 + 1) as usize] = (y4i + y9r);
                    z[(j5) as usize] = (y5r - y8i);
                    z[(j5 + 1) as usize] = (y5i + y8r);
                    z[(j6) as usize] = (y6r - y7i);
                    z[(j6 + 1) as usize] = (y6i + y7r);
                    z[(j7) as usize] = (y6r + y7i);
                    z[(j7 + 1) as usize] = (y6i - y7r);
                    z[(j8) as usize] = (y5r + y8i);
                    z[(j8 + 1) as usize] = (y5i - y8r);
                    z[(j9) as usize] = (y4r + y9i);
                    z[(j9 + 1) as usize] = (y4i - y9r);
                    z[(j10) as usize] = (y3r + y10i);
                    z[(j10 + 1) as usize] = (y3i - y10r);
                    z[(j11) as usize] = (y2r + y11i);
                    z[(j11 + 1) as usize] = (y2i - y11r);
                    z[(j12) as usize] = (y1r + y12i);
                    z[(j12 + 1) as usize] = (y1i - y12r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                    j11 += jstep;
                    j12 += jstep;
                }
                it = (i12 + istep);
                i12 = (i11 + istep);
                i11 = (i10 + istep);
                i10 = (i9 + istep);
                i9 = (i8 + istep);
                i8 = (i7 + istep);
                i7 = (i6 + istep);
                i6 = (i5 + istep);
                i5 = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
        i13 = (i12 + iinc);
        if (i13 >= imax) {
            i13 = (i13 - imax);
        }
        i14 = (i13 + iinc);
        if (i14 >= imax) {
            i14 = (i14 - imax);
        }
        i15 = (i14 + iinc);
        if (i15 >= imax) {
            i15 = (i15 - imax);
        }
        /* if factor is 16 */
        if (ifac == 16) {
            if (mu == 1) {
                c1 = Double(1.0).narrow();
                c2 = P923;
                c3 = P382;
                c4 = P707;
            } else if (mu == 3) {
                c1 = (-Double(1.0)).narrow();
                c2 = P382;
                c3 = P923;
                c4 = (-P707);
            } else if (mu == 5) {
                c1 = Double(1.0).narrow();
                c2 = (-P382);
                c3 = P923;
                c4 = (-P707);
            } else if (mu == 7) {
                c1 = (-Double(1.0)).narrow();
                c2 = (-P923);
                c3 = P382;
                c4 = P707;
            } else if (mu == 9) {
                c1 = Double(1.0).narrow();
                c2 = (-P923);
                c3 = (-P382);
                c4 = P707;
            } else if (mu == 11) {
                c1 = (-Double(1.0)).narrow();
                c2 = (-P382);
                c3 = (-P923);
                c4 = (-P707);
            } else if (mu == 13) {
                c1 = Double(1.0).narrow();
                c2 = P382;
                c3 = (-P923);
                c4 = (-P707);
            } else {
                c1 = (-Double(1.0)).narrow();
                c2 = P923;
                c3 = (-P382);
                c4 = P707;
            }
            c5 = (c1 * c4);
            c6 = (c1 * c3);
            c7 = (c1 * c2);
            for _l in 0..m {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                j11 = i11;
                j12 = i12;
                j13 = i13;
                j14 = i14;
                j15 = i15;
                for _jt in 0..nt {
                    t1r = (z[(j00) as usize] + z[(j8) as usize]);
                    t1i = (z[(j00 + 1) as usize] + z[(j8 + 1) as usize]);
                    t2r = (z[(j4) as usize] + z[(j12) as usize]);
                    t2i = (z[(j4 + 1) as usize] + z[(j12 + 1) as usize]);
                    t3r = (z[(j00) as usize] - z[(j8) as usize]);
                    t3i = (z[(j00 + 1) as usize] - z[(j8 + 1) as usize]);
                    t4r = (c1 * (z[(j4) as usize] - z[(j12) as usize]));
                    t4i = (c1 * (z[(j4 + 1) as usize] - z[(j12 + 1) as usize]));
                    t5r = (t1r + t2r);
                    t5i = (t1i + t2i);
                    t6r = (t1r - t2r);
                    t6i = (t1i - t2i);
                    t7r = (z[(j01) as usize] + z[(j9) as usize]);
                    t7i = (z[(j01 + 1) as usize] + z[(j9 + 1) as usize]);
                    t8r = (z[(j5) as usize] + z[(j13) as usize]);
                    t8i = (z[(j5 + 1) as usize] + z[(j13 + 1) as usize]);
                    t9r = (z[(j01) as usize] - z[(j9) as usize]);
                    t9i = (z[(j01 + 1) as usize] - z[(j9 + 1) as usize]);
                    t10r = (z[(j5) as usize] - z[(j13) as usize]);
                    t10i = (z[(j5 + 1) as usize] - z[(j13 + 1) as usize]);
                    t11r = (t7r + t8r);
                    t11i = (t7i + t8i);
                    t12r = (t7r - t8r);
                    t12i = (t7i - t8i);
                    t13r = (z[(j2) as usize] + z[(j10) as usize]);
                    t13i = (z[(j2 + 1) as usize] + z[(j10 + 1) as usize]);
                    t14r = (z[(j6) as usize] + z[(j14) as usize]);
                    t14i = (z[(j6 + 1) as usize] + z[(j14 + 1) as usize]);
                    t15r = (z[(j2) as usize] - z[(j10) as usize]);
                    t15i = (z[(j2 + 1) as usize] - z[(j10 + 1) as usize]);
                    t16r = (z[(j6) as usize] - z[(j14) as usize]);
                    t16i = (z[(j6 + 1) as usize] - z[(j14 + 1) as usize]);
                    t17r = (t13r + t14r);
                    t17i = (t13i + t14i);
                    t18r = (c4 * (t15r - t16r));
                    t18i = (c4 * (t15i - t16i));
                    t19r = (c5 * (t15r + t16r));
                    t19i = (c5 * (t15i + t16i));
                    t20r = (c1 * (t13r - t14r));
                    t20i = (c1 * (t13i - t14i));
                    t21r = (z[(j3) as usize] + z[(j11) as usize]);
                    t21i = (z[(j3 + 1) as usize] + z[(j11 + 1) as usize]);
                    t22r = (z[(j7) as usize] + z[(j15) as usize]);
                    t22i = (z[(j7 + 1) as usize] + z[(j15 + 1) as usize]);
                    t23r = (z[(j3) as usize] - z[(j11) as usize]);
                    t23i = (z[(j3 + 1) as usize] - z[(j11 + 1) as usize]);
                    t24r = (z[(j7) as usize] - z[(j15) as usize]);
                    t24i = (z[(j7 + 1) as usize] - z[(j15 + 1) as usize]);
                    t25r = (t21r + t22r);
                    t25i = (t21i + t22i);
                    t26r = (t21r - t22r);
                    t26i = (t21i - t22i);
                    t27r = (t9r + t24r);
                    t27i = (t9i + t24i);
                    t28r = (t10r + t23r);
                    t28i = (t10i + t23i);
                    t29r = (t9r - t24r);
                    t29i = (t9i - t24i);
                    t30r = (t10r - t23r);
                    t30i = (t10i - t23i);
                    t31r = (t5r + t17r);
                    t31i = (t5i + t17i);
                    t32r = (t11r + t25r);
                    t32i = (t11i + t25i);
                    t33r = (t3r + t18r);
                    t33i = (t3i + t18i);
                    t34r = ((c2 * t29r) - (c6 * t30r));
                    t34i = ((c2 * t29i) - (c6 * t30i));
                    t35r = (t3r - t18r);
                    t35i = (t3i - t18i);
                    t36r = ((c7 * t27r) - (c3 * t28r));
                    t36i = ((c7 * t27i) - (c3 * t28i));
                    t37r = (t4r + t19r);
                    t37i = (t4i + t19i);
                    t38r = ((c3 * t27r) + (c7 * t28r));
                    t38i = ((c3 * t27i) + (c7 * t28i));
                    t39r = (t4r - t19r);
                    t39i = (t4i - t19i);
                    t40r = ((c6 * t29r) + (c2 * t30r));
                    t40i = ((c6 * t29i) + (c2 * t30i));
                    t41r = (c4 * (t12r - t26r));
                    t41i = (c4 * (t12i - t26i));
                    t42r = (c5 * (t12r + t26r));
                    t42i = (c5 * (t12i + t26i));
                    y1r = (t33r + t34r);
                    y1i = (t33i + t34i);
                    y2r = (t6r + t41r);
                    y2i = (t6i + t41i);
                    y3r = (t35r + t40r);
                    y3i = (t35i + t40i);
                    y4r = (t5r - t17r);
                    y4i = (t5i - t17i);
                    y5r = (t35r - t40r);
                    y5i = (t35i - t40i);
                    y6r = (t6r - t41r);
                    y6i = (t6i - t41i);
                    y7r = (t33r - t34r);
                    y7i = (t33i - t34i);
                    y9r = (t38r - t37r);
                    y9i = (t38i - t37i);
                    y10r = (t42r - t20r);
                    y10i = (t42i - t20i);
                    y11r = (t36r + t39r);
                    y11i = (t36i + t39i);
                    y12r = (c1 * (t11r - t25r));
                    y12i = (c1 * (t11i - t25i));
                    y13r = (t36r - t39r);
                    y13i = (t36i - t39i);
                    y14r = (t42r + t20r);
                    y14i = (t42i + t20i);
                    y15r = (t38r + t37r);
                    y15i = (t38i + t37i);
                    z[(j00) as usize] = (t31r + t32r);
                    z[(j00 + 1) as usize] = (t31i + t32i);
                    z[(j01) as usize] = (y1r - y15i);
                    z[(j01 + 1) as usize] = (y1i + y15r);
                    z[(j2) as usize] = (y2r - y14i);
                    z[(j2 + 1) as usize] = (y2i + y14r);
                    z[(j3) as usize] = (y3r - y13i);
                    z[(j3 + 1) as usize] = (y3i + y13r);
                    z[(j4) as usize] = (y4r - y12i);
                    z[(j4 + 1) as usize] = (y4i + y12r);
                    z[(j5) as usize] = (y5r - y11i);
                    z[(j5 + 1) as usize] = (y5i + y11r);
                    z[(j6) as usize] = (y6r - y10i);
                    z[(j6 + 1) as usize] = (y6i + y10r);
                    z[(j7) as usize] = (y7r - y9i);
                    z[(j7 + 1) as usize] = (y7i + y9r);
                    z[(j8) as usize] = (t31r - t32r);
                    z[(j8 + 1) as usize] = (t31i - t32i);
                    z[(j9) as usize] = (y7r + y9i);
                    z[(j9 + 1) as usize] = (y7i - y9r);
                    z[(j10) as usize] = (y6r + y10i);
                    z[(j10 + 1) as usize] = (y6i - y10r);
                    z[(j11) as usize] = (y5r + y11i);
                    z[(j11 + 1) as usize] = (y5i - y11r);
                    z[(j12) as usize] = (y4r + y12i);
                    z[(j12 + 1) as usize] = (y4i - y12r);
                    z[(j13) as usize] = (y3r + y13i);
                    z[(j13 + 1) as usize] = (y3i - y13r);
                    z[(j14) as usize] = (y2r + y14i);
                    z[(j14 + 1) as usize] = (y2i - y14r);
                    z[(j15) as usize] = (y1r + y15i);
                    z[(j15 + 1) as usize] = (y1i - y15r);
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                    j11 += jstep;
                    j12 += jstep;
                    j13 += jstep;
                    j14 += jstep;
                    j15 += jstep;
                }
                it = (i15 + istep);
                i15 = (i14 + istep);
                i14 = (i13 + istep);
                i13 = (i12 + istep);
                i12 = (i11 + istep);
                i11 = (i10 + istep);
                i10 = (i9 + istep);
                i9 = (i8 + istep);
                i8 = (i7 + istep);
                i7 = (i6 + istep);
                i6 = (i5 + istep);
                i5 = (i4 + istep);
                i4 = (i3 + istep);
                i3 = (i2 + istep);
                i2 = (i1 + istep);
                i1 = (i0 + istep);
                i0 = it;
            }
            continue;
        }
    }
}
