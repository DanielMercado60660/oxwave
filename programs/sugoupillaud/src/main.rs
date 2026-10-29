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

//! SUGOUPILLAUD - calculate the 1D impulse response of a non-absorbing
//! Goupillaud medium. Port of su/main/synthetics_waveforms_testpatterns/
//! sugoupillaud.c (Revision 1.7).
//! Credits: CWP: Albena Mateeva, May 2000, a summer project at Western
//! Geophysical. Polynomial multiplication and division are based on
//! Robinson's Fortran subroutines, Multichannel Time Series Analysis with
//! Digital Computer Programs (1983), Chapter 1.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "  \t\t\t\t\t\t\t\t\t",
    " SUGOUPILLAUD - calculate 1D impulse response of\t \t\t",
    "     non-absorbing Goupillaud medium\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sugoupillaud < stdin > stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    "l=1source layer number; 1 <= l <= tr.ns\t\t\t\t",
    "\t\t\tSource is located at the top of layer l.\t",
    "\tk=1\t\treceiver layer number; 1 <= k\t\t\t",
    "Receiver is located at the top of layer k.\t\t\t\t",
    "tmax  number of output time-samples; default:\t\t\t\t",
    "tmax=NINT((2*tr.ns-(l-1)-(k-1))/2)  if k < tr.ns\t\t\t",
    "\t\t\ttmax=k\t\t\t\tif k >=tr.ns\t",
    "pV=1  flag for vector field seismogram\t\t\t\t\t",
    "\t(displacement, velocity, acceleration);\t\t\t\t",
    "=-1 for pressure seismogram.\t\t\t\t\t\t",
    "verbose=0  silent operation, =1 list warnings\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Input: Reflection coefficient series:\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t impedance[i]-impedance[i+1]\t\t\t\t\t",
    " r[i] = ----------------------------- \t\t\t\t\t",
    "\t impedance[i]+impedance[i+1]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tr[0]= surface refl. coef. (as seen from above)\t\t\t",
    "r[n]= refl. coef. of the deepest interface\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Input file is to be in SU format, i.e., binary floats with a SU header.",
    "\t\t\t\t\t\t\t\t\t",
    " Remarks:\t\t\t\t\t\t\t\t",
    " 1. For vector fields, a buried source produces a spike of amplitude 1\t",
    " propagating downwards and a spike of amplitude -1 propagating upwards.",
    " A buried pressure source produces spikes of amplitude 1 both in the up",
    " and downward directions.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " A surface source induces only a downgoing spike of amplitude 1 at the\t",
    " top of the first layer (both for vector and pressure fields).\t\t",
    " 2. The sampling interval dt in the header of the input reflectivity file",
    " is interpreted as a two-way traveltime thicknes of the layers. The sampling",
    " interval of the output seismogram is the same as that of the input file.",
];

#[allow(non_snake_case, unused_assignments)]
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get subsurface model.
    if gettr(&mut tr) == 0 {
        err!("Can't get reflectivity");
    }
    let n = tr.ns as i32 - 1;
    let (mut k, mut l, mut tmax, mut pV, mut verbose) = (1, 1, 0, 1, 0);
    if getparint("k", &mut k) == 0 {
        k = 1;
    }
    if getparint("l", &mut l) == 0 {
        l = 1;
    }
    if getparint("tmax", &mut tmax) == 0 {
        tmax = if k < n + 1 {
            (2 * n + 2)
                .wrapping_sub(l.wrapping_sub(1))
                .wrapping_sub(k.wrapping_sub(1))
                / 2
        } else {
            k
        };
    }
    if getparint("pV", &mut pV) == 0 {
        pV = 1;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    if n <= 0 {
        err!("The number of subsurface interfaces n must be >=1!");
    }
    if k < 1 {
        err!("Receiver layer k must be >=1 (k=1 corresponds to surface seismogram).");
    }
    if l < 1 {
        err!("Source layer l must be >= 1 (l=1 corresponds to a surface source).");
    }
    if l > n + 1 {
        err!("The current version of the program requires l<=n+1.");
    }
    if tmax < 0 {
        err!("The number of the output time samples tmax cannot be negative.");
    }
    if !(pV == 1 || pV == -1) {
        err!("The field-type flag pV should be either 1 or -1.");
    }

    if verbose != 0 {
        warn!("Number of layers n=%d", n);
        warn!("Source layer l=%d", l);
        warn!("Receiver layer k=%d", k);
        warn!("Output time samples tmax=%d", tmax);
        warn!("Field type flag pV=%d", pV);
    }

    // Helpers return their temporary output arrays, allowing the same in-place
    // calls as C without aliasing mutable Rust references.
    let mut r = vec![0.0f64; (n + 1) as usize];
    let mut x = vec![0.0f64; (2 * tmax) as usize];
    let mut u1 = vec![0.0f64; (2 * tmax + n) as usize];
    let mut d1 = vec![0.0f64; (2 * tmax + n) as usize];
    let mut u1b = vec![0.0f64; (2 * tmax + n) as usize];
    let mut uk = vec![0.0f64; (2 * tmax) as usize];
    let mut dk = vec![0.0f64; (2 * tmax) as usize];
    let mut storx = vec![0.0f64; (2 * tmax) as usize];
    let mut Q = vec![0.0f64; (n + 1) as usize];
    let mut P = vec![0.0f64; (n + 1) as usize];
    let mut qr = vec![0.0f64; (n + 1) as usize];
    let mut pr = vec![0.0f64; (n + 1) as usize];
    let mut cn = vec![0.0f64; (n + 1) as usize];
    let mut t1 = vec![0.0f64; (n + 1) as usize];
    let mut t2 = vec![0.0f64; (n + 1) as usize];
    let mut d0 = vec![0.0f64; (n + 1) as usize];
    let mut td1 = vec![0.0f64; (2 * (n + 1)) as usize];
    let mut mix1 = vec![0.0f64; (2 * tmax + n) as usize];
    let mut mix2 = vec![0.0f64; (2 * tmax + n) as usize];
    let mut mix3 = vec![0.0f64; (2 * (tmax + n)) as usize];
    let mut mix4 = vec![0.0f64; (2 * (tmax + n)) as usize];
    let mut mix = vec![0.0f64; (2 * (tmax + n)) as usize];
    let mut ak: f64;
    let uni = [1.0f64];

    for i in 0..=n {
        r[i as usize] = tr.data[i as usize] as f64;
        if r[i as usize] > 1.0 || r[i as usize] < -1.0 {
            err!("Invalid reflection coefficient encountered.");
        }
    }
    for i in 0..=n {
        r[i as usize] *= pV as f64;
    }

    // Begin computation. The n==0 branch in C is unreachable after validation.
    if n == 0 {
        x.fill(0.0);
        x[0] = 1.0;
        x = pshift(&x, 2 * tmax, k - 1, 2 * tmax);
    } else {
        // Upgoing field in layer 1, sampled in two-way layer times (Z).
        (Q, P) = recurs(n, n, &r);
        cn = pcomb(r[0], -1.0, &Q, &P, n + 1, n + 1, n + 1);
        u1 = pratio(&Q, &cn, n + 1, n + 1, 2 * tmax + n, 2 * tmax + n);

        if l == 1 {
            // I. Surface source.
            if k > n {
                // I.1. Receiver in the homogeneous lower half-space.
                ak = 1.0;
                for i in 1..=n {
                    ak *= 1.0 + r[i as usize];
                }
                ak *= -1.0;
                dk = pratio(&uni, &cn, 1, n + 1, 2 * tmax, 2 * tmax);
                dk = upsample(&dk, 2 * tmax, 2 * tmax);
                dk = pshift(&dk, 2 * tmax, n, 2 * tmax);
                for i in 0..2 * tmax {
                    dk[i as usize] *= ak;
                }
                x = pshift(&dk, 2 * tmax, k - n - 1, 2 * tmax);
            } else {
                // I.2. Receiver in the layers. d1 = -Pn/cn, Z-sampled.
                d1 = pratio(&P, &cn, n + 1, n + 1, 2 * tmax + n, 2 * tmax + n);
                for i in 0..2 * tmax + n {
                    d1[i as usize] *= -1.0;
                }
                if k == 1 {
                    // I.2.1. Surface seismogram, resampled from Z to z.
                    x = pcomb(1.0, 1.0, &u1, &d1, 2 * tmax + n, 2 * tmax + n, 2 * tmax);
                    x = upsample(&x, 2 * tmax, 2 * tmax);
                } else {
                    // I.2.2. Propagate the surface fields to layer k.
                    ak = 1.0;
                    for i in 1..k {
                        ak *= 1.0 - r[i as usize];
                    }
                    (Q, P) = recurs(k - 1, n, &r);
                    mix1 = prod(&Q, &d1, n + 1, 2 * tmax + n, 2 * tmax + n);
                    mix2 = prod(&P, &u1, n + 1, 2 * tmax + n, 2 * tmax + n);
                    mix1 = pcomb(
                        1.0,
                        1.0,
                        &mix1,
                        &mix2,
                        2 * tmax + n,
                        2 * tmax + n,
                        2 * tmax + n,
                    );
                    mix = upsample(&mix1, 2 * tmax + n, 2 * (tmax + n));
                    mix1 = pshift(&mix, 2 * (tmax + n), -(k - 1), 2 * tmax + n);
                    for i in 0..2 * tmax {
                        uk[i as usize] = mix1[i as usize] / ak;
                    }

                    qr = rev(&Q, n + 1);
                    qr = pshift(&qr, n + 1, k - 1 - porder(&Q, n + 1), n + 1);
                    pr = rev(&P, n + 1);
                    pr = pshift(&pr, n + 1, k - 1 - porder(&P, n + 1), n + 1);
                    mix1 = prod(&pr, &d1, n + 1, 2 * tmax + n, 2 * tmax + n);
                    mix2 = prod(&qr, &u1, n + 1, 2 * tmax + n, 2 * tmax + n);
                    mix1 = pcomb(
                        1.0,
                        1.0,
                        &mix1,
                        &mix2,
                        2 * tmax + n,
                        2 * tmax + n,
                        2 * tmax + n,
                    );
                    mix1 = upsample(&mix1, 2 * tmax + n, 2 * tmax + n);
                    mix1 = pshift(&mix1, 2 * tmax + n, -(k - 1), 2 * tmax + n);
                    for i in 0..2 * tmax {
                        dk[i as usize] = mix1[i as usize] / ak;
                    }
                    x = pcomb(1.0, 1.0, &uk, &dk, 2 * tmax, 2 * tmax, 2 * tmax);
                }
            }
        } else {
            // II. Buried source. d0 is Z-sampled; postpone its -(l-1)z shift.
            (Q, P) = recurs(l - 1, n, &r);
            t1 = pcomb(1.0, -r[0], &P, &Q, n + 1, n + 1, n + 1);
            qr = rev(&Q, n + 1);
            qr = pshift(&qr, n + 1, l - 1 - porder(&Q, n + 1), n + 1);
            pr = rev(&P, n + 1);
            pr = pshift(&pr, n + 1, l - 1 - porder(&P, n + 1), n + 1);
            t2 = pcomb(
                -(pV as f64) * 1.0,
                pV as f64 * r[0],
                &qr,
                &pr,
                n + 1,
                n + 1,
                n + 1,
            );
            d0 = pcomb(1.0, 1.0, &t1, &t2, n + 1, n + 1, n + 1);
            ak = 1.0;
            for i in 1..l {
                ak *= 1.0 + r[i as usize];
            }
            for i in 0..n + 1 {
                d0[i as usize] *= 1.0 / ak;
            }
            mix1 = prod(&u1, &d0, 2 * tmax + n, n + 1, 2 * tmax + n);
            mix = upsample(&mix1, 2 * tmax + n, 2 * (tmax + n));
            mix1 = pshift(&mix, 2 * (tmax + n), -(l - 1), 2 * tmax + n);

            if k >= l {
                // II.1. Receiver below the buried source.
                d1 = pcomb(1.0, -r[0], &uni, &u1, 1, 2 * tmax + n, 2 * tmax + n);
                let (kk, delay) = if k <= n + 1 {
                    (k, 0)
                } else {
                    (n + 1, k - n - 1)
                };
                ak = 1.0;
                for i in 1..kk {
                    ak *= 1.0 - r[i as usize];
                }
                (Q, P) = recurs(kk - 1, n, &r);
                if k <= n {
                    mix3 = prod(&Q, &d1, n + 1, 2 * tmax + n, 2 * (tmax + n));
                    mix4 = prod(&P, &u1, n + 1, 2 * tmax + n, 2 * (tmax + n));
                    mix3 = pcomb(
                        1.0 / ak,
                        1.0 / ak,
                        &mix3,
                        &mix4,
                        2 * (tmax + n),
                        2 * (tmax + n),
                        2 * (tmax + n),
                    );
                    mix3 = prod(&d0, &mix3, n + 1, 2 * (tmax + n), 2 * (tmax + n));
                    mix3 = upsample(&mix3, 2 * (tmax + n), 2 * (tmax + n));
                    uk = pshift(&mix3, 2 * (tmax + n), -(l - 1 + kk - 1), 2 * tmax);
                } else {
                    uk.fill(0.0);
                }

                qr = rev(&Q, n + 1);
                qr = pshift(&qr, n + 1, kk - 1 - porder(&Q, n + 1), n + 1);
                pr = rev(&P, n + 1);
                pr = pshift(&pr, n + 1, kk - 1 - porder(&P, n + 1), n + 1);
                mix3 = prod(&pr, &d1, n + 1, 2 * tmax + n, 2 * (tmax + n));
                mix4 = prod(&qr, &u1, n + 1, 2 * tmax + n, 2 * (tmax + n));
                mix3 = pcomb(
                    1.0 / ak,
                    1.0 / ak,
                    &mix3,
                    &mix4,
                    2 * (tmax + n),
                    2 * (tmax + n),
                    2 * (tmax + n),
                );
                mix3 = upsample(&mix3, 2 * (tmax + n), 2 * (tmax + n));
                mix3 = pshift(&mix3, 2 * (tmax + n), -(l - 1 + kk - 1), 2 * (tmax + n));
                mix3 = pshift(&mix3, 2 * (tmax + n), delay, 2 * (tmax + n));
                for i in 0..2 * tmax {
                    dk[i as usize] = mix3[i as usize];
                }
                x = pcomb(1.0, 1.0, &uk, &dk, 2 * tmax, 2 * tmax, 2 * tmax);
                if k == l {
                    for i in 0..2 * tmax {
                        storx[i as usize] = x[i as usize];
                    }
                }
            }
            if k <= l {
                // II.2. Receiver above the buried source.
                (Q, P) = recurs(l - 1, n, &r);
                pr = rev(&P, n + 1);
                pr = pshift(&pr, n + 1, l - 1 - porder(&P, n + 1), n + 1);
                t1 = pcomb(1.0, -(pV as f64) * 1.0, &Q, &pr, n + 1, n + 1, n + 1);
                td1 = upsample(&t1, n + 1, 2 * n + 2);
                td1 = pshift(&td1, 2 * (n + 1), -(l - 1), 2 * (n + 1));
                ak = 1.0;
                for i in 1..l {
                    ak *= 1.0 + r[i as usize];
                }
                for i in 0..2 * (n + 1) {
                    td1[i as usize] *= 1.0 / ak;
                }
                u1b = pcomb(
                    1.0,
                    1.0,
                    &mix1,
                    &td1,
                    2 * tmax + n,
                    2 * (n + 1),
                    2 * tmax + n,
                );
                if k == 1 {
                    for i in 0..2 * tmax {
                        x[i as usize] = u1b[i as usize] * (1.0 - r[0]);
                    }
                } else {
                    // Propagate the fields in the first layer to layer k.
                    (Q, P) = recurs(k - 1, n, &r);
                    t1 = pcomb(1.0, -r[0], &P, &Q, n + 1, n + 1, n + 1);
                    qr = rev(&Q, n + 1);
                    qr = pshift(&qr, n + 1, k - 1 - porder(&Q, n + 1), n + 1);
                    pr = rev(&P, n + 1);
                    pr = pshift(&pr, n + 1, k - 1 - porder(&P, n + 1), n + 1);
                    t2 = pcomb(1.0, -r[0], &qr, &pr, n + 1, n + 1, n + 1);
                    t1 = pcomb(1.0, 1.0, &t1, &t2, n + 1, n + 1, n + 1);
                    td1 = upsample(&t1, n + 1, 2 * n + 2);
                    mix = prod(&td1, &u1b, 2 * n + 2, 2 * tmax + n, 2 * (tmax + n));
                    mix2 = pshift(&mix, 2 * (tmax + n), -(k - 1), 2 * tmax + n);
                    ak = 1.0;
                    for i in 1..k {
                        ak *= 1.0 - r[i as usize];
                    }
                    for i in 0..2 * tmax {
                        x[i as usize] = mix2[i as usize] / ak;
                    }
                    // Average both fields for a colocated buried source/receiver.
                    if k == l {
                        x = pcomb(0.5, 0.5, &storx, &x, 2 * tmax, 2 * tmax, 2 * tmax);
                    }
                }
            }
        }
    }

    // Output synthetic seismogram, resampled to Z.
    let i = (k - l) / 2;
    if 2 * i == k - l {
        for i in 0..tmax {
            tr.data[i as usize] = x[(2 * i) as usize] as f32;
        }
        tr.delrt = 0;
    } else {
        for i in 0..tmax {
            tr.data[i as usize] = x[(2 * i + 1) as usize] as f32;
        }
        tr.delrt = (tr.dt / 2000) as i16;
    }
    tr.ns = tmax as u16;
    tr.trid = 1;
    puttr(&mut tr);
    exit(CWP_Exit());
}

// Local polynomial routines from sugoupillaud.c. All arithmetic is double;
// loop nesting and accumulation order follow the original routines.
fn imin(c1: i32, c2: i32) -> i32 {
    if c1 <= c2 {
        c1
    } else {
        c2
    }
}
fn imax(c1: i32, c2: i32) -> i32 {
    if c1 >= c2 {
        c1
    } else {
        c2
    }
}

fn porder(p: &[f64], l: i32) -> i32 {
    let (mut i, mut k) = (0, 0);
    loop {
        let j = k;
        // For an identically zero polynomial, C reads before its allocation
        // and returns an order below -1 on the reference build. The precise
        // value does not affect division (its numerator loop is empty), but
        // rev then fails its negative-sized allocation. Reproduce that result
        // explicitly rather than reading allocator metadata.
        if l - 1 - i < 0 {
            return -2;
        }
        if p[(l - 1 - i) as usize] == 0.0 {
            k += 1;
        }
        i += 1;
        if k <= j {
            return l - 1 - j;
        }
    }
}

fn pcomb(alfa: f64, beta: f64, p1: &[f64], p2: &[f64], l1: i32, l2: i32, l: i32) -> Vec<f64> {
    let ii = imax(l1, l2);
    // C copies tmp[ii] when l>ii. The extra cleared slot makes this explicit.
    let mut tmp = vec![0.0; (ii + 1) as usize];
    if l1 >= l2 {
        for i in (l2..l1).rev() {
            tmp[i as usize] = alfa * p1[i as usize];
        }
    } else {
        for i in (l1..l2).rev() {
            tmp[i as usize] = beta * p2[i as usize];
        }
    }
    for i in 0..imin(l1, l2) {
        tmp[i as usize] = alfa * p1[i as usize] + beta * p2[i as usize];
    }
    let mut c = vec![0.0; l as usize];
    for i in 0..=imin(ii, l - 1) {
        c[i as usize] = tmp[i as usize];
    }
    c
}

fn rev(p: &[f64], l: i32) -> Vec<f64> {
    let order = porder(p, l);
    if order < -1 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut tmp = vec![0.0; (order + 1) as usize];
    for i in 0..=order {
        tmp[i as usize] = p[(order - i) as usize];
    }
    let mut pr = vec![0.0; l as usize];
    for i in 0..=order {
        pr[i as usize] = tmp[i as usize];
    }
    pr
}

fn pshift(p: &[f64], l: i32, k: i32, ll: i32) -> Vec<f64> {
    let ii = imax(l + k, l);
    let mut tmp = vec![0.0; ii as usize];
    if k >= 0 {
        for i in 0..l {
            tmp[(i + k) as usize] = p[i as usize];
        }
    } else {
        for i in -k..l {
            tmp[(i + k) as usize] = p[i as usize];
        }
    }
    let mut psh = vec![0.0; ll as usize];
    for i in 0..imin(l, ll) {
        psh[i as usize] = tmp[i as usize];
    }
    psh
}

fn prod(p1: &[f64], p2: &[f64], l1: i32, l2: i32, ll: i32) -> Vec<f64> {
    let mut tmp = vec![0.0; (l1 + l2 - 1) as usize];
    for i in 0..l1 {
        for j in 0..l2 {
            let k = i + j;
            tmp[k as usize] = tmp[k as usize] + p1[i as usize] * p2[j as usize];
        }
    }
    let mut pp = vec![0.0; ll as usize];
    for i in 0..imin(l1 + l2 - 1, ll) {
        pp[i as usize] = tmp[i as usize];
    }
    pp
}

#[allow(non_snake_case)]
fn pratio(p1: &[f64], p2: &[f64], l1: i32, l2: i32, L: i32, lq: i32) -> Vec<f64> {
    let (o1, o2) = (porder(p1, l1), porder(p2, l2));
    let mut tmp1 = vec![0.0; (L + o2 + 1) as usize];
    let (mut i, mut k) = (0, 0);
    let jj = loop {
        let jj = k;
        if p2[i as usize] == 0.0 {
            k += 1;
        }
        i += 1;
        if k <= jj {
            break jj;
        }
    };
    let tmp2 = pshift(p2, l2, -jj, o2 + 1);
    let lo = imin(o1, L + o2);
    for i in 0..=lo {
        tmp1[i as usize] = p1[i as usize];
    }
    i = 0;
    loop {
        tmp1[i as usize] = tmp1[i as usize] / tmp2[0];
        if i == L + jj {
            break;
        }
        k = i;
        let ii = imin(o2, L + jj - i);
        for j in 0..ii {
            k += 1;
            tmp1[k as usize] = tmp1[k as usize] - tmp1[i as usize] * tmp2[(j + 1) as usize];
        }
        i += 1;
    }
    tmp1 = pshift(&tmp1, L + o2 + 1, -jj, L + o2 + 1);
    let mut q = vec![0.0; lq as usize];
    for i in 0..imin(L, lq) {
        q[i as usize] = tmp1[i as usize];
    }
    q
}

fn recurs(k: i32, n: i32, r: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut tq = vec![0.0; (n + 1) as usize];
    let mut tp = vec![0.0; (n + 1) as usize];
    tp[0] = 1.0;
    tq[1] = -r[1];
    for j in 2..=k {
        let mut qr = rev(&tq, n + 1);
        qr = pshift(&qr, n + 1, j - 1 - porder(&tq, n + 1), n + 1);
        let mut pr = rev(&tp, n + 1);
        pr = pshift(&pr, n + 1, j - 1 - porder(&tp, n + 1), n + 1);
        let qsh = pshift(&qr, n + 1, 1, n + 1);
        let psh = pshift(&pr, n + 1, 1, n + 1);
        tq = pcomb(1.0, -r[j as usize], &tq, &psh, n + 1, n + 1, n + 1);
        tp = pcomb(1.0, -r[j as usize], &tp, &qsh, n + 1, n + 1, n + 1);
    }
    (tq, tp)
}

fn upsample(p4: &[f64], l: i32, ll: i32) -> Vec<f64> {
    let ii = imax(ll, 2 * l);
    let mut tmp = vec![0.0; ii as usize];
    let mut i = 0;
    loop {
        // The C do/while reads p4[0] even for l=0; ll is then zero, so
        // this read has no observable output and can be omitted safely.
        if l == 0 {
            break;
        }
        tmp[(2 * i) as usize] = p4[i as usize];
        i += 1;
        if i >= l {
            break;
        }
    }
    let mut p2 = vec![0.0; ll as usize];
    for i in 0..ll {
        p2[i as usize] = tmp[i as usize];
    }
    p2
}
