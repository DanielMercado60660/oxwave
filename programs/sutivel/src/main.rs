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

//! SUTIVEL - transversely isotropic velocity table builder.
//! Port of su/main/velocity_analysis/sutivel.c (Revision 1.5).
//! Credits: John E. Anderson, visiting CSM from Mobil, Spring 1994.
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    "  SUTIVEL -  SU Transversely Isotropic velocity table builder\t\t",
    "\tcomputes vnmo or vphase as a function of Thomsen's parameters and",
    "\ttheta and optionally interpolate to constant increments in slowness",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " a=2500.\t\talpha (vertical p velocity)\t\t\t",
    " b=1250.\t\tbeta (vertical sv velocity)\t\t\t",
    " e=.20\t\t\tepsilon (horiz p-wave anisotropy)\t\t",
    " d=.10\t\t\tdelta (strange parameter)\t\t\t",
    " maxangle=90.0\t\tmax angle in degrees\t\t\t\t",
    " nangle=9001\t\tnumber of angles to compute\t\t\t",
    " verbose=0\t\tset to 1 to see full listing\t\t\t",
    " np=8001\t\tnumber of slowness values to output\t\t",
    " option=1\t\t1=output vnmo(p) (result used for TI DMO)\t",
    "\t\t\t2=output vnmo(theta) in degrees\t\t\t",
    "\t\t\t3=output vnmo(theta) in radians\t\t\t",
    "\t\t\t4=output vphase(p)\t\t\t\t",
    "\t\t\t5=output vphase(theta) in degrees\t\t",
    "\t\t\t6=output vphase(theta) in radians\t\t",
    "\t\t\t7=output first derivative vphase(p)\t\t",
    "\t\t\t8=output first derivative vphase(theta) in degrees",
    "\t\t\t9=output first derivative vphase(theta) in radians",
    "\t\t\t10=output second derivative vphase(p)\t\t",
    "\t\t\t11=output second derivative vphase(theta) in degrees",
    "\t\t\t12=output second derivative vphase(theta) in radians",
    "\t\t\t13=( 1/vnmo(0)^2 -1/vnmo(theta)^2 )/p^2 test vs theta",
    "\t\t\t   (result should be zero for all theta for d=e)",
    "\t\t\t14=return vnmo(p) for weak anisotropy\t\t",
    " normalize=0\t\t=1 means scale vnmo by cosine and scale vphase by",
    " \t\t\t    1/sqrt(1+2*e*sin(theta)*sin(theta)\t\t",
    "\t \t\t   (only useful for vphase when d=e for constant",
    "\t\t\t\tresult)\t\t\t\t\t",
    "\t\t\t=0 means output vnmo or vphase unnormalized\t",
    "\t\t\t\t\t\t\t\t\t",
    " Output on standard output is ascii text with:\t\t\t\t",
    " line   1: number of values\t\t\t\t\t\t",
    " line   2: abscissa increment (p or theta increment, always starts at zero)",
    " line 3-n: one value per line\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];
fn vget(a: f32, b: f32, e: f32, d: f32, theta: f32) -> [f32; 5] {
    let sint = (theta as f64).sin() as f32;
    let mut cost = (theta as f64).cos() as f32;
    let sin2t = sint * sint;
    let cos2t = cost * cost;
    let psi = (1.0 - (b * b / (a * a)) as f64) as f32;
    let f = psi * (2.0 * d - e);
    let g = (psi + e) * e;
    let gamma = (0.25 * psi as f64 * psi as f64
        + (f * sin2t * cos2t) as f64
        + (g * sin2t * sin2t) as f64) as f32;
    let dg = f * 2.0 * sint * cost * cos2t + (4.0 * g - 2.0 * f) * sint * sin2t * cost;
    let d2g = f * 2.0 * cos2t * cos2t
        + (g - f) * 12.0 * cos2t * sin2t
        + (2.0 * f - 4.0 * g) * sin2t * sin2t;
    let sq = (gamma as f64).sqrt() as f32;
    let v = (a as f64 * (1.0 + (e * sin2t) as f64 - 0.5 * psi as f64 + sq as f64).sqrt()) as f32;
    let dv =
        (0.5 * a as f64 * a as f64 * ((2.0 * e * sint * cost) as f64 + 0.5 * dg as f64 / sq as f64)
            / v as f64) as f32;
    let d2v = (0.5
        * a as f64
        * a as f64
        * ((2.0 * e * (cos2t - sin2t)) as f64 + 0.5 * d2g as f64 / sq as f64
            - 0.25 * dg as f64 * dg as f64 / (sq * gamma) as f64)
        / v as f64
        - (dv * dv / v) as f64) as f32;
    if cost < 1.0e-7_f32 {
        cost = 1.0e-7_f32;
    }
    let tant = sint / cost;
    let vnmo = (v as f64 * ((1.0f32 + d2v / v) as f64).sqrt()
        / (cost * (1.0 - tant * dv / v)) as f64) as f32;
    [v, dv, d2v, vnmo, sint / v]
}
fn vnmoweak_ti(a: f32, _b: f32, e: f32, d: f32, p: f32) -> f32 {
    let vnmo0 = (a as f64 * ((1.0f32 + 2.0 * d) as f64).sqrt()) as f32;
    let mut y = p * vnmo0;
    y = y * y;
    let dd = p * p * (1.0 + 2.0 * (e - d) * (4.0 * y * y - 9.0 * y + 6.0));
    y = a * a * (1.0 + 2.0 * d);
    y = 1.0 / y;
    if y > dd {
        (1.0 / ((y - dd) as f64).sqrt()) as f32
    } else {
        -1.0
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut a, mut b, mut e, mut d, mut maxangle) =
        (2500.0f32, 1250.0f32, 0.2f32, 0.1f32, 90.0f32);
    let (mut nangle, mut np, mut verbose, mut option, mut normalize) =
        (9001i32, 8001i32, 0i32, 1i32, 0i32);
    getparfloat("a", &mut a);
    getparfloat("b", &mut b);
    getparfloat("e", &mut e);
    getparfloat("d", &mut d);
    getparfloat("maxangle", &mut maxangle);
    getparint("nangle", &mut nangle);
    getparint("np", &mut np);
    getparint("verbose", &mut verbose);
    getparint("option", &mut option);
    getparint("normalize", &mut normalize);
    fprintf!(stderr(), "normalize=%d\n", normalize);
    checkpars();
    if option == 14 {
        let mut pmax = (1.0 / (a as f64 * ((1.0f32 + 2.0 * e) as f64).sqrt())) as f32;
        let mut dp = pmax / (np - 1) as f32;
        for ip in 0..np {
            if vnmoweak_ti(a, b, e, d, ip as f32 * dp) < 0.0 {
                pmax = (ip - 1) as f32 * dp;
                break;
            }
        }
        dp = pmax / (np - 1) as f32;
        printf!("%d\n", np);
        printf!("%e\n", dp);
        for ip in 0..np {
            printf!("%e\n", vnmoweak_ti(a, b, e, d, ip as f32 * dp));
        }
        exit(CWP_Exit());
    }
    if nangle < 2 {
        nangle = 2;
    }
    let mut pin = vec![0.0f32; nangle as usize];
    let mut vin = pin.clone();
    let dangle = maxangle / (nangle - 1) as f32;
    if verbose != 0 {
        printf!(" alpha = a = %f\n", a);
        printf!(" beta = b = %f\n", b);
        printf!(" epsilon = e = %f\n", e);
        printf!(" delta = d = %f\n", d);
        printf!(
            " a*sqrt(1+2*d) = %f\n",
            a as f64 * ((1.0f32 + 2.0 * d) as f64).sqrt()
        );
        printf!(
            " a*sqrt(1+2*e) = %f\n",
            a as f64 * ((1.0f32 + 2.0 * e) as f64).sqrt()
        );
    }
    if [2, 3, 5, 6, 8, 9, 11, 12, 13].contains(&option) {
        if verbose != 0 {
            printf!(" The number of angles =");
        }
        printf!("%d\n", nangle);
        if verbose != 0 {
            if [2, 5, 8, 11, 13].contains(&option) {
                printf!(" The angle increment in degrees =");
            } else {
                printf!(" The angle increment in radians =");
            }
        }
        if [2, 5, 8, 11, 13].contains(&option) {
            printf!("%e\n", dangle);
        } else {
            printf!("%e\n", dangle as f64 * std::f64::consts::PI / 180.0);
        }
        if verbose != 0 {
            printf!(" Following are velocity values for constant increments in angle\n");
        }
    }
    for i in 0..nangle {
        let angle = i as f32 * dangle;
        let theta = (angle as f64 * std::f64::consts::PI / 180.0) as f32;
        let mut vel = vget(a, b, e, d, theta);
        if normalize != 0 {
            vel[0] = (vel[0] as f64
                / (1.0 + (2.0f32 * e) as f64 * (theta as f64).sin() * (theta as f64).sin()).sqrt())
                as f32;
            vel[3] = (vel[3] as f64 * (theta as f64).cos()) as f32;
        }
        let [v, dv, d2v, vnmo, p] = vel;
        pin[i as usize] = p;
        vin[i as usize] = vnmo;
        if option == 4 {
            vin[i as usize] = v;
        }
        if option == 7 {
            vin[i as usize] = dv;
        }
        if option == 10 {
            vin[i as usize] = d2v;
        }
        if option == 2 || option == 3 {
            printf!("%e\n", vnmo);
        }
        if option == 5 || option == 6 {
            printf!("%e\n", v);
        }
        if option == 8 || option == 9 {
            printf!("%e\n", dv);
        }
        if option == 11 || option == 12 {
            printf!("%e\n", d2v);
        }
        if option == 13 {
            let test = ((p * p) as f64 + 1.0 / (vnmo * vnmo) as f64
                - 1.0 / (vin[0] * vin[0]) as f64) as f32;
            printf!("%e\n", test);
        }
        if verbose != 0 {
            printf!(" angle (in degrees) = %6.2f", angle);
            printf!(" vphase = %9.2f", v);
            printf!(" vnmo = %15.6e", vnmo);
            printf!(" p = %15.6e\n", p);
            printf!("\t  dv/dtheta = %15.6e", dv);
            printf!(" d2v/dtheta2 = %15.6e\n", d2v);
        }
    }
    if option == 1 || option == 4 {
        let mut pmax = pin[0];
        for &p in &pin {
            pmax = if pmax > p { pmax } else { p };
        }
        let dp = pmax / (np - 4) as f32;
        let mut pout = vec![0.0f32; np.max(0) as usize];
        let mut vout = pout.clone();
        for ip in 0..np {
            pout[ip as usize] = ip as f32 * dp;
        }
        intlin(
            &pin,
            &vin,
            vin[0],
            vin[nangle as usize - 1],
            &pout,
            &mut vout,
        );
        if verbose != 0 {
            printf!(" np = ");
        }
        printf!("%d\n", np);
        if verbose != 0 {
            printf!(" dp = ");
        }
        printf!("%e\n", dp);
        if verbose != 0 {
            printf!(" Following are velocity values for constant increments in p\n");
        }
        for v in vout {
            printf!("%e\n", v);
        }
    }
    exit(CWP_Exit());
}
