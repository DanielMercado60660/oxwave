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

//! ELASYN - Synthetic seismograms for triangulated elastic media.
//! Port of Trielas/main/elasyn.c (Test Release 1.2, 2011/11/21).
//! Copyright (c) Colorado School of Mines, 2011, 1999. All rights reserved.
//! Authors: Andreas Rueger and Tariq Alkalifah, Colorado School of Mines, 02/02/94.
//! Based on gbbeam.c, Andreas Rueger, 08/12/93, and sdbeam.c, Dave Hale, CSM, 02/26/91.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{abs, sgn, PI};
use oxwave_f_cubicspline::csplin;
use oxwave_inverse_interp::yxtoxy;
use oxwave_pfa_fft::{npfaro, pfacr};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " ELASYN - SYNthetic seismograms for triangulated elastic media\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  elasyn <rayends xg= zg= [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Required Parameters:\t\t\t\t\t\t\t",
    "xg=            x coordinates of receiver surface\t\t\t",
    "zg=            z coordinates of receiver surface\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Optional Parameters:\t\t\t\t\t\t\t",
    "compon=0         horizontal and vertical component seismograms         ",
    "\t\t  =3 vertical component (positive downwards)\t\t",
    "\t\t  =1 horizontal component\t\t\t\t",
    "ng=101           number of receivers (uniform distributed along surface)",
    "krecord=1        integer index of receiver surface (see notes below)\t",
    "nt=251           number of time samples\t\t\t\t",
    "dt=0.004         time sampling interval\t\t\t\t",
    "ft=0.0           first time sample\t\t\t\t\t",
    "inter=0          linear interpolation\t\t\t\t\t",
    "inter=1 (default) cross parabolic interpolation \t\t\t",
    "reftrans=0       =1 complex refl/transm. coefficients considered \t",
    "nameref=-1       all rays recorded at interface <krecord> considered \t",
    "                 =0, only direct hits are considered  \t\t\t",
    "                 >0, only rays reflected from interface <nameref>      ",
    "lscale=          if defined restricts range of extrapolation\t\t",
    "fpeak=0.1/dt     peak frequency of ricker wavelet \t\t\t",
    "infofile         ASCII-file to store useful information \t\t",
    "xfile=x_compon.bin     bin-file to store x_component traces \t\t",
    "zfile=z_compon.bin     bin-file to store z_component traces \t\t",
    "\t\t\t\t\t\t\t\t\t",
    "NOTES:\t\t\t\t\t\t\t\t\t",
    "Only rays that terminate with index krecord will contribute to the\t",
    "synthetic seismograms at the receiver (xg,zg) locations.  The\t\t",
    "receiver locations are determined by cubic spline interpolation\t",
    "of the specified (xg,zg) coordinates.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Warning!!-- This version is not quite complete. There is a bug in the ",
    " interpolation routines that causes a segmentation violation on the last",
    " couple  of traces.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// RayEnd is a native binary record, not a SEG-Y trace. Keep all 27 words,
// including the unused dynamic ray parameters, in their original order.
#[derive(Clone, Default)]
struct RayEnd {
    words: [u32; 27],
}
impl RayEnd {
    fn f(&self, i: usize) -> f32 {
        f32::from_bits(self.words[i])
    }
    fn i(&self, i: usize) -> i32 {
        self.words[i] as i32
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);

    let nxg = countparval("xg");
    let nzg = countparval("zg");
    if nxg == 0 || nzg == 0 {
        err!("must specify both xg and zg");
    }
    if nxg != nzg {
        err!("number of xg must equal number of zg");
    }
    let nxz = nxg as usize;
    let mut ng = 101;
    getparint("ng", &mut ng);
    let mut xg = vec![0.0f32; nxz.max(ng.max(0) as usize)];
    let mut zg = xg.clone();
    getparfloat("xg", &mut xg[..]);
    getparfloat("zg", &mut zg[..]);
    let mut ifile = String::new();
    let ifp = if getparstring("infofile", &mut ifile) != 0 {
        Some(efopen(&ifile, "w"))
    } else {
        None
    };
    let mut lscale = f32::MAX;
    let mut krecord = 1;
    let mut compon = 0;
    let mut nt = 251;
    let mut dt = 0.004f32;
    let mut ft = 0.0f32;
    let mut reftrans = 0;
    getparfloat("lscale", &mut lscale);
    getparint("krecord", &mut krecord);
    getparint("compon", &mut compon);
    getparint("nt", &mut nt);
    getparfloat("dt", &mut dt);
    getparfloat("ft", &mut ft);
    getparint("reftrans", &mut reftrans);
    let mut fpeak = (0.1 / dt as f64) as f32;
    let mut nameref = -1;
    let mut inter = 1;
    getparfloat("fpeak", &mut fpeak);
    getparint("nameref", &mut nameref);
    getparint("inter", &mut inter);
    let mut xfile = String::new();
    let mut zfile = String::new();
    let xfp = if compon == 0 || compon == 1 {
        if getparstring("xfile", &mut xfile) == 0 {
            xfile = "x_compon.bin".into();
        }
        Some(efopen(&xfile, "w"))
    } else {
        None
    };
    let zfp = if compon == 0 || compon == 3 {
        if getparstring("zfile", &mut zfile) == 0 {
            zfile = "z_compon.bin".into();
        }
        Some(efopen(&zfile, "w"))
    } else {
        None
    };
    checkpars();
    if 0.5 / (dt as f64) < 2.0 * fpeak as f64 {
        fprintf!(stderr(), " WARNING: ALIASING POSSIBLE \n");
        fprintf!(stderr(), " decrease dt or reduce fpeak \n");
    }
    if reftrans == 1 {
        warn!("\n WARNING currently no refl/trans coeff supported \n");
    }
    warn!(" \n NOTE This is a very simple, go-bye seismogram builder\n ");
    // Negative dimensions become enormous size_t allocations in the C code.
    if ng < -1 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }

    let mut nre = 0usize;
    let mut nrealloc = 200;
    let mut re = vec![RayEnd::default(); nrealloc];
    // Allocate the arrival grids before reading the endpoints, as in elasyn.c.
    // The extra row is read after makesyn's deliberate one-row shift.
    let nrows = ng.wrapping_add(1).max(0) as usize;
    let mut tg = vec![[0.0f32; 3]; nrows];
    let mut axg = tg.clone();
    let mut azg = tg.clone();
    let mut phg = tg.clone();
    let mut im = vec![0i32; nrows];

    // The C allocation has a spare endpoint, also read by writeinfo (<=nre).
    loop {
        let mut bytes = [0u8; 108];
        let n = fread(&mut bytes, &stdin());
        for i in 0..27 {
            re[nre].words[i] = u32::from_ne_bytes(bytes[4 * i..4 * i + 4].try_into().unwrap());
        }
        if n != bytes.len() {
            break;
        }
        nre += 1;
        if nre == nrealloc {
            nrealloc += 200;
            re.resize(nrealloc, RayEnd::default());
        }
    }
    // C reads x[na], a[na] in interpara even on its last (linear) interval.
    // Retain a zero spare word, matching the original allocation's padding.
    let mut ax = vec![0.0f32; nre + 2];
    let mut az = ax.clone();
    let mut t = ax.clone();
    let mut x = ax.clone();
    let mut z = ax.clone();
    let mut num = vec![0i32; nre + 2];
    let mut px = ax.clone();
    let mut pz = ax.clone();
    let mut ph = ax.clone();
    let dangle = re[0].f(18);
    let mut xmax = -f32::MAX;
    let mut xmin = f32::MAX;
    if let Some(fp) = &ifp {
        writeinfo(&re, fp, nre);
    }
    if ng < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    makexzs(nxz, &mut xg, &mut zg, ng.max(0) as usize);
    let mut nri = 0;
    for r in &re[..nre] {
        if r.i(15) == krecord && (nameref == -1 || nameref == r.i(23)) {
            t[nri] = r.f(5);
            x[nri] = r.f(1);
            z[nri] = r.f(2);
            num[nri] = r.i(24);
            pz[nri] = r.f(4);
            px[nri] = r.f(3);
            ph[nri] = r.f(17);
            let (polarx, polarz) = polar(px[nri], pz[nri], r.f(20), r.f(21), r.f(22), r.i(19));
            ax[nri] = polarx * r.f(16);
            az[nri] = -polarz * r.f(16);
            xmax = if x[nri] > xmax { x[nri] } else { xmax };
            xmin = if x[nri] < xmin { x[nri] } else { xmin };
            nri += 1;
        }
    }
    drop(re);
    if nri == 0 {
        warn!("\n NO RAY CONTRIBUTION. Check <krecord> \n");
        if nt < 0 {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        let zero = vec![0.0f32; nt.max(0) as usize];
        for _ in 0..ng {
            if let Some(fp) = &xfp {
                fwrite(&zero, fp);
            }
            if let Some(fp) = &zfp {
                fwrite(&zero, fp);
            }
        }
        exit(1);
    }
    if lscale == f32::MAX {
        lscale = (xmax - xmin) * 20.0f32 / nre as f32;
    }
    let ng = ng.max(0) as usize;
    // For ng=0 C reads the zero padding word at xg[-1].
    let mut dg = ((if ng == 0 { 0.0 } else { xg[ng - 1] }) - xg[0])
        / if ng > 1 { (ng - 1) as f32 } else { 1.0f32 };
    if ng == 1 {
        dg = (xmax - xmin) * 5.0f32 / nre as f32;
    }
    ray_jacob(dg, &x, &num, &pz, &mut ph, &mut ax, &mut az, nri);
    if inter == 1 {
        interpara(nri, &t, &x, ng, dg, xg[0], lscale, &mut im, &mut tg);
        interpara(nri, &ax, &x, ng, dg, xg[0], lscale, &mut im, &mut axg);
        interpara(nri, &az, &x, ng, dg, xg[0], lscale, &mut im, &mut azg);
        interpara(nri, &ph, &x, ng, dg, xg[0], lscale, &mut im, &mut phg);
    } else {
        interplin(nri, &t, &x, ng, dg, xg[0], lscale, &mut im, &mut tg);
        interplin(nri, &ax, &x, ng, dg, xg[0], lscale, &mut im, &mut axg);
        interplin(nri, &az, &x, ng, dg, xg[0], lscale, &mut im, &mut azg);
        interplin(nri, &ph, &x, ng, dg, xg[0], lscale, &mut im, &mut phg);
    }
    // The source deliberately shifts the interpolation arrays by one.
    makesyn(
        reftrans,
        ng,
        &xg,
        compon,
        &zg,
        fpeak,
        nt,
        dt,
        ft,
        ifp.as_ref(),
        &axg[1..],
        &azg[1..],
        &phg[1..],
        &tg[1..],
        dangle,
        &im[1..],
        xfp.as_ref(),
        zfp.as_ref(),
    );
    // elasyn.c returns 1 after successful synthesis, as well as for no rays.
    exit(1);
}

fn makexzs(nu: usize, xu: &mut [f32], zu: &mut [f32], ns: usize) {
    let mut xud = vec![[0.0f32; 4]; nu];
    let mut zud = xud.clone();
    let mut u = vec![0.0f32; nu];
    for iu in 0..nu {
        u[iu] = iu as f32;
    }
    csplin(&u, &xu[..nu], &mut xud);
    csplin(&u, &zu[..nu], &mut zud);
    let nuu = 20 * nu;
    let duu = (u[nu - 1] - u[0]) / (nuu - 1) as f32;
    let mut s = vec![0.0f32; nuu];
    let mut xlast = xu[0];
    let mut zlast = zu[0];
    let mut uu = 0.0f32;
    for iuu in 0..nuu {
        let mut x = [0.0f32];
        let mut z = [0.0f32];
        intcub(0, &u, &xud, &[uu], &mut x);
        intcub(0, &u, &zud, &[uu], &mut z);
        let dx = x[0] - xlast;
        let dz = z[0] - zlast;
        // C reads s[-1] on iteration zero. Use zero padding for this undefined
        // read; for some curved surfaces the original instead reads a live
        // preceding allocation, so those cases remain heap-layout-dependent.
        s[iuu] = ((if iuu == 0 { 0.0 } else { s[iuu - 1] }) as f64
            + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
        xlast = x[0];
        zlast = z[0];
        uu += duu;
    }
    let mut us = vec![0.0f32; ns];
    let ds = (s[nuu - 1] - s[0]) / if ns > 1 { (ns - 1) as f32 } else { 1.0f32 };
    yxtoxy(duu, 0.0, &s, ds, 0.0, 0.0, (nu - 1) as f32, &mut us);
    intcub(0, &u, &xud, &us, &mut xu[..ns]);
    intcub(0, &u, &zud, &us, &mut zu[..ns]);
}

fn cricker(w: f32, wpeak: f32, delay: f32) -> Complex {
    crmul(
        cwp_cexp(cmplx(-((w / wpeak) as f64).powf(2.0) as f32, delay * w)),
        (4.0 * w as f64 * w as f64 * PI.sqrt() / (wpeak as f64).powf(3.0)) as f32,
    )
}

fn writeinfo(re: &[RayEnd], ifp: &CFile, nre: usize) {
    fprintf!(
        ifp,
        "\n The following information is stored at the rayend \n"
    );
    for r in &re[..=nre] {
        fprintf!(
            ifp,
            "\n Ray stops at (%g,%g) at interface %i\n",
            r.f(1),
            r.f(2),
            r.i(15)
        );
        fprintf!(
            ifp,
            " Ray number %i \tTakeoff angle =%g  \n",
            r.i(24),
            r.f(14)
        );
        fprintf!(ifp, " Slowness components px=%g pz=%g\n", r.f(3), r.f(4));
        fprintf!(ifp, " Traveltime =%g  \t Sigma=%g\n", r.f(5), r.f(0));
        fprintf!(
            ifp,
            " kmah index=%i \t Number of reflections=%i\n",
            r.i(10),
            r.i(11)
        );
        fprintf!(
            ifp,
            " last reflector=%i \t Angle increment=%g \n",
            r.i(23),
            r.f(18) as f64 * 180.0 / PI
        );
        fprintf!(
            ifp,
            " g11=%g \t g33=%g \t g13=%g\n",
            r.f(20),
            r.f(22),
            r.f(21)
        );
        fprintf!(ifp, " vgx=%g \t vgz=%g \t \n", r.f(12), r.f(13));
        fprintf!(
            ifp,
            " Refl/Transm Amplitude =%g Phase=%g\n",
            r.f(16),
            r.f(17)
        );
        fprintf!(
            ifp,
            "\n --------------------------------------------------- \n"
        );
    }
}

#[allow(clippy::too_many_arguments, unused_assignments)]
fn makesyn(
    _reftrans: i32,
    ng: usize,
    _xg: &[f32],
    compon: i32,
    _zg: &[f32],
    fpeak: f32,
    nt: i32,
    dt: f32,
    ft: f32,
    _ifp: Option<&CFile>,
    axg: &[[f32; 3]],
    azg: &[[f32; 3]],
    phg: &[[f32; 3]],
    tg: &[[f32; 3]],
    _dangle: f32,
    im: &[i32],
    xfp: Option<&CFile>,
    zfp: Option<&CFile>,
) {
    let ntfft = npfaro(nt, 2 * nt) as usize;
    let nw = ntfft / 2 + 1;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let wpeak = (2.0 * PI * fpeak as f64) as f32;
    let delay = 0.0f32;
    let mut syn = vec![0.0f32; ntfft];
    let mut cwave = vec![Complex::default(); nw];
    let mut csynx = vec![vec![Complex::default(); nw]; ng];
    let mut csynz = csynx.clone();
    let mut w = 0.0f32;
    for iw in 0..nw {
        cwave[iw] = cricker(w, wpeak, delay);
        cwave[iw] = crmul(cwave[iw], (w as f64).sqrt() as f32);
        w += dw;
    }
    for ig in 0..ng {
        fprintf!(
            stderr(),
            "ig=%i im=%i tg=%f a=%f \n",
            ig,
            im[ig],
            tg[ig][0],
            axg[ig][0]
        );
        for it in 0..im[ig] as usize {
            fprintf!(
                stderr(),
                "ig=%i im=%i tg=%f a=%f phg=%f \n",
                ig,
                im[ig],
                tg[ig][0],
                axg[ig][0],
                arrival(phg, ig, it)
            );
            let cosf = (arrival(phg, ig, it) as f64).cos() as f32;
            let sinf = (arrival(phg, ig, it) as f64).sin() as f32;
            let cosd = ((dw * (arrival(tg, ig, it) - ft)) as f64).cos();
            let sind = ((dw * (arrival(tg, ig, it) - ft)) as f64).sin();
            let mut cosw = 1.0f64;
            let mut sinw = 0.0f64;
            w = 0.0;
            for iw in 0..nw {
                let campr = cosf as f64 * cosw - sinf as f64 * sinw;
                let campi = cosf as f64 * sinw + sinf as f64 * cosw;
                csynx[ig][iw].r =
                    (csynx[ig][iw].r as f64 + campr * arrival(axg, ig, it) as f64) as f32;
                csynx[ig][iw].i =
                    (csynx[ig][iw].i as f64 + campi * arrival(axg, ig, it) as f64) as f32;
                csynz[ig][iw].r =
                    (csynz[ig][iw].r as f64 + campr * arrival(azg, ig, it) as f64) as f32;
                csynz[ig][iw].i =
                    (csynz[ig][iw].i as f64 + campi * arrival(azg, ig, it) as f64) as f32;
                let tempd = cosw * cosd - sinw * sind;
                sinw = cosw * sind + sinw * cosd;
                cosw = tempd;
                w += dw;
            }
        }
        if compon == 1 || compon == 0 {
            for iw in 0..nw {
                csynx[ig][iw] = cmul(csynx[ig][iw], cwave[iw]);
            }
            pfacr(-1, &csynx[ig], &mut syn);
            fwrite(&syn[..nt as usize], xfp.unwrap());
        }
        if compon == 3 || compon == 0 {
            for iw in 0..nw {
                csynz[ig][iw] = cmul(csynz[ig][iw], cwave[iw]);
            }
            pfacr(-1, &csynz[ig], &mut syn);
            fwrite(&syn[..nt as usize], zfp.unwrap());
        }
    }
}

fn polar(px: f32, pz: f32, g11: f32, g13: f32, g33: f32, mode: i32) -> (f32, f32) {
    let (mut polarx, mut polarz);
    if g11 < f32::EPSILON {
        polarx = 0.0;
        polarz = (g33 as f64).sqrt() as f32;
    } else if g33 < f32::EPSILON {
        polarz = 0.0;
        polarx = (g11 as f64).sqrt() as f32;
    } else {
        polarx = ((g11 as f64).sqrt() * sgn(g13)) as f32;
        polarz = (g33 as f64).sqrt() as f32;
    }
    if ((mode == 0 || mode == 3) && px * polarx + pz * polarz < 0.0)
        || ((mode == 4 || mode == 1) && px * polarx < 0.0)
    {
        polarx = -polarx;
        polarz = -polarz;
    }
    (polarx, polarz)
}

fn ray_jacob(
    dg: f32,
    x: &[f32],
    num: &[i32],
    pz: &[f32],
    _ph: &mut [f32],
    ax: &mut [f32],
    az: &mut [f32],
    nrays: usize,
) {
    let mut jac = (x[0] - x[1]) * pz[0];
    let mut diffnum = num[1] - num[0];
    if diffnum != 1 {
        jac /= diffnum as f32;
    }
    // C's ABS preserves negative zero. Its square root and the subsequent
    // division must preserve that sign for rays with zero vertical slowness.
    ax[0] = (ax[0] as f64 / (abs(jac) as f64).sqrt()) as f32;
    az[0] = (az[0] as f64 / (abs(jac) as f64).sqrt()) as f32;
    _ph[0] = 0.0;
    for iray in 1..nrays.saturating_sub(1) {
        let mut jac1 = (x[iray] - x[iray + 1]) * pz[iray];
        diffnum = num[iray + 1] - num[iray];
        if diffnum != 1 {
            jac1 /= diffnum as f32;
        }
        let mut jac2 = (x[iray - 1] - x[iray]) * pz[iray];
        diffnum = num[iray] - num[iray - 1];
        if diffnum != 1 {
            jac2 /= diffnum as f32;
        }
        if ((jac1 * jac1) as f64) < 0.001 * dg as f64 * dg as f64
            && ((jac1 * jac1) as f64) < 0.001 * dg as f64 * dg as f64
        {
            // jac=jac in C (the repeated jac1 test is intentional).
        } else if ((jac1 * jac1) as f64) < 0.001 * dg as f64 * dg as f64 {
            jac = jac2;
        } else if ((jac2 * jac2) as f64) < 0.001 * dg as f64 * dg as f64 {
            jac = jac1;
        } else {
            jac = ((jac1 + jac2) as f64 / 2.0) as f32;
        }
        ax[iray] = (ax[iray] as f64 / (abs(jac) as f64).sqrt()) as f32;
        az[iray] = (az[iray] as f64 / (abs(jac) as f64).sqrt()) as f32;
    }
    // For a single ray C reads x[-1], the zero word preceding this allocation.
    jac = (x[nrays - 1] - if nrays == 1 { 0.0 } else { x[nrays - 2] }) * pz[nrays - 1];
    diffnum = num[1] - num[0];
    if diffnum != 1 {
        jac /= diffnum as f32;
    }
    ax[nrays - 1] = (ax[nrays - 1] as f64 / (abs(jac) as f64).sqrt()) as f32;
    az[nrays - 1] = (az[nrays - 1] as f64 / (abs(jac) as f64).sqrt()) as f32;
}

#[allow(clippy::too_many_arguments)]
fn interplin(
    na: usize,
    a: &[f32],
    x: &[f32],
    ng: usize,
    dg: f32,
    fg: f32,
    tol: f32,
    im: &mut [i32],
    ag: &mut [[f32; 3]],
) {
    let mut itmin = 0i32;
    let mut itmax = 0i32;
    let odg = 1.0f32 / dg;
    for it in 0..ng {
        im[it] = 0;
    }
    let mut x1 = x[0] - fg;
    let mut iti = (x1 * odg) as i32;
    let mut aa1 = a[0];
    for ix in 1..na {
        let mut x2 = x[ix] - fg;
        let mut itf = (x2 * odg) as i32;
        let mut aa2 = a[ix];
        let mut sub = true;
        if itf < iti {
            std::mem::swap(&mut itf, &mut iti);
            std::mem::swap(&mut x1, &mut x2);
            std::mem::swap(&mut aa1, &mut aa2);
            sub = false;
        }
        if x2 - x1 > tol {
            for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                *arrival_mut(ag, it as usize, im[it as usize] as usize) = 0.0;
                im[it as usize] += 1;
            }
        } else {
            for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                itmin = itmin.min(it);
                itmax = itmax.max(it);
                let frac = (it as f32 * dg - x1) / (x2 - x1);
                *arrival_mut(ag, it as usize, im[it as usize] as usize) =
                    ((1.0 - frac as f64) * aa1 as f64 + (frac * aa2) as f64) as f32;
                im[it as usize] += 1;
            }
        }
        if sub {
            iti = itf;
            x1 = x2;
            aa1 = aa2;
        }
    }
    for it in 0..itmin as usize {
        ag[it][0] = 0.0;
        im[it] += 1;
    }
    for it in (itmax + 1) as usize..ng {
        ag[it][0] = 0.0;
        im[it] += 1;
    }
}

#[allow(clippy::too_many_arguments)]
fn interpara(
    na: usize,
    a: &[f32],
    x: &[f32],
    ng: usize,
    dg: f32,
    fg: f32,
    tol: f32,
    im: &mut [i32],
    ag: &mut [[f32; 3]],
) {
    let mut itmin = 0i32;
    let mut itmax = 0i32;
    let odg = 1.0f32 / dg;
    let eps = 0.001f32;
    let epsx = 0.01f32;
    let (mut a1, mut b1, mut c1) = (0.0f32, 0.0f32, 0.0f32);
    for it in 0..ng {
        im[it] = 0;
    }
    let mut x1 = x[0] - fg;
    let mut iti = (x1 * odg) as i32;
    let mut u1 = a[0];
    let mut ixl = -1i32;
    let mut x0 = 0.0f32;
    let mut u0 = 0.0f32;
    for ix in 1..na {
        let x2 = x[ix] - fg;
        let mut itf = (x2 * odg) as i32;
        let u2 = a[ix];
        let mut sub = true;
        let x3 = x[ix + 1] - fg;
        let u3 = a[ix + 1];
        if itf < iti {
            std::mem::swap(&mut itf, &mut iti);
            sub = false;
        }
        if (x2 - x1).abs() > tol {
            for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                *arrival_mut(ag, it as usize, im[it as usize] as usize) = 0.0;
                im[it as usize] += 1;
            }
        } else if ix == 1 || ix == na - 1 {
            let dxx = x2 - x1;
            if dxx.abs() < eps {
                for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                    itmin = itmin.min(it);
                    itmax = itmax.max(it);
                    *arrival_mut(ag, it as usize, im[it as usize] as usize) =
                        (0.5 * (u1 + u2) as f64) as f32;
                    im[it as usize] += 1;
                }
            } else {
                for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                    itmin = itmin.min(it);
                    itmax = itmax.max(it);
                    let frac = (it as f32 * dg - x1) / dxx;
                    *arrival_mut(ag, it as usize, im[it as usize] as usize) =
                        ((1.0 - frac as f64) * u1 as f64 + (frac * u2) as f64) as f32;
                    im[it as usize] += 1;
                }
            }
        } else {
            let dxx = x2 - x1;
            let dxx3 = dxx * (x3 - x2) * (x1 - x0);
            if dxx.abs() < eps {
                for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                    itmin = itmin.min(it);
                    itmax = itmax.max(it);
                    *arrival_mut(ag, it as usize, im[it as usize] as usize) =
                        (0.5 * (u1 + u2) as f64) as f32;
                    im[it as usize] += 1;
                }
            } else if dxx3.abs() < epsx {
                for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                    itmin = itmin.min(it);
                    itmax = itmax.max(it);
                    let frac = (it as f32 * dg - x1) / dxx;
                    *arrival_mut(ag, it as usize, im[it as usize] as usize) =
                        ((1.0 - frac as f64) * u1 as f64 + (frac * u2) as f64) as f32;
                    im[it as usize] += 1;
                }
            } else {
                if ixl != ix as i32 - 1 {
                    let num1 = ((x1 - x0) * (x2 - x0)) as f64;
                    let den1 = 1.0 / num1;
                    let num2 = ((x1 - x0) * (x2 - x1)) as f64;
                    let den2 = 1.0 / num2;
                    let num3 = ((x2 - x1) * (x2 - x0)) as f64;
                    let den3 = 1.0 / num3;
                    a1 = ((u0 * x1 * x2) as f64 * den1 - (u1 * x2 * x0) as f64 * den2
                        + (u2 * x0 * x1) as f64 * den3) as f32;
                    b1 = ((-u0 * (x1 + x2)) as f64 * den1 + (u1 * (x2 + x0)) as f64 * den2
                        - (u2 * (x1 + x0)) as f64 * den3) as f32;
                    c1 = (u0 as f64 * den1 - u1 as f64 * den2 + u2 as f64 * den3) as f32;
                }
                let num1 = ((x2 - x1) * (x3 - x1)) as f64;
                let den1 = 1.0 / num1;
                let num2 = ((x2 - x1) * (x3 - x2)) as f64;
                let den2 = 1.0 / num2;
                let num3 = ((x3 - x2) * (x3 - x1)) as f64;
                let den3 = 1.0 / num3;
                let a2 = ((u1 * x2 * x3) as f64 * den1 - (u2 * x3 * x1) as f64 * den2
                    + (u3 * x1 * x2) as f64 * den3) as f32;
                let b2 = ((-u1 * (x2 + x3)) as f64 * den1 + (u2 * (x3 + x1)) as f64 * den2
                    - (u3 * (x2 + x1)) as f64 * den3) as f32;
                let c2 = (u1 as f64 * den1 - u2 as f64 * den2 + u3 as f64 * den3) as f32;
                for it in (iti + 1).max(0)..=itf.min(ng as i32 - 1) {
                    itmin = itmin.min(it);
                    itmax = itmax.max(it);
                    let xl = it as f32 * dg;
                    let xl2 = xl * xl;
                    let t1 = a1 + b1 * xl + c1 * xl2;
                    let t2 = a2 + b2 * xl + c2 * xl2;
                    *arrival_mut(ag, it as usize, im[it as usize] as usize) =
                        (0.5 * (t1 + t2) as f64) as f32;
                    im[it as usize] += 1;
                }
                a1 = a2;
                b1 = b2;
                c1 = c2;
                ixl = ix as i32;
            }
        }
        if sub {
            iti = itf;
        }
        u0 = u1;
        x0 = x1;
        x1 = x2;
        u1 = u2;
    }
    for it in 0..itmin as usize {
        ag[it][0] = 0.0;
        im[it] += 1;
    }
    for it in (itmax + 1) as usize..ng {
        ag[it][0] = 0.0;
        im[it] += 1;
    }
}

// ealloc2float stores the three-column rows in one contiguous allocation.
// Arrivals beyond column two in C spill into the following row.
fn arrival(ag: &[[f32; 3]], row: usize, slot: usize) -> f32 {
    ag[row + slot / 3][slot % 3]
}
fn arrival_mut(ag: &mut [[f32; 3]], row: usize, slot: usize) -> &mut f32 {
    &mut ag[row + slot / 3][slot % 3]
}
