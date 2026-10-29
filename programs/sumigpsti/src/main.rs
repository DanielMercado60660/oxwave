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

//! SUMIGPSTI - Migration by phase shift for TI media with turning rays.
//! Port of su/main/migration_inversion/sumigpsti.c (Revision 1.4).
//! Credits: CWP: Dave Hale (originally call supsmig.c);
//! modified to TI media by Tariq Alkhalifah.

use oxwave_core::prelude::*;
use oxwave_cwp::{
    complex::Complex,
    cwp::{abs, max, min, nint, PI},
};
use oxwave_index_search::xindex;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfao, npfaro, pfa2cr, pfa2rc, pfacc};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGPSTI - MIGration by Phase Shift for TI media with turning rays\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumigpsti <stdin >stdout [optional parms]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " \tNone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt=from header(dt) or .004\ttime sampling interval\t\t\t",
    " dx=from header(d2) or 1.0\tdistance between sucessive cdp's\t",
    " ffil=0,0,0.5/dt,0.5/dt  trapezoidal window of frequencies to migrate\t",
    " tmig=0.0\ttimes corresponding to interval velocities in vmig\t",
    " vnmig=1500.0\tinterval NMO velocities corresponding to times in tmig\t",
    " vmig=1500.0\tinterval velocities corresponding to times in tmig\t",
    " etamig=0.0\tinterval eta values corresponding to times in tmig\t",
    " vnfile=\tbinary (non-ascii) file containing NMO velocities vn(t)\t",
    " vfile=\tbinary (non-ascii) file containing velocities v(t)\t",
    " etafile=\tbinary (non-ascii) file containing eta values eta(t)\t",
    " nxpad=0\tnumber of cdps to pad with zeros before FFT\t\t",
    " ltaper=0\tlength of linear taper for left and right edges\t\t",
    " verbose=0\t=1 for diagnostic print\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces must be sorted by either increasing or decreasing cdp.\t",
    "\t\t\t\t\t\t\t\t\t",
    " The tmig, vnmig, vmig and etamig arrays specify an interval values\t",
    " function of time. Linear interpolation and constant extrapolation is\t",
    " used to determine interval velocities at times not specified.  Values\t",
    " specified in tmig must increase monotonically.\t\t\t",
    " Alternatively, interval velocities may be stored in a binary file\t",
    " containing one velocity for every time sample.  If vnfile is specified,",
    " then the tmig and vnmig arrays are ignored.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The time of first sample is assumed to be zero, regardless of the value",
    " of the trace header field delrt.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed:  ns and dt\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

struct TATable {
    np: usize,
    p: Vec<f32>,
    taumax: Vec<f32>,
    taumin: Vec<f32>,
    tturn: Vec<f32>,
    t: Vec<Vec<f32>>,
    a: Vec<Vec<f32>>,
}

fn model_file(name: &str) -> CFile {
    if let Some(fp) = fopen(name, "r") {
        return fp;
    }
    // C passes the unchecked fopen result directly to fread. A null FILE
    // terminates with SIGSEGV and no diagnostic, rather than a Rust panic.
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        // Rust installs a SIGSEGV handler for stack overflows; restore SIG_DFL.
        signal(11, 0);
        raise(11);
    }
    exit(139);
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut dt = 0.0f32;
    let mut dx = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f32 as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            warn!("tr.dt not set, assuming dt=0.004");
        }
    }
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }
    let mut ffil = vec![0.0f32; 4];
    let nfil = getparfloat("ffil", &mut ffil);
    if nfil == 0 {
        ffil[2] = (0.5 / dt as f64) as f32;
        ffil[3] = ffil[2];
    } else if nfil != 4 {
        err!("if ffil specified, exactly 4 values must be provided");
    }
    let mut nxpad = 0i32;
    let mut ltaper = 0i32;
    let mut verbose = 0i32;
    let mut np = 200i32;
    let mut ntflag = 1i32;
    getparint("nxpad", &mut nxpad);
    getparint("ltaper", &mut ltaper);
    getparint("verbose", &mut verbose);
    getparint("np", &mut np);
    getparint("ntflag", &mut ntflag);

    /* determine NMO velocity function vnm(t) */
    let mut vnt = vec![0.0f32; nt];
    let mut vnfile = String::new();
    if getparstring("vnfile", &mut vnfile) == 0 {
        let mut ntmig = countparval("tmig");
        if ntmig == 0 {
            ntmig = 1;
        }
        let mut tmig = vec![0.0f32; ntmig as usize];
        getparfloat("tmig", &mut tmig);
        let mut nvnmig = countparval("vnmig");
        if nvnmig == 0 {
            nvnmig = 1;
        }
        if nvnmig != ntmig {
            err!("number of tmig and vnmig must be equal");
        }
        let mut vnmig = vec![0.0f32; nvnmig as usize];
        if getparfloat("vnmig", &mut vnmig) == 0 {
            vnmig[0] = 1500.0;
        }
        for itmig in 1..ntmig as usize {
            if tmig[itmig] <= tmig[itmig - 1] {
                err!("tmig must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..nt {
            intlin(
                &tmig,
                &vnmig,
                vnmig[0],
                vnmig[nvnmig as usize - 1],
                &[t],
                &mut vnt[it..it + 1],
            );
            t += dt;
        }
    } else {
        let fp = model_file(&vnfile);
        if fread(&mut vnt, &fp) != nt {
            err!("cannot read %d velocities from file %s", nt, vnfile);
        }
    }

    /* determine function eta(t) */
    let mut etat = vec![0.0f32; nt];
    let mut etafile = String::new();
    if getparstring("etafile", &mut etafile) == 0 {
        let mut ntmig = countparval("tmig");
        if ntmig == 0 {
            ntmig = 1;
        }
        let mut tmig = vec![0.0f32; ntmig as usize];
        getparfloat("tmig", &mut tmig);
        let mut netamig = countparval("etamig");
        if netamig == 0 {
            netamig = 1;
        }
        if netamig != ntmig {
            err!("number of tmig and etamig must be equal");
        }
        let mut etamig = vec![0.0f32; netamig as usize];
        getparfloat("etamig", &mut etamig);
        for itmig in 1..ntmig as usize {
            if tmig[itmig] <= tmig[itmig - 1] {
                err!("tmig must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..nt {
            intlin(
                &tmig,
                &etamig,
                etamig[0],
                etamig[netamig as usize - 1],
                &[t],
                &mut etat[it..it + 1],
            );
            t += dt;
        }
    } else {
        let fp = model_file(&etafile);
        if fread(&mut etat, &fp) != nt {
            err!("cannot read %d velocities from file %s", nt, etafile);
        }
    }

    /* determine velocity function v(t) */
    let mut vt = vec![0.0f32; nt];
    let mut vfile = String::new();
    if getparstring("vfile", &mut vfile) == 0 {
        let mut ntmig = countparval("tmig");
        if ntmig == 0 {
            ntmig = 1;
        }
        let mut tmig = vec![0.0f32; ntmig as usize];
        getparfloat("tmig", &mut tmig);
        let nvmig = countparval("vmig");
        if nvmig != ntmig {
            for it in 0..nt {
                vt[it] = vnt[it];
            }
        } else {
            let mut vmig = vec![0.0f32; nvmig as usize];
            if getparfloat("vmig", &mut vmig) == 0 {
                vmig[0] = vnt[0];
            }
            for itmig in 1..ntmig as usize {
                if tmig[itmig] <= tmig[itmig - 1] {
                    err!("tmig must increase monotonically");
                }
            }
            let mut t = 0.0f32;
            for it in 0..nt {
                intlin(
                    &tmig,
                    &vmig,
                    vmig[0],
                    vmig[nvmig as usize - 1],
                    &[t],
                    &mut vt[it..it + 1],
                );
                t += dt;
            }
        }
    } else {
        let fp = model_file(&vfile);
        if fread(&mut vt, &fp) != nt {
            err!("cannot read %d velocities from file %s", nt, vfile);
        }
    }
    checkpars();

    // The C table and interpolation require two slopes; smaller values
    // otherwise produce invalid allocations or out-of-bounds accesses.
    if np < 2 {
        err!("np must be at least 2");
    }

    // tableta's extrapolation reads beyond nt for the shear branch of the
    // quadratic. Deliberately back those otherwise undefined reads with zeros;
    // the image table has already reached its last tau sample on these rays.
    let mut a1111 = vec![0.0f32; 2 * nt + 1];
    let mut a3333 = vec![0.0f32; 2 * nt + 1];
    let mut a1313 = vec![0.0f32; 2 * nt + 1];
    let mut a1133 = vec![0.0f32; 2 * nt + 1];
    for it in 0..nt {
        let vnt2 = vnt[it] * vnt[it];
        let vt2 = vt[it] * vt[it];
        a1111[it] = vnt2 * (1.0 + 2.0 * etat[it]);
        a3333[it] = vt2;
        a1313[it] = (0.25 * vt2 as f64) as f32;
        let delta = (0.5 * (vnt2 / vt2 - 1.0) as f64) as f32;
        // sqrt returns double, and subtraction occurs before narrowing.
        a1133[it] = (((2.0 * delta * a3333[it] * (a3333[it] - a1313[it])
            + (a3333[it] - a1313[it]) * (a3333[it] - a1313[it])) as f64)
            .sqrt()
            - a1313[it] as f64) as f32;
    }

    /* copy traces and headers to temporary files */
    let tfp = etmpfile();
    let hfp = etmpfile();
    let mut nx = 0usize;
    loop {
        nx += 1;
        fwrite(tr.hdr.as_bytes(), &hfp);
        fwrite(&tr.data[..nt], &tfp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    efseeko(&hfp, 0, SEEK_SET);
    efseeko(&tfp, 0, SEEK_SET);
    if verbose != 0 {
        fprintf!(stderr(), "\t%d traces input\n", nx);
    }
    // C overflows these bounds for extreme padding, and may select a
    // transform too short for the input when padding is negative.
    let nmin = i32::try_from(nx).ok().and_then(|nx| nx.checked_add(nxpad));
    let Some((nmin, nmax)) = nmin.and_then(|n| n.checked_mul(2).map(|m| (n, m))) else {
        err!("nxpad overflows FFT size");
    };
    let nxfft = npfaro(nmin, nmax);
    if nxfft < 0 || (nxfft as usize) < nx {
        err!("nxpad produces an FFT length smaller than the number of input traces");
    }
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let mut gtk = vec![Complex::default(); nt * nk as usize];
    let mut gtx = vec![0.0f32; nt * nxfft as usize];
    let fftscl = (1.0 / nxfft as f64) as f32;
    for ix in 0..nx {
        efread(&mut gtx[ix * nt..(ix + 1) * nt], &tfp);
        for it in 0..nt {
            gtx[ix * nt + it] *= fftscl;
        }
        if (ix as i32) < ltaper {
            let taper = (ix + 1) as f32 / (ltaper + 1) as f32;
            for it in 0..nt {
                gtx[ix * nt + it] *= taper;
            }
        } else if ix as i32 >= nx as i32 - ltaper {
            let taper = (nx - ix) as f32 / (ltaper + 1) as f32;
            for it in 0..nt {
                gtx[ix * nt + it] *= taper;
            }
        }
    }
    pfa2rc(-1, 2, nt as i32, nxfft, &gtx, &mut gtk);
    if verbose != 0 {
        fprintf!(stderr(), "\tFourier transform done\n");
    }
    let table = tableta(
        np as usize,
        nt,
        dt,
        &mut a3333,
        &mut a1111,
        &mut a1133,
        &mut a1313,
        nt,
        dt,
        verbose,
    );
    if verbose != 0 {
        fprintf!(stderr(), "\tTime/amplitude table built\n");
    }
    let trig = TrigTable::new();
    let mut k = 0.0f32;
    for ik in 0..nk {
        if verbose != 0 && ik % max(nk / 10, 1) == 0 {
            fprintf!(stderr(), "\t%d of %d wavenumbers done\n", ik, nk);
        }
        mig1k(
            &table,
            &trig,
            k,
            &ffil,
            ntflag,
            nt,
            dt,
            &mut gtk[ik as usize * nt..(ik as usize + 1) * nt],
            nt,
            dt,
        );
        k += dk;
    }
    pfa2cr(1, 2, nt as i32, nxfft, &gtk, &mut gtx);
    if verbose != 0 {
        fprintf!(stderr(), "\tinverse Fourier transform done\n");
    }
    for ix in 0..nx {
        efread(tr.hdr.as_bytes_mut(), &hfp);
        tr.data[..nt].copy_from_slice(&gtx[ix * nt..(ix + 1) * nt]);
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

struct TrigTable {
    ctab: Vec<f32>,
    stab: Vec<f32>,
    fntab: f32,
    opi2: f32,
}
impl TrigTable {
    fn new() -> Self {
        let ntab = 1025;
        let dangle = (2.0 * PI / (ntab - 1) as f64) as f32;
        let mut ctab = vec![0.0; ntab];
        let mut stab = vec![0.0; ntab];
        let mut angle = 0.0f32;
        for itab in 0..ntab {
            ctab[itab] = (angle as f64).cos() as f32;
            stab[itab] = (angle as f64).sin() as f32;
            angle += dangle;
        }
        Self {
            ctab,
            stab,
            fntab: ntab as f32,
            opi2: (1.0 / (PI * 2.0)) as f32,
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn mig1k(
    table: &TATable,
    trig: &TrigTable,
    k: f32,
    ffil: &[f32],
    ntflag: i32,
    nt: usize,
    dt: f32,
    cp: &mut [Complex],
    ntau: usize,
    dtau: f32,
) {
    let np = table.np;
    let p = &table.p;
    let ntfft = npfao(2 * nt as i32, 4 * nt as i32);
    let nw = ntfft;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let mut pw = vec![Complex::default(); ntfft as usize];
    let mut qn = vec![Complex::default(); ntau];
    let mut qt = vec![Complex::default(); ntau];
    let wnyq = (PI / dt as f64) as f32;
    let iwnyq = if ntfft % 2 != 0 {
        ntfft / 2
    } else {
        ntfft / 2 + 1
    };
    let iwfil0 = max(0, min(iwnyq, nint(2.0 * PI * ffil[0] as f64 / dw as f64)));
    let iwfil1 = max(0, min(iwnyq, nint(2.0 * PI * ffil[1] as f64 / dw as f64)));
    let iwfil2 = max(0, min(iwnyq, nint(2.0 * PI * ffil[2] as f64 / dw as f64)));
    let iwfil3 = max(0, min(iwnyq, nint(2.0 * PI * ffil[3] as f64 / dw as f64)));
    for it in 0..nt {
        pw[it] = cp[it];
    }
    pfacc(1, &mut pw);
    for iw in 0..nw as usize {
        pw[iw].r = (pw[iw].r as f64 * (1.0 / ntfft as f64)) as f32;
        pw[iw].i = (pw[iw].i as f64 * (1.0 / ntfft as f64)) as f32;
    }
    let wmin = abs(k) / p[np - 1];
    let iwmin = max(max(1, iwfil0), (wmin / dw) as i32);
    let wmax = if p[0] <= abs(k) / wnyq {
        wnyq
    } else {
        abs(k) / p[0]
    };
    let iwmax = min(min(iwnyq - 1, iwfil3), (wmax / dw) as i32);
    let mut w = iwmin as f32 * dw;
    // The original leaves ip1 uninitialised; xindex establishes the bracket.
    let mut ip1 = 0i32;
    for iw in iwmin..=iwmax {
        let kow = abs(k) / w;
        if kow > p[np - 1] {
            w += dw;
            continue;
        }
        let ampf = if iwfil0 <= iw && iw < iwfil1 {
            (iw - iwfil0) as f32 / (iwfil1 - iwfil0) as f32
        } else if iwfil2 < iw && iw <= iwfil3 {
            (iwfil3 - iw) as f32 / (iwfil3 - iwfil2) as f32
        } else {
            1.0
        };
        xindex(p, kow, &mut ip1);
        ip1 = min(ip1, np as i32 - 2);
        let ip2 = ip1 as usize + 1;
        let ip1u = ip1 as usize;
        let p1 = p[ip1u];
        let p2 = p[ip2];
        let a1 = (p2 * p2 - kow * kow) / (p2 * p2 - p1 * p1);
        let a2 = (kow * kow - p1 * p1) / (p2 * p2 - p1 * p1);
        let wa1 = w * a1;
        let wa2 = w * a2;
        let taumaxi = a1 * table.taumax[ip1u] + a2 * table.taumax[ip2];
        let itaumax = min(ntau as i32 - 1, nint((taumaxi / dtau) as f64));
        let taumini = a1 * table.taumin[ip1u] + a2 * table.taumin[ip2];
        let itaumin = max(0, nint((taumini / dtau) as f64));
        let phst = (2.0 * (wa1 * table.tturn[ip1u] + wa2 * table.tturn[ip2]) as f64) as f32;
        let pwsr = ampf * (pw[iw as usize].r + pw[(nw - iw) as usize].r);
        let pwsi = ampf * (pw[iw as usize].i + pw[(nw - iw) as usize].i);
        let pwdr = ampf * (pw[iw as usize].r - pw[(nw - iw) as usize].r);
        let pwdi = ampf * (pw[iw as usize].i - pw[(nw - iw) as usize].i);
        if ntflag & 1 != 0 {
            for itau in 0..itaumax {
                let itau = itau as usize;
                let ampi = a1 * table.a[ip1u][itau] + a2 * table.a[ip2][itau];
                let phsi = wa1 * table.t[ip1u][itau] + wa2 * table.t[ip2][itau];
                let phsn = phsi * trig.opi2;
                let itab = (trig.fntab * (phsn - phsn as i32 as f32)) as usize;
                qn[itau].r += ampi * (pwsr * trig.ctab[itab] + pwdi * trig.stab[itab]);
                qn[itau].i += ampi * (pwsi * trig.ctab[itab] - pwdr * trig.stab[itab]);
            }
        }
        if ntflag & 2 != 0 {
            for itau in itaumin + 1..itaumax {
                let itau = itau as usize;
                let ampi = a1 * table.a[ip1u][itau] + a2 * table.a[ip2][itau];
                let phsi = phst - (wa1 * table.t[ip1u][itau] + wa2 * table.t[ip2][itau]);
                let phsn = phsi * trig.opi2;
                let itab = (trig.fntab * (phsn - phsn as i32 as f32)) as usize;
                qt[itau].r += ampi * (pwsr * trig.stab[itab] - pwdi * trig.ctab[itab]);
                qt[itau].i += ampi * (pwsi * trig.stab[itab] + pwdr * trig.ctab[itab]);
            }
        }
        w += dw;
    }
    for itau in 0..ntau {
        cp[itau].r = qn[itau].r + qt[itau].r;
        cp[itau].i = qn[itau].i + qt[itau].i;
    }
}

#[allow(clippy::too_many_arguments)]
fn tableta(
    np: usize,
    ntau: usize,
    dtau: f32,
    a3333: &mut [f32],
    a1111: &mut [f32],
    a1133: &mut [f32],
    a1313: &mut [f32],
    nt: usize,
    dt: f32,
    _verbose: i32,
) -> TATable {
    // The C allocation is uninitialised; invalid slopes leave some entries
    // unset. Zeroes make those entries deterministic.
    let mut table = TATable {
        np,
        p: vec![0.0; np],
        taumax: vec![0.0; np],
        taumin: vec![0.0; np],
        tturn: vec![0.0; np],
        t: vec![vec![0.0; ntau]; np],
        a: vec![vec![0.0; ntau]; np],
    };
    for it in 0..nt {
        a1111[it] = (0.25 * a1111[it] as f64) as f32;
        a3333[it] = (0.25 * a3333[it] as f64) as f32;
        a1313[it] = (0.25 * a1313[it] as f64) as f32;
        a1133[it] = (0.25 * a1133[it] as f64) as f32;
    }
    for jp in 0..np {
        let pj = (2.0 / (a1111[0] as f64).sqrt()
            * ((jp as f32 / (np as i32 - 1) as f32) as f64).sqrt()) as f32;
        table.p[jp] = pj;
        table.t[jp][0] = 0.0;
        table.a[jp][0] = 1.0;
        let mut ktau = 0usize;
        let mut taui = 0.0f32;
        let mut ti = 0.0f32;
        let mut ai = 1.0f32;
        let mut a1111t = a1111[0];
        let mut a3333t = a3333[0];
        let mut a1313t = a1313[0];
        let mut a1133t = a1133[0];
        let mut f1 = a1111t + a1313t;
        let mut f2 = a3333t + a1313t;
        let mut f3 = a1111t - a1313t;
        let mut f4 = a3333t - a1313t;
        let mut f5 = (2.0 * (a1133t + a1313t) as f64 * (a1133t + a1313t) as f64) as f32;
        let f6 = 2.0f32;
        let eps = 0.0001f32;
        let px2 = pj * pj;
        let mut alpha = f2 * f2 - f4 * f4;
        let mut beta = 2.0 * ((f1 * f2 + f3 * f4 - f5) * px2 - f2 * f6);
        // C evaluates (f1*f1-f3*f3)*px2 in float before the double
        // subtraction introduced by 2.; retain that rounding at both steps.
        let mut gamma = ((f6 * f6) as f64
            - (2.0 * f1 as f64 * f6 as f64 - ((f1 * f1 - f3 * f3) * px2) as f64) * px2 as f64)
            as f32;
        let mut det = ((beta * beta) as f64 - 4.0 * alpha as f64 * gamma as f64) as f32;
        if det < 0.0 {
            continue;
        }
        let mut rad = (det as f64).sqrt() as f32;
        let mut signbeta = if abs(beta) > eps {
            abs(beta) / beta
        } else {
            1.0
        };
        let mut q = (-0.5 * (beta + signbeta * rad) as f64) as f32;
        let mut pz2 = gamma / q;
        if pz2 < 0.0 {
            continue;
        }
        let mut vp = (1.0 / ((px2 + pz2) as f64).sqrt()) as f32;
        let mut angle = min(1.0f64, (pj * vp) as f64).asin() as f32;
        for _it in 1..nt {
            let taul = taui;
            let tl = ti;
            let al = ai;
            let cosa = ((angle as f64).cos() * (a3333t as f64).sqrt() / vp as f64) as f32;
            taui += dt * cosa;
            ti += dt * cosa * cosa;
            ai = 1.0;
            if taui < 0.0 {
                break;
            }
            if taui >= taul {
                table.tturn[jp] = ti;
                table.taumax[jp] = taui;
                table.taumin[jp] = taui;
            } else {
                table.taumin[jp] = taui;
            }
            let itau = (taui / dtau) as i32;
            let ltau = (taul / dtau) as i32;
            for jtau in ltau + 1..=min(itau, ntau as i32 - 1) {
                let tauj = jtau as f32 * dtau;
                let frac = (tauj - taul) / (taui - taul);
                ktau += 1;
                table.t[jp][ktau] = ((1.0 - frac as f64) * tl as f64 + (frac * ti) as f64) as f32;
                table.a[jp][ktau] = ((1.0 - frac as f64) * al as f64 + (frac * ai) as f64) as f32;
            }
            if itau < ntau as i32 - 1 {
                let frac = (taui - itau as f32 * dtau) / dtau;
                let itau = itau as usize;
                a1111t = ((1.0 - frac as f64) * a1111[itau] as f64
                    + (frac * a1111[itau + 1]) as f64) as f32;
                a3333t = ((1.0 - frac as f64) * a3333[itau] as f64
                    + (frac * a3333[itau + 1]) as f64) as f32;
                a1313t = ((1.0 - frac as f64) * a1313[itau] as f64
                    + (frac * a1313[itau + 1]) as f64) as f32;
                a1133t = ((1.0 - frac as f64) * a1133[itau] as f64
                    + (frac * a1133[itau + 1]) as f64) as f32;
            } else {
                let itau = itau as usize;
                a1111t = a1111[ntau - 1]
                    + (taui - (ntau - 1) as f32 * dtau) * (a1111[itau] - a1111[itau - 1]) / dtau;
                a3333t = a3333[ntau - 1]
                    + (taui - (ntau - 1) as f32 * dtau) * (a3333[itau] - a3333[itau - 1]) / dtau;
                a1313t = a1313[ntau - 1]
                    + (taui - (ntau - 1) as f32 * dtau) * (a1313[itau] - a1313[itau - 1]) / dtau;
                a1133t = a1133[ntau - 1]
                    + (taui - (ntau - 1) as f32 * dtau) * (a1133[itau] - a1133[itau - 1]) / dtau;
            }
            f1 = a1111t + a1313t;
            f2 = a3333t + a1313t;
            f3 = a1111t - a1313t;
            f4 = a3333t - a1313t;
            f5 = (2.0 * (a1133t + a1313t) as f64 * (a1133t + a1313t) as f64) as f32;
            alpha = f2 * f2 - f4 * f4;
            beta = 2.0 * ((f1 * f2 + f3 * f4 - f5) * px2 - f2 * f6);
            gamma = ((f6 * f6) as f64
                - (2.0 * f1 as f64 * f6 as f64 - ((f1 * f1 - f3 * f3) * px2) as f64)
                    * px2 as f64) as f32;
            det = ((beta * beta) as f64 - 4.0 * alpha as f64 * gamma as f64) as f32;
            if det < 0.0 {
                continue;
            }
            rad = (det as f64).sqrt() as f32;
            signbeta = if abs(beta) > eps {
                abs(beta) / beta
            } else {
                1.0
            };
            q = (-0.5 * (beta + signbeta * rad) as f64) as f32;
            pz2 = gamma / q;
            if pz2 < 0.0 {
                continue;
            }
            vp = (1.0 / ((px2 + pz2) as f64).sqrt()) as f32;
            angle = min(1.0f64, (pj * vp) as f64).asin() as f32;
        }
        for jtau in ktau + 1..ntau {
            table.t[jp][jtau] = table.t[jp][ktau];
            table.a[jp][jtau] = table.a[jp][ktau];
        }
    }
    table
}
