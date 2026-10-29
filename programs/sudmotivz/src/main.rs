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

//! SUDMOTIVZ - DMO for transversely isotropic V(Z) media, common-offset gathers.
//! Port of su/main/dip_moveout/sudmotivz.c (Revision 1.4).
//! Credits (from the C source): CWP: Tariq Alkhalifah, Craig Artley.

#![allow(non_snake_case, clippy::too_many_arguments)]

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{PI, abs, nint};
use oxwave_inverse_interp::*;
use oxwave_linear_interp::*;
use oxwave_lu::*;
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::*;
mod ray;
use ray::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDMOTIVZ - DMO for Transeversely Isotropic V(Z) media for common-offset",
    "            gathers\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sudmotivz <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " cdpmin         minimum cdp (integer number) for which to apply DMO\t",
    " cdpmax         maximum cdp (integer number) for which to apply DMO\t",
    " dxcdp          distance between adjacent cdp bins (m)\t\t\t",
    " noffmix        number of offsets to mix (see notes)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " vnfile=        binary (non-ascii) file containing NMO interval\t",
    "                  velocities (m/s)\t\t\t\t\t",
    " vfile=         binary (non-ascii) file containing interval velocities\t(m/s)",
    " etafile=       binary (non-ascii) file containing eta interval values (m/s)",
    " tdmo=0.0       times corresponding to interval velocities in vdmo (s)\t",
    " vndmo=1500.0   NMO interval velocities corresponding to times in tdmo (m/s)",
    " vdmo=vndmo    interval velocities corresponding to times in tdmo (m/s)",
    " etadmo=1500.0  eta interval values corresponding to times in tdmo (m/s)",
    " fmax=0.5/dt    maximum frequency in input traces (Hz)\t\t\t",
    " smute=1.5      stretch mute used for NMO correction\t\t\t",
    " speed=1.0      twist this knob for speed (and aliasing)\t\t",
    " verbose=0      =1 for diagnostic print\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces should be sorted into common-offset gathers.  One common-",
    " offset gather ends and another begins when the offset field of the trace",
    " headers changes.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
    " the cdp location expressed in units of meters or feet.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of offsets to mix (noffmix) should typically equal the ratio of",
    " the shotpoint spacing to the cdp spacing.  This choice ensures that every",
    " cdp will be represented in each offset mix.  Traces in each mix will\t",
    " contribute through DMO to other traces in adjacent cdps within that mix.",
    "\t\t\t\t\t\t\t\t\t",
    " vnfile, vfile and etafile should contain the regularly sampled interval",
    " values of NMO velocity, velocity and eta respectivily as a\t\t",
    " function of time.  If, for example, vfile is not supplied, the interval",
    " velocity function is defined by linear interpolation of the values in the",
    " tdmo and vdmo arrays.  The times in tdmo must be monotonically increasing.",
    " If vfile or vdmo are not given it will be equated to vnfile or vndmo. ",
    "\t\t\t\t\t\t\t\t\t",
    " For each offset, the minimum time to process is determined using the\t",
    " smute parameter.  The DMO correction is not computed for samples that\t",
    " have experienced greater stretch during NMO.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed:  nt, dt, delrt, offset, cdp.\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    // The C call omits the %f argument; this build prints zero.
    if ft != 0.0 {
        err!("ft=%f --- nonzero ft not yet implemented!", 0.0f64);
    }
    let mut offset = tr.offset as f32;

    let (mut cdpmin, mut cdpmax, mut noffmix) = (0i32, 0i32, 0i32);
    let (mut dx, mut fmax, mut speed, mut smute) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let mut verbose = 0i32;
    if getparint("cdpmin", &mut cdpmin) == 0 {
        err!("must specify cdpmin");
    }
    if getparint("cdpmax", &mut cdpmax) == 0 {
        err!("must specify cdpmax");
    }
    if cdpmin > cdpmax {
        err!("cdpmin must not be greater than cdpmax");
    }
    if getparfloat("dxcdp", &mut dx) == 0 {
        err!("must specify dxcdp");
    }
    if getparint("noffmix", &mut noffmix) == 0 {
        err!("must specify noffmix");
    }
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (0.5 / dt as f64) as f32;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    if getparfloat("speed", &mut speed) == 0 {
        speed = 1.0;
    }
    if getparfloat("smute", &mut smute) == 0 {
        smute = 1.5;
    }

    let nttab = 2 * nt - 1;
    let mut vt = vec![0.0f32; nttab];
    let mut vnt = vec![0.0f32; nttab];
    let mut etat = vec![0.0f32; nttab];
    let mut a1111 = vec![0.0f32; nttab];
    let mut a3333 = vec![0.0f32; nttab];
    let mut a1313 = vec![0.0f32; nttab];
    let mut a1133 = vec![0.0f32; nttab];
    let mut vrms = vec![0.0f32; nttab];
    let mut vnfile = String::new();
    let mut etafile = String::new();
    let mut vfile = String::new();

    /* determine NMO velocity function vnm(t) */
    if getparstring("vnfile", &mut vnfile) == 0 {
        let mut ntdmo = countparval("tdmo") as usize;
        if ntdmo == 0 {
            ntdmo = 1;
        }
        let mut tdmo = vec![0.0f32; ntdmo];
        if getparfloat("tdmo", &mut tdmo) == 0 {
            tdmo[0] = 0.0;
        }
        let mut nvndmo = countparval("vndmo") as usize;
        if nvndmo == 0 {
            nvndmo = 1;
        }
        if nvndmo != ntdmo {
            err!("number of tdmo and vndmo must be equal");
        }
        let mut vndmo = vec![0.0f32; nvndmo];
        if getparfloat("vndmo", &mut vndmo) == 0 {
            vndmo[0] = 1500.0;
        }
        for itdmo in 1..ntdmo {
            if tdmo[itdmo] <= tdmo[itdmo - 1] {
                err!("tdmo must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..nt {
            intlin(
                &tdmo,
                &vndmo,
                vndmo[0],
                vndmo[ntdmo - 1],
                &[t],
                &mut vnt[it..it + 1],
            );
            t += dt;
        }
    } else {
        read_model(&vnfile, &mut vnt[..nt]);
    }
    for it in nt..nttab {
        vnt[it] = vnt[nt - 1];
    }

    /* determine eta function eta(t) */
    if getparstring("etafile", &mut etafile) == 0 {
        let mut ntdmo = countparval("tdmo") as usize;
        if ntdmo == 0 {
            ntdmo = 1;
        }
        let mut tdmo = vec![0.0f32; ntdmo];
        if getparfloat("tdmo", &mut tdmo) == 0 {
            tdmo[0] = 0.0;
        }
        let mut netadmo = countparval("etadmo") as usize;
        if netadmo == 0 {
            netadmo = 1;
        }
        if netadmo != ntdmo {
            err!("number of tdmo and etadmo must be equal");
        }
        let mut etadmo = vec![0.0f32; netadmo];
        if getparfloat("etadmo", &mut etadmo) == 0 {
            etadmo[0] = 0.0;
        }
        for itdmo in 1..ntdmo {
            if tdmo[itdmo] <= tdmo[itdmo - 1] {
                err!("tdmo must increase monotonically");
            }
        }
        let mut t = 0.0f32;
        for it in 0..nt {
            intlin(
                &tdmo,
                &etadmo,
                etadmo[0],
                etadmo[ntdmo - 1],
                &[t],
                &mut etat[it..it + 1],
            );
            t += dt;
        }
    } else {
        read_model(&etafile, &mut etat[..nt]);
    }
    for it in nt..nttab {
        etat[it] = etat[nt - 1];
    }

    /* determine velocity function v(t) */
    if getparstring("vfile", &mut vfile) == 0 {
        let mut ntdmo = countparval("tdmo") as usize;
        if ntdmo == 0 {
            ntdmo = 1;
        }
        let mut tdmo = vec![0.0f32; ntdmo];
        if getparfloat("tdmo", &mut tdmo) == 0 {
            tdmo[0] = 0.0;
        }
        let nvdmo = countparval("vdmo") as usize;
        if nvdmo != ntdmo {
            for it in 0..nt {
                vt[it] = vnt[it];
            }
        } else {
            let mut vdmo = vec![0.0f32; nvdmo];
            if getparfloat("vdmo", &mut vdmo) == 0 {
                vdmo[0] = vnt[0];
            }
            for itdmo in 1..ntdmo {
                if tdmo[itdmo] <= tdmo[itdmo - 1] {
                    err!("tdmo must increase monotonically");
                }
            }
            let mut t = 0.0f32;
            for it in 0..nt {
                intlin(
                    &tdmo,
                    &vdmo,
                    vdmo[0],
                    vdmo[ntdmo - 1],
                    &[t],
                    &mut vt[it..it + 1],
                );
                t += dt;
            }
        }
    } else {
        read_model(&vfile, &mut vt[..nt]);
    }
    checkpars();
    for it in nt..nttab {
        vt[it] = vt[nt - 1];
    }

    /* compute vrms */
    vrms[0] = vnt[0];
    let mut sum = 0.0f32;
    for it in 1..nttab {
        sum += vnt[it - 1] * vnt[it - 1];
        vrms[it] = ((sum / it as f32) as f64).sqrt() as f32;
    }
    /* convert eta and NMO velocity values to elastic coefficients */
    for it in 0..nttab {
        let vnt2 = vnt[it] * vnt[it];
        let vt2 = vt[it] * vt[it];
        a1111[it] = vnt2 * (1.0f32 + 2.0f32 * etat[it]);
        a3333[it] = vt2;
        a1313[it] = (0.25 * vt2 as f64) as f32;
        let delta = (0.5 * ((vnt2 / vt2 - 1.0f32) as f64)) as f32;
        a1133[it] = (((2.0f32 * delta * a3333[it] * (a3333[it] - a1313[it])
            + (a3333[it] - a1313[it]) * (a3333[it] - a1313[it])) as f64)
            .sqrt()
            - a1313[it] as f64) as f32;
    }
    /* compute half-velocities */
    for it in 0..nttab {
        a1111[it] = (a1111[it] as f64 * 0.25) as f32;
        a3333[it] = (a3333[it] as f64 * 0.25) as f32;
        a1313[it] = (a1313[it] as f64 * 0.25) as f32;
        a1133[it] = (a1133[it] as f64 * 0.25) as f32;
        vt[it] = (vt[it] as f64 * 0.5) as f32;
    }
    let ft = 0.0f32;
    let fp = 0.0f32;
    let lp = (1.0 / (a1111[0] as f64).sqrt()) as f32;
    let nptab = 401usize;
    let dptab = lp / (nptab - 1) as f32;
    let mut tau = vec![vec![0.0f32; nttab]; nptab];
    let mut x = vec![vec![0.0f32; nttab]; nptab];
    let mut a = vec![vec![0.0f32; nttab]; nptab];
    let mut s = vec![vec![0.0f32; nttab]; nptab];
    if verbose != 0 {
        fprintf!(stderr(), "Computing ray tables...\n");
    }
    let mut p = fp;
    for ip in 0..nptab {
        rayvt(
            p,
            nttab as i32,
            dt,
            &a1111,
            &a3333,
            &a1313,
            &a1133,
            &mut tau[ip],
            &mut x[ip],
            &mut a[ip],
            &mut s[ip],
        );
        p += dptab;
    }
    if verbose != 0 {
        fprintf!(stderr(), "Ray tables ready...\n");
    }
    drop((vt, vnt, etat, a1111, a3333, a1133, a1313));

    let nx = (cdpmax - cdpmin + 1) as usize;
    let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
    let nxfft = npfar(nx as i32 + nxpad) as usize;
    let mut q = vec![0.0f32; (nxfft + 2) * nt];
    let mut qmix = vec![vec![0.0f32; nt]; nx];
    let hfp = etmpfile();
    let mut oldoffset = offset;
    let mut gottrace = true;
    let mut ntrace = 0;
    let mut noff = 0;
    loop {
        if gottrace {
            offset = tr.offset as f32;
        }
        if (gottrace && offset != oldoffset) || !gottrace {
            dmooff(
                &x, &a, &tau, &s, nttab, dt, ft, nptab, dptab, fp, lp, nx, dx, oldoffset, nt,
                &vrms, fmax, smute, speed, verbose, &mut q,
            );
            for ix in 0..nx {
                for it in 0..nt {
                    qmix[ix][it] += q[ix * nt + it];
                }
            }
            noff += 1;
            if offset != oldoffset {
                let nxpad = (0.5 * abs(offset / dx) as f64) as i32;
                let nxfft = npfar(nx as i32 + nxpad) as usize;
                q = vec![0.0f32; (nxfft + 2) * nt];
            }
        }
        if noff == noffmix || !gottrace {
            efseeko(&hfp, 0, SEEK_SET);
            for _itrace in 0..ntrace {
                efread(tro.hdr.as_bytes_mut(), &hfp);
                let ix = (tro.cdp - cdpmin) as usize;
                tro.data[..nt].copy_from_slice(&qmix[ix]);
                puttr(&mut tro);
            }
            if verbose != 0 {
                fprintf!(
                    stderr(),
                    "\tCompleted mix of %d offsets with %d traces\n",
                    noff,
                    ntrace
                );
            }
            if !gottrace {
                break;
            }
            efseeko(&hfp, 0, SEEK_SET);
            noff = 0;
            ntrace = 0;
            for ix in 0..nx {
                for it in 0..nt {
                    qmix[ix][it] = 0.0;
                }
            }
        }
        if tr.cdp >= cdpmin && tr.cdp <= cdpmax {
            efwrite(tr.hdr.as_bytes(), &hfp);
            ntrace += 1;
            oldoffset = offset;
            let ix = (tr.cdp - cdpmin) as usize;
            q[ix * nt..(ix + 1) * nt].copy_from_slice(&tr.data[..nt]);
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    exit(CWP_Exit());
}

fn read_model(file: &str, values: &mut [f32]) {
    use oxwave_core::cio::{fopen, fread};
    // The original passes a failed fopen directly to fread and receives SIGSEGV.
    let fp = match fopen(file, "r") {
        Some(fp) => fp,
        None => {
            unsafe extern "C" {
                fn raise(sig: i32) -> i32;
            }
            unsafe {
                raise(11);
            }
            exit(139);
        }
    };
    if fread(values, &fp) != values.len() {
        err!("cannot read %d velocities from file %s", values.len(), file);
    }
}

fn tmute(vrms: &[f32], nt: usize, dt: f32, ft: f32, h: f32, smute: f32) -> f32 {
    let mut t = vec![0.0f32; nt];
    let mut tn = ft;
    for it in 0..nt {
        t[it] = ((tn * tn) as f64 + 4.0 * h as f64 * h as f64 / (vrms[it] * vrms[it]) as f64).sqrt()
            as f32;
        tn += dt;
    }
    let mut it = nt - 1;
    while it > 0 {
        if dt / (t[it] - t[it - 1]) > smute {
            break;
        }
        it -= 1;
    }
    ft + it as f32 * dt
}

fn dmooff(
    x: &[Vec<f32>],
    a: &[Vec<f32>],
    tau: &[Vec<f32>],
    s: &[Vec<f32>],
    nttab: usize,
    dt: f32,
    ft: f32,
    nptab: usize,
    dptab: f32,
    fp: f32,
    lp: f32,
    nx: usize,
    dx: f32,
    offset: f32,
    nt: usize,
    vrms: &[f32],
    fmax: f32,
    smute: f32,
    speed: f32,
    verbose: i32,
    q: &mut [f32],
) {
    let tmin = tmute(vrms, nt, dt, ft, (0.5 * offset as f64) as f32, smute);
    let lt = (nt - 1) as f32 * dt;
    if tmin >= lt {
        err!("tmin=%f >= lt=%f:  no data survive the mute!\n", tmin, lt);
    }
    let dtn = (10.0 * dt as f64) as f32;
    let ftn = tmin;
    let ntn = (nint(((lt - ftn) / dtn) as f64) + 1) as usize;
    let mut dp2 = (tmin as f64 / (fmax as f64 * 0.25 * offset as f64 * offset as f64)) as f32;
    dp2 *= speed;
    let np = nint((lp * lp / dp2) as f64).wrapping_add(1);
    if np < 0 {
        // The signed dimension becomes a huge size_t in ealloc1int.
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let np = np as usize;
    let mut it0min = vec![0i32; np];
    let mut nt0 = vec![0i32; np];
    // One trailing cell preserves the original allocation's look-ahead
    // when a slope has no positive mapping entries (itnmin == ntn).
    let mut t0table = vec![vec![0.0f32; ntn + 1]; np];
    let mut tntable = vec![vec![0.0f32; nt]; np];
    let h = abs((offset as f64 / 2.0) as f32);
    let mut tn = ftn;
    for itn in 0..ntn {
        dmomap(
            x,
            a,
            tau,
            s,
            nttab,
            dt,
            ft,
            nptab,
            dptab,
            np,
            dp2,
            fp,
            lp,
            itn,
            tn,
            h,
            vrms,
            verbose,
            &mut t0table,
        );
        tn += dtn;
    }
    let mut ip = 0;
    while ip < np {
        let mut itnmin = 0usize;
        while itnmin < ntn {
            if t0table[ip][itnmin] > 0.0 {
                break;
            }
            itnmin += 1;
        }
        it0min[ip] = nint((t0table[ip][itnmin] / dt) as f64);
        if it0min[ip] >= nt as i32 {
            it0min[ip] = nt as i32 - 1;
        }
        let mut itnmax = ntn as i32 - 1;
        while itnmax >= 0 {
            if t0table[ip][itnmax as usize] > 0.0 {
                break;
            }
            itnmax -= 1;
        }
        if itnmax <= itnmin as i32 {
            break;
        }
        let itnmax = itnmax as usize;
        let mut it0max = nint((t0table[ip][itnmax] / dt) as f64);
        if it0max >= nt as i32 {
            it0max = nt as i32 - 1;
        }
        nt0[ip] = it0max - it0min[ip] + 1;
        for itn in itnmin + 1..itnmax {
            if t0table[ip][itn] <= t0table[ip][itn - 1] && verbose != 0 {
                fprintf!(
                    stderr(),
                    "t0table not monotonically increasing---extrapolating\n"
                );
            }
        }
        let nout = nt0[ip] as usize;
        yxtoxy(
            dtn,
            ftn + itnmin as f32 * dtn,
            &t0table[ip][itnmin..itnmax + 1],
            dt,
            it0min[ip] as f32 * dt,
            ftn + itnmin as f32 * dtn,
            ftn + itnmax as f32 * dtn,
            &mut tntable[ip][..nout],
        );
        for it in nout..nt {
            // ealloc2float stores rows contiguously. A one-sample inverse
            // mapping reads the previous row's last sample at index -1.
            // Before the first row, the reference allocation has zero bytes.
            let prev = if it >= 1 {
                tntable[ip][it - 1]
            } else if ip > 0 {
                tntable[ip - 1][nt - 1]
            } else {
                0.0
            };
            let prev2 = if it >= 2 {
                tntable[ip][it - 2]
            } else if ip > 0 {
                tntable[ip - 1][nt + it - 2]
            } else {
                0.0
            };
            tntable[ip][it] = (2.0 * prev as f64 - prev2 as f64) as f32;
        }
        nt0[ip] = nt as i32 - it0min[ip];
        ip += 1;
    }
    let npdmo = ip.min(np);
    drop(t0table);
    dmoapply(
        nx, dx, nt, dt, npdmo, dp2, h, &tntable, &it0min, &nt0, tmin, fmax, q,
    );
}

fn dmomap(
    x: &[Vec<f32>],
    a: &[Vec<f32>],
    tau: &[Vec<f32>],
    s: &[Vec<f32>],
    nt: usize,
    dt: f32,
    ft: f32,
    nptab: usize,
    dptab: f32,
    np: usize,
    dp2: f32,
    fp: f32,
    lp: f32,
    itn: usize,
    tn: f32,
    h: f32,
    vrms: &[f32],
    _verbose: i32,
    t0table: &mut [Vec<f32>],
) {
    let lt = ft + (nt - 1) as f32 * dt;
    let mut vars = vec![0.0f32; 5];
    let mut eqns = vec![0.0f32; 5];
    let mut jacob = vec![vec![0.0f32; 5]; 5];
    let mut work = vec![0.0f32; 5];
    let mut ipvt = vec![0i32; 5];
    let vtn = vrms[(tn / dt) as usize];
    let tsg = ((tn * tn + 4.0f32 * h * h / (vtn * vtn)) as f64).sqrt() as f32;
    getX(tsg, h, fp, vtn, &mut vars);
    let mut p2 = fp;
    // C leaves rcond unset if the initial trial converges without a solve.
    let mut rcond = 1.0f32;
    for ip in 0..np {
        let p0 = (p2 as f64).sqrt() as f32;
        for _iter in 0..10 {
            getFJ(
                x,
                a,
                tau,
                s,
                nt as i32,
                nptab as i32,
                dt,
                dptab,
                tsg,
                h,
                p0,
                &vars,
                &mut eqns,
                &mut jacob,
            );
            let mut error = 0.0f32;
            for i in 0..5 {
                error += abs(eqns[i]);
            }
            let tol = 1.0e-4f32;
            if error <= tol {
                break;
            }
            {
                let mut rows: Vec<&mut [f32]> = jacob.iter_mut().map(|r| &mut r[..]).collect();
                sgeco(&mut rows, 5, &mut ipvt, &mut rcond, &mut work);
            }
            let rows: Vec<&[f32]> = jacob.iter().map(|r| &r[..]).collect();
            sgesl(&rows, 5, &ipvt, &mut eqns, 0);
            if rcond != 0.0 {
                for i in 0..5 {
                    vars[i] -= eqns[i];
                }
            } else {
                break;
            }
            if abs(vars[1]) > lp {
                vars[1] = (oxwave_cwp::cwp::sgn(vars[1]) * lp as f64) as f32;
            }
            if abs(vars[2]) > lp {
                vars[2] = (oxwave_cwp::cwp::sgn(vars[2]) * lp as f64) as f32;
            }
            if vars[3] < 0.0 {
                vars[3] = 0.0;
            }
            if vars[3] > tsg {
                vars[3] = tsg;
            }
            if vars[4] < 0.0 {
                vars[4] = 0.0;
            }
            if vars[4] > lt {
                vars[4] = lt;
            }
        }
        let (mut tautp, mut dtaudt, mut dtaudp) = (0.0f32, 0.0f32, 0.0f32);
        intder(
            tau,
            nt as i32,
            dt,
            ft,
            nptab as i32,
            dptab,
            fp,
            vars[4],
            abs(p0),
            &mut tautp,
            &mut dtaudt,
            &mut dtaudp,
        );
        if tautp < dt {
            break;
        }
        if rcond == 0.0 {
            break;
        }
        let x0 = vars[0];
        let t0 = vars[4];
        t0table[ip][itn] = t0 + p0 * x0;
        p2 += dp2;
    }
}

fn dmoapply(
    nx: usize,
    dx: f32,
    nt: usize,
    dt: f32,
    np: usize,
    dp2: f32,
    h: f32,
    tntable: &[Vec<f32>],
    it0min: &[i32],
    nt0: &[i32],
    tmin: f32,
    fmax: f32,
    qx: &mut [f32],
) {
    let nxfft = npfar(nx as i32) as usize;
    let nk = nxfft / 2 + 1;
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    let mut qk = vec![Complex::default(); nk * nt];
    pfa2rc(-1, 2, nt as i32, nxfft as i32, qx, &mut qk);
    let mut k = 0.0f32;
    for ik in 0..nk {
        jakub1k(
            k,
            h,
            np,
            dp2,
            nt,
            dt,
            0.0,
            tntable,
            it0min,
            nt0,
            fmax,
            &mut qk[ik * nt..(ik + 1) * nt],
        );
        k += dk;
    }
    pfa2cr(1, 2, nt as i32, nxfft as i32, &qk, qx);
    let scale = (1.0 / nxfft as f64) as f32;
    let itmin = (nt as i32).min(nint((tmin / dt) as f64)) as usize;
    for ix in 0..nx {
        for it in 0..itmin {
            qx[ix * nt + it] = 0.0;
        }
        for it in itmin..nt {
            qx[ix * nt + it] *= scale;
        }
    }
}

fn jakub1k(
    k: f32,
    _h: f32,
    np: usize,
    dp2: f32,
    nt: usize,
    dt: f32,
    ft: f32,
    tntable: &[Vec<f32>],
    it0min: &[i32],
    nt0: &[i32],
    fmax: f32,
    qq: &mut [Complex],
) {
    let ntfft = npfa(nt as i32) as usize;
    let nw = ntfft as i32;
    let dw = (2.0 * PI / (ntfft as f32 * dt) as f64) as f32;
    let fw = (-PI / dt as f64) as f32;
    let czero = Complex::default();
    let mut qn = vec![czero; ntfft];
    let mut q0 = vec![czero; ntfft];
    let mut p = 0.0f32;
    let mut pl = 0.0f32;
    let mut ph2 = dp2;
    for ip in 0..np {
        let ph = (ph2 as f64).sqrt() as f32;
        let wl = k / ph;
        let wh = if pl == 0.0 {
            (2.0 * PI * fmax as f64) as f32
        } else {
            let l = (k / pl) as f64;
            let r = 2.0 * PI * fmax as f64;
            (if l < r { l } else { r }) as f32
        };
        let iwlm = ((-wl - fw) / dw) as i32;
        let mut iwhm = ((-wh - fw) / dw) as i32;
        if iwhm < 0 {
            iwhm = 0;
        }
        iwhm += 1;
        let mut iwlp = ((wl - fw) / dw) as i32;
        let mut iwhp = ((wh - fw) / dw) as i32;
        if iwhp >= nw {
            iwhp = nw - 1;
        }
        iwlp += 1;
        if iwhm <= iwlm || iwlp <= iwhp {
            let start = it0min[ip] as usize;
            let count = nt0[ip] as usize;
            for it in 0..start {
                qn[it] = czero;
            }
            ints8c(
                dt,
                ft,
                qq,
                czero,
                czero,
                &tntable[ip][..count],
                &mut qn[start..start + count],
            );
            for it in start + count..ntfft {
                qn[it] = czero;
            }
            for it in (1..ntfft).step_by(2) {
                qn[it].r = -qn[it].r;
                qn[it].i = -qn[it].i;
            }
            pfacc(1, &mut qn);
            for iw in iwhm..=iwlm {
                q0[iw as usize].r += qn[iw as usize].r;
                q0[iw as usize].i += qn[iw as usize].i;
            }
            for iw in iwlp..=iwhp {
                q0[iw as usize].r += qn[iw as usize].r;
                q0[iw as usize].i += qn[iw as usize].i;
            }
        }
        pl = p;
        p = ph;
        ph2 += dp2;
    }
    pfacc(-1, &mut q0);
    let scale = (0.5 / ntfft as f64) as f32;
    for it in (0..nt).step_by(2) {
        qq[it].r = q0[it].r * scale;
        qq[it].i = q0[it].i * scale;
    }
    for it in (1..nt).step_by(2) {
        qq[it].r = -q0[it].r * scale;
        qq[it].i = -q0[it].i * scale;
    }
}
