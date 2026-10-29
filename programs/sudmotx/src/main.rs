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

//! SUDMOTX - DMO via T-X domain (Kirchhoff) method for common-offset gathers.
//! Port of su/main/dip_moveout/sudmotx.c (Revision 1.17).
//! Credits: CWP: Dave Hale.
//! Technical reference: A non-aliased integral method for dip-moveout,
//! Dave Hale, submitted to Geophysics, June, 1990.

use oxwave_core::prelude::*;
use oxwave_correlation::xcor;
use oxwave_cwp::cwp::{abs, max, min, nint, PI};
use oxwave_tempfile::temporary_filename;
use oxwave_toeplitz::stoepf;
use std::ffi::CString;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUDMOTX - DMO via T-X domain (Kirchhoff) method for common-offset gathers",
    "\t\t\t\t\t\t\t\t\t",
    " sudmotx <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [optional parms]",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " cdpmin                  minimum cdp (integer number) for which to apply DMO",
    " cdpmax                  maximum cdp (integer number) for which to apply DMO",
    " dxcdp                   distance between successive cdps\t\t",
    " noffmix                 number of offsets to mix (see notes)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " offmax=3000.0           maximum offset\t\t\t\t",
    " tmute=2.0               mute time at maximum offset offmax\t\t",
    " vrms=1500.0             RMS velocity at mute time tmute\t\t",
    " verbose=0               =1 for diagnostic print\t\t\t",
    " tmpdir=\tif non-empty, use the value as a directory path\tprefix\t",
    "\t\tfor storing temporary files; else if the CWP_TMPDIR\t",
    "\t\tenvironment variable is set use\tits value for the path;\t",
    "\t\telse use tmpfile()\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Input traces should be sorted into common-offset gathers.  One common-",
    " offset gather ends and another begins when the offset field of the trace",
    " headers changes.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
    " the cdp location expressed in units of meters or feet.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of offsets to mix (noffmix) should typically equal the ratio of",
    " the shotpoint spacing to the cdp spacing.  This choice ensures that every",
    " cdp will be represented in each offset mix.  Traces in each mix will\t",
    " contribute through DMO to other traces in adjacent cdps within that mix.",
    "\t\t\t\t\t\t\t\t\t",
    " The defaults for offmax and vrms are appropriate only for metric units.",
    " If distances are measured in feet, then these parameters should be\t",
    " specified explicitly.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " offmax, tmute, and vrms need not be specified precisely.\t\t",
    " If these values are unknown, then one should overestimate offmax\t",
    " and underestimate tmute and vrms.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " No muting is actually performed.  The tmute parameter is used only to\t",
    " determine parameters required to perform DMO.\t\t\t\t",
];

// Globals used only by the original program's graceful signal cleanup.
static mut HEADERFP: Option<CFile> = None;
static mut HEADERFILE: String = String::new();

unsafe extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get information from the first header.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;

    let mut cdpmin = 0i32;
    let mut cdpmax = 0i32;
    let mut dxcdp = 0.0f32;
    let mut noffmix = 0i32;
    let mut offmax = 0.0f32;
    let mut tmute = 0.0f32;
    let mut vrms = 0.0f32;
    let mut nsmax = 0i32;
    let mut verbose = 0i32;
    if getparint("cdpmin", &mut cdpmin) == 0 {
        err!("must specify cdpmin");
    }
    if getparint("cdpmax", &mut cdpmax) == 0 {
        err!("must specify cdpmax");
    }
    if cdpmin > cdpmax {
        err!("cdpmin must be less than cdpmax");
    }
    if getparfloat("dxcdp", &mut dxcdp) == 0 {
        err!("must specify dxcdp");
    }
    if getparint("noffmix", &mut noffmix) == 0 {
        err!("must specify noffmix");
    }
    if getparfloat("offmax", &mut offmax) == 0 {
        offmax = 3000.0;
    }
    if getparfloat("tmute", &mut tmute) == 0 {
        tmute = 2.0;
    }
    if getparfloat("vrms", &mut vrms) == 0 {
        vrms = 1500.0;
    }
    if getparint("nsmax", &mut nsmax) == 0 {
        nsmax = 400;
    }
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_str()).unwrap();
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();

    let ncdp = cdpmax.wrapping_sub(cdpmin).wrapping_add(1);
    let mut q = vec![vec![0.0f32; nt]; ncdp as usize];
    let mut p = vec![0.0f32; nt];
    let mut temp = vec![0.0f32; nt];
    // A negative C int converts to a huge size_t and ealloc fails.
    if nsmax < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // malloc(0) in the reference build supplies a readable zeroed slot.
    // makeds and dmotx access as[0] even when no shifts were requested.
    let mut ts = vec![0.0f32; max(nsmax, 1) as usize];
    let mut amps = vec![0.0f32; max(nsmax, 1) as usize];
    let lds = 125;
    let ifds = -100;
    let mut ds = vec![0.0f32; lds];

    let ns = maketa(dxcdp, dt, offmax, tmute, vrms, nsmax, &mut ts, &mut amps);
    if verbose != 0 {
        fprintf!(stderr(), "\tDMO will be performed via %d time shifts\n", ns);
    }
    makeds(ns, &ts, &amps, ifds, &mut ds);
    // xcor reverses its first input before calling convolution. Reverse ds
    // here and adjust its origin to invoke that same library calculation.
    ds.reverse();
    let reversed_ifds = 1 - ifds - lds as i32;

    let headerfile;
    let headerfp;
    if tmpdir.is_empty() {
        headerfile = String::new();
        headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        headerfile = temporary_filename(&tmpdir);
        unsafe {
            HEADERFILE = headerfile.clone();
            for sig in [2, 3, 1, 15] {
                signal(sig, closefiles);
            }
        }
        headerfp = efopen(&headerfile, "w+");
        unsafe {
            HEADERFP = Some(headerfp.clone());
        }
        if verbose != 0 {
            warn!("putting temporary header file in %s", tmpdir);
        }
    }

    let mut oldoffset = tr.offset as f32;
    let mut gottrace = true;
    let mut ntrace = 0;
    let mut noff = 0;
    let mut offset = 0.0f32;
    let mut cdp = 0i32;
    loop {
        if gottrace {
            offset = tr.offset as f32;
            cdp = tr.cdp;
            if offset != oldoffset {
                noff += 1;
            }
            p.copy_from_slice(&tr.data[..nt]);
        }
        if noff == noffmix || !gottrace {
            if !gottrace {
                noff += 1;
            }
            for icdp in 0..ncdp as usize {
                xcor(reversed_ifds, &ds, 0, &q[icdp], 0, &mut temp);
                q[icdp].copy_from_slice(&temp);
            }
            erewind(&headerfp);
            for _itrace in 0..ntrace {
                efread(tro.hdr.as_bytes_mut(), &headerfp);
                let icdp = tro.cdp.wrapping_sub(cdpmin);
                tro.data[..nt].copy_from_slice(&q[icdp as usize]);
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
            erewind(&headerfp);
            noff = 0;
            ntrace = 0;
            for icdp in 0..ncdp as usize {
                for it in 0..nt {
                    q[icdp][it] = 0.0;
                }
            }
        }
        if cdp >= cdpmin && cdp <= cdpmax {
            efwrite(tr.hdr.as_bytes(), &headerfp);
            ntrace += 1;
            let icdp = cdp.wrapping_sub(cdpmin);
            let jcdplo = max(0.0, icdp as f64 - 0.5 * abs(offset / dxcdp) as f64) as i32;
            let jcdphi = min(
                (ncdp - 1) as f64,
                icdp as f64 + 0.5 * abs(offset / dxcdp) as f64,
            ) as i32;
            for jcdp in jcdplo..=jcdphi {
                dmotx(
                    ns,
                    &ts,
                    &amps,
                    offset,
                    (jcdp - icdp) as f32 * dxcdp,
                    dxcdp,
                    0,
                    nt as i32,
                    dt,
                    ft,
                    &p,
                    &mut q[jcdp as usize],
                );
            }
            oldoffset = offset;
        }
        if gettr(&mut tr) == 0 {
            gottrace = false;
        }
    }
    efclose(&headerfp);
    if !headerfile.is_empty() {
        eremove(&headerfile);
    }
    exit(CWP_Exit());
}

/// Time shifts and amplitudes for non-aliased (t,x) domain DMO.
/// Dave Hale, Colorado School of Mines, 11/04/90.
fn maketa(
    dx: f32,
    dt: f32,
    offmax: f32,
    tmute: f32,
    vrms: f32,
    nsmax: i32,
    ts: &mut [f32],
    amps: &mut [f32],
) -> i32 {
    let dx = abs(dx);
    let dt = abs(dt);
    let offmax = abs(offmax);
    let tmute = abs(tmute);
    let vrms = abs(vrms);
    let dalphan = min(
        (2.0 * (tmute / dt) as f64 * ((dx / offmax) as f64).powf(2.0)).sqrt(),
        1.0,
    ) as f32;
    let tsmax = ((tmute / dt) as f64
        * (1.0 - 1.0 / (1.0 + ((offmax / (vrms * tmute)) as f64).powf(2.0)).sqrt()))
        as f32;
    let mut alpha = 0.0f32;
    let mut dalpha = dalphan;
    let mut nsi = 0;
    while alpha * alpha <= tsmax && nsi < nsmax {
        ts[nsi as usize] = alpha * alpha;
        amps[nsi as usize] = dalpha;
        alpha = alpha + dalpha;
        dalpha = (((alpha * alpha) as f64 + 1.0).sqrt() - alpha as f64) as f32;
        if dalpha > dalphan {
            dalpha = dalphan;
        }
        nsi += 1;
    }
    let ltaper = nsi / 3;
    let mut itaper = 1;
    for is in nsi - ltaper..nsi {
        amps[is as usize] = (amps[is as usize] as f64
            * (0.54 + 0.46 * (PI * itaper as f64 / ltaper as f64).cos()))
            as f32;
        itaper += 1;
    }
    nsi
}

/// Shaping filter correcting the DMO horizontal reflection response.
/// Dave Hale, Colorado School of Mines, 11/04/90.
fn makeds(ns: i32, ts: &[f32], amps: &[f32], ifds: i32, ds: &mut [f32]) {
    let ld = 400;
    let ifd = 1 - ld;
    let mut d = [0.0f32; 400];
    let mut dd = [0.0f32; 300];
    let mut di = [0.0f32; 300];
    let mut work = [0.0f32; 300];
    d[ld as usize - 1] = amps[0];
    for is in 1..ns as usize {
        let i = ts[is] as i32;
        d[(ld - i - 1) as usize] = (d[(ld - i - 1) as usize] as f64 + 2.0 * amps[is] as f64) as f32;
    }
    let lds = ds.len();
    xcor(ifd, &d, ifd, &d, 0, &mut dd[..lds]);
    xcor(ifd, &d, 0, &[1.0], ifds, &mut di[..lds]);
    stoepf(&dd[..lds], &di[..lds], ds, &mut work[..lds]);
}

/// Apply DMO to one input trace and accumulate into one output trace.
/// Dave Hale, Colorado School of Mines, 11/04/90.
#[allow(clippy::too_many_arguments)]
fn dmotx(
    ns: i32,
    ts: &[f32],
    amps: &[f32],
    offset: f32,
    x: f32,
    dx: f32,
    itmute: i32,
    nt: i32,
    dt: f32,
    ft: f32,
    p: &[f32],
    q: &mut [f32],
) {
    let h = (offset as f64 / 2.0) as f32;
    let xb = abs(nint((x / dx) as f64) as f32 * dx);
    let ftodt = ft / dt;
    if xb == 0.0 {
        let x1 = (xb as f64 + 0.5 * abs(dx) as f64) as f32;
        let scale1 =
            (1.0 / (1.0 - max(0.0, 1.0 - (x1 * x1 / max(h * h, x1 * x1)) as f64).sqrt())) as f32;
        let asi = amps[0];
        for it in itmute..nt {
            q[it as usize] += asi * p[it as usize];
        }
        for is in 1..ns as usize {
            let tsi = ts[is];
            let asi = amps[is];
            let its = tsi as i32;
            let mut it1 = (scale1 * tsi - ftodt) as i32;
            if it1 < itmute + its {
                it1 = itmute + its;
            }
            if it1 >= nt {
                break;
            }
            for it in it1..nt {
                q[(it - its) as usize] = (q[(it - its) as usize] as f64
                    + 2.0 * asi as f64 * p[it as usize] as f64)
                    as f32;
            }
        }
    } else if ns > 1 {
        let x1 = (xb as f64 + 0.5 * abs(dx) as f64) as f32;
        let x2 = (xb as f64 - 0.5 * abs(dx) as f64) as f32;
        if x2 * x2 > h * h {
            return;
        }
        let scale1 =
            (1.0 / (1.0 - max(0.0, 1.0 - (x1 * x1 / max(h * h, x1 * x1)) as f64).sqrt())) as f32;
        let scale2 =
            (1.0 / (1.0 - max(0.0, 1.0 - (x2 * x2 / max(h * h, x2 * x2)) as f64).sqrt())) as f32;
        for is in 1..ns as usize {
            let tsi = ts[is];
            let asi = amps[is];
            let its = tsi as i32;
            let mut it1 = (scale1 * tsi - ftodt) as i32;
            if it1 < itmute + its {
                it1 = itmute + its;
            }
            let mut it2 = (scale2 * tsi - ftodt) as i32;
            if it2 > nt {
                it2 = nt;
            }
            if it1 >= nt {
                break;
            }
            for it in it1..it2 {
                q[(it - its) as usize] += asi * p[it as usize];
            }
        }
    }
}

// Like the C handler, this runs only on the single main thread and exits.
#[allow(static_mut_refs)]
extern "C" fn closefiles(_sig: i32) {
    unsafe {
        if let Some(fp) = HEADERFP.as_ref() {
            efclose(fp);
        }
        eremove(&HEADERFILE);
    }
    exit(EXIT_FAILURE);
}
