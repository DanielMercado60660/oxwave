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

//! SUINTERP - interpolate traces using automatic event picking.
//! Port of su/main/interp_extrap/suinterp.c (Revision 1.15).
//! Credit: John Anderson (visiting scholar from Mobil), July 1994.
use oxwave_core::cio::as_bytes_mut;
use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_f_vnd::*;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::ints8r;
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, c_int, CStr, CString};
use std::sync::OnceLock;

extern "C" {
    fn access(path: *const c_char, mode: c_int) -> c_int;
    fn signal(sig: c_int, handler: extern "C" fn(c_int)) -> usize;
}

static HEADERFILE: OnceLock<String> = OnceLock::new();

extern "C" fn closefiles(_: c_int) {
    if let Some(file) = HEADERFILE.get() {
        eremove(file);
    }
    // exit flushes the registered CFile streams before process teardown.
    exit(1);
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUINTERP - interpolate traces using automatic event picking\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "           suinterp < stdin > stdout\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " ninterp=1    number of traces to output between each pair of input traces",
    " nxmax=500    maximum number of input traces\t\t\t\t",
    " freq1=4.     starting corner frequency of unaliased range\t\t",
    " freq2=20.    ending corner frequency of unaliased range\t\t",
    " deriv=0      =1 means take vertical derivative on pick section        ",
    "              (useful if interpolating velocities instead of seismic)  ",
    " linear=0     =0 means use 8 point sinc temporal interpolation         ",
    "              =1 means use linear temporal interpolation               ",
    "              (useful if interpolating velocities instead of seismic)  ",
    " lent=5       number of time samples to smooth for dip estimate\t",
    " lenx=1       number of traces to smooth for dip estimate\t\t",
    " lagc=400     number of ms agc for dip estimate\t\t\t",
    " xopt=0       0 compute spatial derivative via FFT\t\t\t",
    "                 (assumes input traces regularly spaced and relatively\t",
    "                  noise-free)\t\t\t\t\t\t",
    "              1 compute spatial derivative via differences\t\t",
    "                 (will work on irregulary spaced data)\t\t\t",
    " iopt=0     0 = interpolate",
    "            1 = output low-pass model: useful for QC if interpolator failing",
    "            2 = output dip picks in units of samples/trace\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0\tverbose = 1 echoes information\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program outputs 'ninterp' interpolated traces between each pair of",
    " input traces.  The values for lagc, freq1, and freq2 are only used for",
    " event tracking. The output data will be full bandwidth with no agc.  The",
    " default parameters typically will do a satisfactory job of interpolation",
    " for dips up to about 12 ms/trace.  Using a larger value for freq2 causes",
    " the algorithm to do a better job on the shallow dips, but to fail on the",
    " steep dips.  Only one dip is assumed at each time sample between each pair",
    " of input traces.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The key assumption used here is that the low frequency data are unaliased",
    " and can be used for event tracking. Those dip picks are used to interpolate",
    " the original full-bandwidth data, giving some measure of interpolation",
    " at higher frequencies which otherwise would be aliased.  Using iopt equal",
    " to 1 allows you to visually check whether the low-pass picking model is",
    " aliased.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Trace headers for interpolated traces are not updated correctly.\t",
    " The output header for an interpolated traces equals that for the preceding",
    " trace in the original input data.  The original input traces are passed",
    " through this module without modification.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The place this code is most likely to fail is on the first breaks.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example run:    suplane | suinterp | suxwigb &\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// VNDtempname owns a VND allocation, which must be freed with VNDfree.
fn v2dop(lwmax: i64, nt: usize, nx: i64) -> VND {
    if nt == 0 || nx <= 0 {
        // Nonpositive dimensions make the C VND geometry undefined.
        err!("invalid VND dimensions");
    }
    let file = VNDtempname("suinterp");
    // SAFETY: VNDtempname returns a live, NUL-terminated allocation.
    let name = unsafe { CStr::from_ptr(file.cast()) }
        .to_string_lossy()
        .into_owned();
    let vnd = V2Dop(2, lwmax, 4, &name, nt as i64, nx);
    unsafe {
        VNDfree(file, "file");
    }
    vnd.unwrap_or_else(|| err!("can't create VND file"))
}

fn runav(a: &[f32], len: i32, b: &mut [f32]) {
    let n = a.len() as i32;
    let lenh = len / 2;
    if len <= 1 {
        b[..a.len()].copy_from_slice(a);
        return;
    }
    let mut sum = 0.0f32;
    for j in 0..len.min(n) {
        sum += a[j as usize];
    }
    for j in 0..(lenh + 1).min(n) {
        b[j as usize] = sum;
    }
    for j in lenh + 1..n - lenh {
        sum = sum + a[(j + lenh) as usize] - a[(j - lenh - 1) as usize];
        b[j as usize] = sum;
    }
    for j in (n - lenh).max(0)..n {
        b[j as usize] = sum;
    }
    sum = (1.0 / len as f64) as f32;
    for j in 0..a.len() {
        b[j] *= sum;
    }
}

fn check_sinc_time(time: &[f32]) {
    // intt8r converts (8 + time) to int, then subtracts 11. Outside
    // this range C has undefined conversions/indexing; reject before the
    // Rust library can panic on an invalid sinc-table index.
    for &t in time {
        let x = (8.0f32 + t) as f64;
        if !x.is_finite() || x < (i32::MIN + 11) as f64 || x > i32::MAX as f64 {
            err!("interpolation time exceeds sinc index range");
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn jea_xinterpolate(
    vndorig: &mut VND,
    vndinterp: &mut VND,
    ninterp: i32,
    nt: usize,
    n2: usize,
    freq1: f32,
    freq2: f32,
    lagc: i32,
    lent: i32,
    len2: i32,
    mut xopt: i32,
    dt: f32,
    iopt: i32,
    deriv: i32,
    linear: i32,
) {
    let lent = 1i32.wrapping_add(2i32.wrapping_mul(lent / 2));
    let len2 = 1i32.wrapping_add(2i32.wrapping_mul(len2 / 2));
    let lagc = (1.0 + lagc as f64 * 0.001 / dt as f64) as i32;
    let ntfft = npfar(nt as i32) as usize;
    let n2fft = npfar(n2 as i32) as usize;
    if ntfft < nt || n2fft < n2 {
        err!("dimensions exceed FFT limits");
    }
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    // C aliases rt with crt. Separate work arrays preserve the transforms
    // while retaining only O(max(nt,n2)) samples in memory.
    let mut rt = vec![0.0f32; ntfft.max(n2fft) + 2];
    let mut crt = vec![Complex::default(); (ntfft / 2 + 1).max(n2fft / 2 + 1)];
    if n2 < 2 || (iopt == 0 && ninterp == 0) {
        for ix in 0..n2 {
            V2Dr0(vndorig, ix as i64, as_bytes_mut(&mut rt), 101);
            V2Dw0(vndinterp, ix as i64, as_bytes_mut(&mut rt), 102);
        }
        // C incorrectly free()s the guarded VND allocation on this path.
        return;
    }
    let mut ccrt = vec![Complex::default(); ntfft / 2 + 1];
    let mut rrt = vec![0.0; ntfft];
    let mut a = vec![0.0; n2.max(nt)];
    let mut b = a.clone();
    let mut p = vec![0.0; nt];
    let mut time = vec![0.0; nt];
    let mut aa = a.clone();
    let mut bb = a.clone();
    let tin: Vec<f32> = if linear != 0 {
        (0..nt).map(|j| j as f32).collect()
    } else {
        Vec::new()
    };
    let mut vnda = v2dop(500000, nt, n2 as i64);
    let mut vndb = v2dop(500000, nt, n2 as i64);

    // Compute filtered data for picking and its time derivative.
    let dff = (2.0 * std::f64::consts::PI / ntfft as f64) as f32;
    for ix in 0..n2 {
        V2Dr0(vndorig, ix as i64, as_bytes_mut(&mut rt), 103);
        for j in 0..nt {
            a[j] = rt[j].abs();
        }
        runav(&a[..nt], lagc, &mut b);
        runav(&b[..nt], lagc, &mut a);
        for j in 0..nt {
            rt[j] = rt[j] / (a[j] + 1.0e-30f32);
        }
        rt[nt..ntfft].fill(0.0);
        if deriv != 0 {
            for j in (1..nt).rev() {
                rt[j] -= rt[j - 1];
            }
            rt[0] = rt[1];
        }
        pfarc(1, &rt[..ntfft], &mut crt[..ntfft / 2 + 1]);
        for j in 0..ntfft / 2 + 1 {
            let mut f = j as f32 * df;
            if freq1 > 0.0 {
                let mut fcl = f / freq1;
                fcl = fcl * fcl * fcl * fcl;
                let mut fch = f / freq2;
                fch = fch * fch * fch * fch;
                f = (fcl as f64 / ((1.0 + fcl as f64) * (1.0 + fch as f64))) as f32;
            } else {
                let fch = f / freq2;
                f = (1.0 / (1.0 + (fch * fch * fch * fch) as f64)) as f32;
            }
            crt[j] = crmul(crt[j], f);
            ccrt[j] = cmul(crt[j], cmplx(0.0, -(j as f32) * dff));
        }
        pfacr(-1, &crt[..ntfft / 2 + 1], &mut rt[..ntfft]);
        V2Dw0(&mut vnda, ix as i64, as_bytes_mut(&mut rt), 104);
        pfacr(-1, &ccrt, &mut rrt);
        V2Dw0(&mut vndb, ix as i64, as_bytes_mut(&mut rrt), 105);
    }
    if iopt == 1 {
        for ix in 0..n2 {
            V2Dr0(&mut vnda, ix as i64, as_bytes_mut(&mut rt), 104);
            V2Dw0(vndinterp, ix as i64, as_bytes_mut(&mut rt), 104);
        }
        VNDcl(vnda, 1);
        VNDcl(vndb, 1);
        return;
    }

    // Compute the spatial derivative of the picking model.
    let dxx = (2.0 * std::f64::consts::PI / (n2fft * n2fft) as f64) as f32;
    if n2 < 16 {
        xopt = 1;
    }
    for it in 0..nt {
        V2Dr1(&mut vnda, it as i64, as_bytes_mut(&mut rt), 106);
        if xopt != 0 {
            for j in 0..n2 - 1 {
                rt[j] = rt[j + 1] - rt[j];
            }
            rt[n2 - 1] = rt[n2 - 2];
        } else {
            rt[n2..n2fft].fill(0.0);
            pfarc(1, &rt[..n2fft], &mut crt[..n2fft / 2 + 1]);
            for j in 0..n2fft / 2 + 1 {
                crt[j] = cmul(crt[j], cmplx(0.0, -(j as f32) * dxx));
            }
            pfacr(-1, &crt[..n2fft / 2 + 1], &mut rt[..n2fft]);
        }
        V2Dw1(&mut vnda, it as i64, as_bytes_mut(&mut rt), 107);
    }
    // Compute dot products and smooth over time.
    for ix in 0..n2 {
        V2Dr0(&mut vnda, ix as i64, as_bytes_mut(&mut a), 108);
        V2Dr0(&mut vndb, ix as i64, as_bytes_mut(&mut b), 109);
        for it in 0..nt {
            aa[it] = a[it] * b[it];
            bb[it] = b[it] * b[it];
        }
        runav(&aa[..nt], lent, &mut a);
        runav(&a[..nt], lent, &mut aa);
        runav(&bb[..nt], lent, &mut b);
        runav(&b[..nt], lent, &mut bb);
        V2Dw0(&mut vnda, ix as i64, as_bytes_mut(&mut aa), 110);
        V2Dw0(&mut vndb, ix as i64, as_bytes_mut(&mut bb), 111);
    }
    if len2 > 1 {
        for it in 0..nt {
            V2Dr1(&mut vnda, it as i64, as_bytes_mut(&mut a), 112);
            V2Dr1(&mut vndb, it as i64, as_bytes_mut(&mut b), 113);
            runav(&a[..n2], len2, &mut aa);
            runav(&aa[..n2], len2, &mut a);
            runav(&b[..n2], len2, &mut bb);
            runav(&bb[..n2], len2, &mut b);
            V2Dw1(&mut vnda, it as i64, as_bytes_mut(&mut a), 114);
            V2Dw1(&mut vndb, it as i64, as_bytes_mut(&mut b), 115);
        }
    }
    // Compute dip picks and interpolate one pair of original traces at a time.
    V2Dr0(vndorig, 0, as_bytes_mut(&mut a), 116);
    for ix in 1..n2 {
        let ixm = (ix - 1) as i64;
        V2Dr0(&mut vnda, ixm, as_bytes_mut(&mut aa), 117);
        V2Dr0(&mut vndb, ixm, as_bytes_mut(&mut bb), 118);
        for it in 0..nt {
            p[it] = -aa[it] / (bb[it] + 1.0e-30f32);
        }
        V2Dr0(vndorig, ix as i64, as_bytes_mut(&mut b), 119);
        if iopt == 2 {
            V2Dw0(vndinterp, ixm, as_bytes_mut(&mut p), 120);
        } else {
            V2Dw0(
                vndinterp,
                ixm * (ninterp as i64 + 1),
                as_bytes_mut(&mut a),
                120,
            );
            for k in 0..ninterp {
                let wa = ((1.0 + k as f64) / (1 + ninterp) as f64) as f32;
                let wb = (1.0 - wa as f64) as f32;
                for it in 0..nt {
                    time[it] = it as f32 - p[it] * wa;
                }
                if linear != 0 {
                    intlin(&tin, &a[..nt], a[0], a[nt - 1], &time, &mut aa[..nt]);
                } else {
                    check_sinc_time(&time);
                    ints8r(1.0, 0.0, &a[..nt], 0.0, 0.0, &time, &mut aa[..nt]);
                }
                for it in 0..nt {
                    time[it] = it as f32 + p[it] * wb;
                }
                if linear != 0 {
                    // The C source deliberately passes a here, too.
                    intlin(&tin, &a[..nt], a[0], a[nt - 1], &time, &mut bb[..nt]);
                } else {
                    check_sinc_time(&time);
                    ints8r(1.0, 0.0, &b[..nt], 0.0, 0.0, &time, &mut bb[..nt]);
                }
                for it in 0..nt {
                    aa[it] = wb * aa[it] + wa * bb[it];
                }
                V2Dw0(
                    vndinterp,
                    k as i64 + 1 + ixm * (ninterp as i64 + 1),
                    as_bytes_mut(&mut aa),
                    121,
                );
            }
        }
        std::mem::swap(&mut a, &mut b);
    }
    if iopt == 2 {
        V2Dw0(vndinterp, (n2 - 1) as i64, as_bytes_mut(&mut p), 122);
    } else {
        V2Dw0(
            vndinterp,
            (n2 - 1) as i64 * (ninterp as i64 + 1),
            as_bytes_mut(&mut a),
            122,
        );
    }
    VNDcl(vnda, 1);
    VNDcl(vndb, 1);
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut lagc = 400;
    let mut freq1 = 4.0;
    let mut freq2 = 20.0;
    let mut lent = 5;
    let mut lenx = 1;
    let mut nxmax = 500;
    let mut xopt = 0;
    let mut ninterp = 1i32;
    let mut iopt = 0;
    let mut deriv = 0;
    let mut linear = 0;
    let mut verbose = 0;
    getparint("lagc", &mut lagc);
    getparfloat("freq1", &mut freq1);
    getparfloat("freq2", &mut freq2);
    getparint("lent", &mut lent);
    getparint("lenx", &mut lenx);
    getparint("nxmax", &mut nxmax);
    getparint("xopt", &mut xopt);
    getparint("ninterp", &mut ninterp);
    getparint("iopt", &mut iopt);
    getparint("deriv", &mut deriv);
    getparint("linear", &mut linear);
    getparint("verbose", &mut verbose);
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.split('\0').next().unwrap()).unwrap();
        // access(WRITE_OK) uses the same credentials and ACLs as C.
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    checkpars();
    let mut vndorig = v2dop(1000000, nt, nxmax as i64);
    let nout = 1i64 + (nxmax as i64 - 1) * (ninterp as i64 + 1);
    // The C expression uses int arithmetic before conversion to long.
    if nout > i32::MAX as i64 || nout <= 0 || ninterp == i32::MAX {
        err!("invalid VND dimensions");
    }
    let mut vndinterp = v2dop(1000000, nt, nout);
    let mut headerfile = String::new();
    let headerfp = if tmpdir.is_empty() {
        let fp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
        fp
    } else {
        headerfile = temporary_filename(&tmpdir);
        let _ = HEADERFILE.set(headerfile.clone());
        // Darwin and POSIX use these signal numbers for HUP, INT, QUIT, TERM.
        unsafe {
            signal(2, closefiles);
            signal(3, closefiles);
            signal(1, closefiles);
            signal(15, closefiles);
        }
        let fp = efopen(&headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
        fp
    };
    let mut nx = 0;
    loop {
        V2Dw0(&mut vndorig, nx, tr.data_bytes_mut(), 1);
        efwrite(tr.hdr.as_bytes(), &headerfp);
        nx += 1;
        if nx >= nxmax as i64 || gettr(&mut tr) == 0 {
            break;
        }
    }
    jea_xinterpolate(
        &mut vndorig,
        &mut vndinterp,
        ninterp,
        nt,
        nx as usize,
        freq1,
        freq2,
        lagc,
        lent,
        lenx,
        xopt,
        dt,
        iopt,
        deriv,
        linear,
    );
    if iopt != 0 {
        ninterp = 0;
    }
    efseeko(&headerfp, 0, SEEK_SET);
    for ix in 0..nx - 1 {
        efread(tr.hdr.as_bytes_mut(), &headerfp);
        for k in 0..=ninterp {
            V2Dr0(
                &mut vndinterp,
                k as i64 + ix * (ninterp as i64 + 1),
                tr.data_bytes_mut(),
                18,
            );
            puttr(&mut tr);
        }
    }
    efread(tr.hdr.as_bytes_mut(), &headerfp);
    V2Dr0(
        &mut vndinterp,
        (nx - 1) * (ninterp as i64 + 1),
        tr.data_bytes_mut(),
        18,
    );
    puttr(&mut tr);
    VNDcl(vndorig, 1);
    VNDcl(vndinterp, 1);
    if VNDtotalmem() != 0 {
        fprintf!(stderr(), "total VND memory at end = %ld\n", VNDtotalmem());
    }
    efclose(&headerfp);
    if !tmpdir.is_empty() {
        eremove(&headerfile);
    }
    exit(CWP_Exit());
}
