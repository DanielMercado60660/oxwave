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

//! SUKTMIG2D - prestack time migration with the DSR operator.
//! Port of su/main/migration_inversion/suktmig2d.c (Revision 1.7).
//! Credits (from C source): CWP: Carlos Pacheco, 2006 (Lanczos filter).

use oxwave_core::prelude::*;
use oxwave_cwp::complex::crmul;
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUKTMIG2D - prestack time migration of a common-offset\t",
    "\tsection with the double-square root (DSR) operator\t",
    "\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    "   suktmig2d < infile vfile= [parameters]  > outfile\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t",
    " vfile=\trms velocity file (units/s) v(t,x) as a function",
    "\t\tof time\t\t\t\t\t\t",
    " dx=\t\tdistance (units) between consecutive traces\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " fcdpdata=tr.cdp\tfirst cdp in data\t\t\t",
    " firstcdp=fcdpdata\tfirst cdp number in velocity file\t",
    " lastcdp=from header\tlast cdp number in velocity file\t",
    " dcdp=from header\tnumber of cdps between consecutive traces",
    " angmax=40\tmaximum aperture angle for migration (degrees)\t",
    " hoffset=.5*tr.offset\t\thalf offset (m)\t\t\t",
    " nfc=16\tnumber of Fourier-coefficients to approximate\t",
    "\t\tlow-pass\t\t\t\t\t",
    "\t\tfilters. The larger nfc the narrower the filter\t",
    " fwidth=5 \thigh-end frequency increment for the low-pass\t",
    " \t\tfilters\t\t\t\t\t\t",
    " \t\tin Hz. The lower this number the more the number",
    "\t\tof lowpass filters to be calculated for each \t",
    "\t\tinput trace.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Caveat: this code may need some work\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " Data must be preprocessed with sufrac to correct for the\t",
    " wave-shaping factor using phasefac=.25 for 2D migration.\t",
    "\t\t\t\t\t\t\t\t",
    " Input traces must be sorted into offset and cdp number.\t",
    " The velocity file consists of rms velocities for all CMPs as a",
    " function of vertical time and horizontal position v(t,x)\t",
    " in C-style binary floating point numbers.  It's easiest to \t",
    " supply v(t,x) that has the same dimensions as the input data to",
    " be migrated. Note that time t is the fast dimension in these  ",
    " the input velocity file.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " The units may be feet or meters, as long as these are\t\t",
    " consistent.\t\t\t\t\t\t\t",
    " Antialias filter is performed using (Gray,1992, Geoph. Prosp), ",
    " using nc low- pass filtered copies of the data. The cutoff\t",
    " frequencies are calculated  as fractions of the Nyquist\t",
    " frequency.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " The maximum allowed angle is 80 degrees(a 10 degree taper is ",
    " applied to the end of the aperture)\t\t\t\t",
];

const LOOKFAC: i32 = 2;
const PFA_MAX: i32 = 720720;

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut oldcdp = 0i32;
    let mut olddeltacdp = 0i32;
    let mut deltacdp: i32;
    let mut angtaper = 0.0f32;
    // C leaves the unused tail untouched between interpolation calls.
    let mut datalo = [0.0f32; 8];
    let mut datahi = [0.0f32; 8];
    let mut check_cdp = false;

    initargs(&argv);
    requestdoc(0, SDOC);
    /* Get info from first trace */
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as usize;
    let dt = intrace.dt as f32 / 1000000.0f32;
    let tmax = (nt as i32 - 1) as f32 * dt;
    let mut dx = 0.0f32;
    let mut vfile = String::new();
    let mut angmax = 40.0f32;
    let mut firstcdp = intrace.cdp;
    let mut fcdpdata = intrace.cdp;
    let mut hoffset = (0.5 * intrace.offset as f64) as f32;
    let mut nfc = 16i32;
    let mut fwidth = 5i32;
    let mut verbose = 0i32;
    mustgetparfloat("dx", &mut dx);
    mustgetparstring("vfile", &mut vfile);
    getparfloat("angmax", &mut angmax);
    getparint("firstcdp", &mut firstcdp);
    getparint("fcdpdata", &mut fcdpdata);
    getparfloat("hoffset", &mut hoffset);
    getparint("nfc", &mut nfc);
    getparint("fwidth", &mut fwidth);
    getparint("verbose", &mut verbose);
    let h = hoffset;

    /* Store traces in tmpfile while getting a count of number of traces */
    let tracefp = etmpfile();
    let hfp = etmpfile();
    let mut ntr = 0i32;
    loop {
        ntr += 1;
        deltacdp = intrace.cdp.wrapping_sub(oldcdp);
        efwrite(intrace.hdr.as_bytes(), &hfp);
        efwrite(&intrace.data[..nt], &tracefp);
        if ntr > 3 && olddeltacdp != deltacdp {
            if verbose != 0 {
                warn!("cdp interval changed in data");
                warn!(
                    "ntr=%d olddeltacdp=%d deltacdp=%d",
                    ntr, olddeltacdp, deltacdp
                );
                check_cdp = true;
            }
        }
        oldcdp = intrace.cdp;
        olddeltacdp = deltacdp;
        if gettr(&mut intrace) == 0 {
            break;
        }
    }
    let mut lastcdp = intrace.cdp;
    let mut dcdp = deltacdp.wrapping_sub(1);
    getparint("lastcdp", &mut lastcdp);
    getparint("dcdp", &mut dcdp);
    checkpars();
    if firstcdp == lastcdp || dcdp == 0 || check_cdp {
        warn!("Check cdp values in data!");
    }
    erewind(&tracefp);
    erewind(&hfp);
    // Invalid C dimensions otherwise become huge usize allocations in Rust.
    let ncdp = lastcdp
        .checked_sub(firstcdp)
        .and_then(|n| n.checked_add(1))
        .filter(|&n| n > 0)
        .unwrap_or_else(|| {
            err!(
                "invalid velocity CDP bounds: firstcdp=%d lastcdp=%d",
                firstcdp,
                lastcdp
            )
        });
    /* Set up FFT parameters */
    let nfft = npfaro(nt as i32, LOOKFAC * nt as i32);
    if nfft >= SU_NFLTS as i32 || nfft >= PFA_MAX {
        err!("Padded nt=%d -- too big", nfft);
    }
    let nf = (nfft / 2 + 1) as usize;
    if dt == 0.0 {
        err!("trace sample interval dt must be positive");
    }
    if fwidth == 0 {
        err!("fwidth must be nonzero");
    }
    let fnyq = (1.0 / (2.0f32 * dt) as f64) as f32;
    let nc = (fnyq / fwidth as f32).ceil() as i32;
    // A negative fwidth can legitimately give nc=0; retain that C behavior.
    let nfilters = nc
        .checked_add(1)
        .filter(|&n| n > 0)
        .unwrap_or_else(|| err!("invalid filter count nc=%d for fwidth=%d", nc, fwidth));
    if verbose != 0 {
        warn!(" The number of filters for antialiasing is nc= %d", nc);
    }
    /* Allocate space */
    // C uses signed int products in its reads and clears. Reject overflow
    // before converting dimensions to usize or requesting any array storage.
    let sample_count = |rows: i32| {
        (nt as i32)
            .checked_mul(rows)
            .filter(|&n| n >= 0)
            .map(|n| n as usize)
            .filter(|&n| {
                n.checked_mul(std::mem::size_of::<f32>())
                    .is_some_and(|bytes| bytes <= isize::MAX as usize)
            })
            .unwrap_or_else(|| err!("array dimensions too large: nt=%d rows=%d", nt, rows))
    };
    let ndata = sample_count(ntr);
    let nlowpass = sample_count(nfilters);
    let nvel = sample_count(ncdp);
    let alloc_floats = |len: usize| {
        let mut array = Vec::new();
        if array.try_reserve_exact(len).is_err() {
            err!("cannot allocate %lu float samples", len);
        }
        array.resize(len, 0.0f32);
        array
    };
    let mut data = Vec::new();
    if data.try_reserve_exact(ntr as usize).is_err() {
        err!("cannot allocate %d trace rows", ntr);
    }
    for _ in 0..ntr {
        data.push(alloc_floats(nt));
    }
    // C clears these arrays through row zero as contiguous blocks.
    let mut lowpass = alloc_floats(nlowpass);
    let mut mig = alloc_floats(ndata);
    let mut vel = alloc_floats(nvel);
    let mut fc = alloc_floats(nfilters as usize);
    let mut rtin = alloc_floats(nfft as usize);
    let mut rtout = alloc_floats(nfft as usize);
    let mut ct = vec![Complex::default(); nf];
    let mut filter = alloc_floats(nf);
    for ix in 0..ntr as usize {
        efread(&mut data[ix], &tracefp);
    }
    let vfp = efopen(&vfile, "r");
    efread(&mut vel, &vfp);
    efclose(&vfp);
    /* Zero all arrays */
    mig.fill(0.0);
    rtin.fill(0.0);
    filter.fill(0.0);
    lowpass.fill(0.0);
    for i in 1..nc + 1 {
        fc[i as usize] = fnyq * i as f32 / nc as f32;
    }

    /* Start migration: loop over input mid-points first */
    if verbose != 0 {
        warn!("Starting migration process...\n");
    }
    for imp in 0..ntr {
        // Preserve the C absolute row expression (it does not subtract
        // firstcdp), but reject its undefined out-of-range accesses.
        let velocity_row = imp as i64 * dcdp as i64 + fcdpdata as i64 - 1;
        if nt != 0 && (velocity_row < 0 || velocity_row >= ncdp as i64) {
            err!(
                "velocity row %ld out of range [0,%d) for input trace %d",
                velocity_row,
                ncdp,
                imp + 1
            );
        }
        let velocity_start = velocity_row as usize * nt;
        let mp = imp as f32 * dx;
        let perc = (imp as f64 * 100.0 / (ntr - 1) as f64) as f32;
        if (imp.wrapping_mul(100) as f64 % (ntr - 1) as f64) == 0.0 && verbose != 0 {
            warn!("migrated %g\n ", perc);
        }
        /* Calculate low-pass filtered versions for antialiasing */
        for it in 0..nt {
            rtin[it] = data[imp as usize][it];
        }
        for ifc in 1..nc + 1 {
            rtout.fill(0.0);
            // C clears nf floats; pfarc subsequently overwrites all ct.
            for it in 0..nf / 2 {
                ct[it] = Complex::default();
            }
            if nf % 2 != 0 {
                ct[nf / 2].r = 0.0;
            }
            lpfilt(nfc, nfft, dt, fc[ifc as usize], &mut filter);
            pfarc(1, &rtin, &mut ct);
            for it in 0..nf {
                ct[it] = crmul(ct[it], filter[it]);
            }
            pfacr(-1, &ct, &mut rtout);
            for it in 0..nt {
                lowpass[(ifc as usize) * nt + (it)] = rtout[it];
            }
        }
        /* Loop over vertical traveltimes */
        for it in 0..nt {
            let mut t0 = it as f32 * dt;
            let v = vel[velocity_start + it];
            let xmax = (((angmax as f64 + 10.0) * PI / 180.0).tan() * v as f64 * t0 as f64) as f32;
            let lxf = imp as f64 - ((xmax / dx) as f64).ceil();
            let lx = (if 0.0 > lxf { 0.0 } else { lxf }) as i32;
            let uxf = imp as f64 + ((xmax / dx) as f64).ceil();
            let ux = (if (ntr as f64) < uxf { ntr as f64 } else { uxf }) as i32;
            /* Output image-points to the left of the midpoint */
            let mut iip = imp;
            while iip > lx {
                let ip = iip as f32 * dx;
                let x = ip - mp;
                let ts = (((t0 / 2.0f32) as f64).powf(2.0) + (((x + h) / v) as f64).powf(2.0))
                    .sqrt() as f32;
                let tr = (((t0 / 2.0f32) as f64).powf(2.0) + (((h - x) / v) as f64).powf(2.0))
                    .sqrt() as f32;
                let t = ts + tr;
                if t >= tmax {
                    break;
                }
                let geoms = ((1.0f32 / (t * v)) as f64).sqrt() as f32;
                let obliq = (0.5
                    * ((1.0f32 + t0 * t0 / (4.0f32 * ts * tr)) as f64
                        - (1.0f32 / (ts * tr)) as f64
                            * ((ts * ts - t0 * t0 / 4.0f32) as f64).sqrt()
                            * ((tr * tr - t0 * t0 / 4.0f32) as f64).sqrt()))
                .sqrt() as f32;
                let ang = (180.0 * ((t0 / t) as f64).acos().abs() / PI) as f32;
                if ang <= angmax {
                    angtaper = 1.0;
                }
                if ang > angmax {
                    angtaper = (((ang - angmax) as f64) * PI / 20.0).cos() as f32;
                }
                /* Evaluate migration operator slowness for antialiasing */
                let pmin = 1.0f32 / (2.0f32 * dx * fnyq);
                let p = ((x + h) as f64 / ((v as f64).powf(2.0) * ts as f64)
                    + (x - h) as f64 / ((v as f64).powf(2.0) * tr as f64))
                    .abs() as f32;
                let mut fplo = 0i32;
                if p > 0.0 {
                    fplo = (nc as f32 * pmin / p).floor() as i32;
                }
                if p == 0.0 {
                    fplo = nc;
                }
                let ref_ = ((nc as f32 * pmin) as f64 % p as f64) as f32;
                let wlo = 1.0f32 - ref_;
                fplo = fplo.wrapping_add(1);
                let fphi = fplo;
                let whi = ref_;
                let itbf = ((t / dt) as f64).ceil() - 3.0;
                let itb = (if itbf > 0.0 { itbf } else { 0.0 }) as i32;
                let ite = (itb + 8).min(nt as i32);
                let firstt = (itb - 1) as f32 * dt;
                /* Move energy from CMP to CIP */
                if fplo >= nc {
                    for k in itb..ite {
                        datalo[(k - itb) as usize] = lowpass[(nc as usize) * nt + (k as usize)];
                    }
                    let mut amplo = [0.0f32];
                    ints8r(dt, firstt, &datalo, 0.0, 0.0, &[t], &mut amplo);
                    mig[(iip as usize) * nt + (it)] += geoms * obliq * angtaper * amplo[0];
                } else if fplo < nc {
                    for k in itb..ite {
                        datalo[(k - itb) as usize] = lowpass[(fplo as usize) * nt + (k as usize)];
                        datahi[(k - itb) as usize] = lowpass[(fphi as usize) * nt + (k as usize)];
                    }
                    let mut amplo = [0.0f32];
                    let mut amphi = [0.0f32];
                    ints8r(dt, firstt, &datalo, 0.0, 0.0, &[t], &mut amplo);
                    ints8r(dt, firstt, &datahi, 0.0, 0.0, &[t], &mut amphi);
                    mig[(iip as usize) * nt + (it)] +=
                        geoms * obliq * angtaper * (wlo * amplo[0] + whi * amphi[0]);
                }
                iip -= 1;
            }
            /* Output image-points to the right of the midpoint */
            let mut iip = imp + 1;
            while iip < ux {
                let ip = iip as f32 * dx;
                let x = ip - mp;
                t0 = it as f32 * dt;
                let ts = (((t0 / 2.0f32) as f64).powf(2.0) + (((x + h) / v) as f64).powf(2.0))
                    .sqrt() as f32;
                let tr = (((t0 / 2.0f32) as f64).powf(2.0) + (((h - x) / v) as f64).powf(2.0))
                    .sqrt() as f32;
                let t = ts + tr;
                if t >= tmax {
                    break;
                }
                let geoms = ((1.0f32 / (t * v)) as f64).sqrt() as f32;
                let obliq = (0.5
                    * ((1.0f32 + t0 * t0 / (4.0f32 * ts * tr)) as f64
                        - (1.0f32 / (ts * tr)) as f64
                            * ((ts * ts - t0 * t0 / 4.0f32) as f64).sqrt()
                            * ((tr * tr - t0 * t0 / 4.0f32) as f64).sqrt()))
                .sqrt() as f32;
                let ang = (180.0 * ((t0 / t) as f64).acos().abs() / PI) as f32;
                if ang <= angmax {
                    angtaper = 1.0;
                }
                if ang > angmax {
                    angtaper = (((ang - angmax) as f64) * PI / 20.0).cos() as f32;
                }
                let pmin = 1.0f32 / (2.0f32 * dx * fnyq);
                let p = ((x + h) as f64 / ((v as f64).powf(2.0) * ts as f64)
                    + (x - h) as f64 / ((v as f64).powf(2.0) * tr as f64))
                    .abs() as f32;
                let mut fplo = 0i32;
                if p > 0.0 {
                    fplo = (nc as f32 * pmin / p).floor() as i32;
                }
                if p == 0.0 {
                    fplo = nc;
                }
                let ref_ = ((nc as f32 * pmin) as f64 % p as f64) as f32;
                let wlo = 1.0f32 - ref_;
                let fphi = fplo.wrapping_add(1);
                let whi = ref_;
                let itbf = ((t / dt) as f64).ceil() - 3.0;
                let itb = (if itbf > 0.0 { itbf } else { 0.0 }) as i32;
                let ite = (itb + 8).min(nt as i32);
                let firstt = (itb - 1) as f32 * dt;
                if fplo >= nc {
                    for k in itb..ite {
                        datalo[(k - itb) as usize] = lowpass[(nc as usize) * nt + (k as usize)];
                    }
                    let mut amplo = [0.0f32];
                    ints8r(dt, firstt, &datalo, 0.0, 0.0, &[t], &mut amplo);
                    mig[(iip as usize) * nt + (it)] += geoms * obliq * angtaper * amplo[0];
                } else if fplo < nc {
                    for k in itb..ite {
                        datalo[(k - itb) as usize] = lowpass[(fplo as usize) * nt + (k as usize)];
                        datahi[(k - itb) as usize] = lowpass[(fphi as usize) * nt + (k as usize)];
                    }
                    let mut amplo = [0.0f32];
                    let mut amphi = [0.0f32];
                    ints8r(dt, firstt, &datalo, 0.0, 0.0, &[t], &mut amplo);
                    ints8r(dt, firstt, &datahi, 0.0, 0.0, &[t], &mut amphi);
                    mig[(iip as usize) * nt + (it)] +=
                        geoms * obliq * angtaper * (wlo * amplo[0] + whi * amphi[0]);
                }
                iip += 1;
            }
        }
    }
    /* Output migrated data */
    erewind(&hfp);
    for ix in 0..ntr as usize {
        efread(outtrace.hdr.as_bytes_mut(), &hfp);
        for it in 0..nt {
            outtrace.data[it] = mig[(ix) * nt + (it)];
        }
        puttr(&mut outtrace);
    }
    efclose(&hfp);
    exit(CWP_Exit());
}

/* Low-pass filter using Lanczos smoothing (R.W. Hamming, 1977).
 * Author: CWP: Carlos Pacheco, 2006. */
fn lpfilt(nfc: i32, nfft: i32, dt: f32, fhi: f32, filter: &mut [f32]) {
    let nf = nfft / 2 + 1;
    let onfft = (1.0 / nfft as f64) as f32;
    let fn_ = (1.0 / (2.0f32 * dt) as f64) as f32;
    let df = onfft / dt;
    let whi = (fhi as f64 * PI / fn_ as f64) as f32;
    let dw = (df as f64 * PI / fn_ as f64) as f32;
    for i in 0..nf {
        filter[i as usize] = (whi as f64 / PI) as f32;
        let w = i as f32 * dw;
        for j in 1..nfc {
            let c = (((whi * j as f32) as f64).sin()
                * (PI * j as f64 / nfc as f64).sin()
                * 2.0
                * nfc as f64
                / (PI * PI * j as f64 * j as f64)) as f32;
            filter[i as usize] =
                (filter[i as usize] as f64 + c as f64 * ((j as f32 * w) as f64).cos()) as f32;
        }
    }
}
