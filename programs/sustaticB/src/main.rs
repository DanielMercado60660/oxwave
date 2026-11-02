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

//! SUSTATICB - Elevation static corrections.
//! Port of su/main/statics/sustaticB.c (Revision 1.1).
//! Credits: CWP: Jamie Burns; Mohammed Alfaraj (header statics and sign, 1992);
//! Timo Tjan (source/receiver statics files, 1995); Thomas Pratt, USGS
//! (suresstat files, 2000); J. W. Neese (hdrs=2,3 logic).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::sgn;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSTATICB - Elevation static corrections, apply corrections from\t",
    "\t      headers or from a source and receiver statics file\t",
    "\t      (beta submitted by J. W. Neese)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "     sustaticB <stdin >stdout  [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "\tv0=v1 or user-defined\tor from header, weathering velocity\t",
    "\tv1=user-defined\t\tor from header, subweathering velocity\t",
    "\thdrs=0\t\t\t=1 to read statics from headers\t\t",
    " \t\t\t\t=2 to read statics from files\t\t",
    "\t\t\t\t=3 to read from output files of suresstat",
    "\tsign=1\t\t\tapply static correction (add tstat values)",
    "\t\t\t\t=-1 apply negative of tstat values\t",
    " Options when hdrs=2 and hdrs=3:\t\t\t\t\t",
    "\tsou_file=\t\tinput file for source statics (ms) \t",
    "\trec_file=\t\tinput file for receiver statics (ms) \t",
    "\tns=240 \t\t(2)number of sources; (3) max fldr\t",
    "\tnr=335 \t\t\tnumber of receivers \t\t\t",
    "\tno=96 \t\t\tnumber of offsets\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " For hdrs=1, statics calculation is not performed, statics correction  ",
    " is applied to the data by reading statics (in ms) from the header.\t",
    "\t\t\t\t\t\t\t\t\t",
    " For hdrs=0, field statics are calculated, and\t\t\t\t",
    " \tinput field sut is assumed measured in ms.\t\t\t",
    " \toutput field sstat = 10^scalel*(sdel - selev + sdepth)/swevel\t",
    " \toutput field gstat = sstat - sut/1000.\t\t\t\t",
    " \toutput field tstat = sstat + gstat + 10^scalel*(selev - gelev)/wevel",
    "\t\t\t\t\t\t\t\t\t",
    " For hdrs=2, statics are surface consistently obtained from the \t",
    " statics files. The geometry should be regular.\t\t\t",
    " The source- and receiver-statics files should be unformated C binary \t",
    " floats and contain the statics (in ms) as a function of surface location.",
    "\t\t\t\t\t\t\t\t\t",
    " For hdrs=3, statics are read from the output files of suresstat, with ",
    " the same options as hdrs=2 (but use no=max traces per shot and assume ",
    " that ns=max fldr number and nr=max receiver number).\t\t\t",
    " For each shot number (trace header fldr) and each receiver number     ",
    " (trace header tracf) the program will look up the appropriate static  ",
    " correction.  The geometry need not be regular as each trace is treated",
    " independently.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat:  The static shifts are computed with the assumption that the  ",
    " desired datum is sea level (elevation=0). You may need to shift the\t",
    " selev and gelev header values via  suchw.\t\t\t\t",
    " Example: subtracting min(selev,gelev)=25094431\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suchw < CR290.su key1=selev,gelev key2=selev,gelev key3=selev,gelev \\ ",
    "            a=-25094431,-25094431 b=1,1 c=0,0 > CR290datum.su\t\t",
];

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut tsd = 0.0f32;
    let mut trd = 0.0f32;

    // Hook up getpar.
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get information from first trace.
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as usize;
    let tmin = (intrace.delrt as f64 / 1000.0) as f32;
    let dt = (intrace.dt as f64 / 1000000.0) as f32;

    // Get parameters. The original does not call checkpars().
    let mut v1 = 0.0f32;
    if getparfloat("v1", &mut v1) == 0 {
        v1 = intrace.swevel as f32;
    }
    let mut v0 = 0.0f32;
    if getparfloat("v0", &mut v0) == 0 {
        v0 = if intrace.wevel != 0 {
            intrace.wevel as f32
        } else {
            v1
        };
    }
    let mut hdrs = 0;
    if getparint("hdrs", &mut hdrs) == 0 {
        hdrs = 0;
    }
    let mut sign = 1;
    if getparint("sign", &mut sign) == 0 {
        sign = 1;
    }

    // Allocate vector of output times.
    let mut t = vec![0.0f32; nt];
    let mut ns = 240i32;
    let mut nr = 335i32;
    let mut no = 96i32;
    let mut sou_statics = Vec::new();
    let mut rec_statics = Vec::new();

    // Read source and receiver statics from files.
    if hdrs == 2 || hdrs == 3 {
        if getparint("ns", &mut ns) == 0 {
            ns = 240;
        }
        if getparint("nr", &mut nr) == 0 {
            nr = 335;
        }
        if getparint("no", &mut no) == 0 {
            no = 96;
        }
        let mut sou_file = String::new();
        let mut rec_file = String::new();
        getparstring("sou_file", &mut sou_file);
        getparstring("rec_file", &mut rec_file);
        // alloc1float returns NULL on failure (not ealloc1float's error).
        // In particular, negative counts must not panic before efopen can
        // report a missing file. Keep failed allocations empty until use.
        if rec_statics.try_reserve_exact(nr as usize).is_ok() {
            rec_statics.resize(nr as usize, 0.0f32);
        }
        if sou_statics.try_reserve_exact(ns as usize).is_ok() {
            sou_statics.resize(ns as usize, 0.0f32);
        }
        let fpr = efopen(&rec_file, "rb");
        // C would read through a NULL allocation here; reject that UB.
        if rec_statics.len() != nr as usize {
            err!("receiver statics allocation failed");
        }
        efread(&mut rec_statics, &fpr);
        efclose(&fpr);
        let fps = efopen(&sou_file, "rb");
        if sou_statics.len() != ns as usize {
            err!("source statics allocation failed");
        }
        efread(&mut sou_statics, &fps);
        efclose(&fps);
    }

    let mut tstat = 0.0f32;
    let mut io = 0i32;
    let mut is = 0i32;
    loop {
        let temp = (sgn(intrace.scalel) * (i32::from(intrace.scalel).abs() as f64).log10()) as i32;
        let scale = 10.0f64.powf(temp as f32 as f64) as f32;

        // Copy and adjust header.
        outtrace.hdr = intrace.hdr;
        if hdrs == 0 {
            tsd = scale
                * (intrace
                    .selev
                    .wrapping_neg()
                    .wrapping_add(intrace.sdel)
                    .wrapping_add(intrace.sdepth)) as f32
                / v1;
            trd = (tsd as f64 - intrace.sut as f64 / 1000.0) as f32;
            tstat = tsd + trd + scale * intrace.selev.wrapping_sub(intrace.gelev) as f32 / v0;
        } else {
            outtrace.sstat = intrace.sstat;
            outtrace.gstat = intrace.gstat;
            outtrace.tstat = (i32::from(intrace.sstat) + i32::from(intrace.gstat)) as i16;
            if hdrs == 1 {
                tstat = (outtrace.tstat as f64 / 1000.0) as f32;
            }
            if hdrs == 2 {
                let ir = is.wrapping_add(io);
                // Use the exclusive bounds of the minimal C repair, retaining
                // the previous static when a lookup is skipped. The nr=1
                // boundary case is byte-identical to the original (UB round 1).
                if is < ns {
                    if is < 0 || is as usize >= sou_statics.len() {
                        err!("source static index %d out of range for ns=%d", is, ns);
                    }
                    tsd = (sou_statics[is as usize] as f64 / 1000.0) as f32;
                }
                if ir > 0 && ir < nr {
                    if ir as usize >= rec_statics.len() {
                        err!("receiver static index %d out of range for nr=%d", ir, nr);
                    }
                    trd = (rec_statics[ir as usize] as f64 / 1000.0) as f32;
                }
                tstat = tsd + trd;
                io = io.wrapping_add(1);
                if io > no.wrapping_sub(1) {
                    io = 0;
                    is = is.wrapping_add(1);
                }
            }
            if hdrs == 3 {
                // Header indices are zero-based; C does not check their bounds.
                if intrace.fldr < 0 || intrace.fldr as usize >= sou_statics.len() {
                    err!(
                        "source static index %d out of range for ns=%d",
                        intrace.fldr,
                        ns
                    );
                }
                tsd = (sou_statics[intrace.fldr as usize] as f64 / 1000.0) as f32;
                if intrace.tracf < 0 || intrace.tracf as usize >= rec_statics.len() {
                    err!(
                        "receiver static index %d out of range for nr=%d",
                        intrace.tracf,
                        nr
                    );
                }
                trd = (rec_statics[intrace.tracf as usize] as f64 / 1000.0) as f32;
                tstat = tsd + trd;
            }
        }

        // Compute output times, then sinc interpolate new data.
        for itime in 0..nt {
            t[itime] = tmin + itime as f32 * dt + sign as f32 * tstat;
        }
        ints8r(
            dt,
            tmin,
            &intrace.data[..nt],
            0.0,
            0.0,
            &t,
            &mut outtrace.data[..nt],
        );

        if hdrs == 0 || hdrs == 2 || hdrs == 3 {
            // Add to existing headers to permit multiple static corrections.
            // C converts the double sum to int before narrowing to short.
            outtrace.sstat = (outtrace.sstat as f64 + 1000.0 * tsd as f64) as i32 as i16;
            outtrace.gstat = (outtrace.gstat as f64 + 1000.0 * trd as f64) as i32 as i16;
            outtrace.tstat = (outtrace.tstat as f64 + 1000.0 * tstat as f64) as i32 as i16;
        }
        puttr(&mut outtrace);
        if gettr(&mut intrace) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
