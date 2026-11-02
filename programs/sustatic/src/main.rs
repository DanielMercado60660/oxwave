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

//! SUSTATIC - elevation static corrections. Port of su/main/statics/sustatic.c (Revision 1.26).
//! Credits (from the C source):
//! CWP: Jamie Burns.
//! CWP: Modified by Mohammed Alfaraj, 11/10/1992, for reading statics
//! from headers and including sign (+-) option.
//! CWP: Modified by Timo Tjan, 29 June 1995, to include input of source
//! and receiver statics from files.
//! Modified by Thomas Pratt, USGS, Feb, 2000 to read statics from the
//! output files of suresstat.

use oxwave_core::prelude::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSTATIC - Elevation static corrections, apply corrections from\t",
    "\t      headers or from a source and receiver statics file\t",
    "\t\t\t\t\t\t\t\t\t",
    "     sustatic <stdin >stdout  [optional parameters]\t \t\t",
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
    "\tns=240 \t\t\tnumber of souces \t\t\t",
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
    " that ns=max shot number and nr=max receiver number).\t\t\t",
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
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as usize;
    let tmin = (intrace.delrt as f64 / 1000.0) as f32;
    let dt = (intrace.dt as f64 / 1000000.0) as f32;
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
    let mut sign = 1;
    getparint("hdrs", &mut hdrs);
    getparint("sign", &mut sign);
    let mut t = vec![0.0f32; nt];
    let mut ns = 240;
    let mut nr = 335;
    let mut no = 96;
    let mut sou_statics = Vec::new();
    let mut rec_statics = Vec::new();
    if hdrs == 2 || hdrs == 3 {
        getparint("ns", &mut ns);
        getparint("nr", &mut nr);
        getparint("no", &mut no);
        let mut sou_file = String::new();
        let mut rec_file = String::new();
        getparstring("sou_file", &mut sou_file);
        getparstring("rec_file", &mut rec_file);
        checkpars();
        // Never convert negative counts into allocation sizes. Defer their
        // errors until the corresponding read to preserve C's file-open order.
        if nr >= 0 {
            rec_statics = vec![0.0f32; nr as usize];
        }
        if ns >= 0 {
            sou_statics = vec![0.0f32; ns as usize];
        }
        let fpr = efopen(&rec_file, "rb");
        if nr < 0 {
            // C's failed allocation makes fread fail with this diagnostic.
            err!("subcalls.c: efread: fread only %d items of %d", 0, nr);
        }
        efread(&mut rec_statics, &fpr);
        efclose(&fpr);
        let fps = efopen(&sou_file, "rb");
        if ns < 0 {
            err!("subcalls.c: efread: fread only %d items of %d", 0, ns);
        }
        efread(&mut sou_statics, &fps);
        efclose(&fps);
    }
    let mut io = 0i32;
    let mut is = 0i32;
    let mut tsd = 0.0f32;
    let mut trd = 0.0f32;
    loop {
        // C converts log10(0) = -infinity to INT_MIN, giving scale = 0.
        let temp = (oxwave_cwp::cwp::sgn(intrace.scalel as f64)
            * ((intrace.scalel as i32).abs() as f64).log10()) as i32;
        let scale = 10.0f64.powf(temp as f32 as f64) as f32;
        outtrace.hdr = intrace.hdr;
        let tstat: f32;
        if hdrs == 0 {
            tsd = scale
                * intrace
                    .selev
                    .wrapping_neg()
                    .wrapping_add(intrace.sdel)
                    .wrapping_add(intrace.sdepth) as f32
                / v1;
            trd = (tsd as f64 - intrace.sut as f64 / 1000.0) as f32;
            tstat = tsd + trd + scale * intrace.selev.wrapping_sub(intrace.gelev) as f32 / v0;
        } else {
            if hdrs == 2 {
                let ir = is.wrapping_add(io);
                // Keep C's conditional lookups (including skipped receivers),
                // but reject its undefined accesses at the inclusive limits.
                if is <= ns {
                    tsd = *sou_statics
                        .get(is as usize)
                        .unwrap_or_else(|| err!("source index %d out of range for ns=%d", is, ns))
                        / 1000.0;
                }
                if ir > 0 && ir <= nr {
                    trd = *rec_statics.get(ir as usize).unwrap_or_else(|| {
                        err!("receiver index %d out of range for nr=%d", ir, nr)
                    }) / 1000.0;
                }
                intrace.tstat = (tsd + trd) as i32 as i16;
                io = io.wrapping_add(1);
                if io > no.wrapping_sub(1) {
                    io = 0;
                    is = is.wrapping_add(1);
                }
            }
            if hdrs == 3 {
                tsd = *sou_statics.get(intrace.fldr as usize).unwrap_or_else(|| {
                    err!(
                        "source index fldr=%d out of range for ns=%d",
                        intrace.fldr,
                        ns
                    )
                });
                trd = *rec_statics.get(intrace.tracf as usize).unwrap_or_else(|| {
                    err!(
                        "receiver index tracf=%d out of range for nr=%d",
                        intrace.tracf,
                        nr
                    )
                });
                intrace.tstat = (tsd + trd) as i32 as i16;
                io = io.wrapping_add(1);
                if io > no.wrapping_sub(1) {
                    io = 0;
                    is = is.wrapping_add(1);
                }
            }
            if intrace.tstat == 0 {
                outtrace.tstat = intrace.sstat.wrapping_add(intrace.gstat);
                tstat = (outtrace.tstat as f64 / 1000.0) as f32;
            } else {
                tstat = (intrace.tstat as f64 / 1000.0) as f32;
            }
        }
        for itime in 0..nt {
            t[itime] = tmin + itime as f32 * dt + sign as f32 * tstat;
        }

        // inttable8.c:177 converts the normalized coordinate to int. Beyond
        // that range C has UB, but a minimal bounds fix still produces zero
        // extrapolation, byte-identical to SU for the huge-shift regression.
        // Interpolate valid spans normally and supply zero for those points.
        let xouts = (1.0f64 / dt as f64) as f32;
        let xoutb = (8.0f64 - (tmin * xouts) as f64) as f32;
        let mut start = 0;
        for itime in 0..nt {
            let xoutn = xoutb + t[itime] * xouts;
            if xoutn as f64 > i32::MAX as f64 || (xoutn as f64) < i32::MIN as f64 + 11.0 {
                ints8r(
                    dt,
                    tmin,
                    &intrace.data[..nt],
                    0.0,
                    0.0,
                    &t[start..itime],
                    &mut outtrace.data[start..itime],
                );
                outtrace.data[itime] = 0.0;
                start = itime + 1;
            }
        }
        ints8r(
            dt,
            tmin,
            &intrace.data[..nt],
            0.0,
            0.0,
            &t[start..],
            &mut outtrace.data[start..nt],
        );
        if hdrs == 0 || hdrs == 2 || hdrs == 3 {
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
