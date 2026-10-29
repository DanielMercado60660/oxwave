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

//! SUSTATICRRS - elevation and residual refraction static corrections.
//! Port of su/main/statics/sustaticrrs.c (Revision 1.23, 2011/11/16).
//! Additional original notice: Copyright (c) Colorado School of Mines, 2011.
//! All rights reserved.
//!
//! Credits (from the C source): CWP: Jamie Burns; Mohammed Alfaraj,
//! 11/10/1992 (headers and sign); Timo Tjan, 29 June 1995 (statics files);
//! Chris Robinson, 11/2000 (residual and average refraction statics).
//! Reference: Cox, M., 1999, Static corrections for seismic reflection
//! surveys: Society of Exploration Geophysicists.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{nint, sgn};
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUSTATICRRS - Elevation STATIC corrections, apply corrections from\t",
    "\t      headers or from a source and receiver statics file,\t",
    "\t      includes application of Residual Refraction Statics\t",
    "\t\t\t\t\t\t\t\t\t",
    "     sustaticrrs <stdin >stdout  [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "\tv0=v1 or user-defined\tor from header, weathering velocity\t",
    "\tv1=user-defined\t\tor from header, subweathering velocity\t",
    "\thdrs=0\t\t\t=1 to read statics from headers\t\t",
    " \t\t\t\t=2 to read statics from files\t\t",
    "\tsign=1\t\t\t=-1 to subtract statics from traces(up shift)",
    " Options when hdrs=2:\t\t\t\t\t\t\t",
    "\tsou_file=\t\tinput file for source statics (ms) \t",
    "\trec_file=\t\tinput file for receiver statics (ms) \t",
    "\tns=240 \t\t\tnumber of sources \t\t\t",
    "\tnr=335 \t\t\tnumber of receivers \t\t\t",
    "\tno=96 \t\t\tnumber of offsets\t\t\t",
    "                                                                       ",
    " Options when hdrs=3:                                                  ",
    "       blvl_file=              base of the near-surface model file (sampled",
    "                                  at CMP locations)                    ",
    "       refr_file=              horizontal reference datum file (sampled at",
    "                                  CMP locations)                       ",
    "       nsamp=                  number of midpoints on line             ",
    "       fx=                     first x location in velocity model      ",
    "       dx=                     midpoint interval                       ",
    "       V_r=                    replacement velocity                    ",
    "       mx=                     number of velocity model samples in     ",
    "                                  lateral direction                    ",
    "       mz=                     number of velocity model samples in     ",
    "                                  vertical direction                   ",
    "       dzv=                    velocity model depth interval           ",
    "       vfile=                  near-surface velocity model             ",
    "                                                                       ",
    " Options when hdrs=4:                                                  ",
    "       nsamp=                  number of midpoints on line             ",
    "       fx=                     first x location in velocity model      ",
    "       dx=                     midpoint interval                       ",
    "                                                                       ",
    " Options when hdrs=5:                                                  ",
    "       none                                                            ",
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
    " For hdrs=3, residual refraction statics and average refraction statics",
    " are computed.  For hdrs=4, residual refraction statics are applied,   ",
    " and for hdrs=5, average refraction statics are applied (Cox, 1999).   ",
    " These three options are coupled in many data processing sequences:    ",
    " before stack residual and average refraction statics are computed but ",
    " only residual refractions statics are applied, and after stack average",
    " refraction statics are applied.  Refraction statics are often split   ",
    " like this to avoid biasing stacking velocities.  The files blvl_file  ",
    " and refr_file are the base of the velocity model defined in vfile and ",
    " the final reference datum, as described by Cox (1999), respectively.  ",
    " Residual refraction statics are stored in the header field gstat, and ",
    " the average statics are stored in the header field tstat.  V_r is the ",
    " replacement velocity as described by Cox (1999).  The velocity file,  ",
    " vfile, is designed to work with a horizontal upper surface defined in ",
    " refr_file.  If the survey has irregular topography, the horizontal    ",
    " upper surface should be above the highest topographic point on the    ",
    " line, and the velocity between this horizontal surface and topography ",
    " should be some very large value, such as 999999999, so that the       ",
    " traveltimes through that region are inconsequential.                  ",
];

// C leaves out-of-range grid accesses undefined. Diagnose those inputs before
// indexing, retaining C's f32 arithmetic and truncation toward zero (including
// fractional positions between -1 and 0, which select the first grid sample).
fn grid_index(coordinate: i32, fx: f32, dx: f32, size: i32, key: &str) -> usize {
    let position = (coordinate as f32 - fx) / dx;
    let index = position as i32;
    if !position.is_finite() || index < 0 || index >= size {
        err!(
            "%s grid index out of range: %g (size %d)",
            key,
            position,
            size
        );
    }
    index as usize
}

fn main() {
    let argv = args();
    let mut intrace = Trace::new();
    let mut outtrace = Trace::new();
    let mut tsd = 0.0f32;
    let mut trd = 0.0f32;
    let mut tstat = 0.0f32;
    let mut res = 0i32;
    let mut bulk = 0i32;

    /* Hook up getpar */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get information from first trace */
    if gettr(&mut intrace) == 0 {
        err!("can't get first trace");
    }
    let nt = intrace.ns as usize;
    let tmin = (intrace.delrt as f64 / 1000.0) as f32;
    let dt = (intrace.dt as f64 / 1000000.0) as f32;

    /* Get parameters */
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
    let mut hdrs = 0i32;
    let mut sign = 1i32;
    getparint("hdrs", &mut hdrs);
    getparint("sign", &mut sign);

    /* Allocate vector of output times */
    let mut t = vec![0.0f32; nt];
    let mut ns = 240i32;
    let mut nr = 335i32;
    let mut no = 96i32;
    // C has no defaults for these variables. Valid hdrs=3/4 input supplies them.
    let mut nsamp = 0i32;
    let mut fx = 0.0f32;
    let mut dx = 0.0f32;
    let mut v_r = 0.0f32;
    let mut mx = 0i32;
    let mut mz = 0i32;
    let mut dzv = 0.0f32;
    let mut sou_statics = Vec::<f32>::new();
    let mut rec_statics = Vec::<f32>::new();
    let mut blvl_arr = Vec::<f32>::new();
    let mut refr_arr = Vec::<f32>::new();
    let mut tcdp_stat = Vec::<f32>::new();
    let mut tcdp_arr = Vec::<f32>::new();
    let mut num = Vec::<f32>::new();
    let mut num_arr = Vec::<f32>::new();
    let mut vel = Vec::<f32>::new();
    let mut fpb = None;
    let mut fpf = None;
    let mut g1 = None;
    let mut g2 = None;
    let mut g3 = None;
    let mut g4 = None;

    /* Reading source and receiver statics from files */
    if hdrs == 2 {
        getparint("ns", &mut ns);
        getparint("nr", &mut nr);
        getparint("no", &mut no);
        let mut sou_file = String::new();
        let mut rec_file = String::new();
        getparstring("sou_file", &mut sou_file);
        getparstring("rec_file", &mut rec_file);
        rec_statics = vec![0.0; nr as usize];
        sou_statics = vec![0.0; ns as usize];
        let fpr = efopen(&rec_file, "rb");
        efread(&mut rec_statics, &fpr);
        efclose(&fpr);
        let fps = efopen(&sou_file, "rb");
        efread(&mut sou_statics, &fps);
        efclose(&fps);
    }

    if hdrs == 3 {
        let mut blvl_file = String::new();
        let mut refr_file = String::new();
        getparstring("blvl_file", &mut blvl_file);
        getparstring("refr_file", &mut refr_file);
        getparint("nsamp", &mut nsamp);
        getparfloat("fx", &mut fx);
        getparfloat("dx", &mut dx);
        getparfloat("V_r", &mut v_r);
        if nsamp <= 0 {
            err!("nsamp must be positive");
        }
        blvl_arr = vec![0.0; nsamp as usize];
        fpb = Some(efopen(&blvl_file, "rb"));
        efread(&mut blvl_arr, fpb.as_ref().unwrap());
        refr_arr = vec![0.0; nsamp as usize];
        fpf = Some(efopen(&refr_file, "rb"));
        efread(&mut refr_arr, fpf.as_ref().unwrap());
        g1 = fopen("tcdp.b", "wb");
        g2 = fopen("num.b", "wb");
        let mut vfile = String::new();
        if getparstring("vfile", &mut vfile) == 0 {
            err!("must specify vfile");
        }
        if getparint("mx", &mut mx) == 0 {
            err!("must specify mx");
        }
        if getparint("mz", &mut mz) == 0 {
            err!("must specify mz");
        }
        if getparfloat("dzv", &mut dzv) == 0 {
            err!("must specify dzv");
        }
        if mx <= 0 || mz < 0 {
            err!("mx must be positive and mz must be nonnegative");
        }
        let vfp = fopen(&vfile, "r");
        // alloc2float(mz,mx) is contiguous, with depth the inner dimension.
        vel = vec![0.0; (mx as usize) * (mz as usize)];
        if let Some(vfp) = vfp {
            fread(&mut vel, &vfp);
            fclose(&vfp);
        } else {
            // C calls fread on NULL here. Preserve its shell exit status
            // without a C signal call; malformed input has no defined output.
            exit(139);
        }
        tcdp_stat = vec![0.0; nsamp as usize];
        num = vec![0.0; nsamp as usize];
        for i in 0..nsamp as usize {
            tcdp_stat[i] = 0.0;
            num[i] = 0.0;
        }
    }

    if hdrs == 4 {
        getparint("nsamp", &mut nsamp);
        getparfloat("fx", &mut fx);
        getparfloat("dx", &mut dx);
        if nsamp <= 0 {
            err!("nsamp must be positive");
        }
        tcdp_arr = vec![0.0; nsamp as usize];
        g3 = Some(efopen("tcdp.b", "rb"));
        efread(&mut tcdp_arr, g3.as_ref().unwrap());
        num_arr = vec![0.0; nsamp as usize];
        g4 = Some(efopen("num.b", "rb"));
        efread(&mut num_arr, g4.as_ref().unwrap());
    }
    checkpars();
    if hdrs == 3 || hdrs == 4 {
        if !fx.is_finite() {
            err!("fx must be finite");
        }
        // Negative spacing is valid; zero or NaN cannot produce a C int index.
        // getpar already diagnoses infinite parameters as eatof overflow.
        if dx == 0.0 || !dx.is_finite() {
            err!("dx must be finite and nonzero");
        }
    }

    /* Loop on traces */
    let mut io = 0i32;
    let mut is = 0i32;
    loop {
        let temp = (sgn(intrace.scalel) * ((intrace.scalel as i32).abs() as f64).log10()) as i32;
        let scale = 10.0f64.powf(temp as f32 as f64) as f32;
        /* Copy and adjust header */
        outtrace.hdr = intrace.hdr;

        /* Compute static correction if necessary */
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
            /* Else, read statics from headers or files */
            if hdrs == 2 {
                let ir = is.wrapping_add(io);
                if is < ns {
                    tsd = (sou_statics[is as usize] as f64 / 1000.0) as f32;
                }
                if ir >= 0 && ir < nr {
                    trd = (rec_statics[ir as usize] as f64 / 1000.0) as f32;
                }
                intrace.tstat = nint(tsd as f64 * 1000.0 + trd as f64 * 1000.0) as i16;
                io = io.wrapping_add(1);
                if io > no.wrapping_sub(1) {
                    io = 0;
                    is = is.wrapping_add(1);
                }
            }

            if hdrs == 3 {
                let is_c = grid_index(intrace.sx, fx, dx, nsamp, "sx");
                let z_blvl_s = blvl_arr[is_c];
                let z_refr_s = refr_arr[is_c];
                let ir_c = grid_index(intrace.gx, fx, dx, nsamp, "gx");
                let z_blvl_r = blvl_arr[ir_c];
                let z_refr_r = refr_arr[ir_c];
                tsd = 0.0;
                trd = 0.0;
                // The velocity grid may have fewer columns than nsamp. It is
                // accessed only when the depth integration loop is nonempty.
                let columns = if mz > 1 { mx } else { nsamp };
                let is_c = grid_index(intrace.sx, fx, dx, columns, "sx");
                let ir_c = grid_index(intrace.gx, fx, dx, columns, "gx");
                for j in 0..mz.wrapping_sub(1) {
                    tsd = tsd + dzv / vel[is_c * mz as usize + j as usize];
                    trd = trd + dzv / vel[ir_c * mz as usize + j as usize];
                }
                tsd = tsd - (z_blvl_s - z_refr_s) / v_r;
                trd = trd - (z_blvl_r - z_refr_r) / v_r;
                intrace.tstat = nint(tsd as f64 * 1000.0 + trd as f64 * 1000.0) as i16;
                let ic_c = grid_index(intrace.cdp, fx, dx, nsamp, "cdp");
                tcdp_stat[ic_c] = tcdp_stat[ic_c] + tsd + trd;
                num[ic_c] = (num[ic_c] as f64 + 1.0) as f32;
            }

            if hdrs == 4 {
                let ic_c = grid_index(intrace.cdp, fx, dx, nsamp, "cdp");
                if num_arr[ic_c] as f64 >= 0.5 {
                    bulk = nint(1000.0 * tcdp_arr[ic_c] as f64 / num_arr[ic_c] as f64);
                    res = (intrace.tstat as i32).wrapping_sub(bulk);
                }
                if num_arr[ic_c] as f64 <= 0.5 {
                    bulk = 0;
                    res = 0;
                }
            }

            /* If total statics not supplied, calculate it */
            if intrace.tstat == 0 && hdrs != 3 && hdrs != 4 {
                outtrace.tstat = (intrace.sstat as i32 + intrace.gstat as i32) as i16;
                tstat = (outtrace.tstat as f64 / 1000.0) as f32;
            }
            /* For hdrs=5, here is where tstat is calculated */
            if intrace.tstat != 0 && hdrs != 3 && hdrs != 4 {
                tstat = (intrace.tstat as f64 / 1000.0) as f32;
            }
            if hdrs == 3 {
                tstat = 0.0;
            }
            if hdrs == 4 {
                tstat = (res as f64 / 1000.0) as f32;
            }
        }

        /* Compute output times */
        for itime in 0..nt {
            t[itime] = tmin + itime as f32 * dt + sign as f32 * tstat;
        }

        /* Sinc interpolate new data */
        ints8r(
            dt,
            tmin,
            &intrace.data[..nt],
            0.0,
            0.0,
            &t,
            &mut outtrace.data[..nt],
        );

        /* Set header field for output trace */
        if hdrs != 1 && hdrs != 3 && hdrs != 4 && hdrs != 5 {
            /* Add to existing header values to permit multiple corrections. */
            outtrace.sstat = (outtrace.sstat as f64 + 1000.0 * tsd as f64) as i32 as i16;
            outtrace.gstat = (outtrace.gstat as f64 + 1000.0 * trd as f64) as i32 as i16;
            outtrace.tstat = (outtrace.tstat as f64 + 1000.0 * tstat as f64) as i32 as i16;
        }
        if hdrs == 3 {
            outtrace.sstat = 0;
            outtrace.gstat = 0;
            outtrace.tstat = intrace.tstat;
        }
        if hdrs == 4 {
            outtrace.sstat = 0;
            outtrace.gstat = res as i16;
            outtrace.tstat = bulk as i16;
        }
        if hdrs == 5 {
            outtrace.tstat = intrace.tstat;
            outtrace.sstat = intrace.sstat;
            outtrace.gstat = intrace.gstat;
        }
        puttr(&mut outtrace);
        if gettr(&mut intrace) == 0 {
            break;
        }
    }

    if hdrs == 3 {
        fwrite(&tcdp_stat, g1.as_ref().unwrap());
        fwrite(&num, g2.as_ref().unwrap());
        fclose(fpb.as_ref().unwrap());
        fclose(fpf.as_ref().unwrap());
        fclose(g1.as_ref().unwrap());
        fclose(g2.as_ref().unwrap());
    }
    if hdrs == 4 {
        fclose(g3.as_ref().unwrap());
        fclose(g4.as_ref().unwrap());
    }
    exit(CWP_Exit());
}
