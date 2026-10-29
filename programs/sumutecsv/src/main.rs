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
//! SUMUTECSV port of su/main/windowing_sorting_muting/sumutecsv.c.
//! Credits: SEP: Shuki Ronen, Chuck Sword; CWP: Shuki Ronen, Jack K. Cohen,
//! Dave Hale, Bjoern Rommel, Carlos E. Theodoro, Sang-yong Suh, John Stockwell;
//! DELPHI: Alexander Koek; USBG: Florian Bleibinhaus; Andre Latour.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUMUTECSV - MUTE above (or below) bilinearly interpolated polygonal curves ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sumutecsv <stdin >stdout [required parameters] [optional parameters]      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:                                                       ",
    "                                                                            ",
    " qin=             Mute functions can be input via this file.                ",
    "                  This file is optional, but if you do not input it,        ",
    "                  you must use parameters cdp=, offs=, tims=.               ",
    "                  See external document Q_FILE_STANDARDS.                   ",
    "                                                                            ",
    "                  The following 3 parameters cannot be specified if         ",
    "                  you input mute functions via the qin= file.               ",
    "                                                                            ",
    " cdp=             CDPs for which offs & tims are specified.                 ",
    " offs=            offsets corresponding to times in tims.                   ",
    " tims=            times corresponding to offsets in offs.                   ",
    "                  If any tims value is greater than 100 all times are       ",
    "                  assumed to be milliseconds. Otherwise seconds.            ",
    "            ***   If qin= is not specified, all 3 previous parameters       ",
    "                  must be specified. There must be at least 1 number        ",
    "                  in the cdp= list. There must be the same number of        ",
    "                  tims= parameters as numbers in the cdp= list.             ",
    "                  There must be the same number of offs= parameters         ",
    "                  as numbers in the cdp= list, or, there can be just        ",
    "                  one offs= parameter provided there are the same           ",
    "                  number of tims in all tims= lists.                        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rfile=           If set, read a K-file containing 3D grid definition.      ",
    "                  Assume 2D if no K-file and no grid definition is          ",
    "                  supplied via command line parameters. The required        ",
    "                  command line parameters are: grid_xa,grid_ya,             ",
    "                  grid_xb,grid_yb, grid_xc,grid_yc, grid_wb, grid_wc.       ",
    "                  (See program SUBINCSV for 3D grid details).               ",
    "                  If this is a 3D then the input CDPs for the mute          ",
    "                  locations and the seismic trace CDP key should be         ",
    "                  3D grid cell numbers as produced by SUBINCSV.             ",
    "                  A 3D also forces restrictions on the locations of         ",
    "                  the input mute locations. Their CDP locations must        ",
    "                  form aligned rectangles (see Notes).                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " offkey=offset    Key header word specifying trace offset                   ",
    " abs=1            use the absolute value of offkey.                         ",
    "               =0 do not use absolute value of offkey.                      ",
    " ntaper=0         number of samples to taper (sine-squared) from            ",
    "                  computed full-mute sample                                 ",
    " mode=0           mute ABOVE the polygonal curves                           ",
    "               =1 to zero BELOW the polygonal curves                        ",
    "               =2 to mute below AND above a straight line. In this case     ",
    "                  offs,tims describe the total time length of the muted     ",
    "                  zone as a function of offs. the slope of the line is      ",
    "                  given by vel=                                             ",
    "               =3 to mute below AND above a constant velocity hyperbola     ",
    "                  as in mode=2 offs,tims describe the total time length of  ",
    "                  mute zone as a function of offs, the velocity is given by ",
    "                  the value of vel=                                         ",
    " vel=330          constant velocity for linear or hyperbolic mute           ",
    " tzero=0          time shift (ms.) of linear or hyperbolic mute at the      ",
    "                  offkey value of 0. Note: If any tims value is greater     ",
    "                  than 100, tzero must be MILLISECONDS. Otherwise seconds.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " extrapi=0        do not extrapolate at ends in igi direction.              ",
    "                  (Mute times beyond ends are held constant).               ",
    "               =1 extrapolate at both igi ends                              ",
    "               =2 extrapolate only at lower igi end                         ",
    "               =3 extrapolate only at higher igi end                        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " extrapc=0        do not extrapolate at ends in igc direction.              ",
    "                  (Mute times beyond ends are held constant).               ",
    "               =1 extrapolate at both igc ends                              ",
    "               =2 extrapolate only at lower igc end                         ",
    "               =3 extrapolate only at higher igc end                        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " extrapt=0        do not extrapolate at ends in offset direction.           ",
    "                  (Mute times beyond ends are held constant).               ",
    "               =1 extrapolate at both offset ends                           ",
    "               =2 extrapolate only at lower offset end                      ",
    "               =3 extrapolate only at higher offset end                     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " check=0          Do not print grid checking and function locations.        ",
    "               =1 If grid definition input, run some grid functions on      ",
    "                  the 4 grid corner points and print results. Also print    ",
    "                  input mute function location information by using the     ",
    "                  grid definition and the input cdp= values.                ",
    "                  The output values are:                                    ",
    "                     G,cdp,igi,igc,xgrid,ygrid,xworld,yworld.               ",
    "                  (The format is intended to make it easy to review and     ",
    "                  also to make it easy to load to spreadsheets).            ",
    "                  This information can be written to a file by              ",
    "                  putting 2>yourfile on the command line.                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " print=0          Do not print INPUT mute functions.                        ",
    "               =1 Print INPUT mute functions. The cdp number and its        ",
    "                  input mute function values are printed. This print is     ",
    "                  not intended to be pretty. It just allows easy checking   ",
    "                  that the q-file contains the expected values.             ",
    "                  It also allows programmers and others to confirm that     ",
    "                  command line input via cdp=,offs=,tims= produces the      ",
    "                  same results as identical input via the qin= file.        ",
    "                                                                            ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The offsets specified in the offs array must be monotonically increasing.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For 3D, user needs to input mute locations (cdp numbers) which form aligned",
    " rectangles. That is, how-ever-many grid inlines the user chooses to put    ",
    " mute locations on, there must always be the same number of mute locations  ",
    " on each inline and those functions must be located at same grid crosslines.",
    " For instance, if user inputs locations for inline 7 at crosslines 15,25,40 ",
    " then the user must input locations at crosslines 15,25,40 for any other    ",
    " inlines that the user wants to supply locations for. (If user is lucky     ",
    " enough that the grid has 100 inlines, then the input locations could be at ",
    " CDPs 715,725,740 and 1115,1125,1140 and 2015,2025,2040 and 2515,2525,2540. ",
    " Note that the CDP numbers specified on the cdp= parameter and also in the  ",
    " input seismic trace cdp key are translated to inline and crossline numbers ",
    " using the input 3D grid definition - so those cdp numbers need to          ",
    " correspond to the input 3D grid definition.                                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For trace CDPs not listed in cdp= parameter or qin= file, bilinear         ",
    " interpolation is done if the trace CDP location is surrounded by 4 mute    ",
    " functions specified in the cdp= list. If the trace CDP is not surrounded   ",
    " by 4 input mute functions, the result depends on the extrapi and extrapc   ",
    " options. The result can be any combination of linear interpolation and     ",
    " linear extrapolation and constant extrapolation. If input mute functions   ",
    " are only located along 1 grid inline or 1 grid crossline, result is linear ",
    " interpolation in that direction (and linear or constant extrapolation      ",
    " the outer ending functions).                                               ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The interpolation related to cdp is actually done after the interpolation  ",
    " related to offset. That is, first the trace offset is used to compute times",
    " from the offs,tims arrays for 4, 2, or 1 mute functions and then weighting ",
    " based on the cdp locations of those mute functions is applied. Note also   ",
    " that restricting the mute to the earliest and latest trace time is done    ",
    " just before application to each trace. A consequence of this is that both  ",
    " negative offs= and negative tims= values are allowed and honored even if   ",
    " the abs= option is 1 (the default).                                        ",
    "\t\t\t\t\t\t\t\t\t     ",
];

use oxwave_bilinear_interp::*;
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;
use oxwave_q_file::*;

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut qname = String::new();
    let hasq = getparstring("qin", &mut qname) != 0;
    let mut mode = 0;
    getparint("mode", &mut mode);
    if !(0..=3).contains(&mode) {
        err!("error: mode must be 0,1,2, or 3");
    }
    let mut ntaper = 0;
    getparint("ntaper", &mut ntaper);
    if ntaper < 0 {
        err!("error: ntaper cannot be negative.");
    }
    let mut abs = 1;
    getparint("abs", &mut abs);
    let mut ei = 0;
    getparint("extrapi", &mut ei);
    if !(0..=3).contains(&ei) {
        err!("error: extrapi option not in range ");
    }
    let mut ec = 0;
    getparint("extrapc", &mut ec);
    if !(0..=3).contains(&ec) {
        err!("error: extrapc option not in range ");
    }
    let mut et = 0;
    getparint("extrapt", &mut et);
    let taper: Vec<f32> = (0..ntaper)
        .map(|k| {
            let s = (((k + 1) as f64) * std::f64::consts::PI / (2 * ntaper) as f64).sin() as f32;
            s * s
        })
        .collect();
    let mut vel = 330f32;
    getparfloat("vel", &mut vel);
    if vel == 0. {
        err!("error: vel cannot be 0");
    }
    let mut tm0 = 0f32;
    getparfloat("tzero", &mut tm0);
    let mut key = String::new();
    if getparstring("offkey", &mut key) == 0 {
        key = "offset".into();
    }
    let typ = hdtype(&key);
    let index = getindex(&key);
    let mut rname = String::new();
    let hasr = getparstring("rfile", &mut rname) != 0;
    let mut maygrid = 0;
    gridcommand(&mut maygrid);
    if maygrid == 1 && hasr {
        err!("error: input k-file not allowed when full grid on command line.");
    }
    if maygrid == -1 && !hasr {
        err!("error: input k-file required when partial grid on command line.");
    }
    let is3d = maygrid != 0 || hasr;
    let mut check = 0;
    getparint("check", &mut check);
    let mut print = 0;
    getparint("print", &mut print);
    let mut gv = [0f64; 999];
    if is3d {
        let fp = if maygrid != 1 {
            Some(
                fopen(&rname, "r")
                    .unwrap_or_else(|| err!("error: input K-file did not open correctly.")),
            )
        } else {
            None
        };
        let mut ew = 1;
        gridread(fp.as_ref(), &mut gv, &mut ew);
        if ew > 0 {
            err!("error reading grid (from K-file or command line parameters)");
        }
        gridset(&mut gv, &mut ew);
        match ew {
            1 => err!("gridset error: grid_wb cell width must be positive."),
            2 => err!("gridset error: grid_wc cell width must be positive."),
            3 => err!("gridset error: corner B is within grid_wb cell width of corner A."),
            x if x > 0 => err!("gridset error: returned with some unrecognized error code."),
            -1 => warn!("gridset warning: corner C is near A and is reset to A."),
            _ => (),
        }
        gridcheck(&gv, check, &mut ew);
        if ew > 0 {
            err!("gridcheck error: returned with some unrecognized error code.");
        }
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("error: can't get first trace");
    }
    let seismic = matches!(tr.trid, 0..=17 | 109);
    if seismic {
        if tr.dt == 0 {
            err!("dt header field must be set");
        }
    } else if tr.trid == 130 {
        if tr.d1 == 0. {
            err!("d1 header field must be set");
        }
    } else {
        err!("tr.trid = %d, unsupported trace id", tr.trid);
    }
    let mut names = vec!["cdp".to_string(), "offs".to_string(), "tims".to_string()];
    let result = if hasq {
        let fp = fopen(&qname, "r")
            .unwrap_or_else(|| err!("error: input Q-file did not open correctly."));
        getviaqfile(&fp, &names, 3, 0)
    } else {
        getviacommand(&mut names, 1, 0)
    };
    let mut q = match result {
        Ok(q) => q,
        Err(e) => {
            if hasq {
                match e {1=>err!("getqinfo error: extra C_SU_NAMES record in q-file"),2=>err!("getqinfo error: extra C_SU_NDIMS record in q-file"),3=>err!("getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),11=>err!("readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"),12=>err!("readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."),22=>err!("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),23=>err!("getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),24=>err!("getviaqfile error: C_SU_NDIMS blank where valid number expected"),25=>err!("getviaqfile error: C_SU_NDIMS non-number where valid number expected"),26=>err!("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),27=>err!("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),x if x>100=>err!("getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",x-99),_=>err!("getviaqfile error: unrecognized error code %d",e)}
            } else {
                match e {1=>err!("getviacommand error: no non-tuple name passed in."),2=>err!("getviacommand error: non-tuple names have different amounts of values."),3=>err!("getviacommand error: independent dimension parameter is empty."),4=>err!("getviacommand error: an independent dimension parameter is empty."),5=>err!("getviacommand error: members of tuple have different amounts at same location."),_=>err!("getviacommand error: returned unrecognized error number = %d",e)}
            }
        }
    };
    if q.ifixd == 2 || q.iztuple != 1 || q.ktuple != 2 {
        if hasq {
            err!("error: q-file does not contain cdp,offs,tims names (or they are in non-standard order)");
        } else {
            err!("error: command line does not contain cdp,offs,tims (or they are incorrect combination). ");
        }
    }
    checkpars();
    if q.names.first().map(|x| x.as_str()) != Some("cdp") {
        err!("error: input must have cdp (before offs,tims).");
    }
    for r in &mut q.records {
        let c = r.dlots[0].round_ties_even() as i32;
        let mut gi = c;
        let mut gc = 1;
        if is3d {
            gridcdpic(&gv, c, &mut gi, &mut gc);
            if gi < -2147483644 {
                err!("error: input cdp %d is not in grid", c);
            }
            if check > 0 {
                let (mut xg, mut yg, mut xw, mut yw) = (0., 0., 0., 0.);
                gridicgridxy(&gv, gi, gc, &mut xg, &mut yg);
                gridicrawxy(&gv, gi, gc, &mut xw, &mut yw);
                warn!(
                    "G,%12d,%12d,%12d,  %.20g,%.20g,  %.20g,%.20g",
                    c, gi, gc, xg, yg, xw, yw
                );
            }
        }
        r.kinf = vec![c, gi, gc];
    }
    if print > 0 {
        if q.ifixd == 0 {
            for r in &q.records {
                warn!(" cpd= %d   Number of offs,tims pairs= %d", r.kinf[0], r.nto);
                for i in 0..r.nto as usize {
                    warn!(" %f %f ", r.dlots[1 + r.nto as usize + i], r.dlots[1 + i]);
                }
            }
        } else {
            warn!(" Number of offs= %d", q.records[0].nto);
            for v in &q.pindepa {
                warn!(" %f ", *v);
            }
            for r in &q.records {
                warn!(" cpd= %d   Number of tims= %d", r.kinf[0], r.nto);
                for i in 0..q.records[0].nto as usize {
                    warn!(" %f ", r.dlots[1 + i]);
                }
            }
        }
    }
    for r in &q.records {
        let n = r.nto as usize;
        let offs = if q.ifixd == 0 {
            &r.dlots[1 + n..1 + 2 * n]
        } else {
            &q.pindepa[..n]
        };
        for i in 1..n {
            if offs[i - 1] >= offs[i] {
                if q.ifixd == 0 {
                    err!(
                        "error: offs values are not increasing (input cdp = %d)",
                        r.kinf[0]
                    );
                } else {
                    err!("error: offs values are not increasing.");
                }
            }
        }
    }
    let scale = if q
        .records
        .iter()
        .any(|r| r.dlots[1..1 + r.nto as usize].iter().any(|x| *x > 100.))
    {
        0.001
    } else {
        1.
    };
    tm0 = (tm0 as f64 * scale) as f32;
    q.records.sort_by_key(|r| (r.kinf[2], r.kinf[1]));
    for p in q.records.windows(2) {
        if p[0].kinf[1..3] == p[1].kinf[1..3] {
            let r = &p[1];
            err!(
                "error: Two sets of values input for cdp,igi,igc = %d %d %d",
                r.kinf[0],
                r.kinf[1],
                r.kinf[2]
            );
        }
    }
    let (mgi, mgc) = qsplit(&q.records).unwrap_or_else(|_| {
        err!("error: Input cdps do not form enclosed-rectangles (needed for bilinear).")
    });
    loop {
        let nt = tr.ns as i32;
        let mut tmin = (tr.delrt as f64 / 1000.) as f32;
        let mut dt = (tr.dt as f64 / 1000000.) as f32;
        if !seismic {
            tmin = 0.;
            dt = tr.d1;
        }
        let mut val = Value::default();
        gethval(&tr, index, &mut val);
        let mut fval = vtof(typ, val);
        if abs == 1 {
            fval = fval.abs();
        }
        let (mut wi, mut wc) = (0., 0.);
        let (mut ix, mut cx) = (1, 1);
        let mut gi = tr.cdp;
        let mut gc = 1;
        if is3d && q.records.len() > 1 {
            gridcdpic(&gv, tr.cdp, &mut gi, &mut gc);
            if gi < -2147483644 {
                err!("Error: input cdp= %d not in grid.", tr.cdp);
            }
        }
        if q.records.len() > 1 {
            binterpfind(
                gi, &mgi, ei, gc, &mgc, ec, &mut ix, &mut cx, &mut wi, &mut wc,
            );
        }
        let ni = mgi.len() as i32;
        let nc = mgc.len() as i32;
        let nd = if q.records.len() == 1 {
            1
        } else {
            ix + ni * (cx - 1)
        };
        let ndc = nd + if ni == 1 || nc == 1 { 0 } else { ni };
        let indices = if q.records.len() == 1 {
            [0, 0, 0, 0]
        } else {
            [
                (nd - 1) as usize,
                nd as usize,
                (ndc - 1) as usize,
                ndc as usize,
            ]
        };
        let arrays = indices.map(|i| {
            let r = &q.records[i];
            let n = r.nto as usize;
            let o = if q.ifixd == 0 {
                &r.dlots[1 + n..1 + 2 * n]
            } else {
                &q.pindepa[..n]
            };
            (o, &r.dlots[1..1 + n])
        });
        let mut wt = 0.;
        binterpvalue(
            fval as f64,
            et,
            arrays[0].0,
            arrays[0].1,
            arrays[1].0,
            arrays[1].1,
            ni,
            wi,
            arrays[2].0,
            arrays[2].1,
            arrays[3].0,
            arrays[3].1,
            nc,
            wc,
            &mut wt,
        );
        let t = (wt * scale) as f32;
        if mode == 0 {
            let nm = nint(((t - tmin) / dt) as f64).min(nt);
            if nm > 0 {
                tr.data[..nm as usize].fill(0.);
            }
            for i in 0..ntaper {
                if i + nm > 0 {
                    tr.data[(i + nm) as usize] *= taper[i as usize];
                }
            }
            tr.muts = if seismic {
                nint((t * 1000.) as f64) as i16
            } else {
                nint(t as f64) as i16
            };
        } else if mode == 1 {
            let nm = nint(((tmin + nt as f32 * dt - t) / dt) as f64).max(0);
            if nm > nt {
                // In C the memset starts before data, overwriting the tail of
                // the adjacent trace header as well as all trace samples.
                let preceding = (nm - nt) as usize;
                if preceding > HDRBYTES / 4 {
                    err!("error: mute extends before trace header");
                }
                tr.hdr.as_bytes_mut()[HDRBYTES - preceding * 4..].fill(0);
                tr.data[..nt as usize].fill(0.);
            } else {
                tr.data[(nt - nm) as usize..nt as usize].fill(0.);
            }
            for i in 0..ntaper {
                if nt > nm + i && nm + i > 0 {
                    tr.data[(nt - nm - 1 - i) as usize] *= taper[i as usize];
                }
            }
            tr.mute = if seismic {
                nint((t * 1000.) as f64) as i16
            } else {
                nint(t as f64) as i16
            };
        } else {
            let nm = nint(((tmin + t) / dt) as f64);
            let a = tm0 / dt;
            let b = fval / vel / dt;
            let air = if mode == 2 {
                nint((a + b) as f64)
            } else {
                nint(((a * a + b * b) as f64).sqrt())
            };
            let top = (air - nm / 2).max(0).min(nt);
            let bot = (air + nm / 2).min(nt);
            if bot >= top && top >= 0 {
                tr.data[top as usize..bot as usize].fill(0.);
            }
            for i in 0..ntaper {
                let it = air - nm / 2 - i;
                if it > 0 && (it as usize) < tr.data.len() {
                    tr.data[it as usize] *= taper[i as usize];
                }
            }
            for i in 0..ntaper {
                let it = air + nm / 2 + i;
                if it >= 0 && it < nt {
                    tr.data[it as usize] *= taper[i as usize];
                }
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
