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

//! SUNMOCSV - NMO for an arbitrary velocity function of time and 3D or 2D CDP.
//! Port of su/main/stretching_moveout_resamp/sunmocsv.c (Revision 1.2).
//!
//! Credits (from the C source):
//! SEP: Shuki Ronen, Chuck Sword
//! CWP: Shuki Ronen, Jack, Dave Hale, Bjoern Rommel
//! Modified 08/08/98: Carlos E. Theodoro, lateral offset option.
//! Modified 07/11/02: Sang-yong Suh, upward option.
//! CWP Sept 2010: John Stockwell, lateral offset and anisotropy changes.
//! Modified Sept 2021, Jan 2022, Apr 2022: Andre Latour, 3D grids,
//! library velocity readers and bilinear interpolation, millisecond times.

use oxwave_bilinear_interp::*;
use oxwave_core::cio::fopen;
use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;
use oxwave_inverse_interp::*;
use oxwave_linear_interp::*;
use oxwave_q_file::*;
use oxwave_sinc_interp::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUNMOCSV - NMO for an arbitrary velocity function of time and 3D or 2D CDP ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sunmocsv <stdin >stdout [optional parameters]\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qin=                  Velocity functions can be input via this file.       ",
    "                       This file is optional, but if you do not input it,   ",
    "                       you must use parameters cdp=, tnmo=, vnmo=.          ",
    "                       See external document Q_FILE_STANDARDS.              ",
    "\t\t\t\t\t\t\t\t\t     ",
    "                       The following 3 parameters cannot be specified if    ",
    "                       you input velocity functions via the qin= file.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cdp=\t\t        CDPs for which vnmo & tnmo are specified.            ",
    " tnmo=                 NMO times corresponding to velocities in vnmo.       ",
    "                       If any tnmo value is greater than 100 they are all   ",
    "                       assumed to be milliseconds. Otherwise seconds.       ",
    " vnmo=                 NMO velocities corresponding to times in tnmo.\t     ",
    "                       If qin= is not specified, all 3 of these parameters  ",
    "                       must be specified. There must be at least 1 number   ",
    "                       in the cdp= list. There must be the same number of   ",
    "                       vnmo= parameters as numbers in the cdp= list.        ",
    "                       There must be the same number of tnmo= parameters    ",
    "                       as numbers in the cdp= list, or, there can be just   ",
    "                       one tnmo= parameter provided there are the same      ",
    "                       number of velocities in all vnmo= lists.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " smute=1.5             samples with NMO stretch exceeding smute are zeroed  ",
    " lmute=25\t\tlength (in samples) of linear ramp for stretch mute  ",
    " sscale=1\t\t=1 to divide output samples by NMO stretch factor    ",
    " invert=0\t\t=1 to perform (approximate) inverse NMO\t\t     ",
    " upward=0\t\t=1 to scan upward to find first sample to kill\t     ",
    " voutfile=\t\tif set, interplolated velocity function v[cdp][t] is ",
    "\t\t\toutput to named file.\t\t\t     \t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rfile=                If set, read a K-file containing 3D grid definition. ",
    "                       Assume 2D if no K-file and no grid definition is     ",
    "                       supplied via command line parameters. The required   ",
    "                       command line parameters are: grid_xa,grid_ya,        ",
    "                       grid_xb,grid_yb, grid_xc,grid_yc, grid_wb, grid_wc.  ",
    "                       (See program SUBINCSV for 3D grid details).          ",
    "                       If this is a 3D then the input CDPs for the velocity ",
    "                       functions and the seismic trace CDP key should be    ",
    "                       3D grid cell numbers as produced by SUBINCSV.        ",
    "                       A 3D also forces restrictions on the locations of the",
    "                       input velocity functions. Their CDP locations must   ",
    "                       form aligned rectangles (see Notes).                 ",
    "\t\t\t\t\t\t\t\t\t     ",
    " check=0               Do not print grid checking and function locations.   ",
    "                       =1 If grid definition input, run some grid functions ",
    "                          on the 4 grid corner points and print results.    ",
    "                          Also print input velocity function location       ",
    "                          information by using the grid definition and      ",
    "                          the input cdp= values. The format and values are: ",
    "                            G,cdp,igi,igc,xgrid,ygrid,xworld,yworld.        ",
    "                          (The format is intended to make it easy to see    ",
    "                          and also to make it easy to load to spreadsheets).",
    "                          This information can be written to a file by      ",
    "                          putting 2>yourfile on the command line.           ",
    "\t\t\t\t\t\t\t\t\t     ",
    " print=0               Do not print INPUT velocity functions.               ",
    "                       =1 Print INPUT velocity functions. The cdp number    ",
    "                          and its input velocity function are printed.      ",
    "                          This print is not intended to be pretty. It just  ",
    "                          allows easy checking that the q-file contains     ",
    "                          the expected values (or just use a text editor).  ",
    "                          It also allows programmers and others to confirm  ",
    "                          that command line input via cdp=,tnmo=,vnmo=      ",
    "                          is the same as identical input via the qin= file. ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " NMO interpolation error is less than 1% for frequencies less than 60% of   ",
    " the Nyquist frequency.\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Exact inverse NMO is impossible, particularly for early times at large     ",
    " offsets and for frequencies near Nyquist with large interpolation errors.  ",
    " \t\t\t\t\t\t\t\t     \t     ",
    " The \"offset\" header field must be set.\t\t\t\t     ",
    " Use suazimuth to set offset header field when sx,sy,gx,gy are all\t     ",
    " nonzero. \t\t\t\t\t\t\t   \t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The times specified in tnmo= parameters must be monotonically increasing.  ",
    " Velocities can also be input via the qin= file (see C_QFILE_STANDARDS).    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Velocities at each input cdp are interpolated to the sample times and then ",
    " spatial interpolation is done, as explained next:                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For 3D, the user needs to input velocity functions which form aligned      ",
    " rectangles. That is, how-ever-many grid inlines the user chooses to put    ",
    " velocity functions on, there must always be the same number of functions   ",
    " on each inline and those functions must be located at same grid crosslines.",
    " For instance, if user inputs functions for inline 7 at crosslines 15,25,40 ",
    " then the user must input the functions at crosslines 15,25,40 for any other",
    " inlines that the user wants to supply functions for. (If user is lucky     ",
    " enough that the grid has 100 inlines, then the input functions could be at ",
    " CDPs 715,725,740 and 1115,1125,1140 and 2015,2025,2040 and 2515,2525,2540. ",
    " Note that the CDP numbers specified on the cdp= parameter and also in the  ",
    " input seismic trace cdp key are translated to inline and crossline numbers ",
    " using the input 3D grid definition - so those cdp numbers need to          ",
    " correspond to the input 3D grid definition.                                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For CDPs not listed in cdp= parameter or qin= file, bilinear interpolation ",
    " of 1/velocity_squared is done if CDP is surrounded by 4 input functions,   ",
    " linear interpolation of 1/velocity_squared is done if the CDP is outside   ",
    " the sides of fully surrounded area, constant extrapolation if the CDP is   ",
    " outside the corners of fully surrounded area. If functions are only input  ",
    " along 1 inline (or 1 crossline) then the result is linear interpolation    ",
    " of 1/velocity_squared along that direction (and constant beyond ends).     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The format of the output interpolated velocity file is unformatted C floats",
    " with vout[cdp][t], with time as the fast dimension and may be used as an   ",
    " input velocity file for further processing.\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut qname = String::new();
    let mut rname = String::new();
    let mut voutfile = String::new();
    let has_qname = getparstring("qin", &mut qname) != 0;
    let has_rname = getparstring("rfile", &mut rname) != 0;
    getparstring("voutfile", &mut voutfile);
    let mut smute = 1.5f32;
    getparfloat("smute", &mut smute);
    if smute <= 0.0 {
        err!("error: smute must be greater than 0.0");
    }
    let mut lmute = 25i32;
    let mut sscale = 1i32;
    let mut invert = 0i32;
    let mut upward = 0i32;
    let mut icheck = 0i32;
    let mut iprint = 0i32;
    getparint("lmute", &mut lmute);
    getparint("sscale", &mut sscale);
    getparint("invert", &mut invert);
    getparint("upward", &mut upward);
    getparint("check", &mut icheck);
    getparint("print", &mut iprint);

    /* Process and set the grid definition values. */
    let mut maygrid = 0;
    let mut gvals = [0.0f64; 999];
    gridcommand(&mut maygrid);
    if maygrid == 1 && has_rname {
        err!("error: input k-file not allowed when full grid on command line.");
    }
    if maygrid == -1 && !has_rname {
        err!("error: input k-file required when partial grid on command line.");
    }
    let is3d = !(maygrid == 0 && !has_rname);
    if is3d {
        let fpr = if maygrid != 1 {
            match fopen(&rname, "r") {
                Some(fp) => Some(fp),
                None => err!("error: input K-file did not open correctly."),
            }
        } else {
            None
        };
        let mut errwarn = 1;
        gridread(fpr.as_ref(), &mut gvals, &mut errwarn);
        if errwarn > 0 {
            err!("error reading grid (from K-file or command line parameters)");
        }
        gridset(&mut gvals, &mut errwarn);
        match errwarn {
            1 => err!("gridset error: grid_wb cell width must be positive."),
            2 => err!("gridset error: grid_wc cell width must be positive."),
            3 => err!("gridset error: corner B is within grid_wb cell width of corner A."),
            n if n > 0 => err!("gridset error: returned with some unrecognized error code."),
            -1 => warn!("gridset warning: corner C is near A and is reset to A."),
            _ => {}
        }
        gridcheck(&gvals, icheck, &mut errwarn);
        if errwarn > 0 {
            err!("gridcheck error: returned with some unrecognized error code.");
        }
    }

    /* Get information from the first header. */
    if gettr(&mut tr) == 0 {
        err!("error: can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut ovvt = vec![0.0f32; nt];
    // The C reads ttn[1] even at nt=1. A zero spare allocation element
    // reproduces the original allocator's observed value for that case.
    let mut ttn = vec![0.0f32; nt.max(2)];
    let mut atn = vec![0.0f32; nt.max(2)];
    let mut qtn = vec![0.0f32; nt];
    // Inverse NMO sets itmute=-1 when nt=1. Keep an explicit leading
    // allocation word for the C accesses at [-1], observed as zero in SU.
    // These three arrays use logical index i at Rust index i+1.
    let mut tnt = vec![0.0f32; nt + 1];
    let mut at = vec![0.0f32; nt + 1];
    let mut qt = vec![0.0f32; nt + 1];
    let mut dind = vec![0.0f64; nt];
    let mut dswap = vec![0.0f64; nt];
    let voutfp = if !voutfile.is_empty() {
        match fopen(&voutfile, "w") {
            Some(fp) => Some(fp),
            None => err!("error: cannot open voutfile=%s\n", voutfile),
        }
    } else {
        None
    };
    for i in 0..nt {
        dind[i] = (ft + i as f32 * dt) as f64;
    }

    /* Set parameters and names of the three values needed from the readers. */
    let mut pname = vec!["cdp".to_owned(), "tnmo".to_owned(), "vnmo".to_owned()];
    let mut input = if !has_qname {
        match getviacommand(&mut pname, 1, nt as i32) {
            Ok(v) => v,
            Err(1) => err!("getviacommand error: no non-tuple name passed in."),
            Err(2) => {
                err!("getviacommand error: non-tuple names have different amounts of values.")
            }
            Err(3) => err!("getviacommand error: independent dimension parameter is empty."),
            Err(4) => err!("getviacommand error: an independent dimension parameter is empty."),
            Err(5) => err!(
                "getviacommand error: members of tuple have different amounts at same location."
            ),
            Err(n) => err!(
                "getviacommand error: returned unrecognized error number = %d",
                n
            ),
        }
    } else {
        let fpq = match fopen(&qname, "r") {
            Some(fp) => fp,
            None => err!("error: input Q-file did not open correctly."),
        };
        match getviaqfile(&fpq, &pname, 3, nt as i32) {
            Ok(v) => v,
            Err(1) => err!("getqinfo error: extra C_SU_NAMES record in q-file"),
            Err(2) => err!("getqinfo error: extra C_SU_NDIMS record in q-file"),
            Err(3) => err!(
                "getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."
            ),
            Err(11) => err!(
                "readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"
            ),
            Err(12) => err!(
                "readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."
            ),
            Err(22) => {
                err!("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.")
            }
            Err(23) => err!("getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
            Err(24) => err!("getviaqfile error: C_SU_NDIMS blank where valid number expected"),
            Err(25) => err!("getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
            Err(26) => {
                err!("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple")
            }
            Err(27) => {
                err!("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.")
            }
            Err(n) if n > 100 => err!(
                "getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",
                n - 99
            ),
            Err(n) => err!("getviaqfile error: unrecognized error code %d", n),
        }
    };
    if input.ifixd == 2 || input.iztuple != 1 || input.ktuple != 2 {
        if !has_qname {
            err!("error: command line does not contain cdp,tnmo,vnmo (or they are incorrect combination). ");
        } else {
            err!("error: q-file does not contain cdp,tnmo,vnmo names (or they are in non-standard order)");
        }
    }
    checkpars();
    let iztuple = input.iztuple as usize;
    let ifixd = input.ifixd;
    let mut pindepa = std::mem::take(&mut input.pindepa);
    let mut vinfo = std::mem::take(&mut input.records);
    let ncdp = vinfo.len();
    let jnloc1 = input.names[..iztuple].iter().position(|n| n == "cdp");
    if jnloc1 != Some(0) {
        err!("error: input must have cdp (before tnmo,vnmo).");
    }

    /* Sort on inline/crossline for 3D and on CDP/1 for 2D. */
    for rec in &mut vinfo {
        rec.kinf = vec![0; 3];
        rec.kinf[0] = rec.dlots[0].round_ties_even() as i64 as i32;
        if is3d {
            let (mut igi, mut igc) = (0, 0);
            gridcdpic(&gvals, rec.kinf[0], &mut igi, &mut igc);
            rec.kinf[1] = igi;
            rec.kinf[2] = igc;
            if rec.kinf[1] < -2147483644 {
                err!("error: input cdp %d is not in grid", rec.kinf[0]);
            }
        } else {
            rec.kinf[1] = rec.kinf[0];
            rec.kinf[2] = 1;
        }
    }

    if iprint > 0 {
        if ifixd == 0 {
            for rec in &vinfo {
                warn!(
                    " cpd= %d   Number of time,vel pairs= %d",
                    rec.kinf[0], rec.nto
                );
                let times = &rec.dlots[iztuple + rec.nto as usize..];
                for i in 0..rec.nto as usize {
                    warn!(" %f %f ", times[i], rec.dlots[iztuple + i]);
                }
            }
        } else if ifixd == 1 {
            warn!(" Number of times= %d", vinfo[0].nto);
            for i in 0..vinfo[0].nto as usize {
                warn!(" %f ", pindepa[i]);
            }
            for rec in &vinfo {
                warn!(" cpd= %d   Number of vels= %d", rec.kinf[0], rec.nto);
                for i in 0..vinfo[0].nto as usize {
                    warn!(" %f ", rec.dlots[iztuple + i]);
                }
            }
        }
    }

    /* Check increasing times, then convert all times if any exceed 100. */
    for rec in &vinfo {
        // Empty tuples make the C last-time indexing undefined.
        if rec.nto <= 0 {
            err!(
                "error: no time/velocity pairs (input cdp = %d)",
                rec.kinf[0]
            );
        }
    }
    let mut milliseconds = false;
    if ifixd == 0 {
        for rec in &vinfo {
            let times = &rec.dlots[iztuple + rec.nto as usize..iztuple + 2 * rec.nto as usize];
            for i in 1..rec.nto as usize {
                if times[i - 1] >= times[i] {
                    err!(
                        "error: tnmo values are not increasing (input cdp = %d)",
                        rec.kinf[0]
                    );
                }
            }
            if times[rec.nto as usize - 1] > 100.0 {
                milliseconds = true;
            }
        }
    } else if ifixd == 1 {
        for i in 1..vinfo[0].nto as usize {
            if pindepa[i - 1] >= pindepa[i] {
                err!("error: tnmo values are not increasing.");
            }
        }
        if pindepa[vinfo[0].nto as usize - 1] > 100.0 {
            milliseconds = true;
        }
    }
    if milliseconds {
        if ifixd == 0 {
            for rec in &mut vinfo {
                for i in 0..rec.nto as usize {
                    rec.dlots[iztuple + rec.nto as usize + i] *= 0.001;
                }
            }
        } else if ifixd == 1 {
            for i in 0..vinfo[0].nto as usize {
                pindepa[i] *= 0.001;
            }
        }
    }

    vinfo.sort_by(comp_sort2);
    for jcdp in 1..ncdp {
        if vinfo[jcdp - 1].kinf[1] == vinfo[jcdp].kinf[1]
            && vinfo[jcdp - 1].kinf[2] == vinfo[jcdp].kinf[2]
        {
            err!(
                "error: Two sets of values input for cdp,igi,igc = %d %d %d",
                vinfo[jcdp].kinf[0],
                vinfo[jcdp].kinf[1],
                vinfo[jcdp].kinf[2]
            );
        }
    }
    let (mgi, mgc) = match qsplit(&vinfo) {
        Ok(v) => v,
        Err(_) => err!("error: Input cdps do not form enclosed-rectangles (needed for bilinear)."),
    };
    let mgi_tot = mgi.len() as i32;
    let mgc_tot = mgc.len() as i32;
    let mgi_totdeg = if mgi_tot == 1 || mgc_tot == 1 {
        0
    } else {
        mgi_tot
    };

    /* Interpolate in time and square-and-inverse the velocities. */
    for rec in &mut vinfo {
        let times = if ifixd == 0 {
            &rec.dlots[iztuple + rec.nto as usize..iztuple + 2 * rec.nto as usize]
        } else {
            &pindepa[..]
        };
        for it in 0..nt {
            linterpd(
                dind[it],
                times,
                &rec.dlots[iztuple..iztuple + rec.nto as usize],
                0,
                &mut dswap[it],
            );
            if dswap[it] > 1.0e-98 {
                dswap[it] = 1.0 / (dswap[it] * dswap[it]);
            } else if dswap[it] < -1.0e-98 {
                dswap[it] = -1.0 / (dswap[it] * dswap[it]);
            }
        }
        for it in 0..nt {
            rec.dlots[it + iztuple] = dswap[it];
        }
    }

    let mut oldcdp = tr.cdp as i64 - 1;
    let mut oldoffset = tr.offset as i64 - 1;
    let mut newsloth = false;
    let mut itmute = 0i32;
    let mut kigc = 1;
    if ncdp < 2 {
        for it in 0..nt {
            ovvt[it] = vinfo[0].dlots[it + iztuple] as f32;
        }
        if let Some(fp) = &voutfp {
            for it in 0..nt {
                let veltemp = ovvt[it];
                let voutt = [(1.0 / veltemp as f64).sqrt() as f32];
                efwrite(&voutt, fp);
            }
        }
    }

    loop {
        if tr.cdp as i64 != oldcdp {
            if ncdp > 1 {
                let kigi;
                if is3d {
                    let mut igi = 0;
                    gridcdpic(&gvals, tr.cdp, &mut igi, &mut kigc);
                    if igi < -2147483644 {
                        err!("Error: input cdp= %d not in grid.", tr.cdp);
                    }
                    kigi = igi;
                } else {
                    kigi = tr.cdp;
                }
                let (mut mgix, mut mgcx, mut wi, mut wc) = (0, 0, 0.0, 0.0);
                binterpfind(
                    kigi, &mgi, 0, kigc, &mgc, 0, &mut mgix, &mut mgcx, &mut wi, &mut wc,
                );
                let ndxi = mgix + mgi_tot * (mgcx - 1);
                let ndxc = ndxi + mgi_totdeg;
                let mdxi = ndxi - 1;
                let mdxc = ndxc - 1;
                binterpapply(
                    &vinfo[mdxi as usize].dlots[iztuple..],
                    &vinfo[ndxi as usize].dlots[iztuple..],
                    mgi_tot,
                    wi,
                    &vinfo[mdxc as usize].dlots[iztuple..],
                    &vinfo[ndxc as usize].dlots[iztuple..],
                    mgc_tot,
                    wc,
                    &mut dswap,
                );
                for it in 0..nt {
                    ovvt[it] = dswap[it] as f32;
                }
                newsloth = true;
                if let Some(fp) = &voutfp {
                    for it in 0..nt {
                        let veltemp = ovvt[it];
                        let voutt = [(1.0 / veltemp as f64).sqrt() as f32];
                        efwrite(&voutt, fp);
                    }
                }
            }
        } else {
            newsloth = false;
        }

        if newsloth || tr.offset as i64 != oldoffset {
            let offset = tr.offset as f32;
            let temp = (offset * offset) / (dt * dt);
            let mut tn = ft / dt;
            for it in 0..nt {
                let tsq = temp * ovvt[it];
                ttn[it] = ((tn * tn + tsq) as f64).sqrt() as f32;
                tn = (tn as f64 + 1.0) as f32;
            }
            atn[0] = ttn[1] - ttn[0];
            for it in 1..nt {
                atn[it] = ttn[it] - ttn[it - 1];
            }
            let osmute = (1.0 / smute as f64) as f32;
            let mut it;
            if upward == 0 {
                it = 0i32;
                while it < nt as i32 - 1 && atn[it as usize] < osmute {
                    it += 1;
                }
            } else {
                it = nt as i32 - 1;
                while it > 0 && atn[it as usize] >= osmute {
                    it -= 1;
                }
            }
            itmute = it;
            if invert != 0 {
                let start = itmute as usize;
                yxtoxy(
                    1.0,
                    ft / dt + itmute as f32,
                    &ttn[start..nt],
                    1.0,
                    ft / dt + itmute as f32,
                    ft / dt - nt as f32,
                    ft / dt + nt as f32,
                    &mut tnt[start + 1..nt + 1],
                );
                itmute = (1.0 + ttn[start] as f64 - (ft / dt) as f64) as i32;
                itmute = if (nt as i32 - 2) < itmute {
                    nt as i32 - 2
                } else {
                    itmute
                };
                if sscale != 0 {
                    for it in itmute + 1..nt as i32 {
                        at[(it + 1) as usize] = tnt[(it + 1) as usize] - tnt[it as usize];
                    }
                    at[(itmute + 1) as usize] = at[(itmute + 2) as usize];
                }
            }
        }
        if invert == 0 {
            let start = itmute as usize;
            ints8r(
                1.0,
                ft / dt,
                &tr.data[..nt],
                0.0,
                0.0,
                &ttn[start..nt],
                &mut qtn[start..nt],
            );
            for it in 0..start {
                qtn[it] = 0.0;
            }
            let mut it = itmute;
            while it < itmute.wrapping_add(lmute) && it < nt as i32 {
                qtn[it as usize] *= (it - itmute + 1) as f32 / lmute as f32;
                it += 1;
            }
            if sscale != 0 {
                for it in start..nt {
                    qtn[it] *= atn[it];
                }
            }
            tr.data[..nt].copy_from_slice(&qtn);
        } else {
            let start = (itmute + 1) as usize;
            ints8r(
                1.0,
                ft / dt,
                &tr.data[..nt],
                0.0,
                0.0,
                &tnt[start..nt + 1],
                &mut qt[start..nt + 1],
            );
            for it in 0..itmute {
                qt[(it + 1) as usize] = 0.0;
            }
            if sscale != 0 {
                for it in start..nt + 1 {
                    qt[it] *= at[it];
                }
            }
            tr.data[..nt].copy_from_slice(&qt[1..nt + 1]);
        }
        puttr(&mut tr);
        oldoffset = tr.offset as i64;
        oldcdp = tr.cdp as i64;
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    if let Some(fp) = voutfp {
        efclose(&fp);
    }
    exit(CWP_Exit());
}

/* Compare crossline first, then inline, just as compSort2 in the C. */
fn comp_sort2(p1: &QInfo, p2: &QInfo) -> std::cmp::Ordering {
    if p1.kinf[2] < p2.kinf[2] {
        return std::cmp::Ordering::Less;
    }
    if p1.kinf[2] > p2.kinf[2] {
        return std::cmp::Ordering::Greater;
    }
    if p1.kinf[1] < p2.kinf[1] {
        return std::cmp::Ordering::Less;
    }
    if p1.kinf[1] > p2.kinf[1] {
        return std::cmp::Ordering::Greater;
    }
    std::cmp::Ordering::Equal
}
