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

//! SUNMO - NMO for an arbitrary velocity function of time and CDP.
//! Port of su/main/stretching_moveout_resamp/sunmo.c (Revision 1.31).
//! Credits: SEP: Shuki Ronen, Chuck Sword; CWP: Shuki Ronen, Jack, Dave Hale, Bjoern Rommel; Carlos E. Theodoro, Sang-yong Suh, John Stockwell.

use oxwave_core::prelude::*;
use oxwave_index_search::xindex;
use oxwave_inverse_interp::yxtoxy;
use oxwave_linear_interp::intlin;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUNMO - NMO for an arbitrary velocity function of time and CDP\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sunmo <stdin >stdout [optional parameters]\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " tnmo=0,...\t\tNMO times corresponding to velocities in vnmo\t     ",
    " vnmo=1500,...\t\tNMO velocities corresponding to times in tnmo\t     ",
    " cdp=\t\t\tCDPs for which vnmo & tnmo are specified (see Notes) ",
    " smute=1.5\t\tsamples with NMO stretch exceeding smute are zeroed  ",
    " lmute=25\t\tlength (in samples) of linear ramp for stretch mute  ",
    " sscale=1\t\t=1 to divide output samples by NMO stretch factor    ",
    " invert=0\t\t=1 to perform (approximate) inverse NMO\t\t     ",
    " upward=0\t\t=1 to scan upward to find first sample to kill\t     ",
    " voutfile=\t\tif set, interplolated velocity function v[cdp][t] is ",
    "\t\t\toutput to named file.\t\t\t     \t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    " For constant-velocity NMO, specify only one vnmo=constant and omit tnmo.   ",
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
    " For NMO with a velocity function of time only, specify the arrays\t     ",
    "\t   vnmo=v1,v2,... tnmo=t1,t2,...\t\t\t\t     ",
    " where v1 is the velocity at time t1, v2 is the velocity at time t2, ...    ",
    " The times specified in the tnmo array must be monotonically increasing.    ",
    " Linear interpolation and constant extrapolation of the specified velocities",
    " is used to compute the velocities at times not specified.\t\t     ",
    " The same holds for the anisotropy coefficients as a function of time only. ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For NMO with a velocity function of time and CDP, specify the array\t     ",
    "\t   cdp=cdp1,cdp2,...\t\t\t\t\t\t     ",
    " and, for each CDP specified, specify the vnmo and tnmo arrays as described ",
    " above. The first (vnmo,tnmo) pair corresponds to the first cdp, and so on. ",
    " Linear interpolation and constant extrapolation of 1/velocity^2 is used    ",
    " to compute velocities at CDPs not specified.\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The format of the output interpolated velocity file is unformatted C floats",
    " with vout[cdp][t], with time as the fast dimension and may be used as an   ",
    " input velocity file for further processing.\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Note that this version of sunmo does not attempt to deal with\tanisotropy.  ",
    " The version of sunmo with experimental anisotropy support is \"sunmo_a\"   ",
    "\t\t\t\t\t\t\t\t\t     ",
];

fn interpovv(nt: usize, cdp: &[f32], ovv: &[Vec<f32>], cdpt: f32, ovvt: &mut [f32], index: &mut i32) {
    if cdpt <= cdp[0] { ovvt[..nt].copy_from_slice(&ovv[0][..nt]); }
    else if cdpt >= cdp[cdp.len()-1] { ovvt[..nt].copy_from_slice(&ovv[cdp.len()-1][..nt]); }
    else {
        xindex(cdp, cdpt, index);
        let i = *index as usize;
        let a1 = (cdp[i+1]-cdpt)/(cdp[i+1]-cdp[i]);
        let a2 = (cdpt-cdp[i])/(cdp[i+1]-cdp[i]);
        for it in 0..nt { ovvt[it] = a1*ovv[i][it]+a2*ovv[i+1][it]; }
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut voutfile = String::new();
    getparstring("voutfile", &mut voutfile);
    let voutfp = if !voutfile.is_empty() {
        match fopen(&voutfile, "w") { Some(fp) => Some(fp), None => { err!("cannot open voutfile=%s\n", voutfile); } }
    } else { None };
    let mut ncdp = countparval("cdp");
    if ncdp > 0 {
        if countparname("vnmo") != ncdp { err!("a vnmo array must be specified for each cdp"); }
        if countparname("tnmo") != ncdp { err!("a tnmo array must be specified for each cdp"); }
    } else {
        ncdp = 1;
        if countparname("vnmo") > 1 { err!("only one (or no) vnmo array must be specified"); }
        if countparname("tnmo") > 1 { err!("only one (or no) tnmo array must be specified"); }
    }
    let mut cdp = vec![0.0; ncdp as usize];
    if getparfloat("cdp", &mut cdp) == 0 { cdp[0] = tr.cdp as f32; }
    let mut ovv = vec![vec![0.0; nt]; ncdp as usize];
    for icdp in 0..ncdp as usize {
        let mut nvnmo = countnparval(icdp as i32 + 1, "vnmo");
        let mut ntnmo = countnparval(icdp as i32 + 1, "tnmo");
        if nvnmo != ntnmo && !(ncdp == 1 && nvnmo == 1 && ntnmo == 0) { err!("number of vnmo and tnmo values must be equal"); }
        if nvnmo == 0 { nvnmo = 1; }
        if ntnmo == 0 { ntnmo = nvnmo; }
        let mut vnmo = vec![0.0; nvnmo as usize];
        let mut tnmo = vec![0.0; nvnmo as usize];
        if getnparfloat(icdp as i32 + 1, "vnmo", &mut vnmo) == 0 { vnmo[0] = 1500.0; }
        if getnparfloat(icdp as i32 + 1, "tnmo", &mut tnmo) == 0 { tnmo[0] = 0.0; }
        for it in 1..ntnmo as usize { if tnmo[it] <= tnmo[it-1] { err!("tnmo values must increase monotonically"); } }
        let mut tn = ft;
        for it in 0..nt {
            let mut v = [0.0];
            intlin(&tnmo, &vnmo, vnmo[0], vnmo[nvnmo as usize-1], &[tn], &mut v);
            ovv[icdp][it] = (1.0f64 / (v[0]*v[0]) as f64) as f32;
            tn += dt;
        }
    }
    for j in 1..ncdp as usize {
        let mut i = j;
        while i > 0 && cdp[i-1] > cdp[i] { cdp.swap(i-1,i); ovv.swap(i-1,i); i -= 1; }
    }
    let mut smute = 1.5f32; getparfloat("smute", &mut smute);
    if smute <= 0.0 { err!("smute must be greater than 0.0"); }
    let mut lmute = 25; getparint("lmute", &mut lmute);
    let mut sscale = 1; getparint("sscale", &mut sscale);
    let mut invert = 0; getparint("invert", &mut invert);
    let mut upward = 0; getparint("upward", &mut upward);
    checkpars();
    let mut ovvt = vec![0.0; nt];
    let mut ttn = vec![0.0; nt.max(2)];
    let mut atn = vec![0.0; nt.max(2)];
    let mut qtn = vec![0.0; nt];
    let mut tnt = vec![0.0; nt.max(2)];
    let mut at = vec![0.0; nt.max(2)];
    let mut qt = vec![0.0; nt];
    let mut index = 0;
    interpovv(nt, &cdp, &ovv, tr.cdp as f32, &mut ovvt, &mut index);
    if let Some(ref fp) = voutfp { for &v in &ovvt { efwrite(&[((1.0f64 / v as f64).sqrt()) as f32], fp); } }
    let mut oldcdp = tr.cdp;
    let mut oldoffset = (tr.offset as i64) - 1;
    let mut itmute = 0usize;
    loop {
        let newsloth = tr.cdp != oldcdp && ncdp > 1;
        if newsloth {
            interpovv(nt, &cdp, &ovv, tr.cdp as f32, &mut ovvt, &mut index);
            if let Some(ref fp) = voutfp { for &v in &ovvt { efwrite(&[((1.0f64 / v as f64).sqrt()) as f32], fp); } }
        }
        if newsloth || tr.offset as i64 != oldoffset {
            let offset = tr.offset as f32;
            let temp = offset*offset/(dt*dt);
            let mut tn = ft/dt;
            for it in 0..nt {
                let tsq = temp*ovvt[it];
                ttn[it] = ((tn*tn + tsq) as f64).sqrt() as f32;
                tn = (tn as f64 + 1.0) as f32;
            }
            atn[0] = ttn[1]-ttn[0];
            for it in 1..nt { atn[it] = ttn[it]-ttn[it-1]; }
            let osmute = (1.0f64/smute as f64) as f32;
            let mut it = if upward == 0 { 0 } else { nt.saturating_sub(1) };
            if upward == 0 { while it < nt.saturating_sub(1) && atn[it] < osmute { it += 1; } }
            else { while it > 0 && atn[it] >= osmute { it -= 1; } }
            itmute = it;
            if invert != 0 && nt >= 2 {
                let base = ft/dt;
                let start = (base as f64 + itmute as f64) as f32;
                yxtoxy(1.0, start, &ttn[itmute..nt], 1.0, start, base - nt as f32, base + nt as f32, &mut tnt[itmute..nt]);
                let m = (1.0f64 + ttn[itmute] as f64 - (ft/dt) as f64) as i32;
                itmute = m.min(nt as i32 - 2).max(0) as usize;
                if sscale != 0 { for it in itmute+1..nt { at[it] = tnt[it]-tnt[it-1]; } at[itmute] = at[itmute+1]; }
            }
        }
        if invert == 0 {
            ints8r(1.0, ft/dt, &tr.data[..nt], 0.0, 0.0, &ttn[itmute..nt], &mut qtn[itmute..nt]);
            qtn[..itmute].fill(0.0);
            // Compare the ramp length directly so a large lmute cannot overflow
            // when added to the mute index.
            for it in itmute..nt { if ((it - itmute) as i64) < lmute as i64 { qtn[it] *= (it-itmute+1) as f32 / lmute as f32; } }
            if sscale != 0 { for it in itmute..nt { qtn[it] *= atn[it]; } }
            tr.data[..nt].copy_from_slice(&qtn);
        } else {
            ints8r(1.0, ft/dt, &tr.data[..nt], 0.0, 0.0, &tnt[itmute..nt], &mut qt[itmute..nt]);
            qt[..itmute].fill(0.0);
            if sscale != 0 { for it in itmute..nt { qt[it] *= at[it]; } }
            tr.data[..nt].copy_from_slice(&qt);
        }
        puttr(&mut tr);
        oldoffset = tr.offset as i64; oldcdp = tr.cdp;
        if gettr(&mut tr) == 0 { break; }
    }
    if let Some(fp) = voutfp { efclose(&fp); }
    exit(CWP_Exit());
}
