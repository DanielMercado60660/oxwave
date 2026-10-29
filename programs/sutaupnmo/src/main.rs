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
//! SUTAUPNMO - NMO for an arbitrary velocity function of tau and CDP.
//! Port of su/main/stretching_moveout_resamp/sutaupnmo.c, revision 1.3.
//! Credits: Durham, Richard Hobbs; SEP: Shuki Ronen, Chuck Sword;
//! CWP: Shuki Ronen, Jack K. Cohen, Dave Hale.

use oxwave_core::prelude::*;
use oxwave_index_search::xindex;
use oxwave_linear_interp::intlin;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUTAUPNMO - NMO for an arbitrary velocity function of tau and CDP\t",
    "\t\t\t\t\t\t\t\t\t",
    "  sutaupnmo <stdin >stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " tnmo=0,...\t\tNMO times corresponding to velocities in vnmo\t",
    " vnmo=1500,...\t\tNMO velocities corresponding to times in tnmo\t",
    " cdp=\t\t\tCDPs for which vnmo & tnmo are specified (see Notes) ",
    " smute=1.5\t\tsamples with NMO stretch exceeding smute are zeroed  ",
    " lmute=25\t\tlength (in samples) of linear ramp for stretch mute  ",
    " sscale=1\t\t=1 to divide output samples by NMO stretch factor    ",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For constant-velocity NMO, specify only one vnmo=constant and omit tnmo.",
    "\t\t\t\t\t\t\t\t\t",
    " For NMO with a velocity function of tau only, specify the arrays\t",
    "\t   vnmo=v1,v2,... tnmo=t1,t2,...\t\t\t\t",
    " where v1 is the velocity at tau t1, v2 is the velocity at tau t2, ...    ",
    " The taus specified in the tnmo array must be monotonically increasing.    ",
    " Linear interpolation and constant extrapolation of the specified velocities",
    " is used to compute the velocities at taus not specified.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For NMO with a velocity function of tau and CDP, specify the array\t",
    "\t   cdp=cdp1,cdp2,...\t\t\t\t\t\t",
    " and, for each CDP specified, specify the vnmo and tnmo arrays as described ",
    " above. The first (vnmo,tnmo) pair corresponds to the first cdp, and so on. ",
    " Linear interpolation and constant extrapolation of velocity^2 is used\t ",
    " to compute velocities at CDPs not specified.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Moveout is defined by\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  tau^2 + tau^2.p^2.vel^2\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: In general, the user should set the cdp parameter.  The default is   ",
    "\tto use tr.cdp from the first trace and assume only one cdp.\t ",
    " Caveat:\t\t\t\t\t\t\t\t",
    " Taunmo should handle triplication\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " NMO interpolation error is less than 1% for frequencies less than 60% of   ",
    " the Nyquist frequency.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Exact inverse NMO is not implemented, nor has anisotropy\t\t",
    " Example implementation:\t\t\t\t\t\t",
    "   sutaup dx=25 option=2 pmin=0 pmax=0.0007025 < cmpgather.su |\t",
    "   supef minlag=0.2 maxlag=0.8 |\t\t\t\t\t",
    "   sutaupnmo tnmo=0.5,2,4 vnmo=1500,2000,3200 smute=1.5 |\t\t",
    "   sumute key=tracr mode=1 ntaper=20 xmute=1,30,40,50,85,15  \t\t",
    "\t\t\t\t tmute=7.8,7.8,4.5,3.5,2.0,0.35 |\t",
    "   sustack key=cdp | ... [...]\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn interpvv(nt: usize, cdp: &[f32], vv: &[Vec<f32>], cdpt: f32, vvt: &mut [f32], indx: &mut i32) {
    let ncdp = cdp.len();
    if cdpt <= cdp[0] { vvt[..nt].copy_from_slice(&vv[0]); }
    else if cdpt >= cdp[ncdp-1] { vvt[..nt].copy_from_slice(&vv[ncdp-1]); }
    else {
        xindex(cdp, cdpt, indx);
        let i = *indx as usize;
        let a1 = (cdp[i+1]-cdpt)/(cdp[i+1]-cdp[i]);
        let a2 = (cdpt-cdp[i])/(cdp[i+1]-cdp[i]);
        for it in 0..nt { vvt[it] = a1*vv[i][it]+a2*vv[i+1][it]; }
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
    let mut vv = vec![vec![0.0; nt]; ncdp as usize];
    for icdp in 0..ncdp as usize {
        let mut nvnmo = countnparval(icdp as i32+1, "vnmo");
        let mut ntnmo = countnparval(icdp as i32+1, "tnmo");
        if nvnmo != ntnmo && !(ncdp == 1 && nvnmo == 1 && ntnmo == 0) { err!("number of vnmo and tnmo values must be equal"); }
        if nvnmo == 0 { nvnmo = 1; }
        if ntnmo == 0 { ntnmo = nvnmo; }
        let mut vnmo = vec![0.0; nvnmo as usize];
        let mut tnmo = vec![0.0; nvnmo as usize];
        if getnparfloat(icdp as i32+1, "vnmo", &mut vnmo) == 0 { vnmo[0] = 1500.0; }
        if getnparfloat(icdp as i32+1, "tnmo", &mut tnmo) == 0 { tnmo[0] = 0.0; }
        for it in 1..ntnmo as usize { if tnmo[it] <= tnmo[it-1] { err!("tnmo values must increase monotonically"); } }
        let mut tn = ft;
        for it in 0..nt {
            let mut v = [0.0];
            intlin(&tnmo[..ntnmo as usize], &vnmo[..ntnmo as usize], vnmo[0], vnmo[nvnmo as usize-1], &[tn], &mut v);
            vv[icdp][it] = v[0]*v[0];
            tn += dt;
        }
    }
    for j in 1..ncdp as usize {
        let mut i = j;
        while i > 0 && cdp[i-1] > cdp[i] { cdp.swap(i-1,i); vv.swap(i-1,i); i -= 1; }
    }
    let mut smute = 1.5f32;
    getparfloat("smute", &mut smute);
    if smute <= 0.0 { err!("smute must be greater than 0.0"); }
    let mut lmute = 25i32;
    getparint("lmute", &mut lmute);
    let mut sscale = 1i32;
    getparint("sscale", &mut sscale);
    checkpars();
    let mut vvt = vec![0.0; nt];
    // The original reads ttn[1] even when nt=1; reserve a zeroed extra slot.
    let mut ttn = vec![0.0; nt.max(2)];
    let mut atn = vec![0.0; nt.max(1)];
    let mut qtn = vec![0.0; nt];
    let mut indx = 0;
    interpvv(nt, &cdp, &vv, tr.cdp as f32, &mut vvt, &mut indx);
    let mut oldcdp = tr.cdp;
    let mut oldp = -99999.0f32;
    let mut itmute = 0usize;
    loop {
        let newvel2 = tr.cdp != oldcdp && ncdp > 1;
        if newvel2 { interpvv(nt, &cdp, &vv, tr.cdp as f32, &mut vvt, &mut indx); }
        let newp = (tr.f2 as f64 + (tr.tracr as f32 as f64 - 1.0)*tr.d2 as f64) as f32;
        if newvel2 || newp != oldp {
            let mut tn = ft/dt;
            for it in 0..nt {
                let tsq = tn*tn - newp*newp*tn*tn*vvt[it];
                ttn[it] = if tsq < 0.0 { 0.0 } else { (tsq as f64).sqrt() as f32 };
                tn = (tn as f64 + 1.0) as f32;
            }
            atn[0] = ttn[1]-ttn[0];
            for it in 1..nt { atn[it] = ttn[it]-ttn[it-1]; }
            let osmute = (1.0/smute as f64) as f32;
            itmute = 0;
            while itmute+1 < nt && atn[itmute] > osmute { itmute += 1; }
        }
        ints8r(1.0, ft/dt, &tr.data[..nt], 0.0, 0.0, &ttn[..itmute], &mut qtn[..itmute]);
        for it in itmute..nt { qtn[it] = 0.0; }
        // Negative indices in the original ramp address memory outside qtn.
        let start = (itmute as i64-lmute as i64).max(0) as usize;
        for it in start..itmute.min(nt) { qtn[it] *= (itmute-it) as f32/lmute as f32; }
        if sscale != 0 { for it in itmute..nt { qtn[it] *= atn[it]; } }
        tr.data[..nt].copy_from_slice(&qtn);
        puttr(&mut tr);
        oldp = newp;
        oldcdp = tr.cdp;
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
