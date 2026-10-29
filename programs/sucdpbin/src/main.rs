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

//! SUCDPBIN - Compute CDP bin number. Port of su/main/headers/sucdpbin.c.
//! Credits: Dominique Rousset and Mohamed Hamza, Université de Pau et des Pays de l'Adour.
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUCDPBIN - Compute CDP bin number\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " sucdpbin <stdin >stdout xline= yline= dcdp=\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " xline=\t\tarray of X defining the CDP line\t\t",
    " yline=\t\tarray of Y defining the CDP line\t\t",
    " dcdp=\t\t\tdistance between bin centers\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters\t\t\t\t\t\t\t",
    " verbose=0\t\t<>0 output informations\t\t\t\t",
    " cdpmin=1001\t\tmin cdp bin number\t\t\t\t",
    " distmax=dcdp\t\tsearch radius\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " xline,yline defines the CDP line made of continuous straight lines. \t",
    " If a smoother line is required, use unisam to interpolate.\t\t",
    " Bin centers are located at dcdp constant interval on this line. \t",
    " Each trace will be numbered with the number of the closest bin. If no  ",
    " bin center is found within the search radius. cdp is set to 0\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

use oxwave_core::prelude::*;

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv); requestdoc(1, SDOC);
    let nxline = countparval("xline");
    if nxline == 0 { err!("must give xline= vector"); }
    let nyline = countparval("yline");
    if nyline == 0 { err!("must give yline= vector"); }
    if nxline != nyline { err!("lengths of xline, yline must be the same"); }
    let mut dcdp = 0f32;
    if getparfloat("dcdp", &mut dcdp) == 0 { err!("must give dcdp"); }
    // The original reads one point past each allocation on the final segment.
    let mut xline = vec![0f32; nxline as usize + 1];
    let mut yline = vec![0f32; nyline as usize + 1];
    getparfloat("xline", &mut xline); getparfloat("yline", &mut yline);
    for i in 1..nxline as usize {
        if xline[i] <= xline[i-1] { err!("xline values must increase monotonically"); }
    }
    let mut verbose = 0i32; getparint("verbose", &mut verbose);
    if verbose != 0 {
        warn!("%d points on the line", nxline);
        for i in 0..nxline as usize { warn!("point # %d x= %f y= %f", i as i32, xline[i], yline[i]); }
    }
    let mut cdpmin = 1001i32; getparint("cdpmin", &mut cdpmin);
    let mut distmax = dcdp; getparfloat("distmax", &mut distmax);
    checkpars();
    let mut length = 0f32;
    for i in 1..nxline as usize {
        length = (length as f64 + (((xline[i]-xline[i-1]) as f64).powf(2.0) + ((yline[i]-yline[i-1]) as f64).powf(2.0)).sqrt()) as f32;
    }
    let maxbin = (length / dcdp + 1.0) as i32;
    let mut xbin = vec![0f32; (maxbin.wrapping_add(1) as usize).max(2)];
    let mut ybin = vec![0f32; xbin.len()];
    let mut ibin = 0usize; xbin[0] = xline[0]; ybin[0] = yline[0];
    let mut iseg = 0usize;
    let mut sl = (yline[1]-yline[0])/(xline[1]-xline[0]);
    let mut ord = (yline[0]*xline[1]-yline[1]*xline[0])/(xline[1]-xline[0]);
    let mut dx = (dcdp as f64 * ((1.0f32 / (1.0f32 + sl*sl)) as f64).sqrt()) as f32;
    loop {
        if verbose != 0 { warn!("iseg=%d", iseg as i32); }
        loop {
            ibin += 1;
            if ibin >= xbin.len() { xbin.push(0.0); ybin.push(0.0); }
            xbin[ibin] = xbin[ibin-1]+dx; ybin[ibin] = sl*xbin[ibin]+ord;
            if verbose != 0 { warn!("ibin=%d x=%f y=%f ", ibin as i32, xbin[ibin], ybin[ibin]); }
            if !(xbin[ibin] <= xline[iseg+1]) { break; }
        }
        iseg += 1; ibin -= 1;
        sl = (yline[iseg+1]-yline[iseg])/(xline[iseg+1]-xline[iseg]);
        dx = (dcdp as f64 * ((1.0f32 / (1.0f32 + sl*sl)) as f64).sqrt()) as f32;
        ord = (yline[iseg]*xline[iseg+1]-yline[iseg+1]*xline[iseg])/(xline[iseg+1]-xline[iseg]);
        if verbose != 0 { warn!("sl=%f ord=%f", sl, ord); }
        let a = 1f32+sl*sl;
        let b = -2f32*xbin[ibin]+2f32*sl*(ord-ybin[ibin]);
        let c = ((xbin[ibin]*xbin[ibin]) as f64 + ((ord-ybin[ibin]) as f64).powf(2.0) - (dcdp*dcdp) as f64) as f32;
        let delta = b*b-4f32*a*c;
        let x2 = ((-(b as f64)+(delta as f64).sqrt())/2.0/a as f64) as f32;
        ibin += 1;
        if ibin >= xbin.len() { xbin.push(0.0); ybin.push(0.0); }
        xbin[ibin] = x2; ybin[ibin] = sl*xbin[ibin]+ord;
        if verbose != 0 {
            warn!("delta=%f x2=%f ibin=%d xbin=%f", delta, x2, ibin as i32, xbin[ibin]);
            warn!("ibin=%d x=%f y=%f ", ibin as i32, xbin[ibin], ybin[ibin]);
        }
        if !(xbin[ibin] <= xline[nxline as usize-1]) { break; }
    }
    let nbin = ibin;
    if verbose != 0 {
        warn!("length of the line :%f ", length);
        warn!("maximum number of bin:%d", maxbin);
        warn!("actual number of bin:%d", nbin as i32);
    }
    if gettr(&mut tr) == 0 { err!("can't read first trace"); }
    loop {
        let scale: f32;
        if tr.scalco < 0 { scale = (1.0/(tr.scalco as i32).abs() as f64) as f32; }
        else if tr.scalco > 0 { scale = tr.scalco as f32; }
        else { warn!("scalco = 0 ; 1 assumed"); scale = 1.0; }
        let xmp = (tr.gx.wrapping_add(tr.sx) as f64 / 2.0 * scale as f64) as f32;
        let ymp = (tr.gy.wrapping_add(tr.sy) as f64 / 2.0 * scale as f64) as f32;
        tr.cdp = cdpmin;
        let mut distmin = (((xmp-xbin[0]) as f64).powf(2.0)+((ymp-ybin[0]) as f64).powf(2.0)).sqrt() as f32;
        for i in 1..nbin {
            let dist = (((xmp-xbin[i]) as f64).powf(2.0)+((ymp-ybin[i]) as f64).powf(2.0)).sqrt() as f32;
            if dist < distmin { distmin = dist; tr.cdp = (i as i32).wrapping_add(cdpmin); }
        }
        if distmin > distmax { tr.cdp = 0; }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
