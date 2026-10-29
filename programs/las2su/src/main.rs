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
//! LAS2SU - LAS well log curves to SU traces.
//! Port of su/main/well_logs/las2su.c (Revision 1.9).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: Chris Liner, based on code by Ferhan Ahmed and a2b.c; Saudi Aramco
//! granted permission to release code developed at EXPEC-ARC. CWP: John Stockwell;
//! CENPET lasstrip: Werner Heigl; Rob Hardy.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " LAS2SU - convert las2 format well log curves to su traces\t",
    "\t\t\t\t\t\t\t\t",
    "  las2su <stdin nskip= ncurve= >stdout [optional params]\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " ncurve=automatic\tnumber of log curves in the las file\t",
    " dz=0.5\t\tinput depth sampling (ft)\t\t",
    " m=0\t\t\toutput (d1,f1) in feet\t\t\t",
    "\t\t\t=1 output (d1,f1) in meters\t\t",
    " ss=0\t\t\tdo not subsample (unless nz > 32767 )\t",
    "\t\t\t=1 pass every other sample\t\t",
    " verbose=0\t\t=1 to echo las header lines to screen\t",
    " outhdr=las_header.asc\tname of file for las headers\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " 1. It is recommended to run LAS_CERTIFY available from CWLS\t",
    "    at http://cwls.org.\t\t\t\t\t",
    " 2. First log curve MUST BE depth.\t\t\t\t",
    " 3. If number of depth levels > 32767 (segy NT limit)\t\t",
    "    then input is subsampled by factor of 2 or 4 as needed\t",
    " 4. Logs may be isolated using tracl header word (1,2,...,ncurve) ",
    "    tracl=1 is depth\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " If the input LAS file contains sonic data as delta t or interval",
    " transit time and you plan to use these data for generating a ",
    " reflection coefficient time series in suwellrf, convert the sonic",
    " trace to velocity using suop with op=s2v (sonic to velocity) ",
    " before input of the velocity trace to suwellrf.\t\t",
    "\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t",
    " No trap is made for the commonly used null value in LAS files ",
    " (-999.25). The null value will be output as ?999.25, which\t",
    " really messes up a suxwigb display of density data because the",
    " ?999.25 skews the more or less 2.5 values of density.\t\t",
    " The user needs to edit out null values (-999.25) before running",
    " other programs, such as \"suwellrf\".\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut verbose, mut m, mut ss, mut ncurve) = (0, 0, 0, 0);
    let mut dz = 0.5f32;
    let mut outhdr = String::new();
    getparint("verbose", &mut verbose);
    getparint("m", &mut m);
    getparint("ss", &mut ss);
    getparfloat("dz", &mut dz);
    if getparstring("outhdr", &mut outhdr) == 0 { outhdr = "las_header.asc".into(); }
    let outhdrfp = efopen(&outhdr, "w");
    let is_ncurve_set = getparint("ncurve", &mut ncurve) == 0;
    checkpars();
    let mut discard = 1;
    let mut in_tilde_c_block = false;
    let input = stdin();
    loop {
        let mut line = Vec::new();
        while line.len() < 999 {
            let c = fgetc(&input);
            if c < 0 { break; }
            line.push(c as u8);
            if c == 10 { break; }
        }
        for &c in &line { fputc(c as i32, &outhdrfp); }
        fputc(10, &outhdrfp);
        if line.starts_with(b"~C") { in_tilde_c_block = true; }
        else if line.starts_with(b"~P") { in_tilde_c_block = false; }
        if line.starts_with(b"#") && in_tilde_c_block { discard += 1; }
        if in_tilde_c_block && is_ncurve_set { ncurve += 1; }
        if line.starts_with(b"~A") { break; }
    }
    if is_ncurve_set { ncurve -= discard; }
    efclose(&outhdrfp);
    let nzmax = SU_NFLTS as i32;
    if verbose != 0 { warn!("ncurve %d nzmax %d", ncurve, nzmax); }
    // C converts a negative ncurve to size_t and fails in ealloc2float.
    // Report that failure before converting it to a Rust allocation size.
    if ncurve < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    // With no curves C dereferences x[0] outside its allocation.
    if ncurve == 0 { err!("las2su.c: ncurve must be positive"); }
    let mut x = vec![vec![0.0f32; nzmax as usize]; ncurve as usize];
    let mut nz = 0usize;
    while let Some(string) = fgets(400, &input) {
        // C writes past the allocated columns before reaching subsampling.
        // Keep storage bounded and fail cleanly on this undefined input.
        if nz == nzmax as usize {
            err!("las2su.c: number of depth levels exceeds %d", nzmax);
        }
        let mut tokens = string.split(|c| c == ' ' || c == '\n' || c == '\t').filter(|s| !s.is_empty());
        for curve in 0..ncurve as usize {
            let token = tokens.next().unwrap_or("");
            x[curve][nz] = atof(token) as f32;
        }
        nz += 1;
    }
    warn!("nz=%i", nz as i32);
    for _ in 0..2 {
        if nz > nzmax as usize || ss == 1 {
            nz /= 2;
            dz = (2.0 * dz as f64) as f32;
            for curve in 0..ncurve as usize {
                for iz in 0..nz { x[curve][iz] = x[curve][2 * iz]; }
            }
            if verbose != 0 { warn!("New: nz=%i dz=%g ft\n", nz as i32, dz); }
        }
    }
    let mut tr = Trace::new();
    tr.trid = 1;
    tr.ns = nz as u16;
    tr.dt = (1000.0f32 * dz) as u32 as u16;
    if m == 0 { tr.d1 = dz; tr.f1 = x[0][0]; }
    else { tr.d1 = (dz as f64 / 3.28084) as f32; tr.f1 = (x[0][0] as f64 / 3.28084) as f32; }
    for curve in 0..ncurve as usize {
        tr.tracl = curve as i32 + 1;
        for iz in 0..nz { tr.data[iz] = x[curve][iz]; }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
