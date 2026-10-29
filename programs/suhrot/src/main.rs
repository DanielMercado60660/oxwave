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

//! SUHROT - horizontal rotation of three-component data.
//! Port of su/main/multicomponent/suhrot.c (Revision 1.7).
//! Credits: Nils Maercklin, Geophysics, Kiel University, Germany, 1999.
use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUHROT - Horizontal ROTation of three-component data\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suhrot <stdin >stdout [optional parameters]\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " angle=rad\tunit of angles, choose \"rad\", \"deg\", or \"gon\"\t",
    " inv=0\t\t1 = inverse rotation (counter-clockwise)\t\t",
    " verbose=0\t1 = echo angle for each 3-C station\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " a=...\t\tarray of user-supplied rotation angles\t\t\t",
    " x=0.0,...\tarray of corresponding header value(s)\t\t\t",
    " key=tracf\theader word defining 3-C station (\"x\")\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " ... or input angles from files:\t\t\t\t\t",
    " n=0\t\t number of x and a values in input files\t\t",
    " xfile=...   file containing the x values as specified by the\t\t",
    " \t\t\t\t\"key\" parameter\t\t\t",
    " afile=...   file containing the a values\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Three adjacent traces are considered as one three-component\t\t",
    " dataset.\t\t\t\t\t\t\t\t",
    " By default, the data will be rotated from the Z-North-East (Z,N,E)\t",
    " coordinate system into Z-Radial-Transverse (Z,R,T).\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tIf one of the parameters \"a=\" or \"afile=\" is set, the data\t",
    "\tare rotated by these user-supplied angles. Specified x values\t",
    "\tmust be monotonically increasing or decreasing, and afile and\t",
    "\txfile are files of binary (C-style) floats.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut angle = String::new();
    let mut key = String::new();
    let mut inv = 0;
    let mut verbose = 0;
    if getparstring("angle", &mut angle) == 0 { angle = "rad".into(); }
    if getparstring("key", &mut key) == 0 { key = "tracf".into(); }
    getparint("inv", &mut inv);
    getparint("verbose", &mut verbose);
    let fangle = match angle.as_str() {
        "rad" => 1.0f32,
        "deg" => (180.0 / std::f64::consts::PI) as f32,
        "gon" => (200.0 / std::f64::consts::PI) as f32,
        _ => err!("unknown angle=%s", angle),
    };
    let mut x = Vec::<f32>::new();
    let mut a = Vec::<f32>::new();
    let nav = countparval("a");
    let user;
    if nav != 0 {
        let nxv = countparval("x");
        if nxv == 0 { x.push(0.0); } else { getparfloat("x", &mut x); }
        if nav as usize != x.len() { err!("number of values in \"x=\" and \"a=\" must be equal"); }
        getparfloat("a", &mut a);
        user = true;
    } else {
        let mut afile = String::new();
        if getparstring("afile", &mut afile) != 0 {
            let mut xfile = String::new();
            let mut n = 0;
            mustgetparstring("xfile", &mut xfile);
            mustgetparint("n", &mut n);
            if n < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
            if n == 0 { err!("n must be positive"); }
            let read_values = |path: &str, name: &str| -> Vec<f32> {
                let bytes = match std::fs::read(path) { Ok(b) => b, Err(_) => err!("cannot open %s=%s", name, path) };
                if bytes.len() / 4 < n as usize { err!("error reading %s=%s", name, path); }
                bytes.chunks_exact(4).take(n as usize).map(|b| f32::from_ne_bytes([b[0], b[1], b[2], b[3]])).collect()
            };
            a = read_values(&afile, "afile");
            x = read_values(&xfile, "xfile");
            user = true;
        } else { user = false; }
    }
    checkpars();
    if fangle != 1.0 { for value in &mut a { *value /= fangle; } }
    let typ = hdtype(&key);
    let index = getindex(&key);
    let mut headers = [tr.hdr; 3];
    let mut indata = vec![vec![0.0f32; nt]; 3];
    let mut outdata = vec![vec![0.0f32; nt]; 2];
    let mut icomp = 0usize;
    let mut nstat = 0;
    let mut phi = 0.0f32;
    loop {
        headers[icomp] = tr.hdr;
        indata[icomp].copy_from_slice(&tr.data[..nt]);
        if user && icomp == 0 {
            let mut val = Value::default();
            gethval(&tr, index, &mut val);
            let fval = vtof(typ, val);
            let mut result = [phi];
            intlin(&x, &a, a[0], a[a.len()-1], &[fval], &mut result);
            phi = result[0];
        }
        icomp += 1;
        if icomp == 3 {
            icomp = 0;
            nstat += 1;
            if !user {
                let dx = tr.gx as f32 - tr.sx as f32;
                let dy = tr.gy as f32 - tr.sy as f32;
                if dy != 0.0 {
                    phi = (dx as f64).atan2(dy as f64) as f32;
                    if phi < 0.0 { phi = (phi as f64 + 2.0 * std::f64::consts::PI) as f32; }
                } else {
                    phi = (if dx > 0.0 { 0.5 } else { 1.5 } * std::f64::consts::PI) as f32;
                }
                if dy == 0.0 || dx == 0.0 { phi = 0.0; }
            }
            if inv != 0 { phi = -phi; }
            if verbose != 0 { warn!("angle %g %s at station %d", phi * fangle, angle, nstat); }
            for i in 0..nt {
                outdata[0][i] = ((-(phi as f64)).cos() * indata[1][i] as f64 + (-(phi as f64)).sin() * indata[2][i] as f64) as f32;
                outdata[1][i] = (-( -(phi as f64)).sin() * indata[1][i] as f64 + (-(phi as f64)).cos() * indata[2][i] as f64) as f32;
            }
            for i in 0..3 {
                tr.hdr = headers[i];
                if i == 0 { tr.data[..nt].copy_from_slice(&indata[0]); tr.trid = 15; }
                else { tr.data[..nt].copy_from_slice(&outdata[i-1]); tr.trid = if i == 1 {17} else {16}; }
                puttr(&mut tr);
            }
        }
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
