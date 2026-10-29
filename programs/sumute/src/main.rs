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
//! SUMUTE - mute traces along a polygonal curve.
//! Port of su/main/windowing_sorting_muting/sumute.c (Revision 1.41).
//! Credits: SEP: Shuki Ronen; CWP: Jack K. Cohen, Dave Hale, John Stockwell;
//! DELPHI: Alexander Koek; USBG: Florian Bleibinhaus.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t   \t\t\t\t\t\t\t\t",
    " SUMUTE - MUTE above (or below) a user-defined polygonal curve with\t",
    "\t   the distance along the curve specified by key header word \t",
    " \t   \t\t\t\t\t\t\t\t",
    " sumute <stdin >stdout xmute= tmute= [optional parameters]\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " xmute=\t\tarray of position values as specified by\t",
    " \t\t\tthe `key' parameter\t\t\t\t",
    " tmute=\t\tarray of corresponding time values (sec)\t",
    " \t\t\tin case of air wave muting, correspond to \t",
    " \t\t\tair blast duration\t\t\t\t",
    "  ... or input via files:\t\t\t\t\t\t",
    " nmute=\t\tnumber of x,t values defining mute\t\t",
    " xfile=\t\tfile containing position values as specified by\t",
    " \t\t\tthe `key' parameter\t\t\t\t",
    " tfile=\t\tfile containing corresponding time values (sec)\t",
    "  ... or via header:\t\t\t\t\t\t\t",
    " hmute=\t\tkey header word specifying mute time\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " key=offset\t\tKey header word specifying trace offset \t",
    " \t\t\t\t=tracl  use trace number instead\t",
    " ntaper=0\t\tnumber of points to taper before hard\t\t",
    "\t\t\tmute (sine squared taper)\t\t\t",
    " mode=0\t   mute ABOVE the polygonal curve\t\t\t",
    "\t\t=1 to zero BELOW the polygonal curve\t\t\t",
    "\t\t=2 to mute below AND above a straight line. In this case",
    "\t\t \txmute,tmute describe the total time length of   ",
    "\t\t\tthe muted zone as a function of xmute the slope ",
    "\t\t\tof the line is given by the velocity linvel=\t",
    "\t \t=3 to mute below AND above a constant velocity hyperbola",
    "\t\t\tas in mode=2 xmute,tmute describe the total time",
    "\t\t\tlength of the mute zone as a function of xmute,  ",
    "\t\t\tthe velocity is given by the value of linvel=\t",
    " \t\t=4 to mute below AND above a user defined polygonal line",
    "\t\t\tgiven by xmute, tmute pairs. The widths in time ",
    "\t\t\tof the muted zone are given by the twindow vector",
    " linvel=330   \t\tconstant velocity for linear or hyperbolic mute\t",
    " tm0=0   \t\ttime shift of linear or hyperbolic mute at\t",
    "\t\t\t 'key'=0\t\t\t\t\t",
    " twindow=\tvector of mute zone widths in time, operative only in mode=4",
    "  ... or input via file:\t\t\t\t\t\t",
    " twfile= \t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " The tmute interpolant is extrapolated to the left by the smallest time",
    " sample on the trace and to the right by the last value given in the\t",
    " tmute array.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The files tfile and xfile are files of binary (C-style) floats.\t",
    "\t\t\t\t\t\t\t\t\t",
    " In the context of this program \"above\" means earlier time and\t",
    " \"below\" means later time (above and below as seen on a seismic section.",
    "\t\t\t\t\t\t\t\t\t",
    " The mode=2 option is intended for removing air waves. The mute is\t",
    " is over a narrow window above and below the line specified by the\t",
    " the line specified by the velocity \"linvel\". Here the values of     ",
    " tmute, xmute or tfile and xfile define the total time width of the mute.",
    "\t\t\t\t\t\t\t\t\t",
    " If data are spatial, such as the (z-x) output of a migration, then    ",
    " depth values are used in place of times in tmute and tfile. The value ",
    " of the depth sampling interval is given by the d1 header field\t",
    " You must use the option key=tracl in sumute in this case.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat: if data are seismic time sections, then tr.dt must be set. If ",
    " data are seismic depth sections, then tr.trid must be set to the value",
    " of TRID_DEPTH and the tr.d1 header field must be set.\t\t\t",
    " To find the value of TRID_DEPTH:  \t\t\t\t\t",
    " type: \t\t\t\t\t\t\t\t",
    "     sukeyword trid\t\t\t\t\t\t\t",
    "\tand look for the entry for \"Depth-Range (z-x) traces\"\t\t",
    "  \t\t\t\t\t\t\t\t\t",
];

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_linear_interp::intlin;
use std::io::Read;

fn read_values(name: &str, path: &str, n: usize) -> Vec<f32> {
    let mut file = match std::fs::File::open(path) { Ok(f) => f, Err(_) => err!("cannot open %s=%s\n", name, path) };
    let mut values = Vec::new();
    for _ in 0..n {
        let mut bytes = [0; 4];
        if file.read_exact(&mut bytes).is_err() { err!("error reading %s=%s\n", name, path); }
        values.push(f32::from_ne_bytes(bytes));
    }
    values
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut mode = 0;
    getparint("mode", &mut mode);
    let mut key = String::new();
    let hmute = getparstring("hmute", &mut key) != 0;
    let mut xmute = Vec::<f32>::new();
    let mut tmute = Vec::<f32>::new();
    let mut twindow = Vec::<f32>::new();
    let mut xfile = String::new();
    let mut tfile = String::new();
    if !hmute {
        if getparstring("xfile", &mut xfile) == 0 || getparstring("tfile", &mut tfile) == 0 {
            let nxmute = countparval("xmute");
            if nxmute == 0 { err!("must give xmute= vector"); }
            let ntmute = countparval("tmute");
            if ntmute == 0 { err!("must give tmute= vector"); }
            if nxmute != ntmute { err!("lengths of xmute, tmute must be the same"); }
            getparfloat("xmute", &mut xmute);
            getparfloat("tmute", &mut tmute);
            if mode == 4 {
                let ntwindow = countparval("twindow");
                if ntwindow == 0 { err!("must give twindow= vector"); }
                if nxmute != ntwindow { err!("lengths of xmute, twindow must be the same"); }
                getparfloat("twindow", &mut twindow);
            }
        } else {
            let mut n = 0;
            mustgetparint("nmute", &mut n);
            if n < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
            xmute = read_values("xfile", &xfile, n as usize);
            tmute = read_values("tfile", &tfile, n as usize);
            if mode == 4 { err!("cannot open twfile=%s\n", ""); }
        }
    }
    let mut ntaper = 0;
    getparint("ntaper", &mut ntaper);
    let mut below = 0;
    if getparint("below", &mut below) != 0 { mode = below; warn!("use of below parameter is obsolete. mode value set to %d \n", mode); }
    let mut absolute = 1;
    getparint("absolute", &mut absolute);
    if !hmute && getparstring("key", &mut key) == 0 { key = "offset".to_string(); }
    let mut linvel = 330.0f32;
    let mut tm0 = 0.0f32;
    getparfloat("linvel", &mut linvel);
    getparfloat("tm0", &mut tm0);
    checkpars();
    if linvel == 0.0 { err!("linear velocity can't be 0"); }
    let type_ = hdtype(&key);
    let index = getindex(&key);
    if ntaper < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut taper = vec![0.0f32; ntaper as usize];
    for k in 0..ntaper {
        let s = (((k + 1) as f64 * std::f64::consts::PI / (2 * ntaper) as f64).sin()) as f32;
        taper[k as usize] = s * s;
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't read first trace"); }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    if seismic { if tr.dt == 0 { err!("dt header field must be set"); } }
    else if tr.trid == 130 { if tr.d1 == 0.0 { err!("d1 header field must be set"); } }
    else { err!("tr.trid = %d, unsupported trace id", tr.trid as i32); }
    loop {
        let nt = tr.ns as i32;
        let mut tmin = (tr.delrt as f64 / 1000.0) as f32;
        let mut dt = (tr.dt as f64 / 1000000.0) as f32;
        if !seismic { tmin = 0.0; dt = tr.d1; }
        let mut val = Value::default();
        gethval(&tr, index, &mut val);
        let mut fval = vtof(type_, val);
        let mut t = 0.0f32;
        if hmute { t = (fval as f64 / 1000.0) as f32; }
        else {
            if tmute.is_empty() { err!("nmute must be positive"); }
            intlin(&xmute, &tmute, tmin, tmute[tmute.len()-1], &[fval], std::slice::from_mut(&mut t));
        }
        if absolute != 0 { fval = fval.abs(); }
        if mode == 0 {
            let nmute = nint(((t-tmin)/dt) as f64).min(nt);
            if nmute > 0 { tr.data[..nmute as usize].fill(0.0); }
            for i in 0..ntaper { if i + nmute > 0 { tr.data[(i+nmute) as usize] *= taper[i as usize]; } }
            tr.muts = if seismic { nint((t * 1000.0) as f64) as i16 } else { nint(t as f64) as i16 };
        } else if mode == 1 {
            let nmute = 0.max(nint(((tmin + nt as f32 * dt - t)/dt) as f64));
            if nmute > nt {
                // C's memset also clears the adjacent header before the samples.
                let preceding = (nmute - nt) as usize;
                if preceding > HDRBYTES / 4 { err!("mute extends before trace header"); }
                tr.hdr.as_bytes_mut()[HDRBYTES - preceding * 4..].fill(0);
                tr.data[..nt as usize].fill(0.0);
            } else {
                tr.data[(nt-nmute) as usize..nt as usize].fill(0.0);
            }
            for i in 0..ntaper { if nt > nmute+i && nmute+i > 0 { tr.data[(nt-nmute-1-i) as usize] *= taper[i as usize]; } }
            tr.mute = if seismic { nint((t * 1000.0) as f64) as i16 } else { nint(t as f64) as i16 };
        } else if (2..=4).contains(&mode) {
            let nmute;
            let ntair;
            if mode == 4 {
                // The original passes an uninitialized twfval to intlin.
                let mut tw = 0.0;
                intlin(&xmute, &twindow, twindow[0], twindow[twindow.len()-1], &[0.0], std::slice::from_mut(&mut tw));
                nmute = nint((tw/dt) as f64);
                ntair = nint((t/dt) as f64);
            } else {
                nmute = nint(((tmin+t)/dt) as f64);
                if mode == 2 { ntair = nint((tm0/dt+fval/linvel/dt) as f64); }
                else { let a = tm0/dt; let b = fval/linvel/dt; ntair = nint(((a*a+b*b) as f64).sqrt()); }
            }
            let topmute = 0.max(ntair-nmute/2).min(nt);
            let botmute = nt.min(ntair+nmute/2);
            if botmute >= topmute && topmute >= 0 { tr.data[topmute as usize..botmute as usize].fill(0.0); }
            for i in 0..ntaper { let itaper = ntair-nmute/2-i; if itaper > 0 && (itaper as usize) < tr.data.len() { tr.data[itaper as usize] *= taper[i as usize]; } }
            for i in 0..ntaper { let itaper = ntair+nmute/2+i; if itaper < nt && itaper >= 0 { tr.data[itaper as usize] *= taper[i as usize]; } }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
