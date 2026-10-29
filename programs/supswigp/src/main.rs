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
// Original supswigp.c: Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.

//! SUPSWIGP - port of su/graphics/psplot/supswigp.c (Revision 1.37).
//! Credits: CWP: Dave Hale, Zhiming Li, Jack Cohen and John Stockwell;
//! Delphi: Alexander Koek (irregular spacing); Brian Zook, Southwest Research
//! Institute (scale factors and vsp style).

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;
extern "C" {
    fn access(path: *const i8, mode: i32) -> i32;
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUPSWIGP - PostScript Polygon-filled WIGgle plot of a segy data set\t",
    "\t\t\t\t\t\t\t\t\t",
    " supswigp <stdin [optional parameters] | ...\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t \t",
    " key=(keyword)\t\tif set, values of x2 are set from header field\t",
    "\t\t\tspecified by keyword\t\t\t\t",
    " n2=tr.ntr or number of traces in the data set\t(ntr is an alias for n2)",
    " d1=tr.d1 or tr.dt/10^6\tsampling interval in the fast dimension\t",
    "   =.004 for seismic \t\t(if not set)\t\t\t\t",
    "   =1.0 for nonseismic\t\t(if not set)\t\t\t\t",
    " d2=tr.d2\t\t\tsampling interval in the slow dimension\t",
    "   =1.0 \t\t\t(if not set)\t\t\t\t",
    " f1=tr.f1 or tr.delrt/10^3 or 0.0  first sample in the fast dimension\t",
    " f2=tr.f2 or tr.tracr or tr.tracl  first sample in the slow dimension\t",
    "   =1.0 for seismic\t\t    (if not set)\t\t\t",
    "   =d2 for nonseismic\t\t    (if not set)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " style=seismic\t\t normal (axis 1 horizontal, axis 2 vertical) or ",
    "\t\t\t vsp (same as normal with axis 2 reversed)\t",
    "\t\t\t Note: vsp requires use of a keyword\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0              =1 to print some useful information\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir=\t \tif non-empty, use the value as a directory path\t",
    "\t\t \tprefix for storing temporary files; else if the\t",
    "\t         \tthe CWP_TMPDIR environment variable is set use\t",
    "\t         \tits value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note that for seismic time domain data, the \"fast dimension\" is\t",
    " time and the \"slow dimension\" is usually trace number or range.\t",
    " Also note that \"foreign\" data tapes may have something unexpected\t",
    " in the d2,f2 fields, use segyclean to clear these if you can afford\t",
    " the processing time or use d2= f2= to override the header values if\t",
    " not.\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If key=keyword is set, then the values of x2 are taken from the header",
    " field represented by the keyword (for example key=offset, will show   ",
    " traces in true offset). This permit unequally spaced traces to be plotted.",
    " Type   sukeyword -o   to see the complete list of SU keywords.\t",
    "\t\t\t\t\t\t\t\t\t",
    " This program is really just a wrapper for the plotting program: pswigp",
    " See the pswigp selfdoc for the remaining parameters.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " On NeXT:\tsupswigp < infile [optional parameters]  | open\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Trace header fields accessed: ns, ntr, tracr, tracl, delrt, trid,     ",
    "\tdt, d1, d2, f1, f2, key specified by key\t\t\t",
];

fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let nt = tr.ns as usize;
    let mut ntr = tr.ntr; let mut have_ntr = ntr != 0;
    let mut verbose = 0; getparint("verbose", &mut verbose);
    let mut d1 = 0.0f32;
    if getparfloat("d1", &mut d1) == 0 {
        if tr.d1 != 0.0 { d1 = tr.d1; }
        else if tr.dt != 0 { d1 = (tr.dt as f64 / 1000000.0) as f32; }
        else if seismic { d1 = 0.004; warn!("tr.dt not set, assuming dt=0.004"); }
        else { d1 = 1.0; warn!("tr.d1 not set, assuming d1=1.0"); }
    }
    let mut f1 = 0.0f32;
    if getparfloat("f1", &mut f1) == 0 {
        if tr.f1 != 0.0 { f1 = tr.f1; }
        else if tr.delrt != 0 { f1 = (tr.delrt as f64 / 1000.0) as f32; }
    }
    let mut fill = 1; getparint("fill", &mut fill);
    if getparint("n2", &mut ntr) != 0 || getparint("ntr", &mut ntr) != 0 { have_ntr = true; }
    let mut d2 = 0.0f32;
    if getparfloat("d2", &mut d2) == 0 { d2 = if tr.d2 != 0.0 { tr.d2 } else { 1.0 }; }
    let mut f2 = 0.0f32;
    if getparfloat("f2", &mut f2) == 0 {
        f2 = if tr.f2 != 0.0 { tr.f2 } else if tr.tracr != 0 { tr.tracr as f32 } else if tr.tracl != 0 { tr.tracl as f32 } else if seismic { 1.0 } else { 0.0 };
    }
    let mut style = String::new(); if getparstring("style", &mut style) == 0 { style = "seismic".into(); }
    let mut key = String::new(); let mut irregular = getparstring("key", &mut key) != 0;
    let mut index1 = 0; let mut type1 = ""; let mut index2 = 0; let mut type2 = "";
    let depth = ["gelev", "selev", "sdepth", "gdel", "sdel", "swdep", "gwdep"].contains(&key.as_str());
    let coord = ["sx", "sy", "gx", "gy"].contains(&key.as_str());
    if irregular {
        type1 = hdtype(&key); index1 = getindex(&key);
        if index1 == -1 { err!("%s: keyword not in segy.h: '%s'", "supswigp.c", key); }
        if depth || coord { let k = if depth { "scalel" } else { "scalco" }; index2 = getindex(k); type2 = hdtype(k); }
    }
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 { tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default(); }
    if !tmpdir.is_empty() && unsafe { access(CString::new(tmpdir.as_str()).unwrap().as_ptr(), 2) } != 0 { err!("you can't write in %s (or it doesn't exist)", tmpdir); }
    let bindir = match std::env::var("CWPBIN") { Ok(s) => s, Err(_) => {
        let root = std::env::var("CWPROOT").unwrap_or_default();
        if root.is_empty() { warn!("CWPROOT environment variable is not set! "); err!("Set CWPROOT in shell environment as per instructions in CWP/SU Installation README files"); }
        format!("{}/bin", root)
    }};
    let mut datafp = None; let mut datafile = String::new(); let mut x2string = String::from("x2=");
    if !have_ntr || irregular {
        if verbose != 0 {
            if irregular { warn!("trace spacing from key %s", key); warn!("....  getting positions"); }
            else { warn!("n2 not getparred and ntr header field not set"); warn!("....  counting traces"); }
        }
        let fp = if tmpdir.is_empty() {
            if verbose != 0 { warn!("using tmpfile() call"); } etmpfile()
        } else {
            datafile = temporary_filename(&tmpdir); let fp = efopen(&datafile, "w+");
            if verbose != 0 { warn!("putting temporary files in %s", tmpdir); } fp
        };
        ntr = 0; let mut xmin = f32::MAX; let mut xmax = -f32::MAX;
        loop {
            if irregular {
                if ntr != 0 { x2string.push(','); }
                let mut val = Value::default(); gethval(&tr, index1, &mut val);
                let x = if depth || coord {
                    let mut scale = Value::default(); gethval(&tr, index2, &mut scale);
                    (vtod(type1, val) * 10.0f64.powf(vtod(type2, scale))) as f32
                } else { vtof(type1, val) };
                x2string.push_str(&sprintf!("%g", x));
                xmin = if xmin < x { xmin } else { x }; xmax = if xmax > x { xmax } else { x };
                if depth && style == "vsp" { for i in 0..nt { tr.data[i] *= -1.0; } fill = -1; }
            }
            ntr += 1; efwrite(&tr.data[..nt], &fp); if gettr(&mut tr) == 0 { break; }
        }
        if irregular {
            if depth && style == "vsp" { x2string.push_str(&sprintf!(" x2beg=%g x2end=%g", xmax, xmin)); style = "normal".into(); }
            if xmin == xmax { warn!("values in header %s all equal,", key); warn!("using f2=%f d2=%f", f2, d2); irregular = false; have_ntr = false; }
        } else if style == "vsp" { style = "normal".into(); warn!("style=vsp requires key= to be set"); }
        datafp = Some(fp);
    }
    if style == "vsp" { style = "normal".into(); }
    let mut plotcmd = if irregular {
        sprintf!("%s/pswigp n1=%d d1=%f f1=%f %s style=%s fill=%d", bindir, nt as i32, d1, f1, x2string, style, fill)
    } else {
        sprintf!("%s/pswigp n1=%d n2=%d d1=%f d2=%f f1=%f f2=%f style=%s fill=%d", bindir, nt as i32, ntr, d1, d2, f1, f2, style, fill)
    };
    for arg in argv.iter().skip(1) {
        if !["d1=", "d2=", "f1=", "f2=", "fill=", "style="].iter().any(|p| arg.starts_with(p)) {
            plotcmd.push_str(" \""); plotcmd.push_str(arg); plotcmd.push('"');
        }
    }
    let plotfp = epopen(&plotcmd, "w");
    if !have_ntr || irregular {
        let fp = datafp.as_ref().unwrap(); erewind(fp); let mut trbuf = vec![0.0f32; nt];
        for _ in 0..ntr { efread(&mut trbuf, fp); efwrite(&trbuf, &plotfp); }
    } else {
        loop { efwrite(&tr.data[..nt], &plotfp); if gettr(&mut tr) == 0 { break; } }
    }
    epclose(&plotfp);
    if let Some(fp) = datafp { efclose(&fp); if !datafile.is_empty() { eremove(&datafile); } }
    exit(CWP_Exit());
}
