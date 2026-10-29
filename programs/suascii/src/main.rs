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

//! SUASCII - print nonzero headers and data in various formats.
//! Port of su/main/data_conversion/suascii.c (Revision 1.19).
//! Credits: CWP: Jack K. Cohen; CENPET: Werner M. Heigl; RISSC: Nils Maercklin.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUASCII - print non zero header values and data in various formats    ",
    "                                                                       ",
    " suascii <stdin >ascii_file                                            ",
    "                                                                       ",
    " Optional parameter:                                                   ",
    "    bare=0     print headers and data                                  ",
    "        =1     print only data                                         ",
    "        =2     print headers only                                      ",
    "        =3     print data in print data in .csv format, e.g. for Excel ",
    "        =4     print data as tab delimited .txt file, e.g. for GnuPlot ",
    "        =5     print data as .xyz file, e.g. for plotting with GMT     ",
    "                                                                       ",
    "    ntr=50     maximum number of output traces (bare=3 or bare=4 only) ",
    "    index=0    don't include time/depth index in ascii file (bare=4)   ",
    "         =1    include time/depth index in ascii file                  ",
    "                                                                       ",
    "    key=       if set, name of keyword containing x-value              ",
    "               in .xyz output (bare=5 only)                            ",
    "    sep=       if set, string separating traces in .xyz output         ",
    "               (bare=5; default is no separation)                      ",
    "                                                                       ",
    "    verbose=0  =1 for detailed information                             ",
    "                                                                       ",
    " Notes:                                                                ",
    "    The programs suwind and suresamp provide trace selection and       ",
    "    subsampling, respectively.                                         ",
    "    With bare=0 and bare=1 traces are separated by a blank line.       ",
    "                                                                       ",
    "    With bare=3 a maximum of ntr traces are output in .csv format      ",
    "    (\"comma-separated value\"), e.g. for import into spreadsheet      ",
    "    applications like Excel.                                           ",
    "                                                                       ",
    "    With bare=4 a maximum of ntr traces are output in as tab delimited ",
    "    columns. Use bare=4 for plotting in GnuPlot.                       ",
    "                                                                       ",
    "    With bare=5 traces are written as \"x y z\" triples as required    ",
    "    by certain plotting programs such as the Generic Mapping Tools     ",
    "    (GMT). If sep= is set, traces are separated by a line containing   ",
    "    the string provided, e.g. sep=\">\" for GMT multisegment files.    ",
    "                                                                       ",
    "    \"option=\" is an acceptable alias for \"bare=\".                  ",
    "                                                                       ",
    " Related programs: sugethw, sudumptrace                                ",
    "                                                                       ",
];

fn seismic(id: i16) -> bool { (0..=17).contains(&id) || id == 109 }

fn print_xyz(tr: &Trace, key: Option<&str>, itr: &mut i32, d2: &mut f32, f2: &mut f32) {
    let x: f32;
    if let Some(key) = key {
        let typ = hdtype(key);
        let index = getindex(key);
        let mut val = Value::default();
        gethval(tr, index, &mut val);
        x = vtof(typ, val);
    } else {
        if *itr == 0 {
            *d2 = if tr.d2 != 0.0 { tr.d2 } else { 1.0 };
            *f2 = tr.f2;
        }
        x = (*itr as f32) * *d2 + *f2;
        *itr += 1;
    }
    let (d1, f1) = if seismic(tr.trid) && tr.dt != 0 {
        ((tr.dt as f64 / 1000000.0) as f32, (tr.delrt as f64 / 1000.0) as f32)
    } else {
        (if tr.d1 != 0.0 { tr.d1 } else { 1.0 }, tr.f1)
    };
    for i in 0..tr.ns as usize {
        printf!("%g %g %g\n", x, (i as f32) * d1 + f1, tr.data[i]);
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut bare = 0;
    let mut index = 0;
    let mut ntr = 50;
    let mut key = String::new();
    let mut sep = String::new();
    let mut verbose = 0;
    if getparint("bare", &mut bare) == 0 && getparint("option", &mut bare) == 0 { bare = 0; }
    getparint("index", &mut index);
    getparint("ntr", &mut ntr);
    let has_key = getparstring("key", &mut key) != 0;
    let has_sep = getparstring("sep", &mut sep) != 0;
    getparint("verbose", &mut verbose);
    if bare < 0 || bare > 5 { err!("bare=%d must be in range 0...5", bare); }
    if ntr < 1 { err!("ntr=%d must be positive", ntr); }
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let mut dt: f32;
    let tmin: f32;
    if !seismic(tr.trid) {
        if verbose != 0 { warn!("input is not seismic data, trid=%d", tr.trid); }
        dt = tr.d1;
        if dt == 0.0 { getparfloat("dt", &mut dt); }
        if dt == 0.0 { err!("d1 field is zero and not getparred"); }
        tmin = tr.f1;
        if tmin == 0.0 && verbose != 0 { warn!("f1 field is zero or not set"); }
    } else {
        if verbose != 0 { warn!("input is seismic data, trid=%d", tr.trid); }
        dt = (tr.dt as f64 / 1000000.0) as f32;
        if dt == 0.0 { getparfloat("dt", &mut dt); }
        if dt == 0.0 { err!("dt field is zero and not getparred"); }
        tmin = (tr.delrt as f64 / 1000.0) as f32;
        if tmin == 0.0 && verbose != 0 { warn!("delrt field is zero or not set"); }
    }
    checkpars();
    let nt = tr.ns as usize;
    if verbose != 0 { warn!("dt=%g  tmin=%g  nt=%d", dt, tmin, nt as i32); }
    let mut data = if bare == 3 || bare == 4 { vec![vec![0.0f32; nt]; ntr as usize] } else { Vec::new() };
    let mut j = 0usize;
    let mut itr = 0;
    let mut d2 = 0.0;
    let mut f2 = 0.0;
    loop {
        match bare {
            0 => {
                printheader(&tr);
                for i in 0..nt { printf!("%5d ", (i + 1) as i32); printf!("%11.4e\n", tr.data[i]); }
                putchar(b'\n' as i32);
            }
            1 => {
                for i in 0..nt { printf!("%15.8f\n", tr.data[i]); }
                putchar(b'\n' as i32);
            }
            2 => printheader(&tr),
            3 | 4 => data[j].copy_from_slice(&tr.data[..nt]),
            5 => {
                print_xyz(&tr, if has_key { Some(&key) } else { None }, &mut itr, &mut d2, &mut f2);
                if has_sep { printf!("%s\n", sep); }
            }
            _ => {}
        }
        j += 1;
        if (bare == 3 || bare == 4) && j == ntr as usize { break; }
        if gettr(&mut tr) == 0 { break; }
    }
    if verbose != 0 {
        if bare == 3 || bare == 4 { warn!("bare=%d, j=%d traces read, maximum number is ntr=%d", bare, j as i32, ntr); }
        else { warn!("bare=%d, j=%d traces read", bare, j as i32); }
    }
    if bare == 3 {
        for i in 0..tr.ns as usize {
            for k in 0..j-1 { printf!("%g,", data[k][i]); }
            printf!("%g\n", data[j-1][i]);
        }
    }
    if bare == 4 {
        for k in 0..nt {
            if index != 0 { printf!("%15.8f\t", tmin + (k as f32) * dt); }
            if j == 1 {
                for i in 0..j { printf!("%15.8f\n", data[i][k]); }
            } else {
                for i in 0..j-1 { printf!("%15.8f\t", data[i][k]); }
                printf!("%15.8f\n", data[j-1][k]);
            }
        }
    }
    exit(CWP_Exit());
}
