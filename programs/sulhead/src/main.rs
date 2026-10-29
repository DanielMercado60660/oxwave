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
//! SULHEAD - load column data into trace headers. Port of su/main/headers/sulhead.c (Revision 1.6).
//! Credits: Balasz Nemeth, Potash Corporation, Saskatoon Saskatchewan; given to CWP in 2008.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "SULHEAD - Load information from an ascii column file into HEADERS based",
    "\t   on the value of the user specified header field\t\t",
    "  sulhead < inflie > outfile cf=Column_file key=..  [ optional parameters]",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " cf=Name of column file\t\t\t\t\t\t",
    " key=key1,key2,...Number of column entires\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " mc=1\t\tColumn number to use to match rows to traces\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Notes:\t\t\t\t\t\t\t\t\t",
    " Caveat: This is not simple trace header setting, but conditional\t",
    " setting.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " This utility reads the column file and loads the values into the\t",
    " specified header locations. Each column represents one set of header  ",
    " words, one of them (#mc) is used to match the rows to the traces\t",
    " using header tr.key[mc].\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " key=cdp,ep,sx   mc=1\tcf=file\t\t\t\t\t\t",
    " file contains:\t\t\t\t\t\t\t",
    "\t1  2  3\t\t\t\t\t\t\t\t",
    "\t2  3  4\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " if tr.cdp = 1 then tr.ep and tr.sx will be set to 2 and 3\t\t",
    " if tr.cdp = 2 then tr.ep and tr.sx will be set to 3 and 4\t\t",
    " if tr.cdp=other than tr.trid=3\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Caveat: the user has to make it sure that number of entires in key=\t",
    "\t is equal the number of columns stored in the file.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For simple mass setting of header words, see selfdoc of:  sushw\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn setval_f_i(type_: &str, a: i32, val: &mut Value) {
    match type_.as_bytes()[0] {
        b's' => err!("can't set char header word"),
        b'h' => val.set_h(a as i16),
        b'u' => val.set_u(a as u16),
        b'l' => val.set_l(a as i64),
        b'v' => val.set_v(a as u64),
        b'i' => val.set_i(a),
        b'p' => val.set_p(a as u32),
        b'f' => val.set_f(a as f32),
        b'd' => { val.set_d(a as f64); err!("unknown type %s", type_); }
        _ => err!("unknown type %s", type_),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut cf = String::new();
    mustgetparstring("cf", &mut cf);
    let nc = countparval("key");
    let mut key = vec![String::new(); SU_NKEYS];
    if nc != 0 { getparstringarray("key", &mut key); } else { key[0] = "cdp".to_string(); }
    let mut mc = 1;
    let mut verbose = 0;
    getparint("mc", &mut mc);
    getparint("verbose", &mut verbose);
    checkpars();
    // The C check uses && and leaves these cases undefined; reject them before indexing.
    if mc < 1 || mc > nc { err!("Master column index has to be in the range [1..nc].\n"); }
    if nc as usize > SU_NKEYS { err!("Too many header keys (maximum %d).", SU_NKEYS as i32); }
    let mut kind = vec![0; SU_NKEYS];
    let mut ktype = vec![""; SU_NKEYS];
    for ikey in 0..nc as usize {
        ktype[ikey] = hdtype(&key[ikey]);
        kind[ikey] = getindex(&key[ikey]);
    }
    let tmkind = kind[(mc - 1) as usize];
    let tmtype = ktype[(mc - 1) as usize];
    let fp = match fopen(&cf, "r") { Some(fp) => fp, None => err!("No file %s was found!", cf) };
    let mut rc = 0;
    while !feof(&fp) { if fgetc(&fp) == 0x0a { rc += 1; } }
    if verbose != 0 {
        warn!("Number of rows in file: %d", rc);
        warn!("Number of columns in file: %d", nc);
    }
    erewind(&fp);
    let mut tab = vec![vec![0i32; nc as usize]; rc as usize];
    let mut tmp = 0.0f32;
    for ir in 0..rc as usize {
        for ic in 0..nc as usize {
            if let Some(x) = fscanf_f32(&fp) { tmp = x; }
            tab[ir][ic] = (tmp as f64 + 0.5) as i32;
        }
    }
    fclose(&fp);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let mut val = Value::default();
    loop {
        let mut tmval = Value::default();
        gethval(&tr, tmkind, &mut tmval);
        let mut flag = 0;
        for ir in 0..rc as usize {
            if vtoi(tmtype, tmval) == tab[ir][(mc - 1) as usize] {
                for ic in 0..nc as usize {
                    setval_f_i(ktype[ic], tab[ir][ic], &mut val);
                    puthval(&mut tr, kind[ic], &val);
                }
                flag = 1;
            }
        }
        if flag == 0 {
            tr.trid = 3;
            if verbose != 0 {
                warn!("No data for trace with: ");
                warn!(" %s ", key[(mc - 1) as usize]);
                warn!("=");
                fprintfval(&stderr(), tmtype, tmval);
                warn!("\n");
            }
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
