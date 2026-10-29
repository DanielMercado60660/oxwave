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
//! SUSHW - Set header words from trace numbers or a file.
//! Port of su/main/headers/sushw.c (Revision 1.25).
//! Credits: SEP: Einar Kajartansson; CWP: Jack K. Cohen;
//! CWP: John Stockwell, added multiple fields and infile options.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUSHW - Set one or more Header Words using trace number, mod and\t",
    "\t integer divide to compute the header word values or input\t",
    "\t the header word values from a file\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " ... compute header fields\t\t\t\t\t\t",
    "   sushw <stdin >stdout key=cdp,.. a=0,..  b=0,.. c=0,.. d=0,.. j=..,..",
    " \t\t\t\t\t\t\t\t\t",
    " ... or read headers from a binary file\t\t\t\t",
    "   sushw <stdin > stdout  key=key1,..    infile=binary_file\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters for setting headers from infile:\t\t\t",
    " key=key1,key2 ... is the list of header fields as they appear in infile",
    " infile= \tbinary file of values for field specified by\t\t",
    " \t\tkey1,key2,...\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters ():\t\t\t\t\t\t",
    " key=cdp,...\t\t\theader key word(s) to set \t\t",
    " a=0,...\t\t\tvalue(s) on first trace\t\t\t",
    " b=0,...\t\t\tincrement(s) within group\t\t",
    " c=0,...\t\t\tgroup increment(s)\t \t\t",
    " d=0,...\t\t\ttrace number shift(s)\t\t\t",
    " j=ULONG_MAX,ULONG_MAX,...\tnumber of elements in group\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Fields that are getparred must have the same number of entries as key\t",
    " words being set. Any field that is not getparred is set to the default",
    " value(s) above. Explicitly setting j=0 will set j to ULONG_MAX.\t",
    " \t\t\t\t\t\t\t\t\t",
    " The value of each header word key is computed using the formula:\t",
    " \ti = itr + d\t\t\t\t\t\t\t",
    " \tval(key) = a + b * (i % j) + c * (int(i / j))\t\t\t",
    " where itr is the trace number (first trace has itr=0, NOT 1)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " 1. set every dt field to 4ms\t\t\t\t\t\t",
    " \tsushw <indata key=dt a=4000 |...\t\t\t\t",
    " 2. set the sx field of the first 32 traces to 6400, the second 32 traces",
    "    to 6300, decrementing by -100 for each 32 trace groups\t\t",
    "   ...| sushw key=sx a=6400 c=-100 j=32 |...\t\t\t\t",
    " 3. set the offset fields of each group of 32 traces to 200,400,...,6400",
    "   ...| sushw key=offset a=200 b=200 j=32 |...\t\t\t\t",
    " 4. perform operations 1., 2., and 3. in one call\t\t\t",
    "  ..| sushw key=dt,sx,offset a=4000,6400,200 b=0,0,200 c=0,-100,0 j=0,32,32 |",
    " \t\t\t\t\t\t\t\t\t",
    " In this example, we set every dt field to 4ms.  Then we set the first\t",
    " 32 shotpoint fields to 6400, the second 32 shotpoint fields to 6300 and",
    " so forth.  Next we set each group of 32 offset fields to 200, 400, ...,",
    " 6400.\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Example of a typical processing sequence using suchw:\t\t\t",
    "  sushw <indata key=dt a=4000 |\t\t\t\t\t",
    "  sushw key=sx a=6400 c=-100 j=32 |\t\t\t\t\t",
    "  sushw key=offset a=200 b=200 j=32 |\t\t\t     \t\t",
    "  suchw key1=gx key2=offset key3=sx b=1 c=1 |\t\t     \t\t",
    "  suchw key1=cdp key2=gx key3=sx b=1 c=1 d=2 >outdata\t     \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Again, it is possible to eliminate the multiple calls to both sushw and",
    " sushw, as in Example 4.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Reading header values from a binary file:\t\t\t\t",
    " If the parameter infile=binary_file is set, then the values that are to",
    " be set for the fields specified by key=key1,key2,... are read from that",
    " file. The values are read sequentially from the file and assigned trace",
    " by trace to the input SU data. The infile consists of C (unformated)\t",
    " binary floats in the form of an array of size (nkeys)*(ntraces) where\t",
    " nkeys is the number of floats in the first (fast) dimension and ntraces",
    " is the number of traces.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Comment: \t\t\t\t\t\t\t\t",
    " Users wishing to edit one or more header fields (as in geometry setting)",
    " may do this via the following sequence:\t\t\t\t",
    "     sugethw < sudata output=geom  key=key1,key2 ... > hdrfile \t",
    " Now edit the ASCII file hdrfile with any editor, setting the fields\t",
    " appropriately. Convert hdrfile to a binary format via:\t\t",
    "     a2b < hdrfile n1=nfields > binary_file\t\t\t\t",
    " Then set the header fields via:\t\t\t\t\t",
    "     sushw < sudata infile=binary_file key=key1,key2,... > sudata.edited",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat: \t\t\t\t\t\t\t\t",
    " If the (number of traces)*(number of key words) exceeds the number of\t",
    " values in the infile then the user may still set a single header field",
    " on the remaining traces via the parameters key=keyword a,b,c,d,j.\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    "    sushw < sudata=key1,key2 ... infile=binary_file [Optional Parameters]",
];

fn setval(type_: &str, val: &mut Value, a: f64, b: f64, c: f64, i: f64, j: f64) {
    let m = if j == 0.0 { i } else { i - j * (i / j).floor() };
    let x = a + b * m + c * ((i / j) as i32 as f64);
    match type_.as_bytes()[0] {
        b's' => err!("can't set char header word"),
        b'h' => val.set_h(x as i32 as i16),
        b'u' => val.set_u(x as i32 as u16),
        b'l' => val.set_l(x as i64),
        b'v' => val.set_v(x as u64),
        b'i' => val.set_i(x as i32),
        b'p' => val.set_p(x as i64 as u32),
        b'f' => val.set_f(x as f32),
        b'd' => { val.set_d(x); err!("unknown type %s", type_); }
        _ => err!("unknown type %s", type_),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let nkeys = countparval("key") as usize;
    let mut key = vec![String::new(); SU_NKEYS];
    if nkeys != 0 { getparstringarray("key", &mut key); } else { key[0] = "cdp".into(); }
    let mut types = Vec::new();
    let mut indices = Vec::new();
    for k in &key[..nkeys] {
        types.push(hdtype(k));
        indices.push(getindex(k));
    }
    let mut infile = String::new();
    getparstring("infile", &mut infile);
    let from_file = !infile.is_empty();
    let infp = if from_file { Some(efopen(&infile, "r")) } else { None };
    let mut params = vec![vec![0.0f64; nkeys]; 5];
    let names = ["a", "b", "c", "d", "j"];
    if !from_file {
        for p in 0..5 {
            let n = countparval(names[p]);
            if n != 0 {
                if n as usize != nkeys { err!("number of %s values not equal to number of keys", names[p]); }
                getpardouble(names[p], &mut params[p]);
                if p == 4 { for v in &mut params[p] { if *v == 0.0 { *v = u64::MAX as f64; } } }
            } else if p == 4 {
                for v in &mut params[p] { *v = u64::MAX as f64; }
            }
        }
    }
    let mut afile = vec![0.0f32; nkeys];
    let mut count: i32 = 0;
    let mut itr: i32 = 0;
    let mut val = Value::default();
    checkpars();
    while gettr(&mut tr) != 0 {
        if let Some(ref fp) = infp {
            if efread(&mut afile, fp) != 0 {
                for k in 0..nkeys {
                    setval(types[k], &mut val, afile[k] as f64, 0.0, 0.0, 0.0, u64::MAX as f64);
                    puthval(&mut tr, indices[k], &val);
                    count = count.wrapping_add(1);
                }
            }
        } else {
            for k in 0..nkeys {
                let i = itr as f64 + params[3][k];
                setval(types[k], &mut val, params[0][k], params[1][k], params[2][k], i, params[4][k]);
                puthval(&mut tr, indices[k], &val);
            }
        }
        itr = itr.wrapping_add(1);
        puttr(&mut tr);
    }
    if let Some(ref fp) = infp {
        efclose(fp);
        if count < itr.wrapping_mul(nkeys as i32) {
            warn!("itr=%d > count=%d %s", itr.wrapping_mul(count), count);
            warn!("n traces=%d > data count =%d", itr.wrapping_mul(nkeys as i32), count);
        }
    }
    exit(CWP_Exit());
}
