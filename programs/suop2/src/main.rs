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

//! SUOP2 - binary operations on two data sets.
//! Port of su/main/operations/suop2.c (Revision 1.25).
//! Credits: SEP: Shuki Ronen; CWP: Jack K. Cohen; CWP: John Stockwell;
//! Fernando M. Roxo da Motta.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUOP2 - do a binary operation on two data sets\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suop2 data1 data2 op=diff [trid=111] >stdout\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t\t",
    " \top=diff\t\tdifference of two panels of seismic data\t",
    " \t\t\t=sum  sum of two panels of seismic data\t\t",
    " \t\t\t=prod product of two panels of seismic data\t",
    " \t\t\t=quo quotient of two panels of seismic data\t",
    " \t\t\t=ptdiff differences of a panel and single trace\t",
    " \t\t\t=ptsum sum of a panel and single trace\t\t",
    " \t\t\t=ptprod product of a panel and single trace\t",
    " \t\t\t=ptquo quotient of a panel and single trace\t",
    " \t\t\t=zipper do \"zipper\" merge of two panels\t",
    "\t\t\t=zippol convert polar to rectangular and then zip",
    " \t\t\t\t\t\t\t\t\t",
    "  trid=FUNPACKNYQ\toutput trace identification code. (This option  ",
    " \t\t\tis active only for op=zipper)\t\t\t",
    "\t\t\tFor SU version 39-43 FUNPACNYQ=111\t\t",
    " \t\t\t(See: sukeyword trid     for current value)\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note1: Output = data1 \"op\" data2 with the header of data1\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note2: For convenience and backward compatibility, this\t\t",
    " \tprogram may be called without an op code as:\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For:  panel \"op\" panel  operations: \t\t\t\t",
    " \tsusum  file1 file2 == suop2 file1 file2 op=sum\t\t\t",
    " \tsudiff file1 file2 == suop2 file1 file2 op=diff\t\t\t",
    " \tsuprod file1 file2 == suop2 file1 file2 op=prod\t\t\t",
    " \tsuquo  file1 file2 == suop2 file1 file2 op=quo\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For:  panel \"op\" trace  operations: \t\t\t\t",
    " \tsuptsum  file1 file2 == suop2 file1 file2 op=ptsum\t\t",
    " \tsuptdiff file1 file2 == suop2 file1 file2 op=ptdiff\t\t",
    " \tsuptprod file1 file2 == suop2 file1 file2 op=ptprod\t\t",
    " \tsuptquo  file1 file2 == suop2 file1 file2 op=ptquo\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note3: If an explicit op code is used it must FOLLOW the\t\t",
    "\tfilenames.\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note4: With op=quo and op=ptquo, divide by 0 is trapped and 0 is returned.",
    " \t\t\t\t\t\t\t\t\t",
    " Note5: Weighted operations can be specified by setting weighting\t",
    "\tcoefficients for the two datasets:\t\t\t\t",
    "\tw1=1.0\t\t\t\t\t\t\t\t",
    "\tw2=1.0\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note6: With op=zipper, it is possible to set the output trace id code ",
    " \t\t(See: sukeyword trid)\t\t\t\t\t",
    "  This option processes the traces from two files interleaving its samples.",
    "  Both files must have the same trace length and must not longer than\t",
    "  SU_NFLTS/2  (as in SU 39-42  SU_NFLTS=32768).\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    "  Being \"tr1\" a trace of data1 and \"tr2\" the corresponding trace of",
    "  data2, The merged trace will be :\t\t\t\t\t",
    "\t\t\t  \t\t\t\t\t\t",
    "  tr[2*i]= tr1[i]\t\t\t\t\t\t\t",
    "  tr[2*i+1] = tr2[i]\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "  The default value of output tr.trid is that used by sufft and suifft,",
    "  which is the trace id reserved for the complex traces obtained through",
    "  the application of sufft. See also, suamp.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note 7: op=zippol is like op=zipper, but the input samples are polar\t",
    "\t(amplitude and phase) and are converted to cartesian (real, imag)",
    "\tbefore interleaving them.\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For operations on non-SU binary files  use:farith \t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(2, SDOC);
    let fp1 = efopen(&argv[1], "r");
    let fp2 = efopen(&argv[2], "r");
    let mut tr1 = Trace::new();
    let mut tr2 = Trace::new();
    let mut w1 = 0.0f32;
    let mut w2 = 0.0f32;
    getparfloat("w1", &mut w1);
    getparfloat("w2", &mut w2);
    let mut op = String::from("diff");
    getparstring("op", &mut op);
    let iop = match op.as_str() {
        "sum" => 1, "diff" => 2, "prod" => 3, "quo" => 4,
        "zipper" => 5, "zippol" => 6, "ptsum" => 7,
        "ptdiff" => 8, "ptprod" => 9, "ptquo" => 10,
        _ => err!("unknown operation=\"%s\", see self-doc", op),
    };
    let mut trid = 0i32;
    if iop == 5 || iop == 6 {
        if getparint("trid", &mut trid) == 0 { trid = 111; }
    }
    checkpars();
    let mut nt = 0usize;
    let mut itr = 0i32;
    if iop <= 7 {
        let mut nbytes = 0;
        while fgettr(&fp1, &mut tr1) != 0 {
            nbytes = fgettr(&fp2, &mut tr2);
            if nbytes == 0 { break; }
            nt = tr1.ns as usize;
            if nt != tr2.ns as usize {
                warn!("trace %d:", itr);
                err!("%s and %s have different ns (%d vs %d)", argv[1], argv[2], tr1.ns, tr2.ns);
            }
            let complex = tr1.trid == 111 && tr2.trid == 111;
            if complex { warn!("performing complex arithmetic"); }
            if w1 != 0.0 && w2 != 0.0 {
                for i in 0..nt { tr1.data[i] *= w1; tr2.data[i] *= w2; }
            }
            // C writes past the trace array for oversized interleaved output.
            // Reject that undefined case before either expansion loop writes.
            if (iop == 5 || iop == 6) && nt > SU_NFLTS / 2 {
                err!("op=%s requires ns <= %d (got %d)", op, SU_NFLTS / 2, nt);
            }
            match iop {
                1 => for i in 0..nt { tr1.data[i] += tr2.data[i]; },
                2 => for i in 0..nt { tr1.data[i] -= tr2.data[i]; },
                3 if complex => { let mut i = 2; while i + 2 < nt { complex_mul(&mut tr1, &tr2, i); i += 2; } },
                3 => for i in 0..nt { tr1.data[i] *= tr2.data[i]; },
                4 if complex => { let mut i = 2; while i < nt { complex_div(&mut tr1, &tr2, i); i += 2; } },
                4 => for i in 0..nt { let denom = tr2.data[i]; if denom == 0.0 { tr1.data[i] = 0.0; } else { tr1.data[i] /= denom; } },
                5 => { for i in (0..nt).rev() { tr1.data[2*i] = tr1.data[i]; tr1.data[2*i+1] = tr2.data[i]; } tr1.trid = trid as i16; tr1.ns = (2*nt) as u16; },
                6 => { for i in (0..nt).rev() { tr1.data[2*i+1] = (tr1.data[i] as f64 * (tr2.data[i] as f64).sin()) as f32; tr1.data[2*i] = (tr1.data[i] as f64 * (tr2.data[i] as f64).cos()) as f32; } tr1.trid = trid as i16; tr1.ns = (2*nt) as u16; },
                _ => err!("mysterious operation=\"%s\"", op),
            }
            puttr(&mut tr1);
            itr += 1;
        }
        if nbytes == 0 {
            warn!("%s still had traces when %s was exhausted", argv[1], argv[2]);
            warn!("processed %d pairs of traces before EOF", itr);
        } else if fgettr(&fp2, &mut tr2) != 0 {
            warn!("%s still had traces when %s was exhausted", argv[2], argv[1]);
            warn!("processed %d pairs of traces before EOF", itr);
        }
    } else {
        fgettr(&fp2, &mut tr2);
        if w2 != 0.0 { for i in 0..nt { tr2.data[i] *= w2; } }
        while fgettr(&fp1, &mut tr1) != 0 {
            if w1 != 0.0 { for i in 0..nt { tr1.data[i] *= w1; } }
            nt = tr1.ns as usize;
            if nt != tr2.ns as usize {
                warn!("trace %d:", itr);
                err!("%s and %s have different ns (%d vs %d)", argv[1], argv[2], tr1.ns, tr2.ns);
            }
            let complex = tr1.trid == 111 && tr2.trid == 111;
            match iop {
                8 => for i in 0..nt { tr1.data[i] -= tr2.data[i]; },
                9 if complex => { let mut i = 0; while i < nt { complex_mul(&mut tr1, &tr2, i); i += 2; } },
                9 => for i in 0..nt { tr1.data[i] *= tr2.data[i]; },
                10 if complex => { let mut i = 0; while i < nt { complex_div(&mut tr1, &tr2, i); i += 2; } },
                10 => for i in 0..nt { let denom = tr2.data[i]; if denom == 0.0 { tr1.data[i] = 0.0; } else { tr1.data[i] /= denom; } },
                _ => err!("mysterious operation=\"%s\"", op),
            }
            puttr(&mut tr1);
            itr += 1;
        }
    }
    exit(CWP_Exit());
}

fn complex_mul(a: &mut Trace, b: &Trace, i: usize) {
    let (x, y, c, d) = (a.data[i], a.data[i+1], b.data[i], b.data[i+1]);
    a.data[i] = x*c - y*d;
    a.data[i+1] = x*d + y*c;
}
fn complex_div(a: &mut Trace, b: &Trace, i: usize) {
    let (x, y, c, d) = (a.data[i], a.data[i+1], b.data[i], b.data[i+1]);
    if c == 0.0 && d == 0.0 { a.data[i] = 0.0; a.data[i+1] = 0.0; }
    else if c == 0.0 { a.data[i] = -y/d; }
    else if d == 0.0 { a.data[i+1] = y/c; }
    else { a.data[i] = x/c - y/d; a.data[i+1] = x/d + y/c; }
}
