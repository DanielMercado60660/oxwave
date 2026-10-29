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

//! SUANALYTIC - generate an analytic trace using the Hilbert transform.
//! Port of su/main/transforms/suanalytic.c (Revision 1.2).
//! Credits: CWP: John Stockwell, based on suhilb by Jack K. Cohen.

use oxwave_core::prelude::*;
use oxwave_hilbert::hilbert;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUANALYTIC - use the Hilbert transform to generate an ANALYTIC\t",
    "\t\t(complex) trace\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suanalytic <stdin >sdout \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t       \t\t",
    " Optional Parameter:\t\t\t\t\t\t\t",
    " phaserot=\t\tphase rotation in degrees of complex trace\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The output are complex valued traces. The analytic trace is defined as",
    "   ctr[ i ] = indata[i] + i hilb[indata[t]]\t\t\t\t",
    " where the imaginary part is the hilbert tranform of the original trace",
    " \t\t\t\t\t\t\t\t\t",
    " The Hilbert transform is computed in the direct (time) domain\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " If phaserot is set, then a phase rotated complex trace is produced\t",
    "   ctr[ i ] = cos[phaserot]*indata[i] + i sin[phaserot]* hilb[indata[t]]",
    " \t\t\t\t\t\t\t\t\t",
    " Use \"suamp\" to extract real, imaginary, amplitude (modulus), etc \t",
    " Exmple:\t\t\t\t\t\t\t\t",
    " suanalytic < sudata | suamp mode=amp | suxgraph \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Use \"suattributes\" for instantaneous phase, frequency, etc.\t\t",
    "\t\t\t\t\t\t\t       \t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let mut phaserot = 0.0f32;
    let is_phaserot = getparfloat("phaserot", &mut phaserot) != 0;
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    if !seismic { warn!("input is not seismic data, trid=%d", tr.trid); }
    // C overruns tr.data at nt >= 32768. Reject before either output branch.
    if nt > tr.data.len() / 2 {
        err!("nt=%d is too large for analytic output (maximum %d)", nt, tr.data.len() / 2);
    }
    let mut data = vec![0.0f32; nt];
    let mut hdata = vec![0.0f32; nt];
    loop {
        data.fill(0.0);
        hdata.fill(0.0);
        for i in 0..nt { data[i] = tr.data[i]; }
        hilbert(&data, &mut hdata);
        if !is_phaserot {
            for i in 0..nt {
                tr.data[2 * i] = data[i];
                tr.data[2 * i + 1] = hdata[i];
            }
        } else {
            for i in 0..nt {
                tr.data[2 * i] = (data[i] as f64 * (3.141592653589793_f64 * phaserot as f64 / 180.0).cos()) as f32;
                tr.data[2 * i + 1] = (hdata[i] as f64 * (3.141592653589793_f64 * phaserot as f64 / 180.0).sin()) as f32;
            }
        }
        tr.trid = 111;
        tr.ns = (2 * nt) as u16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
