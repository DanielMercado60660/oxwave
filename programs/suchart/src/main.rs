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
//! SUCHART - prepare data for x vs. y plot.
//! Port of su/main/headers/suchart.c (Revision 1.19).
//! Credits: SEP: Einar Kjartansson; CWP: Jack K. Cohen.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUCHART - prepare data for x vs. y plot\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " suchart <stdin >stdout key1=sx key2=gx\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tkey1=sx  \tabscissa \t\t\t\t",
    " \tkey2=gx\t\tordinate\t\t\t\t",
    "\toutpar=null\tname of parameter file\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " The output is the (x, y) pairs of binary floats\t\t",
    " \t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t",
    " suchart < sudata outpar=pfile >plot_data\t\t\t",
    " psgraph <plot_data par=pfile title=\"CMG\" \\\t\t\t",
    "\tlinewidth=0 marksize=2 mark=8 | ...\t\t\t",
    " rm plot_data \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " suchart < sudata | psgraph n=1024 d1=.004 \\\t\t\t",
    "\tlinewidth=0 marksize=2 mark=8 | ...\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " fold chart: \t\t\t\t\t\t\t",
    " suchart < stacked_data key1=cdp key2=nhs |\t\t\t",
    "            psgraph n=NUMBER_OF_TRACES d1=.004 \\\t\t",
    "\tlinewidth=0 marksize=2 mark=8 > chart.ps\t\t",
    " \t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if isatty(fileno(&stdout())) {
        err!("must redirect or pipe binary output");
    }

    let mut key1 = String::new();
    let mut key2 = String::new();
    if getparstring("key1", &mut key1) == 0 { key1 = "sx".to_string(); }
    if getparstring("key2", &mut key2) == 0 { key2 = "gx".to_string(); }
    let type1 = hdtype(&key1);
    let type2 = hdtype(&key2);
    let index1 = getindex(&key1);
    let index2 = getindex(&key2);

    let mut npairs: i32 = 0;
    while gettr(&mut tr) != 0 {
        let mut val1 = Value::default();
        let mut val2 = Value::default();
        gethval(&tr, index1, &mut val1);
        gethval(&tr, index2, &mut val2);
        let x = vtof(type1, val1);
        let y = vtof(type2, val2);
        efwrite(&[x], &stdout());
        efwrite(&[y], &stdout());
        npairs = npairs.wrapping_add(1);
    }

    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) != 0 {
        let fp = efopen(&outpar, "w");
        fprintf!(&fp, "n=%d label1=%s label2=%s\n", npairs, key1, key2);
    }
    checkpars();
    exit(CWP_Exit());
}
