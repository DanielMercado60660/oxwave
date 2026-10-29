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

//! SUMATRIX: port of su/main/windowing_sorting_muting/sumatrix.c, revision 1.0.
//! Copyright (c) Colorado School of Mines, 2023. All rights reserved.
//! Author Andre Latour, August 2024; started from surescsv.c for its shot and receiver identifier code.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_header_case::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUMATRIX - Assemble a Matrix of Trace Sample Zones.                   ",
    "                                                                       ",
    " Notes:                                                                ",
    " Typically this program creates output traces which can be displayed   ",
    " to highlight surface-consistent differences and/or geometry-errors.   ",
    " Each output trace has time zones assembled from multiple input traces.",
    "                                                                       ",
    " ----------------------------------------------------------------------",
    "                                                                       ",
    " sumatrix  <stdin                                                      ",
    "                                                                       ",
    " skeyloc=fldr    List of keys that identify source locations.          ",
    "           Note: For 3d surveys you can use fldr (but possibly better  ",
    "                 to use grnofr,grnlof if sugeomcsv used SPS2 files).   ",
    " \t\t\t\t\t\t\t\t\t",
    " sdivider=1.0    List of dividers for skeyloc values. The values in    ",
    "                 this list can be set to -1 to reverse the order of    ",
    "                 sources in the output matrix.                         ",
    " \t\t\t\t\t\t\t\t\t",
    " rkeyloc=gaps    List of keys that identify receiver locations.        ",
    "           Note: For 3d surveys you often MUST use 2 keys              ",
    "                 (such as grnors,gaps if sugeomcsv used SPS2 files)    ",
    " \t\t\t\t\t\t\t\t\t",
    " rdivider=1.0    List of dividers for rkeyloc values. The values in    ",
    "                 this list can be set to -1 to reverse the order of    ",
    "                 receivers in the output matrix.                       ",
    " \t\t\t\t\t\t\t\t\t",
    " rfill=1         Output zeroed traces when rkeyloc combination         ",
    "                 has no input trace within a numzone output range.     ",
    "             =0  Do not output zeroed traces when rkeyloc combination  ",
    "                 has no input trace within a numzone output range.     ",
    "                 This option might be better for 3d surveys.           ",
    "             =2  Output zeroed traces for all rkeyloc combinations.    ",
    "                 Note this differs from option 1 in that it ALSO makes ",
    "                 traces OUTSIDE each numzone output range.             ",
    "           Note: None of these options creates an output trace for a   ",
    "                 rkeyloc combination THAT DOES NOT EXIST in input.     ",
    " \t\t\t\t\t\t\t\t\t",
    " lenzone=400.0   Length of Zones (ms., rounded to nearest sample.)     ",
    " \t\t\t\t\t\t\t\t\t",
    " numzone=20      Number of Zones To Put In Each Output Trace.          ",
    "           Note: Output traces have lenzone*numzone time length.       ",
    "                 Do not exceed the maximum trace samples of SU or      ",
    "                 other software that needs to handle these traces.     ",
    " \t\t\t\t\t\t\t\t\t",
    " maxsources=10000 Maximum Amount of Sources. This is also the          ",
    "                  maximum unique combinations of skeyloc values.       ",
    "                                                                       ",
    " maxreceivers=20000 Maximum Amount of Receivers. This is also the      ",
    "                    maximum unique combinations of rkeyloc values.     ",
    "                                                                       ",
    " maxtraces=1000000 Maximum Amount of Traces. Note that most memory for ",
    "                   the traces is only allocated for the actual amount  ",
    "                   of traces read (and just the lenzone, no headers).  ",
    "                                                                       ",
    "                                                                       ",
    " NOTE 1:  Output is sets of traces with numzone zones per each trace.  ",
    "          Each zone is made from 1 combination of skeyloc values.      ",
    "          Completely missing skeyloc combinations are not output.      ",
    "          A set of output traces is made for each increment of numzone.",
    "          Each set makes an output trace for each rkeyloc combination  ",
    "          that has an input trace in the numzone range of that set.    ",
    "                                                                       ",
    " NOTE 2:  Output trace headers are COPIES of the FIRST input header    ",
    "          with some values reset. The output skeyloc values are set    ",
    "          to the skeyloc value combination of the TOP ZONE. The output ",
    "          rkeyloc keys are set to the rkeyloc value combination for    ",
    "          each trace. And key nhs is set to the number of input trace  ",
    "          zones copied to each output trace.                           ",
    "                                                                       ",
    " NOTE 3:  For 3D surveys you may want to reverse only the line or      ",
    "          point ordering. For example, a divider list like 1,-1        ",
    "                                                                       ",
    " NOTE 4:  Usually, specify source keys in skeyloc and receiver keys    ",
    "          in rkeyloc. But you can also specify them vice-versa in      ",
    "          order to exchange the matrix orientation.                    ",
    "                                                                       ",
    " NOTE 5:  For some unusual situations you might want to set the        ",
    "          sdivider and/or rdivider to scale the input values.          ",
    "          The internal computation is: int(floor(header_value/divider))",
    "          which essentially means scale and round down to near integer.",
    "                                                                       ",
];

fn bhighii(a: &[Vec<i32>], g: &[i32]) -> usize {
    let (mut l, mut h) = (0, a.len());
    while l < h {
        let m = l + (h - l) / 2;
        if g >= a[m].as_slice() {
            l = m + 1
        } else {
            h = m
        }
    }
    l
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut sk = Vec::new();
    if countparval("skeyloc") > 0 {
        getparstringarray("skeyloc", &mut sk);
    } else {
        sk.push("fldr".to_string());
    }
    let mut sc = Vec::new();
    for k in &sk {
        let c = GetCase(k);
        if c < 1 {
            err!(
                "**** Error: Specified skeyloc name %s is not recognized.",
                k
            );
        }
        sc.push(c);
    }
    let mut ms = 10000;
    if countparval("maxsources") > 0 {
        getparint("maxsources", &mut ms);
    }
    if ms < 1 {
        err!("**** Error: maxsources must be greater than 0");
    }
    let mut sd = vec![1.0f64; sk.len()];
    if countparval("sdivider") > 0 {
        if countparval("sdivider") != sk.len() as i32 {
            err!("**** Error: sdivider list must have same amount as skeyloc ");
        }
        getpardouble("sdivider", &mut sd);
        for v in &sd {
            if *v == 0.0 {
                err!("**** Error: sdivider values cannot be 0.0 ");
            }
        }
    }
    let mut rk = Vec::new();
    if countparval("rkeyloc") > 0 {
        getparstringarray("rkeyloc", &mut rk);
    } else {
        rk.push("gaps".to_string());
    }
    let mut rc = Vec::new();
    for k in &rk {
        let c = GetCase(k);
        if c < 1 {
            err!(
                "**** Error: Specified rkeyloc name %s is not recognized.",
                k
            );
        }
        rc.push(c);
    }
    let mut mr = 20000;
    if countparval("maxreceivers") > 0 {
        getparint("maxreceivers", &mut mr);
    }
    if mr < 1 {
        err!("**** Error: maxreceivers must be greater than 0");
    }
    let mut rd = vec![1.0f64; rk.len()];
    if countparval("rdivider") > 0 {
        if countparval("rdivider") != rk.len() as i32 {
            err!("**** Error: rdivider list must have same amount as rkeyloc ");
        }
        getpardouble("rdivider", &mut rd);
        for v in &rd {
            if *v == 0.0 {
                err!("**** Error: rdivider values cannot be 0.0 ");
            }
        }
    }
    let mut fill = 1;
    if countparval("rfill") > 0 {
        getparint("rfill", &mut fill);
    }
    if fill < 0 || fill > 2 {
        err!("**** Error: rfill parameter out of range.");
    }
    let mut mt = 1000000;
    if countparval("maxtraces") > 0 {
        getparint("maxtraces", &mut mt);
    }
    if mt < 1 {
        err!("**** Error: maxtraces must be greater than 0");
    }
    let mut fl = 400.0f32;
    if countparval("lenzone") > 0 {
        getparfloat("lenzone", &mut fl);
    }
    if fl < 1.0 {
        err!("**** Error: lenzone must be greater than 0");
    }
    let mut nz = 20;
    if countparval("numzone") > 0 {
        getparint("numzone", &mut nz);
    }
    if nz < 1 {
        err!("**** Error: numzone must be greater than 0");
    }
    checkpars();
    let (mut tr, mut out) = (Trace::new(), Trace::new());
    let (mut su, mut ru): (Vec<Vec<i32>>, Vec<Vec<i32>>) = (Vec::new(), Vec::new());
    let (mut st, mut rt): (Vec<Vec<i32>>, Vec<Vec<i32>>) = (Vec::new(), Vec::new());
    let mut sam: Vec<Vec<f32>> = Vec::new();
    let (mut len, mut lensam) = (0usize, 0usize);
    while gettr(&mut tr) != 0 {
        let n = sam.len();
        if n == 0 {
            out.hdr = tr.hdr;
            let zone_samples = fl as f64 / (tr.dt as f64 / 1000.0);
            // The C conversion and subsequent writes are undefined for these sizes.
            if !zone_samples.is_finite() || zone_samples > i32::MAX as f64 - 0.5 {
                err!("**** Error: output trace exceeds sample capacity or has invalid dt");
            }
            len = nint(zone_samples) as usize;
            let output_samples = match len.checked_mul(nz as usize) {
                Some(count) if count <= out.data.len() && count <= u16::MAX as usize => count,
                _ => err!("**** Error: output trace exceeds sample capacity or has invalid dt"),
            };
            out.ns = output_samples as u16;
            out.dt = tr.dt;
            lensam = (tr.ns as usize).min(len);
        }
        if n >= mt as usize {
            err!(
                "error: Number of input traces is greater than maxtraces=%d",
                mt
            );
        }
        let mut sg = vec![0; sc.len()];
        for i in 0..sg.len() {
            sg[i] = (fromhead(&tr, sc[i]) / sd[i]).floor() as i32;
        }
        let j = bhighii(&su, &sg);
        if j == 0 || sg != su[j - 1] {
            if su.len() >= ms as usize {
                err!(
                    "error: At input trace=%d  number of unique skeyloc combinations is greater than maxsources=%d",
                    n as i32,
                    ms
                );
            }
            su.insert(j, sg.clone());
        }
        let mut rg = vec![0; rc.len()];
        for i in 0..rg.len() {
            rg[i] = (fromhead(&tr, rc[i]) / rd[i]).floor() as i32;
        }
        let k = bhighii(&ru, &rg);
        if k == 0 || rg != ru[k - 1] {
            if ru.len() >= mr as usize {
                err!(
                    "error: At input trace=%d  number of unique rkeyloc combinations is greater than maxreceivers=%d",
                    n as i32,
                    mr
                );
            }
            ru.insert(k, rg.clone());
        }
        sam.push(tr.data[..lensam].to_vec());
        st.push(sg);
        rt.push(rg);
    }
    let (ns, nr) = (su.len(), ru.len());
    let mut jk = vec![vec![-1i32; nr]; ns];
    let mut dup = 0;
    for i in 0..sam.len() {
        let j = bhighii(&su, &st[i]) - 1;
        let k = bhighii(&ru, &rt[i]) - 1;
        if jk[j][k] != -1 {
            dup += 1;
        }
        jk[j][k] = i as i32;
    }
    let mut n = 0;
    while n < ns {
        let end = (n + nz as usize).min(ns);
        let (mut kb, mut ke) = (0, nr);
        if fill < 2 {
            kb = nr;
            ke = 0;
            for k in 0..nr {
                for j in n..end {
                    if jk[j][k] > -1 {
                        if kb == nr {
                            kb = k;
                        }
                        ke = k + 1;
                    }
                }
            }
        }
        for i in 0..sc.len() {
            tohead(&mut out, sc[i], su[n][i] as f64 * sd[i]);
        }
        for k in kb..ke {
            let mut fold = 0;
            for i in 0..len * nz as usize {
                out.data[i] = 0.0;
            }
            for j in n..end {
                if jk[j][k] > -1 {
                    fold += 1;
                    for i in 0..lensam {
                        out.data[i + (j - n) * len] = sam[jk[j][k] as usize][i];
                    }
                }
            }
            if fold > 0 || fill > 0 {
                for i in 0..rc.len() {
                    tohead(&mut out, rc[i], ru[k][i] as f64 * rd[i]);
                }
                out.nhs = fold;
                puttr(&mut out);
            }
        }
        n += nz as usize;
    }
    warn!(
        "Unique combinations of:  skeyloc values=%d   rkeyloc values=%d",
        ns as i32, nr as i32
    );
    if dup > 0 {
        warn!(
            "There were %d traces that REPLACED another trace in the SAME matrix location.",
            dup
        );
    }
    exit(CWP_Exit());
}
