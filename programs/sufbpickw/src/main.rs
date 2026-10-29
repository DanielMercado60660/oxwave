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

//! SUFBPICKW - First break auto picker, port of su/main/picking/sufbpickw.c (Revision 1.3).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_f_smooth::smooth_segmented_array;
use oxwave_gather_io::*;
use oxwave_linear_interp::intlin;
use oxwave_sblas::{isamax, sasum};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " SUFBPICKW - First break auto picker\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    "   sufbpickw < infile >outfile\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "  none\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " keyg=ep\t\t\t\t\t\t \t",
    " window=.03\tLength of forward and backward windows (s)\t",
    " test=1\tOutput the characteristic function\t \t",
    "\t\tThis can be used for testing window size\t",
    " Template\t\t\t\t\t\t\t",
    " o=\t\toffset...\t\t\t\t  \t",
    " t=\t\ttime pairs for defining first break search\t",
    "\t\t\twindow centre\t\t\t\t",
    " tdv=.05\tHalf length of the search window\t\t",
    "\t\t\t\t\t\t\t\t",
    " If the template is specified the maximum value of the\t\t",
    " characteristic function is searched in the window\t\t",
    "  defined by the template only.Default is the whole trace.\t",
    "\t\t\t\t\t\t\t  \t",
    "  The time of the pick is stored in header word unscale \t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let nowindow = countparval("t") == 0;
    let np = if nowindow {
        2u32
    } else {
        countparval("t") as u32
    };
    let mut t = Vec::<f32>::new();
    let mut o = Vec::<f32>::new();
    if !nowindow {
        t.resize(np as usize, 0.0);
        o.resize(np as usize, 0.0);
        if np == countparval("o") as u32 {
            getparfloat("t", &mut t);
            getparfloat("o", &mut o);
        } else {
            err!(" t and o has different number of elements\n");
        }
    }
    let mut keyg = String::new();
    if getparstring("keyg", &mut keyg) == 0 {
        keyg = "ep".into();
    }
    let mut window = 0.0f32;
    let mut tdv = 0.0f32;
    let mut test = 0i32;
    if getparfloat("window", &mut window) == 0 {
        window = 0.02;
    }
    if getparfloat("tdv", &mut tdv) == 0 {
        tdv = -1.0;
    }
    if getparint("test", &mut test) == 0 {
        test = 1;
    }
    checkpars();
    let mut state = GatherState::default();
    let mut typeg = "";
    let mut valg = Value::default();
    let (mut nt, mut ntr, mut first) = (0, 0, 0);
    let mut dt = 0.0f32;
    let mut rec_o = get_gather(
        &mut state, &keyg, &mut typeg, &mut valg, &mut nt, &mut ntr, &mut dt, &mut first,
    );
    let iwindow = nint((window / dt) as f64);
    let itdv = if tdv == -1.0 {
        nt
    } else {
        nint(2.0 * tdv as f64 / dt as f64)
    };
    if ntr == 0 {
        err!("Can't get first record\n");
    }
    while let Some(mut rec) = rec_o {
        // C allocates cf once per gather and clears it for each trace.
        // Its peak search may read past nt; zero padding avoids a Rust panic,
        // but cannot reproduce undefined contents beyond the C allocation.
        let extra = (itdv.max(0) as usize)
            .saturating_add(iwindow.max(0) as usize)
            .saturating_add(1);
        let mut cf = vec![0.0f32; nt as usize + extra];
        let mut find = vec![0.0f32; nt as usize];
        for it in 0..nt {
            find[it as usize] = it as f32;
        }
        for trace in &mut rec {
            cf[..nt as usize].fill(0.0);
            let ifbt = if nowindow {
                0
            } else {
                let mut fbt = [0.0];
                intlin(
                    &o,
                    &t,
                    t[0],
                    t[np as usize - 1],
                    &[trace.offset as f32],
                    &mut fbt,
                );
                nint((fbt[0] / dt) as f64)
            };
            let iw = iwindow.max(0) as usize;
            // The gather crate stores only nt samples. C allocates a full
            // segy for each trace, so retain padding for oversized windows.
            let mut data = trace.data.clone();
            data.resize(data.len().max(iw.saturating_mul(2).saturating_add(1)), 0.0);
            let mut ampb = sasum(iwindow, &data, 1) + f32::EPSILON;
            let mut ampf = sasum(iwindow, &data[iw..], 1) + f32::EPSILON;
            for it in iwindow..nt - iwindow - 1 {
                let it = it as usize;
                cf[it] = ampf / ampb;
                ampb = (ampb as f64 - (data[it - iw] as f64).abs()) as f32;
                ampf = (ampf as f64 - (data[it] as f64).abs()) as f32;
                ampb = (ampb as f64 + (data[it] as f64).abs()) as f32;
                ampf = (ampf as f64 + (data[it + iw] as f64).abs()) as f32;
            }
            // Smooth the characteristic function, including when test=0.
            // The C count is nt-iwindow-1, so the final sample is excluded.
            // A negative count becomes a huge size_t in C's ealloc1float.
            // The slice-based library API cannot express that failing call.
            if nt - iwindow - 1 < 0 {
                err!("ealloc.c: malloc failed (Cannot allocate memory)");
            }
            if nt - iwindow - 1 > 0 {
                smooth_segmented_array(
                    &find[iw..nt as usize - 1],
                    &mut cf[iw..nt as usize - 1],
                    iw,
                    1,
                    5,
                );
            }
            let it = (ifbt - itdv / 2).max(0).min(nt - 1) as usize;
            let mut itime = isamax(itdv, &cf[it..], 1) + it as i32;
            if itime == ifbt - itdv / 2 || itime == ifbt + itdv / 2 {
                itime = 0;
            }
            trace.unscale = dt * itime as f32;
            if test != 0 {
                trace.data[..nt as usize].copy_from_slice(&cf[..nt as usize]);
            }
        }
        put_gather(&rec, nt);
        rec_o = get_gather(
            &mut state, &keyg, &mut typeg, &mut valg, &mut nt, &mut ntr, &mut dt, &mut first,
        );
    }
    exit(CWP_Exit());
}
