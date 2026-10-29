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

//! Port of su/lib/tabplot.c (revision 1.5).
//! Authors: CWP: Brian Sumner, Jack K. Cohen.

/// Plot selected sample points from a trace to standard output.
///
/// C prototype: `void tabplot(segy *tp, int itmin, int itmax);`
/// `data` is the trace's sample array, and both bounds are inclusive.
/// The bounds must refer to valid samples, as required by the C function.
pub fn tabplot(data: &[f32], itmin: i32, itmax: i32) {
    const PLOTWIDTH: i32 = 29;
    let abs = |x: f32| if x < 0.0 { -x } else { x };
    let mut amp = abs(data[itmin as usize]);
    for i in (itmin + 1)..=itmax {
        let value = abs(data[i as usize]);
        amp = if amp > value { amp } else { value };
    }

    for i in itmin..=itmax {
        let (val, plt) = if amp == 0.0 {
            (0.0, PLOTWIDTH)
        } else {
            let val = data[i as usize];
            let plt = (PLOTWIDTH as f64 * ((val / amp) as f64 + 1.0)) as i32;
            (val, plt)
        };
        let plot = if plt < PLOTWIDTH {
            format!("{}{}|", " ".repeat(plt as usize), "-".repeat((PLOTWIDTH - plt) as usize))
        } else if plt == PLOTWIDTH {
            format!("{}*", " ".repeat(PLOTWIDTH as usize))
        } else {
            format!("{}|{}", " ".repeat(PLOTWIDTH as usize), "+".repeat((plt - PLOTWIDTH) as usize))
        };
        oxwave_core::printf!("%5d %11.4e%s\n", i.wrapping_add(1), val, plot);
    }
}
