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

//! SUAMP - extract amplitude, phase, real or imaginary frequency traces.
//! Port of su/main/transforms/suamp.c (Revision 1.24).
//! Credits: CWP: Shuki Ronen, Jack K. Cohen c.1986.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::closeto;
use oxwave_f_smooth::smooth_1;
use oxwave_phase_unwrap::{oppenheim_unwrap_phase, simple_unwrap_phase};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUAMP - output amp, phase, real or imag trace from\t\t\t",
    " \t(frequency, x) domain data\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suamp <stdin >stdout mode=amp\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t\t",
    " mode=amp\toutput flag\t\t \t\t\t\t",
    " \t\t=amp\toutput amplitude traces\t\t\t\t",
    " \t\t=logamp\toutput log(amplitude) traces\t\t\t",
    " \t\t\t=phase\toutput phase traces\t\t\t",
    " \t\t\t=ouphase output unwrapped phase traces (oppenheim)",
    " \t\t\t=suphase output unwrapped phase traces (simple)\t",
    " \t\t\t=real\toutput real parts\t\t\t",
    " \t     \t=imag\toutput imag parts\t \t\t\t",
    " jack=0\t=1  divide value at zero frequency by 2   \t\t",
    "\t\t(operative only for mode=amp)\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " .... phase unwrapping options\t..... \t\t\t\t\t",
    " unwrap=1\t |dphase| > pi/unwrap constitutes a phase wrapping\t",
    "\t\t\t(operative only for mode=suphase)\t\t",
    " trend=1\tremove linear trend from the unwrapped phase\t\t",
    " zeromean=0\tassume phase(0)=0.0, else assume phase is zero mean\t",
    " smooth=0\tapply damped least squares smoothing to unwrapped phase ",
    " r=10.0\t    ... damping coefficient, only active when smooth=1\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " \tThe trace returned is half length from 0 to Nyquist. \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " \tsufft <data | suamp >amp_traces\t\t\t\t\t",
    " Example: \t\t\t\t\t\t\t\t",
    "\tsufft < data > complex_traces\t\t\t\t\t",
    " \t suamp < complex_traces mode=real > real_traces\t\t\t",
    " \t suamp < complex_traces mode=imag > imag_traces\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Note: the inverse of the above operation is: \t\t\t\t",
    "\tsuop2 real_traces imag_traces op=zipper > complex_traces\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Note: Explanation of jack=1 \t\t\t\t\t\t",
    " The amplitude spectrum is the modulus of the complex output of\t",
    " the fft. f(0) is thus the average over the range of integration\t",
    " of the transform. For causal functions, or equivalently, half\t\t",
    " transforms, f(0) is 1/2 of the average over the full range.\t\t",
    " Most oscillatory functions encountered in wave applications are\t",
    " zero mean, so this is usually not an issue.\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Note: Phase unwrapping: \t\t\t\t\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " The mode=ouphase uses the phase unwrapping method of Oppenheim and\t",
    " Schaffer, 1975. \t\t\t\t\t\t\t",
    " The mode=suphase generates unwrapped phase assuming that jumps\t",
    " in phase larger than pi/unwrap constitute a phase wrapping.\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    // Get information from first trace.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    if tr.trid != 111 {
        err!("input not complex freq data, trid=%d", tr.trid);
    }
    let nf = tr.ns as usize;
    let nfby2 = nf / 2;
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut df = tr.d1;
    if df == 0.0 {
        let mut dt = (tr.dt as f64 / 1000000.0) as f32;
        if dt == 0.0 {
            dt = 0.004;
            if verbose != 0 {
                warn!("dt or d1 not set, assumed to be .004");
            }
        }
        df = (1.0 / nf as f64 * dt as f64) as f32;
    }
    let mut mode = String::new();
    if getparstring("mode", &mut mode) == 0 {
        mode = "amp".into();
    }
    let mut jack = 0;
    let mut smooth = 0;
    let mut r = 10.0f32;
    getparint("jack", &mut jack);
    getparint("smooth", &mut smooth);
    getparfloat("r", &mut r);
    let imode = match mode.as_str() {
        "real" => 1,
        "imag" => 2,
        "amp" => 3,
        "phase" => 4,
        "suphase" => 5,
        "ouphase" => 6,
        "logamp" => 7,
        _ => err!("unknown mode=\"%s\", see self-doc", mode),
    };
    let mut trend = 1;
    let mut zeromean = 0;
    let mut unwrap = 1.0f32;
    if imode == 5 || imode == 6 {
        getparint("trend", &mut trend);
        getparint("zeromean", &mut zeromean);
        getparfloat("unwrap", &mut unwrap);
    }
    checkpars();
    // Allocate arrays for real and imaginary parts.
    let mut xr = vec![0.0f32; nfby2];
    let mut xi = vec![0.0f32; nfby2];
    loop {
        for i in 0..nfby2 {
            xr[i] = tr.data[2 * i];
            xi[i] = tr.data[2 * i + 1];
        }
        match imode {
            1 | 2 => {
                for i in 0..nfby2 {
                    tr.data[i] = if imode == 1 { xr[i] } else { xi[i] };
                }
                tr.trid = if imode == 1 { 116 } else { 117 };
            }
            3 | 7 => {
                for i in 0..nfby2 {
                    let re = xr[i];
                    let im = xi[i];
                    let amp = ((re * re + im * im) as f64).sqrt() as f32;
                    // C leaves the zero-frequency sample as amplitude even
                    // in LOGAMP mode, including the optional jack division.
                    tr.data[i] = if i == 0 && jack != 0 {
                        (amp as f64 / 2.0) as f32
                    } else if i > 0 && imode == 7 {
                        if closeto(amp, 0.0) {
                            0.0
                        } else {
                            (amp as f64).ln() as f32
                        }
                    } else {
                        amp
                    };
                }
                tr.trid = if imode == 3 { 118 } else { 126 };
            }
            4 => {
                for i in 0..nfby2 {
                    let re = xr[i];
                    let im = xi[i];
                    tr.data[i] = if re * re + im * im != 0.0 {
                        (im as f64).atan2(re as f64) as f32
                    } else {
                        0.0
                    };
                }
                tr.trid = 119;
            }
            5 => {
                // Allocate space for the phase.
                let mut phase = vec![0.0f32; nfby2];
                let mut xout = vec![0.0f32; nfby2];
                for i in 0..nfby2 {
                    let re = xr[i];
                    let im = xi[i];
                    xout[i] = i as f32;
                    phase[i] = if re * re + im * im != 0.0 {
                        (im as f64).atan2(re as f64) as f32
                    } else {
                        0.0
                    };
                }
                // Unwrap the phase, then smooth even if unwrapping is disabled.
                if unwrap != 0.0 && nfby2 != 0 {
                    simple_unwrap_phase(trend, zeromean, unwrap, &mut phase);
                }
                if smooth != 0 {
                    smooth_1(&xout, &mut phase, r);
                }
                // Write unwrapped phase.
                for i in 0..nfby2 {
                    tr.data[i] = phase[i];
                }
                tr.trid = 119;
            }
            6 => {
                let mut xout = vec![0.0f32; nfby2];
                let mut phase = vec![0.0f32; nfby2];
                // Oppenheim starts with zero phase, including when unwrap=0.
                if unwrap != 0.0 && nfby2 != 0 {
                    oppenheim_unwrap_phase(trend, zeromean, df, &xr, &xi, &mut phase);
                }
                for i in 0..nfby2 {
                    xout[i] = i as f32;
                }
                if smooth != 0 {
                    smooth_1(&xout, &mut phase, r);
                }
                for i in 0..nfby2 {
                    tr.data[i] = phase[i];
                }
                tr.trid = 119;
            }
            _ => unreachable!(),
        }
        // Output the half-length trace.
        tr.ns = nfby2 as u16;
        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
