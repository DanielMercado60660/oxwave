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

//! SUOP - unary arithmetic operations on traces.
//! Port of su/main/operations/suop.c (revision 1.32).
//! Credits: CWP: Shuki Ronen, Jack K Cohen; Toralf Foerster,
//! Reg Beardsley, Chris Liner, and others.
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUOP - do unary arithmetic operation on segys \t\t",
    " \t\t\t\t\t\t\t\t",
    " suop <stdin >stdout op=abs\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    "\top=abs\t\toperation flag\t\t\t\t",
    "\t\t\tabs   : absolute value\t\t\t",
    "\t\t\tavg   : remove average value\t\t",
    "\t\t\tssqrt : signed square root\t\t",
    "\t\t\tsqr   : square\t\t\t\t",
    "\t\t\tssqr  : signed square\t\t\t",
    "\t\t\tsgn   : signum function\t\t\t",
    "\t\t\texp   : exponentiate\t\t\t",
    "\t\t\tsexp  : signed exponentiate\t\t",
    "\t\t\tslog  : signed natural log\t\t",
    "\t\t\tslog2 : signed log base 2\t\t",
    "\t\t\tslog10: signed common log\t\t",
    "\t\t\tcos   : cosine\t\t\t\t",
    "\t\t\tsin   : sine\t\t\t\t",
    "\t\t\ttan   : tangent\t\t\t\t",
    "\t\t\tcosh  : hyperbolic cosine\t\t",
    "\t\t\tsinh  : hyperbolic sine\t\t\t",
    "\t\t\ttanh  : hyperbolic tangent\t\t",
    "\t\t\tcnorm : norm complex samples by modulus ",
    "\t\t\tnorm  : divide trace by Max. Value\t",
    "\t\t\tdb    : 20 * slog10 (data)\t\t",
    "\t\t\tneg   : negate value\t\t\t",
    "\t\t\tposonly : pass only positive values\t",
    "\t\t\tnegonly : pass only negative values\t",
    "                       sum   : running sum trace integration   ",
    "                       diff  : running diff trace differentiation",
    "                       refl  : (v[i+1] - v[i])/(v[i+1] + v[i]) ",
    "\t\t\tmod2pi : modulo 2 pi\t\t\t",
    "\t\t\tinv   : inverse\t\t\t\t",
    "\t\t\trmsamp : rms amplitude\t\t\t",
    "                       s2v   : sonic to velocity (ft/s) conversion     ",
    "                       s2vm  : sonic to velocity (m/s) conversion     ",
    "                       d2m   : density (g/cc) to metric (kg/m^3) conversion ",
    "                       drv2  : 2nd order vertical derivative ",
    "                       drv4  : 4th order vertical derivative ",
    "                       integ : top-down integration            ",
    "                       spike : local extrema to spikes         ",
    "                       saf   : spike and fill to next spike    ",
    "                       freq  : local dominant freqeuncy        ",
    "                       lnza  : preserve least non-zero amps    ",
    "                       --------- window operations ----------- ",
    "                       mean  : arithmetic mean                 ",
    "                       despike  : despiking based on median filter",
    "                       std   : standard deviation              ",
    "                       var   : variance                        ",
    "       nw=21           number of time samples in window        ",
    "                       --------------------------------------- ",
    "\t\t\tnop   : no operation\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Note:\tBinary ops are provided by suop2.\t\t\t",
    " Operations inv, slog, slog2, and slog10 are \"punctuated\",\t",
    " meaning that if, the input contains 0 values,\t\t\t",
    " 0 values are returned.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " For file operations on non-SU format binary files use:  farith",
];
use oxwave_core::prelude::*;
fn sign(x: f32) -> f64 {
    if x < 0.0 {
        -1.0
    } else {
        1.0
    }
}
fn abs(x: f32) -> f32 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let two = 2.0f32 * dt;
    let mut op = String::from("abs");
    getparstring("op", &mut op);
    let valid = [
        "abs", "ssqrt", "sqr", "ssqr", "sgn", "exp", "sexp", "slog", "slog2", "slog10", "cos",
        "sin", "tan", "cosh", "sinh", "tanh", "norm", "cnorm", "db", "neg", "nop", "posonly",
        "negonly", "sum", "refl", "diff", "mod2pi", "inv", "avg", "rmsamp", "std", "var", "mean",
        "despike", "s2v", "s2vm", "d2m", "drv2", "drv4", "integ", "spike", "freq", "lnza", "saf",
    ];
    if !valid.contains(&op.as_str()) {
        err!("unknown operation=\"%s\", see self-doc", op);
    }
    let mut nw = 21i32;
    getparint("nw", &mut nw);
    checkpars();
    if nw % 2 == 0 {
        nw += 1;
    }
    // C allocates the median buffer for every operation. Negative nw becomes
    // an impossible size_t allocation; report ealloc1float's ENOMEM diagnostic
    // even when Rust rejects the capacity before calling the allocator.
    let mut med = Vec::<f32>::new();
    if nw < 0 || med.try_reserve_exact(nw as usize).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // Unused median storage stays uninitialized, like malloc in C. In
    // particular, a very large positive nw must not force us to touch it.
    if op == "despike" {
        // C reads outside tmp for an oversized window. Reject that undefined
        // input before indexing; the other window operations handle it safely.
        if nw as usize > nt {
            err!("despike: nw=%d exceeds trace length nt=%d", nw, nt);
        }
        med.resize(nw as usize, 0.0);
    }
    let mut tmp = vec![0f32; nt + 1];
    loop {
        match op.as_str() {
            "abs" | "ssqrt" | "sqr" | "ssqr" | "sgn" | "exp" | "sexp" | "slog" | "slog2"
            | "slog10" | "cos" | "sin" | "tan" | "cosh" | "sinh" | "tanh" | "db" | "neg"
            | "posonly" | "negonly" | "inv" | "s2v" | "s2vm" | "d2m" => {
                for i in 0..nt {
                    let x = tr.data[i];
                    let a = abs(x) as f64;
                    let s = sign(x);
                    tr.data[i] = match op.as_str() {
                        "abs" => abs(x),
                        "ssqrt" => (s * a.sqrt()) as f32,
                        "sqr" => x * x,
                        "ssqr" => (s * x as f64 * x as f64) as f32,
                        "sgn" => s as f32,
                        "exp" => (x as f64).exp() as f32,
                        "sexp" => (s * a.exp()) as f32,
                        "slog" => {
                            if a > 0.0 {
                                (s * a.ln()) as f32
                            } else {
                                0.0
                            }
                        }
                        "slog2" => {
                            let y = if a > 0.0 {
                                (s * a.ln() / 2f64.ln()) as f32
                            } else {
                                0.0
                            };
                            if abs(y) > 0.0 {
                                (sign(y) * (abs(y) as f64).log10()) as f32
                            } else {
                                0.0
                            }
                        }
                        "slog10" => {
                            if a > 0.0 {
                                (s * a.log10()) as f32
                            } else {
                                0.0
                            }
                        }
                        "cos" => (x as f64).cos() as f32,
                        "sin" => (x as f64).sin() as f32,
                        "tan" => (x as f64).tan() as f32,
                        "cosh" => (x as f64).cosh() as f32,
                        "sinh" => (x as f64).sinh() as f32,
                        "tanh" => (x as f64).tanh() as f32,
                        "db" => {
                            if a > 0.0 {
                                (20.0 * s * a.log10()) as f32
                            } else {
                                0.0
                            }
                        }
                        "neg" => -x,
                        "posonly" => {
                            if x > 0.0 {
                                x
                            } else {
                                0.0
                            }
                        }
                        "negonly" => {
                            if x < 0.0 {
                                x
                            } else {
                                0.0
                            }
                        }
                        "inv" => {
                            if x != 0.0 {
                                (1.0 / x as f64) as f32
                            } else {
                                x
                            }
                        }
                        "s2v" => (1000000.0 / x as f64) as f32,
                        "s2vm" => (304800.0 / x as f64) as f32,
                        "d2m" => (1000.0 * x as f64) as f32,
                        _ => x,
                    };
                }
            }
            "norm" => {
                let mut max = 0.0;
                for i in 0..nt {
                    let x = abs(tr.data[i]);
                    if max < x {
                        max = x;
                    }
                }
                if max != 0.0 {
                    for i in 0..nt {
                        tr.data[i] /= max;
                    }
                }
            }
            "cnorm" => {
                if tr.trid != 111 {
                    let error = std::io::Error::last_os_error().to_string();
                    let message = error.split(" (os error ").next().unwrap_or(&error);
                    eprintln!("cnorm option only works for trid=111: {}", message);
                }
                // For odd nt, C also normalizes the sample just past ns. Keep
                // that behavior while it is inside the trace's sample array,
                // but reject the maximum odd ns instead of reading past it.
                if nt % 2 != 0 && nt == tr.data.len() {
                    err!("cnorm: ns=%d leaves an unpaired complex sample", nt);
                }
                for i in (0..nt).step_by(2) {
                    let x = ((tr.data[i] * tr.data[i] + tr.data[i + 1] * tr.data[i + 1]) as f64)
                        .sqrt() as f32;
                    tr.data[i] /= x;
                    tr.data[i + 1] /= x;
                }
            }
            "sum" => {
                for i in 1..nt {
                    tr.data[i] += tr.data[i - 1];
                }
            }
            "integ" => {
                let mut val = 0.0f32;
                for i in 0..nt {
                    val += tr.data[i];
                    tr.data[i] = val;
                }
            }
            "refl" => {
                for i in (1..nt).rev() {
                    let n = tr.data[i] - tr.data[i - 1];
                    let mut d = tr.data[i] + tr.data[i - 1];
                    if d == 0.0 {
                        d = 1.0;
                    }
                    tr.data[i] = n / d;
                }
                tr.data[0] = 0.0;
            }
            "diff" => {
                let mut d = vec![0f32; nt.max(3)];
                d[..nt].copy_from_slice(&tr.data[..nt]);
                for i in 2..nt.saturating_sub(2) {
                    tr.data[i] = (d[i + 1] - d[i - 1]) / two;
                }
                tr.data[0] = (d[1] - d[0]) / dt;
                if nt >= 2 {
                    tr.data[nt - 1] = (d[nt - 1] - d[nt - 2]) / dt;
                }
                tr.data[1] = (d[2] - d[0]) / two;
                if nt >= 3 {
                    tr.data[nt - 2] = (d[nt - 1] - d[nt - 3]) / two;
                }
            }
            "avg" | "rmsamp" => {
                let mut sum = 0f32;
                for i in 0..nt {
                    sum += if op == "avg" {
                        tr.data[i]
                    } else {
                        tr.data[i] * tr.data[i]
                    };
                }
                if op == "avg" {
                    let a = sum / nt as f32;
                    for i in 0..nt {
                        tr.data[i] -= a;
                    }
                } else {
                    let a = ((sum / nt as f32) as f64).sqrt() as f32;
                    tr.data[..nt].fill(0.0);
                    tr.data[0] = a;
                }
            }
            "mod2pi" => {
                for i in 0..nt {
                    while tr.data[i] < 0.0 {
                        tr.data[i] = (tr.data[i] as f64 + 2.0 * 3.141592653589793) as f32;
                    }
                    while tr.data[i] as f64 >= 2.0 * 3.141592653589793 {
                        tr.data[i] = (tr.data[i] as f64 - 2.0 * 3.141592653589793) as f32;
                    }
                }
            }
            "lnza" => {
                tmp[..nt].copy_from_slice(&tr.data[..nt]);
                for i in 1..nt.saturating_sub(1) {
                    let x = tr.data[i - 1];
                    let y = tr.data[i];
                    let z = tr.data[i + 1];
                    if sign(x) == sign(y) && sign(y) == sign(z) {
                        tmp[i] = 0.0;
                    }
                    if sign(x) == sign(y) && sign(y) != sign(z) {
                        tmp[i] = if y.abs() < z.abs() { y } else { z };
                    }
                    if sign(x) != sign(y) && sign(y) == sign(z) {
                        tmp[i] = if x.abs() < y.abs() { x } else { y };
                    }
                }
                tmp[0] = 0.0;
                tmp[nt] = 0.0;
                tr.data[..nt].copy_from_slice(&tmp[..nt]);
            }
            "spike" | "saf" | "freq" => {
                let mut old = 0;
                let mut vold = 0.0;
                for i in 1..nt.saturating_sub(1) {
                    let x = tr.data[i - 1];
                    let y = tr.data[i];
                    let z = tr.data[i + 1];
                    let high = y > x && y > z;
                    let low = y < x && y < z;
                    tmp[i] = if op == "freq" {
                        if high {
                            1.0
                        } else if low {
                            -1.0
                        } else {
                            0.0
                        }
                    } else if high || low {
                        y
                    } else {
                        0.0
                    };
                    if op == "saf" && (high || low) {
                        for j in old..i {
                            tmp[j] = vold;
                        }
                        old = i;
                        vold = y;
                    }
                }
                tmp[0] = 0.0;
                tmp[nt] = 0.0;
                if op == "freq" {
                    // The original dereferences its unallocated tmp1 buffer.
                    let _ = std::process::Command::new("/bin/kill")
                        .args(["-s", "SEGV", &std::process::id().to_string()])
                        .status();
                } else {
                    tr.data[..nt].copy_from_slice(&tmp[..nt]);
                }
            }
            "drv2" => {
                for i in 1..nt.saturating_sub(1) {
                    let x = tr.data[i - 1] - tr.data[i];
                    tmp[i] = (x as f64 / (2.0 * dt as f64)) as f32;
                }
                if nt > 1 {
                    tmp[0] = tmp[1];
                    tmp[nt] = tmp[nt - 1];
                }
                for i in 1..nt.saturating_sub(1) {
                    tr.data[i] = tmp[i];
                }
            }
            "drv4" => {
                for i in 2..nt.saturating_sub(2) {
                    let x = (tr.data[i - 2] as f64 - 8.0 * tr.data[i - 1] as f64
                        + 8.0 * tr.data[i + 1] as f64
                        - tr.data[i + 2] as f64) as f32;
                    tmp[i] = -x / (12.0f32 * dt);
                }
                if nt >= 3 {
                    tmp[0] = tmp[2];
                    tmp[1] = tmp[2];
                    tmp[nt - 1] = tmp[nt - 2];
                    tmp[nt] = tmp[nt - 2];
                }
                for i in 2..nt.saturating_sub(2) {
                    tr.data[i] = tmp[i];
                }
            }
            "mean" | "var" | "std" | "despike" => {
                tmp[..nt].copy_from_slice(&tr.data[..nt]);
                let mt = (nw - 1) / 2;
                for i in 0..nt {
                    let ii = i as i32;
                    if op == "despike" {
                        for j in 0..nw {
                            let k = if ii - mt > 0 && ii + mt < nt as i32 {
                                ii - mt + j
                            } else if ii - mt < 0 {
                                j
                            } else {
                                nt as i32 - nw + j
                            };
                            med[j as usize] = tmp[k as usize];
                        }
                        med.sort_by(|a, b| {
                            if a < b {
                                std::cmp::Ordering::Less
                            } else if a > b {
                                std::cmp::Ordering::Greater
                            } else {
                                std::cmp::Ordering::Equal
                            }
                        });
                        tr.data[i] = med[(nw / 2) as usize];
                    } else if ii - mt > 0 && ii + mt < nt as i32 {
                        let mut val = 0f32;
                        for j in ii - mt..=ii + mt {
                            val += tmp[j as usize];
                        }
                        if op == "mean" {
                            tr.data[i] = val / nw as f32;
                        } else {
                            let mean = val / nw as f32;
                            for j in ii - mt..=ii + mt {
                                let d = tmp[j as usize] - mean;
                                val += d * d;
                            }
                            let v = val / nw as f32;
                            tr.data[i] = if op == "std" {
                                (v as f64).sqrt() as f32
                            } else {
                                v
                            };
                        }
                    } else {
                        tr.data[i] = 0.0;
                    }
                }
            }
            _ => {}
        }
        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
