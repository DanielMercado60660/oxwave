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
//! FARITH - file arithmetic. Port of par/main/data_conversion/farith.c (Revision 1.21).
//! Credits: Dave Hale, Zhaobo Meng, Tony Kocurko, John Stockwell.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " FARITH - File ARITHmetic -- perform simple arithmetic with binary files",
    " \t\t\t\t\t\t\t\t\t",
    " farith <infile >outfile [optional parameters]\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " in=stdin\tinput file\t\t\t\t\t\t",
    " out=stdout\toutput file\t\t\t\t\t\t",
    " in2=\t   second input file (required for binary operations)\t\t",
    "\t\t   if it can't be opened as a file, it might be a scalar",
    " n=size_of_in,  fastest dimension (used only for op=cartprod is set)\t",
    " isig=\t\t index at which signum function acts (used only for \t",
    "\t\t\top=signum)\t\t\t\t\t",
    " scale=\tvalue to scale in by, used only for op=scale)\t\t",
    " bias=\t\tvalue to bias in by, used only for op=bias)\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " op=noop   noop for out = in\t\t\t\t\t\t",
    "\t   neg  for out = -in\t\t\t\t\t\t",
    "\t   abs  for out = abs(in)\t\t\t\t\t",
    "\t   scale for out = in *scale\t\t\t\t\t",
    "\t   bias for out = in + bias \t\t\t\t\t",
    "\t   exp  for out = exp(in)\t\t\t\t\t",
    "\t   sin  for out = sin(in)\t\t\t\t\t",
    "\t   cos  for out = cos(in)\t\t\t\t\t",
    "\t   log  for out = log(in)\t\t\t\t\t",
    "\t   sqrt for out = (signed) sqrt(in)\t\t\t\t",
    "\t   sqr  for out = in*in\t\t\t\t\t\t",
    "\t   degrad  for out = in*PI/180\t\t\t\t\t",
    "\t   raddeg  for out = in*180/PI\t\t\t\t\t",
    "\t   pinv  for out = (punctuated) 1 / in   \t\t\t",
    "\t   pinvsqr  for out = (punctuated) 1 /in*in \t\t\t",
    "\t   pinvsqrt for out = (punctuated signed) 1 /sqrt(in) \t\t",
    "\t   add  for out = in + in2\t\t\t\t\t",
    "\t   sub  for out = in - in2\t\t\t\t\t",
    "\t   mul  for out = in * in2\t\t\t\t\t",
    "\t   div  for out = in / in2\t\t\t\t\t",
    "\t\tcartprod for out = in x in2\t\t\t\t\t",
    "\t\trequires: n=size_of_in, fastest dimension in output\t",
    "\t\tsignum for out[i] = in[i] for i< isig  and\t\t\t",
    "\t\t\t\t= -in[i] for i>= isig\t\t\t",
    "\t\trequires: isig=point where signum function acts\t\t",
    " Seismic operations:\t\t\t\t\t\t\t",
    "\t   sloth   for  out =  1/in^2\t\tSloth from velocity in\t",
    "\t   slowp   for  out =  1/in - 1/in2\tSlowness perturbation\t",
    "\t   slothp  for  out =  1/in^2 - 1/in2^2   Sloth perturbation\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " op=sqrt takes sqrt(x) for x>=0 and -sqrt(ABS(x)) for x<0 (signed sqrt)",
    "\t\t\t\t\t\t\t\t\t",
    " op=pinv takes y=1/x for x!=0,  if x=0 then y=0. (punctuated inverse)\t",
    "\t\t\t\t\t\t\t\t\t",
    " The seismic operations assume that in and in2 are wavespeed profiles.\t",
    " \"Slowness\" is 1/wavespeed and \"sloth\" is  1/wavespeed^2.\t\t",
    " Use \"suop\" and \"suop2\" to perform unary and binary operations on\t",
    " data in the SU (SEGY trace) format.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The options \"pinvsq\" and \"pinvsqrt\" are also useful for seismic\t",
    " computations involving converting velocity to sloth and vice versa.\t",
    "\t\t\t\t\t\t\t\t\t",
    " The option \"cartprod\" (cartesian product) requires also that the\t",
    " parameter n=size_of_in be set. This will be the fastest dimension\t",
    " of the rectangular array that is output.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The option \"signum\" causes a flip in sign for all values with index\t",
    " greater than \"isig\"\t(really -1*signum(index)).\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For file operations on SU format files, please use:  suop, suop2\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut name = String::new();
    let infp = if getparstring("in", &mut name) != 0 {
        fopen(&name, "r").unwrap_or_else(|| err!("Cannot open in=%s\n", name))
    } else { stdin() };
    let mut op = String::from("noop");
    getparstring("op", &mut op);
    let mut in2fp = None;
    let mut scalar = 0.0f32;
    let mut isscalar = false;
    if matches!(op.as_str(), "add" | "sub" | "mul" | "div" | "cartprod" | "slowp" | "slothp") {
        let mut in2 = String::new();
        if getparstring("in2", &mut in2) == 0 { err!("Must specify in2 for op=%s\n", op); }
        in2fp = fopen(&in2, "r");
        if in2fp.is_none() {
            let (value, end) = strtod(&in2);
            scalar = value as f32;
            if end == 0 { err!("Can neither open in2=%s nor convert it to a scalar", in2); }
            else if op == "div" && scalar == 0.0 { err!("Cannot divide by in2=%s", in2); }
            else { isscalar = true; }
        }
    }
    let mut out = String::new();
    let outfp = if getparstring("out", &mut out) != 0 {
        fopen(&out, "w").unwrap_or_else(|| err!("Cannot open out=%s", out))
    } else { stdout() };
    let mut scale = 1.0f32;
    let mut bias = 0.0f32;
    let mut isig = 0i32;
    if op == "scale" { getparfloat("scale", &mut scale); }
    if op == "bias" { getparfloat("bias", &mut bias); }
    if op == "signum" && getparint("isig", &mut isig) == 0 { err!("Must specify isig"); }
    if op == "cartprod" {
        let mut n = 0i32;
        if getparint("n", &mut n) == 0 { err!("Must specify n=size_of_in\n"); }
        let mut y1 = vec![0.0f32; n as usize];
        let nread = efread(&mut y1, &infp);
        if nread != n as usize { err!("Number of samples n_read=%d != n=%d", nread as i32, n); }
        let fp = in2fp.as_ref().unwrap_or_else(|| unsafe { raise(11); std::process::abort() });
        let mut x = [0.0f32];
        while efread(&mut x, fp) != 0 {
            for value in &y1 { efwrite(&[x[0] * *value], &outfp); }
        }
    } else {
        if !matches!(op.as_str(), "noop" | "neg" | "abs" | "scale" | "bias" | "signum" | "exp" | "sin" | "cos" | "log" | "sqrt" | "sqr" | "sloth" | "degrad" | "raddeg" | "pinv" | "pinvsqr" | "pinvsqrt" | "add" | "sub" | "mul" | "div" | "slowp" | "slothp") {
            err!("op=%s is not a valid operation", op);
        }
        let mut x = [0.0f32];
        let mut i1 = 0i32;
        while efread(&mut x, &infp) != 0 {
            let x = x[0];
            let y = match op.as_str() {
                "noop" => x, "neg" => -x, "abs" => (x as f64).abs() as f32,
                "scale" => scale*x, "bias" => x+bias,
                "signum" => { let y = if i1 >= isig { -x } else { x }; i1 = i1.wrapping_add(1); y },
                "exp" => (x as f64).exp() as f32, "sin" => (x as f64).sin() as f32,
                "cos" => (x as f64).cos() as f32, "log" => (x as f64).ln() as f32,
                "sqrt" => if x < 0.0 { -(-(x as f64)).sqrt() as f32 } else { (x as f64).sqrt() as f32 },
                "sqr" => x*x,
                "sloth" => if x != 0.0 { (1.0f64 / x as f64 * x as f64) as f32 } else { warn!("zero value detected in input!"); 0.0 },
                "degrad" => (x as f64 * std::f64::consts::PI / 180.0) as f32,
                "raddeg" => (x as f64 * 180.0 / std::f64::consts::PI) as f32,
                "pinv" => if x != 0.0 { 1.0f32/x } else { 0.0 },
                "pinvsqr" => if x != 0.0 { 1.0f32/(x*x) } else { 0.0 },
                "pinvsqrt" => if x < 0.0 { (-1.0/(-(x as f64)).sqrt()) as f32 } else if x > 0.0 { (1.0/(x as f64).sqrt()) as f32 } else { 0.0 },
                _ => {
                    let x2 = if isscalar && matches!(op.as_str(), "add" | "sub" | "mul" | "div") { scalar } else {
                        let fp = in2fp.as_ref().unwrap_or_else(|| unsafe { raise(11); std::process::abort() });
                        let mut x2 = [0.0f32];
                        if efread(&mut x2, fp) == 0 { err!("error or end of file reading in2 file"); }
                        x2[0]
                    };
                    match op.as_str() {
                        "add" => x+x2, "sub" => x-x2, "mul" => x*x2,
                        "div" => if !isscalar && x2 == 0.0 { warn!("zero values detected in denominator!"); warn!("outputting numerator*FLT_MAX"); x*f32::MAX } else { x/x2 },
                        "slowp" | "slothp" => {
                            if x <= 0.0 || x2 <= 0.0 { err!("zero value detected!"); }
                            if op == "slowp" { 1.0f32/x - 1.0f32/x2 } else { (1.0f32/x)*x - (1.0f32/x2)*x2 }
                        }, _ => unreachable!()
                    }
                }
            };
            efwrite(&[y], &outfp);
        }
    }
    checkpars();
    exit(CWP_Exit());
}

unsafe extern "C" { fn raise(sig: i32) -> i32; }
