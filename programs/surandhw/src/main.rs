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

// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SURANDHW - set a header word to a random variable.
//! Port of su/main/headers/surandhw.c (Revision 1.3).
//! Credits (from the C source):
//! Eric Verschuur Delft University of Technology
//! eric@delphi.tn.tudelft.nl
//!
//! Based on sushw written by:
//! SEP: Einar Kjartannson
//! CWP: Jack K. Cohen
//!
//! Caveat: All constants are cast to doubles.

use oxwave_core::prelude::*;
use oxwave_f_frannor::*;
use oxwave_f_franuni::*;
use std::time::{SystemTime, UNIX_EPOCH};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " surandhw - set header word to random variable \t\t",
    " \t\t\t\t\t\t\t\t",
    " surandhw <stdin >stdout key=tstat a=0 min=0 max=1\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " \tnone (no op)\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " \tkey=tstat\theader key word to set\t\t\t",
    " \ta=0\t\t=1 flag to add original value to final key",
    " \tnoise=gauss\tnoise probability distribution\t\t",
    " \t\t\t=flat for uniform; default Gaussian\t",
    " \tseed=from_clock\trandom number seed (integer)\t\t",
    " \tmin=0\t\tminimum random number\t\t\t",
    " \tmax=1\t\tmaximum radnom number\t\t \t",
    " \t\t\t\t\t\t\t\t",
    " NOTES:\t\t\t\t\t\t\t",
    " The value of header word key is computed using the formula:\t",
    " \tval(key) = a * val(key) + rand\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t",
    "  \tsurandhw <indata key=tstat a=0 min=0 max=10  > outdata\t",
    " \t\t\t\t\t\t\t\t",
];

// Noise probability distributions.
const GAUSS: i32 = 0;
const FLAT: i32 = 1;

fn setrandval(type_: &str, val: &mut Value, a: f64, noise: f64) {
    match type_.as_bytes()[0] {
        b's' => err!("can't set char header word"),
        // The C build converts through a 32-bit register before storing
        // the low 16 bits; casting directly to a Rust short saturates early.
        b'h' => val.set_h((a * val.h() as f64 + noise) as i32 as i16),
        b'u' => val.set_u((a * val.u() as f64 + noise) as u32 as u16),
        b'l' => val.set_l((a * val.l() as f64 + noise) as i64),
        b'v' => val.set_v((a * val.v() as f64 + noise) as u64),
        b'i' => val.set_i((a * val.i() as f64 + noise) as i32),
        b'p' => val.set_p((a * val.p() as f64 + noise) as u32),
        b'f' => val.set_f((a * val.f() as f64 + noise) as f32),
        b'd' => {
            val.set_d(a * val.d() + noise);
            err!("unknown type %s", type_);
        }
        _ => err!("unknown type %s", type_),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 {
        key = "tstat".into();
    }
    let mut a = 1.0f64;
    let mut min = 0.0f64;
    let mut max = 1.0f64;
    getpardouble("a", &mut a);
    getpardouble("min", &mut min);
    getpardouble("max", &mut max);
    let mut stype = String::new();
    if getparstring("noise", &mut stype) == 0 {
        stype = "gauss".into();
    }
    let itype = if stype == "gauss" {
        GAUSS
    } else if stype == "flat" {
        FLAT
    } else {
        err!("noise=\"%s\", must be gauss or flat", stype);
    };
    let mut seed = 0u32;
    if getparuint("seed", &mut seed) == 0 {
        seed = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|x| x.as_secs() as u32)
            .unwrap_or(u32::MAX);
        if seed == u32::MAX {
            err!("time() failed to set seed");
        }
    }
    checkpars();

    // C passes the unsigned seed to functions taking a signed int.
    if itype == GAUSS {
        srannor(seed as i32);
    }
    sranuni(seed as i32);

    let type_ = hdtype(&key);
    let index = getindex(&key);
    while gettr(&mut tr) != 0 {
        // Get new noise value, including when the endpoints are equal.
        let noise = match itype {
            // C evaluates 1 + frannor() in float before promoting to double.
            GAUSS => min + ((max - min) * (1.0f32 + frannor()) as f64 / 2.0),
            FLAT => (max - min) * franuni() as f64 + min,
            _ => err!("%d: mysterious itype = %d", 126, itype),
        };

        // Get desired key from header and add noise.
        let mut val = Value::default();
        gethval(&tr, index, &mut val);
        setrandval(type_, &mut val, a, noise);
        puthval(&mut tr, index, &val);

        // Write updated trace.
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
