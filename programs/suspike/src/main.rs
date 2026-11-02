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

//! SUSPIKE - make a small spike data set.
//! Port of su/main/synthetics_waveforms_testpatterns/suspike.c (Revision 1.15).
//! Credits: CWP: Shuki Ronen, Chris Liner.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUSPIKE - make a small spike data set \t\t\t",
    "\t\t\t\t\t\t\t\t",
    " suspike [optional parameters] > out_data_file  \t\t",
    "\t\t\t\t\t\t\t\t",
    " Creates a common offset su data file with up to four spikes\t",
    " for impulse response studies\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\tnt=64 \t\tnumber of time samples\t\t\t",
    "\tntr=32\t\tnumber of traces\t\t\t",
    " \tdt=0.004 \ttime sampling interval in seconds\t",
    " \toffset=400 \toffset\t\t\t\t\t",
    "\tnspk=4\t\tnumber of spikes\t\t\t",
    "\tix1= ntr/4\ttrace number (from left) for spike #1\t",
    "\tit1= nt/4 \ttime sample to spike #1\t\t\t",
    "\tix2 = ntr/4\ttrace for spike #2\t\t\t",
    "\tit2 = 3*nt/4 \ttime for spike #2\t\t\t",
    "\tix3 = 3*ntr/4;\ttrace for spike #3\t\t\t",
    "\tit3 = nt/4;\ttime for spike #3\t\t\t",
    "\tix4 = 3*ntr/4;\ttrace for spike #4\t\t\t",
    "\tit4 = 3*nt/4;\ttime for spike #4\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut nt = 64;
    getparint("nt", &mut nt);
    CHECK_NT("nt", nt);
    tr.ns = nt as u16;
    let mut ntr = 32;
    getparint("ntr", &mut ntr);
    let mut dt: f32 = 0.004;
    getparfloat("dt", &mut dt);
    tr.dt = (dt * 1000000.0) as u32 as u16;
    let mut offset: f32 = 400.0;
    getparfloat("offset", &mut offset);
    tr.offset = offset as i32;
    let mut nspk = 4;
    getparint("nspk", &mut nspk);
    let mut ix1 = ntr / 4;
    getparint("ix1", &mut ix1);
    let mut it1 = nt / 4;
    getparint("it1", &mut it1);
    let mut ix2 = ntr / 4;
    getparint("ix2", &mut ix2);
    let mut it2 = 3i32.wrapping_mul(nt) / 4;
    getparint("it2", &mut it2);
    let mut ix3 = 3i32.wrapping_mul(ntr) / 4;
    getparint("ix3", &mut ix3);
    let mut it3 = nt / 4;
    getparint("it3", &mut it3);
    let mut ix4 = 3i32.wrapping_mul(ntr) / 4;
    getparint("ix4", &mut ix4);
    let mut it4 = 3i32.wrapping_mul(nt) / 4;
    getparint("it4", &mut it4);

    for itr in 0..ntr {
        // A negative count makes the C memset exceed the trace object.
        if nt < 0 {
            err!("nt=%d must be nonnegative", nt);
        }
        tr.data[..nt as usize].fill(0.0);
        if itr == ix1.wrapping_sub(1) {
            set_spike(&mut tr, it1);
        }
        if nspk > 1 && itr == ix2.wrapping_sub(1) {
            set_spike(&mut tr, it2);
        }
        if nspk > 2 && itr == ix3.wrapping_sub(1) {
            set_spike(&mut tr, it3);
        }
        if nspk > 3 && itr == ix4.wrapping_sub(1) {
            set_spike(&mut tr, it4);
        }
        tr.tracl = itr + 1;
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

// Reproduce the C assignment relative to tr.data, including positions beyond
// nt and writes into the preceding header (e.g. the default it1=0 for nt=1).
// Discard writes outside the whole C segy object. In UB round 1, guarding
// it=65536 and it=-60 in the original C left stdout, stderr and status unchanged;
// those cases must therefore generate the original trace instead of an error.
fn set_spike(tr: &mut Trace, it: i32) {
    let index = it.wrapping_sub(1);
    if index >= 0 && (index as usize) < SU_NFLTS {
        tr.data[index as usize] = 1.0;
    } else if index < 0 && index >= -(HDRBYTES as i32 / 4) {
        let offset = (HDRBYTES as i32 + index * 4) as usize;
        tr.hdr.as_bytes_mut()[offset..offset + 4].copy_from_slice(&1.0f32.to_ne_bytes());
    }
}
