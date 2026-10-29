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

//! SUGOUPILLAUDPO - primaries-only impulse response of a Goupillaud medium.
//! Port of su/main/synthetics_waveforms_testpatterns/sugoupillaudpo.c.
//! Credits (from the C source): CWP: Albena Mateeva, April 2001.

use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t  ",
    " SUGOUPILLAUDPO - calculate Primaries-Only impulse response of a lossless",
    "\t      GOUPILLAUD medium for plane waves at normal incidence\t",
    "\t\t\t\t\t\t\t\t\t  ",
    " sugoupillaudpo < stdin > stdout [optional parameters]\t\t      ",
    "\t\t\t\t\t\t\t\t\t    ",
    " Required parameters:\t\t\t\t\t\t\t     ",
    "\tnone\t\t\t\t\t\t\t\t",
    "      \t\t\t\t\t\t\t\t\t     ",
    " Optional parameters:\t\t\t\t\t\t       ",
    "\tl=1\t   source layer number; 1 <= l <= tr.ns\t\t  ",
    "\t\t      Source is located at the top of layer l.\t\t     ",
    "\tk=1\t   receiver layer number; 1 <= k\t\t\t ",
    "\t\t      Receiver is located at the top of layer k.\t    ",
    "\ttmax\t  number of output time-samples;\t\t\t",
    "\t\t      default: long enough to capture all primaries\t ",
    "\tpV=1\t  flag for vector field seismogram\t\t      ",
    "\t\t      (displacement, velocity, acceleration);\t       ",
    "\t\t      =-1 for pressure seismogram.\t\t\t  ",
    "\tverbose=0     silent operation, =1 list warnings\t\t    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Input: Reflection coefficient series:\t\t\t\t      ",
    "\t\t\t\t\t\t\t\t\t    ",
    "\t\t\t       impedance[i]-impedance[i+1]\t\t   ",
    "\t\t       r[i] = -----------------------------\t\t  ",
    "\t\t\t       impedance[i]+impedance[i+1]\t\t   ",
    "\t\t\t\t\t\t\t\t\t    ",
    "\tr[0]= surface refl. coef. (as seen from above)\t\t      ",
    "\tr[n]= refl. coef. of the deepest interface\t\t\t  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Input file is to be in SU format, i.e., binary floats with a SU header.    ",
    "\t\t\t\t\t\t\t\t\t    ",
    " Remarks:\t\t\t\t\t\t\t\t   ",
    " 1. For vector fields, a buried source produces a spike of amplitude 1      ",
    " propagating downwards and a spike of amplitude -1 propagating upwards.     ",
    " A buried pressure source produces spikes of amplitude 1 both in the up-    ",
    " and downward directions.\t\t\t\t\t\t   ",
    "    A surface source induces only a downgoing spike of amplitude 1 at the   ",
    " top of the first layer (both for vector and pressure fields).\t      ",
    " 2. The sampling interval dt in the header of the input reflectivity file   ",
    " is interpreted as a two-way traveltime thicknes of the layers. The sampling",
    " interval of the output seismogram is the same as that of the input file.   ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    if gettr(&mut tr) == 0 { err!("Can't get reflectivity"); }
    let n = tr.ns as i32 - 1;
    let mut k = 1i32;
    let mut l = 1i32;
    getparint("k", &mut k);
    getparint("l", &mut l);
    let skl = if k > l { 1 } else if k < l { -1 } else { 0 };
    let n1 = if k <= l { k } else { l };
    let n2 = if k >= l { k } else { l };
    let mut tmax = 0i32;
    if getparint("tmax", &mut tmax) == 0 {
        let extent = n1.wrapping_sub(1)
            .max(n.wrapping_add(1).wrapping_sub(n2))
            .max(n.wrapping_add(1).wrapping_sub(n1));
        tmax = n2.wrapping_sub(n1).wrapping_add(2)
            .wrapping_add(2i32.wrapping_mul(extent)) / 2;
    }
    let mut p_v = 1i32;
    let mut verbose = 0i32;
    getparint("pV", &mut p_v);
    getparint("verbose", &mut verbose);
    if n < 0 { err!("The number of subsurface interfaces n cannot be negative."); }
    if n == 0 { warn!("WARNING: model without subsurface reflectors!"); }
    if k < 1 { err!("Receiver layer k must be >=1 (k=1 corresponds to surface seismogram)."); }
    if l < 1 { err!("Source layer l must be >= 1 (l=1 corresponds to a surface source)."); }
    if l > n + 1 { err!("The current version of the program requires l<=n+1."); }
    if tmax < 0 { err!("The number of the output time samples tmax cannot be negative."); }
    if p_v != 1 && p_v != -1 { err!("The field-type flag pV should be either 1 or -1."); }
    if verbose != 0 {
        warn!("Source layer l=%d", l);
        warn!("Number of layers n=%d", n);
        warn!("Receiver layer k=%d", k);
        warn!("Output time samples tmax=%d", tmax);
        warn!("Field type flag pV=%d", p_v);
    }
    let mut r = vec![0.0f32; (n + 1) as usize];
    for i in 0..=n as usize {
        r[i] = tr.data[i];
        if r[i] > 1.0 || r[i] < -1.0 { err!("Invalid reflection coefficient encountered."); }
    }
    for v in &mut r { *v *= p_v as f32; }
    if (n2 - n1) as i64 >= 2 * tmax as i64 { err!("Seismogram length too small -- cannot observe any signal."); }
    // These inputs would access beyond the C reflectivity or trace buffers.
    if n2 > n + 1 { err!("Receiver layer k exceeds the reflectivity model."); }
    if tmax as usize > tr.data.len() { err!("The number of output time samples exceeds the trace buffer."); }
    let mut x = vec![0.0f32; (2 * tmax) as usize];
    let mut transm1 = 1.0f32;
    for i in n1..n2 { transm1 *= 1.0f32 + skl as f32 * r[i as usize]; }
    x[(n2 - n1) as usize] = if p_v == -1 || l == 1 { transm1 } else { skl as f32 * transm1 };
    let mut effsd = 1.0f32;
    if l != 1 { effsd = 1.0f32 + p_v as f32 * r[(l - 1) as usize]; }
    let rmax = n.min(tmax - 1 + (n2 + n1 - 1) / 2);
    let mut transm2 = 1.0f32;
    for i in n2..=rmax {
        x[(2 * (i + 1) - n2 - n1) as usize] = effsd * transm1 * r[i as usize] * transm2;
        transm2 *= 1.0f32 - r[i as usize] * r[i as usize];
    }
    if l != 1 {
        transm2 = 1.0;
        let mut i = n1 - 1;
        while i >= 0 && n2 - n1 + 2 * (n1 - 1 - i) < 2 * tmax {
            x[(n2 - n1 + 2 * (n1 - 1 - i)) as usize] += p_v as f32 * transm1 * r[i as usize] * transm2;
            transm2 *= 1.0f32 - r[i as usize] * r[i as usize];
            i -= 1;
        }
    }
    if 2 * ((k - l) / 2) == k - l {
        for i in 0..tmax as usize { tr.data[i] = x[2 * i]; }
        tr.delrt = 0;
    } else {
        for i in 0..tmax as usize { tr.data[i] = x[2 * i + 1]; }
        tr.delrt = (tr.dt / 2000) as i16;
    }
    tr.ns = tmax as u16;
    tr.trid = 1;
    puttr(&mut tr);
    exit(CWP_Exit());
}
