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

//! SUWAVEFORM - generate a seismic wavelet.
//! Port of su/main/synthetics_waveforms_testpatterns/suwaveform.c (Revision 1.9).
//! Author: Nils Maercklin; credits: CWP and the authors of waveforms.c.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{nint, PI};
use oxwave_waveforms::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUWAVEFORM - generate a seismic wavelet\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " suwaveform <stdin >stdout [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t  \t",
    "\tone of the optional parameters listed below\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t  \t",
    "\ttype=akb\twavelet type\t\t\t\t\t",
    "\t\t   akb:\tAKB wavelet defined by max frequency fpeak\t",
    "\t\t   berlage: Berlage wavelet\t\t\t\t",
    "\t\t   gauss:   Gaussian wavelet defined by frequency fpeak\t",
    "\t\t   gaussd:  Gaussian first derivative wavelet\t\t",
    "\t\t   ricker1: Ricker wavelet defined by frequency fpeak\t",
    "\t\t   ricker2: Ricker wavelet defined by half and period\t",
    "\t\t   spike:   spike wavelet, shifted by time tspike\t",
    "\t\t   unit:\tunit wavelet, i.e. amplitude = 1 = const.",
    "\t\t\t\t\t\t\t\t\t",
    "\tdt=0.004\ttime sampling interval in seconds\t\t",
    "\tns=\t\tif set, number of samples in  output trace\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tfpeak=20.0\tpeak frequency of a Berlage, Ricker, or Gaussian,",
    "\t\t   and maximum frequency of an AKB wavelet in Hz\t",
    "\t\t\t\t\t\t\t\t\t",
    "\thalf=1/fpeak   Ricker wavelet \"ricker2\": half-length\t\t",
    "\tperiod=c*half  Ricker wavelet \"ricker2\": period (c=sqrt(6)/pi)",
    "\tdistort=0.0\tRicker wavelet \"ricker2\": distortion factor\t",
    "\tdecay=4*fpeak  Berlage wavelet: exponential decay factor in 1/sec",
    "\ttn=2\t   Berlage wavelet: time exponent\t\t\t",
    "\tipa=-90\tBerlage wavelet: initial phase angle in degrees\t\t",
    "\ttspike=0.0\tSpike wavelet: time at spike in seconds\t\t",
    "\tverbose=0\t1: echo output wavelet length\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    "\tIf ns is not defined, the program determines the trace length\t",
    "\tdepending on the dominant signal period.\t\t\t   ",
    "\t\t\t\t\t\t\t\t\t",
    "\tThe Ricker wavelet \"ricker1\" and the Gaussian wavelet \"gauss\"  ",
    "\tare zero-phase. For these two wavelets, the trace header word\t",
    "\tdelrt is set such that the peak amplitude is at t=0 seconds.\t",
    "\tIf this is not acceptable, use \"sushw key=delrt a=0\".\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tThe Ricker wavelets can be defined either by the peak frequency\t",
    "\tfpeak (\"ricker1\") or by its half-length, the period, and a\t",
    "\tdistortion factor (\"ricker2\"). \"ricker\" is an acceptable\t",
    "\talias for \"ricker1\".\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tThe Berlage wavelet is defined by the peak frequency fpeak, a time ",
    "\ttime exponent tn describing the wavelet shape at its beginning,\t",
    "\tand an exponential decay factor describing the amplitude decay\t",
    "\ttowards later times. The parameters tn and decay are non-negative, ",
    "\treal numbers; tn is typically a small integer number and decay a   ",
    "\tmultiple of the dominant signal period 1/fpeak. Additionally, an   ",
    "\tinitial phase angle can be given; use -90 or 90 degrees for\t",
    "\tzero-amplitude at the beginning.\t\t\t\t   ",
    "\t\t\t\t\t\t\t\t\t",
    "\tFor an AKB wavelet, fpeak is the maximum frequency; the peak\t",
    "\tfrequency is about 1/3 of the fpeak value.\t\t\t ",
    "\t\t\t\t\t\t\t\t\t",
    "\tThe output wavelet can be normalized or scaled with \"sugain\".\t",
    "\tUse \"suvibro\" to generate a Vibroseis sweep.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " A normalized, zero-phase Ricker wavelet with a peak frequency\t\t",
    " of 15 Hz is generated and convolved with a spike dataset:\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsuwaveform type=ricker1 fpeak=15 | sugain pbal=1 > wavelet.su\t",
    "\tsuplane npl=1 | suconv sufile=wavelet.su | suxwigb\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Gaussian and derivatives of Gaussians:\t\t\t\t",
    " Use \"sudgwaveform\" to generate these\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Technical references:\t\t\t\t\t\t\t",
    " Aldridge, D. F., 1990, The Berlage wavelet: Geophysics, 55,1508-1511.\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    "\tThis program does not check for aliasing.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut dt = 0.004f32;
    getparfloat("dt", &mut dt);
    if dt <= 0.0 { err!("dt=%g must be positive", dt); }
    let mut fpeak = 20.0f32;
    getparfloat("fpeak", &mut fpeak);
    if fpeak <= 0.0 { err!("fpeak=%g must be positive", fpeak); }
    let mut half = 0.0f32;
    if getparfloat("half", &mut half) == 0 { half = (1.0 / fpeak as f64) as f32; }
    if half <= 0.0 { err!("half=%g must be positive", half); }
    let mut period = 0.0f32;
    if getparfloat("period", &mut period) == 0 { period = (half as f64 * 0.77969680123) as f32; }
    if period <= 0.0 { err!("period=%g must be positive", period); }
    let mut tspike = 0.0f32;
    getparfloat("tspike", &mut tspike);
    if tspike < 0.0 { err!("tspike=%g must be non-negative", tspike); }
    let mut ampl = 1.0f32;
    getparfloat("ampl", &mut ampl);
    let mut distort = 0.0f32;
    getparfloat("distort", &mut distort);
    let mut tn = 2.0f32;
    getparfloat("tn", &mut tn);
    if tn < 0.0 { err!("tn=%g must be non-negative", tn); }
    let mut decay = 0.0f32;
    if getparfloat("decay", &mut decay) == 0 { decay = (4.0 * fpeak as f64) as f32; }
    if decay < 0.0 { err!("decay=%g must be non-negative", decay); }
    let mut ipa = -90.0f32;
    getparfloat("ipa", &mut ipa);
    ipa = (ipa as f64 * (PI / 180.0)) as f32;
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);
    let mut wtype = String::new();
    if getparstring("type", &mut wtype) == 0 { wtype = "akb".into(); }
    if wtype == "ricker" { wtype = "ricker1".into(); }
    let mut nt = 0i32;
    if getparint("ns", &mut nt) == 0 && getparint("nt", &mut nt) == 0 {
        nt = match wtype.as_str() {
            "berlage" if decay != 0.0 => nint((fpeak as f64 * (8.0 + 2.0 * tn as f64) / decay as f64).floor() / (dt * fpeak) as f64).wrapping_add(1),
            "berlage" | "ricker1" | "gauss" | "gaussd" => nint(2.0 / (dt * fpeak) as f64).wrapping_add(1),
            "ricker2" => nint(2.0 * half as f64 / dt as f64).wrapping_add(1),
            "akb" => nint(4.0 / (fpeak * dt) as f64).wrapping_add(1),
            "spike" => nint((tspike / dt) as f64).wrapping_add(5),
            "unit" => nint((half / dt) as f64).wrapping_add(1),
            _ => err!("unknown wavelet type=%s", wtype),
        };
    }
    if nt < 1 { warn!("ns=%d too small, using ns=1", nt); nt = 1; }
    if nt > SU_NFLTS as i32 { warn!("trace too long, using maximum ns=%d", SU_NFLTS as i32); nt = SU_NFLTS as i32; }
    let hlw = nint((half / dt) as f64).max(1) as usize;
    let spike_index = nint((tspike / dt) as f64) as usize;
    // C may write beyond nt; allocate room while copying only nt samples.
    let len = match wtype.as_str() {
        "ricker2" => (nt as usize).max(hlw.saturating_mul(2).saturating_sub(1)),
        "spike" => (nt as usize).max(spike_index.saturating_add(1)),
        _ => nt as usize,
    };
    let mut wavelet = vec![0.0f32; len];
    match wtype.as_str() {
        "berlage" => berlage_wavelet(dt, fpeak, ampl, tn, decay, ipa, &mut wavelet),
        "ricker1" => { ricker1_wavelet(dt, fpeak, &mut wavelet); tr.delrt = (-nint(1000.0 / fpeak as f64)) as i16; }
        "ricker2" => { ricker2_wavelet(hlw, dt, period, ampl, distort, &mut wavelet); if distort == 0.0 { tr.delrt = (-nint(1000.0 * half as f64)) as i16; } }
        "gauss" => { gaussian_wavelet(dt, fpeak, &mut wavelet); tr.delrt = (-nint(1000.0 / fpeak as f64)) as i16; }
        "gaussd" => gaussderiv_wavelet(dt, fpeak, &mut wavelet),
        "akb" => akb_wavelet(dt, fpeak, &mut wavelet),
        "spike" => spike_wavelet(spike_index, &mut wavelet),
        "unit" => unit_wavelet(&mut wavelet),
        _ => err!("unknown wavelet type=%s", wtype),
    }
    tr.data[..nt as usize].copy_from_slice(&wavelet[..nt as usize]);
    tr.tracl = 1;
    tr.dt = nint(dt as f64 * 1000000.0) as u16;
    tr.ns = nt as u16;
    tr.trid = 1;
    puttr(&mut tr);
    if verbose != 0 { warn!("%s wavelet with %d samples (%g seconds)", wtype, nt, dt * (nt - 1) as f32); }
    exit(CWP_Exit());
}
