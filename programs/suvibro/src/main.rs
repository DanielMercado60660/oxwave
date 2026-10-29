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

//! SUVIBRO - generate a Vibroseis sweep.
//! Port of su/main/synthetics_waveforms_testpatterns/suvibro.c (Revision 1.15).
//! Credits: CWP: Michel Dietrich; rewrite: Tagir Galikeev, CWP, 7 October 1994.
use oxwave_core::prelude::*;
use std::f64::consts::PI;
use oxwave_cwp::cwp::nint;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " SUVIBRO - Generates a Vibroseis sweep (linear, linear-segment,",
    "\t\t\tdB per Octave, dB per Hertz, T-power)\t",
    "\t\t\t\t\t\t\t\t",
    " suvibro [optional parameters] > out_data_file\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t",
    " dt=0.004\t\ttime sampling interval\t\t\t",
    " sweep=1\t  \tlinear sweep\t\t\t  \t",
    " \t\t  \t=2 linear-segment\t\t\t",
    " \t\t  \t=3 decibel per octave\t \t\t",
    " \t\t  \t=4 decibel per hertz\t  \t\t",
    " \t\t  \t=5 t-power\t\t\t\t",
    " swconst=0.0\t\tsweep constant (see note)\t\t",
    " f1=10.0\t\tsweep frequency at start\t\t",
    " f2=60.0\t\tsweep frequency at end\t\t\t",
    " tv=10.0\t\tsweep length\t\t\t\t",
    " phz=0.0\t\tinitial phase (radians=1 default)\t",
    " radians=1\t\t=0 degrees\t\t\t\t",
    " fseg=10.0,60.0\tfrequency segments (see notes)\t\t",
    " tseg=0.0,10.0\t\ttime segments (see notes)\t\t",
    " t1=1.0\t\tlength of taper at start (see notes)\t",
    " t2=1.0\t\tlength of taper at end (see notes)\t",
    " taper=1\t\tlinear\t\t\t\t\t",
    "\t\t  \t=2 sine\t\t\t\t\t",
    "\t\t\t=3 cosine\t\t\t\t",
    "\t\t\t=4 gaussian (+/-3.8)\t\t\t",
    "\t\t\t=5 gaussian (+/-2.0)\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " The default tapers are linear envelopes. To eliminate the\t",
    " taper, choose t1=t2=0.0.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " \"swconst\" is active only with nonlinear sweeps, i.e. when\t",
    " sweep=3,4,5.\t\t\t\t\t\t\t",
    " \"tseg\" and \"fseg\" arrays are used when only sweep=2\t",
    "\t\t\t\t\t\t\t\t",
    " Sweep is a modulated cosine function.\t\t\t\t",
];

fn main() {
    let argv = args(); initargs(&argv); requestdoc(0, SDOC);
    let (mut sweep, mut swconst, mut f1, mut f2, mut tv) = (1, 0.0f32, 10.0f32, 60.0f32, 10.0f32);
    let (mut radians, mut phz, mut t1, mut t2, mut taper, mut dt) = (1, 0.0f32, 1.0f32, 1.0f32, 1, 0.004f32);
    getparint("sweep", &mut sweep); getparfloat("swconst", &mut swconst);
    getparfloat("f1", &mut f1); getparfloat("f2", &mut f2); getparfloat("tv", &mut tv);
    getparint("radians", &mut radians); getparfloat("phz", &mut phz);
    getparfloat("t1", &mut t1); getparfloat("t2", &mut t2); getparint("taper", &mut taper); getparfloat("dt", &mut dt);
    let mut tr = Trace::new();
    tr.dt = nint(1000000.0 * dt as f64) as u16;
    let nt = (tv / dt + 1.0f32) as i32;
    if nt >= SU_NFLTS as i32 { err!("nt=tv/dt=%d -- too big", nt); }
    if t1 + t2 > tv { err!("sum of tapers t1=%f, t2=%f exceeds tv=%f", t1, t2, tv); }
    if radians == 0 { phz = (phz as f64 * (2.0 * PI) / 360.0) as f32; }
    let mut freq = Vec::new(); let mut time = Vec::new();
    if sweep == 2 {
        let mut ntseg = countparval("tseg"); if ntseg == 0 { ntseg = 2; }
        let mut tseg = vec![0.0f32; ntseg as usize];
        if getparfloat("tseg", &mut tseg) == 0 { tseg[0] = 0.0; tseg[1] = tv; }
        let mut nfseg = countparval("fseg"); if nfseg == 0 { nfseg = 2; }
        if nfseg != ntseg { err!("number of tseg and fseg must be equal"); }
        freq = vec![0.0f32; nfseg as usize];
        if getparfloat("fseg", &mut freq) == 0 { freq[0] = f1; freq[1] = f2; }
        if ntseg > 2 { for i in 1..ntseg as usize { if tseg[i] <= tseg[i-1] { err!("tseg must increase monotonically"); } } }
        let isegm = ntseg as usize - 1; tv = tseg[isegm];
        for i in 0..isegm { time.push(tseg[i+1] - tseg[i]); }
    }
    tr.ns = nt as u16; tr.tracl = 1; tr.sfs = f1 as i32 as i16; tr.sfe = f2 as i32 as i16;
    tr.slen = (1000.0 * tv as f64) as i32 as i16; tr.styp = sweep as i16;
    // tseg replaces tv after the original size check. Keep tr.ns based on
    // the original tv, and validate the actual segment ranges below using
    // this recomputed count. Exactly SU_NFLTS generated samples still fit.
    let n = (tv / dt + 1.0f32) as i32;
    if sweep == 1 {
        let rate = (f2 - f1) / tv;
        for i in 0..n { let t = i as f32 * dt; tr.data[i as usize] = (2.0 * PI * (f1 + rate / 2.0 * t) as f64 * t as f64 + phz as f64).cos() as f32; }
    }
    if sweep == 2 {
        let (mut m, mut k, mut phase) = (0i32, 0i32, 0.0f32);
        for i in 0..time.len() {
            k = k.wrapping_add(m);
            m = (n as f32 / tv * time[i]) as i32;
            // A negative first tseg can make the segment span exceed n.
            // Check the exclusive end in i64 before writing any samples.
            // Validate actual writes: n alone can exceed capacity when no
            // segment generates samples, which is safe in the reference.
            let segment_end = k as i64 + m as i64;
            if m > 0 && (k < 0 || segment_end > SU_NFLTS as i64) {
                if n > SU_NFLTS as i32 { err!("nt=tv/dt=%d -- too big", n); }
                err!("segment samples exceed trace buffer");
            }
            let aa = (2.0 * PI * freq[i] as f64 * dt as f64) as f32;
            let ab = (PI * (freq[i+1] - freq[i]) as f64 * (dt * dt) as f64 / time[i] as f64) as f32;
            for j in 0..m { let phi = j as f32 * (ab * j as f32 + aa) + phase + phz; tr.data[(j+k) as usize] = (phi as f64).cos() as f32; }
            phase = (ab * m as f32 + aa) * m as f32 + phase;
        }
    }
    if sweep == 3 {
        if swconst == -6.0 { swconst = -5.999; }
        let swcon = swconst / 6.0 + 1.0; let s = (swcon + 1.0) / swcon;
        let k1 = (f1 as f64).powf(swcon as f64) as f32;
        let k2 = ((f2 as f64).powf(swcon as f64) as f32 - (f1 as f64).powf(swcon as f64) as f32) / tv;
        for i in 0..n { let t = i as f32 * dt; tr.data[i as usize] = (2.0 * PI / (s*k2) as f64 * ((k1+k2*t) as f64).powf(s as f64) + phz as f64).cos() as f32; }
    }
    if sweep == 4 || sweep == 5 {
        let k1 = (20.0 / (swconst as f64 * 10.0f64.ln())) as f32;
        let k2 = (((swconst as f64 * 10.0f64.ln() * (f2-f1) as f64 / 20.0).exp() - 1.0) as f32) / tv;
        for i in 0..n {
            let t = i as f32 * dt;
            if swconst == 0.0 { tr.data[i as usize] = f1; }
            else { let v = 1.0 + (t*k2) as f64; tr.data[i as usize] = (2.0*PI*((f1*t) as f64 + (k1/k2) as f64 * (v*v.ln()-v)) + phz as f64).cos() as f32; }
        }
    }
    if t1 != 0.0 || t2 != 0.0 { sweep_taper(&mut tr, t1, t2, taper, tv, dt); }
    puttr(&mut tr); exit(CWP_Exit());
}
fn sweep_taper(tr: &mut Trace, t1: f32, t2: f32, tap_type: i32, tv: f32, dt: f32) {
    let nt = (tv/dt + 1.0f32) as i32;
    let nt1 = (t1/dt + 1.0f32) as i32;
    let nt2 = (t2/dt + 1.0f32) as i32;
    for (count, end) in [(nt1, false), (nt2, true)] {
        if count > 1 { for i in 0..count {
            let f = i as f32 / count as f32;
            let env = match tap_type as i8 {
                1 => f, 2 => (PI*f as f64/2.0).sin() as f32,
                3 => (0.5*(1.0-(PI*f as f64).cos())) as f32,
                4 => { let x = (3.8090232 * (1.0f32-f) as f64) as f32; (-(x*x) as f64).exp() as f32 },
                5 => { let x = (2.0 * (1.0f32-f) as f64) as f32; (-(x*x) as f64).exp() as f32 },
                _ => { err!(" taper ?!"); }
            };
            let idx = if end { nt.wrapping_sub(i) } else { i };
            // C starts the end taper at nt, which is valid only while that
            // index (and every following taper index) fits the static buffer.
            let Some(sample) = tr.data.get_mut(idx as usize) else {
                err!("taper sample index=%d outside trace buffer", idx);
            };
            *sample *= env;
        } }
    }
}
