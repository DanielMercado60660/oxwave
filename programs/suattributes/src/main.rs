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
//! SUATTRIBUTES - instantaneous trace attributes, port of suattributes.c.
//! Credits: Jack Cohen, John Stockwell, Agung Wiyono, Kylee Brown and Steven Rennolet.

use oxwave_core::prelude::*;
use oxwave_hilbert::hilbert;
const PI: f64 = 3.141592653589793;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUATTRIBUTES - instantaneous trace ATTRIBUTES \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suattributes <stdin >stdout mode=amp\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t\t",
    " \tmode=amp\toutput flag \t\t\t\t\t",
    " \t       \t\t=amp envelope traces\t\t\t\t",
    " \t       \t\t=phase phase traces\t\t\t\t",
    " \t       \t\t=uphase unwrapped phase traces\t\t\t",
    " \t       \t\t=freq frequency traces\t\t\t\t",
    "\t\t\t=freqw Frequency Weighted Envelope\t\t",
    "\t\t\t=thin  Thin-Bed (inst. freq - average freq)\t",
    "\t\t\t=bandwith Instantaneous bandwidth\t\t",
    "\t\t\t=normamp Normalized Phase (Cosine Phase)\t",
    " \t       \t\t=fdenv 1st envelope traces derivative\t\t",
    " \t       \t\t=sdenv 2nd envelope traces derivative\t\t",
    " \t       \t\t=q Ins. Q Factor\t\t\t\t",
    "\tunwrap=\t\tdefault unwrap=0 for mode=phase\t\t\t",
    " \t\t\tdefault unwrap=1 for freq, uphase, freqw, Q\t",
    " \t\t\tdphase_min=PI/unwrap\t\t\t\t",
    "\twint=\t\twindowing for freqw\t\t\t\t",
    "\t\t\twindowing for thin\t\t\t\t",
    "\t\t\tdefault=1 \t\t\t\t\t",
    " \t\t\to--------o--------o\t\t\t\t",
    " \t\t\tdata-1\tdata\tdata+1\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program performs complex trace attribute analysis. The first three",
    " attributes, amp,phase,freq are the classical Taner, Kohler, and\t",
    " Sheriff, 1979 attributes.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The unwrap parameter is active only for mode=freq and mode=phase. The\t",
    " quantity dphase_min is the minimum change in the phase angle taken to be",
    " the result of phase wrapping, rather than natural phase variation in the",
    " data. Setting unwrap=0 turns off phase-unwrapping altogether. Choosing",
    " unwrap > 1 makes the unwrapping function more sensitive to phase changes.",
    " Setting unwrap > 1 may be necessary to resolve higher frequencies in\t",
    " data (or sample data more finely). The phase unwrapping is crude. The ",
    " differentiation needed to compute the instantaneous frequency\t\t",
    " freq(t)= d(phase)/dt is a simple centered difference.\t\t\t",
    "\t \t\t\t\t\t       \t\t\t",
    " The mode=uphase generates uwrapped phase traces by integrating the\t",
    " instantaneous amplitude traces.\t\t       \t\t\t",
    "\t \t\t\t\t\t       \t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " suvibro f1=10 f2=50 t1=0 t2=0 tv=1 | suattributes2 mode=amp | ...\t",
    " suvibro f1=10 f2=50 t1=0 t2=0 tv=1 | suattributes2 mode=phase | ...\t",
    " suvibro f1=10 f2=50 t1=0 t2=0 tv=1 | suattributes2 mode=freq | ...\t",
    " suplane | suattributes mode=... | supswigb |...       \t\t",
];

fn differentiate(f: &mut [f32], h: f32) {
    let n = f.len();
    if n == 1 { f[0] = (0.0 - f[0]) / h; return; }
    let mut temp = vec![0.0; n];
    let h2 = 2.0 * h;
    temp[0] = (f[1] - f[0]) / h;
    for i in 1..n-1 { temp[i] = (f[i+1] - f[i-1]) / h2; }
    temp[n-1] = (f[n-1] - f[n-2]) / h;
    f.copy_from_slice(&temp);
}
fn unwrap_phase(phase: &mut [f32], w: f32) {
    if w == 0.0 { err!("wrapping parameter is zero"); }
    let pibyw = (PI / w as f64) as f32;
    let n = phase.len();
    let mut dphase = vec![0.0; n];
    let mut temp = vec![0.0; n];
    temp[0] = phase[0];
    for i in 1..n {
        dphase[i] = (phase[i] - phase[i-1]).abs();
        if (dphase[i] - dphase[i-1]).abs() >= pibyw { dphase[i] = dphase[i-1]; }
        temp[i] = temp[i-1] + dphase[i];
    }
    phase.copy_from_slice(&temp);
}
fn twindow(data: &mut [f32], wtime: i32) {
    let nt = data.len() as i32;
    // C reads outside its allocation for invalid windows; pad those reads.
    let offset = wtime.abs() as usize + 1;
    let mut padded = vec![0.0; nt as usize + 2*offset+2];
    padded[offset..offset+nt as usize].copy_from_slice(data);
    let get = |i: i32| -> f32 { padded[(i + offset as i32) as usize] };
    let mut temp = vec![0.0; padded.len()];
    let mut nwin = 2*wtime+1;
    let mut sum = 0.0f32;
    for i in 0..wtime+1 { sum += get(i); }
    temp[offset] = sum / nwin as f32;
    for i in 1..wtime { sum += get(i+wtime); nwin += 1; temp[(i+offset as i32) as usize] = sum/nwin as f32; }
    for i in wtime..nt-wtime { sum += get(i+wtime); sum -= get(i-wtime); temp[(i+offset as i32) as usize] = sum/nwin as f32; }
    for i in nt-wtime..nt { sum -= get(i-wtime); nwin -= 1; temp[(i+offset as i32) as usize] = sum/nwin as f32; }
    data.copy_from_slice(&temp[offset..offset+nt as usize]);
}
fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    if !((0..=17).contains(&tr.trid) || tr.trid == 109) { warn!("input is not seismic data, trid=%d", tr.trid); }
    let mut mode = String::new();
    if getparstring("mode", &mut mode) == 0 { mode = "amp".into(); }
    if !["amp", "phase", "freq", "uphase", "bandwidth", "normamp", "freqw", "thin", "fdenv", "sdenv", "q"].contains(&mode.as_str()) { err!("unknown mode=\"%s\", see self-doc", mode); }
    let mut unwrap = if mode == "phase" { 0.0 } else { 1.0 };
    let mut wint = 3;
    if ["phase", "freq", "uphase", "q", "freqw", "thin"].contains(&mode.as_str()) { getparfloat("unwrap", &mut unwrap); }
    if mode == "freqw" || mode == "thin" { getparint("wint", &mut wint); }
    let mut data = vec![0.0; nt]; let mut hdata = vec![0.0; nt];
    loop {
        data.copy_from_slice(&tr.data[..nt]); hilbert(&data, &mut hdata);
        let mut amp = vec![0.0; nt]; let mut phase = vec![0.0; nt];
        for i in 0..nt {
            let re = data[i]; let im = hdata[i]; let mag = re*re + im*im;
            amp[i] = (mag as f64).sqrt() as f32;
            if mag != 0.0 { phase[i] = (im as f64).atan2(re as f64) as f32; }
        }
        let h = (2.0*PI*dt as f64) as f32;
        match mode.as_str() {
            "amp" => { tr.data[..nt].copy_from_slice(&amp); tr.trid = 123; }
            "phase" => { if unwrap != 0.0 { unwrap_phase(&mut phase, unwrap); } tr.data[..nt].copy_from_slice(&phase); tr.trid = 124; }
            "normamp" => { for i in 0..nt { tr.data[i] = (phase[i] as f64).cos() as f32; } tr.trid = 124; }
            "fdenv" | "sdenv" => { differentiate(&mut amp, h); if mode == "sdenv" { differentiate(&mut amp, h); } tr.data[..nt].copy_from_slice(&amp); tr.trid = 123; }
            "bandwidth" => { let amp2 = amp.clone(); differentiate(&mut amp, dt); for i in 0..nt { let den = 2.0*PI*amp2[i] as f64; tr.data[i] = if den != 0.0 { (amp[i] as f64 / den).abs() as f32 } else { 0.0 }; } tr.trid = 123; }
            _ => {
                let amp2 = amp.clone(); if mode == "q" { differentiate(&mut amp, dt); }
                if unwrap != 0.0 { unwrap_phase(&mut phase, unwrap); }
                differentiate(&mut phase, h); let fnyq = (0.5 / dt as f64) as f32;
                for v in &mut phase { if *v > fnyq { *v = 2.0*fnyq - *v; } }
                match mode.as_str() {
                    "freq" => { tr.data[..nt].copy_from_slice(&phase); tr.trid = 125; }
                    "uphase" => { for i in 1..nt { tr.data[0] = phase[0]; tr.data[i] += phase[i-1]; } tr.trid = 124; }
                    "freqw" => { let mut weighted = vec![0.0; nt]; for i in 0..nt { weighted[i] = amp[i]*phase[i]; } twindow(&mut amp, wint); twindow(&mut weighted, wint); for i in 0..nt { tr.data[i] = if amp[i] == 0.0 { 0.0 } else { weighted[i]/amp[i] }; } tr.trid = 125; }
                    "thin" => { let mut phase2 = vec![0.0; nt]; for i in 0..nt { phase2[i] = 2.0*fnyq - phase[i]; } twindow(&mut phase2, wint); for i in 0..nt { tr.data[i] = phase[i] - phase2[i]; } tr.trid = 125; }
                    "q" => { for i in 0..nt { tr.data[i] = if amp[i] != 0.0 { (-PI*phase[i] as f64*amp2[i] as f64/amp[i] as f64) as f32 } else { 0.0 }; } tr.trid = 125; }
                    _ => unreachable!(),
                }
            }
        }
        tr.d1 = dt; puttr(&mut tr); if gettr(&mut tr) == 0 { break; }
    }
    exit(CWP_Exit());
}
