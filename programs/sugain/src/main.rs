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

//! SUGAIN - port of su/main/amplitudes/sugain.c (Revision 1.63).
//! Credits: SEP: Jon Claerbout; CWP: Jack K. Cohen, Brian Sumner, Dave Hale.
//! Technical reference: Claerbout, Imaging the Earth, pages 233-236.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{closeto, nint};
use oxwave_tempfile::{temporary_filename, temporary_stream};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUGAIN - apply various types of gain\t\t\t\t  \t",
    "\t\t\t\t\t\t\t\t\t",
    " sugain <stdin >stdout [optional parameters]\t\t\t   \t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t  \t",
    "\tnone (no-op)\t\t\t\t\t\t    \t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t  \t",
    "\tpanel=0\t        =1  gain whole data set (vs. trace by trace)\t",
    "\ttpow=0.0\tmultiply data by t^tpow\t\t\t \t",
    "\tepow=0.0\tmultiply data by exp(epow*t)\t\t    \t",
    "\tetpow=1.0\tmultiply data by exp(epow*t^etpow)\t    \t",
    "\tgpow=1.0\ttake signed gpowth power of scaled data\t \t",
    "\tagc=0\t   flag; 1 = do automatic gain control\t     \t\t",
    "\tgagc=0\t  flag; 1 = ... with gaussian taper\t\t\t",
    "\twagc=0.5\tagc window in seconds (use if agc=1 or gagc=1)  ",
    "\ttrap=none\tzero any value whose magnitude exceeds trapval  ",
    "\tclip=none\tclip any value whose magnitude exceeds clipval  ",
    "\tpclip=none\tclip any value greater than clipval  \t\t",
    "\tnclip=none\tclip any value less than  clipval \t\t",
    "\tqclip=1.0\tclip by quantile on absolute values on trace    ",
    "\tqbal=0\t  flag; 1 = balance traces by qclip and scale     \t",
    "\tpbal=0\t  flag; 1 = bal traces by dividing by rms value   \t",
    "\tmbal=0\t  flag; 1 = bal traces by subtracting the mean    \t",
    "\tmaxbal=0\tflag; 1 = balance traces by subtracting the max ",
    "\tscale=1.0\tmultiply data by overall scale factor\t   \t",
    "\tnorm=0.0\tdivide data by overall scale factor\t     \t",
    "\tbias=0.0\tbias data by adding an overall bias value\t",
    "\tjon=0\t   \tflag; 1 means tpow=2, gpow=.5, qclip=.95\t",
    "\tverbose=0\tverbose = 1 echoes info\t\t\t\t",
    "\tmark=0\t\tapply gain only to traces with tr.mark=0\t",
    "\t\t\t=1 apply gain only to traces with tr.mark!=0    ",
    "\tvred=0\t  reducing velocity of data to use with tpow\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " \ttmpdir=\t\tif non-empty, use the value as a directory path\t",
    "\t\t\tprefix for storing temporary files; else if the ",
    "\t\t\tthe CWP_TMPDIR environment variable is set use  ",
    "\t\t\tits value for the path; else use tmpfile()\t",
    "\t\t\t\t\t\t\t\t\t",
    " Operation order:\t\t\t\t\t\t\t",
    " if (norm) scale/norm\t\t\t\t\t\t  \t",
    "\t\t\t\t\t\t\t\t\t",
    " out(t) = scale * BAL{CLIP[AGC{[t^tpow * exp(epow * t^tpow) * ( in(t)-bias )]^gpow}]}",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    "\tThe jon flag selects the parameter choices discussed in\t\t",
    "\tClaerbout's Imaging the Earth, pp 233-236.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tExtremely large/small values may be lost during agc. Windowing  ",
    "\tthese off and applying a scale in a preliminary pass through\t",
    "\tsugain may help.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tSugain only applies gain to traces with tr.mark=0. Use sushw,\t",
    "\tsuchw, suedit, or suxedit to mark traces you do not want gained.",
    "\tSee the selfdocs of sushw, suchw, suedit, and suxedit for more\t",
    "\tinformation about setting header fields. Use \"sukeyword mark\" ",
    "\tfor more information about the mark header field.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "      debias data by using mbal=1\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "      option etpow only becomes active if epow is nonzero\t\t",
];

struct Params {
    tpow: f32,
    epow: f32,
    etpow: f32,
    gpow: f32,
    vred: f32,
    trap: f32,
    clip: f32,
    pclip: f32,
    nclip: f32,
    qclip: f32,
    scale: f32,
    bias: f32,
    agc: i32,
    gagc: i32,
    pbal: i32,
    qbal: i32,
    mbal: i32,
    maxbal: i32,
    iw: usize,
    tmin: f32,
    dt: f32,
}
fn quant(a: &mut [f32], k: usize) -> f32 {
    let (mut low, mut hi) = (0i32, a.len() as i32 - 1);
    while low < hi {
        let ak = a[k];
        let (mut i, mut j) = (low, hi);
        loop {
            while a[i as usize] < ak {
                i += 1
            }
            while a[j as usize] > ak {
                j -= 1
            }
            if i <= j {
                a.swap(i as usize, j as usize);
                i += 1;
                j -= 1
            }
            if i > j {
                break;
            }
        }
        if j < (k as i32) {
            low = i
        }
        if (k as i32) < i {
            hi = j
        }
    }
    a[k]
}
fn clip(d: &mut [f32], c: f32) {
    for x in d {
        if *x > c {
            *x = c
        } else if *x < -c {
            *x = -c
        }
    }
}
fn gain(d: &mut [f32], p: &Params, offset: i32) {
    let n = d.len();
    if n == 0 {
        return;
    }
    if p.bias != 0.0 {
        for x in d.iter_mut() {
            *x += p.bias
        }
    }
    if p.tpow != 0.0 {
        for (i, x) in d.iter_mut().enumerate() {
            let mut t = p.tmin + (i as f32) * p.dt;
            if i > 0 && p.vred > 0.0 {
                t = p.tmin + ((offset as f32) / p.vred).abs() + (i as f32) * p.dt
            }
            let f = if i == 0 && p.tmin == 0.0 {
                0.0
            } else {
                (t as f64).powf(p.tpow as f64) as f32
            };
            *x *= f
        }
    }
    if p.epow != 0.0 {
        for (i, x) in d.iter_mut().enumerate() {
            let t = p.tmin + (i as f32) * p.dt;
            let f = if i == 0 && p.tmin == 0.0 {
                0.0
            } else {
                (t as f64).powf(p.etpow as f64) as f32
            };
            *x *= ((p.epow * f) as f64).exp() as f32
        }
    }
    if !closeto(p.gpow, 1.0) {
        for x in d.iter_mut() {
            let v = *x;
            *x = if closeto(p.gpow, 0.5) {
                if v >= 0.0 {
                    (v as f64).sqrt() as f32
                } else {
                    -((-v) as f64).sqrt() as f32
                }
            } else if closeto(p.gpow, 2.0) {
                v * v.abs()
            } else if v >= 0.0 {
                (v as f64).powf(p.gpow as f64) as f32
            } else {
                -((-v) as f64).powf(p.gpow as f64) as f32
            }
        }
    }
    if p.agc != 0 {
        let w = p.iw;
        let d2: Vec<f32> = d.iter().map(|v| v * v).collect();
        let mut out = vec![0.0; n];
        let mut sum = 0.0f32;
        for v in d2.iter().take(w) {
            sum += v
        }
        let mut nw = w as i32;
        let rms = sum / nw as f32;
        out[0] = if rms == 0.0 {
            0.0
        } else {
            (d[w] as f64 / (rms as f64).sqrt()) as f32
        };
        for i in 1..w {
            sum += d2[i + w - 1];
            nw += 1;
            let rms = sum / nw as f32;
            out[i] = if rms == 0.0 {
                0.0
            } else {
                (d[i] as f64 / (rms as f64).sqrt()) as f32
            }
        }
        nw += 1;
        for i in w..n - w {
            sum = 0.0;
            for v in &d2[i - w..i + w] {
                sum += v
            }
            let rms = sum / nw as f32;
            out[i] = if rms == 0.0 {
                0.0
            } else {
                (d[i] as f64 / (rms as f64).sqrt()) as f32
            }
        }
        for i in n - w..n {
            sum = 0.0;
            for v in &d2[i - w..n] {
                sum += v
            }
            nw -= 1;
            let rms = sum / nw as f32;
            out[i] = if d[i] == 0.0 {
                0.0
            } else {
                (d[i] as f64 / (rms as f64).sqrt()) as f32
            }
        }
        d.copy_from_slice(&out);
    }
    if p.gagc != 0 {
        let w = p.iw;
        let u = (3.8090232f64 / w as f64) as f32;
        let usq = u * u;
        let d2: Vec<f32> = d.iter().map(|v| v * v).collect();
        let mut s = d2.clone();
        for j in 1..w {
            let f = j as f32;
            let weight = (-(usq * f * f) as f64).exp() as f32;
            for i in j..n {
                s[i] += weight * d2[i - j]
            }
            for i in 0..n - j {
                s[i] += weight * d2[i + j]
            }
        }
        for i in 0..n {
            d[i] = if s[i] == 0.0 {
                0.0
            } else {
                (d[i] as f64 / (s[i] as f64).sqrt()) as f32
            }
        }
    }
    if p.trap > 0.0 {
        for x in d.iter_mut() {
            if x.abs() > p.trap {
                *x = 0.0
            }
        }
    }
    if p.clip > 0.0 {
        clip(d, p.clip)
    }
    if p.pclip < f32::MAX {
        for x in d.iter_mut() {
            if *x > p.pclip {
                *x = p.pclip
            }
        }
    }
    if p.nclip > -f32::MAX {
        for x in d.iter_mut() {
            if *x < p.nclip {
                *x = p.nclip
            }
        }
    }
    if p.qclip < 1.0 || p.qbal != 0 {
        let bal = if p.qclip == 1.0 {
            let mut b = d[0].abs();
            for x in &d[1..] {
                let v = x.abs();
                b = if b > v { b } else { v }
            }
            b
        } else {
            let mut a: Vec<f32> = d.iter().map(|x| x.abs()).collect();
            let k = ((p.qclip * n as f32) as f64 - 0.5) as usize;
            quant(&mut a, k)
        };
        if p.qbal != 0 {
            if bal != 0.0 {
                for x in d.iter_mut() {
                    *x /= bal
                }
                if p.qclip < 1.0 {
                    clip(d, 1.0)
                }
            }
        } else {
            clip(d, bal)
        }
    }
    if p.pbal != 0 {
        let mut s = 0.0f32;
        for x in d.iter() {
            s += x * x
        }
        let rms = ((s / n as f32) as f64).sqrt() as f32;
        if rms != 0.0 {
            for x in d.iter_mut() {
                *x /= rms
            }
        }
    }
    if p.mbal != 0 {
        let mut mean = 0.0f32;
        for x in d.iter() {
            mean += x
        }
        mean /= n as f32;
        if mean != 0.0 {
            for x in d.iter_mut() {
                *x -= mean
            }
        }
    }
    if p.maxbal != 0 {
        let mut m = d[0];
        for x in d.iter() {
            if *x > m {
                m = *x
            }
        }
        for x in d.iter_mut() {
            *x -= m
        }
    }
    if !closeto(p.scale, 1.0) {
        for x in d.iter_mut() {
            *x *= p.scale
        }
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default()
    }
    if !tmpdir.is_empty() {
        let writable = match std::fs::metadata(&tmpdir) {
            Ok(m) if m.is_dir() => temporary_stream(Some(&tmpdir)).is_some(),
            Ok(_) => std::fs::OpenOptions::new().write(true).open(&tmpdir).is_ok(),
            Err(_) => false,
        };
        if !writable {
            err!("you can't write in %s (or it doesn't exist)", tmpdir)
        }
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace")
    };
    let nt = tr.ns as usize;
    let mut dt = (tr.dt as f64 / 1000000.0) as f32;
    let tmin = (tr.delrt as f64 / 1000.0) as f32;
    if dt == 0.0 {
        getparfloat("dt", &mut dt);
        if dt == 0.0 {
            mustgetparfloat("dt", &mut dt);
        }
    }
    macro_rules! f {
        ($name:literal,$default:expr) => {{
            let mut v: f32 = $default;
            getparfloat($name, &mut v);
            v
        }};
    }
    macro_rules! i {
        ($name:literal) => {{
            let mut v = 0i32;
            getparint($name, &mut v);
            v
        }};
    }
    let mut p = Params {
        tpow: f!("tpow", 0.0),
        epow: f!("epow", 0.0),
        etpow: f!("etpow", 1.0),
        gpow: f!("gpow", 1.0),
        vred: f!("vred", 0.0),
        trap: f!("trap", 0.0),
        clip: f!("clip", 0.0),
        pclip: f!("pclip", f32::MAX),
        nclip: f!("nclip", -f32::MAX),
        qclip: f!("qclip", 1.0),
        scale: f!("scale", 1.0),
        bias: 0.0,
        agc: 0,
        gagc: 0,
        pbal: 0,
        qbal: 0,
        mbal: 0,
        maxbal: 0,
        iw: 0,
        tmin,
        dt,
    };
    let norm = f!("norm", 0.0);
    p.bias = f!("bias", 0.0);
    let wagc = f!("wagc", 0.5);
    p.agc = i!("agc");
    p.gagc = i!("gagc");
    p.pbal = i!("pbal");
    p.qbal = i!("qbal");
    p.mbal = i!("mbal");
    let panel = i!("panel");
    let jon = i!("jon");
    let verbose = i!("verbose");
    p.maxbal = i!("maxbal");
    let mark = i!("mark");
    checkpars();
    if p.vred < 0.0 {
        err!("vred = %f, must be positive", p.vred)
    }
    if p.trap < 0.0 {
        err!("trap = %f, must be positive", p.trap)
    }
    if p.clip < 0.0 {
        err!("clip = %f, must be positive", p.clip)
    }
    if p.qclip < 0.0 || p.qclip > 1.0 {
        err!("qclip = %f, must be between 0 and 1", p.qclip)
    }
    if p.agc != 0 || p.gagc != 0 {
        let iw = nint((wagc / dt) as f64);
        if iw < 1 {
            err!("wagc=%g must be positive", wagc)
        }
        if iw > nt as i32 {
            err!("wagc=%g too long for trace", wagc)
        }
        p.iw = (iw >> 1) as usize;
        if verbose != 0 {
            warn!("half window: %d samples", p.iw as i32)
        }
    }
    if jon != 0 {
        p.tpow = 2.0;
        p.gpow = 0.5;
        p.qclip = 0.95
    }
    if norm != 0.0 {
        p.scale /= norm
    }
    if panel == 0 {
        loop {
            if (tr.mark == 0 && mark == 0) || (tr.mark != 0 && mark != 0) {
                let offset = tr.offset;
                gain(&mut tr.data[..nt], &p, offset)
            }
            puttr(&mut tr);
            if gettr(&mut tr) == 0 {
                break;
            }
        }
    } else {
        let (tracefp, headerfp, tracefile, headerfile) = if tmpdir.is_empty() {
            let tracefp = etmpfile();
            let headerfp = etmpfile();
            if verbose != 0 {
                warn!("using tmpfile() call");
            }
            (tracefp, headerfp, String::new(), String::new())
        } else {
            let tracefile = temporary_filename(&tmpdir);
            let headerfile = temporary_filename(&tmpdir);
            let tracefp = efopen(&tracefile, "w+");
            let headerfp = efopen(&headerfile, "w+");
            if verbose != 0 {
                warn!("putting temporary files in %s", tmpdir);
            }
            (tracefp, headerfp, tracefile, headerfile)
        };
        let mut ntr = 0;
        loop {
            ntr += 1;
            efwrite(tr.hdr.as_bytes(), &headerfp);
            efwrite(&tr.data[..nt], &tracefp);
            if gettr(&mut tr) == 0 {
                break;
            }
        }
        erewind(&tracefp);
        erewind(&headerfp);
        let mut data = vec![0.0; nt * ntr];
        efread(&mut data, &tracefp);
        efclose(&tracefp);
        if !tmpdir.is_empty() {
            eremove(&tracefile);
        }
        gain(&mut data, &p, tr.offset);
        for j in 0..ntr {
            tr.data[..nt].copy_from_slice(&data[j * nt..(j + 1) * nt]);
            efread(tr.hdr.as_bytes_mut(), &headerfp);
            puttr(&mut tr)
        }
        efclose(&headerfp);
        if !tmpdir.is_empty() {
            eremove(&headerfile);
        }
    }
    exit(CWP_Exit());
}
