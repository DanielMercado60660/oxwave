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

//! SUGROLL - Ground roll supression using Karhunen-Loeve transform.
//! Port of su/main/noise/sugroll.c (Revision 1.1).
//! Original copyright (c) Colorado School of Mines, 1999. All rights reserved.
//! Credits: BRGM: Adnand Bitri, 1999.
//! Reference: Xuewei Liu, F., 1999, Ground roll supresion using the
//! Karhunen-Loeve transform, Geophysics vol. 64 No. 2 pp 564-566.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUGROLL - Ground roll supression using Karhunen-Loeve transform\t",
    "                                                                       ",
    "    sugroll <infile >outfile  [optional parameters]                 \t",
    "                                                                       ",
    " Optional Parameters:                                                  ",
    " dt=tr.dt (from header) \ttime sampling interval (secs)           ",
    " nx=ntr   (counted from data)\tnumber of horizontal samples (traces)\t",
    " sb=0      1 for the graund-roll                                       ",
    " verbose=0\tverbose = 1 echoes information\t\t\t",
    " nrot=3        the principal components for the m largest eigenvalues ",
    " tmpdir= \t if non-empty, use the value as a directory path",
    "\t\t prefix for storing temporary files; else if the",
    "\t         the CWP_TMPDIR environment variable is set use\t",
    "\t         its value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t",
    " Notes:                                                                ",
    " The method developed here is to extract the ground-roll from the\t",
    " common-shot gathers using Karhunen-Loeve transform, and then to substract  ",
    " it from the original data. The advantage is the ground-roll is suppresed  ",
    " with negligible distortion of the signal.\t",
    "                                                                       ",
];

// The C program keeps these files global for its interrupt handler. The
// frozen Rust core has no signal-handler API; keep filenames local here.
struct TempFiles {
    tracefile: String,
    headerfile: String,
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        // Rust's permission query covers missing paths and read-only mode
        // bits. A faithful access(WRITE_OK) including ACLs needs a core API.
        if !std::fs::metadata(&tmpdir).is_ok_and(|m| !m.permissions().readonly()) {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let dtheader = (tr.dt as f32 as f64 / 1000000.0) as f32;
    let mut sb = 0;
    if getparint("sb", &mut sb) == 0 {
        sb = 0;
    }
    let mut nrot = 3;
    if getparint("nrot", &mut nrot) == 0 {
        nrot = 3;
    }

    let istmpdir = !tmpdir.is_empty();
    let mut files = TempFiles {
        tracefile: String::new(),
        headerfile: String::new(),
    };
    let (tracefp, headerfp) = if !istmpdir {
        let tracefp = etmpfile();
        let headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
        (tracefp, headerfp)
    } else {
        files.tracefile = temporary_filename(&tmpdir);
        files.headerfile = temporary_filename(&tmpdir);
        // The library returns a fixed fallback when its /tmp mkstemp is
        // denied by the sandbox. Keep the two storage files distinct then.
        if files.headerfile == files.tracefile {
            files.headerfile.push_str(".header");
        }
        let tracefp = efopen(&files.tracefile, "w+");
        let headerfp = efopen(&files.headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
        (tracefp, headerfp)
    };
    let mut ntr = 0usize;
    loop {
        ntr += 1;
        efwrite(tr.hdr.as_bytes(), &headerfp);
        efwrite(&tr.data[..nt], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    let mut nx = ntr as i32;
    if getparint("nx", &mut nx) == 0 {
        nx = ntr as i32;
    }
    let mut dt = dtheader;
    if getparfloat("dt", &mut dt) == 0 {
        dt = dtheader;
    }
    // There is no checkpars() in the original program.
    if dt == 0.0 {
        err!("header field dt not set, must be getparred");
    }

    let mut in_traces = vec![vec![0.0f32; nt]; ntr];
    let mut out_traces = vec![vec![0.0f32; nt]; ntr];
    let mut xt = vec![vec![0.0f32; ntr]; nt];
    let mut g = vec![vec![0.0f32; ntr]; ntr];
    let mut v = vec![vec![0.0f32; ntr]; ntr];
    let mut u = vec![vec![0.0f32; nt]; ntr];
    let mut w = vec![0.0f32; ntr];
    erewind(&headerfp);
    erewind(&tracefp);
    for ix in 0..ntr {
        fread(&mut in_traces[ix], &tracefp);
    }
    efclose(&tracefp);

    // Transpose the data; construct G = X * XT in the C accumulation order.
    for ix in 0..ntr {
        for it in 0..nt {
            xt[it][ix] = in_traces[ix][it];
        }
    }
    for i in 0..ntr {
        for j in 0..ntr {
            g[i][j] = 0.0;
            for k in 0..nt {
                g[i][j] += in_traces[i][k] * xt[k][j];
            }
        }
    }
    for i in 0..ntr {
        for j in 0..ntr {
            if g[i][j] != g[j][i] {
                warn!(" G n'est pas symetrique ");
            }
        }
    }
    let mut nrots = 0;
    jacobi(&mut g, ntr, &mut w, &mut v, &mut nrots);
    eigsrt(&mut w, &mut v, ntr);
    for i in 0..ntr {
        let singval = w[i] / w[0];
        warn!(" Singular value number singval=%g", singval);
    }

    // Psi = transpose(v) * X; retain the selected rows; X' = v * Psi'.
    for i in 0..ntr {
        for j in 0..nt {
            u[i][j] = 0.0;
            for k in 0..ntr {
                u[i][j] += v[k][i] * in_traces[k][j];
            }
        }
    }
    for i in nrot..nx {
        // Invalid component bounds dereference an invalid row pointer in C.
        // Match its shell-visible status (128 + SIGSEGV) without a Rust panic.
        // Delivering SIGSEGV itself requires a signal API in the frozen core.
        if nt != 0 && (i < 0 || i as usize >= ntr) {
            exit(139);
        }
        for j in 0..nt {
            u[i as usize][j] = 0.0;
        }
    }
    for i in 0..ntr {
        for j in 0..nt {
            out_traces[i][j] = 0.0;
            for k in 0..ntr {
                out_traces[i][j] += v[i][k] * u[k][j];
            }
        }
    }
    if istmpdir {
        eremove(&files.tracefile);
    }
    erewind(&headerfp);
    for itr in 0..ntr {
        efread(tr.hdr.as_bytes_mut(), &headerfp);
        for it in 0..nt {
            tr.data[it] = if sb == 0 {
                in_traces[itr][it] - out_traces[itr][it]
            } else {
                out_traces[itr][it]
            };
        }
        puttr(&mut tr);
    }
    efclose(&headerfp);
    if istmpdir {
        eremove(&files.headerfile);
    }
    exit(CWP_Exit());
}

fn eigsrt(d: &mut [f32], v: &mut [Vec<f32>], n: usize) {
    for i in 0..n - 1 {
        let mut k = i;
        let mut p = d[i];
        for j in i + 1..n {
            if d[j] >= p {
                k = j;
                p = d[j];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 0..n {
                p = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }
}

fn jacobi(a: &mut [Vec<f32>], n: usize, d: &mut [f32], v: &mut [Vec<f32>], nrots: &mut i32) {
    let mut b = vec![0.0f32; n];
    let mut z = vec![0.0f32; n];
    for ip in 0..n {
        for iq in 0..n {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    for ip in 0..n {
        d[ip] = a[ip][ip];
        b[ip] = d[ip];
        z[ip] = 0.0;
    }
    *nrots = 0;
    for i in 1..=50 {
        let mut sm = 0.0f32;
        for ip in 0..n - 1 {
            for iq in ip + 1..n {
                sm = (sm as f64 + (a[ip][iq] as f64).abs()) as f32;
            }
        }
        if sm == 0.0 {
            return;
        }
        let tresh = if i < 4 {
            (0.2 * sm as f64 / (n * n) as f64) as f32
        } else {
            0.0
        };
        for ip in 0..n - 1 {
            for iq in ip + 1..n {
                let g = (100.0 * (a[ip][iq] as f64).abs()) as f32;
                if i > 4
                    && ((d[ip] as f64).abs() + g as f64) as f32 == (d[ip] as f64).abs() as f32
                    && ((d[iq] as f64).abs() + g as f64) as f32 == (d[iq] as f64).abs() as f32
                {
                    a[ip][iq] = 0.0;
                } else if (a[ip][iq] as f64).abs() > tresh as f64 {
                    let mut h = d[iq] - d[ip];
                    let mut t;
                    if ((h as f64).abs() + g as f64) as f32 == (h as f64).abs() as f32 {
                        t = a[ip][iq] / h;
                    } else {
                        let theta = (0.5 * h as f64 / a[ip][iq] as f64) as f32;
                        t = (1.0 / ((theta as f64).abs() + (1.0 + (theta * theta) as f64).sqrt()))
                            as f32;
                        if theta < 0.0 {
                            t = -t;
                        }
                    }
                    // In C, 1+t*t is evaluated entirely in float before sqrt.
                    let c = (1.0 / ((1.0f32 + t * t) as f64).sqrt()) as f32;
                    let s = t * c;
                    let tau = (s as f64 / (1.0 + c as f64)) as f32;
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    // Preserve the original ip-1 / iq-1 bounds, including empty ranges.
                    for j in 0..(ip as i32 - 1) {
                        rotate(a, j as usize, ip, j as usize, iq, s, tau);
                    }
                    for j in ip + 1..iq - 1 {
                        rotate(a, ip, j, j, iq, s, tau);
                    }
                    for j in iq + 1..n {
                        rotate(a, ip, j, iq, j, s, tau);
                    }
                    for j in 0..n {
                        rotate(v, j, ip, j, iq, s, tau);
                    }
                    *nrots += 1;
                }
            }
        }
        for ip in 0..n {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    fprintf!(stderr(), "To many iterations in routine jacobi\n");
    exit(-2);
}

// The ROTATE macro from sugroll.c; all operations are float.
fn rotate(a: &mut [Vec<f32>], i: usize, j: usize, k: usize, l: usize, s: f32, tau: f32) {
    let g = a[i][j];
    let h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}
