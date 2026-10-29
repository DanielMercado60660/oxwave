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

//! SUMAX - trace and global extrema, port of sumax.c revision 1.23.
//! Credits: John Stockwell, Garry Perratt, Reginald Beardsley, Shuki Ronen, Jack K. Cohen, Gerald Klein.
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::abs;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMAX - get trace by trace local/global maxima, minima, or absolute maximum",
    "\t\t\t\t\t\t\t\t\t",
    " sumax <stdin >stdout [optional parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    "\tnone\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters: \t\t\t\t\t\t",
    "\toutput=ascii \t\twrite ascii data to outpar\t\t",
    "\t\t\t\t=binary for binary floats to stdout\t",
    "\t\t\t\t=segy for SEGY traces to stdout\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tmode=maxmin\t\toutput both minima and maxima\t\t",
    "\t\t\t\t=max maxima only\t\t\t",
    "\t\t\t\t=min minima only\t\t\t",
    "\t\t\t\t=abs absolute maxima only      \t\t",
    "\t\t\t\t=rms RMS \t\t      \t\t",
    "\t\t\t\t=thd search first max above threshold\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tthreshamp=0\t\tthreshold amplitude value\t\t",
    "\tthreshtime=0\t\ttmin to start search for threshold \t",
    "\t\t\t\t\t\t\t\t\t",
    "\tverbose=0 \t\twrites global quantities to outpar\t",
    "\t\t\t\t=1 trace number, values, sample location",
    "\t\t\t\t=2 key1 & key2 instead of trace number  ",
    "\tkey1=fldr\t\tkey for verbose=2                       ",
    "\tkey2=ep\t\t\tkey for verbose=2                       ",
    "\t\t\t\t\t\t\t\t\t",
    "\toutpar=/dev/tty\t\toutput parameter file; contains output\t",
    "\t\t\t\t\tfrom verbose\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples: \t\t\t\t\t\t\t\t",
    " For global max and min values:  sumax < segy_data\t\t\t",
    " For local and global max and min values:  sumax < segy_data verbose=1\t",
    " To plot values specified by mode:\t\t\t\t\t",
    "    sumax < segy_data output=binary mode=modeval | xgraph n=npairs\t",
    " To plot seismic data with the only values nonzero being those specified",
    " by mode=modeval:\t\t\t\t\t\t\t",
    "    sumax < segy_data output=segy mode=modeval | suxwigb\t\t",
    "  \t\t\t\t\t\t\t\t\t",
    " Note:\twhile traces are counted from 1, sample values are counted from 0.",
    "\tAlso, if multiple min, max, or abs max values exist on a trace,\t",
    "       only the first one is captured.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " See also: suxmax, supsmax\t\t\t\t\t\t",
];
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    let mut output = String::new();
    if getparstring("output", &mut output) == 0 {
        output = "ascii".into();
    }
    if !["ascii", "binary", "segy"].contains(&output.as_str()) {
        err!("%s unknown output selection", output);
    }
    let mut key1 = String::new();
    let mut key2 = String::new();
    if getparstring("key1", &mut key1) == 0 {
        key1 = "fldr".into();
    }
    if getparstring("key2", &mut key2) == 0 {
        key2 = "ep".into();
    }
    let mut mode = String::new();
    if getparstring("mode", &mut mode) == 0 {
        mode = "maxmin".into();
    }
    if !["maxmin", "max", "min", "thd", "abs", "rms"].contains(&mode.as_str()) {
        err!("%s unknown mode", mode);
    }
    if mode == "rms" && output == "segy" {
        err!("%s output incompatible with %s mode", output, mode);
    }
    if output == "binary" && mode == "maxmin" {
        mode = "max".into();
        warn!("outputting max values");
    }
    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 {
        outpar = "/dev/tty".into();
    }
    let fp = efopen(&outpar, "w");
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut dt = 0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    let mut gmax = tr.data[0];
    let mut gmin = gmax;
    let mut gabsmax = abs(gmax);
    let mut grmssumsq = gmax * gmax;
    let mut thamp = 0f32;
    let mut thtime = 0f32;
    if mode == "thd" {
        getparfloat("threshamp", &mut thamp);
        getparfloat("threshtime", &mut thtime);
    }
    checkpars();
    let mut itr = 0i32;
    let (mut maxit, mut minit, mut absit, mut thit) = (0usize, 0usize, 0usize, 0usize);
    let (mut gmaxit, mut gminit, mut gabsit, mut thgit) = (0usize, 0usize, 0usize, 0usize);
    let (mut gmaxtr, mut gmintr, mut gabstr, mut thgtr) = (1i32, 1i32, 1i32, 0i32);
    let mut thgmax = 0f32;
    loop {
        itr += 1;
        let (mut id1, mut id2) = (0, 0);
        if verbose == 2 {
            let mut v = Value::default();
            gethdval(&tr, &key1, &mut v);
            id1 = vtoi(hdtype(&key1), v);
            gethdval(&tr, &key2, &mut v);
            id2 = vtoi(hdtype(&key2), v);
        }
        let mut max = tr.data[0];
        let mut min = max;
        let mut absmax = abs(max);
        let mut sum = max * max;
        grmssumsq += sum;
        let mut thmax = 0f32;
        let mut thfit = (thtime / dt) as i32;
        for it in 1..nt {
            if tr.data[it] > max {
                max = tr.data[it];
                maxit = it;
            }
            if tr.data[it] < min {
                min = tr.data[it];
                minit = it;
            }
            if abs(tr.data[it]) > absmax {
                absmax = abs(tr.data[it]);
                absit = it;
            }
            if it as i32 >= thfit {
                if tr.data[it] > thamp {
                    if tr.data[it] > thmax {
                        thmax = tr.data[it];
                        thit = it;
                    }
                } else if thmax > thamp {
                    thfit = nt as i32 + 1;
                }
            }
            if mode == "rms" {
                let val = tr.data[it];
                sum += val * val;
                grmssumsq += val * val;
            }
            tr.data[it] = 0.0;
        }
        let mut rms = 0f32;
        if mode == "rms" {
            rms = ((sum / nt as f32) as f64).sqrt() as f32;
        }
        if output == "segy" {
            if mode == "maxmin" {
                tr.data[maxit] = max;
                tr.data[minit] = min;
            }
            if mode == "max" {
                tr.data[maxit] = max;
            }
            if mode == "min" {
                tr.data[minit] = min;
            }
            if mode == "abs" {
                tr.data[absit] = absmax;
            }
            if mode == "thd" {
                tr.data[thit] = thmax;
            }
        }
        if max > gmax {
            gmax = max;
            gmaxtr = itr;
            gmaxit = maxit;
        }
        if min < gmin {
            gmin = min;
            gmintr = itr;
            gminit = minit;
        }
        if absmax > gabsmax {
            gabsmax = absmax;
            gabstr = itr;
            gabsit = absit;
        }
        if thmax > thgmax {
            thgmax = thmax;
            thgtr = itr;
            thgit = thit;
        }
        if output == "ascii" && verbose != 0 {
            let ident = if verbose == 2 {
                sprintf!("%d %d", id1, id2)
            } else {
                sprintf!("%d", itr)
            };
            if mode == "maxmin" {
                fprintf!(
                    fp,
                    "%s %e at point %d time %e\n",
                    ident,
                    max,
                    maxit as i32,
                    dt * maxit as f32
                );
                fprintf!(
                    fp,
                    "   %e at point %d time %e\n",
                    min,
                    minit as i32,
                    dt * minit as f32
                );
            }
            let selected = match mode.as_str() {
                "max" => Some((max, maxit)),
                "min" => Some((min, minit)),
                "abs" => Some((absmax, absit)),
                "thd" => Some((thmax, thit)),
                _ => None,
            };
            if let Some((v, i)) = selected {
                fprintf!(
                    fp,
                    "%s %e at point %d time %e\n",
                    ident,
                    v,
                    i as i32,
                    dt * i as f32
                );
            }
            if mode == "rms" {
                fprintf!(fp, "%s %e\n", ident, rms);
            }
        }
        if output == "binary" {
            let v = match mode.as_str() {
                "max" => max,
                "min" => min,
                "abs" => absmax,
                "thd" => thmax,
                "rms" => rms,
                _ => 0.0,
            };
            efwrite(&[itr as f32], &stdout());
            efwrite(&[v], &stdout());
        }
        if output == "segy" {
            puttr(&mut tr);
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    if mode == "maxmin" {
        if verbose != 0 {
            fprintf!(
                fp,
                "global max = %e, at trace = %d, at sample = %d\n",
                gmax,
                gmaxtr,
                gmaxit as i32
            );
            fprintf!(
                fp,
                "global min = %e, at trace = %d, at sample = %d\n",
                gmin,
                gmintr,
                gminit as i32
            );
        } else {
            fprintf!(fp, "%e %e\n", gmax, gmin);
        }
    }
    let global = match mode.as_str() {
        "max" => Some((
            "global max = %e, at trace = %d, at sample = %d\n",
            gmax,
            gmaxtr,
            gmaxit,
        )),
        "min" => Some((
            "global min = %e, at trace = %d, at sample = %d\n",
            gmin,
            gmintr,
            gminit,
        )),
        "abs" => Some((
            "global abs max = %e, at trace = %d, at sample = %d\n",
            gabsmax,
            gabstr,
            gabsit,
        )),
        "thd" => Some((
            "global max = %e, at trace = %d, at sample = %d\n",
            thgmax,
            thgtr,
            thgit,
        )),
        _ => None,
    };
    if let Some((fmt, v, t, i)) = global {
        if verbose != 0 {
            fprintf!(fp, fmt, v, t, i as i32);
        } else {
            fprintf!(fp, "%e\n", v);
        }
    }
    if mode == "rms" {
        let rms = ((grmssumsq / itr.wrapping_mul(nt as i32) as f32) as f64).sqrt() as f32;
        if verbose != 0 {
            fprintf!(fp, "global rms = %e\n", rms);
        } else {
            fprintf!(fp, "%e\n", rms);
        }
    }
    if output == "binary" {
        fprintf!(fp, "npairs = %d\n", itr);
    }
    exit(CWP_Exit());
}
