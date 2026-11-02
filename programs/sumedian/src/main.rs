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

//! SUMEDIAN - median filter along a polygonal curve. Port of su/main/filters/sumedian.c (Revision 1.16).
//! Credits: John Stockwell, Werner M. Heigl, Richard Hobbs.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_linear_interp::intlin;
use oxwave_tempfile::temporary_filename;
use std::ffi::{c_char, c_int, CString};
use std::io::{BufReader, Read};
use std::sync::atomic::{AtomicPtr, Ordering};

extern "C" {
    fn access(path: *const c_char, mode: c_int) -> c_int;
    fn signal(sig: c_int, handler: extern "C" fn(c_int)) -> usize;
}

// Globals in C, so closefiles can remove the named files on interruption.
#[derive(Default)]
struct TempFiles {
    tracefile: String,
    headerfile: String,
    tracefp: Option<CFile>,
    headerfp: Option<CFile>,
}

static FILES: AtomicPtr<TempFiles> = AtomicPtr::new(std::ptr::null_mut());

extern "C" fn closefiles(_sig: c_int) {
    // main retains this allocation for the lifetime of the handlers.
    let files = unsafe { &*FILES.load(Ordering::Relaxed) };
    if let Some(fp) = &files.headerfp {
        efclose(fp);
    }
    if let Some(fp) = &files.tracefp {
        efclose(fp);
    }
    eremove(&files.headerfile);
    eremove(&files.tracefile);
    exit(EXIT_FAILURE);
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t   \t\t\t\t\t\t\t\t",
    " SUMEDIAN - MEDIAN filter about a user-defined polygonal curve with\t",
    "\t   the distance along the curve specified by key header word \t",
    " \t   \t\t\t\t\t\t\t\t",
    " sumedian <stdin >stdout xshift= tshift= [optional parameters]\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " xshift=               array of position values as specified by\t",
    "                       the `key' parameter\t\t\t\t",
    " tshift=               array of corresponding time values (sec)\t",
    "  ... or input via files:\t\t\t\t\t\t",
    " nshift=               number of x,t values defining median times\t",
    " xfile=                file containing position values as specified by\t",
    "                       the `key' parameter\t\t\t\t",
    " tfile=                file containing corresponding time values (sec)\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " key=tracl             Key header word specifying trace number \t",
    "                       =offset  use trace offset instead\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " mix=.6,1,1,1,.6       array of weights for mix (weighted moving average)",
    " median=0              =0  for mix\t\t\t\t\t",
    "                       =1  for median filter\t\t\t\t",
    " nmed=5                odd no. of traces to median filter\t\t",
    " sign=-1               =-1  for upward shift\t\t\t\t",
    "                       =+1  for downward shift\t\t\t\t",
    " subtract=1            =1  subtract filtered data from input\t\t",
    "                       =0  don't subtract\t\t\t\t",
    " verbose=0             =1  echoes information\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " ------\t\t\t\t\t\t\t\t",
    " Median filtering is a process for suppressing a particular moveout on ",
    " seismic sections. Its advantage over traditional dip filtering is that",
    " events with arbitrary moveout may be suppressed. Median filtering is\t",
    " commonly used in up/down wavefield separation of VSP data.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The process generally consists of 3 steps:\t\t\t\t",
    " 1. A copy of the data panel is shifted such that the polygon in x,t\t",
    "    specifying moveout is flattened to horizontal. The x,t pairs are \t",
    "    specified either by the vector xshift,tshift or by the values in\t",
    "    the datafiles xfile,tfile.\tFor fractional shift, the shifted data\t",
    "    is interpolated.\t\t\t\t\t\t\t",
    " 2. Then a mix (weighted moving average) is performed over the shifted\t",
    "    panel to emphasize events with the specified moveout and suppress\t",
    "    events with other moveouts.\t\t\t\t\t",
    " 3. The panel is then shifted back (and interpolated) to its original\t",
    "    moveout, and subtracted from the original data. Thus all events\t",
    "    with the user-specified moveout are removed.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For VSP data the following modifications are provided:\t\t",
    " 1. The moveout polygon in x,t is usually the first break times for\t",
    "    each trace. The parameter sign allows for downward shift in order\t",
    "    align upgoing events.\t\t\t\t\t\t",
    " 2. Alternative to a mix, a median filter can be applied by setting\t",
    "    the parameter median=1 and nmed= to the number of traces filtered.\t",
    " 3. By setting subtract=0 the filtered panel is only shifted back but\t",
    "    not subtracted from the original data.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The values of tshift are linearly interpolated for traces falling\t",
    " between given xshift values. The tshift interpolant is extrapolated\t",
    " to the left by the smallest time sample on the trace and to the right\t",
    " by the last value given in the tshift array. \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The files tfile and xfile are files of binary (C-style) floats.\t",
    " \t\t\t\t\t\t\t\t\t",
    " The number of values defined by mix=val1,val2,... determines the\t",
    " number of traces to be averaged, the values determine the weights.\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " The median filter may perform poorly on the edges of a section.\t",
    " Choosing larger beginning and ending mix values may help, but may\t",
    " also introduce additional artifacts.\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples:\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn read_file(name: &str, n: usize, kind: &str) -> Vec<f32> {
    let file = match std::fs::File::open(name) {
        Ok(file) => file,
        Err(_) => err!("cannot open %s=%s\n", kind, name),
    };
    let mut file = BufReader::new(file);
    let mut values = Vec::new();
    // Read only the requested floats, including for pipes and device files.
    for _ in 0..n {
        let mut bytes = [0; 4];
        if file.read_exact(&mut bytes).is_err() {
            err!("error reading %s=%s\n", kind, name);
        }
        if values.len() == values.capacity() && values.try_reserve(1).is_err() {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        values.push(f32::from_ne_bytes(bytes));
    }
    values
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut xf = String::new();
    let mut tf = String::new();
    let (xs, ts) = if getparstring("xfile", &mut xf) != 0 && getparstring("tfile", &mut tf) != 0 {
        let mut n = 0;
        mustgetparint("nshift", &mut n);
        // A negative C allocation size fails before either file is opened.
        if n < 0 {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        (
            read_file(&xf, n as usize, "xfile"),
            read_file(&tf, n as usize, "tfile"),
        )
    } else {
        let n = countparval("xshift");
        if n == 0 {
            err!("must give xshift= vector");
        }
        let m = countparval("tshift");
        if m == 0 {
            err!("must give tshift= vector");
        }
        if n != m {
            err!("lengths of xshift, tshift must be the same");
        }
        let mut x = vec![0.0; n as usize];
        let mut t = vec![0.0; n as usize];
        getparfloat("xshift", &mut x);
        getparfloat("tshift", &mut t);
        (x, t)
    };
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 {
        key = "tracl".into();
    }
    let typ = hdtype(&key);
    let index = getindex(&key);
    let mut nmix = countparval("mix") as usize;
    let mut mix = if nmix != 0 {
        let mut v = vec![0.0; nmix];
        getparfloat("mix", &mut v);
        if nmix % 2 == 0 {
            err!("number of mixing coefficients must be odd");
        }
        v
    } else {
        nmix = 5;
        vec![0.6, 1.0, 1.0, 1.0, 0.6]
    };
    let mut median = 0;
    let mut nmed = 0;
    let mut sign = -1;
    let mut subtract = 1;
    let mut verbose = 0;
    getparint("median", &mut median);
    if getparint("nmed", &mut nmed) == 0 && median != 0 {
        nmed = 5;
    }
    getparint("sign", &mut sign);
    getparint("subtract", &mut subtract);
    getparint("verbose", &mut verbose);
    if median != 0 && nmed % 2 == 0 {
        nmed += 1;
        warn!("increased nmed by 1 to ensure it is odd");
    }
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_str()).unwrap();
        // WRITE_OK is W_OK (2): test access before reading the first trace,
        // but defer opening the temporary files until after that read.
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }
    if median == 0 {
        for w in &mut mix {
            *w /= nmix as f32;
        }
    } else {
        nmix = nmed as usize;
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't read first trace");
    }
    if tr.dt == 0 {
        err!("dt header field must be set");
    }
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let nt = tr.ns as usize;
    let databytes = 4 * nt;
    let mut tracefile = String::new();
    let mut headerfile = String::new();
    let (tracefp, headerfp) = if tmpdir.is_empty() {
        (etmpfile(), etmpfile())
    } else {
        tracefile = temporary_filename(&tmpdir);
        headerfile = temporary_filename(&tmpdir);
        let files = Box::leak(Box::new(TempFiles {
            tracefile: tracefile.clone(),
            headerfile: headerfile.clone(),
            ..TempFiles::default()
        }));
        FILES.store(files, Ordering::Relaxed);
        // SIGINT, SIGQUIT, SIGHUP and SIGTERM, as in C.
        for sig in [2, 3, 1, 15] {
            unsafe {
                signal(sig, closefiles);
            }
        }
        let tracefp = efopen(&tracefile, "w+");
        files.tracefp = Some(tracefp.clone());
        let headerfp = efopen(&headerfile, "w+");
        files.headerfp = Some(headerfp.clone());
        (tracefp, headerfp)
    };
    if verbose != 0 {
        if tmpdir.is_empty() {
            warn!("using tmpfile() call");
        } else {
            warn!("putting temporary files in %s", tmpdir);
        }
    }
    // As in C, spool just the header and the nt actual samples. Memory use
    // depends on the mixing window and one integer shift per input trace.
    let mut ntr = 0;
    loop {
        ntr += 1;
        efwrite(tr.hdr.as_bytes(), &headerfp);
        efwrite(&tr.data_bytes()[..databytes], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    erewind(&headerfp);
    erewind(&tracefp);
    let mut shifts = vec![0i32; ntr];
    let mut shiftmax = 0;
    let mut shiftmin = 0;
    // With nshift=0 C dereferences tshift[-1] during interpolation.
    // Keep file-open, temporary-file and input errors in their original order.
    if xs.is_empty() {
        err!("nshift must be positive");
    }
    for i in 0..ntr {
        // C computes tmin before reading this header, including on trace 0.
        let tmin = (tr.delrt as f64 / 1000.0) as f32;
        efread(tr.hdr.as_bytes_mut(), &headerfp);
        let mut val = Value::default();
        gethval(&tr, index, &mut val);
        let x = vtof(typ, val);
        let mut t = [0.0];
        intlin(&xs, &ts, tmin, ts[ts.len() - 1], &[x], &mut t);
        shifts[i] = nint(((t[0] - tmin) / dt) as f64);
        if i == 0 {
            shiftmax = shifts[i];
            shiftmin = shifts[i];
        } else {
            shiftmax = shiftmax.max(shifts[i]);
            shiftmin = shiftmin.min(shifts[i]);
        }
    }
    erewind(&headerfp);
    erewind(&tracefp);
    if verbose != 0 {
        for i in 0..ntr {
            warn!("inshift[%d]=%d", i, shifts[i]);
        }
    }
    let ntdshift = (nt as i32).wrapping_add(shiftmax);
    let shiftbytes = (ntdshift as usize).wrapping_mul(4);
    let mixbytes = shiftbytes.wrapping_mul(nmix);
    if verbose != 0 {
        warn!("nt=%d  shiftmax=%d  shiftmin=%d", nt, shiftmax, shiftmin);
        warn!(
            "ntdshift=%d  shiftbytes=%d  mixbytes=%d",
            ntdshift, shiftbytes, mixbytes
        );
    }
    // Negative dimensions reach impossible allocations in C. Preserve its
    // allocation error, after printing the sizes and before allocating.
    if ntdshift < 0 || (median != 0 && nmed < 0) {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // The C loop processes nt samples even when every shift is negative
    // and nt + shiftmax is smaller than nt. Allocate the full loop extent:
    // the minimal C allocation fix leaves the negative-shift null panel's
    // output unchanged (ub-2026-11, round 1), so it must remain supported.
    let ntd = nt.max(ntdshift as usize);
    let mut data = vec![vec![0.0f32; ntd]; nmix];
    let mut temp = vec![0.0f32; ntd];
    let mut dtemp = vec![0.0f32; nt];
    let mut stemp = vec![0.0f32; nmix];
    // Padding retains the previous port's zero extension for C's accesses
    // past the inshift vector while preloading and draining short panels.
    shifts.resize(ntr + nmix + 2, 0);
    let middle = (nmix - 1) / 2;
    let halfwidth = middle + 1;
    let mut trcount = 0;
    for i in 0..halfwidth {
        efread(&mut tr.data_bytes_mut()[..databytes], &tracefp);
        trcount = i;
        for it in 0..nt {
            let p = (it as i32)
                .wrapping_add(shiftmax)
                .wrapping_add(sign.wrapping_mul(shifts[i])) as usize;
            if p < ntd {
                data[middle - i][p] = tr.data[it];
                data[middle + i][p] = tr.data[it];
            }
        }
    }
    for i in 0..ntr {
        efread(tr.hdr.as_bytes_mut(), &headerfp);
        temp.fill(0.0);
        dtemp.fill(0.0);
        for it in 0..nt {
            if median == 0 {
                for j in 0..nmix {
                    temp[it] += data[j][it] * mix[j];
                }
            } else {
                for j in 0..nmix {
                    stemp[j] = data[j][it];
                }
                for j in 0..nmix - 1 {
                    for k in 0..nmix - 1 - j {
                        if stemp[k + 1] < stemp[k] {
                            stemp.swap(k, k + 1);
                        }
                    }
                }
                temp[it] = stemp[middle];
            }
            let p = (it as i32)
                .wrapping_sub(shiftmax)
                .wrapping_sub(sign.wrapping_mul(shifts[i]));
            if ((subtract != 0 && p >= 0) || (subtract == 0 && it as i32 - shiftmax >= 0))
                && p >= 0
                && (p as usize) < nt
            {
                dtemp[p as usize] = if subtract != 0 {
                    data[middle][it] - temp[it]
                } else {
                    temp[it]
                };
            }
        }
        tr.data[..nt].copy_from_slice(&dtemp);
        for j in (1..nmix).rev() {
            let (left, right) = data.split_at_mut(j);
            right[0].copy_from_slice(&left[j - 1]);
        }
        tr.ns = nt as u16;
        puttr(&mut tr);
        if trcount < ntr {
            efread(&mut tr.data_bytes_mut()[..databytes], &tracefp);
            trcount += 1;
            for it in 0..nt {
                let p = (it as i32)
                    .wrapping_add(shiftmax)
                    .wrapping_add(sign.wrapping_mul(shifts[trcount]))
                    as usize;
                if p < ntd {
                    data[0][p] = tr.data[it];
                }
            }
        } else {
            let r = 2 * (trcount - ntr);
            if r < nmix {
                data[0] = data[r].clone();
            }
            trcount += 1;
        }
    }
    if verbose != 0 && subtract != 0 {
        warn!("filtered data subtracted from input");
    }
    efclose(&headerfp);
    if !tmpdir.is_empty() {
        eremove(&headerfile);
    }
    efclose(&tracefp);
    if !tmpdir.is_empty() {
        eremove(&tracefile);
    }
    exit(CWP_Exit());
}
