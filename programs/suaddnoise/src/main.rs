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

//! SUADDNOISE - add noise to traces.
//! Port of su/main/noise/suaddnoise.c (Revision 1.48).
//! Credits: CWP: Jack Cohen, Brian Sumner, Ken Larner;
//! John Stockwell (fixed filtered noise option).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, max};
use oxwave_f_frannor::{frannor, srannor};
use oxwave_f_franuni::{franuni, sranuni};
use oxwave_tempfile::temporary_filename;
use std::ffi::CString;
use std::sync::atomic::{AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUADDNOISE - add noise to traces\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suaddnoise <stdin >stdout  sn=20  noise=gauss  seed=from_clock\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tif any of f=f1,f2,... and amp=a1,a2,... are specified by the user",
    "\tand if dt is not set in header, then dt is mandatory\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " \tsn=20\t\t\tsignal to noise ratio\t\t\t",
    " \tnoise=gauss\t\tnoise probability distribution\t\t",
    " \t\t\t\t=flat for uniform; default Gaussian\t",
    " \tseed=from_clock\t\trandom number seed (integer)\t\t",
    "\tf=f1,f2,...\t\tarray of filter frequencies (as in sufilter)",
    "\tamps=a1,a2,...\t\tarray of filter amplitudes\t\t",
    " \tdt= (from header)\ttime sampling interval (sec)\t\t",
    "\tverbose=0\t\t=1 for echoing useful information\t",
    " \t\t\t\t\t\t\t\t\t",
    " \ttmpdir=\t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Output = Signal +  scale * Noise\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " scale = (1/sn) * (absmax_signal/sqrt(2))/sqrt(energy_per_sample)\t",
    " \t\t\t\t\t\t\t\t\t",
    " If the signal is already band-limited, f=f1,f2,... and amps=a1,a2,...\t",
    " can be used, as in sufilter, to bandlimit the noise traces to match\t",
    " the signal band prior to computing the scale defined above.\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Examples of noise bandlimiting:\t\t\t\t\t",
    " low freqency:    suaddnoise < data f=40,50 amps=1,0 | ...\t\t",
    " high freqency:   suaddnoise < data f=40,50 amps=0,1 | ...\t\t",
    " near monochromatic: suaddnoise < data f=30,40,50 amps=0,1,0 | ...\t",
    " with a notch:    suaddnoise < data f=30,40,50 amps=1,0,1 | ...\t",
    " bandlimited:     suaddnoise < data f=20,30,40,50 amps=0,1,1,0 | ...\t",
    " \t\t\t\t\t\t\t\t\t",
];

const GAUSS: i32 = 0;
const FLAT: i32 = 1;

extern "C" {
    fn access(path: *const std::ffi::c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

#[derive(Default)]
struct TempFiles {
    tracefile: String,
    headerfile: String,
    bandoutfile: String,
    tracefp: Option<CFile>,
    headerfp: Option<CFile>,
    bandoutfp: Option<CFile>,
}

static FILES: AtomicPtr<TempFiles> = AtomicPtr::new(std::ptr::null_mut());

fn trap_signals() {
    // POSIX SIGHUP, SIGINT, SIGQUIT, SIGTERM, as in the C program.
    for sig in [2, 3, 1, 15] {
        unsafe {
            signal(sig, closefiles);
        }
    }
}

// The C helper also uses stdio and remove from its signal handler.
extern "C" fn closefiles(_sig: i32) {
    // The pointee is leaked by main, before any handlers are installed.
    let files = unsafe { &*FILES.load(Ordering::Relaxed) };
    if let Some(fp) = &files.headerfp {
        efclose(fp);
    }
    if let Some(fp) = &files.tracefp {
        efclose(fp);
    }
    eremove(&files.headerfile);
    eremove(&files.tracefile);
    if let Some(fp) = &files.bandoutfp {
        efclose(fp);
    }
    eremove(&files.bandoutfile);
    exit(EXIT_FAILURE);
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get parameters */
    let mut verbose = 0;
    getparint("verbose", &mut verbose);

    /* Look for user-supplied tmpdir */
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() {
        let path = CString::new(tmpdir.as_str()).unwrap();
        if unsafe { access(path.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }

    /* Get noise type; recall itype initialized as GAUSS */
    let mut stype = String::new();
    if getparstring("noise", &mut stype) == 0 {
        stype = "gauss".into();
    }
    let mut itype = GAUSS;
    if stype == "flat" {
        itype = FLAT;
    } else if stype != "gauss" {
        err!("noise=\"%s\", must be gauss or flat", stype);
    }

    /* Get signal to noise ratio */
    let mut sn = 20.0f32;
    getparfloat("sn", &mut sn);
    if sn <= 0.0 {
        // C erroneously formats a promoted double with %d. Apple arm64 varargs
        // read its low 32 bits; preserve even the negative fractional cases.
        err!("sn=%d must be positive", (sn as f64).to_bits() as i32);
    }

    /* Set seed */
    let mut seed = 0u32;
    if getparuint("seed", &mut seed) == 0 {
        match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
            Ok(seed_time) => seed = seed_time.as_secs() as u32,
            Err(_) => err!("time() failed to set seed"),
        }
    }
    if itype == GAUSS {
        srannor(seed as i32);
    } else {
        sranuni(seed as i32);
    }

    /* Prepare temporary files to hold headers and data */
    let files = Box::leak(Box::new(TempFiles::default()));
    FILES.store(files as *mut TempFiles, Ordering::Relaxed);
    let mut istmpdir = false;
    if tmpdir.is_empty() {
        files.tracefp = Some(etmpfile());
        files.headerfp = Some(etmpfile());
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
    } else {
        files.tracefile = temporary_filename(&tmpdir);
        files.headerfile = temporary_filename(&tmpdir);
        trap_signals();
        files.tracefp = Some(efopen(&files.tracefile, "w+"));
        files.headerfp = Some(efopen(&files.headerfile, "w+"));
        istmpdir = true;
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
    }
    let tracefp = files.tracefp.as_ref().unwrap().clone();
    let headerfp = files.headerfp.as_ref().unwrap().clone();

    /* Get info from first trace */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let databytes = nt as usize * std::mem::size_of::<f32>();

    /* Loop over input traces & write headers and data to tmp files */
    let mut ntr = 0i32;
    loop {
        ntr = ntr.wrapping_add(1);
        efwrite(tr.hdr.as_bytes(), &headerfp);
        efwrite(&tr.data_bytes()[..databytes], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    let nfloats = ntr.wrapping_mul(nt);

    /* Compute absmax of signal over entire data set */
    erewind(&tracefp);
    let mut absmaxsig = 0.0f32;
    for _i in 0..nfloats {
        let mut sigval = [0.0f32];
        efread(&mut sigval, &tracefp);
        absmaxsig = max(absmaxsig, abs(sigval[0]));
    }

    /* Compute noise vector elements in [-1, 1] */
    let mut noise = vec![0.0f32; nfloats as usize];
    match itype {
        GAUSS => {
            /* frannor gives elements in N(0,1)--ie. pos & negs */
            for i in 0..nfloats {
                noise[i as usize] = frannor();
            }
        }
        FLAT => {
            /* franuni gives elements in [0, 1] */
            for i in 0..nfloats {
                noise[i as usize] = (2.0 * franuni() as f64 - 1.0) as f32;
            }
        }
        _ => err!("%d: mysterious itype = %d", 214, itype),
    }

    /* Band limit noise traces if user getpars any of the f's */
    let mut f = String::new();
    let mut amps = String::new();
    getparstring("f", &mut f);
    getparstring("amps", &mut amps);
    if !f.is_empty() || !amps.is_empty() {
        let nsegy = HDRBYTES + databytes;
        let mut segybuf = vec![0u8; nsegy];

        trap_signals();

        /* Prepare temporary files to hold traces */
        let bandinfp = etmpfile();
        // temporary_filename returns a new name without modifying its prefix.
        // Preserve the C bug: bandoutfile remains empty for command and remove.
        let bandoutfp = efopen(&temporary_filename(&files.bandoutfile), "w+");
        files.bandoutfp = Some(bandoutfp.clone());

        /* Paste headers on noise traces and put in tmpfile */
        erewind(&headerfp);
        for itr in 0..ntr {
            efread(tr.hdr.as_bytes_mut(), &headerfp);
            let trshift = itr.wrapping_mul(nt) as usize;
            tr.data[..nt as usize].copy_from_slice(&noise[trshift..trshift + nt as usize]);
            fputtr(&bandinfp, &mut tr);
        }

        /* Build cmdbuf; append sufilter command */
        let mut cmdbuf = sprintf!("sufilter >%s", files.bandoutfile);
        for arg in &argv[1..] {
            /* Don't pass sn=, seed=, noise= */
            if !arg.starts_with("sn=") && !arg.starts_with("seed=") && !arg.starts_with("noise=") {
                cmdbuf.push(' ');
                cmdbuf.push_str(arg);
            }
        }

        let fp = epopen(&cmdbuf, "w");
        erewind(&bandinfp);
        for _itr in 0..ntr {
            efread(&mut segybuf, &bandinfp);
            efwrite(&segybuf, &fp);
        }
        efclose(&bandinfp);
        epclose(&fp);

        /* Load bandlimited traces back into noise vector */
        erewind(&bandoutfp);
        for itr in 0..ntr {
            fgettr(&bandoutfp, &mut tr);
            let trshift = itr.wrapping_mul(nt) as usize;
            noise[trshift..trshift + nt as usize].copy_from_slice(&tr.data[..nt as usize]);
        }
        efclose(&bandoutfp);
        eremove(&files.bandoutfile);
    }

    /* Compute noise power */
    let mut noipow = 0.0f32;
    for i in 0..nfloats {
        let noiseval = noise[i as usize];
        noipow += noiseval * noiseval;
    }

    /* Compute noise scale for desired noise/signal ratio */
    absmaxsig = (absmaxsig as f64 / 2.0f64.sqrt()) as f32;
    noipow /= nfloats as f32;
    let noiscale = if absmaxsig != 0.0 {
        (absmaxsig as f64 / (sn as f64 * (noipow as f64).sqrt())) as f32
    } else {
        1.0f32
    };

    /* Add scaled noise to trace and output sum */
    erewind(&headerfp);
    erewind(&tracefp);
    let mut output_bytesper = 4usize;
    for itr in 0..ntr {
        let trshift = itr.wrapping_mul(nt) as usize;
        efread(tr.hdr.as_bytes_mut(), &headerfp);
        efread(&mut tr.data_bytes_mut()[..databytes], &tracefp);
        for i in 0..nt as usize {
            tr.data[i] += noiscale * noise[trshift + i];
        }
        if itr == 0 {
            output_bytesper = match tr.trid {
                201 => 1,
                202 => 2,
                _ => 4,
            };
        }
        if output_bytesper != 4 && tr.trid != 201 && tr.trid != 202 {
            // fputtr.c:123-130 encodes all floats, but writes only the byte
            // count saved from the first (packed) trace. The shared Rust
            // writer writes all encoded bytes and then reports an error.
            // Keep SU's truncated output here, including a partial float.
            let fp = stdout();
            efwrite(&tr.hdr.to_xdr(), &fp);
            let mut encoded = Vec::with_capacity(databytes);
            for i in 0..nt as usize {
                encoded.extend_from_slice(&tr.data[i].to_bits().to_be_bytes());
            }
            efwrite(&encoded[..output_bytesper * nt as usize], &fp);
        } else {
            puttr(&mut tr);
        }
    }

    /* Clean up */
    efclose(&headerfp);
    if istmpdir {
        eremove(&files.headerfile);
    }
    efclose(&tracefp);
    if istmpdir {
        eremove(&files.tracefile);
    }
    exit(CWP_Exit());
}
