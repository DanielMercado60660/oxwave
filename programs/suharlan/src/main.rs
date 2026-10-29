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

//! SUHARLAN - signal-noise separation by the invertible linear transformation
//! method of Harlan, 1984. Port of su/main/noise/suharlan.c (Revision 1.11).
//! Credits: Gabriel Alvarez CWP (1995). Some subroutines are direct
//! translations to C from Fortran versions written by Dr. Bill Harlan (1984).
//! The original source also carries: Copyright (c) Colorado School of Mines,
//! 2011. All rights reserved.

#![deny(unsafe_code)]

use oxwave_core::prelude::*;
use oxwave_correlation::xcor;
use oxwave_cwp::cwp::abs;
use oxwave_f_convolution::convolve_cwp;
use oxwave_f_frannor::{frannor, srannor};
use oxwave_f_franuni::{franuni, sranuni};
use oxwave_f_smooth::{dlsq_smoothing, smooth_histogram};
use oxwave_f_taup::{forward_p_transform, inverse_p_transform};
use oxwave_hilbert::hilbert;
use oxwave_index_search::xindex;
use oxwave_tempfile::temporary_filename;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUHARLAN - signal-noise separation by the invertible linear\t\t",
    "\t    transformation method of Harlan, 1984\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   suharlan <infile >outfile  [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t \t",
    " <none>\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " FLAGS:\t\t\t\t\t\t\t\t",
    " niter=1\tnumber of requested iterations\t\t\t\t",
    " anenv=1\t=1 for positive analytic envelopes\t\t\t",
    "\t\t=0 for no analytic envelopes (not recommended)\t\t",
    " scl=0\t\t=1 to scale output traces (not recommended)\t\t",
    " plot=3\t=0 for no plots. =1 for 1-D plots only\t\t\t",
    "\t\t=2 for 2-D plots only. =3 for all plots\t\t\t",
    " norm=1\t=0 not to normalize reliability values\t\t\t",
    " verbose=1\t=0 not to print processing information\t\t\t",
    " rgt=2\t\t=1 for uniform random generator\t\t\t\t",
    "\t\t=2 for gaussian random generator\t\t\t",
    " sts=1\t\t=0 for no smoothing (not recommended)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " General Parameters:\t\t\t\t\t\t\t",
    " dx=20\t\toffset sampling interval (m)\t\t\t\t",
    " fx=0\t  \toffset on first trace (m)\t\t\t\t",
    " dt=0.004\ttime sampling interval (s)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Tau-P Transform Parameters:\t\t\t\t\t\t",
    " gopt=1\t=1 for parabolic transform. =2 for Foster/Mosher\t",
    "\t\t=3 for linear. =4 for absolute value of linear\t\t",
    " pmin1=-400\tminimum moveout at farthest offset for fwd transf(ms)\t",
    " pmax1=400\tmaximum moveout at farthest offset for fwd transf(ms)\t",
    " pmin2=pmin1\tminimum moveout at farthest offset for inv transf(ms)\t",
    " pmax2=pmax1\tmaximum moveout at farthest offset for inv transf(ms)\t",
    " np=100\tnumber of p-values for taup transform\t\t\t",
    " prewhite=0.01\tprewhitening value (suggested between 0.1 and 0,01)\t",
    " offref=2000\treference offset for p-values (m)\t\t\t",
    " depthref=500\treference depth for Foster/Mosher taup (if gopt=4)\t",
    " pmula=pmax1\tmaximum p-value preserved in the data (ms)\t\t",
    " pmulb=pmax1\tminimum p-value muted on the data (ms)\t\t\t",
    " ninterp=0\tnumber of traces to interpolate in input data\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Extraction Parameters:\t\t\t\t\t\t",
    " nintlh=50\tnumber of intervals (bins) in histograms\t\t",
    " sditer=5\tnumber of steepest descent iterations to compute ps\t",
    " c=0.04\tmaximum noise allowed in a sample of signal(%)\t\t",
    " rel1=0.5\treliability value for first pass of the extraction\t",
    " rel2=0.75\treliability value for second pass of the extraction\t",
    "\t\t\t\t\t\t\t\t   \t",
    " Smoothing Parameters:\t\t\t\t\t\t\t",
    " r1=10\t\tnumber of points for damped lsq vertical smoothing\t",
    " r2=2\t\tnumber of points for damped lsq horizontal smoothing\t",
    "\t\t\t\t\t\t\t\t   \t",
    "\t\t\t\t\t\t\t\t   \t",
    " Output Files:\t\t\t\t\t\t\t\t",
    " signal=out_signal \tname of output file for extracted signal\t",
    " noise=out_noise \tname of output file for extracted noise\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The signal-noise separation algorithm was developed by Dr. Bill Harlan",
    " in 1984. It can be used to separate events that can be focused by a\t",
    " linear transformation (signal) from events that can't (noise). The\t",
    " linear transform is whatever is well siuted for the application at\t",
    " hand. Here, only the discrete Radon transform is used, so the program\t",
    " is capable of separating events focused by that transform (linear,\t",
    " parabolic or time-invariantly hyperbolic). Should other transform be\t",
    " required, the changes to the program will be relatively\t\t",
    " straightforward.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The reliability parameter is the most critical one to determine what\t",
    " to extract as signal and what to reject as noise. It should be tested\t",
    " for every dataset. The way to test it is to start with a small value,\t",
    " say 0.1 or 0.01. If too much noise is present in the extracted noise,\t",
    " it is too low. If too much signal was extracted, that is, part of the\t",
    " signal was lost, it is too big. All other parameters have good default",
    " values and should perhaps not be changed in a first encounter with the",
    " program. The transform parameters are also critical. They should be\t",
    " chosen such that no aliasing is present and such that the range of\t",
    " interesting slopes is spanned by the transform but not much more. The ",
    " program suradon.c has more documentation on the transform paramters.\t",
    "\t\t\t\t\t\t\t\t\t",
];

#[derive(Clone)]
struct Parameters {
    anenv: i32,
    scl: i32,
    plot: i32,
    norm: i32,
    verbose: i32,
    rgt: i32,
    sts: i32,
    gopt: i32,
    dx: f32,
    fx: f32,
    nintlh: i32,
    sditer: i32,
    c: f32,
    rel1: f32,
    rel2: f32,
    niter: i32,
    dt: f32,
    r1: f32,
    r2: f32,
    pmin1: f32,
    pmax1: f32,
    pmin2: f32,
    pmax2: f32,
    np: i32,
    prewhite: f32,
    offref: f32,
    interoff: f32,
    depthref: f32,
    pmula: f32,
    pmulb: f32,
    ninterp: i32,
    seed: i32,
}

fn main() {
    let argv = args();
    let mut tro1 = Trace::new();
    let mut tro2 = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    if gettr(&mut tro1) == 0 {
        err!("can't get first trace");
    }
    let nt = tro1.ns as usize;
    let mut dt: f32 = ((tro1.dt as f32) as f64 / 1000000.0) as f32;
    let mut anenv: i32 = 0;
    if getparint("anenv", &mut anenv) == 0 {
        anenv = 1;
    }
    let mut scl: i32 = 0;
    if getparint("scl", &mut scl) == 0 {
        scl = 0;
    }
    let mut plot: i32 = 0;
    if getparint("plot", &mut plot) == 0 {
        plot = 3;
    }
    let mut norm: i32 = 0;
    if getparint("norm", &mut norm) == 0 {
        norm = 1;
    }
    let mut verbose: i32 = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 1;
    }
    let mut rgt: i32 = 0;
    if getparint("rgt", &mut rgt) == 0 {
        rgt = 2;
    }
    let mut sts: i32 = 0;
    if getparint("sts", &mut sts) == 0 {
        sts = 1;
    }
    let mut gopt: i32 = 0;
    if getparint("gopt", &mut gopt) == 0 {
        gopt = 1;
    }
    let mut dx: f32 = 0.0;
    if getparfloat("dx", &mut dx) == 0 {
        dx = 20.0;
    }
    let mut fx: f32 = 0.0;
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.0;
    }
    let mut nintlh: i32 = 0;
    if getparint("nintlh", &mut nintlh) == 0 {
        nintlh = 50;
    }
    let mut sditer: i32 = 0;
    if getparint("sditer", &mut sditer) == 0 {
        sditer = 5;
    }
    let mut c: f32 = 0.0;
    if getparfloat("c", &mut c) == 0 {
        c = 0.05;
    }
    let mut rel1: f32 = 0.0;
    if getparfloat("rel1", &mut rel1) == 0 {
        rel1 = 0.5;
    }
    let mut rel2: f32 = 0.0;
    if getparfloat("rel2", &mut rel2) == 0 {
        rel2 = 0.75;
    }
    let mut niter: i32 = 0;
    if getparint("niter", &mut niter) == 0 {
        niter = 1;
    }
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.004;
    }
    let mut r1: f32 = 0.0;
    if getparfloat("r1", &mut r1) == 0 {
        r1 = 10.0;
    }
    let mut r2: f32 = 0.0;
    if getparfloat("r2", &mut r2) == 0 {
        r2 = 2.0;
    }
    let mut pmin1: f32 = 0.0;
    if getparfloat("pmin1", &mut pmin1) == 0 {
        pmin1 = -400.0;
    }
    let mut pmax1: f32 = 0.0;
    if getparfloat("pmax1", &mut pmax1) == 0 {
        pmax1 = 400.0;
    }
    let mut pmin2: f32 = 0.0;
    if getparfloat("pmin2", &mut pmin2) == 0 {
        pmin2 = pmin1;
    }
    let mut pmax2: f32 = 0.0;
    if getparfloat("pmax2", &mut pmax2) == 0 {
        pmax2 = pmax1;
    }
    let mut np: i32 = 0;
    if getparint("np", &mut np) == 0 {
        np = 100;
    }
    let mut prewhite: f32 = 0.0;
    if getparfloat("prewhite", &mut prewhite) == 0 {
        prewhite = 0.01;
    }
    let mut offref: f32 = 0.0;
    if getparfloat("offref", &mut offref) == 0 {
        offref = 2000.0;
    }
    let mut interoff: f32 = 0.0;
    if getparfloat("interoff", &mut interoff) == 0 {
        interoff = 0.0;
    }
    let mut depthref: f32 = 0.0;
    if getparfloat("depthref", &mut depthref) == 0 {
        depthref = 0.0;
    }
    let mut pmula: f32 = 0.0;
    if getparfloat("pmula", &mut pmula) == 0 {
        pmula = pmax1;
    }
    let mut pmulb: f32 = 0.0;
    if getparfloat("pmulb", &mut pmulb) == 0 {
        pmulb = pmax1;
    }
    let mut ninterp: i32 = 0;
    if getparint("ninterp", &mut ninterp) == 0 {
        ninterp = 0;
    }

    let mut signalfile = String::new();
    let mut noisefile = String::new();
    if getparstring("signal", &mut signalfile) == 0 {
        signalfile = "out_signal".into();
    }
    if getparstring("noise", &mut noisefile) == 0 {
        noisefile = "out_noise".into();
    }
    let mut seed = 0;
    if getparint("seed", &mut seed) == 0 {
        seed = match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
            Ok(now) => now.as_secs() as u32 as i32,
            Err(_) => -1,
        };
        if seed == -1 {
            seed = 1;
            warn!("time() failed to set seed, setting it to one");
        }
    }
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    if !tmpdir.is_empty() && !has_write_access(&tmpdir) {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }

    let (tracefp, headerfp, mut temporary) = if tmpdir.is_empty() {
        let tracefp = etmpfile();
        let headerfp = etmpfile();
        if verbose != 0 {
            warn!("using tmpfile() call");
        }
        (tracefp, headerfp, None)
    } else {
        let tracefile = temporary_filename(&tmpdir);
        let headerfile = temporary_filename(&tmpdir);
        let tracefp = efopen(&tracefile, "w+");
        let headerfp = efopen(&headerfile, "w+");
        if verbose != 0 {
            warn!("putting temporary files in %s", tmpdir);
        }
        (tracefp, headerfp, Some((tracefile, headerfile)))
    };

    let mut ntr = 0usize;
    loop {
        ntr += 1;
        fwrite(tro1.hdr.as_bytes(), &headerfp);
        fwrite(&tro1.data[..nt], &tracefp);
        if gettr(&mut tro1) == 0 {
            break;
        }
    }
    let mut traces = Matrix::new(nt, ntr);
    let mut out_signal = Matrix::new(nt, ntr);
    let mut out_noise = Matrix::new(nt, ntr);
    erewind(&tracefp);
    for ix in 0..ntr {
        fread(&mut traces[ix], &tracefp);
    }
    efclose(&tracefp);

    let p = Parameters {
        anenv,
        scl,
        plot,
        norm,
        verbose,
        rgt,
        sts,
        gopt,
        dx,
        fx,
        nintlh,
        sditer,
        c,
        rel1,
        rel2,
        niter,
        dt,
        r1,
        r2,
        pmin1,
        pmax1,
        pmin2,
        pmax2,
        np,
        prewhite,
        offref,
        interoff,
        depthref,
        pmula,
        pmulb,
        ninterp,
        seed,
    };
    let mut state = ExtractionState::default();
    separate_signal_noise(
        &p,
        nt,
        ntr,
        &traces,
        &mut out_signal,
        &mut out_noise,
        &mut state,
    );

    if !signalfile.is_empty() {
        let signal_file = match fopen(&signalfile, "w") {
            Some(fp) => fp,
            None => err!("cannot open signal file=%s\n", signalfile),
        };
        erewind(&headerfp);
        for itr in 0..ntr {
            fread(tro1.hdr.as_bytes_mut(), &headerfp);
            for it in 0..nt {
                tro1.data[it] = out_signal[itr][it];
            }
            fputtr(&signal_file, &mut tro1);
        }
        efclose(&signal_file);
    }
    if !noisefile.is_empty() {
        let noise_file = efopen(&noisefile, "w");
        erewind(&headerfp);
        for itr in 0..ntr {
            efread(tro2.hdr.as_bytes_mut(), &headerfp);
            for it in 0..nt {
                tro2.data[it] = out_noise[itr][it];
            }
            fputtr(&noise_file, &mut tro2);
        }
        efclose(&noise_file);
    }
    efclose(&headerfp);
    if let Some((tracefile, headerfile)) = temporary.take() {
        eremove(&headerfile);
        eremove(&tracefile);
    }
    exit(CWP_Exit());
}

// Call the same POSIX interface as main in suharlan.c. Permission mode bits
// alone cannot reproduce access(WRITE_OK): the real IDs and ACLs also matter.
// This binds the operating-system function; it does not implement an SU library.
#[allow(unsafe_code)]
fn has_write_access(path: &str) -> bool {
    unsafe extern "C" {
        fn access(path: *const std::ffi::c_char, mode: std::ffi::c_int) -> std::ffi::c_int;
    }
    let path = std::ffi::CString::new(path.split('\0').next().unwrap()).unwrap();
    // SAFETY: path is NUL-terminated and remains alive for the synchronous call.
    // WRITE_OK is 2 in su.h (POSIX W_OK).
    unsafe { access(path.as_ptr(), 2) == 0 }
}

// The original registers closefiles for INT, QUIT, HUP and TERM here.
// Rust std has no signal-handler registration API. Normal-exit cleanup is
// preserved above; signal cleanup requires a supported core API.

// The C golden_search and extract_signal functions have process-local static
// state. Carry it through both extraction passes.
#[derive(Default)]
struct ExtractionState {
    count: i32,
    golden: GoldenSearch,
}
#[derive(Default)]
struct GoldenSearch {
    ifill: usize,
    rmag: f32,
    x: [f32; 4],
    f: [f32; 4],
}

// Library matrix APIs use rows; plot_two_d writes them in C allocation order.
#[derive(Clone)]
struct Matrix {
    data: Vec<Vec<f32>>,
}
impl Matrix {
    fn new(n1: usize, n2: usize) -> Self {
        Self {
            data: vec![vec![0.0; n1]; n2],
        }
    }
}
impl std::ops::Index<usize> for Matrix {
    type Output = [f32];
    fn index(&self, i: usize) -> &[f32] {
        &self.data[i]
    }
}
impl std::ops::IndexMut<usize> for Matrix {
    fn index_mut(&mut self, i: usize) -> &mut [f32] {
        &mut self.data[i]
    }
}

fn separate_signal_noise(
    p: &Parameters,
    nt: usize,
    nx: usize,
    in_traces: &Matrix,
    signal: &mut Matrix,
    noise: &mut Matrix,
    state: &mut ExtractionState,
) {
    let np = p.np as usize;
    let mut rand = 0.0f32;
    let mut rel = 0.0f32;
    let mut plotname = "";
    if p.plot < 0 || p.plot > 3 {
        err!("plot flag has to be 0<=plot<=3");
    }
    let dp = (p.pmax1 - p.pmin1) / (p.np - 1) as f32;
    if p.verbose == 1 {
        warn!("nx=%d nt=%d\n", nx as i32, nt as i32);
    }
    let mut temp_traces_taup = Matrix::new(nt, np);
    let mut traces_taup = Matrix::new(nt, np);
    let mut noise_taup = Matrix::new(nt, np);

    forward_p_transform(
        nt,
        p.dt,
        p.pmax1,
        p.pmin1,
        dp,
        p.depthref,
        60.0,
        80.0,
        4.0,
        20.0,
        400,
        5,
        1,
        0,
        p.ninterp as usize,
        1,
        p.prewhite,
        p.interoff,
        p.offref,
        p.gopt,
        p.dx,
        p.fx,
        p.pmula,
        p.pmulb,
        &in_traces.data,
        &mut traces_taup.data,
    );
    if p.plot == 2 || p.plot == 3 {
        plotname = "input";
        plot_two_d(nx * nt, &in_traces.data, plotname);
        plotname = "input_taup";
        plot_two_d(nx * nt, &traces_taup.data, plotname);
    }
    for ix in 0..nx {
        for it in 0..nt {
            noise[ix][it] = in_traces[ix][it];
        }
    }
    for ip in 0..np {
        for it in 0..nt {
            temp_traces_taup[ip][it] = traces_taup[ip][it];
        }
    }
    for iter in 0..p.niter {
        if p.rgt == 1 {
            sranuni(p.seed);
        } else if p.rgt == 2 {
            srannor(p.seed);
        }
        for ix in 0..nx {
            if p.rgt == 1 {
                rand = (2.0 * franuni() as f64 - 1.0) as f32;
            } else if p.rgt == 2 {
                rand = frannor();
            }
            if rand < 0.0 {
                for it in 0..nt {
                    noise[ix][it] = (noise[ix][it] as f64 * -1.0) as f32;
                }
            }
        }
        forward_p_transform(
            nt,
            p.dt,
            p.pmax1,
            p.pmin1,
            dp,
            p.depthref,
            60.0,
            80.0,
            4.0,
            20.0,
            400,
            5,
            1,
            0,
            p.ninterp as usize,
            1,
            p.prewhite,
            p.interoff,
            p.offref,
            p.gopt,
            p.dx,
            p.fx,
            p.pmula,
            p.pmulb,
            &noise.data,
            &mut noise_taup.data,
        );
        if p.plot == 2 || p.plot == 3 {
            if iter == 0 {
                plotname = "rand1";
            } else if iter == 1 {
                plotname = "rand2";
            }
            plot_two_d(nx * nt, &noise.data, plotname);
            if iter == 0 {
                plotname = "rtaup1";
            } else if iter == 1 {
                plotname = "rtaup2";
            }
            plot_two_d(nx * nt, &noise_taup.data, plotname);
        }
        if iter == 0 {
            rel = p.rel1;
        } else if iter == 1 {
            rel = p.rel2;
        }
        extract_signal(p, rel, np, nt, &mut traces_taup, &noise_taup, state);
        inverse_p_transform(
            nt,
            p.dt,
            p.pmax2,
            p.pmin2,
            dp,
            p.depthref,
            60.0,
            80.0,
            p.interoff,
            p.offref,
            p.gopt,
            p.dx,
            p.fx,
            &traces_taup.data,
            &mut signal.data,
        );
        if p.plot == 2 || p.plot == 3 {
            if iter == 0 {
                plotname = "signal1";
            } else if iter == 1 {
                plotname = "signal2";
            }
            plot_two_d(nx * nt, &signal.data, plotname);
            if iter == 0 {
                plotname = "staup1";
            } else if iter == 1 {
                plotname = "staup2";
            }
            plot_two_d(nx * nt, &traces_taup.data, plotname);
        }
        if p.scl != 0 {
            scale_traces(p.scl, nx, nt, in_traces, signal);
            if p.plot == 2 {
                if iter == 0 {
                    plotname = "s1scaled";
                } else if iter == 1 {
                    plotname = "s2scaled";
                }
                plot_two_d(nx * nt, &signal.data, plotname);
            }
        }
        for ix in 0..nx {
            for it in 0..nt {
                if p.scl != -1 {
                    noise[ix][it] = in_traces[ix][it] - signal[ix][it];
                }
            }
        }
        for ip in 0..np {
            for it in 0..nt {
                traces_taup[ip][it] = temp_traces_taup[ip][it];
            }
        }
        if p.plot == 2 || p.plot == 3 {
            if iter == 0 {
                plotname = "noise1";
            } else if iter == 1 {
                plotname = "noise2";
            }
            plot_two_d(nx * nt, &noise.data, plotname);
        }
    }
}

fn extract_signal(
    p: &Parameters,
    r: f32,
    nx: usize,
    nt: usize,
    traces: &mut Matrix,
    rand_traces: &Matrix,
    state: &mut ExtractionState,
) {
    let nintlh = p.nintlh as usize;
    let mut rindex = 0;
    if p.verbose == 1 {
        warn!("PARAMETERS: c=%g rel=%g r1=%g r2=%g\n", p.c, r, p.r1, p.r2);
        warn!("FLAGS: anenv=%d plot=%d sts=%d\n", p.anenv, p.plot, p.sts);
    }
    let mut pd = vec![0.0f32; nintlh];
    let mut ps = vec![0.0f32; nintlh];
    let mut pn = vec![0.0f32; nintlh];
    let mut amps = vec![0.0f32; nintlh];
    let mut xamps = vec![0.0f32; nintlh];
    let mut pspn = vec![0.0f32; nintlh];
    let mut esd = vec![0.0f32; nintlh];
    let mut reliability = vec![0.0f32; nintlh];
    let mut anenv_traces = Matrix::new(nt, nx);
    if p.anenv != 0 {
        compute_analytic_envelopes(p.anenv, nx, nt, traces, &mut anenv_traces);
    }
    compute_histogram_stuff(
        p.verbose,
        nx,
        nt,
        nintlh,
        traces,
        rand_traces,
        &mut rindex,
        &mut amps,
        &mut xamps,
    );
    make_histogram(nintlh, nt, nx, traces, &amps, &mut pd);
    make_histogram(nintlh, nt, nx, rand_traces, &amps, &mut pn);
    smooth_histogram(&mut pn);
    deconvolve_histograms(
        nintlh,
        rindex,
        p.sditer,
        &pd,
        &pn,
        &mut ps,
        &mut state.golden,
    );
    convolve_cwp(-rindex, &ps, -rindex, &pn, -rindex, &mut pspn);
    expected_value(nintlh, rindex, &pspn, &pn, &ps, &xamps, &mut esd);
    compute_reliability(
        p.norm,
        nintlh,
        rindex,
        p.c,
        &ps,
        &pn,
        &pspn,
        &mut reliability,
    );
    zero_noisy_samples(
        p.anenv,
        p.sts,
        nx,
        nt,
        p.r1,
        p.r2,
        r,
        nintlh,
        &amps,
        traces,
        &anenv_traces,
        &reliability,
    );
    if p.verbose == 1 {
        let (mut spd, mut spn, mut sps, mut spspn) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
        for ix in 0..nintlh {
            spd += pd[ix];
            spn += pn[ix];
            sps += ps[ix];
            spspn += pspn[ix];
        }
        warn!("Area under the curve for pd=%g\n", spd);
        warn!("Area under the curve for pn=%g\n", spn);
        warn!("Area under the curve for ps=%g\n", sps);
        warn!("Area under the curve for ps*pn=%g\n", spspn);
    }
    state.count += 1;
    if p.plot == 1 || p.plot == 3 {
        let mut plotname = "";
        if state.count == 1 {
            plotname = "pd1";
        } else if state.count == 2 {
            plotname = "pd2";
        }
        plot_one_d(nintlh, &xamps, &pd, plotname);
        if state.count == 1 {
            plotname = "pn1";
        } else if state.count == 2 {
            plotname = "pn2";
        }
        plot_one_d(nintlh, &xamps, &pn, plotname);
        if state.count == 1 {
            plotname = "ps1";
        } else if state.count == 2 {
            plotname = "ps2";
        }
        plot_one_d(nintlh, &xamps, &ps, plotname);
        if state.count == 1 {
            plotname = "pspn1";
        } else if state.count == 2 {
            plotname = "pspn2";
        }
        plot_one_d(nintlh, &xamps, &pspn, plotname);
        if state.count == 1 {
            plotname = "Esd1";
        } else if state.count == 2 {
            plotname = "Esd2";
        }
        plot_one_d(nintlh, &xamps, &esd, plotname);
        if state.count == 1 {
            plotname = "rel1";
        } else if state.count == 2 {
            plotname = "rel2";
        }
        plot_one_d(nintlh, &xamps, &reliability, plotname);
    }
}

fn expected_value(
    nintlh: usize,
    rindex: i32,
    pspn: &[f32],
    pn: &[f32],
    ps: &[f32],
    xamps: &[f32],
    esd: &mut [f32],
) {
    let max = xamps[nintlh - 1] + (xamps[1] - xamps[0]);
    let mut xps = vec![0.0f32; nintlh];
    for ix in 0..nintlh {
        xps[ix] = xamps[ix] * ps[ix];
    }
    convolve_cwp(-rindex, &xps, -rindex, pn, -rindex, esd);
    for ix in 0..nintlh {
        esd[ix] = divide(max, esd[ix], pspn[ix]);
    }
}

fn compute_reliability(
    norm: i32,
    nintlh: usize,
    rindex: i32,
    c: f32,
    ps: &[f32],
    pn: &[f32],
    pspn: &[f32],
    reliability: &mut [f32],
) {
    let mut max = -999999999i32 as f32;
    let mut numerator = vec![0.0f32; nintlh];
    conv1(-rindex, ps, -rindex, pn, -rindex, &mut numerator, -1, c);
    for i in 0..nintlh {
        reliability[i] = divide(1.0, numerator[i], pspn[i]);
    }
    if norm == 1 {
        for i in 0..nintlh {
            if max < reliability[i] {
                max = reliability[i];
            }
        }
        for i in 0..nintlh {
            reliability[i] /= max;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn zero_noisy_samples(
    anenv: i32,
    sts: i32,
    nx: usize,
    nt: usize,
    r1: f32,
    r2: f32,
    r: f32,
    nintlh: usize,
    amps: &[f32],
    traces: &mut Matrix,
    anenv_traces: &Matrix,
    rel: &[f32],
) {
    let mut index = 0i32;
    let mut oz = Matrix::new(nt, nx);
    for ix in 0..nx {
        for it in 0..nt {
            let d = if anenv == 0 {
                traces[ix][it]
            } else {
                anenv_traces[ix][it]
            };
            xindex(&amps[..nintlh], d, &mut index);
            oz[ix][it] = if rel[index as usize] < r { 0.0 } else { 1.0 };
        }
    }
    if sts != 0 {
        dlsq_smoothing(0, nt, 0, nx, r1, r2, 0.0, &mut oz.data);
    }
    for ix in 0..nx {
        for it in 0..nt {
            traces[ix][it] *= oz[ix][it];
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn compute_histogram_stuff(
    verbose: i32,
    nx: usize,
    nt: usize,
    nintlh: usize,
    traces: &Matrix,
    rand_traces: &Matrix,
    rindex: &mut i32,
    amps: &mut [f32],
    xamps: &mut [f32],
) {
    let (min1, max1, sum1) = compute_max_min_sum(nx, nt, traces);
    let (min2, max2, _sum2) = compute_max_min_sum(nx, nt, rand_traces);
    let mut min = (min1 as f64 - 0.5 * (max2 - min2) as f64) as f32;
    let mut max = (max1 as f64 + 0.5 * (max2 - min2) as f64) as f32;
    let del_int = (max - min) / (nintlh as i32 - 1) as f32;
    let fint = min;
    min = (min as f64 - del_int as f64 / 2.0) as f32;
    max = (max as f64 + del_int as f64 / 2.0) as f32;
    for ix in 0..nintlh {
        amps[ix] = min + del_int * ix as f32;
        xamps[ix] = fint + del_int * ix as f32;
    }
    xindex(&amps[..nintlh], sum1 / (nt * nx) as f32, rindex);
    if verbose == 1 {
        warn!("for histograms: min=%g max=%g", min, max);
        warn!("first bin=%g binwidth=%g\n", fint, del_int);
    }
}

fn make_histogram(
    nintlh: usize,
    nt: usize,
    nx: usize,
    traces: &Matrix,
    amps: &[f32],
    pdf: &mut [f32],
) {
    let mut index = 0i32;
    let ns = (nt * nx) as i32;
    for i in 0..nintlh {
        pdf[i] = 0.0;
    }
    for ix in 0..nx {
        for it in 0..nt {
            xindex(&amps[..nintlh], traces[ix][it], &mut index);
            pdf[index as usize] += 1.0;
        }
    }
    for i in 0..nintlh {
        pdf[i] /= ns as f32;
    }
}

fn compute_analytic_envelopes(
    sgn: i32,
    nx: usize,
    nt: usize,
    traces: &Matrix,
    analytic_envelopes: &mut Matrix,
) {
    let mut hilbert_trace = vec![0.0f32; nt];
    for ix in 0..nx {
        hilbert(&traces[ix], &mut hilbert_trace);
        for it in 0..nt {
            let amp1 = hilbert_trace[it];
            let amp2 = traces[ix][it];
            let mut amp3 = amp1 * amp1 + amp2 * amp2;
            if sgn == 1 || sgn == 2 || sgn == 3 {
                amp3 = (amp3 as f64).sqrt() as f32;
            }
            if sgn == 3 || sgn == 6 {
                analytic_envelopes[ix][it] = -amp3;
            } else if sgn == 1 || sgn == 4 {
                if amp2 < 0.0 {
                    analytic_envelopes[ix][it] = -amp3;
                } else {
                    analytic_envelopes[ix][it] = amp3;
                }
            } else {
                analytic_envelopes[ix][it] = amp3;
            }
        }
    }
}

fn scale_traces(scl: i32, nx: usize, nt: usize, in_traces: &Matrix, out_traces: &mut Matrix) {
    if scl == 1 || scl == 3 {
        for ix in 0..nx {
            scale_one_trace(nt, &in_traces[ix], &mut out_traces[ix]);
        }
    }
    if scl == 2 || scl == 3 {
        let mut temp_in = Matrix::new(nx, nt);
        let mut temp_out = Matrix::new(nx, nt);
        matrix_transpose(nx, nt, in_traces, &mut temp_in);
        matrix_transpose(nx, nt, out_traces, &mut temp_out);
        for it in 0..nt {
            scale_one_trace(nx, &temp_in[it], &mut temp_out[it]);
        }
        matrix_transpose(nt, nx, &temp_out, out_traces);
    }
}

fn scale_one_trace(ns: usize, in_trace: &[f32], out_trace: &mut [f32]) {
    let num = dot_product(ns, in_trace, in_trace);
    let den = dot_product(ns, out_trace, out_trace);
    if den == 0.0 {
        return;
    }
    let scale = ((num / den) as f64).sqrt() as f32;
    for is in 0..ns {
        out_trace[is] *= scale;
    }
}

fn matrix_transpose(n1: usize, n2: usize, matrix: &Matrix, tr_matrix: &mut Matrix) {
    for i1 in 0..n1 {
        for i2 in 0..n2 {
            tr_matrix[i2][i1] = matrix[i1][i2];
        }
    }
}

fn dot_product(ns: usize, vector1: &[f32], vector2: &[f32]) -> f32 {
    let mut sum = 0.0f32;
    for is in 0..ns {
        sum += vector1[is] * vector2[is];
    }
    sum
}

fn compute_max_min_sum(nx: usize, nt: usize, data: &Matrix) -> (f32, f32, f32) {
    let mut min = 9999999.0f32;
    let mut max = -9999999.0f32;
    let mut sum = 0.0f32;
    for ix in 0..nx {
        for it in 0..nt {
            let d = data[ix][it];
            if d < min {
                min = d;
            }
            if d > max {
                max = d;
            }
            sum += d;
        }
    }
    (min, max, sum)
}

fn plot_one_d(npoints: usize, xamps: &[f32], data: &[f32], plotname: &str) {
    let out = fopen(plotname, "w").expect("plot_one_d fopen returned NULL");
    fwrite(&xamps[..npoints], &out);
    fclose(&out);
    let out = fopen(plotname, "a").expect("plot_one_d fopen returned NULL");
    fwrite(&data[..npoints], &out);
    fclose(&out);
}
fn plot_two_d(npoints: usize, data: &[Vec<f32>], plotname: &str) {
    // The caller uses nx*nt even for tau-p arrays (whose size is np*nt).
    // Do not silently shorten a C out-of-bounds plot when np < nx.
    assert!(data.iter().map(Vec::len).sum::<usize>() >= npoints);
    let out = fopen(plotname, "w").expect("plot_two_d fopen returned NULL");
    let mut remaining = npoints;
    for row in data {
        let count = remaining.min(row.len());
        fwrite(&row[..count], &out);
        remaining -= count;
        if remaining == 0 {
            break;
        }
    }
    fclose(&out);
}

// This helper belongs to suharlan.c (it is distinct from convolve_cwp).
// Direct C translation of Dr. Bill Harlan's Fortran version, 1984.
#[allow(clippy::too_many_arguments)]
fn conv1(fx: i32, x: &[f32], fy: i32, y: &[f32], fz: i32, z: &mut [f32], flag: i32, perc: f32) {
    let nx = x.len() as i32;
    let ny = y.len() as i32;
    let nz = z.len() as i32;
    for iz in 0..nz {
        z[iz as usize] = 0.0;
    }
    for iz in 1..=nz {
        let jz = iz + fz - 1;
        for ix in 1..=nx {
            let jx = ix + fx - 1;
            let mut jy = jz - jx;
            if flag == 1 {
                jy = jz + jx;
            }
            let iy = jy - fy + 1;
            if iy < 1 || iy > ny {
                continue;
            }
            let mut s = 1.0f32;
            if flag == -1 {
                let mut r = jz as f32;
                r = abs(r) * perc;
                let mut rl = (jy as f64 - 0.5) as f32;
                let mut rr = (jy as f64 + 0.5) as f32;
                if rl > r {
                    rl = r;
                }
                if rl < -r {
                    rl = -r;
                }
                if rr > r {
                    rr = r;
                }
                if rr < -r {
                    rr = -r;
                }
                s = rr - rl;
            }
            z[(iz - 1) as usize] += s * x[(ix - 1) as usize] * y[(iy - 1) as usize];
        }
    }
}

// Gabriel Alvarez (1995), translated from Dr. Bill Harlan's FORTRAN IV (1982).
fn deconvolve_histograms(
    nsamples: usize,
    mean_index: i32,
    niter: i32,
    pd: &[f32],
    pn: &[f32],
    ps: &mut [f32],
    golden: &mut GoldenSearch,
) {
    let ndiv = 20;
    let rmax = 1000.0f32;
    let mut alpha = 0.5f32;
    let mut xvalue = f32::MIN_POSITIVE;
    let mut fvalue = 0.0f32;
    let mut agold = 0.0f32;
    let mut psp = vec![0.0f32; nsamples];
    let mut grad = vec![0.0f32; nsamples];
    for j in 0..nsamples {
        ps[j] = 0.0;
    }
    ps[mean_index as usize] = 1.0;
    for _iter in 0..niter {
        let mut scale = 2.0f32 * alpha;
        if alpha == 0.0 {
            scale = 2.0f32 * xvalue;
        }
        if scale > 1.0 {
            scale = 1.0;
        }
        gradient(nsamples, mean_index, rmax, pd, pn, ps, &mut grad);
        let mut idiv = 1;
        while idiv > 0 {
            let xgold = golden_search(fvalue, &mut idiv, &mut agold, golden);
            xvalue = xgold * scale;
            alpha = agold * scale;
            if idiv >= ndiv {
                break;
            }
            let r1 = (1.0 - xvalue as f64) as f32;
            let r2 = xvalue;
            let mut sum = 0.0f32;
            for j in 0..nsamples {
                psp[j] = r1 * ps[j] + r2 * grad[j];
                sum += psp[j];
            }
            for j in 0..nsamples {
                psp[j] /= abs(sum);
            }
            fvalue = cross_entropy(nsamples, mean_index, rmax, pd, pn, &psp);
        }
        let r1 = (1.0 - alpha as f64) as f32;
        let r2 = alpha;
        let mut sum = 0.0f32;
        for j in 0..nsamples {
            ps[j] = r1 * ps[j] + r2 * grad[j];
            sum += ps[j];
        }
        for j in 0..nsamples {
            ps[j] /= abs(sum);
        }
    }
}

// Gabriel Alvarez (1995), translated from Dr. Bill Harlan's Fortran IV (1983).
fn gradient(ns: usize, si: i32, rmax: f32, pd: &[f32], pn: &[f32], ps: &[f32], grad: &mut [f32]) {
    let mut pspn = vec![0.0f32; ns];
    convolve_cwp(-si, pn, -si, ps, -si, &mut pspn);
    for j in 0..ns {
        pspn[j] = divide(rmax, pd[j], pspn[j]);
    }
    xcor(-si, pn, -si, &pspn, -si, grad);
    let mut sum = 0.0f32;
    for j in 0..ns {
        sum += grad[j];
    }
    for j in 0..ns {
        grad[j] /= abs(sum);
    }
}

// Gabriel Alvarez (1995), translated from Dr. Bill Harlan's Fortran IV (1983).
fn cross_entropy(ns: usize, si: i32, rmax: f32, pd: &[f32], pn: &[f32], ps: &[f32]) -> f32 {
    let mut sum = 0.0f32;
    let mut pspn = vec![0.0f32; ns];
    convolve_cwp(-si, ps, -si, pn, -si, &mut pspn);
    for j in 0..ns {
        pspn[j] = divide(rmax, pspn[j], pd[j]);
    }
    for j in 0..ns {
        pspn[j] = divide(rmax, 1.0, pspn[j]);
        sum = (sum as f64 + pd[j] as f64 * (pspn[j] as f64).ln()) as f32;
    }
    sum
}

// Gabriel Alvarez (1995), translated from Dr. Bill Harlan's Fortran IV (1983).
fn divide(mut rmax: f32, mut a: f32, mut b: f32) -> f32 {
    if abs(a) as f64 > 1000000000.0 && abs(b) as f64 > 1000000000.0 {
        a = (a as f64 / 1000000000.0) as f32;
        b = (b as f64 / 1000000000.0) as f32;
    }
    if abs(a) as f64 > 1000000000.0 && abs(b) as f64 > 1000000000.0 {
        a = (a as f64 / 1000000000.0) as f32;
        b = (b as f64 / 1000000000.0) as f32;
    }
    rmax = abs(rmax);
    let mut quotient;
    if abs(a) < rmax * abs(b) {
        quotient = a / b;
        if abs(quotient) > rmax {
            quotient = 0.0;
        }
    } else {
        quotient = rmax;
        if (a > 0.0 && b < 0.0) || (a < 0.0 && b > 0.0) {
            quotient = (quotient as f64 * -1.0) as f32;
        }
    }
    if a == 0.0 {
        quotient = 0.0;
    }
    quotient
}

// Gabriel Alvarez (1995), translated from Dr. Bill Harlan's Fortran IV (1983).
fn golden_search(fvalue: f32, iter: &mut i32, alpha: &mut f32, s: &mut GoldenSearch) -> f32 {
    if *iter < 1 {
        *iter = 1;
    }
    if *iter == 1 {
        *alpha = 0.0;
        s.rmag = (0.5 * (5.0f64.sqrt() - 1.0)) as f32;
        s.x[0] = 0.0;
        s.x[1] = (1.0 - s.rmag as f64) as f32;
        s.x[2] = s.rmag;
        s.x[3] = 1.0;
        s.f[s.ifill] = fvalue;
        s.ifill = (*iter - 1) as usize;
    } else if *iter <= 4 {
        s.f[s.ifill] = fvalue;
        s.ifill = (*iter - 1) as usize;
    } else {
        s.f[s.ifill] = fvalue;
        let mut im = 0;
        let mut fm = s.f[0];
        for j in 1..4 {
            if s.f[j] >= fm {
                continue;
            }
            im = j;
            fm = s.f[j];
        }
        *alpha = s.x[im];
        if im == 0 || im == 1 {
            s.x[3] = s.x[2];
            s.f[3] = s.f[2];
            s.x[2] = s.x[1];
            s.f[2] = s.f[1];
            s.x[1] = (s.x[2] - s.x[0]) * s.rmag + s.x[0];
            s.ifill = 1;
        } else {
            s.x[0] = s.x[1];
            s.f[0] = s.f[1];
            s.x[1] = s.x[2];
            s.f[1] = s.f[2];
            s.x[2] = s.x[3] - (s.x[3] - s.x[1]) * s.rmag;
            s.ifill = 2;
        }
    }
    let xvalue = s.x[s.ifill];
    *iter += 1;
    xvalue
}
