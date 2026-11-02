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

//! SUTIFOWLER - VTI constant velocity prestack time migration.
//! Port of su/main/migration_inversion/sutifowler.c (Revision 1.5).
//! Credits: CWP: John Anderson (visitor to CSM from Mobil) Spring 1993.
//! Taper based on code by David Hale.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
use oxwave_cwp::Complex;

use oxwave_f_vnd::*;
use oxwave_linear_interp::*;
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::*;
use std::ffi::CStr;

const NP: i32 = 10003;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUTIFOWLER   VTI constant velocity prestack time migration\t\t",
    "\t      velocity analysis via Fowler's method\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " sutifowler ncdps=250 vmin=1500 vmax=6000 dx=12.5\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameter:\t\t\t\t\t\t\t",
    " ncdps=\t\tnumber of input cdp's\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " choose=1\t1 do full prestack time migration\t\t\t",
    "\t\t2 do only DMO\t\t\t\t\t\t",
    "\t\t3 do only post-stack migrations\t\t\t\t",
    "\t\t4 do only stacking velocity analysis\t\t\t",
    " getcvstacks=0\tflag to set to 1 if inputting precomputed cvstacks\t",
    "\t\t(vmin, nvstack, and ncdps must match SUCVS4FOWLER job)\t",
    " vminstack=vmin\tminimum velocity panel in m/s in input cvstacks\t",
    " etamin=0.\t\tminimum eta (see paper by Tariq Alkhalifah)\t",
    " etamax=0.5\tmaximum eta (see paper by Tariq Alkhalifah)\t\t",
    " neta=1\tnumber of eta values to image\t\t\t\t",
    " d=0.\t\tThomsen's delta\t\t\t\t\t\t",
    " vpvs=0.5\tassumed vp/vs ratio (not critical -- default almost always ok)",
    " dx=25.\tcdp x increment\t\t\t\t\t\t",
    " vmin=1500.\tminimum velocity panel in m/s to output\t\t\t",
    " vmax=8000.\tmaximum velocity panel in m/s to output\t\t\t",
    " nv=75\t number of velocity panels to output\t\t\t\t",
    " nvstack=180\tnumber of stacking velocity panels to compute\t\t",
    "\t\t     ( Let offmax be the maximum offset, fmax be\t",
    "\t\t     the maximum freq to preserve, and tmute be\t\t",
    "\t\t     the starting mute time in sec on offmax, then\t",
    "\t\t     the recommended value for nvstack would be\t\t",
    "\t\t     nvstack = 4 +(offmax*offmax*fmax)/(0.6*vmin*vmin*tmute)",
    "\t\t     ---you may want to make do with less---)\t\t",
    " nxpad=0\t  number of traces to padd for spatial fft\t\t",
    "\t\t     Ideally nxpad = (0.5*tmax*vmax+0.5*offmax)/dx\t",
    " lmute=24\t length of mute taper in ms\t\t\t\t",
    " nonhyp=1\t  1 if do mute at 2*offset/vmin to avoid non-hyperbolic ",
    "\t\t\t\tmoveout, 0 otherwise\t\t\t",
    " lbtaper=0\tlength of bottom taper in ms\t\t\t\t",
    " lstaper=0\tlength of side taper in traces\t\t\t\t",
    " dtout=1.5*dt\toutput sample rate in s,   note: typically\t\t",
    "\t\t\t\tfmax=salias*0.5/dtout\t\t\t",
    " mxfold=120\tmaximum number of offsets/input cmp\t\t\t",
    " salias=0.8\tfraction of output frequencies to force within sloth\t",
    "\t\t     antialias limit.  This controls muting by offset of",
    "\t\t     the input data prior to computing the cv stacks\t",
    "\t\t     for values of choose=1 or choose=2.\t\t",
    " file=sutifowler\troot name for temporary files\t\t\t",
    " p=not\t\tEnter a path name where to put temporary files.\t\t",
    "\t  \tspecified  Can enter multiple times to use multiple disk",
    "\t\tsystems.\t\t\t\t\t\t",
    "\t\t     The default uses information from the .VND file\t",
    "\t\t     in the current directory if it exists, or puts \t",
    "\t\t     unique temporary files in the current directory.\t",
    " ngroup=20\tNumber of cmps per velocity analysis group.\t\t",
    " printfile=stderr    The output file for printout from this program.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required trace header words on input are ns, dt, cdp, offset.\t\t",
    " On output, trace headers are rebuilt from scratch with\t\t",
    " ns - number of samples\t\t\t\t\t\t",
    " dt - sample rate in usec\t\t\t\t\t\t",
    " cdp - the output cmp number (0 based)\t\t\t\t\t",
    " offset - the output velocity\t\t\t\t\t\t",
    " tracf\t- the output velocity index (0 based)\t\t\t\t",
    " fldr - index for velocity analysis group (0 based, groups of ngroup cdps)",
    " ep - central cmp for velocity analysis group\t\t\t\t",
    " igc - index for choice of eta (0 based)\t\t\t\t",
    " igi - eta*100\t\t\t\t\t\t\t\t",
    " sx=gx\t- x coordinate as icmp*dx\t\t\t\t\t",
    " tracl=tracr -sequential trace count (1 based)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: Due to aliasing considerations, the small offset-to-depth\t",
    " ratio assumption inherent in the TI DMO derivation, and the\t\t",
    " poor stacking of some large-offset events associated with TI non-hyperbolic",
    " moveout, a fairly stiff initial mute is recommended for the\t\t",
    " long offsets.  As a result, this method may not work well\t\t",
    " where you have multiple reflections to remove via stacking.\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Note: The temporary files can be split over multiple disks by building",
    " a .VND file in your working directory.  The .VND file is ascii text\t",
    " with the first line giving the number of directories followed by\t",
    " successive lines with one line per directory name.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: The output data order has primary key equal to cdp, secondary\t",
    " key equal to eta, and tertiary key equal to velocity.\t\t\t",
];

// Keep parameter parsing and numerical operations in the original C order.
#[allow(unused_variables, unused_assignments)]
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get first trace and extract critical info from header */
    if gettr(&mut tr) == 0 {
        err!("Can't get first trace \n");
    }
    let nt = tr.ns as i64;
    let dt = (0.000001 * tr.dt as f64) as f32;
    let mut oldcmp = tr.cdp;
    let mut printfile = String::new();
    let fpl = if getparstring("printfile", &mut printfile) == 0 {
        Some(stderr())
    } else {
        fopen(&printfile, "w")
    };

    let mut salias = 0.8f32;
    getparfloat("salias", &mut salias);
    if salias > 1.0 {
        salias = 1.0;
    }
    let mut dtout = (1.5 * dt as f64) as f32;
    getparfloat("dtout", &mut dtout);
    let mut ntout = (1.0f32 + nt as f32 * dt / dtout) as i64;
    let mut getcvstacks = 0i32;
    getparint("getcvstacks", &mut getcvstacks);
    if getcvstacks != 0 {
        dtout = dt;
        ntout = nt;
    }
    let fmax = (salias as f64 * 0.5 / dtout as f64) as f32;
    let fpl = match fpl {
        Some(fp) => fp,
        // The C program dereferences an unchecked fopen result at the first
        // fprintf, after parsing salias, dtout, and getcvstacks.
        None => null_file(),
    };
    fprintf!(fpl, "sutifowler: ntin=%ld dtin=%f\n", nt, dt);
    fprintf!(fpl, "sutifowler: ntout=%ld dtout=%f\n", ntout, dtout);
    let mut file = "sutifowler".to_string();
    getparstring("file", &mut file);
    let mut dx = 25.0f32;
    getparfloat("dx", &mut dx);
    let mut vmin = 1500.0f32;
    getparfloat("vmin", &mut vmin);
    let mut vmax = 8000.0f32;
    getparfloat("vmax", &mut vmax);
    let mut vminstack = vmin;
    getparfloat("vminstack", &mut vminstack);
    let mut d = 0.0f32;
    getparfloat("d", &mut d);
    let mut etamin = 0.0f32;
    getparfloat("etamin", &mut etamin);
    let mut etamax = 0.5f32;
    getparfloat("etamax", &mut etamax);
    let mut vsvp = 0.5f32;
    getparfloat("vsvp", &mut vsvp);
    let mut neta = 1i32;
    getparint("neta", &mut neta);
    if ((etamax - etamin) as f64).abs() < 1.0e-7 {
        neta = 1;
    }
    if neta < 1 {
        neta = 1;
    }
    let mut ichoose = 1i32;
    getparint("choose", &mut ichoose);
    let mut ncmps = 0i32;
    if getparint("ncdps", &mut ncmps) == 0 {
        err!("sutifowler: must enter ncdps");
    }
    let mut nv = 75i32;
    getparint("nv", &mut nv);
    let mut nvstack = 180i32;
    getparint("nvstack", &mut nvstack);
    let mut ntpad = (0.1 * ntout as f64) as i32;
    getparint("ntpad", &mut ntpad);
    let mut nxpad = 0i32;
    getparint("nxpad", &mut nxpad);
    let mut lmute = 24i32;
    getparint("lmute", &mut lmute);
    lmute = (1.0 + 0.001 * lmute as f64 / dtout as f64) as i32;
    let mut lbtaper = 0i32;
    getparint("lbtaper", &mut lbtaper);
    let mut lstaper = 0i32;
    getparint("lstaper", &mut lstaper);
    let mut mxfold = 120i32;
    getparint("mxfold", &mut mxfold);
    let mut nonhyp = 1i32;
    getparint("nonhyp", &mut nonhyp);
    let mut ngroup = 20i32;
    getparint("ngroup", &mut ngroup);
    let mut ndir = countparname("p");
    let mut dir = Vec::new();
    if ndir == 0 {
        ndir = -1;
    } else {
        // The C program leaves this VND allocation outstanding at exit.
        let _dir = VNDemalloc(ndir as usize * size_of::<*mut u8>(), "dir");
        for k in 0..ndir {
            let mut path = String::new();
            getnparstring(k + 1, "p", &mut path);
            dir.push(path);
        }
    }
    checkpars();

    lbtaper = (lbtaper as f64 / (1000.0 * dt as f64)) as i32;
    let ti = d.abs() > 0.0 || etamin.abs() > 0.0 || neta > 1;
    if ti {
        fprintf!(fpl, "sutifowler: operation in TI mode\n");
    }
    let mut deta = 0.0f32;
    if neta > 1 {
        deta = (etamax - etamin) / (neta - 1) as f32;
    }
    let mut dp = (1.0 / (vmin * (NP - 5) as f32) as f64) as f32;
    if ti {
        dp = (dp as f64 * (1.0 + 2.0 * etamin.abs() as f64).sqrt()) as f32;
    }
    if ichoose > 2 {
        nvstack = nv;
    }
    let mut ntfft = if ichoose == 1 || ichoose == 2 || ichoose == 3 {
        (ntout + ntpad as i64) as i32
    } else {
        1
    };
    let mut ntffts = if ichoose == 1 || ichoose == 3 {
        (ntout.wrapping_mul(2) as f64 / 0.6) as i32
    } else {
        1
    };

    ntfft = npfao(ntfft, 2 * ntfft);
    ntffts = npfao(ntffts, 2 * ntffts);
    let mut dw = (2.0 * PI / (ntfft as f32 * dtout) as f64) as f32;
    let nxfft = npfar(ncmps + nxpad);
    let dk = (2.0 * PI / (nxfft as f32 * dx) as f64) as f32;
    fprintf!(
        fpl,
        "sutifowler: ntfft=%d ntffts=%d nxfft=%d\n",
        ntfft,
        ntffts,
        nxfft
    );
    let czero = Complex { r: 0.0, i: 0.0 };
    let mut scale = 1.0f32;
    if ichoose < 5 {
        scale = (1.0 / nxfft as f64) as f32;
    }
    if ichoose == 1 || ichoose == 2 {
        scale = (scale as f64 * (1.0 / ntfft as f64)) as f32;
    }
    if ichoose == 1 || ichoose == 3 {
        scale = (scale as f64 * (1.0 / ntffts as f64)) as f32;
    }
    let nxfftny = nxfft / 2 + 1;
    let ntfftny = ntfft / 2 + 1;
    let mut nmax = nxfftny;
    if ntfft > nmax {
        nmax = ntfft;
    }
    if NP / 2 + 1 > nmax {
        nmax = NP / 2 + 1;
    }
    if nvstack > nmax {
        nmax = nvstack;
    }
    if nv * neta > nmax {
        nmax = nv * neta;
    }
    let mut ctemp = vec![czero; nmax as usize];
    let mut rindex = vec![0.0f32; nmax as usize];
    if ntffts > nmax {
        nmax = ntffts;
    }
    let mut crt = vec![czero; nmax as usize];
    // C aliases crt as float *rt and char *ccrt. real/cbytes below expose
    // the same initialized storage, with each borrow ending before the next.
    fprintf!(fpl, "sutifowler: nv=%d nvstack=%d\n", nv, nvstack);
    // Negative counts drive empty C loops. A one-element allocation removes
    // VNDemalloc's size wraparound without changing either negative-count case's
    // output (including the existing nvstack validation below).
    let mut v = vec![0.0f32; nv.max(1) as usize];
    let mut p2stack = vec![0.0f32; nvstack.max(1) as usize];
    // C stores one extra offset/mute before checking noff > mxfold. Reserve
    // that slot so the original error is reached without a Rust bounds panic.
    let mut mute = vec![0.0f32; mxfold as usize + 1];
    let mut off = vec![0.0f32; mxfold as usize + 1];
    fprintf!(fpl, "sutifowler: allocating and filling w array\n");
    let mut w = vec![0.0f32; ntfft as usize];
    for iw in 0..ntfft {
        if iw < ntfftny {
            w[iw as usize] = iw as f32 * dw;
        } else {
            w[iw as usize] = (iw - ntfft) as f32 * dw;
        }
        if iw == 0 {
            w[0] = (0.1 * dw as f64) as f32;
        }
    }

    /* build function for stacking velocity analysis */
    fprintf!(
        fpl,
        "sutifowler: building function for stacking velocity analysis\n"
    );
    let dv = (vmax - vmin) / (nv - 1).max(1) as f32;
    for iv in 0..nv {
        v[iv as usize] = vmin + iv as f32 * dv;
    }
    let mut dp2 = 0.0f32;
    if ichoose >= 3 {
        for iv in 0..nvstack {
            p2stack[iv as usize] = (1.0 / (v[iv as usize] * v[iv as usize]) as f64) as f32;
            fprintf!(
                fpl,
                "\t    stacking velocity %ld %f\n",
                iv as i64,
                v[iv as usize]
            );
        }
    } else {
        if nvstack < 6 {
            err!("sutifowler: nvstack must be 6 or more");
        }
        dp2 = (1.0 / (vminstack * vminstack * (nvstack - 5) as f32) as f64) as f32;
        for iv in 0..nvstack {
            p2stack[iv as usize] = iv as f32 * dp2;
            if iv > 0 {
                let factor = (1.0 / (p2stack[iv as usize] as f64).sqrt()) as f32;
                fprintf!(fpl, "\t    stacking velocity %ld %f\n", iv as i64, factor);
            } else {
                fprintf!(fpl, "\t    stacking velocity %ld infinity\n", iv as i64);
            }
        }
    }

    fprintf!(
        fpl,
        "sutifowler: Opening and zeroing large block matrix disk file\n"
    );
    fprintf!(
        fpl,
        "\t    This can take a while, but all is fruitless if the \n"
    );
    fprintf!(fpl, "\t    necessary disk space is not there...\n");
    let n = [(nxfft + 2) as i64, ntout * (nv * neta).max(nvstack) as i64];
    let fname = VNDtempname(&file);
    let mut vnd = VNDop(
        2,
        0,
        &n,
        1,
        size_of::<f32>() as i64,
        unsafe { CStr::from_ptr(fname.cast()) }.to_str().unwrap(),
        ndir,
        &dir.iter().map(String::as_str).collect::<Vec<_>>(),
        1,
    )
    .unwrap_or_else(|| null_file());
    unsafe { VNDfree(fname, "main: freeing fname 1") };
    fprintf!(
        fpl,
        "sutifowler: large file RAM mem buf = %ld bytes\n",
        vnd.num_bytes_mem_buf
    );
    fprintf!(
        fpl,
        "sutifowler: large file disk area = %ld bytes\n",
        vnd.num_bytes_per_block * vnd.num_blocks_per_panel * vnd.num_panels
    );

    if getcvstacks != 0 {
        fprintf!(fpl, "sutifowler: reading input cvstacks\n");
        for icmp in 0..ncmps {
            let key = [icmp as i64, 0];
            for iv in 0..nvstack {
                VNDrw(
                    b'w',
                    0,
                    &mut vnd,
                    1,
                    &key,
                    0,
                    fbytes(&mut tr.data),
                    iv as i64 * ntout,
                    1,
                    ntout,
                    1,
                    "writing cvstacks to disk",
                );
                if gettr(&mut tr) == 0 {
                    if icmp == ncmps - 1 && iv == nvstack - 1 {
                        /* all ok, read all the input data */
                    } else {
                        err!("sutifowler: error reading input cvstacks");
                    }
                }
            }
        }
        // Original goto xffts skips gather setup and leaves input headers intact.
    } else {
        fprintf!(
            fpl,
            "sutifowler: beginning constant velocity stacks of the input cmp gathers\n"
        );
        let fname = VNDtempname(&file);
        let mut vnda = V2Dop(
            2,
            1000000,
            size_of::<f32>() as i64,
            unsafe { CStr::from_ptr(fname.cast()) }.to_str().unwrap(),
            nt,
            mxfold as i64,
        )
        .unwrap_or_else(|| null_file());
        unsafe { VNDfree(fname, "main: freeing fname 2") };
        fprintf!(
            fpl,
            "sutifowler: cmp gather RAM mem buf = %ld bytes\n",
            vnda.num_bytes_mem_buf
        );
        let mut icmp = 0;
        let mut noff = 0;
        // The labelled break implements goto done_with_input, including the
        // skipped final cvstack and icmp increment when ncdps is reached.
        'done_with_input: {
            loop {
                if tr.cdp != oldcmp {
                    cvstack(
                        &mut vnda, &mut vnd, icmp, noff, &off, &mute, lmute, nvstack, &p2stack, dt,
                        dtout, &fpl,
                    );
                    icmp += 1;
                    if icmp == ncmps {
                        fprintf!(fpl, "sutifowler: more input cdps than ncdps parameter\n");
                        fprintf!(fpl, "\t    Will only process ncdps gathers.\n");
                        break 'done_with_input;
                    }
                    oldcmp = tr.cdp;
                    noff = 0;
                }
                if lbtaper > 0 || lstaper > 0 {
                    taper(lstaper, lbtaper, ncmps, icmp, nt as i32, &mut tr.data);
                }
                let factor = scale;
                for it in 0..nt {
                    tr.data[it as usize] *= factor;
                }
                V2Dw0(&mut vnda, noff as i64, fbytes(&mut tr.data), 1);
                off[noff as usize] = tr.offset as f32;
                if ichoose == 1 || ichoose == 2 {
                    mute[noff as usize] = fmax * off[noff as usize] * off[noff as usize] * dp2;
                } else {
                    mute[noff as usize] = 0.0;
                }
                if nonhyp != 0 {
                    let other = 2.0f32 * off[noff as usize] / vmin;
                    mute[noff as usize] = if mute[noff as usize] > other {
                        mute[noff as usize]
                    } else {
                        other
                    };
                }
                noff += 1;
                if noff > mxfold {
                    err!("tifowler: input cdp has more traces than mxfold");
                }
                if gettr(&mut tr) == 0 {
                    break;
                }
            }
            cvstack(
                &mut vnda, &mut vnd, icmp, noff, &off, &mute, lmute, nvstack, &p2stack, dt, dtout,
                &fpl,
            );
            icmp += 1;
        }
        ncmps = icmp;
        fprintf!(fpl, "sutifowler: read and stacked %d cmp gathers\n", ncmps);
        VNDcl(vnda, 1);
    }
    /* xffts: */
    VNDflush(&mut vnd);
    if ichoose < 5 {
        fprintf!(fpl, "sutifowler: doing forward x -> k spatial fft's\n");
        for it in 0..ntout * nvstack as i64 {
            V2Dr0(&mut vnd, it, cbytes(&mut crt), 21);
            for k in ncmps..nxfft + 2 {
                real(&mut crt)[k as usize] = 0.0;
            }
            // The Rust FFT API uses disjoint slices. Snapshot the C input
            // view before overwriting the aliased output buffer.
            let rt = real(&mut crt)[..nxfft as usize].to_vec();
            pfarc(1, &rt, &mut crt[..nxfftny as usize]);
            V2Dw0(&mut vnd, it, cbytes(&mut crt), 22);
        }
        VNDr2c(&mut vnd);
    }

    if ichoose <= 3 {
        fprintf!(fpl, "sutifowler: looping over k\n");
        let mut vndvnmo = None;
        let mut vndvphase = None;
        let mut vndb = None;
        if ti && (ichoose == 1 || ichoose == 2) {
            let table = ptabledmo(nv, &v, etamin, deta, neta, d, vsvp, NP, dp, dp2, &file);
            fprintf!(
                fpl,
                "sutifowler: dmo index(p) RAM mem buf = %ld bytes\n",
                table.num_bytes_mem_buf
            );
            vndvnmo = Some(table);
        }
        if ti && (ichoose == 1 || ichoose == 3) {
            let table = ptablemig(nv, &v, etamin, deta, neta, d, vsvp, NP, &file);
            fprintf!(
                fpl,
                "sutifowler: migration scaler(p) RAM mem buf = %ld bytes\n",
                table.num_bytes_mem_buf
            );
            vndvphase = Some(table);
        }
        if ichoose == 1 || ichoose == 2 {
            let iv = (nv * neta).max(nvstack) as i64;
            let fname = VNDtempname(&file);
            let matrix = V2Dop(
                2,
                750000,
                size_of::<Complex>() as i64,
                unsafe { CStr::from_ptr(fname.cast()) }.to_str().unwrap(),
                ntfft as i64,
                iv,
            )
            .unwrap_or_else(|| null_file());
            fprintf!(
                fpl,
                "sutifowler: (w,v) RAM mem buf = %ld bytes\n",
                matrix.num_bytes_mem_buf
            );
            unsafe { VNDfree(fname, "main: freeing fname 3") };
            vndb = Some(matrix);
        }

        /* loop over spatial wavenumbers */
        for k in 0..nxfftny {
            if k == 20 * (k / 20) {
                fprintf!(fpl, "sutifowler: k index = %d out of %d\n", k, nxfftny);
            }
            let ak = k as f32 * dk;
            let key = [k as i64, 0];
            if ichoose == 1 || ichoose == 2 {
                /* do Fowler DMO */
                let vndb = vndb.as_mut().unwrap();
                for iv in 0..nvstack {
                    VNDrw(
                        b'r',
                        0,
                        &mut vnd,
                        1,
                        &key,
                        0,
                        cbytes(&mut crt),
                        iv as i64 * ntout,
                        1,
                        ntout,
                        31,
                        "Fowler DMO t -> w fft read",
                    );
                    for it in ntout..ntfft as i64 {
                        crt[it as usize] = czero;
                    }
                    pfacc(-1, &mut crt[..ntfft as usize]);
                    V2Dw0(vndb, iv as i64, cbytes(&mut crt), 32);
                }
                for iw in 0..ntfft {
                    let p = (0.5 * ak as f64 / (w[iw as usize] as f64).abs()) as f32;
                    if ti {
                        let ip = (p / dp) as i64;
                        if ip < NP as i64 {
                            V2Dr0(vndvnmo.as_mut().unwrap(), ip, fbytes(&mut rindex), 40);
                        } else {
                            for iv in 0..nv * neta {
                                rindex[iv as usize] = -1.0;
                            }
                        }
                    } else {
                        let p2 = p * p;
                        for iv in 0..nv {
                            let v2 = v[iv as usize] * v[iv as usize];
                            rindex[iv as usize] = (1.0f32 - v2 * p2) / (v2 * dp2);
                        }
                    }
                    V2Dr1(vndb, iw as i64, cbytes(&mut crt), 41);
                    for iv in 0..nvstack {
                        ctemp[iv as usize] = crt[iv as usize];
                    }
                    ints8c(
                        1.0,
                        0.0,
                        &ctemp[..nvstack as usize],
                        czero,
                        czero,
                        &rindex[..(nv * neta).max(0) as usize],
                        &mut crt[..(nv * neta).max(0) as usize],
                    );
                    V2Dw1(vndb, iw as i64, cbytes(&mut crt), 42);
                }
                for iv in 0..nv * neta {
                    V2Dr0(vndb, iv as i64, cbytes(&mut crt), 51);
                    pfacc(1, &mut crt[..ntfft as usize]);
                    VNDrw(
                        b'w',
                        0,
                        &mut vnd,
                        1,
                        &key,
                        0,
                        cbytes(&mut crt),
                        iv as i64 * ntout,
                        1,
                        ntout,
                        52,
                        "Fowler DMO w -> t fft write",
                    );
                }
            }
            if ichoose == 3 && neta > 1 {
                /* fix up disk order if only doing TI migrations */
                for iv in 0..nv {
                    VNDrw(
                        b'r',
                        0,
                        &mut vnd,
                        1,
                        &key,
                        0,
                        cbytes(&mut crt),
                        iv as i64 * ntout,
                        1,
                        ntout,
                        57,
                        "option 3 fixup for multiple eta read",
                    );
                    for ieta in 1..neta {
                        VNDrw(
                            b'w',
                            0,
                            &mut vnd,
                            1,
                            &key,
                            0,
                            cbytes(&mut crt),
                            iv as i64 * ntout + ieta as i64 * nv as i64 * ntout,
                            1,
                            ntout,
                            58,
                            "option 3 fixup for multiple eta write",
                        );
                    }
                }
            }
            if ichoose == 1 || ichoose == 3 {
                /* do Stolt migration */
                for iv in 0..nv * neta {
                    if ti {
                        V2Dr0(vndvphase.as_mut().unwrap(), iv as i64, cbytes(&mut crt), 50);
                        let rt = real(&mut crt);
                        let dpm = rt[0];
                        dw = (2.0 * PI / (ntfft as f32 * dtout) as f64) as f32;
                        let iwmin = (0.5 * ak as f64 / ((NP - 3) as f32 * dpm * dw) as f64) as i32;
                        for iw in iwmin + 1..ntfftny {
                            let p = (0.5 * ak as f64 / (w[iw as usize] as f64).abs()) as f32;
                            let rp = (1.0 + (p / dpm) as f64) as f32;
                            let ip = rp as i64;
                            let wgt = rp - ip as f32;
                            let factor = ((wgt * rt[(ip + 1) as usize]) as f64
                                + (1.0 - wgt as f64) * rt[ip as usize] as f64)
                                as f32;
                            rindex[iw as usize] = w[iw as usize] * factor;
                            rindex[(ntfft - iw) as usize] = w[(ntfft - iw) as usize] * factor;
                        }
                        let fw = (-2.0 * PI / dtout as f64) as f32;
                        rindex[0] = fw;
                        // The repaired C stops before ntfft: larger iw writes
                        // rindex[ntfft-iw] before its allocation and corrupts
                        // VND's allocation guard (sutifowler.c:565-567).
                        for iw in 1..(iwmin + 1).min(ntfft) {
                            rindex[iw as usize] = fw;
                            rindex[(ntfft - iw) as usize] = fw;
                        }
                    } else {
                        scale = (0.5 * v[iv as usize] as f64 * ak as f64) as f32;
                        for iw in 0..ntfft {
                            if (w[iw as usize] as f64).abs() > scale as f64 {
                                let mut factor = scale / w[iw as usize];
                                factor = ((1.0f32 + factor * factor) as f64).sqrt() as f32;
                                rindex[iw as usize] = w[iw as usize] * factor;
                            } else {
                                rindex[iw as usize] = (-2.0 * PI / dtout as f64) as f32;
                            }
                        }
                    }
                    VNDrw(
                        b'r',
                        0,
                        &mut vnd,
                        1,
                        &key,
                        0,
                        cbytes(&mut crt),
                        iv as i64 * ntout,
                        1,
                        ntout,
                        61,
                        "Stolt t -> w fft read",
                    );
                    for it in (1..ntout).step_by(2) {
                        crt[it as usize].r = -crt[it as usize].r;
                        crt[it as usize].i = -crt[it as usize].i;
                    }
                    for it in ntout..ntffts as i64 {
                        crt[it as usize] = czero;
                    }
                    pfacc(1, &mut crt[..ntffts as usize]);
                    dw = (2.0 * PI / (ntffts as f32 * dtout) as f64) as f32;
                    let fw = (-PI / dtout as f64) as f32;
                    ints8c(
                        dw,
                        fw,
                        &crt[..ntffts as usize],
                        czero,
                        czero,
                        &rindex[..ntfft as usize],
                        &mut ctemp[..ntfft as usize],
                    );
                    /* obliquity factor code */
                    for iw in 0..ntfft {
                        let factor = ((w[iw as usize] / rindex[iw as usize]) as f64).abs() as f32;
                        crt[iw as usize].r = factor * ctemp[iw as usize].r;
                        crt[iw as usize].i = factor * ctemp[iw as usize].i;
                    }
                    pfacc(-1, &mut crt[..ntfft as usize]);
                    VNDrw(
                        b'w',
                        0,
                        &mut vnd,
                        1,
                        &key,
                        0,
                        cbytes(&mut crt),
                        iv as i64 * ntout,
                        1,
                        ntout,
                        62,
                        "Stolt w->t fft write",
                    );
                }
            }
        }
        fprintf!(fpl, "sutifowler: completed loop over wavenumbers\n");
        if let Some(matrix) = vndb {
            VNDcl(matrix, 1);
        }
        if let Some(table) = vndvnmo {
            VNDcl(table, 1);
        }
        if let Some(table) = vndvphase {
            VNDcl(table, 1);
        }
    }
    if ichoose < 5 {
        fprintf!(fpl, "sutifowler: doing inverse spatial fft's k->x\n");
        for it in 0..ntout * nv as i64 * neta as i64 {
            V2Dr0(&mut vnd, it, cbytes(&mut crt), 71);
            let input = crt[..nxfftny as usize].to_vec();
            pfacr(-1, &input, &mut real(&mut crt)[..nxfft as usize]);
            V2Dw0(&mut vnd, it, cbytes(&mut crt), 72);
        }
        VNDc2r(&mut vnd);
    }

    fprintf!(fpl, "sutifowler: outputting results\n");
    let mut it = 0i64;
    for icmp in 0..ncmps {
        let key = [icmp as i64, 0];
        for ieta in 0..neta {
            for iv in 0..nv {
                VNDrw(
                    b'r',
                    0,
                    &mut vnd,
                    1,
                    &key,
                    0,
                    fbytes(&mut tr.data),
                    iv as i64 * ntout + ieta as i64 * nv as i64 * ntout,
                    1,
                    ntout,
                    82,
                    "outputting all velocities for each cmp",
                );
                tr.ns = ntout as u16;
                tr.dt = (1000000.0f32 * dtout) as i32 as u16;
                tr.cdp = icmp;
                tr.tracf = iv;
                tr.offset = v[iv as usize] as i32;
                tr.cdpt = iv;
                tr.sx = (icmp as f32 * dx) as i32;
                tr.gx = (icmp as f32 * dx) as i32;
                it += 1;
                tr.tracl = it as i32;
                tr.tracr = it as i32;
                // ARM64 integer division by zero returns zero in the C build.
                tr.fldr = icmp.checked_div(ngroup).unwrap_or(0);
                tr.ep = 10 + tr.fldr * ngroup;
                tr.igc = ieta as i16;
                tr.igi = (100.0f32 * (etamin + ieta as f32 * deta)) as i32 as i16;
                tr.d1 = dtout;
                tr.f1 = 0.0;
                tr.d2 = 1.0;
                tr.f2 = 0.0;
                puttr(&mut tr);
            }
        }
    }
    /* close files and return */
    VNDcl(vnd, 1);
    drop(crt);
    drop(ctemp);
    drop(v);
    drop(p2stack);
    drop(mute);
    drop(off);
    drop(rindex);
    drop(w);
    if VNDtotalmem() != 0 {
        fprintf!(stderr(), "total VND memory at end = %ld\n", VNDtotalmem());
    }
    exit(CWP_Exit());
}

/// Compute constant velocity stacks for this cmp and load out to disk.
/// John Anderson (visitor to CSM from Mobil), Spring 1993.
#[allow(clippy::too_many_arguments)]
fn cvstack(
    vnda: &mut VND,
    vnd: &mut VND,
    icmp: i32,
    noff: i32,
    off: &[f32],
    mute: &[f32],
    lmute: i32,
    nv: i32,
    p2: &[f32],
    dt: f32,
    dtout: f32,
    fpl: &CFile,
) {
    let nt = vnda.n[0];
    // ARM64 integer division by zero returns zero in the C build.
    let ntout = vnd.n[1].checked_div(nv as i64).unwrap_or(0);
    let mut rt = vec![0.0f32; ntout as usize];
    let mut buf = vec![0.0f32; nt.max(ntout) as usize];
    let mut ttn = vec![0.0f32; ntout as usize];
    let mut qtn = vec![0.0f32; ntout as usize];
    let mut fold = vec![0.0f32; ntout as usize];
    let mut big = vec![0.0f32; noff as usize];
    /* check for completely dead traces so can skip them */
    for ioff in 0..noff {
        big[ioff as usize] = 0.0;
        V2Dr0(vnda, ioff as i64, fbytes(&mut buf), 3);
        for it in 0..nt {
            let value = (buf[it as usize] as f64).abs();
            big[ioff as usize] = if big[ioff as usize] as f64 > value {
                big[ioff as usize]
            } else {
                value as f32
            };
        }
    }
    for iv in 0..nv {
        for it in 0..ntout {
            rt[it as usize] = 0.0;
        }
        for it in 0..ntout {
            fold[it as usize] = 0.0;
        }
        for ioff in 0..noff {
            if big[ioff as usize] > 0.0 {
                V2Dr0(vnda, ioff as i64, fbytes(&mut buf), 4);
                let factor = off[ioff as usize] * off[ioff as usize] * p2[iv as usize];
                let itmute = (mute[ioff as usize] / dtout) as i32;
                for it in itmute as i64..ntout {
                    let t = it as f32 * dtout;
                    ttn[it as usize] = ((t * t + factor) as f64).sqrt() as f32;
                }
                /* do nmo via 8-point sinc interpolation */
                // C passes ntout-itmute to ints8r. A nonpositive count does
                // no interpolation; the subsequent taper/sum loops also skip.
                if (itmute as i64) < ntout {
                    ints8r(
                        dt,
                        0.0,
                        &buf[..nt as usize],
                        0.0,
                        0.0,
                        &ttn[itmute as usize..ntout as usize],
                        &mut qtn[itmute as usize..ntout as usize],
                    );
                }
                /* apply linear ramp to taper mute */
                let mut it = itmute;
                while it < itmute + lmute && (it as i64) < ntout {
                    qtn[it as usize] *= (it - itmute + 1) as f32 / lmute as f32;
                    it += 1;
                }
                /* sum NMO corrected trace to stacked trace */
                for it in itmute as i64..ntout {
                    rt[it as usize] += qtn[it as usize];
                }
                /* count fold information */
                for it in itmute as i64..ntout {
                    fold[it as usize] += 1.0;
                }
            }
        }
        for it in 0..ntout {
            if fold[it as usize] == 0.0 {
                fold[it as usize] = 1.0;
            }
        }
        for it in 0..ntout {
            rt[it as usize] /= fold[it as usize];
        }
        let key = [icmp as i64, 0];
        VNDrw(
            b'w',
            0,
            vnd,
            1,
            &key,
            0,
            fbytes(&mut rt),
            iv as i64 * ntout,
            1,
            ntout,
            1,
            "cv stacking",
        );
    }
    drop(rt);
    drop(buf);
    drop(ttn);
    drop(qtn);
    drop(fold);
    drop(big);
    if icmp == 20 * (icmp / 20) {
        fprintf!(
            fpl,
            "sutifowler: completed stacking velocity analysis for cmp %d\n",
            icmp
        );
    }
}

/// Build the TI DMO table, dimension 0 along velocity and 1 along slowness.
/// Values are 1/(vnmo(p)*vnmo(p)*dp2), the stacking velocity analysis index.
/// John Anderson (visitor to CSM from Mobil), Spring 1993.
#[allow(clippy::too_many_arguments)]
fn ptabledmo(
    nv: i32,
    v: &[f32],
    etamin: f32,
    deta: f32,
    neta: i32,
    d: f32,
    vsvp: f32,
    np: i32,
    dp: f32,
    dp2: f32,
    file: &str,
) -> VND {
    let newnv = (nv * neta) as i64;
    let newnp = np as i64;
    let fname = VNDtempname(file);
    let dscale = (1.0 / ((1.0f32 + 2.0f32 * d) as f64).sqrt()) as f32;
    let mut vndvel = V2Dop(
        2,
        125000,
        size_of::<f32>() as i64,
        unsafe { CStr::from_ptr(fname.cast()) }.to_str().unwrap(),
        newnv,
        newnp,
    )
    .unwrap_or_else(|| null_file());
    unsafe { VNDfree(fname, "ptabledmo: freeing fname") };
    let nangle = 10001;
    let dangle = (89.0 / (nangle - 1) as f64) as f32;
    let mut pin = vec![0.0f32; nangle as usize];
    let mut vin = vec![0.0f32; nangle as usize];
    let mut pout = vec![0.0f32; np as usize];
    let mut vout = vec![0.0f32; np as usize];
    let mut vel = [0.0f32; 5];
    for ip in 0..np {
        pout[ip as usize] = ip as f32 * dp;
    }
    for ieta in 0..neta {
        let e = (etamin + ieta as f32 * deta) * (1.0f32 + 2.0f32 * d) + d;
        for iv in 0..nv {
            let a = v[iv as usize] * dscale; /* vertical p-wave phase velocity */
            let b = vsvp * a; /* vertical s-wave phase velocity */
            for iangle in 0..nangle {
                let angle = iangle as f32 * dangle;
                let theta = (angle as f64 * PI / 180.0) as f32;
                vget(a, b, e, d, theta, &mut vel);
                pin[iangle as usize] = vel[4];
                vin[iangle as usize] = vel[3];
            }
            intlin(
                &pin,
                &vin,
                vin[0],
                vin[(nangle - 1) as usize],
                &pout,
                &mut vout,
            );
            for ip in 0..np {
                vout[ip as usize] =
                    (1.0 / (vout[ip as usize] * vout[ip as usize] * dp2) as f64) as f32;
            }
            V2Dw1(&mut vndvel, (iv + ieta * nv) as i64, fbytes(&mut vout), 101);
        }
    }
    drop(pin);
    drop(vin);
    drop(pout);
    drop(vout);
    vndvel
}

/// Build the TI Stolt migration table. Its first dimension contains dp
/// followed by vphase/(a*cos(theta)), sampled in tan(theta)/a; dimension 1
/// indexes output velocity and eta. John Anderson, Spring 1993.
#[allow(clippy::too_many_arguments)]
fn ptablemig(
    nv: i32,
    v: &[f32],
    etamin: f32,
    deta: f32,
    neta: i32,
    d: f32,
    vsvp: f32,
    np: i32,
    file: &str,
) -> VND {
    let newnv = (nv * neta) as i64;
    let newnp = np as i64;
    let fname = VNDtempname(file);
    let mut vndvel = V2Dop(
        2,
        125000,
        size_of::<f32>() as i64,
        unsafe { CStr::from_ptr(fname.cast()) }.to_str().unwrap(),
        newnp,
        newnv,
    )
    .unwrap_or_else(|| null_file());
    unsafe { VNDfree(fname, "pmigtable: freeing fname") };
    let dscale = (1.0 / ((1.0f32 + 2.0f32 * d) as f64).sqrt()) as f32;
    let nangle = 10001;
    let dangle = (89.0 / (nangle - 1) as f64) as f32;
    let mut pin = vec![0.0f32; nangle as usize];
    let mut vin = vec![0.0f32; nangle as usize];
    let mut pout = vec![0.0f32; np as usize];
    let mut vout = vec![0.0f32; np as usize];
    let mut vel = [0.0f32; 5];
    for ieta in 0..neta {
        let e = (etamin + ieta as f32 * deta) * (1.0f32 + 2.0f32 * d) + d;
        for iv in 0..nv {
            let a = v[iv as usize] * dscale; /* vertical p-wave phase velocity */
            let b = vsvp * a; /* vertical s-wave phase velocity */
            for iangle in 0..nangle {
                let angle = iangle as f32 * dangle;
                let theta = (angle as f64 * PI / 180.0) as f32;
                vget(a, b, e, d, theta, &mut vel);
                pin[iangle as usize] = ((theta as f64).tan() / a as f64) as f32;
                vin[iangle as usize] = (vel[0] as f64 / (a as f64 * (theta as f64).cos())) as f32;
            }
            let pmax = pin[(nangle - 1) as usize];
            let dp = pmax / (np - 2) as f32;
            vout[0] = dp;
            for ip in 0..np - 1 {
                pout[ip as usize] = ip as f32 * dp;
            }
            intlin(
                &pin,
                &vin,
                vin[0],
                vin[(nangle - 1) as usize],
                &pout[..(np - 1) as usize],
                &mut vout[1..],
            );
            V2Dw0(&mut vndvel, (iv + ieta * nv) as i64, fbytes(&mut vout), 101);
        }
    }
    drop(pin);
    drop(vin);
    drop(pout);
    drop(vout);
    vndvel
}

// An unchecked null FILE*/VND* in the original terminates with SIGSEGV.
// Raise the same signal without introducing an invalid Rust dereference.
fn null_file() -> ! {
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(139);
}

// These views translate the C pointer casts; they perform no arithmetic or I/O.
// Complex is repr(C), exactly two adjacent f32s without padding. f32 accepts
// every bit pattern. The returned views have exclusive, bounded lifetimes.
fn real(values: &mut [Complex]) -> &mut [f32] {
    unsafe { std::slice::from_raw_parts_mut(values.as_mut_ptr().cast(), values.len() * 2) }
}

fn cbytes(values: &mut [Complex]) -> &mut [u8] {
    unsafe { std::slice::from_raw_parts_mut(values.as_mut_ptr().cast(), size_of_val(values)) }
}

fn fbytes(values: &mut [f32]) -> &mut [u8] {
    unsafe { std::slice::from_raw_parts_mut(values.as_mut_ptr().cast(), size_of_val(values)) }
}

/// Phase and NMO velocities for a TI medium. Program-local C helper;
/// John Anderson (visitor to CSM from Mobil), Spring 1993.
/// vel: phase velocity, its first and second angular derivatives, NMO
/// velocity, and ray parameter sin(theta)/vphase.
fn vget(a: f32, b: f32, e: f32, d: f32, theta: f32, vel: &mut [f32; 5]) {
    let eps = 1.0e-7f32;
    let sint = (theta as f64).sin() as f32;
    let mut cost = (theta as f64).cos() as f32;
    let sin2t = sint * sint;
    let cos2t = cost * cost;
    let psi = (1.0 - ((b * b) / (a * a)) as f64) as f32;
    let gamma = (0.25 * psi as f64 * psi as f64
        + ((2.0f32 * d - e) * psi * sin2t * cos2t) as f64
        + ((psi + e) * e * sin2t * sin2t) as f64) as f32;
    let dgamma = (2.0f32 * d - e) * psi * 2.0f32 * sint * cost * cos2t
        + (4.0f32 * (psi + e) * e - 2.0f32 * (2.0f32 * d - e) * psi) * sint * sin2t * cost;
    let d2gamma = (((2.0f32 * d - e) * psi * 2.0f32) as f64
        * ((cos2t * cos2t) as f64 - 3.0 * cos2t as f64 * sin2t as f64)
        + (4.0f32 * (psi + e) * e - 2.0f32 * (2.0f32 * d - e) * psi) as f64
            * (3.0 * cos2t as f64 * sin2t as f64 - (sin2t * sin2t) as f64))
        as f32;
    let sqgamma = (gamma as f64).sqrt() as f32;
    let v =
        (a as f64 * (1.0 + (e * sin2t) as f64 - 0.5 * psi as f64 + sqgamma as f64).sqrt()) as f32;
    let dv = (0.5
        * a as f64
        * a as f64
        * ((2.0f32 * e * sint * cost) as f64 + 0.5 * dgamma as f64 / sqgamma as f64)
        / v as f64) as f32;
    let d2v = (0.5
        * a as f64
        * a as f64
        * ((2.0f32 * e * (cos2t - sin2t)) as f64 + 0.5 * d2gamma as f64 / sqgamma as f64
            - 0.25 * dgamma as f64 / (sqgamma * gamma) as f64)
        / v as f64
        - 0.5
            * a as f64
            * a as f64
            * dv as f64
            * ((2.0f32 * e * sint * cost) as f64 + 0.5 * dgamma as f64 / sqgamma as f64)
            / (v * v) as f64) as f32;
    if cost < eps {
        cost = eps;
    }
    let tant = sint / cost;
    let vnmo = (v as f64 * ((1.0f32 + d2v / v) as f64).sqrt()
        / (cost * (1.0f32 - tant * dv / v)) as f64) as f32;
    let p = sint / v;
    vel[0] = v;
    vel[1] = dv;
    vel[2] = d2v;
    vel[3] = vnmo;
    vel[4] = p;
}

/// Taper traces near the sides and bottom of the trace window.
/// Program-local helper based on code by David Hale.
fn taper(lxtaper: i32, lbtaper: i32, nx: i32, ix: i32, nt: i32, trace: &mut [f32]) {
    let xtaper: f32;
    /* if near left side */
    if ix < lxtaper {
        xtaper =
            (0.54 + 0.46 * (PI * lxtaper.wrapping_sub(ix) as f64 / lxtaper as f64).cos()) as f32;
    /* else if near right side */
    } else if ix >= nx.wrapping_sub(lxtaper) {
        xtaper = (0.54
            + 0.46
                * (PI * lxtaper.wrapping_add(ix).wrapping_add(1).wrapping_sub(nx) as f64
                    / lxtaper as f64)
                    .cos()) as f32;
    /* else x tapering is unnecessary */
    } else {
        xtaper = 1.0;
    }
    /* if x tapering is necessary, apply it */
    if xtaper != 1.0 {
        for it in 0..nt {
            trace[it as usize] *= xtaper;
        }
    }
    /* if requested, apply t tapering */
    for it in 0.max(nt.wrapping_sub(lbtaper))..nt {
        trace[it as usize] = (trace[it as usize] as f64
            * (0.54
                + 0.46
                    * (PI * lbtaper.wrapping_add(it).wrapping_add(1).wrapping_sub(nt) as f64
                        / lbtaper as f64)
                        .cos())) as f32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotropic_velocities_obey_angle_dependence() {
        // In an isotropic medium phase velocity is independent of angle,
        // both angular derivatives vanish, and vnmo = a/cos(theta).
        for i in 0..90 {
            let theta = (i as f64 * PI / 180.0) as f32;
            let mut vel = [0.0; 5];
            vget(2000.0, 1000.0, 0.0, 0.0, theta, &mut vel);
            assert_eq!(vel[0], 2000.0);
            assert_eq!(vel[1], 0.0);
            assert_eq!(vel[2], 0.0);
            let cost = (theta as f64).cos() as f32;
            let expected = (2000.0 / cost as f64) as f32;
            assert_eq!(vel[3], expected);
            assert_eq!(vel[4], (theta as f64).sin() as f32 / 2000.0);
        }
    }

    #[test]
    fn vertical_ti_velocity_matches_thomsen_delta() {
        for d in [-0.1f32, 0.0, 0.05, 0.2] {
            let mut vel = [0.0; 5];
            vget(2000.0, 1000.0, 0.3, d, 0.0, &mut vel);
            assert_eq!(vel[0], 2000.0);
            assert_eq!(vel[1], 0.0);
            assert_eq!(vel[4], 0.0);
            let expected = 2000.0 * (1.0 + 2.0 * d as f64).sqrt();
            assert!((vel[3] as f64 - expected).abs() < 0.001);
        }
    }

    #[test]
    fn ti_velocities_match_original_static_helper() {
        // Reference words were produced by compiling the original unmodified
        // sutifowler.c into a test-only driver and calling its static vget.
        // Cover general TI, negative anisotropy, near-horizontal incidence,
        // and the cost<eps branch; these are not identities of this port.
        let cases = [
            (
                [2000.0, 1000.0, 0.3, 0.05, 0.4],
                [1157471064, 1129008489, 1149589705, 1160877126, 961102720],
            ),
            (
                [3200.0, 1440.0, 0.2, 0.1, 1.2],
                [1164362997, 1139772076, 3299414652, 1178483214, 964974641],
            ),
            (
                [1500.0, 750.0, -0.1, -0.05, 0.8],
                [1152667890, 3273141723, 1125474542, 1156922990, 973233304],
            ),
            (
                [6000.0, 2400.0, 0.5, 0.2, 1.55],
                [1174704243, 1122025599, 3316251350, 1227115824, 955720161],
            ),
            (
                [2000.0, 1000.0, 0.3, 0.05, 1.5707964],
                [1159601448, 3095295026, 3299703918, 1347176558, 969883190],
            ),
        ];
        for (x, expected) in cases {
            let mut vel = [0.0; 5];
            vget(x[0], x[1], x[2], x[3], x[4], &mut vel);
            for j in 0..5 {
                let reference = f32::from_bits(expected[j]);
                assert!(
                    (vel[j] - reference).abs() <= 1.0e-5 * reference.abs(),
                    "input={x:?} component={j} actual={} reference={reference}",
                    vel[j]
                );
            }
        }
    }

    #[test]
    fn taper_endpoints_and_untapered_samples() {
        let mut left = [1.0f32; 6];
        let mut right = left;
        taper(2, 0, 7, 0, 6, &mut left);
        taper(2, 0, 7, 6, 6, &mut right);
        assert_eq!(left, [0.08f32; 6]);
        assert_eq!(left, right);
        let mut middle = [1.0f32; 6];
        taper(2, 2, 7, 3, 6, &mut middle);
        assert_eq!(middle, [1.0, 1.0, 1.0, 1.0, 0.54, 0.08]);
        let mut untouched = [1.0f32; 6];
        taper(0, 0, 7, 3, 6, &mut untouched);
        assert_eq!(untouched, [1.0f32; 6]);
        taper(-1, -1, 7, 3, 6, &mut untouched);
        assert_eq!(untouched, [1.0f32; 6]);
    }
}
