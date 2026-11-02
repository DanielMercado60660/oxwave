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

//! SUFDMOD1 - Finite difference modelling (1-D first order) for the acoustic wave equation.
//! Port of su/main/synthetics_waveforms_testpatterns/sufdmod1.c (Revision 1.4).
//! The C source contains no separate author credits.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUFDMOD1 - Finite difference modelling (1-D 1rst order) for the\t",
    " acoustic wave equation\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t",
    " sufdmod1 <vfile >sfile nz= tmax= sz= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters :\t\t\t\t\t\t\t",
    " <vfile or vfile=\tbinary file containing velocities[nz]\t\t",
    " >sfile or sfile=\tSU file containing seimogram[nt]\t\t",
    " nz=\t\t number of z samples\t\t\t\t   \t",
    " tmax=\t\tmaximum propagation time\t\t\t\t",
    " sz=\t\t z coordinate of source\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters :\t\t\t\t\t\t\t",
    " dz=1\t   z sampling interval\t\t\t\t\t\t",
    " fz=0.0\t first depth sample\t\t\t\t\t",
    " rz=1\t   coordinate of receiver\t\t\t\t\t",
    " sz=1\t   coordinate of source\t\t\t\t\t\t",
    " dfile=\t binary input file containing density[nz]\t\t",
    " wfile=\t output file for wave field (snapshots in a SU trace panel)",
    " abs=0,1\tabsorbing conditions on top and bottom\t\t\t",
    " styp=0\t source type (0: gauss, 1: ricker 1, 2: ricker 2)\t",
    " freq=15.0\tapproximate source center frequency (Hz)\t\t",
    " nt=1+tmax/dt   number od time samples (dt determined for numerical\t",
    " stability)\t\t\t\t\t\t\t\t",
    " zt=1\t   trace undersampling factor for trace and snapshots\t \t",
    " zd=1\t   depth undersampling factor for snapshots\t\t   \t",
    " press=1\tto record the pressure field; 0 records the particle\t",
    "\t\tvelocity\t\t\t\t\t\t",
    " verbose=0\t=1 for diagnostic messages\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes :\t\t\t\t\t\t\t\t",
    "  This program uses a first order explicit velocity/pressure  finite\t",
    "  difference equation.\t\t\t\t\t\t\t",
    "  The source function is applied on the pressure component.\t\t",
    "  If no density file is given, constant density is assumed\t \t",
    "  Wavefield  can be easily viewed with suximage, user must provide f2=0",
    "  to the ximage program in order to  get correct time labelling\t",
    "  Seismic trace is shifted in order to get a zero phase source\t\t",
    "  Source begins and stop when it's amplitude is 10^-4 its maximum\t",
    "  Time and depth undersampling only modify the output trace and snapshots.",
    "  These parameters are useful for keeping snapshot file small and\t",
    "  the number of samples under SU_NFLTS.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut snapsh = Trace::new();
    let mut sismo = Trace::new();
    let mut seisfp = stdout();
    let mut wavefp = None;
    let mut velocityfp = stdin();
    let mut dfile = String::new();
    let mut wfile = String::new();
    let mut velfile = String::new();
    let mut sfile = String::new();

    /* hook up getpar to handle the parameters */
    initargs(&argv);
    requestdoc(0, SDOC);

    let mut verbose = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }

    /* get required parameters */
    let mut nz = 0i32;
    if getparint("nz", &mut nz) == 0 {
        err!("must specify nz ! ");
    }
    if verbose != 0 {
        warn!("nz= %d", nz);
    }
    let mut tmax = 0.0f32;
    if getparfloat("tmax", &mut tmax) == 0 {
        err!("must specify tmax ! ");
    }
    if verbose != 0 {
        warn!("tmax= %f", tmax);
    }
    let mut sz = 0.0f32;
    if getparfloat("sz", &mut sz) == 0 {
        err!("must specify sz ! ");
    }
    if verbose != 0 {
        warn!("sz= %f", sz);
    }

    /* get optional parameters */
    let mut nt = 0i32;
    if getparint("nt", &mut nt) == 0 {
        nt = 0;
    }
    if verbose != 0 {
        warn!("nt= %d", nt);
    }
    let mut styp = 0i32;
    if getparint("styp", &mut styp) == 0 {
        styp = 0;
    }
    if verbose != 0 {
        warn!("styp= %d ", styp);
    }
    let mut dz = 1.0f32;
    if getparfloat("dz", &mut dz) == 0 {
        dz = 1.0;
    }
    let mut fz = 0.0f32;
    if getparfloat("fz", &mut fz) == 0 {
        fz = 0.0;
    }

    /* source coordinates to samples */
    let isz = nint(((sz - fz) / dz) as f64);
    if verbose != 0 {
        warn!("source on knot number %d ", isz);
    }
    let mut rz = 0.0f32;
    if getparfloat("rz", &mut rz) == 0 {
        rz = 0.0;
    }
    let irz = nint(((rz - fz) / dz) as f64);
    if verbose != 0 {
        warn!("receiver depth : %f on knot # %d ", rz, irz);
    }
    let mut freq = 0.0f32;
    if getparfloat("freq", &mut freq) == 0 {
        freq = 15.0;
    }
    if verbose != 0 {
        warn!("frequency : %f  Hz", freq);
    }

    getparstring("velfile", &mut velfile);
    if verbose != 0 {
        if !velfile.is_empty() {
            warn!("Velocity file : %s ", velfile);
        } else {
            warn!("Velocity file supplied via stdin");
        }
    }
    getparstring("sfile", &mut sfile);
    if verbose != 0 {
        if !sfile.is_empty() {
            warn!("Output trace file : %s ", sfile);
        } else {
            warn!("Output trace via stdout");
        }
    }
    getparstring("dfile", &mut dfile);
    if verbose != 0 {
        if !dfile.is_empty() {
            warn!("Density file : %s ", dfile);
        } else {
            warn!("No density file supplied ");
        }
    }
    getparstring("wfile", &mut wfile);
    if verbose != 0 {
        if !wfile.is_empty() {
            warn!("Wave file : %s ", wfile);
        } else {
            warn!("No wave file requested ");
        }
    }
    // zd is still its initial value, 1, at this point in the C program.
    let mut zd = 1i32;
    if nint((nz as f32 / zd as f32) as f64).wrapping_add(1) > SU_NFLTS as i32 {
        warn!("Too many depth points : impossible to output wave field. Increase zd ?");
        wfile.clear();
    }

    /* get absorbing conditions */
    let mut abs = vec![0i32; 2];
    if getparint("abs", &mut abs) == 0 {
        abs[0] = 0;
        abs[1] = 1;
    }
    if verbose != 0 {
        if abs[0] == 1 {
            warn!("absorbing condition on top ");
        }
        if abs[1] == 1 {
            warn!("absorbing condition on bottom ");
        }
    }
    /* get decimation coefficients */
    let mut td = 1i32;
    if getparint("td", &mut td) == 0 {
        td = 1;
    }
    if verbose != 0 {
        warn!("time decimation ccoefficent: %d ", td);
    }
    if getparint("zd", &mut zd) == 0 {
        zd = 1;
    }
    if verbose != 0 {
        warn!("depth decimation ccoefficent: %d ", zd);
    }
    let mut press = 1i32;
    if getparint("press", &mut press) == 0 {
        press = 1;
    }
    if press != 0 && press != 1 {
        err!("press must equal 0 or 1");
    }
    if verbose != 0 {
        if press == 1 {
            warn!("program will output pressure values");
        } else if press == 0 {
            warn!("program will output particle velocity values");
        }
    }

    // Reject coordinates outside the model before the C source injection
    // and receiver recording would access beyond their allocated arrays.
    if isz < 0 || isz >= nz {
        err!("source index %d outside [0, %d)", isz, nz);
    }
    if irz < 0 || irz >= nz {
        err!("receiver index %d outside [0, %d)", irz, nz);
    }

    /* allocate space */
    let mut p = vec![0.0f32; nz as usize];
    let mut v = vec![0.0f32; nz as usize];
    let mut rv = vec![0.0f32; nz as usize];
    let mut rd = vec![0.0f32; nz as usize];
    let mut rd1_5 = vec![0.0f32; nz as usize];

    /* read velocity file */
    if !velfile.is_empty() {
        velocityfp =
            fopen(&velfile, "r").unwrap_or_else(|| err!("cannot open velfile=%s ", velfile));
    }
    if efread(&mut rv, &velocityfp) != nz as usize {
        err!("cannot read %d velocity values ", nz);
    }

    /* read density file and linearly interpolate on correct location */
    if !dfile.is_empty() {
        let densityfp = fopen(&dfile, "r").unwrap_or_else(|| err!("cannot open dfile=%s ", dfile));
        if fread(&mut rd, &densityfp) != nz as usize {
            err!("error reading dfile %s", dfile);
        }
        fclose(&densityfp);
    } else {
        for iz in 0..nz as usize {
            rd[iz] = 2500.0;
        }
    }
    for iz in 0..(nz - 1) as usize {
        rd1_5[iz] = (rd[iz] + rd[iz + 1]) / 2.0f32;
    }
    rd1_5[(nz - 1) as usize] = rd[(nz - 1) as usize];

    /* time step computation */
    let mut vmax = 0.0f32;
    for iz in 0..nz as usize {
        if rv[iz] > vmax {
            vmax = rv[iz];
        }
    }
    if verbose != 0 {
        warn!("vmax= %f ", vmax);
    }
    let dt = (dz as f64 / 1.414 / vmax as f64 / 2.0) as f32;
    if verbose != 0 {
        warn!("time step dt= %f ", dt);
    }

    /* maximum number of iterations */
    if nt == 0 {
        nt = (1.0f32 + tmax / dt) as i32;
    }
    if verbose != 0 {
        warn!("number of time steps nt= %d ", nt);
    }
    if nint((nt as f32 / td as f32) as f64).wrapping_add(1) > SU_NFLTS as i32 {
        err!("too many time steps. Increase td ?");
    }

    /* source parameter computation */
    let alpha = (2.0 * 9.8696 * freq as f64 * freq as f64) as f32;
    let mut epst0 = 0.0f32;
    if styp == 0 || styp == 2 {
        epst0 = (source(0.0, styp, dt, dz, 0.0, alpha) as f64 / 1e4).abs() as f32;
    } else if styp == 1 {
        epst0 = ((source(
            (1.0 / ((2.0f32 * alpha) as f64).sqrt()) as f32,
            styp,
            dt,
            dz,
            0.0,
            alpha,
        ) as f64)
            .abs()
            / 1e4) as f32;
    }
    if verbose != 0 {
        warn!("epst0 = %f ", epst0);
    }

    /* time shift to get a t0 centered source */
    let mut t = tmax + dt;
    loop {
        t = t - dt;
        if !((source(t, styp, dt, dz, 0.0, alpha) as f64).abs() < epst0 as f64) {
            break;
        }
    }
    let t0 = t;
    let ies = (2.0f32 * t / dt) as i32;
    if verbose != 0 {
        warn!("time shift t0 = %f s", t0);
    }

    /* array initialization */
    for iz in 0..nz as usize {
        v[iz] = 0.0;
        p[iz] = 0.0;
    }
    if !wfile.is_empty() {
        // The C code passes a null FILE pointer to fputtr/fclose on failure.
        // Report this undefined-input path through SU's error handler instead.
        wavefp = Some(fopen(&wfile, "w").unwrap_or_else(|| err!("cannot open wfile=%s ", wfile)));
        snapsh.d1 = dz * zd as f32;
        snapsh.f1 = fz;
        // arm64 C integer division by zero returns zero.
        snapsh.ns = nz.checked_div(zd).unwrap_or(0).wrapping_add(1) as u16;
        snapsh.d2 = dt * td as f32;
        snapsh.f2 = 0.0;
    }

    /* propagation computation */
    let mut itsis = 0usize;
    for it in 0..=nt {
        t = it as f32 * dt;
        if abs[0] == 1 {
            p[0] = (p[0] * (1.0f32 - rv[0] * dt / dz)
                + 2.0f32 * rd[0] * rv[0] * rv[0] * dt / dz * v[0])
                / (1.0f32 + rv[0] * dt / dz);
        } else {
            p[0] = 0.0;
        }
        for iz in 1..nz as usize {
            p[iz] = p[iz] + rd[iz] * rv[iz] * rv[iz] * dt / dz * (v[iz] - v[iz - 1]);
        }
        if abs[1] != 1 {
            p[(nz - 1) as usize] = 0.0;
        }
        if it < ies {
            p[isz as usize] = p[isz as usize] + source(t, styp, dt, dz, t0, alpha);
        }
        for iz in 0..(nz - 1) as usize {
            v[iz] = v[iz] + dt / rd1_5[iz] / dz * (p[iz + 1] - p[iz]);
        }
        let iz = (nz - 1) as usize;
        if abs[1] != 1 {
            v[iz] = 0.0;
        } else {
            v[iz] = ((rd1_5[iz] * dz - dt * rd[iz] * rv[iz]) * v[iz] - 2.0f32 * dt * p[iz])
                / (rd1_5[iz] * dz + dt * rd[iz] * rv[iz]);
        }
        // arm64 C remainder uses a zero quotient when td is zero.
        let rem = it.wrapping_sub(it.checked_div(td).unwrap_or(0).wrapping_mul(td));
        if rem == 0 {
            if press == 1 {
                sismo.data[itsis] = p[irz as usize];
            } else {
                sismo.data[itsis] = v[irz as usize];
            }
            itsis += 1;
        }
        if let Some(fp) = wavefp.as_ref().filter(|_| rem == 0) {
            for iz in 0..nz.checked_div(zd).unwrap_or(0) {
                if press == 1 {
                    snapsh.data[iz as usize] = p[(iz * zd) as usize];
                } else {
                    snapsh.data[iz as usize] = v[(iz * zd) as usize];
                }
            }
            fputtr(fp, &mut snapsh);
        }
    }
    if let Some(fp) = wavefp.as_ref() {
        fclose(fp);
    }
    // The C conversions use a 32-bit register before storing the short.
    sismo.dt = ((td as f32 * dt) as f64 * 1e6) as u32 as u16;
    sismo.ns = itsis as u16;
    sismo.delrt = (-t0 * 1000.0f32) as i32 as i16;
    sismo.trid = 1; // TREAL
    sismo.tracl = 1;
    if !sfile.is_empty() {
        seisfp = efopen(&sfile, "w");
    }
    if !sfile.is_empty() && !wfile.is_empty() && itsis != 0 && snapsh.ns != sismo.ns {
        // The C fclose(wavefp)/efopen(sfile) reuses the FILE pointer. SU's
        // fputtr keeps nsfirst for that pointer, so the receiver header says
        // itsis while the payload still contains snapsh.ns samples. CFile
        // does not reuse identities. Seed the core writer on a scratch stream
        // with the snapshot count, then overwrite it with the receiver trace.
        // Copy the resulting bytes so sfile can also be a pipe or /dev/stdout.
        // All header encoding and sample-count handling remain in core fputtr.
        let fp = etmpfile();
        let ns = sismo.ns;
        sismo.ns = snapsh.ns;
        fputtr(&fp, &mut sismo);
        erewind(&fp);
        sismo.ns = ns;
        fputtr(&fp, &mut sismo);
        erewind(&fp);
        let mut bytes = vec![0u8; HDRBYTES + 4 * snapsh.ns as usize];
        efread(&mut bytes, &fp);
        efwrite(&bytes, &seisfp);
        efclose(&fp);
    } else {
        fputtr(&seisfp, &mut sismo);
    }
    exit(CWP_Exit());
}

fn source(t: f32, styp: i32, dt: f32, dz: f32, t0: f32, alpha: f32) -> f32 {
    let mut sou = 0.0f32;
    let x = -alpha * (t - t0) * (t - t0);
    if x > -40.0f32 {
        if styp == 0 {
            sou = (x as f64).exp() as f32;
        } else if styp == 1 {
            sou = ((-2.0f32 * alpha * (t - t0)) as f64 * (x as f64).exp()) as f32;
        } else if styp == 2 {
            sou = ((2.0f32 * alpha * (1.0f32 + 2.0f32 * x)) as f64 * (x as f64).exp()) as f32;
        }
    } else {
        sou = 0.0;
    }
    ((sou / dz * dt) as f64 * 1e8) as f32
}
