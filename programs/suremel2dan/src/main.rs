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

//! SUREMEL2DAN - Elastic anisotropic 2D Fourier modeling with REM integration.
//! Port of su/main/synthetics_waveforms_testpatterns/suremel2dan.c, Revision 1.4.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: University of Hamburg: Ekkehart Tessmer, July 2013.
//! Reference: Kosloff, D., Fihlo A.Q., Tessmer, E. and Behle, A., 1989,
//! Numerical solution of the acoustic and elastic wave equations by a new
//! rapid expansion method, Geophysical Prospecting, 37, 383-394.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_pfa_fft::{Complex, npfa, pfacc};
use oxwave_sinc_interp::ints8r;
const NXB: i32 = 21;
const NZB: i32 = 21;
fn intpar(name: &str, default: i32) -> i32 {
    let mut x = default;
    getparint(name, &mut x);
    x
}
fn floatpar(name: &str, default: f32) -> f32 {
    let mut x = default;
    getparfloat(name, &mut x);
    x
}
fn strpar(name: &str, default: &str) -> String {
    let mut x = String::new();
    if getparstring(name, &mut x) == 0 {
        x = default.into();
    }
    x
}
fn open_model(name: &str, message: &str) -> CFile {
    fopen(name, "r").unwrap_or_else(|| err!("%s", message))
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    macro_rules! reqi {
        ($name:literal) => {{
            let mut x = 0;
            if getparint($name, &mut x) == 0 {
                err!(concat!($name, " required!\n"));
            }
            x
        }};
    }
    macro_rules! reqf {
        ($name:literal) => {{
            let mut x = 0.0f32;
            if getparfloat($name, &mut x) == 0 {
                err!(concat!($name, " required!\n"));
            }
            x
        }};
    }
    let amode = reqi!("amode");
    let nx = reqi!("nx");
    let nz = reqi!("nz");
    let nt = reqi!("nt");
    let dx = reqf!("dx");
    let dz = reqf!("dz");
    let dt = reqf!("dt");
    let _nzpad = nz;
    let nsourc = countparval("isx");
    if nsourc != countparval("isz") {
        err!("isx and isz must have the same number of coordinates.\n");
    }
    if nsourc != countparval("styp") {
        err!("isx and styp arrays must be of same size.\n");
    }
    if nsourc != countparval("samp") {
        err!("isx and samp arrays must be of same size.\n");
    }
    let mut isx = vec![0; nsourc as usize];
    let mut isz = isx.clone();
    let mut styp = isx.clone();
    let mut samp = vec![1.0f32; nsourc as usize];
    if getparint("isx", &mut isx) == 0 {
        err!("isx required!\n");
    }
    if getparint("isz", &mut isz) == 0 {
        err!("isz required!\n");
    }
    if getparint("styp", &mut styp) == 0 {
        err!("styp required!\n");
    }
    if getparfloat("samp", &mut samp) == 0 {
        err!("samp required!\n");
    }
    let mut vmax0 = 0.0f32;
    let mut vmin0 = 0.0f32;
    if amode != 0 {
        vmax0 = reqf!("vmax");
        vmin0 = reqf!("vmin");
    }
    let dtsnap = floatpar("dtsnap", 0.);
    let mut nsntyp = 0;
    let mut sntyp = Vec::<i32>::new();
    if dtsnap != 0. {
        nsntyp = countparval("sntyp");
        sntyp = vec![0; nsntyp as usize];
        if getparint("sntyp", &mut sntyp) == 0 {
            err!("sntyp required if dtsnap<>0!\n");
        }
        getparint("sntyp", &mut sntyp);
    }
    let sflag = intpar("sflag", 2);
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 && sflag == 2 {
        err!("fmax required if Ricker wavelet selected!\n");
    }
    let prec = intpar("prec", 0);
    let nsectx = countparval("irz");
    let mut irz = vec![0; nsectx as usize];
    let mut rxtyp = vec![0; nsectx as usize];
    if nsectx != 0 {
        getparint("irz", &mut irz);
        if countparval("rxtyp") == nsectx {
            getparint("rxtyp", &mut rxtyp);
        } else {
            err!("irz and rxtyp arrays must be of same size.\n");
        }
    }
    let nsectz = countparval("irx");
    let mut irx = vec![0; nsectz as usize];
    let mut rztyp = vec![0; nsectz as usize];
    if nsectz != 0 {
        getparint("irx", &mut irx);
        if countparval("rztyp") == nsectz {
            getparint("rztyp", &mut rztyp);
        } else {
            err!("irx and rztyp must be of same size.\n");
        }
    }
    let iabso = intpar("iabso", 1);
    let abso = floatpar("abso", 0.1);
    let verbose = intpar("verbose", 0);
    let nbwx = intpar("nbwx", 20);
    let nbwz = intpar("nbwz", 20);
    let w = floatpar("w", 0.1);
    let fx = floatpar("fx", 0.);
    let fz = floatpar("fz", 0.);
    let vmaxu = floatpar("vmaxu", 0.);
    getparfloat("samp", &mut samp);
    let mut jpfile = String::new();
    let jpfp_file = if getparstring("jpfile", &mut jpfile) == 0 {
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };
    // C leaves a failed fopen unchecked. Preserve errors reached before the
    // first use of jpfp, but report a clean error instead of using a NULL FILE.
    let jpfp = || {
        jpfp_file
            .as_ref()
            .unwrap_or_else(|| err!("Could not open diagnostic file: %s", jpfile))
    };
    let c11file = strpar("c11file", "c11");
    let c13file = strpar("c13file", "c13");
    let c15file = strpar("c15file", "c15");
    let c33file = strpar("c33file", "c33");
    let c35file = strpar("c35file", "c35");
    let c55file = strpar("c55file", "c55");
    let vpfile = strpar("vpfile", "vp");
    let vsfile = strpar("vsfile", "vs");
    let densfile = strpar("densfile", "dens");
    let xsectb = strpar("xsect", "xsect_");
    let zsectb = strpar("zsect", "zsect_");
    let snapb = strpar("snap", "snap_");
    let sname = strpar("sname", "wavelet.su");
    let typ = ["p", "s", "ux", "uz"];
    if verbose != 0 {
        fprintf!(jpfp(), "\namode=%d ", amode);
        if amode != 0 {
            fprintf!(jpfp(), "(anisotropic elastic modeling)\n");
        } else {
            fprintf!(jpfp(), "(isotropic elastic modeling)\n");
        }
        fprintf!(jpfp(), "\nnx=%d  dx=%f\n", nx, dx);
        fprintf!(jpfp(), "nz=%d  dz=%f\n", nz, dz);
        fprintf!(jpfp(), "nt=%d  dt=%f\n", nt, dt);
        fprintf!(jpfp(), "fx=%f  fz=%f\n", fx, fz);
        fprintf!(jpfp(), "\nw=%.3f\n", w);
        fprintf!(jpfp(), "sources: %d\n", nsourc);
        fprintf!(jpfp(), "  x-coordinates:");
        for i in 0..nsourc as usize {
            fprintf!(jpfp(), "%7d", isx[i]);
        }
        fprintf!(jpfp(), "\n");
        fprintf!(jpfp(), "  z-coordinates:");
        for i in 0..nsourc as usize {
            fprintf!(jpfp(), "%7d", isz[i]);
        }
        fprintf!(jpfp(), "\n");
        fprintf!(jpfp(), "  types        :");
        for i in 0..nsourc as usize {
            fprintf!(jpfp(), "%7s", typ[styp[i] as usize]);
        }
        fprintf!(jpfp(), "\n");
        fprintf!(jpfp(), "  amplitudes   :");
        for i in 0..nsourc as usize {
            fprintf!(jpfp(), "%7.3f", samp[i]);
        }
        fprintf!(jpfp(), "\n");
        fprintf!(jpfp(), "\nsflag=%d\n", sflag);
        fprintf!(jpfp(), "fmax=%e\n", fmax);
        fprintf!(jpfp(), "prec=%d\n", prec);
        if vmaxu > 0. {
            fprintf!(jpfp(), "vmaxu=%f\n", vmaxu);
        }
    }
    if npfa(nx) != nx {
        err!("Error: Invalid nx !\n");
    }
    // The original spells this !nx % 2, which only rejects zero.
    if nx == 0 {
        err!("Error: nx must be odd!\n");
    }
    if npfa(nz) != nz {
        err!("Error: Invalid nz !");
    }
    if nz == 0 {
        err!("Error: nz must be odd!\n");
    }
    for l in 0..nsourc as usize {
        if isx[l] - 11 < 0 || isx[l] + 11 > nx - 1 {
            err!(
                "Error: (%d) source box beyond model's boundary in x-direction!\n",
                l as i32
            );
        }
        if isz[l] - 11 < 0 || isz[l] + 11 > nz - 1 {
            err!(
                "Error: (%d) source box beyond model's boundary in z-direction!\n",
                l as i32
            );
        }
    }
    let tmax = nt as f32 * dt;
    if dtsnap < 0. {
        err!("stopped: dtsnap must be >= 0!\n");
    }
    let nsnap = if dtsnap > 0. {
        nint((tmax / dtsnap) as f64)
    } else {
        0
    };
    if verbose != 0 {
        if nsectx != 0 {
            fprintf!(jpfp(), "\nx-receiver lines: %d\n", nsectx);
            fprintf!(jpfp(), "  z-coordinates:");
            for il in 0..nsectx as usize {
                fprintf!(jpfp(), "%5d", irz[il]);
            }
            fprintf!(jpfp(), "\n");
            fprintf!(jpfp(), "  section types:");
            for il in 0..nsectx as usize {
                fprintf!(jpfp(), "%5s", typ[rxtyp[il] as usize]);
            }
            fprintf!(jpfp(), "\n");
        }
        if nsectz != 0 {
            fprintf!(jpfp(), "z-receiver lines: %d\n", nsectz);
            fprintf!(jpfp(), "  x-coordinates:");
            for il in 0..nsectz as usize {
                fprintf!(jpfp(), "%5d", irx[il]);
            }
            fprintf!(jpfp(), "\n");
            fprintf!(jpfp(), "  section types:");
            for il in 0..nsectz as usize {
                fprintf!(jpfp(), "%5s", typ[rztyp[il] as usize]);
            }
            fprintf!(jpfp(), "\n");
        }
        if nsntyp != 0 {
            fprintf!(jpfp(), "number of snapshot types: %d\n", nsntyp);
            fprintf!(jpfp(), "  nsnap=%d\n", nsnap);
            fprintf!(jpfp(), "  snapshot types:");
            for il in 0..nsntyp as usize {
                fprintf!(jpfp(), "%5s", typ[sntyp[il] as usize]);
            }
            fprintf!(jpfp(), "\n");
        }
    }
    if prec == 2 {
        fprintf!(
            jpfp(),
            "sflag and fmax are ignored if precomputed Bessel coefficients are used (prec=2).\n"
        );
    }
    if nt < 1 {
        err!("stopped: nt must be > 0!\n");
    }
    let mut bufz = vec![0.0f32; nz as usize];
    let mut bwx = vec![0.0f32; nbwx as usize];
    let mut bwz = vec![0.0f32; nbwz as usize];
    let mut qux = grid(nx, nz);
    let mut quz = grid(nx, nz);
    let mut q2ux = grid(nx, nz);
    let mut q2uz = grid(nx, nz);
    let mut a1 = grid(nx, nz);
    let mut a2 = grid(nx, nz);
    let mut a3 = grid(nx, nz);
    let mut rhoinv = grid(nx, nz);
    let mut c11 = if amode != 0 { grid(nx, nz) } else { Vec::new() };
    let mut c13 = if amode != 0 { grid(nx, nz) } else { Vec::new() };
    let mut c15 = if amode != 0 { grid(nx, nz) } else { Vec::new() };
    let mut c33 = if amode != 0 { grid(nx, nz) } else { Vec::new() };
    let mut c35 = if amode != 0 { grid(nx, nz) } else { Vec::new() };
    let mut c55 = if amode != 0 { grid(nx, nz) } else { Vec::new() };
    let mut lambda = if amode == 0 { grid(nx, nz) } else { Vec::new() };
    let mut twmy = if amode == 0 { grid(nx, nz) } else { Vec::new() };
    let mut gbox = grid(21, 21);
    let mut snap = vec![vec![grid(nx, nz); nsnap.max(0) as usize]; nsntyp as usize];
    let mut xsect = vec![grid(nx, nt); nsectx as usize];
    let mut zsect = vec![grid(nz, nt); nsectz as usize];
    let mut dtwave = 0.0f32;
    let mut nwav = 0;
    let mut wave = Vec::<f32>::new();
    if sflag == 0 {
        let sfp = open_model(&sname, "Could not open user wavelet file\n");
        let mut tri = Trace::new();
        fgettr(&sfp, &mut tri);
        if tri.ns == 0 {
            err!("user wavelet: tri.ns must be > 0!\n");
        }
        if tri.dt == 0 {
            err!("user wavelet: tri.dt must be > 0!\n");
        }
        nwav = tri.ns as i32;
        dtwave = (tri.dt as f64 * 1.0e-6) as f32;
        wave = tri.data[..nwav as usize].to_vec();
        efclose(&sfp);
        fprintf!(jpfp(), "nwav=%5d\n", nwav);
    }
    if sflag == 2 {
        dtwave = dt;
        nwav = (8.0f64 / (dt * fmax) as f64) as i32;
        wave = vec![0.0; nwav as usize];
        for it in 0..nwav {
            wave[it as usize] = fwave((it + 1) as f32 * dt, fmax);
        }
    }
    let mut vmin = 1.0e30f32;
    let mut vmax = 0.0f32;
    let fpd = open_model(&densfile, "Could not open density file.");
    for i in 0..nx as usize {
        if fread(&mut bufz, &fpd) != nz as usize {
            err!("densfile has not enough items");
        }
        for k in 0..nz as usize {
            rhoinv[k][i] = (1.0f64 / bufz[k] as f64) as f32;
        }
    }
    if amode != 0 {
        let fp11 = open_model(&c11file, "Could not open c11 file.");
        let fp13 = open_model(&c13file, "Could not open c13 file.");
        let fp15 = open_model(&c15file, "Could not open c15 file.");
        let fp33 = open_model(&c33file, "Could not open c33 file.");
        let fp35 = open_model(&c35file, "Could not open c35 file.");
        let fp55 = open_model(&c55file, "Could not open c55 file.");
        for i in 0..nx as usize {
            macro_rules! readc {
                ($fp:expr,$c:ident,$msg:literal) => {
                    if fread(&mut bufz, &$fp) != nz as usize {
                        err!($msg);
                    }
                    for k in 0..nz as usize {
                        $c[k][i] = bufz[k];
                    }
                };
            }
            readc!(fp11, c11, "c11file has not enough items");
            readc!(fp13, c13, "c13file has not enough items");
            readc!(fp15, c15, "c15file has not enough items");
            readc!(fp33, c33, "c33file has not enough items");
            readc!(fp35, c35, "c35file has not enough items");
            readc!(fp55, c55, "c55file has not enough items");
        }
        for fp in [&fp11, &fp13, &fp15, &fp33, &fp35, &fp55] {
            fclose(fp);
        }
        vmax = vmax0;
        vmin = vmin0;
    } else {
        let fps = open_model(&vsfile, "Could not open S-velocity file.");
        let fpp = open_model(&vpfile, "Could not open P-velocity file.");
        for i in 0..nx as usize {
            if fread(&mut bufz, &fps) != nz as usize {
                err!("vsfile has not enough items");
            }
            for k in 0..nz as usize {
                twmy[k][i] = (2.0 * bufz[k] as f64 * bufz[k] as f64 / rhoinv[k][i] as f64) as f32;
                vmin = if bufz[k] < vmin { bufz[k] } else { vmin };
            }
            if fread(&mut bufz, &fpp) != nz as usize {
                err!("vpfile has not enough items");
            }
            for k in 0..nz as usize {
                lambda[k][i] = bufz[k] * bufz[k] / rhoinv[k][i] - twmy[k][i];
                vmax = if bufz[k] > vmax { bufz[k] } else { vmax };
            }
        }
        fclose(&fps);
        fclose(&fpp);
    }
    fclose(&fpd);
    if verbose != 0 {
        fprintf!(jpfp(), "max., min. velocities: %e, %e\n", vmax, vmin);
    }
    if vmax == 0. {
        err!("stopped: vmax is zero!\n");
    }
    if (if dx > dz { dx } else { dz }) as f64 > vmin as f64 / (2. * fmax as f64) {
        err!("stopped: fmax is too high for spatial sampling!\n");
    }
    if iabso != 0 {
        for i in 0..nx as usize {
            for k in (nz - nbwz)..nz {
                let k = k as usize;
                rhoinv[k][i] = rhoinv[0][i];
                if amode != 0 {
                    c11[k][i] = c11[0][i];
                    c13[k][i] = c13[0][i];
                    c15[k][i] = c15[0][i];
                    c33[k][i] = c33[0][i];
                    c35[k][i] = c35[0][i];
                    c55[k][i] = c55[0][i];
                } else {
                    twmy[k][i] = twmy[0][i];
                    lambda[k][i] = lambda[0][i];
                }
            }
        }
        damp(nbwx, abso, &mut bwx);
        damp(nbwz, abso, &mut bwz);
    }
    if vmaxu > 0. {
        if vmax > vmaxu {
            err!("vmaxu is smaller than vmax!\n");
        }
        vmax = vmaxu;
    }
    let pi = (4. * 1.0f64.atan()) as f32;
    let r =
        (1.1 * pi as f64 * vmax as f64 * (1.0 / (dx * dx) as f64 + 1.0 / (dz * dz) as f64).sqrt())
            as f32;
    let m = (tmax * r) as i32;
    let m2 = m / 2;
    if verbose != 0 {
        fprintf!(jpfp(), "r=%e  tmax*r=%e  m=%d\n", r, tmax * r, m);
    }
    let mut bessn = vec![vec![0.0f64; m2 as usize]; nsnap.max(0) as usize];
    let mut bestr = vec![vec![0.0f64; m2 as usize]; nt as usize];
    let mut btemp = vec![0.0f64; (m + 50) as usize];
    let start_bk = cpu_time();
    match prec {
        0 | 1 => {
            if verbose != 0 {
                fprintf!(jpfp(), "computing Bessel coefficients ...");
            }
            for it in 0..nsnap {
                let t = (it + 1) as f32 * dtsnap;
                if sflag == 1 {
                    bessel_jn((t * r) as f64, m + 50, &mut btemp);
                    for l in 0..m2 {
                        bessn[it as usize][l as usize] = btemp[(2 * l + 1) as usize];
                    }
                } else {
                    lintgr(t, r, dtwave, nwav, m, &wave, &mut bessn[it as usize]);
                }
            }
            for it in 0..nt {
                let t = (it + 1) as f32 * dt;
                if sflag == 1 {
                    bessel_jn((t * r) as f64, m + 50, &mut btemp);
                    for l in 0..m2 {
                        bestr[it as usize][l as usize] = btemp[(2 * l + 1) as usize];
                    }
                } else {
                    lintgr(t, r, dtwave, nwav, m, &wave, &mut bestr[it as usize]);
                }
            }
            if verbose != 0 {
                fprintf!(jpfp(), " done.\n");
            }
            if prec == 1 {
                let fpbes = efopen("b_k", "w");
                efwrite(&[vmax], &fpbes);
                efwrite(&[r], &fpbes);
                efwrite(&[m], &fpbes);
                efwrite(&[nsnap], &fpbes);
                efwrite(&[dtsnap], &fpbes);
                efwrite(&[nt], &fpbes);
                efwrite(&[dt], &fpbes);
                for it in 0..nsnap {
                    efwrite(&bessn[it as usize], &fpbes);
                }
                for it in 0..nt {
                    efwrite(&bestr[it as usize], &fpbes);
                }
                fclose(&fpbes);
            }
        }
        2 => {
            let fpbes = open_model("b_k", "Bessel coefficients must be precomputed (prec=1).\n");
            let mut v0 = [0.0f32];
            let mut r0 = [0.0f32];
            let mut m0 = [0i32];
            let mut ns0 = [0i32];
            let mut ds0 = [0.0f32];
            let mut nt0 = [0i32];
            let mut dt0 = [0.0f32];
            efread(&mut v0, &fpbes);
            efread(&mut r0, &fpbes);
            efread(&mut m0, &fpbes);
            efread(&mut ns0, &fpbes);
            efread(&mut ds0, &fpbes);
            efread(&mut nt0, &fpbes);
            efread(&mut dt0, &fpbes);
            if vmax != v0[0]
                || r != r0[0]
                || m != m0[0]
                || nsnap != ns0[0]
                || dtsnap != ds0[0]
                || nt != nt0[0]
                || dt != dt0[0]
            {
                err!("Wrong set of Bessel coefficients in file b_k!\n");
            }
            for it in 0..nsnap {
                efread(&mut bessn[it as usize], &fpbes);
            }
            for it in 0..nt {
                efread(&mut bestr[it as usize], &fpbes);
            }
            fclose(&fpbes);
        }
        _ => {}
    }
    let start_mod = cpu_time();
    if verbose != 0 {
        fprintf!(
            jpfp(),
            "CPU time for Bessel coefficients: %.3f seconds\n",
            ((start_mod - start_bk) as f32) as f64 / 1000000.
        );
        fprintf!(jpfp(), "starting modeling.\n");
    }
    let mut deriv = Derivatives::new(nx, nz, dx, dz);
    // tbc2d initializes its static taper geometry only when it is called.
    let boundary = if iabso != 0 {
        Boundary::new(nbwx, nbwz, nx, nz)
    } else {
        Boundary::new(0, 0, nx, nz)
    };
    let (mut indx, mut indz): (i32, i32);
    let (mut tmp1, mut tmp2, mut tmp3, mut theta, mut amx, mut amz): (
        f32,
        f32,
        f32,
        f32,
        f32,
        f32,
    );
    for k in 0i32..nz {
        for i in 0i32..nx {
            qux[k as usize][i as usize] = 0.0;
            quz[k as usize][i as usize] = 0.0;
        }
        for i in 0i32..nx {
            rhoinv[k as usize][i as usize] =
                ((rhoinv[k as usize][i as usize] as f64) * (4.0f64 / ((r * r) as f64))) as f32;
        }
    }
    g2d(w, 21, 21, &mut gbox);
    // Initialize the first odd expansion term Q_1 from the sources.
    for l in 0i32..nsourc {
        if styp[l as usize] == 0i32 {
            for k in 0i32..NZB {
                indz = (isz[l as usize] + k) - (NZB / 2i32);
                for i in 0i32..NXB {
                    indx = (isx[l as usize] + i) - (NXB / 2i32);
                    a1[indz as usize][indx as usize] =
                        (((-0.5f64) * (gbox[k as usize][i as usize] as f64))
                            * (samp[l as usize] as f64)) as f32;
                }
            }
            deriv.difx(&a1, &mut qux, None, 1i32);
            deriv.difz(&a1, &mut quz, None, 1i32);
        } else if styp[l as usize] == 1i32 {
            for k in 0i32..NZB {
                indz = (isz[l as usize] + k) - (NZB / 2i32);
                for i in 0i32..NXB {
                    indx = (isx[l as usize] + i) - (NXB / 2i32);
                    a1[indz as usize][indx as usize] =
                        (((-0.5f64) * (gbox[k as usize][i as usize] as f64))
                            * (samp[l as usize] as f64)) as f32;
                    a2[indz as usize][indx as usize] = -a1[indz as usize][indx as usize];
                }
            }
            deriv.difx(&a2, &mut quz, None, 1i32);
            deriv.difz(&a1, &mut qux, None, 1i32);
        } else if styp[l as usize] == 2i32 {
            for k in 0i32..NZB {
                indz = (isz[l as usize] + k) - (NZB / 2i32);
                for i in 0i32..NXB {
                    indx = (isx[l as usize] + i) - (NXB / 2i32);
                    qux[indz as usize][indx as usize] = qux[indz as usize][indx as usize]
                        + (gbox[k as usize][i as usize] * samp[l as usize]);
                }
            }
        } else if styp[l as usize] == 3i32 {
            for k in 0i32..NZB {
                indz = (isz[l as usize] + k) - (NZB / 2i32);
                for i in 0i32..NXB {
                    indx = (isx[l as usize] + i) - (NXB / 2i32);
                    quz[indz as usize][indx as usize] = quz[indz as usize][indx as usize]
                        + (gbox[k as usize][i as usize] * samp[l as usize]);
                }
            }
        }
    }
    for k in 0i32..nz {
        for i in 0i32..nx {
            qux[k as usize][i as usize] =
                ((qux[k as usize][i as usize] as f64) * (1.0f64 / (((r * dx) * dz) as f64))) as f32;
            quz[k as usize][i as usize] =
                ((quz[k as usize][i as usize] as f64) * (1.0f64 / (((r * dx) * dz) as f64))) as f32;
        }
    }
    // Strains for Q_3; accumulate shear outputs before converting to stress.
    deriv.difx(&qux, &mut a1, None, 0i32);
    deriv.difz(&quz, &mut a2, None, 0i32);
    deriv.difx(&quz, &mut a3, None, 0i32);
    for il in 0i32..nsntyp {
        if sntyp[il as usize] == 1i32 {
            for it in 0i32..nsnap {
                for k in 0i32..nz {
                    for i in 0i32..nx {
                        snap[il as usize][it as usize][k as usize][i as usize] = ((2.0f64
                            * (a3[k as usize][i as usize] as f64))
                            * bessn[it as usize][0i32 as usize])
                            as f32;
                    }
                }
            }
        }
    }
    for il in 0i32..nsectx {
        if rxtyp[il as usize] == 1i32 {
            for it in 0i32..nt {
                for i in 0i32..nx {
                    xsect[il as usize][it as usize][i as usize] =
                        ((2.0f64 * (a3[irz[il as usize] as usize][i as usize] as f64))
                            * bestr[it as usize][0i32 as usize]) as f32;
                }
            }
        }
    }
    for il in 0i32..nsectz {
        if rztyp[il as usize] == 1i32 {
            for it in 0i32..nt {
                for k in 0i32..nz {
                    zsect[il as usize][it as usize][k as usize] =
                        ((2.0f64 * (a3[k as usize][irx[il as usize] as usize] as f64))
                            * bestr[it as usize][0i32 as usize]) as f32;
                }
            }
        }
    }
    deriv.difz(&qux, &mut a3, None, 1i32);
    for il in 0i32..nsntyp {
        if sntyp[il as usize] == 1i32 {
            for it in 0i32..nsnap {
                for k in 0i32..nz {
                    for i in 0i32..nx {
                        snap[il as usize][it as usize][k as usize][i as usize] =
                            ((snap[il as usize][it as usize][k as usize][i as usize] as f64)
                                - ((a3[k as usize][i as usize] as f64)
                                    * bessn[it as usize][0i32 as usize]))
                                as f32;
                    }
                }
            }
        }
    }
    for il in 0i32..nsectx {
        if rxtyp[il as usize] == 1i32 {
            for it in 0i32..nt {
                for i in 0i32..nx {
                    xsect[il as usize][it as usize][i as usize] =
                        ((xsect[il as usize][it as usize][i as usize] as f64)
                            - ((a3[irz[il as usize] as usize][i as usize] as f64)
                                * bestr[it as usize][0i32 as usize]))
                            as f32;
                }
            }
        }
    }
    for il in 0i32..nsectz {
        if rztyp[il as usize] == 1i32 {
            for it in 0i32..nt {
                for k in 0i32..nz {
                    zsect[il as usize][it as usize][k as usize] =
                        ((zsect[il as usize][it as usize][k as usize] as f64)
                            - ((a3[k as usize][irx[il as usize] as usize] as f64)
                                * bestr[it as usize][0i32 as usize]))
                            as f32;
                }
            }
        }
    }
    if amode != 0 {
        for k in 0i32..nz {
            for i in 0i32..nx {
                tmp1 = ((a1[k as usize][i as usize] * c11[k as usize][i as usize])
                    + (a2[k as usize][i as usize] * c13[k as usize][i as usize]))
                    + (a3[k as usize][i as usize] * c15[k as usize][i as usize]);
                tmp2 = ((a1[k as usize][i as usize] * c13[k as usize][i as usize])
                    + (a2[k as usize][i as usize] * c33[k as usize][i as usize]))
                    + (a3[k as usize][i as usize] * c35[k as usize][i as usize]);
                tmp3 = ((a1[k as usize][i as usize] * c15[k as usize][i as usize])
                    + (a2[k as usize][i as usize] * c35[k as usize][i as usize]))
                    + (a3[k as usize][i as usize] * c55[k as usize][i as usize]);
                a1[k as usize][i as usize] = tmp1;
                a2[k as usize][i as usize] = tmp2;
                a3[k as usize][i as usize] = tmp3;
            }
        }
    } else {
        for k in 0i32..nz {
            for i in 0i32..nx {
                theta = a1[k as usize][i as usize] + a2[k as usize][i as usize];
                a1[k as usize][i as usize] = (lambda[k as usize][i as usize] * theta)
                    + (twmy[k as usize][i as usize] * a1[k as usize][i as usize]);
                a2[k as usize][i as usize] = (lambda[k as usize][i as usize] * theta)
                    + (twmy[k as usize][i as usize] * a2[k as usize][i as usize]);
                a3[k as usize][i as usize] =
                    (((twmy[k as usize][i as usize] * a3[k as usize][i as usize]) as f64) * 0.5f64)
                        as f32;
            }
        }
    }
    for il in 0i32..nsntyp {
        if sntyp[il as usize] == 0i32 {
            for it in 0i32..nsnap {
                for k in 0i32..nz {
                    for i in 0i32..nx {
                        snap[il as usize][it as usize][k as usize][i as usize] = ((0.5f64
                            * ((a1[k as usize][i as usize] + a2[k as usize][i as usize]) as f64))
                            * bessn[it as usize][0i32 as usize])
                            as f32;
                    }
                }
            }
        }
    }
    for il in 0i32..nsectx {
        if rxtyp[il as usize] == 0i32 {
            for it in 0i32..nt {
                for i in 0i32..nx {
                    xsect[il as usize][it as usize][i as usize] = ((0.5f64
                        * ((a1[irz[il as usize] as usize][i as usize]
                            + a2[irz[il as usize] as usize][i as usize])
                            as f64))
                        * bestr[it as usize][0i32 as usize])
                        as f32;
                }
            }
        }
    }
    for il in 0i32..nsectz {
        if rztyp[il as usize] == 1i32 {
            for it in 0i32..nt {
                for k in 0i32..nz {
                    zsect[il as usize][it as usize][k as usize] = ((0.5f64
                        * ((a1[k as usize][irx[il as usize] as usize]
                            + a2[k as usize][irx[il as usize] as usize])
                            as f64))
                        * bestr[it as usize][0i32 as usize])
                        as f32;
                }
            }
        }
    }
    deriv.difx(&a1, &mut q2ux, Some(&rhoinv), 2i32);
    deriv.difz(&a3, &mut q2ux, Some(&rhoinv), 3i32);
    deriv.difz(&a2, &mut q2uz, Some(&rhoinv), 2i32);
    deriv.difx(&a3, &mut q2uz, Some(&rhoinv), 3i32);
    if iabso != 0 {
        boundary.apply(&mut qux, &bwx, &bwz);
        boundary.apply(&mut quz, &bwx, &bwz);
        boundary.apply(&mut q2ux, &bwx, &bwz);
        boundary.apply(&mut q2uz, &bwx, &bwz);
    }
    for k in 0i32..nz {
        for i in 0i32..nx {
            q2ux[k as usize][i as usize] = ((q2ux[k as usize][i as usize] as f64)
                + (3.0f64 * (qux[k as usize][i as usize] as f64)))
                as f32;
            q2uz[k as usize][i as usize] = ((q2uz[k as usize][i as usize] as f64)
                + (3.0f64 * (quz[k as usize][i as usize] as f64)))
                as f32;
        }
    }
    for il in 0i32..nsntyp {
        if sntyp[il as usize] == 2i32 {
            for it in 0i32..nsnap {
                for k in 0i32..nz {
                    for i in 0i32..nx {
                        snap[il as usize][it as usize][k as usize][i as usize] =
                            ((qux[k as usize][i as usize] as f64)
                                * bessn[it as usize][0i32 as usize])
                                as f32;
                    }
                }
            }
        }
        if sntyp[il as usize] == 3i32 {
            for it in 0i32..nsnap {
                for k in 0i32..nz {
                    for i in 0i32..nx {
                        snap[il as usize][it as usize][k as usize][i as usize] =
                            ((quz[k as usize][i as usize] as f64)
                                * bessn[it as usize][0i32 as usize])
                                as f32;
                    }
                }
            }
        }
    }
    for il in 0i32..nsectx {
        if rxtyp[il as usize] == 2i32 {
            for it in 0i32..nt {
                for i in 0i32..nx {
                    xsect[il as usize][it as usize][i as usize] =
                        ((qux[irz[il as usize] as usize][i as usize] as f64)
                            * bestr[it as usize][0i32 as usize]) as f32;
                }
            }
        }
        if rxtyp[il as usize] == 3i32 {
            for it in 0i32..nt {
                for i in 0i32..nx {
                    xsect[il as usize][it as usize][i as usize] =
                        ((quz[irz[il as usize] as usize][i as usize] as f64)
                            * bestr[it as usize][0i32 as usize]) as f32;
                }
            }
        }
    }
    for il in 0i32..nsectz {
        if rztyp[il as usize] == 2i32 {
            for it in 0i32..nt {
                for k in 0i32..nz {
                    zsect[il as usize][it as usize][k as usize] =
                        ((qux[k as usize][irx[il as usize] as usize] as f64)
                            * bestr[it as usize][0i32 as usize]) as f32;
                }
            }
        }
        if rztyp[il as usize] == 3i32 {
            for it in 0i32..nt {
                for k in 0i32..nz {
                    zsect[il as usize][it as usize][k as usize] =
                        ((quz[k as usize][irx[il as usize] as usize] as f64)
                            * bestr[it as usize][0i32 as usize]) as f32;
                }
            }
        }
    }
    // Recurrence for all remaining odd expansion terms Q_(2n+1).
    for l in 1i32..m2 {
        if verbose == 2i32 {
            amx = 0.0f64 as f32;
            amz = 0.0f64 as f32;
            for k in 0i32..nz {
                for i in 0i32..nx {
                    amx = (if (q2ux[k as usize][i as usize] as f64).abs() > (amx as f64) {
                        (q2ux[k as usize][i as usize] as f64).abs()
                    } else {
                        amx as f64
                    }) as f32;
                    amz = (if (q2uz[k as usize][i as usize] as f64).abs() > (amz as f64) {
                        (q2uz[k as usize][i as usize] as f64).abs()
                    } else {
                        amz as f64
                    }) as f32;
                }
            }
            fprintf!(&stderr(), "l=%5d    amx=%e  amz=%e\n", l, amx, amz);
        }
        for k in 0i32..nz {
            for i in 0i32..nx {
                a1[k as usize][i as usize] = q2ux[k as usize][i as usize];
                a2[k as usize][i as usize] = q2uz[k as usize][i as usize];
                q2ux[k as usize][i as usize] = ((2.0f64 * (q2ux[k as usize][i as usize] as f64))
                    - (qux[k as usize][i as usize] as f64))
                    as f32;
                q2uz[k as usize][i as usize] = ((2.0f64 * (q2uz[k as usize][i as usize] as f64))
                    - (quz[k as usize][i as usize] as f64))
                    as f32;
                qux[k as usize][i as usize] = a1[k as usize][i as usize];
                quz[k as usize][i as usize] = a2[k as usize][i as usize];
            }
        }
        deriv.difx(&qux, &mut a1, None, 0i32);
        deriv.difz(&quz, &mut a2, None, 0i32);
        deriv.difx(&quz, &mut a3, None, 0i32);
        for il in 0i32..nsntyp {
            if sntyp[il as usize] == 1i32 {
                for it in 0i32..nsnap {
                    for k in 0i32..nz {
                        for i in 0i32..nx {
                            snap[il as usize][it as usize][k as usize][i as usize] =
                                ((snap[il as usize][it as usize][k as usize][i as usize] as f64)
                                    + ((2.0f64 * (a3[k as usize][i as usize] as f64))
                                        * bessn[it as usize][l as usize]))
                                    as f32;
                        }
                    }
                }
            }
        }
        for il in 0i32..nsectx {
            if rxtyp[il as usize] == 1i32 {
                for it in 0i32..nt {
                    for i in 0i32..nx {
                        xsect[il as usize][it as usize][i as usize] =
                            ((xsect[il as usize][it as usize][i as usize] as f64)
                                + ((2.0f64 * (a3[irz[il as usize] as usize][i as usize] as f64))
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        for il in 0i32..nsectz {
            if rztyp[il as usize] == 1i32 {
                for it in 0i32..nt {
                    for k in 0i32..nz {
                        zsect[il as usize][it as usize][k as usize] =
                            ((zsect[il as usize][it as usize][k as usize] as f64)
                                + ((2.0f64 * (a3[k as usize][irx[il as usize] as usize] as f64))
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        deriv.difz(&qux, &mut a3, None, 1i32);
        for il in 0i32..nsntyp {
            if sntyp[il as usize] == 1i32 {
                for it in 0i32..nsnap {
                    for k in 0i32..nz {
                        for i in 0i32..nx {
                            snap[il as usize][it as usize][k as usize][i as usize] =
                                ((snap[il as usize][it as usize][k as usize][i as usize] as f64)
                                    - ((a3[k as usize][i as usize] as f64)
                                        * bessn[it as usize][l as usize]))
                                    as f32;
                        }
                    }
                }
            }
        }
        for il in 0i32..nsectx {
            if rxtyp[il as usize] == 1i32 {
                for it in 0i32..nt {
                    for i in 0i32..nx {
                        xsect[il as usize][it as usize][i as usize] =
                            ((xsect[il as usize][it as usize][i as usize] as f64)
                                - ((a3[irz[il as usize] as usize][i as usize] as f64)
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        for il in 0i32..nsectz {
            if rztyp[il as usize] == 1i32 {
                for it in 0i32..nt {
                    for k in 0i32..nz {
                        zsect[il as usize][it as usize][k as usize] =
                            ((zsect[il as usize][it as usize][k as usize] as f64)
                                - ((a3[k as usize][irx[il as usize] as usize] as f64)
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        if amode != 0 {
            for k in 0i32..nz {
                for i in 0i32..nx {
                    tmp1 = ((a1[k as usize][i as usize] * c11[k as usize][i as usize])
                        + (a2[k as usize][i as usize] * c13[k as usize][i as usize]))
                        + (a3[k as usize][i as usize] * c15[k as usize][i as usize]);
                    tmp2 = ((a1[k as usize][i as usize] * c13[k as usize][i as usize])
                        + (a2[k as usize][i as usize] * c33[k as usize][i as usize]))
                        + (a3[k as usize][i as usize] * c35[k as usize][i as usize]);
                    tmp3 = ((a1[k as usize][i as usize] * c15[k as usize][i as usize])
                        + (a2[k as usize][i as usize] * c35[k as usize][i as usize]))
                        + (a3[k as usize][i as usize] * c55[k as usize][i as usize]);
                    a1[k as usize][i as usize] = tmp1;
                    a2[k as usize][i as usize] = tmp2;
                    a3[k as usize][i as usize] = tmp3;
                }
            }
        } else {
            for k in 0i32..nz {
                for i in 0i32..nx {
                    theta = a1[k as usize][i as usize] + a2[k as usize][i as usize];
                    a1[k as usize][i as usize] = (lambda[k as usize][i as usize] * theta)
                        + (twmy[k as usize][i as usize] * a1[k as usize][i as usize]);
                    a2[k as usize][i as usize] = (lambda[k as usize][i as usize] * theta)
                        + (twmy[k as usize][i as usize] * a2[k as usize][i as usize]);
                    a3[k as usize][i as usize] =
                        (((twmy[k as usize][i as usize] * a3[k as usize][i as usize]) as f64)
                            * 0.5f64) as f32;
                }
            }
        }
        for il in 0i32..nsntyp {
            if sntyp[il as usize] == 0i32 {
                for it in 0i32..nsnap {
                    for k in 0i32..nz {
                        for i in 0i32..nx {
                            snap[il as usize][it as usize][k as usize][i as usize] =
                                ((snap[il as usize][it as usize][k as usize][i as usize] as f64)
                                    + ((0.5f64
                                        * ((a1[k as usize][i as usize] + a2[k as usize][i as usize])
                                            as f64))
                                        * bessn[it as usize][l as usize]))
                                    as f32;
                        }
                    }
                }
            }
        }
        for il in 0i32..nsectx {
            if rxtyp[il as usize] == 0i32 {
                for it in 0i32..nt {
                    for i in 0i32..nx {
                        xsect[il as usize][it as usize][i as usize] =
                            ((xsect[il as usize][it as usize][i as usize] as f64)
                                + ((0.5f64
                                    * ((a1[irz[il as usize] as usize][i as usize]
                                        + a2[irz[il as usize] as usize][i as usize])
                                        as f64))
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        for il in 0i32..nsectz {
            if rztyp[il as usize] == 1i32 {
                for it in 0i32..nt {
                    for k in 0i32..nz {
                        zsect[il as usize][it as usize][k as usize] =
                            ((zsect[il as usize][it as usize][k as usize] as f64)
                                + ((0.5f64
                                    * ((a1[k as usize][irx[il as usize] as usize]
                                        + a2[k as usize][irx[il as usize] as usize])
                                        as f64))
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        deriv.difx(&a1, &mut q2ux, Some(&rhoinv), 3i32);
        deriv.difz(&a3, &mut q2ux, Some(&rhoinv), 3i32);
        deriv.difz(&a2, &mut q2uz, Some(&rhoinv), 3i32);
        deriv.difx(&a3, &mut q2uz, Some(&rhoinv), 3i32);
        if iabso != 0 {
            boundary.apply(&mut qux, &bwx, &bwz);
            boundary.apply(&mut quz, &bwx, &bwz);
            boundary.apply(&mut q2ux, &bwx, &bwz);
            boundary.apply(&mut q2uz, &bwx, &bwz);
        }
        for il in 0i32..nsntyp {
            if sntyp[il as usize] == 2i32 {
                for it in 0i32..nsnap {
                    for k in 0i32..nz {
                        for i in 0i32..nx {
                            snap[il as usize][it as usize][k as usize][i as usize] =
                                ((snap[il as usize][it as usize][k as usize][i as usize] as f64)
                                    + ((qux[k as usize][i as usize] as f64)
                                        * bessn[it as usize][l as usize]))
                                    as f32;
                        }
                    }
                }
            }
            if sntyp[il as usize] == 3i32 {
                for it in 0i32..nsnap {
                    for k in 0i32..nz {
                        for i in 0i32..nx {
                            snap[il as usize][it as usize][k as usize][i as usize] =
                                ((snap[il as usize][it as usize][k as usize][i as usize] as f64)
                                    + ((quz[k as usize][i as usize] as f64)
                                        * bessn[it as usize][l as usize]))
                                    as f32;
                        }
                    }
                }
            }
        }
        for il in 0i32..nsectx {
            if rxtyp[il as usize] == 2i32 {
                for it in 0i32..nt {
                    for i in 0i32..nx {
                        xsect[il as usize][it as usize][i as usize] =
                            ((xsect[il as usize][it as usize][i as usize] as f64)
                                + ((qux[irz[il as usize] as usize][i as usize] as f64)
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
            if rxtyp[il as usize] == 3i32 {
                for it in 0i32..nt {
                    for i in 0i32..nx {
                        xsect[il as usize][it as usize][i as usize] =
                            ((xsect[il as usize][it as usize][i as usize] as f64)
                                + ((quz[irz[il as usize] as usize][i as usize] as f64)
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
        for il in 0i32..nsectz {
            if rztyp[il as usize] == 2i32 {
                for it in 0i32..nt {
                    for k in 0i32..nz {
                        zsect[il as usize][it as usize][k as usize] =
                            ((zsect[il as usize][it as usize][k as usize] as f64)
                                + ((qux[k as usize][irx[il as usize] as usize] as f64)
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
            if rztyp[il as usize] == 3i32 {
                for it in 0i32..nt {
                    for k in 0i32..nz {
                        zsect[il as usize][it as usize][k as usize] =
                            ((zsect[il as usize][it as usize][k as usize] as f64)
                                + ((quz[k as usize][irx[il as usize] as usize] as f64)
                                    * bestr[it as usize][l as usize]))
                                as f32;
                    }
                }
            }
        }
    }

    if verbose != 0 {
        let user_mod = cpu_time() - start_mod;
        fprintf!(
            jpfp(),
            "CPU time for modeling: %.3f seconds\n",
            user_mod as f32 as f64 / 1000000.
        );
    }
    // In the supplied macOS binary each filename suffix is formatted by
    // __sprintf_chk into char dum[4]. Even "_0.su" exceeds that buffer and
    // the fortified libc raises SIGTRAP before opening any output file.
    // Preserve the observed binary behavior without an unsafe buffer write.
    #[cfg(target_os = "macos")]
    if nsectx != 0 || nsectz != 0 || nsntyp != 0 {
        unsafe extern "C" {
            fn raise(signal: i32) -> i32;
        }
        unsafe {
            raise(5);
        }
    }
    let mut tro = Trace::new();
    tro.trid = 1;
    tro.scalco = -10;
    tro.scalel = -10;
    tro.sx = nint(((fx + isx[0] as f32 * dx) as f64) * 10.);
    tro.ns = nt as u16;
    tro.dt = nint(dt as f64 * 1000000.) as u16;
    tro.ntr = nx;
    tro.d2 = dx;
    for il in 0..nsectx as usize {
        let filename = format!("{}{}_{}.su", xsectb, typ[rxtyp[il] as usize], il);
        let fp = fopen(&filename, "w").unwrap();
        tro.fldr = il as i32 + 1;
        for i in 0..nx as usize {
            tro.tracl = i as i32 + 1;
            tro.tracr = tro.tracl;
            tro.gx = nint((fx + i as f32 * dx) as f64 * 10.);
            tro.sdepth = nint((fz + isz[0] as f32 * dz) as f64 * 10.);
            tro.offset = nint((tro.gx - tro.sx) as f64 / 10.);
            for it in 0..nt as usize {
                tro.data[it] = xsect[il][it][i];
            }
            fvputtr(&fp, &mut tro);
        }
        efclose(&fp);
    }
    tro.ntr = nz;
    tro.d2 = dz;
    for il in 0..nsectz as usize {
        let filename = format!("{}{}_{}.su", zsectb, typ[rztyp[il] as usize], il);
        let fp = fopen(&filename, "w").unwrap();
        tro.fldr = il as i32 + 1;
        for k in 0..nz as usize {
            tro.tracl = k as i32 + 1;
            tro.tracr = tro.tracl;
            tro.gx = nint((fz + k as f32 * dz) as f64 * 10.);
            tro.sdepth = nint((fz + isz[0] as f32 * dz) as f64 * 10.);
            tro.offset = nint((tro.gx - tro.sdepth) as f64 / 10.);
            for it in 0..nt as usize {
                tro.data[it] = zsect[il][it][k];
            }
            fvputtr(&fp, &mut tro);
        }
        efclose(&fp);
    }
    if nsntyp != 0 {
        let mut sno = Trace::new();
        sno.ns = nz as u16;
        sno.ntr = nx;
        sno.scalco = -10;
        sno.trid = 130;
        sno.scalco = -10;
        sno.scalel = -10;
        sno.sx = nint((fx + isx[0] as f32 * dx) as f64 * 10.);
        sno.sdepth = nint((fz + isz[0] as f32 * dz) as f64 * 10.);
        sno.dt = nint(dtsnap as f64 * 1000000.) as u16;
        sno.d1 = dz;
        sno.d2 = dx;
        for is in 0..nsntyp as usize {
            let filename = format!("{}{}_{}.su", snapb, typ[sntyp[is] as usize], is);
            let fp = fopen(&filename, "w").unwrap();
            let mut sntracl = 0;
            for l in 0..nsnap as usize {
                sno.fldr = l as i32 + 1;
                for i in 0..nx as usize {
                    sntracl += 1;
                    sno.tracl = sntracl;
                    sno.tracr = sntracl;
                    sno.gx = nint((fx + i as f32 * dx) as f64 * 10.);
                    sno.offset = nint((sno.gx - sno.sx) as f64 / 10.);
                    for k in 0..nz as usize {
                        sno.data[k] = snap[is][l][k][i];
                    }
                    fvputtr(&fp, &mut sno);
                }
            }
            efclose(&fp);
        }
    }
    if verbose != 0 {
        fprintf!(jpfp(), "finished modeling.\n");
    }
    if let Some(fp) = jpfp_file.as_ref() {
        efclose(fp);
    } else {
        // This libc returns EOF for fclose(NULL); efclose then reports it.
        err!("%s: efclose: fclose failed", "subcalls.c");
    }
    exit(CWP_Exit());
}

// g2d: Gaussian weights for the spatial source distribution.
fn g2d(mut w: f32, nx: i32, nz: i32, g: &mut Grid) {
    let pi = (4. * 1.0f64.atan()) as f32;
    let i0 = nx / 2;
    let k0 = nz / 2;
    if (w as f64) < 0.1 {
        w = 0.1;
    }
    let sigma = (((-w * w) as f64 / (2. * 0.1f64.ln())).sqrt()) as f32;
    let fac = (1. / (sigma as f64 * (2. * pi as f64).sqrt())) as f32;
    let mut sum = 0.0f32;
    for k in 0..nz {
        for i in 0..nx {
            let r2 = ((i - i0) * (i - i0) + (k - k0) * (k - k0)) as f32;
            let arg = (r2 as f64 / (2. * sigma as f64 * sigma as f64)) as f32;
            g[k as usize][i as usize] = 0.;
            if arg < 80. {
                g[k as usize][i as usize] = (fac as f64 * (-arg as f64).exp()) as f32;
            }
            sum += g[k as usize][i as usize];
        }
    }
    for k in 0..nz as usize {
        for i in 0..nx as usize {
            g[k][i] /= sum;
        }
    }
}

// rk: wavenumbers, including the original 1/n normalization.
fn rk(n: i32, d: f32, ind: i32) -> Vec<f32> {
    let pi = (4. * 1.0f64.atan()) as f32;
    let isign = -(-1.0f64).powf(ind as f32 as f64) as i32;
    let n2 = n / 2;
    let dn = (2. * pi as f64 / (n as f32 * d) as f64) as f32;
    let mut ak = vec![0.; n as usize];
    for i in 0..n {
        let mut c = dn * i as f32;
        if i > n2 {
            c = -dn * (n - i) as f32;
        }
        ak[i as usize] = ((c as f64).powf(ind as f32 as f64) * isign as f64 / n as f64) as f32;
    }
    ak
}

// The C difx/difz static buffers are owned by this modeling run.
struct Derivatives {
    nx: usize,
    nz: usize,
    ax: Vec<Complex>,
    az: Vec<Complex>,
    rkx: Vec<f32>,
    rkz: Vec<f32>,
}
impl Derivatives {
    fn new(nx: i32, nz: i32, dx: f32, dz: f32) -> Self {
        Self {
            nx: nx as usize,
            nz: nz as usize,
            ax: vec![Complex::default(); nx as usize],
            az: vec![Complex::default(); nz as usize],
            rkx: rk(nx, dx, 1),
            rkz: rk(nz, dz, 1),
        }
    }
    fn difx(&mut self, vin: &Grid, vout: &mut Grid, bulk: Option<&Grid>, iadd: i32) {
        let nx = self.nx;
        let nz = self.nz;
        let ieo = nz % 2;
        let num = ((nz + 1) / 2) * 2;
        for k in (0..num).step_by(2) {
            for i in 0..nx {
                self.ax[i].r = vin[k][i];
                self.ax[i].i = if ieo == 0 || k != num - 2 {
                    vin[k + 1][i]
                } else {
                    0.
                };
            }
            pfacc(-1, &mut self.ax);
            for i in 0..nx {
                let tmp = self.ax[i].r * self.rkx[i];
                self.ax[i].r = -self.ax[i].i * self.rkx[i];
                self.ax[i].i = tmp;
            }
            pfacc(1, &mut self.ax);
            for i in 0..nx {
                match iadd {
                    0 => {
                        vout[k][i] = self.ax[i].r;
                        if ieo == 0 || k != num - 2 {
                            vout[k + 1][i] = self.ax[i].i;
                        }
                    }
                    1 => {
                        vout[k][i] += self.ax[i].r;
                        if ieo == 0 || k != num - 2 {
                            vout[k + 1][i] += self.ax[i].i;
                        }
                    }
                    2 => {
                        let b = bulk.unwrap();
                        vout[k][i] = self.ax[i].r * b[k][i];
                        if ieo == 0 || k != num - 2 {
                            vout[k + 1][i] = self.ax[i].i * b[k + 1][i];
                        }
                    }
                    3 => {
                        let b = bulk.unwrap();
                        vout[k][i] += self.ax[i].r * b[k][i];
                        if ieo == 0 || k != num - 2 {
                            vout[k + 1][i] += self.ax[i].i * b[k + 1][i];
                        }
                    }
                    _ => {}
                }
            }
        }
    }
    fn difz(&mut self, vin: &Grid, vout: &mut Grid, bulk: Option<&Grid>, iadd: i32) {
        // nzpad == nz and iload == 0 at every call in main.
        let nx = self.nx;
        let nz = self.nz;
        let ieo = nx % 2;
        let num = ((nx + 1) / 2) * 2;
        for i in (0..num).step_by(2) {
            for k in 0..nz {
                self.az[k].r = vin[k][i];
                self.az[k].i = if ieo == 0 || i != num - 2 {
                    vin[k][i + 1]
                } else {
                    0.
                };
            }
            pfacc(-1, &mut self.az);
            for k in 0..nz {
                let tmp = self.az[k].r * self.rkz[k];
                self.az[k].r = -self.az[k].i * self.rkz[k];
                self.az[k].i = tmp;
            }
            pfacc(1, &mut self.az);
            for k in 0..nz {
                match iadd {
                    0 => {
                        vout[k][i] = self.az[k].r;
                        if ieo == 0 || i != num - 2 {
                            vout[k][i + 1] = self.az[k].i;
                        }
                    }
                    1 => {
                        vout[k][i] += self.az[k].r;
                        if ieo == 0 || i != num - 2 {
                            vout[k][i + 1] += self.az[k].i;
                        }
                    }
                    2 => {
                        let b = bulk.unwrap();
                        vout[k][i] = self.az[k].r * b[k][i];
                        if ieo == 0 || i != num - 2 {
                            vout[k][i + 1] = self.az[k].i * b[k][i + 1];
                        }
                    }
                    3 => {
                        let b = bulk.unwrap();
                        vout[k][i] += self.az[k].r * b[k][i];
                        if ieo == 0 || i != num - 2 {
                            vout[k][i + 1] += self.az[k].i * b[k][i + 1];
                        }
                    }
                    _ => {}
                }
            }
        }
    }
}

fn bessel_jn(x: f64, n: i32, bes: &mut [f64]) {
    if x > n as f64 {
        err!("bessel_jn: n(%d) must be larger than x(%e)!\n", n, x);
    }
    let n = n as usize;
    if x.abs() > 1.0e-130 {
        bes[n - 1] = 0.;
        bes[n - 2] = 1.;
        let mut sum = 0.;
        for i in (0..=n - 3).rev() {
            bes[i] = (2 * (i + 1)) as f64 / x * bes[i + 1] - bes[i + 2];
            if i % 2 == 0 {
                sum = sum + 2. * bes[i];
                if sum.abs() > 1.0e30 {
                    for j in i..n {
                        bes[j] /= sum;
                    }
                    sum = 1.;
                }
            }
        }
        sum = bes[0];
        for i in (2..n).step_by(2) {
            sum += 2. * bes[i];
        }
        for i in 0..n {
            bes[i] /= sum;
        }
    } else {
        bes[0] = 1.;
        for i in 1..n {
            bes[i] = 0.;
        }
    }
}
fn legpol(n: i32, x: f64) -> (f64, f64, f64) {
    let mut p0 = 1.;
    let mut p1 = x;
    for i in 1..n {
        let tmp = ((2 * i + 1) as f64 * x * p1 - i as f64 * p0) / (i + 1) as f64;
        p0 = p1;
        p1 = tmp;
    }
    (p1, p0, n as f64 * (x * p1 - p0) / (x * x - 1.))
}
fn null(n: i32, x0: &mut f64) -> (f64, f64) {
    let eps = f64::EPSILON;
    let mut iter = 0;
    let mut d = 2. * eps;
    let mut p = 0.;
    let mut pm = 0.;
    while d > eps {
        iter += 1;
        if iter > 10 {
            fprintf!(&stderr(), "stop in null: too many iterations!\n");
            exit(1);
        }
        let ps;
        (p, pm, ps) = legpol(n, *x0);
        let x1 = *x0 - p / ps;
        d = (x1 - *x0).abs();
        *x0 = x1;
    }
    (p, pm)
}
fn legcoeffs(n: i32, x: &mut [f64], a: &mut [f64]) {
    let n2 = n / 2;
    let con = 1. - 1. / (8. * n as f64 * n as f64) + 1. / (8. * n as f64 * n as f64 * n as f64);
    for i in 0..n2 {
        x[i as usize] =
            con * (std::f64::consts::PI * (4. * (i + 1) as f64 - 1.) / (4. * n as f64 + 2.)).cos();
    }
    for i in 0..n2 as usize {
        let (_, pm) = null(n, &mut x[i]);
        a[i] = 2. * (1. - x[i] * x[i]) / ((n * n) as f64 * pm * pm);
    }
    if n % 2 != 0 {
        let i = n2 as usize;
        x[i] = 0.;
        let (_, pm, _) = legpol(n, x[i]);
        a[i] = 2. / ((n * n) as f64 * pm * pm);
    }
    for i in 0..n2 as usize {
        x[n as usize - 1 - i] = -x[i];
        a[n as usize - 1 - i] = a[i];
    }
}
fn lintgr(t: f32, r: f32, dt: f32, nt: i32, m: i32, wave: &[f32], sum: &mut [f64]) {
    let mut n = (t * r) as i32;
    n /= 2;
    if n < 20 {
        n = 20;
    }
    let m2 = m / 2;
    let mm = m + 50;
    let mut bes = vec![0.0f64; (mm + 1) as usize];
    let mut x = vec![0.0f64; n as usize];
    let mut a = x.clone();
    let mut u = vec![0.0f32; n as usize];
    let mut xout = u.clone();
    let mut yout = u.clone();
    legcoeffs(n, &mut x, &mut a);
    for j in 0..m2 as usize {
        sum[j] = 0.;
    }
    for i in 0..n as usize {
        u[i] = (0.5 * t as f64 * (x[i] + 1.)) as f32;
        xout[i] = t - u[i];
    }
    let mut wmax = 0.0f64;
    for i in 0..nt as usize {
        if (wave[i] as f64).abs() > wmax {
            wmax = (wave[i] as f64).abs();
        }
    }
    wmax *= 1.0e-6;
    ints8r(dt, dt, wave, 0., 0., &xout, &mut yout);
    for i in 0..n as usize {
        if (yout[i] as f64).abs() < wmax {
            continue;
        }
        let f = yout[i] as f64 * a[i];
        let arg = (u[i] * r) as f64;
        bessel_jn(arg, mm, &mut bes);
        for j in 0..m2 as usize {
            sum[j] += bes[2 * j + 1] * f;
        }
    }
    for j in 0..m2 as usize {
        sum[j] *= 0.5 * t as f64;
    }
}
fn fwave(t: f32, fmax: f32) -> f32 {
    let pi = (4. * 1.0f64.atan()) as f32;
    let pi2 = (2. * pi as f64) as f32;
    let agauss = (0.5 * fmax as f64) as f32;
    let tcut = (1.5 / agauss as f64) as f32;
    let s = (t - tcut) * agauss;
    if (s as f64).abs() < 4. {
        ((-2. * s as f64 * s as f64).exp() * ((pi2 * s) as f64).cos()) as f32
    } else {
        0.
    }
}
fn damp(nbw: i32, abso: f32, bw: &mut [f32]) {
    let pi = (4. * 1.0f64.atan()) as f32;
    let delta = pi / nbw as f32;
    for i in 0..nbw as usize {
        bw[i] = (1. - abso as f64 * (1. + ((i as f32 * delta) as f64).cos()) * 0.5) as f32;
    }
}
struct Boundary {
    ist: Vec<i32>,
    nx: i32,
    nz: i32,
    nwbx: i32,
    nwbz: i32,
}
impl Boundary {
    fn new(nwbx: i32, nwbz: i32, nx: i32, nz: i32) -> Self {
        let mut ist = vec![0; nz as usize];
        let mut ng1 = nwbx;
        let mut ng2 = nwbz;
        let mut ihv = 0;
        if ng1 > ng2 {
            ng1 = nwbz;
            ng2 = nwbx;
            ihv = 1;
        }
        let mut ifl = -ng2;
        let mut istart = nwbx;
        let mut iz = 0;
        if ihv == 1 {
            ist[0] = nwbx;
        }
        for i in 0..ng2 {
            if ifl >= 0 {
                ifl -= ng2;
                iz += 1;
                if ihv == 1 {
                    ist[iz] = nwbx - i;
                }
                istart -= 1;
            }
            ifl += ng1;
            if ihv == 0 {
                ist[i as usize] = istart;
            }
        }
        Self {
            ist,
            nx,
            nz,
            nwbx,
            nwbz,
        }
    }
    fn apply(&self, a: &mut Grid, wbx: &[f32], wbz: &[f32]) {
        let Self {
            nx, nz, nwbx, nwbz, ..
        } = *self;
        let nz0 = nz - 2 * nwbz;
        for k in 0..nz {
            let mut nwb = nwbx;
            if k >= nz0 + nwbz {
                nwb = self.ist[(nz - k - 1) as usize];
            }
            if k >= nz0 && k < nz0 + nwbz {
                nwb = self.ist[(k - nz0) as usize];
            }
            for i in 0..nwb {
                let wb = wbx[i as usize];
                a[k as usize][i as usize] *= wb;
                a[k as usize][(nx - i - 1) as usize] *= wb;
            }
        }
        for k in 0..nwbz {
            let is = self.ist[k as usize];
            let wb = wbz[(nwbz - k - 1) as usize];
            for i in is..nx - is {
                a[(k + nz0) as usize][i as usize] *= wb;
                a[(nz - k - 1) as usize][i as usize] *= wb;
            }
        }
    }
}

// Read only ru_utime from the platform's getrusage structure, as in the C.
fn cpu_time() -> i64 {
    unsafe extern "C" {
        fn getrusage(who: i32, usage: *mut i64) -> i32;
    }
    let mut usage = [0i64; 64];
    unsafe {
        getrusage(0, usage.as_mut_ptr());
    }
    usage[0] * 1000000 + usage[1]
}
type Grid = Vec<Vec<f32>>;
fn grid(nx: i32, nz: i32) -> Grid {
    vec![vec![0.; nx as usize]; nz as usize]
}
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                        ",
    " SUREMEL2DAN - Elastic anisotropic 2D Fourier method modeling with high ",
    "               accuracy Rapid Expansion Method (REM) time integration   ",
    "                                                                        ",
    " suremel2dan [parameters]                                               ",
    "                                                                        ",
    " Required parameters:                                                   ",
    "                                                                        ",
    " nx=         number of grid points in horizontal direction              ",
    " nz=         number of grid points in vertical direction                ",
    " nt=         number of time samples                                     ",
    " dx=         spatial increment in horizontal direction                  ",
    " dz=         spatial increment in vertical direction                    ",
    " dt=         time sample interval in seconds                            ",
    " isx=        grid point # of horizontal source positions                ",
    " isz=        grid point # of vertical source positions                  ",
    " styp=       source types (pressure, shear, single forces)              ",
    " samp=       amplitudes of sources                                      ",
    " amode=      0: isotropic,  1: anisotropic                              ",
    " vmax=       global maximum velocity (only if amode=1)                  ",
    " vmin=       global minimum velocity (only if amode=1)                  ",
    "                                                                        ",
    " Optional parameters:                                                   ",
    " fx=0.0      first horizontal coordinate                                ",
    " fz=0.0      first vertical coordinate                                  ",
    " irx=        horizontal grid point # of vertical receiver lines         ",
    " irz=        vertical grid point # of horizontal receiver lines         ",
    " rxtyp=      types of horizontal receiver lines                         ",
    " rztyp=      types of vertical receivers lines                          ",
    " sntyp=      types of snapshots                                         ",
    "             0: P,  1: S,  2: UX,  3: UZ                                ",
    " w=0.1       width of spatial source distribution (see notes)           ",
    " sflag=2     source time function                                       ",
    "             0: user supplied source function                           ",
    "             1: impulse (spike at t=0)                                  ",
    "             2: Ricker wavelet                                          ",
    " fmax=       maximum frequency of Ricker (default) wavelet              ",
    " amps=1.0    amplitudes of sources                                      ",
    " prec=0      1: precompute Bessel coefficients b_k (see notes)          ",
    "             2: use precomputed Bessel coefficients b_k                 ",
    " vmaxu=      user-defined maximum velocity                              ",
    " dtsnap=0.0  time interval in seconds of wave field snapshots           ",
    " iabso=1     apply absorbing boundary conditions (0: none)              ",
    " abso=0.1    damping parameter for absorbing boundaries                 ",
    " nbwx=20     horizontal width of absorbing boundary                     ",
    " nbwz=20     vertical width of absorbing boundary                       ",
    " verbose=0   1: show parameters used                                    ",
    "             2: print maximum amplitude at every expansion term         ",
    "                                                                        ",
    " c11file=c11       c11 filename                                         ",
    " c13file=c13       c13 filename                                         ",
    " c15file=c15       c15 filename                                         ",
    " c33file=c33       c33 filename                                         ",
    " c35file=c35       c35 filename                                         ",
    " c55file=c55       c55 filename                                         ",
    " vpfile=vp         P-velocity filename                                  ",
    " vsfile=vs         S-velocity filename                                  ",
    " densfile=dens     density filename                                     ",
    "                                                                        ",
    " sname=wavelet.su  user supplied source time function filename          ",
    "                                                                        ",
    " Basenames of seismogram and snapshot files:                            ",
    " xsect=xsect_     x-direction section files basename                    ",
    " zsect=zsect_     z-direction section files basename                    ",
    " snap=snap_       snapshot files basename                               ",
    "                                                                        ",
    " jpfile=stderr        diagnostic output                                 ",
    "                                                                        ",
    " Notes:                                                                 ",
    "  0. The combination of the Fourier method with REM time integration    ",
    "     allows the computation of synthetic seismograms which are free     ",
    "     of numerical grid dispersion. REM has no restriction on the        ",
    "     time step size dt. The Fourier method requires at least two        ",
    "     grid points per shortest wavelength.                               ",
    "  1. nx and nz must be valid numbers for pfafft transform lengths.      ",
    "     nx and nz must be odd numbers. For valid numbers see e.g.          ",
    "     numbers in structure 'nctab' in source file                        ",
    "     $CWPROOT/src/cwp/lib/pfafft.c.                                     ",
    "  2. Velocities and densities are stored as plain C style files         ",
    "     of floats where the fast dimension is along the z-direction.       ",
    "  3. Units must be consistent, e.g. m, s and m/s.                       ",
    "  4. A 20 grid points wide border at the sides and the bottom of        ",
    "     the modeling grid is used for sponge boundary conditions           ",
    "     (default: iabso=1).                                                ",
    "     Source and receiver lines should be placed some (e.g. 10) grid     ",
    "     points away from the absorbing boundaries in order to reduce       ",
    "     reflections due to obliquely incident wavefronts.                  ",
    "  5. Dominant frequency is about fmax/2 (sflag=2), absolute maximum     ",
    "     is delayed by 3/fmax from beginning of wavelet.                    ",
    "  6. If source is not single force (i.e. pressure or shear source)      ",
    "     it should be not a spike in space; the parameter w determines      ",
    "     at which distance (in grid points) from the source's center        ",
    "     the Gaussian weight decays to 10 percent of its maximum.           ",
    "     w=2 may be a reasonable choice; however, the waveform will be      ",
    "     distorted.                                                         ",
    "  7. Horizontal and vertical receiver line sections are written to      ",
    "     separate files. Each file can hold more than one line.             ",
    "  8. Parameter vmaxu may need to be chosen larger than the highest      ",
    "     propagation velocity if the modeling run becomes unstable.         ",
    "     This happens if the largest eigenvalue of the modeling             ",
    "     operator L is larger than estimated from the largest velocity      ",
    "     due to variations of the density.                                  ",
    "  9. Bessel coefficients can be precomputed (prec=1) and stored on      ",
    "     disk to save CPU time when several shots need to be run.           ",
    "     In this case computation of Bessel coefficients can be skipped     ",
    "     and read from disk file for reuse (prec=2).                        ",
    "     For reuse of Bessel coefficients the user may need to define       ",
    "     the overall maximum velocity (vmaxu).                              ",
    " 10. If snapshots are not required, a spike source (sflag=1) may be     ",
    "     applied and the resulting impulse response seismograms can be      ",
    "     convolved later with a desired wavelet.                            ",
    " 11. Output is written to SU style files.                               ",
    "     Basenames of seismogram and snapshot output files will be          ",
    "     extended by the type of the data (p, s, ux, or uz).                ",
    "     Additionally seismogram files will be consecutively numbered.      ",
    "                                                                        ",
    " Caveat:                                                                ",
    "     Time sections and snapshots are kept entirely in memory during     ",
    "     run time. Therefore, lots of time section and snapshots may        ",
    "     eat up a large amount of memory.                                   ",
];
