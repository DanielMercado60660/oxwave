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

//! SUREMAC2D - Acoustic 2D Fourier modeling with REM time integration.
//! Port of su/main/synthetics_waveforms_testpatterns/suremac2d.c (Revision 1.7).
//! Credits: University of Hamburg: Ekkehart Tessmer, October 2012.
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_pfa_fft::npfa;
mod numerics;
use numerics::*;

// getrusage is a system call; retain the original diagnostic CPU timing.
#[repr(C)]
#[derive(Default)]
struct Timeval {
    sec: i64,
    usec: i32,
}
#[repr(C)]
#[derive(Default)]
struct Rusage {
    user: Timeval,
    system: Timeval,
    counters: [i64; 14],
}
unsafe extern "C" {
    fn getrusage(who: i32, usage: *mut Rusage) -> i32;
}
fn user_time() -> i64 {
    let mut cput = Rusage::default();
    unsafe {
        getrusage(0, &mut cput);
    }
    cput.user.sec * 1_000_000 + cput.user.usec as i64
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                        ",
    " SUREMAC2D - Acoustic 2D Fourier method modeling with high accuracy     ",
    "             Rapid Expansion Method (REM) time integration              ",
    "                                                                        ",
    " suremac2d [parameters]                                                 ",
    "                                                                        ",
    " Required parameters:                                                   ",
    "                                                                        ",
    " opflag=     0: variable density wave equation                          ",
    "             1: constant density wave equation                          ",
    "             2: non-reflecting wave equation                            ",
    "                                                                        ",
    " nx=         number of grid points in horizontal direction              ",
    " nz=         number of grid points in vertical direction                ",
    " nt=         number of time samples                                     ",
    " dx=         spatial increment in horizontal direction                  ",
    " dz=         spatial increment in vertical direction                    ",
    " dt=         time sample interval in seconds                            ",
    " isx=        grid point # of horizontal source positions                ",
    " isz=        grid point # of vertical source positions                  ",
    "                                                                        ",
    " Optional parameters:                                                   ",
    " fx=0.0      first horizontal coordinate                                ",
    " fz=0.0      first vertical coordinate                                  ",
    " irx=        horizontal grid point # of vertical receiver lines         ",
    " irz=        vertical grid point # of horizontal receiver lines         ",
    " w=0.1       width of spatial source distribution (see notes)           ",
    " sflag=2     source time function                                       ",
    "             0: user supplied source function                           ",
    "             1: impulse (spike at t=0)                                  ",
    "             2: Ricker wavelet                                          ",
    " fmax=       maximum frequency of Ricker (default) wavelet              ",
    " amps=1.0    amplitudes of sources                                      ",
    " prec=0      1: precompute Bessel coefficients b_k (see notes)          ",
    "             2: use precomputed Bessel coefficients b_k                 ",
    " fsflag=0    1: perform run with free surface b.c.                      ",
    " vmaxu=      user-defined maximum velocity                              ",
    " dtsnap=0.0  time interval in seconds of wave field snapshots           ",
    " iabso=1     apply absorbing boundary conditions (0: none)              ",
    " abso=0.1    damping parameter for absorbing boundaries                 ",
    " nbwx=20     horizontal width of absorbing boundary                     ",
    " nbwz=20     vertical width of absorbing boundary                       ",
    " verbose=0   1: show parameters used                                    ",
    "             2: print maximum amplitude at every expansion term         ",
    "                                                                        ",
    " velfile=vel          velocity filename                                 ",
    " densfile=dens        density filename                                  ",
    " sname=wavelet.su     user supplied source time function filename       ",
    " sepxname=sectx.su    x-direction pressure sections filename            ",
    " sepzname=sectz.su    z-direction pressure sections filename            ",
    " snpname=snap.su      pressure snapshot filename                        ",
    " jpfile=stderr        diagnostic output                                 ",
    "                                                                        ",
    " Notes:                                                                 ",
    "  0. The combination of the Fourier method with REM time integration    ",
    "     allows the computation of synthetic seismograms which are free     ",
    "     of numerical grid dispersion. REM has no restriction on the        ",
    "     time step size dt. The Fourier method requires at least two        ",
    "     grid points per shortest wavelength.                               ",
    "  1. nx and nz must be valid numbers for pfafft transform lengths.      ",
    "     nx and nz must be odd numbers (unless opflag=1). For valid         ",
    "     numbers see e.g. numbers in structure 'nctab' in source file       ",
    "     $CWPROOT/src/cwp/lib/pfafft.c.                                     ",
    "  2. Velocities (and densities) are stored as plain C style files       ",
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
    "  6. If opflag!=1 the source should be not a spike in space; the        ",
    "     parameter w determines at which distance (in grid points) from     ",
    "     the source's center the Gaussian weight decays to 10 percent       ",
    "     of its maximum. w=2 may be a reasonable choice; however, the       ",
    "     waveform will be distorted.                                        ",
    "  7. Horizontal and vertical receiver line sections are written to      ",
    "     separate files. Each file can hold more than one line.             ",
    "  8. Parameter vmaxu may need to be chosen larger than the highest      ",
    "     propagation velocity if the modeling run becomes unstable.         ",
    "     This happens if the largest eigenvalue of the modeling             ",
    "     operator L is larger than estimated from the largest velocity      ",
    "     due to variations of the density.                                  ",
    "     In particular if using the variable density acoustic wave          ",
    "     equation the eigenvalues depend also on the density and it is      ",
    "     impossible to estimated the largest eigenvalue analytically.       ",
    "  9. Bessel coefficients can be precomputed (prec=1) and stored on      ",
    "     disk to save CPU time when several shots need to be run.           ",
    "     In this case computation of Bessel coefficients can be skipped     ",
    "     and read from disk file for reuse (prec=2).                        ",
    "     For reuse of Bessel coefficients the user may need to define       ",
    "     the overall maximum velocity (vmaxu).                              ",
    " 10. If snapshots are not required, a spike source (sflag=1) may be     ",
    "     applied and the resulting impulse response seismograms can be      ",
    "     convolved later with a desired wavelet.                            ",
    " 11. The free surface (fsflag=1) does not coincide with the first       ",
    "     vertical grid index (0). It appears to be half a grid spacing      ",
    "     above that position.                                               ",
    "                                                                        ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    /* get parameters */
    let mut opflag: i32 = 0;
    if getparint("opflag", &mut opflag) == 0 {
        err!("opflag required!\n");
    }
    let mut nx: i32 = 0;
    if getparint("nx", &mut nx) == 0 {
        err!("nx required!\n");
    }
    let mut nz: i32 = 0;
    if getparint("nz", &mut nz) == 0 {
        err!("nz required!\n");
    }
    let mut nt: i32 = 0;
    if getparint("nt", &mut nt) == 0 {
        err!("nt required!\n");
    }
    let mut dx: f32 = 0.0;
    if getparfloat("dx", &mut dx) == 0 {
        err!("dx required!\n");
    }
    let mut dz: f32 = 0.0;
    if getparfloat("dz", &mut dz) == 0 {
        err!("dz required!\n");
    }
    let mut dt: f32 = 0.0;
    if getparfloat("dt", &mut dt) == 0 {
        err!("dt required!\n");
    }

    let nsourc = countparval("isx") as usize;
    if nsourc != countparval("isz") as usize {
        err!("isx and isz must have the same number of coordinates.\n");
    }
    let mut isx = vec![0i32; nsourc];
    let mut isz = vec![0i32; nsourc];
    let mut amps = vec![1.0f32; nsourc];
    if getparint("isx", &mut isx) == 0 {
        err!("isx required!\n");
    }
    if getparint("isz", &mut isz) == 0 {
        err!("isz required!\n");
    }
    let mut dtsnap = 0.0f32;
    if getparfloat("dtsnap", &mut dtsnap) == 0 {
        dtsnap = 0.0;
    }
    let mut sflag = 2i32;
    if getparint("sflag", &mut sflag) == 0 {
        sflag = 2;
    }

    // C leaves fmax indeterminate if omitted with a non-Ricker source.
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 && sflag == 2 {
        err!("fmax required if Ricker wavelet selected!\n");
    }
    let mut prec = 0i32;
    if getparint("prec", &mut prec) == 0 {
        prec = 0;
    }
    let mut fsflag = 0i32;
    if getparint("fsflag", &mut fsflag) == 0 {
        fsflag = 0;
    }

    let nsectx = countparval("irz") as usize;
    let mut irz = vec![0i32; nsectx];
    if nsectx != 0 {
        getparint("irz", &mut irz);
        // Reject the unchecked row access in suremac2d.c:659, as in fixed C.
        for &row in &irz {
            if row < 0 || row >= nz {
                err!("Error: irz must be >= 0 and < nz!\n");
            }
        }
    }
    let nsectz = countparval("irx") as usize;
    let mut irx = vec![0i32; nsectz];
    if nsectz != 0 {
        getparint("irx", &mut irx);
    }
    let mut iabso = 1i32;
    if getparint("iabso", &mut iabso) == 0 {
        iabso = 1;
    }
    let mut abso = 0.1f32;
    if getparfloat("abso", &mut abso) == 0 {
        abso = 0.1;
    }
    let mut verbose = 0i32;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let mut nbwx = 20i32;
    if getparint("nbwx", &mut nbwx) == 0 {
        nbwx = 20;
    }
    let mut nbwz = 20i32;
    if getparint("nbwz", &mut nbwz) == 0 {
        nbwz = 20;
    }
    let mut w = 0.1f32;
    if getparfloat("w", &mut w) == 0 {
        w = 0.1;
    }
    let mut fx = 0.0f32;
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.0;
    }
    let mut fz = 0.0f32;
    if getparfloat("fz", &mut fz) == 0 {
        fz = 0.0;
    }
    let mut vmaxu = 0.0f32;
    if getparfloat("vmaxu", &mut vmaxu) == 0 {
        vmaxu = 0.0;
    }

    getparfloat("amps", &mut amps);
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        stderr()
    } else {
        fopen(&jpfile, "w").unwrap_or_else(|| exit(139))
    };
    let mut velfile = String::new();
    if getparstring("velfile", &mut velfile) == 0 {
        velfile = "vel".to_string();
    }
    let mut densfile = String::new();
    if getparstring("densfile", &mut densfile) == 0 {
        densfile = "dens".to_string();
    }
    let mut snpname = String::new();
    if getparstring("snpname", &mut snpname) == 0 {
        snpname = "snap.su".to_string();
    }
    let mut sepxname = String::new();
    if getparstring("sepxname", &mut sepxname) == 0 {
        sepxname = "sectx.su".to_string();
    }
    let mut sepzname = String::new();
    if getparstring("sepzname", &mut sepzname) == 0 {
        sepzname = "sectz.su".to_string();
    }
    let mut sname = String::new();
    if getparstring("sname", &mut sname) == 0 {
        sname = "wavelet.su".to_string();
    }

    let mut nzpad = nz;
    if fsflag == 1 {
        nzpad = npfa((nz as f64 * 1.5) as i32);
        if opflag == 1 {
            nzpad = npfa(nzpad);
        } else {
            while nzpad % 2 == 0 {
                nzpad = npfa(nzpad + 1);
            }
        }
    }
    if verbose != 0 {
        fprintf!(jpfp, "opflag=%d\n", opflag);
        fprintf!(jpfp, "nx=%d  dx=%f\n", nx, dx);
        fprintf!(jpfp, "nz=%d  dz=%f\n", nz, dz);
        fprintf!(jpfp, "fsflag=%d\n", fsflag);
        if fsflag == 1 {
            fprintf!(jpfp, "nzpad=%d\n", nzpad);
        }
        fprintf!(jpfp, "nt=%d  dt=%f\n", nt, dt);
        fprintf!(jpfp, "fx=%f  fz=%f\n", fx, fz);
        fprintf!(jpfp, "w=%f\n", w);
        if nsourc == 1 {
            fprintf!(jpfp, "isx=%d  isz=%d\n", isx[0], isz[0]);
        } else {
            fprintf!(jpfp, "sources: nsourc=%d\n", nsourc);
            fprintf!(jpfp, "  coordinates  amplitudes\n");
            for il in 0..nsourc {
                fprintf!(jpfp, "  (%4d,%4d)  %e\n", isx[il], isz[il], amps[il]);
            }
        }
        fprintf!(jpfp, "sflag=%d\n", sflag);
        fprintf!(jpfp, "fmax=%e\n", fmax);
        fprintf!(jpfp, "prec=%d\n", prec);
        if nsectx != 0 {
            fprintf!(jpfp, "x-receiver lines: nsectx=%d\n", nsectx);
            fprintf!(jpfp, "  z-coordinates:");
            for il in 0..nsectx {
                fprintf!(jpfp, "  %d", irz[il]);
            }
            fprintf!(jpfp, "\n");
        }
        if nsectz != 0 {
            fprintf!(jpfp, "z-receiver lines: nsectx=%d\n", nsectx);
            fprintf!(jpfp, "  x-coordinates:");
            for il in 0..nsectz {
                fprintf!(jpfp, "  %d", irx[il]);
            }
            fprintf!(jpfp, "\n");
        }
        if vmaxu > 0.0 {
            fprintf!(jpfp, "vmaxu=%f\n", vmaxu);
        }
    }
    if npfa(nx) != nx {
        err!("Error: Invalid nx !\n");
    }
    // Preserve C precedence: !nx % 2 tests zero, rather than evenness.
    if nx == 0 && opflag != 1 {
        err!("Error: nx must be odd!\n");
    }
    if npfa(nz) != nz {
        err!("Error: Invalid nz !");
    }
    if nz == 0 && opflag != 1 {
        err!("Error: nz must be odd!\n");
    }
    for l in 0..nsourc {
        if isx[l] - 11 < 0 || isx[l] + 11 > nx - 1 {
            err!(
                "Error: (%d) source box beyond model's boundary in x-direction!\n",
                l
            );
        }
        if isz[l] - 11 < 0 || isz[l] + 11 > nz - 1 {
            err!(
                "Error: (%d) source box beyond model's boundary in z-direction!\n",
                l
            );
        }
    }
    let tmax = nt as f32 * dt;
    if dtsnap < 0.0 {
        err!("stopped: dtsnap must be >= 0!\n");
    }
    let nsnap = if dtsnap > 0.0 {
        nint((tmax / dtsnap) as f64)
    } else {
        0
    };
    if verbose != 0 {
        if nsectx != 0 {
            fprintf!(jpfp, "sepxname=%s\n", sepxname);
        }
        if nsectz != 0 {
            fprintf!(jpfp, "sepzname=%s\n", sepzname);
        }
        fprintf!(jpfp, "nsnap=%d\n", nsnap);
        if nsnap != 0 {
            fprintf!(jpfp, "dtsnap=%e\n", dtsnap);
            fprintf!(jpfp, "snpname=%s\n", snpname);
        }
    }
    if prec == 2 {
        fprintf!(
            jpfp,
            "sflag and fmax are ignored if precomputed Bessel coefficients are used (prec=2).\n"
        );
    }
    if nt < 1 {
        err!("stopped: nt must be > 0!\n");
    }
    let (nx, nz, nzpad, nt, nsnap) = (
        nx as usize,
        nz as usize,
        nzpad as usize,
        nt as usize,
        nsnap as usize,
    );
    let mut bwx = vec![0.0f32; nbwx as usize];
    let mut bwz = vec![0.0f32; nbwz as usize];
    let mut a1 = vec![0.0f32; nx * nz];
    let mut a2 = vec![0.0f32; nx * nz];
    let mut c2rho = vec![0.0f32; nx * nz];
    let mut aux1 = vec![0.0f32; nx * nzpad];
    let mut aux2 = if opflag != 1 {
        vec![0.0f32; nx * nzpad]
    } else {
        Vec::new()
    };
    let mut rhoinv = if opflag == 0 {
        vec![0.0f32; nx * nz]
    } else {
        Vec::new()
    };
    let mut snap = vec![vec![0.0f32; nx * nz]; nsnap];
    let mut sectx = vec![vec![vec![0.0f32; nx]; nt]; nsectx];
    let mut sectz = vec![vec![vec![0.0f32; nz]; nt]; nsectz];
    let mut dtwave = 0.0f32;
    let mut nwav = 0usize;
    let mut wave = Vec::new();
    if sflag == 0 {
        let sfp = fopen(&sname, "r").unwrap_or_else(|| err!("Could not open user wavelet file\n"));
        let mut tri = Trace::new();
        fgettr(&sfp, &mut tri);
        if tri.ns == 0 {
            err!("user wavelet: tri.ns must be > 0!\n");
        }
        if tri.dt == 0 {
            err!("user wavelet: tri.dt must be > 0!\n");
        }
        nwav = tri.ns as usize;
        dtwave = (tri.dt as f64 * 1.0e-6) as f32;
        wave = vec![0.0f32; nwav];
        for it in 0..nwav {
            wave[it] = tri.data[it];
        }
        efclose(&sfp);
        fprintf!(jpfp, "nwav=%5d\n", nwav);
    }
    if sflag == 2 {
        dtwave = dt;
        nwav = (8.0 / (dt * fmax) as f64) as usize;
        wave = vec![0.0f32; nwav];
        for it in 0..nwav {
            wave[it] = fwave((it + 1) as f32 * dt, fmax);
        }
    }
    let vel = fopen(&velfile, "r").unwrap_or_else(|| err!("Could not open velocity file\n"));
    let dens = if opflag == 0 {
        Some(fopen(&densfile, "r").unwrap_or_else(|| err!("Could not open density file\n")))
    } else {
        None
    };
    let (vmin, mut vmax) = stru(nx, nz, &mut c2rho, &mut rhoinv, opflag, &vel, dens.as_ref());
    if iabso != 0 {
        for k in 0..nbwz as usize {
            for i in 0..nx {
                c2rho[(nz - 1 - k) * nx + i] = c2rho[i];
            }
            if opflag == 0 {
                for i in 0..nx {
                    rhoinv[(nz - 1 - k) * nx + i] = rhoinv[i];
                }
            }
        }
    }
    efclose(&vel);
    if let Some(dens) = dens {
        efclose(&dens);
    }
    if verbose != 0 {
        fprintf!(jpfp, "max., min. velocities: %e, %e\n", vmax, vmin);
    }
    if iabso != 0 {
        bwx = damp(nbwx as usize, abso);
        bwz = damp(nbwz as usize, abso);
    }
    if vmax == 0.0 {
        err!("stopped: vmax is zero!\n");
    }
    if (if dx > dz { dx } else { dz }) as f64 > vmin as f64 / (2.0 * fmax as f64) {
        err!("stopped: fmax is too high for spatial sampling!\n");
    }
    if vmaxu > 0.0 {
        if vmax > vmaxu {
            err!("vmaxu is smaller than vmax!\n");
        }
        vmax = vmaxu;
    }
    let pi = (4.0 * 1.0f64.atan()) as f32;
    let r =
        (1.1 * pi as f64 * vmax as f64 * (1.0 / (dx * dx) as f64 + 1.0 / (dz * dz) as f64).sqrt())
            as f32;
    let r2 = r * r;
    let m = (tmax * r) as i32;
    let m2 = m / 2;
    if verbose != 0 {
        fprintf!(jpfp, "r=%e  tmax*r=%e  m=%d\n", r, tmax * r, m);
    }
    let (m, m2) = (m as usize, m2 as usize);
    // The C allocation is contiguous. Extra slots retain zero-filled heap
    // contents on tiny expansions which read the first two terms regardless
    // of m2; row indexing below preserves reads into the following row.
    let mut bessn = vec![0.0f64; nsnap * m2 + 2];
    let mut bestr = vec![0.0f64; nt * m2 + 2];
    let mut btemp = vec![0.0f64; m + 50];
    let user_bk_start = if verbose != 0 { user_time() } else { 0 };
    match prec {
        0 | 1 => {
            if verbose != 0 {
                fprintf!(jpfp, "computing Bessel coefficients ...");
            }
            for it in 0..nsnap {
                let t = (it + 1) as f32 * dtsnap;
                if sflag == 1 {
                    bessel_jn((t * r) as f64, m + 50, &mut btemp);
                    for l in 0..m2 {
                        bessn[it * m2 + l] = btemp[2 * l + 1];
                    }
                } else {
                    lintgr(
                        t,
                        r,
                        dtwave,
                        nwav,
                        m,
                        &wave,
                        &mut bessn[it * m2..(it + 1) * m2],
                    );
                }
            }
            for it in 0..nt {
                let t = (it + 1) as f32 * dt;
                if sflag == 1 {
                    bessel_jn((t * r) as f64, m + 50, &mut btemp);
                    for l in 0..m2 {
                        bestr[it * m2 + l] = btemp[2 * l + 1];
                    }
                } else {
                    lintgr(
                        t,
                        r,
                        dtwave,
                        nwav,
                        m,
                        &wave,
                        &mut bestr[it * m2..(it + 1) * m2],
                    );
                }
            }
            // C scans coefficient extrema here but does not use them.
            for k in 0..m2 {
                let mut amx = 0.0f32;
                let mut amn = 0.0f32;
                for it in 0..nt {
                    amx = if amx as f64 > bestr[it * m2 + k] {
                        amx
                    } else {
                        bestr[it * m2 + k] as f32
                    };
                    amn = if (amn as f64) < bestr[it * m2 + k] {
                        amn
                    } else {
                        bestr[it * m2 + k] as f32
                    };
                }
            }
            if verbose != 0 {
                fprintf!(jpfp, " done.\n");
            }
            if prec == 1 {
                let fpbes = efopen("b_k", "w");
                efwrite(&[vmax], &fpbes);
                efwrite(&[r], &fpbes);
                efwrite(&[m as i32], &fpbes);
                efwrite(&[nsnap as i32], &fpbes);
                efwrite(&[dtsnap], &fpbes);
                efwrite(&[nt as i32], &fpbes);
                efwrite(&[dt], &fpbes);
                for it in 0..nsnap {
                    efwrite(&bessn[it * m2..(it + 1) * m2], &fpbes);
                }
                for it in 0..nt {
                    efwrite(&bestr[it * m2..(it + 1) * m2], &fpbes);
                }
                fclose(&fpbes);
            }
        }
        2 => {
            let fpbes = fopen("b_k", "r")
                .unwrap_or_else(|| err!("Could not open Bessel coefficients file (b_k)!\n"));
            let mut vmax0 = [0.0f32];
            let mut r0 = [0.0f32];
            let mut m0 = [0i32];
            let mut nsnap0 = [0i32];
            let mut dtsnap0 = [0.0f32];
            let mut nt0 = [0i32];
            let mut dt0 = [0.0f32];
            efread(&mut vmax0, &fpbes);
            efread(&mut r0, &fpbes);
            efread(&mut m0, &fpbes);
            efread(&mut nsnap0, &fpbes);
            efread(&mut dtsnap0, &fpbes);
            efread(&mut nt0, &fpbes);
            efread(&mut dt0, &fpbes);
            if vmax != vmax0[0]
                || r != r0[0]
                || m as i32 != m0[0]
                || nsnap as i32 != nsnap0[0]
                || dtsnap != dtsnap0[0]
                || nt as i32 != nt0[0]
                || dt != dt0[0]
            {
                err!("Wrong set of Bessel coefficients in file b_k!\n");
            }
            for it in 0..nsnap {
                efread(&mut bessn[it * m2..(it + 1) * m2], &fpbes);
            }
            for it in 0..nt {
                efread(&mut bestr[it * m2..(it + 1) * m2], &fpbes);
            }
            fclose(&fpbes);
        }
        _ => {}
    }
    let mut user_mod_start = 0;
    if verbose != 0 {
        let current = user_time();
        let user_bk = current - user_bk_start;
        user_mod_start = current;
        fprintf!(
            jpfp,
            "CPU time for Bessel coefficients: %.3f seconds\n",
            user_bk as f32 as f64 / 1000000.0
        );
    }
    if verbose != 0 {
        fprintf!(jpfp, "starting modeling.\n");
    }
    for k in 0..nz {
        for i in 0..nx {
            a1[k * nx + i] = 0.0;
            a2[k * nx + i] = 0.0;
        }
    }
    let gbox = g2d(w, 21, 21);
    for l in 0..nsourc {
        for k in 0..21 {
            for i in 0..21 {
                let indz = isz[l] + k as i32 - 10;
                let indx = isx[l] + i as i32 - 10;
                if gbox[k * 21 + i] as f64 > 1.0e-2 {
                    a1[indz as usize * nx + indx as usize] = gbox[k * 21 + i] * amps[l];
                }
            }
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            a1[k * nx + i] = (a1[k * nx + i] as f64 * (1.0 / (r * dx * dz) as f64)) as f32;
        }
    }
    let order = if opflag == 1 { 2 } else { 1 };
    let mut workx = Derivative::new(nx, dx, order);
    let mut workz = Derivative::new(nzpad, dz, order);

    /* initialize first two expansion terms Q_1 and Q_3 */
    match opflag {
        0 | 2 => {
            let bulk = if opflag == 0 { &rhoinv } else { &c2rho };
            difx(Some(&a1), &mut aux1, bulk, nx, nz, &mut workx, 2);
            difx(None, &mut aux1, &[], nx, nz, &mut workx, 0);
            difz(&a1, &mut aux2, bulk, nx, nz, nzpad, &mut workz, 2, 0);
            difz(&aux2, &mut aux1, &[], nx, nz, nzpad, &mut workz, 1, 1);
            for k in 0..nz {
                for i in 0..nx {
                    aux1[k * nx + i] = c2rho[k * nx + i] * aux1[k * nx + i] / r2;
                    a2[k * nx + i] =
                        (4.0 * aux1[k * nx + i] as f64 + 3.0 * a1[k * nx + i] as f64) as f32;
                }
            }
        }
        1 => {
            difx2(&a1, &mut aux1, &[], nx, nz, &mut workx, 0);
            difz2(&a1, &mut aux1, &[], nx, nz, nzpad, &mut workz, 1);
            for k in 0..nz {
                for i in 0..nx {
                    aux1[k * nx + i] = c2rho[k * nx + i] * aux1[k * nx + i] / r2;
                    a2[k * nx + i] =
                        (4.0 * aux1[k * nx + i] as f64 + 3.0 * a1[k * nx + i] as f64) as f32;
                }
            }
        }
        _ => err!("no such opflag!\n"),
    }
    for it in 0..nsnap {
        for k in 0..nz {
            for i in 0..nx {
                snap[it][k * nx + i] = (2.0 * a1[k * nx + i] as f64 * bessn[it * m2]
                    + 2.0 * a2[k * nx + i] as f64 * bessn[it * m2 + 1])
                    as f32;
            }
        }
    }
    for il in 0..nsectx {
        for it in 0..nt {
            for i in 0..nx {
                sectx[il][it][i] = (2.0 * a1[irz[il] as usize * nx + i] as f64 * bestr[it * m2]
                    + 2.0 * a2[irz[il] as usize * nx + i] as f64 * bestr[it * m2 + 1])
                    as f32;
            }
        }
    }
    for il in 0..nsectz {
        for it in 0..nt {
            for k in 0..nz {
                sectz[il][it][k] = (2.0 * a1[k * nx + irx[il] as usize] as f64 * bestr[it * m2]
                    + 2.0 * a2[k * nx + irx[il] as usize] as f64 * bestr[it * m2 + 1])
                    as f32;
            }
        }
    }
    if iabso != 0 {
        tbc2d(&mut a1, &bwx, &bwz, nbwx, nbwz, nx as i32, nz as i32);
        tbc2d(&mut a2, &bwx, &bwz, nbwx, nbwz, nx as i32, nz as i32);
    }
    /* calculate all other odd expansion terms Q_2n+1 */
    for l in 2..m2 {
        if verbose == 2 {
            let mut amx = 0.0f32;
            for k in 0..nz {
                for i in 0..nx {
                    let tmp = (a1[k * nx + i] as f64).abs() as f32;
                    amx = if amx > tmp { amx } else { tmp };
                }
            }
            fprintf!(jpfp, "l=%4d  max= %e\n", l, amx);
        }
        match opflag {
            0 | 2 => {
                let bulk = if opflag == 0 { &rhoinv } else { &c2rho };
                difx(Some(&a2), &mut aux1, bulk, nx, nz, &mut workx, 2);
                difx(None, &mut aux1, &[], nx, nz, &mut workx, 0);
                difz(&a2, &mut aux2, bulk, nx, nz, nzpad, &mut workz, 2, 0);
                difz(&aux2, &mut aux1, &[], nx, nz, nzpad, &mut workz, 1, 1);
                for k in 0..nz {
                    for i in 0..nx {
                        aux1[k * nx + i] = c2rho[k * nx + i] * aux1[k * nx + i] / r2;
                        aux1[k * nx + i] =
                            (4.0 * aux1[k * nx + i] as f64 + 2.0 * a2[k * nx + i] as f64
                                - a1[k * nx + i] as f64) as f32;
                        a1[k * nx + i] = a2[k * nx + i];
                        a2[k * nx + i] = aux1[k * nx + i];
                    }
                }
            }
            1 => {
                difx2(&a2, &mut aux1, &[], nx, nz, &mut workx, 0);
                difz2(&a2, &mut aux1, &[], nx, nz, nzpad, &mut workz, 1);
                for k in 0..nz {
                    for i in 0..nx {
                        aux1[k * nx + i] = c2rho[k * nx + i] * aux1[k * nx + i] / r2;
                        aux1[k * nx + i] =
                            (4.0 * aux1[k * nx + i] as f64 + 2.0 * a2[k * nx + i] as f64
                                - a1[k * nx + i] as f64) as f32;
                        a1[k * nx + i] = a2[k * nx + i];
                        a2[k * nx + i] = aux1[k * nx + i];
                    }
                }
            }
            _ => err!("no such opflag!\n"),
        }
        for it in 0..nsnap {
            for k in 0..nz {
                for i in 0..nx {
                    snap[it][k * nx + i] = (snap[it][k * nx + i] as f64
                        + 2.0 * a2[k * nx + i] as f64 * bessn[it * m2 + l])
                        as f32;
                }
            }
        }
        for il in 0..nsectx {
            for it in 0..nt {
                for i in 0..nx {
                    sectx[il][it][i] = (sectx[il][it][i] as f64
                        + 2.0 * a2[irz[il] as usize * nx + i] as f64 * bestr[it * m2 + l])
                        as f32;
                }
            }
        }
        for il in 0..nsectz {
            for it in 0..nt {
                for k in 0..nz {
                    sectz[il][it][k] = (sectz[il][it][k] as f64
                        + 2.0 * a2[k * nx + irx[il] as usize] as f64 * bestr[it * m2 + l])
                        as f32;
                }
            }
        }
        if iabso != 0 {
            tbc2d(&mut a1, &bwx, &bwz, nbwx, nbwz, nx as i32, nz as i32);
            tbc2d(&mut a2, &bwx, &bwz, nbwx, nbwz, nx as i32, nz as i32);
        }
    }
    if verbose != 0 {
        let user_mod = user_time() - user_mod_start;
        fprintf!(
            jpfp,
            "CPU time for modeling: %.3f seconds\n",
            user_mod as f32 as f64 / 1000000.0
        );
    }
    /* output results to files: common trace headers */
    let mut tro = Trace::new();
    tro.trid = 1;
    tro.scalco = -10;
    tro.scalel = -10;
    tro.sx = nint(((fx + isx[0] as f32 * dx) as f64) * 10.0);
    tro.sdepth = nint(((fz + isz[0] as f32 * dz) as f64) * 10.0);
    tro.ns = nt as u16;
    tro.dt = nint(dt as f64 * 1000000.0) as u16;
    tro.d2 = dx;
    if nsectx != 0 {
        let sepx = efopen(&sepxname, "w");
        for il in 0..nsectx {
            tro.fldr = (il + 1) as i32;
            for i in 0..nx {
                tro.tracl = i as i32;
                tro.tracf = i as i32;
                tro.gx = nint((fx + i as f32 * dx) as f64 * 10.0);
                for it in 0..nt {
                    tro.data[it] = sectx[il][it][i];
                }
                fvputtr(&sepx, &mut tro);
            }
        }
        efclose(&sepx);
    }
    if nsectz != 0 {
        let sepz = efopen(&sepzname, "w");
        for il in 0..nsectz {
            tro.fldr = (il + 1) as i32;
            for k in 0..nz {
                tro.tracl = k as i32;
                tro.tracf = k as i32;
                tro.gx = nint((fz + k as f32 * dz) as f64 * 10.0);
                for it in 0..nt {
                    tro.data[it] = sectz[il][it][k];
                }
                fvputtr(&sepz, &mut tro);
            }
        }
        efclose(&sepz);
    }
    if nsnap != 0 {
        let snp = efopen(&snpname, "w");
        let mut sno = Trace::new();
        sno.ns = nz as u16;
        sno.scalco = -10;
        sno.trid = 130;
        sno.scalco = -10;
        sno.scalel = -10;
        sno.sx = nint((fx + isx[0] as f32 * dx) as f64 * 10.0);
        sno.sdepth = nint((fz + isz[0] as f32 * dz) as f64 * 10.0);
        sno.dt = nint(dtsnap as f64 * 1000000.0) as u16;
        sno.d1 = dz;
        sno.d2 = dx;
        for l in 0..nsnap {
            sno.fldr = (l + 1) as i32;
            for i in 0..nx {
                sno.tracl = i as i32;
                sno.tracf = i as i32;
                sno.gx = nint((fx + i as f32 * dx) as f64 * 10.0);
                sno.offset = nint((sno.gx - sno.sx) as f64 / 10.0);
                for k in 0..nz {
                    sno.data[k] = snap[l][k * nx + i];
                }
                fvputtr(&snp, &mut sno);
            }
        }
        efclose(&snp);
    }
    if verbose != 0 {
        fprintf!(jpfp, "finished modeling.\n");
    }
    efclose(&jpfp);
    exit(CWP_Exit());
}
