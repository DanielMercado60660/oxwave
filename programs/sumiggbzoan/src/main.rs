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

//! Anisotropic Gaussian beam migration. Port of su/main/migration_inversion/sumiggbzoan.c (Revision 1.5).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: Tariq Alkhalifah, based on MIGGBZO by Dave Hale.
//! CWP: repackaged as an SU program by John Stockwell, April 2006.

#![allow(non_snake_case, unused_parens)]
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUMIGGBZOAN - MIGration via Gaussian beams ANisotropic media (P-wave)\t",
    "\t\t\t\t\t\t\t\t\t",
    " sumiggbzoan <infile >outfile vfile= nt= nx= nz= [optional parameters]\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " a3333file=\t\tname of file containing a3333(x,z)\t\t",
    " nx=                    number of inline samples (traces)\t\t",
    " nz=                    number of depth samples\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt=tr.dt               time sampling interval\t\t\t\t",
    " dx=tr.d2               inline sampling interval (trace spacing)\t",
    " dz=1.0                 depth sampling interval\t\t\t",
    " fmin=0.025/dt          minimum frequency\t\t\t\t",
    " fmax=10*fmin           maximum frequency\t\t\t\t",
    " amin=-amax             minimum emergence angle; must be > -90 degrees\t",
    " amax=60                maximum emergence angle; must be < 90 degrees\t",
    " bwh=0.5*vavg/fmin      beam half-width; vavg denotes average velocity\t",
    " verbose=0\t\t silent, =1 chatty \t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Files for general anisotropic parameters confined to a vertical plane:",
    " a1111file=\t\tname of file containing a1111(x,z)\t\t",
    " a1133file=          \tname of file containing a1133(x,z)\t\t",
    " a1313file=          \tname of file containing a1313(x,z)\t\t",
    " a1113file=          \tname of file containing a1113(x,z)\t\t",
    " a3313file=          \tname of file containing a3313(x,z)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " For transversely isotropic media Thomsen's parameters could be used:\t",
    " deltafile=\t\tname of file containing delta(x,z)\t\t",
    " epsilonfile=\t\tname of file containing epsilon(x,z)\t\t",
    " a1313file=          \tname of file containing a1313(x,z)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " if anisotropy parameters are not given the program considers\t\t",
    " the medium to be isotropic.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_cwp::cwp::{nint, PI};
use oxwave_pfa_fft::*;
use oxwave_sinc_interp::*;

mod beam;
mod cells;
mod migration;
mod ray;
mod vel;
use beam::{beamData, formBeams};
use cells::{cellBeam, cellTimeAmp};
use migration::miggbzoaniso;
use ray::makeRay;
use vel::Vel2;

fn read_model(file: &str, key: &str, nx: usize, nz: usize) -> Vec<Vec<f32>> {
    let fp = match fopen(file, "r") {
        Some(fp) => fp,
        None => model_error(key, file, false),
    };
    let mut a = vec![0.0f32; nx * nz];
    if fread(&mut a, &fp) != nx * nz {
        model_error(key, file, true);
    }
    fclose(&fp);
    (0..nx)
        .map(|ix| a[ix * nz..(ix + 1) * nz].to_vec())
        .collect()
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut a3333file = String::new();
    let mut nx = 0i32;
    let mut nz = 0i32;
    mustgetparstring("a3333file", &mut a3333file);
    mustgetparint("nx", &mut nx);
    mustgetparint("nz", &mut nz);
    let mut verbose = 0i32;
    getparint("verbose", &mut verbose);
    let mut dt = 0.0f32;
    let mut dx = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        if tr.dt != 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        } else {
            dt = 0.004;
            warn!("tr.dt not set, assuming dt=0.004");
        }
    }
    if getparfloat("dx", &mut dx) == 0 {
        if tr.d2 != 0.0 {
            dx = tr.d2;
        } else {
            dx = 1.0;
            warn!("tr.d2 not set, assuming d2=1.0");
        }
    }
    let mut dz = 1.0f32;
    getparfloat("dz", &mut dz);
    let mut fmin = 0.0f32;
    if getparfloat("fmin", &mut fmin) == 0 {
        fmin = (0.025 / dt as f64) as f32;
    }
    let mut fmax = 0.0f32;
    if getparfloat("fmax", &mut fmax) == 0 {
        fmax = (10.0 * fmin as f64) as f32;
    }
    let mut amax = 60.0f32;
    getparfloat("amax", &mut amax);
    let mut amin = 0.0f32;
    if getparfloat("amin", &mut amin) == 0 {
        amin = -amax;
    }
    if verbose != 0 {
        warn!(" dz = %g", dz);
        warn!(" dx = %g", dx);
        warn!(" dt = %g", dt);
        warn!(" fmin = %g", fmin);
        warn!(" fmax = %g", fmax);
        warn!(" amin = %g", amin);
        warn!(" amax = %g", amax);
    }
    // Store traces in a temporary file while getting the actual trace count.
    let tracefp = etmpfile();
    nx = 0;
    loop {
        nx += 1;
        efwrite(&tr.data[..nt], &tracefp);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    let nxu = nx as usize;
    let nzu = nz as usize;
    let mut fflat = vec![0.0f32; nt * nxu];
    erewind(&tracefp);
    efread(&mut fflat, &tracefp);
    efclose(&tracefp);
    let f: Vec<Vec<f32>> = (0..nxu)
        .map(|ix| fflat[ix * nt..(ix + 1) * nt].to_vec())
        .collect();
    drop(fflat);
    let mut files = vec![String::new(); 7];
    let mut supplied = [false; 7];
    let keys = [
        "a1111file",
        "a1133file",
        "a1313file",
        "a1113file",
        "a3313file",
        "deltafile",
        "epsilonfile",
    ];
    for i in 0..7 {
        supplied[i] = getparstring(keys[i], &mut files[i]) != 0;
    }
    let mut a3333 = read_model(&a3333file, "a3333file", nxu, nzu);
    let mut a1111 = if supplied[0] {
        read_model(&files[0], keys[0], nxu, nzu)
    } else {
        a3333.clone()
    };
    let mut a1133 = if supplied[1] {
        read_model(&files[1], keys[1], nxu, nzu)
    } else {
        let mut a = a3333.clone();
        for ix in 0..nxu {
            for iz in 0..nzu {
                a[ix][iz] = (0.4 * a3333[ix][iz] as f64) as f32;
            }
        }
        a
    };
    let mut a1313 = if supplied[2] {
        read_model(&files[2], keys[2], nxu, nzu)
    } else {
        let mut a = a3333.clone();
        for ix in 0..nxu {
            for iz in 0..nzu {
                a[ix][iz] = (0.3 * a3333[ix][iz] as f64) as f32;
            }
        }
        a
    };
    let mut a1113 = if supplied[3] {
        read_model(&files[3], keys[3], nxu, nzu)
    } else {
        vec![vec![0.0; nzu]; nxu]
    };
    let mut a3313 = if supplied[4] {
        read_model(&files[4], keys[4], nxu, nzu)
    } else {
        vec![vec![0.0; nzu]; nxu]
    };
    if supplied[5] || supplied[6] {
        let delta = if supplied[5] {
            read_model(&files[5], keys[5], nxu, nzu)
        } else {
            vec![vec![0.0; nzu]; nxu]
        };
        let epsilon = if supplied[6] {
            read_model(&files[6], keys[6], nxu, nzu)
        } else {
            vec![vec![0.0; nzu]; nxu]
        };
        for ix in 0..nxu {
            for iz in 0..nzu {
                a1111[ix][iz] = a3333[ix][iz] + 2.0f32 * epsilon[ix][iz] * a3333[ix][iz];
                // sqrt returns double; the subtraction is performed before narrowing.
                a1133[ix][iz] =
                    (((2.0f32 * delta[ix][iz] * a3333[ix][iz] * (a3333[ix][iz] - a1313[ix][iz])
                        + (a3333[ix][iz] - a1313[ix][iz]) * (a3333[ix][iz] - a1313[ix][iz]))
                        as f64)
                        .sqrt()
                        - a1313[ix][iz] as f64) as f32;
            }
        }
    }
    let mut vavg1 = 0.0f32;
    let mut vavg3 = 0.0f32;
    for ix in 0..nxu {
        for iz in 0..nzu {
            a1111[ix][iz] = (a1111[ix][iz] as f64 * 0.25) as f32;
            vavg1 = (vavg1 as f64 + (a1111[ix][iz] as f64).sqrt()) as f32;
            a3333[ix][iz] = (a3333[ix][iz] as f64 * 0.25) as f32;
            vavg3 = (vavg3 as f64 + (a3333[ix][iz] as f64).sqrt()) as f32;
            a1133[ix][iz] = (a1133[ix][iz] as f64 * 0.25) as f32;
            a1313[ix][iz] = (a1313[ix][iz] as f64 * 0.25) as f32;
            a1113[ix][iz] = (a1113[ix][iz] as f64 * 0.25) as f32;
            a3313[ix][iz] = (a3313[ix][iz] as f64 * 0.25) as f32;
        }
    }
    let mut vavg = (vavg1 + vavg3) / 2.0f32;
    vavg /= (nx * nz) as f32;
    let mut bwh = 0.0f32;
    if getparfloat("bwh", &mut bwh) == 0 {
        bwh = vavg / fmin;
    }
    checkpars();
    if verbose != 0 {
        warn!(" vavg = %g", vavg);
        warn!(" bwh = %g", bwh);
    }
    let mut g = vec![vec![0.0f32; nzu]; nxu];
    miggbzoaniso(
        bwh, fmin, fmax, amin, amax, nt as i32, dt, nx, dx, nz, dz, &a1111, &a3333, &a1133, &a1313,
        &a1113, &a3313, &f, &mut g, verbose,
    );
    tr.ns = nz as u16;
    tr.trid = 130;
    tr.d1 = dz;
    tr.d2 = dx;
    for ix in 0..nxu {
        tr.tracl = ix as i32 + 1;
        for iz in 0..nzu {
            tr.data[iz] = g[ix][iz];
        }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}

#[derive(Clone, Copy, Default)]
struct RayStep {
    t: f32,
    x: f32,
    z: f32,
    q1: f32,
    p1: f32,
    q2: f32,
    p2: f32,
    kmah: i32,
    c: f32,
    s: f32,
    v: f32,
    dvdx: f32,
    dvdz: f32,
}
struct Ray {
    rs: Vec<RayStep>,
    ic: i32,
    circles: Vec<Circle>,
}
#[derive(Clone, Copy, Default)]
struct Circle {
    irsf: i32,
    irsl: i32,
    x: f32,
    z: f32,
    r: f32,
}

fn makeCircles(nc: i32, nrs: i32, rs: &[RayStep]) -> Vec<Circle> {
    let mut c = vec![Circle::default(); nc as usize];
    // ARM64 C integer division returns zero if the circle count is zero.
    let nrsc = 1 + (nrs - 1).checked_div(nc).unwrap_or(0);
    for ic in 0..nc {
        let mut irsf = ic * nrsc;
        let mut irsl = irsf + nrsc - 1;
        if irsf >= nrs {
            irsf = nrs - 1;
        }
        if irsl >= nrs {
            irsl = nrs - 1;
        }
        let mut xmin = rs[irsf as usize].x;
        let mut xmax = xmin;
        let mut zmin = rs[irsf as usize].z;
        let mut zmax = zmin;
        for irs in irsf + 1..=irsl {
            if rs[irs as usize].x < xmin {
                xmin = rs[irs as usize].x;
            }
            if rs[irs as usize].x > xmax {
                xmax = rs[irs as usize].x;
            }
            if rs[irs as usize].z < zmin {
                zmin = rs[irs as usize].z;
            }
            if rs[irs as usize].z > zmax {
                zmax = rs[irs as usize].z;
            }
        }
        let x = (0.5 * (xmin + xmax) as f64) as f32;
        let z = (0.5 * (zmin + zmax) as f64) as f32;
        let r = (((x - xmin) * (x - xmin) + (z - zmin) * (z - zmin)) as f64).sqrt() as f32;
        c[ic as usize] = Circle {
            irsf,
            irsl,
            x,
            z,
            r,
        };
    }
    c
}
fn nearestRayStep(ray: &mut Ray, x: f32, z: f32) -> i32 {
    let nrs = ray.rs.len() as i32;
    if ray.circles.is_empty() {
        ray.ic = 0;
        ray.circles = makeCircles((nrs as f64).sqrt() as i32, nrs, &ray.rs);
    }
    let c = &ray.circles;
    let nc = c.len() as i32;
    let mut ic = ray.ic;
    let mut dx = x - c[ic as usize].x;
    let mut dz = z - c[ic as usize].z;
    let mut dmin = (2.0 * (((dx * dx + dz * dz) as f64).sqrt() + c[ic as usize].r as f64)) as f32;
    let mut dsmin = dmin * dmin;
    let mut irsmin = 0;
    let mut jc = ic;
    let mut js = 0;
    for _ in 0..nc {
        dx = x - c[jc as usize].x;
        dz = z - c[jc as usize].z;
        let ds = dx * dx + dz * dz;
        let rdmin = c[jc as usize].r + dmin;
        if ds <= rdmin * rdmin {
            let mut update = false;
            for irs in c[jc as usize].irsf..=c[jc as usize].irsl {
                let xrs = ray.rs[irs as usize].x;
                let zrs = ray.rs[irs as usize].z;
                dx = x - xrs;
                dz = z - zrs;
                let ds = dx * dx + dz * dz;
                if ds < dsmin {
                    dsmin = ds;
                    irsmin = irs;
                    update = true;
                }
            }
            if update {
                dmin = (dsmin as f64).sqrt() as f32;
                ic = jc;
            }
        }
        js = if js > 0 { -js - 1 } else { -js + 1 };
        jc += js;
        if jc < 0 || jc >= nc {
            js = if js > 0 { -js - 1 } else { -js + 1 };
            jc += js;
        }
    }
    ray.ic = ic;
    if irsmin < 0 || irsmin >= nrs {
        fprintf!(stderr(), "irsmin=%d\n", irsmin);
    }
    irsmin
}

struct BeamData {
    ntr: i32,
    dtr: f32,
    ftr: f32,
    nti: i32,
    dti: f32,
    fti: f32,
    cf: Vec<Vec<Complex>>,
}
#[derive(Clone, Copy, Default)]
struct Cell {
    live: i32,
    dead: i32,
    tr: f32,
    ti: f32,
    ar: f32,
    ai: f32,
}
struct Cells<'a> {
    lx: i32,
    mx: i32,
    nx: i32,
    dx: f32,
    fx: f32,
    lz: i32,
    mz: i32,
    nz: i32,
    dz: f32,
    fz: f32,
    live: i32,
    dead: i32,
    wmin: f32,
    lmin: f32,
    cell: Vec<Vec<Cell>>,
    ray: Ray,
    bd: BeamData,
    g: &'a mut [Vec<f32>],
}

fn accBeam(
    ray: Ray,
    fmin: f32,
    lmin: f32,
    nt: i32,
    dt: f32,
    ft: f32,
    f: &[f32],
    nx: i32,
    dx: f32,
    fx: f32,
    nz: i32,
    dz: f32,
    fz: f32,
    g: &mut [Vec<f32>],
) {
    let wmin = (2.0 * PI * fmin as f64) as f32;
    let lx = 8;
    let lz = 8;
    let mx = 2 + (nx - 1) / lx;
    let mz = 2 + (nz - 1) / lz;
    let bd = beamData(wmin, nt, dt, ft, f);
    let jx = nint(((ray.rs[0].x - fx) / dx / lx as f32) as f64);
    let jz = nint(((ray.rs[0].z - fz) / dz / lz as f32) as f64);
    // franuni is only used for opaque visited-cell labels in the C source.
    // Each Rust cell allocation is initialized, so fixed nonzero labels give
    // the same traversal without depending on uninitialized allocator memory.
    let mut cells = Cells {
        lx,
        mx,
        nx,
        dx,
        fx,
        lz,
        mz,
        nz,
        dz,
        fz,
        live: 1,
        dead: 2,
        wmin,
        lmin,
        cell: vec![vec![Cell::default(); mz as usize]; mx as usize],
        ray,
        bd,
        g,
    };
    setCell(&mut cells, jx, jz);
    accCell(&mut cells, jx, jz);
}
fn setCell(cells: &mut Cells, jx: i32, jz: i32) {
    if jx < 0 || jx >= cells.mx || jz < 0 || jz >= cells.mz {
        return;
    }
    if cells.cell[jx as usize][jz as usize].live == cells.live {
        return;
    }
    cells.cell[jx as usize][jz as usize].live = cells.live;
    if cellTimeAmp(cells, jx, jz) != 0 {
        setCell(cells, jx + 1, jz);
        setCell(cells, jx - 1, jz);
        setCell(cells, jx, jz + 1);
        setCell(cells, jx, jz - 1);
    }
}
fn accCell(cells: &mut Cells, jx: i32, jz: i32) {
    if jx < 0 || jx >= cells.mx - 1 || jz < 0 || jz >= cells.mz - 1 {
        return;
    }
    if cells.cell[jx as usize][jz as usize].dead == cells.dead {
        return;
    }
    cells.cell[jx as usize][jz as usize].dead = cells.dead;
    if cells.cell[jx as usize][jz as usize].live != cells.live {
        return;
    }
    if cells.cell[jx as usize + 1][jz as usize].live == cells.live
        && cells.cell[jx as usize][jz as usize + 1].live == cells.live
        && cells.cell[jx as usize + 1][jz as usize + 1].live == cells.live
    {
        cellBeam(cells, jx, jz);
    }
    accCell(cells, jx + 1, jz);
    accCell(cells, jx - 1, jz);
    accCell(cells, jx, jz + 1);
    accCell(cells, jx, jz - 1);
}
// The C caller dereferences the NULL ray returned for an invalid takeoff point.
fn segv() -> ! {
    unsafe extern "C" {
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        raise(11);
    }
    exit(139)
}

fn model_error(key: &str, file: &str, reading: bool) -> ! {
    match (key, reading) {
        ("a3333file", false) => err!("error opening a3333file=%s!\n", file),
        ("a3333file", true) => err!("error reading a3333file=%s!\n", file),
        ("a1111file", false) => err!("error opening a1111file=%s!\n", file),
        ("a1111file", true) => err!("error reading a1111file=%s!\n", file),
        ("a1133file", false) => err!("error opening a1133file=%s!\n", file),
        ("a1133file", true) => err!("error reading a1133file=%s!\n", file),
        ("a1313file", false) => err!("error opening a1313file=%s!\n", file),
        ("a1313file", true) => err!("error reading a1313file=%s!\n", file),
        ("a1113file", false) => err!("error opening a1113file=%s!\n", file),
        ("a1113file", true) => err!("error reading a1113file=%s!\n", file),
        ("a3313file", false) => err!("error opening a3313file=%s!\n", file),
        ("a3313file", true) => err!("error reading a3313file=%s!\n", file),
        ("deltafile", false) => err!("error opening deltafile=%s!\n", file),
        ("deltafile", true) => err!("error reading deltafile=%s!\n", file),
        ("epsilonfile", false) => err!("error opening epsilonfile=%s!\n", file),
        ("epsilonfile", true) => err!("error reading epsilonfile=%s!\n", file),
        _ => unreachable!(),
    }
}
