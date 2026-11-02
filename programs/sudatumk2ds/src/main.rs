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

//! SUDATUMK2DS - Kirchhoff datuming of sources for 2D prestack data.
//! Port of su/main/datuming/sudatumk2ds.c (Revision 1.8).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Trino Salinas, 05/01/96, Colorado School of Mines.
//! Based on sukzmig2d.c by Zhenyue Liu, 03/01/95.
//! Subroutines adapted from Dave Hale's modeling library to define
//! topography using cubic splines (Dave Hale, CSM, 09/17/91).

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_cwp::{
    complex::Complex,
    cwp::{PI, abs as cabs, max as cmax, min as cmin},
};
use oxwave_f_cubicspline::csplin;
use oxwave_inverse_interp::yxtoxy;
use oxwave_pfa_fft::{npfaro, pfacr, pfarc};
use std::ops::{Index, IndexMut};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    "SUDATUMK2DS - Kirchhoff datuming of sources for 2D prestack data\t",
    " \t\t(input data are receiver gathers) \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "    sudatumk2ds  infile=  outfile=  [parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " infile=stdin\t\tfile for input seismic traces\t\t\t",
    " outfile=stdout\tfile for common offset migration output  \t",
    " ttfile=\t\tfile for input traveltime tables\t\t",
    "   The following 9 parameters describe traveltime tables:\t\t",
    " fzt=\t\t\tfirst depth sample in traveltime table\t\t",
    " nzt= \t\t\tnumber of depth samples in traveltime table\t",
    " dzt=\t\t\tdepth interval in traveltime table\t\t",
    " fxt=\t\t\tfirst lateral sample in traveltime table\t",
    " nxt=\t\t\tnumber of lateral samples in traveltime table\t",
    " dxt=\t\t\tlateral interval in traveltime table\t\t",
    " fs= \t\t\tx-coordinate of first source\t\t\t",
    " ns= \t\t\tnumber of sources\t\t\t\t",
    " ds= \t\t\tx-coordinate increment of sources\t\t",
    "                                                                       ",
    " fxi=                   x-coordinate of the first surface location      ",
    " dxi=                   horizontal spacing on surface                   ",
    " nxi=                   number of input surface locations               ",
    " sgn=                   Sign of the datuming process (up=-1 or down=1)  ",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt= or from header (dt) \ttime sampling interval of input data\t",
    " ft= or from header (ft) \tfirst time sample of input data\t\t",
    " surf=\"0,0;99999,0\"  The first surface defined the recording surface ",
    " surf=\"0,0;99999,0\"  and the second one, the new datum.              ",
    "                       \"x1,z1;x2,z2;x3,z3;...\"                       ",
    " fzo=fzt\t\tz-coordinate of first point in output trace \t",
    " dzo=0.2*dzt\t\tvertical spacing of output trace \t\t",
    " nzo=5*(nzt-1)+1 \tnumber of points in output trace\t\t",
    " fxso=fxt\t\tx-coordinate of first shot\t \t\t",
    " dxso=0.5*dxt\t\tshot horizontal spacing\t\t \t\t",
    " nxso=2*(nxt-1)+1  \tnumber of shots \t\t\t\t",
    " fxgo=fxt\t\tx-coordinate of first receiver\t\t\t",
    " exgo=fxgo+(nxgo-1)*dxgo\tx-coordinate of the last receiver\t",
    " dxgo=0.5*dxt\t\treceiver horizontal spacing\t\t\t",
    " nxgo=nxso\t\tnumber of receivers per shot\t\t\t",
    " fmax=0.25/dt\t\tfrequency-highcut for input traces\t\t",
    " offmax=99999\t\tmaximum absolute offset allowed in migration \t",
    " aperx=nxt*dxt/2  \tmigration lateral aperature \t\t\t",
    " angmax=60\t\tmigration angle aperature from vertical \t",
    " v0=1500(m/s)\t\treference velocity value at surface\t\t",
    " dvz=0.0  \t\treference velocity vertical gradient\t\t",
    " antiali=1             Antialiase filter (no-filter = 0)               ",
    " jpfile=stderr\t\tjob print file name \t\t\t\t",
    " mtr=100  \t\tprint verbal information at every mtr traces\t",
    " ntr=100000\t\tmaximum number of input traces to be migrated\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. Traveltime tables were generated by program rayt2d (or other ones)\t",
    "    on relatively coarse grids, with dimension ns*nxt*nzt. In the\t",
    "    migration process, traveltimes are interpolated into shot/gephone \t",
    "    positions and output grids.\t\t\t\t\t",
    " 2. Input traces must be SU format and organized in common shot gathers",
    " 3. If the offset value of an input trace is not in the offset array \t",
    "    of output, the nearest one in the array is chosen. \t\t",
    " 4. Amplitudes are computed using the reference velocity profile, v(z),",
    "    specified by the parameters v0= and dvz=.\t\t\t\t",
    " 5. Input traces must specify source and receiver positions via the header",
    "    fields tr.sx and tr.gx. Offset is computed automatically.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// All table dimensions retain the C order: allocation(nz,nx), t[ix][iz].
#[derive(Clone)]
struct Grid {
    nz: usize,
    data: Vec<f32>,
}
impl Grid {
    fn new(nz: i32, nx: i32) -> Self {
        Self {
            nz: nz as usize,
            data: vec![0.0; nz as usize * nx as usize],
        }
    }
}
impl Index<(i32, i32)> for Grid {
    type Output = f32;
    fn index(&self, (ix, iz): (i32, i32)) -> &f32 {
        &self.data[ix as usize * self.nz + iz as usize]
    }
}
impl IndexMut<(i32, i32)> for Grid {
    fn index_mut(&mut self, (ix, iz): (i32, i32)) -> &mut f32 {
        &mut self.data[ix as usize * self.nz + iz as usize]
    }
}
macro_rules! reqi {
    ($name:ident) => {
        let mut $name: i32 = 0;
        if getparint(stringify!($name), &mut $name) == 0 {
            err!(concat!("must specify ", stringify!($name), "!\n"));
        }
    };
}
macro_rules! reqf {
    ($name:ident) => {
        let mut $name: f32 = 0.0;
        if getparfloat(stringify!($name), &mut $name) == 0 {
            err!(concat!("must specify ", stringify!($name), "!\n"));
        }
    };
}
macro_rules! opti {
    ($name:ident, $default:expr) => {
        let mut $name: i32 = 0;
        if getparint(stringify!($name), &mut $name) == 0 {
            $name = $default;
        }
    };
}
macro_rules! optf {
    ($name:ident, $default:expr) => {
        let mut $name: f32 = 0.0;
        if getparfloat(stringify!($name), &mut $name) == 0 {
            $name = $default;
        }
    };
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut infile = "stdin".to_string();
    let mut outfile = "stdout".to_string();
    let infp = if getparstring("infile", &mut infile) == 0 {
        stdin()
    } else {
        fopen(&infile, "r").unwrap_or_else(|| err!("cannot open infile=%s\n", infile))
    };
    let outfp = if getparstring("outfile", &mut outfile) == 0 {
        stdout()
    } else {
        fopen(&outfile, "w").unwrap_or_else(|| exit(139))
    };
    efseeko(&infp, 0, SEEK_CUR);
    efseeko(&outfp, 0, SEEK_CUR);

    let mut ttfile = String::new();
    if getparstring("ttfile", &mut ttfile) == 0 {
        err!("must specify ttfile!\n");
    }
    let ttfp = fopen(&ttfile, "r").unwrap_or_else(|| err!("cannot open ttfile=%s\n", ttfile));
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        stderr()
    } else {
        fopen(&jpfile, "w").unwrap_or_else(|| exit(139))
    };
    let mut tr = Trace::new();
    if fgettr(&infp, &mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    optf!(dt, (tr.dt as f64 / 1000000.0) as f32);
    if (dt as f64) < 0.0000001 {
        err!("dt must be positive!\n");
    }
    optf!(ft, (tr.delrt as f32 as f64 / 1000.0) as f32);
    reqi!(nxt);
    reqf!(fxt);
    reqf!(dxt);
    reqi!(nzt);
    reqf!(fzt);
    reqf!(dzt);
    reqi!(ns);
    reqf!(fs);
    reqf!(ds);
    let ext = fxt + (nxt - 1) as f32 * dxt;
    let ezt = fzt + (nzt - 1) as f32 * dzt;
    let es = fs + (ns - 1) as f32 * ds;
    reqi!(nxi);
    reqf!(fxi);
    reqf!(dxi);
    reqi!(sgn);

    opti!(nxso, (nxt - 1) * 2 + 1);
    optf!(fxso, fxt);
    optf!(dxso, (dxt as f64 * 0.5) as f32);
    opti!(nxgo, (nxt - 1) * 2 + 1);
    optf!(fxgo, fxt);
    optf!(dxgo, (dxt as f64 * 0.5) as f32);
    optf!(exgo, fxgo + (nxgo - 1) as f32 * dxgo);
    opti!(antiali, 1);
    opti!(nzo, (nzt - 1) * 5 + 1);
    optf!(fzo, fzt);
    optf!(dzo, (dzt as f64 * 0.2) as f32);
    fzo *= sgn as f32;
    dzo *= sgn as f32;
    let exso = fxso + (nxso - 1) as f32 * dxso;
    let ezo = fzo + (nzo - 1) as f32 * dzo;
    if fxt > fxso || fxt > fxgo || ext < exso || ext < exgo || fzt > fzo || ezt < ezo {
        err!(" migration output range is out of traveltime table!\n");
    }
    optf!(v0, 1500.0);
    optf!(dvz, 0.0);
    optf!(angmax, 60.0);
    if (angmax as f64) < 0.00001 {
        err!("angmax must be positive!\n");
    }
    optf!(aperx, (0.5 * nxt as f64 * dxt as f64) as f32);
    optf!(offmax, 999999.0);
    optf!(fmax, (0.25 / dt as f64) as f32);
    opti!(ntr, 100000);
    opti!(mtr, 100);
    checkpars();
    let surfaces = decode_surfaces(sgn);
    let srf = makesurf((0.025 * dxso as f64) as f32, &surfaces);

    fprintf!(&jpfp, "\n");
    fprintf!(&jpfp, " Datuming parameters\n");
    fprintf!(&jpfp, " ================\n");
    fprintf!(&jpfp, " infile=%s \n", infile);
    fprintf!(&jpfp, " outfile=%s \n", outfile);
    fprintf!(&jpfp, " ttfile=%s \n", ttfile);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzt=%d fzt=%g dzt=%g\n", nzt, fzt, dzt);
    fprintf!(&jpfp, " nxt=%d fxt=%g dxt=%g\n", nxt, fxt, dxt);
    fprintf!(&jpfp, " ns=%d fs=%g ds=%g\n", ns, fs, ds);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nxi=%d fxi=%g dxi=%g sgn=%d\n", nxi, fxi, dxi, sgn);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzo=%d fzo=%g dzo=%g\n", nzo, fzo, dzo);
    fprintf!(&jpfp, " nxso=%d fxso=%g dxso=%g\n", nxso, fxso, dxso);
    fprintf!(&jpfp, " nxgo=%d fxgo=%g dxgo=%g\n", nxgo, fxgo, dxgo);
    fprintf!(&jpfp, " \n");

    let mut szif = Grid::new(2, nxi);
    let mut sz = vec![0.0; ns as usize];
    let mut nangl = vec![0.0; nxi as usize];
    zcoor_surfaces(
        fxi, dxi, nxi, fs, ds, ns, &srf, &mut szif, &mut sz, &mut nangl,
    );
    let mut rmax = cmax(es - fxt, ext - fs);
    rmax = cmin(rmax as f64, 0.5 * offmax as f64 + aperx as f64) as f32;
    let nr = 2 + (rmax / dxgo) as i32;
    let mut tb = vec![Grid::new(nzt, nr); ns as usize];
    let mut pb = tb.clone();
    let mut sigb = tb.clone();
    let mut angb = tb.clone();
    for is in 0..ns as usize {
        let v00 = v0 + sz[is] * dvz;
        timeb(
            nr,
            nzt,
            dxso,
            dzt,
            fzt,
            sz[is],
            dvz,
            v00,
            &mut tb[is],
            &mut pb[is],
            &mut sigb[is],
            &mut angb[is],
        );
    }
    fprintf!(&jpfp, " nt=%d ft=%g dt=%g \n", nt, ft, dt);
    fprintf!(&jpfp, " fmax=%g\n", fmax);
    fprintf!(&jpfp, " v0=%g dvz=%g \n", v0, dvz);
    fprintf!(
        &jpfp,
        " aperx=%g offmax=%g angmax=%g\n",
        aperx,
        offmax,
        angmax
    );
    fprintf!(&jpfp, " ntr=%d mtr=%d \n", ntr, mtr);
    fprintf!(&jpfp, " ================\n");
    fflush(&jpfp);

    let nxgt = ((exgo - fxgo) / dxgo + 1.0) as i32;
    let mut ng = vec![0i32; nxgt as usize];
    let mut ttab = vec![Grid::new(nzt, nxt); ns as usize];
    let mut tsum = Grid::new(nzt, nxt);
    let mut tbs = Grid::new(nzt, nr);
    let mut pbs = tbs.clone();
    let mut sigbs = tbs.clone();
    let mut angbs = tbs.clone();
    let ndat = if nxso < nxgo { nxso } else { nxgo };
    // ealloc3float keeps both the samples and row pointers contiguous.
    // When nxso > nxgo, dat2d can advance into the next receiver's rows;
    // retain that layout instead of allocating each gather separately.
    let mut datg = vec![0.0f32; nt as usize * ndat as usize * nxgt as usize];
    fprintf!(&jpfp, " input traveltime tables \n");
    for is in 0..ns as usize {
        let nseek = nxt as i64 * nzt as i64 * is as i64;
        efseeko(&ttfp, nseek * 4, SEEK_SET);
        efread(&mut ttab[is].data, &ttfp);
        let s = fs + is as f32 * ds;
        resit(nxt, fxt, dxt, nzt, nr, dxgo, &tb[is], &mut ttab[is], s);
    }
    fprintf!(&jpfp, " start datuming sources... \n");
    fprintf!(&jpfp, " \n");
    fflush(&jpfp);
    let mut mtmax =
        ((2.0f32 * dxso) as f64 * (angmax as f64 * PI / 180.0).sin() / (v0 * dt) as f64) as i32;
    if mtmax < 1 {
        mtmax = 1;
    }
    let mut jtr: i32 = 1;
    let mut ktr = 0;
    let mut fxin = fxso;
    let mut s = 0.0;
    let hdrfp = etmpfile();
    let mut filter = Filter::default();
    loop {
        let sx = tr.sx as f32;
        let gx = tr.gx as f32;
        if gx != s {
            fxin = sx;
            s = gx;
        }
        let io = ((gx - fxgo) / dxgo) as i32;
        if cmin(sx, gx) >= fs && cmax(sx, gx) <= es && cmax(gx - sx, sx - gx) <= offmax {
            ng[io as usize] += 1;
            efwrite(&[tr.hdr], &hdrfp);
            let a = (sx - fs) / ds;
            let mut is = a as i32;
            if is == ns - 1 {
                is = ns - 2;
            }
            let mut res = a - is as f32;
            if res as f64 <= 0.01 {
                res = 0.0;
            }
            if res as f64 >= 0.99 {
                res = 1.0;
            }
            let is = is as usize;
            sum2(
                nxt,
                nzt,
                1.0 - res,
                res,
                &ttab[is],
                &ttab[is + 1],
                &mut tsum,
            );
            sum2(nr, nzt, 1.0 - res, res, &pb[is], &pb[is + 1], &mut pbs);
            sum2(nr, nzt, 1.0 - res, res, &tb[is], &tb[is + 1], &mut tbs);
            sum2(
                nr,
                nzt,
                1.0 - res,
                res,
                &sigb[is],
                &sigb[is + 1],
                &mut sigbs,
            );
            sum2(
                nr,
                nzt,
                1.0 - res,
                res,
                &angb[is],
                &angb[is + 1],
                &mut angbs,
            );
            let a = (sx - fxi) / dxi;
            let mut is = a as i32;
            if is == nxi - 1 {
                is = nxi - 2;
            }
            let mut res = a - is as f32;
            if res as f64 <= 0.01 {
                res = 0.0;
            }
            if res as f64 >= 0.99 {
                res = 1.0;
            }
            let nangls = (1.0 - res) * nangl[is as usize] + res * nangl[is as usize + 1];
            // For sgn=-1 the C interpolation can read trace[-1]. In segy
            // this is the final four native bytes of the input header.
            let header = oxwave_core::cio::as_bytes(std::slice::from_ref(&tr.hdr));
            let trace_left = f32::from_ne_bytes(header[HDRBYTES - 4..HDRBYTES].try_into().unwrap());
            dat2d(
                &mut tr.data,
                nt,
                ft,
                dt,
                sx,
                gx,
                &mut datg[io as usize * ndat as usize * nt as usize..],
                aperx,
                nxso,
                fxin,
                dxso,
                nzo as f32,
                fzo,
                dzo,
                mtmax,
                sx,
                fmax,
                nxi,
                fxi,
                dxi,
                angmax,
                &tbs,
                &pbs,
                &angbs,
                nr,
                &tsum,
                nzt,
                fzt,
                dzt,
                nxt,
                fxt,
                dxt,
                antiali,
                sgn,
                &szif,
                nangls,
                &sigbs,
                &mut filter,
                trace_left,
            );
            ktr += 1;
            // arm64 C division by zero returns zero, including this remainder.
            if (jtr - 1).checked_rem(mtr).unwrap_or(0) == 0 {
                fprintf!(&jpfp, " Datumed source %d\n", jtr);
                fflush(&jpfp);
            }
        }
        jtr += 1;
        // C reads the next trace before checking its trace limit.
        if fgettr(&infp, &mut tr) == 0 || jtr > ntr {
            break;
        }
    }
    fprintf!(&jpfp, " Datumed %d sources in total\n", ktr);
    erewind(&hdrfp);
    let scal = (1.0 / (2.0 * PI).sqrt() * dxso as f64) as f32;
    for io in 0..nxgt {
        for ixo in 0..ng[io as usize] {
            efread(std::slice::from_mut(&mut tr.hdr), &hdrfp);
            for it in 0..nt {
                tr.data[it as usize] = datg[((io * ndat + ixo) * nt + it) as usize];
            }
            for it in 0..nt {
                tr.data[it as usize] *= scal;
            }
            fputtr(&outfp, &mut tr);
        }
    }
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " output done\n");
    fflush(&jpfp);
    efclose(&jpfp);
    efclose(&outfp);
    efclose(&hdrfp);
    exit(CWP_Exit());
}

struct SurfaceSegment {
    x: f32,
    z: f32,
    s: f32,
}
struct Surface {
    ss: Vec<SurfaceSegment>,
}

fn decode_surfaces(sgn: i32) -> Vec<(Vec<f32>, Vec<f32>)> {
    let mut nr = countparname("surf");
    if nr == 0 {
        nr = 1;
    }
    let mut surfaces = Vec::new();
    for ir in 0..nr {
        let mut s = String::new();
        if getnparstring(ir + 1, "surf", &mut s) == 0 {
            s = "0,0;99999,0".into();
        }
        // C counts with strtok, then walks the modified string using
        // strlen+1. Preserve its treatment of consecutive delimiters too.
        let mut bytes = s.into_bytes();
        bytes.push(0);
        let mut pos = 0;
        let mut first = None;
        let mut ntokens = 0;
        while bytes[pos] != 0 {
            while bytes[pos] == b',' || bytes[pos] == b';' {
                pos += 1;
            }
            if bytes[pos] == 0 {
                break;
            }
            if first.is_none() {
                first = Some(pos);
            }
            ntokens += 1;
            while bytes[pos] != 0 && bytes[pos] != b',' && bytes[pos] != b';' {
                pos += 1;
            }
            if bytes[pos] != 0 {
                bytes[pos] = 0;
                pos += 1;
            }
        }
        if ntokens % 2 != 0 || ntokens < 4 {
            err!("Surface number %d specified incorrectly!\n", ir + 1);
        }
        let mut x = Vec::new();
        let mut z = Vec::new();
        let mut pos = first.unwrap();
        for _ in 0..ntokens / 2 {
            let len = bytes[pos..].iter().position(|&b| b == 0).unwrap();
            x.push(atof(&String::from_utf8_lossy(&bytes[pos..pos + len])) as f32);
            pos += len + 1;
            let len = bytes[pos..].iter().position(|&b| b == 0).unwrap();
            z.push((atof(&String::from_utf8_lossy(&bytes[pos..pos + len])) * sgn as f64) as f32);
            pos += len + 1;
        }
        surfaces.push((x, z));
    }
    surfaces
}

fn cub(derivative: i32, u: &[f32], coefficients: &[[f32; 4]], point: f32) -> f32 {
    let mut value = [0.0];
    intcub(derivative, u, coefficients, &[point], &mut value);
    value[0]
}

fn makesurf(dsmax: f32, surfaces: &[(Vec<f32>, Vec<f32>)]) -> Vec<Surface> {
    let mut rr = Vec::new();
    for (xu, zu) in surfaces {
        let nu = xu.len();
        let mut u = vec![0.0; nu];
        for iu in 0..nu {
            u[iu] = iu as f32;
        }
        let mut xud = vec![[0.0; 4]; nu];
        csplin(&u, xu, &mut xud);
        let mut zud = vec![[0.0; 4]; nu];
        csplin(&u, zu, &mut zud);
        let nuu = 20 * nu;
        let duu = (u[nu - 1] - u[0]) / (nuu - 1) as f32;
        let mut s = vec![0.0; nuu];
        let mut xlast = xu[0];
        let mut zlast = zu[0];
        let mut uu = duu;
        for iuu in 1..nuu {
            let x = cub(0, &u, &xud, uu);
            let z = cub(0, &u, &zud, uu);
            let dx = x - xlast;
            let dz = z - zlast;
            s[iuu] = (s[iuu - 1] as f64 + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
            xlast = x;
            zlast = z;
            uu += duu;
        }
        let ns = (1.0f32 + s[nuu - 1] / dsmax) as i32;
        let ds = s[nuu - 1] / ns as f32;
        let fs = (0.5 * ds as f64) as f32;
        let mut us = vec![0.0; ns as usize];
        yxtoxy(duu, 0.0, &s, ds, fs, 0.0, (nu - 1) as f32, &mut us);
        let mut ss = Vec::with_capacity(ns as usize);
        for is in 0..ns as usize {
            let rsx = cub(0, &u, &xud, us[is]);
            let rsz = cub(0, &u, &zud, us[is]);
            let rsxd = cub(1, &u, &xud, us[is]);
            let rszd = cub(1, &u, &zud, us[is]);
            // The C also stores c=rsxd/sqrt(...); that field is never read.
            let s = (-rszd as f64 / ((rsxd * rsxd + rszd * rszd) as f64).sqrt()) as f32;
            ss.push(SurfaceSegment { x: rsx, z: rsz, s });
        }
        rr.push(Surface { ss });
    }
    rr
}

#[allow(clippy::too_many_arguments)]
fn zcoor_surfaces(
    fx: f32,
    dx: f32,
    nx: i32,
    fxs: f32,
    dxs: f32,
    nxs: i32,
    srf: &[Surface],
    szif: &mut Grid,
    sz: &mut [f32],
    nangl: &mut [f32],
) {
    // decodeSurfaces supplies only one surface when surf is omitted. C
    // nevertheless reads srf[1]; depending on adjacent allocator contents,
    // it either faults or skips the missing surface's loops, leaving its
    // freshly allocated depths at zero. Preserve that successful path
    // explicitly, without an invalid read or an unconditional signal.
    for (is, surface) in srf.iter().take(2).enumerate() {
        let ss = &surface.ss;
        let nss = ss.len();
        let mut ixi = 0;
        for ix in 0..nx {
            let x = fx + ix as f32 * dx;
            for ik in ixi..nss - 1 {
                if ss[ik].x <= x && ss[ik + 1].x >= x {
                    let az = ss[ik].z - ss[ik + 1].z;
                    let ax0 = ss[ik + 1].x - x;
                    let ax = ss[ik + 1].x - ss[ik].x;
                    szif[(ix, is as i32)] = ax0 * az / ax + ss[ik + 1].z;
                    if is == 0 {
                        let az = ss[ik].s - ss[ik + 1].s;
                        let temp = ax0 * az / ax + ss[ik + 1].s;
                        nangl[ix as usize] = (temp as f64).asin() as f32;
                    }
                    ixi = ik;
                }
            }
        }
        if is == 0 {
            let mut ixi = 0;
            for ix in 0..nxs {
                let x = fxs + ix as f32 * dxs;
                for ik in ixi..nss - 1 {
                    if ss[ik].x <= x && ss[ik + 1].x >= x {
                        let az = ss[ik].z - ss[ik + 1].z;
                        let ax0 = ss[ik + 1].x - x;
                        let ax = ss[ik + 1].x - ss[ik].x;
                        sz[ix as usize] = ax0 * az / ax + ss[ik + 1].z;
                        ixi = ik;
                    }
                }
            }
        }
    }
    for is in 0..2 {
        szif[(0, is)] = 2.0 * szif[(1, is)] - szif[(2, is)];
        szif[(nx - 1, is)] = 2.0 * szif[(nx - 2, is)] - szif[(nx - 3, is)];
    }
    // C lines 751/753 read beyond sz when nxs=2. Retaining the depths
    // already interpolated above matches both original and bounds-fixed
    // C for the round 1 case (byte-identical output).
    if nxs >= 3 {
        sz[0] = 2.0 * sz[1] - sz[2];
    }
    nangl[0] = 2.0 * nangl[1] - nangl[2];
    if nxs >= 3 {
        sz[nxs as usize - 1] = 2.0 * sz[nxs as usize - 2] - sz[nxs as usize - 3];
    }
    nangl[nx as usize - 1] = 2.0 * nangl[nx as usize - 2] - nangl[nx as usize - 3];
}

#[allow(clippy::too_many_arguments)]
fn resit(nx: i32, fx: f32, dx: f32, nz: i32, nr: i32, dr: f32, tb: &Grid, t: &mut Grid, x0: f32) {
    for ix in 0..nx {
        let xi = fx + ix as f32 * dx - x0;
        let ar = xi.abs() / dr;
        let mut jr = ar as i32;
        let sr = ar - jr as f32;
        let sr0 = (1.0 - sr as f64) as f32;
        if jr > nr - 2 {
            jr = nr - 2;
        }
        for iz in 0..nz {
            t[(ix, iz)] -= sr0 * tb[(jr, iz)] + sr * tb[(jr + 1, iz)];
        }
    }
}

fn sum2(nx: i32, nz: i32, a1: f32, a2: f32, t1: &Grid, t2: &Grid, t: &mut Grid) {
    for ix in 0..nx {
        for iz in 0..nz {
            t[(ix, iz)] = a1 * t1[(ix, iz)] + a2 * t2[(ix, iz)];
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn timeb(
    nr: i32,
    nz: i32,
    dr: f32,
    dz: f32,
    fz: f32,
    z0: f32,
    a: f32,
    v0: f32,
    t: &mut Grid,
    p: &mut Grid,
    sig: &mut Grid,
    ang: &mut Grid,
) {
    if a == 0.0 {
        let mut r = 0.0f32;
        for ir in 0..nr {
            let mut z = fz - z0;
            for iz in 0..nz {
                let mut rou = ((r * r + z * z) as f64).sqrt() as f32;
                if rou < dz {
                    rou = dz;
                }
                t[(ir, iz)] = rou / v0;
                p[(ir, iz)] = r / (rou * v0);
                sig[(ir, iz)] = v0 * rou;
                ang[(ir, iz)] = ((r / rou) as f64).asin() as f32;
                z += dz;
            }
            r += dr;
        }
    } else {
        let oa = (1.0 / a as f64) as f32;
        let zc = v0 * oa;
        let mut r = 0.0f32;
        for ir in 0..nr {
            let mut z = fz + zc - z0;
            for iz in 0..nz {
                let _rou = ((r * r + z * z) as f64).sqrt() as f32;
                let v = v0 + a * (z - zc + z0);
                if ir == 0 {
                    t[(ir, iz)] = (((v / v0) as f64).ln() * oa as f64) as f32;
                    p[(ir, iz)] = 0.0;
                    ang[(ir, iz)] = 0.0;
                    sig[(ir, iz)] = (0.5
                        * (z as f64 + 0.1 * dz as f64 - zc as f64 + z0 as f64)
                        * (v0 + v) as f64) as f32;
                } else {
                    let rc = ((r * r + z * z - zc * zc) as f64 / (2.0 * r as f64)) as f32;
                    let rou = ((zc * zc + rc * rc) as f64).sqrt() as f32;
                    t[(ir, iz)] =
                        (((v * (rou + rc) / (v0 * (rou + rc - r))) as f64).ln() * oa as f64) as f32;
                    p[(ir, iz)] =
                        (((rou * rou - rc * rc) as f64).sqrt() / (rou * v0) as f64) as f32;
                    let mut temp = v0 * p[(ir, iz)];
                    if temp > 1.0 {
                        temp = 1.0;
                    }
                    ang[(ir, iz)] = (temp as f64).asin() as f32;
                    sig[(ir, iz)] = a * rou * r;
                }
                z += dz;
            }
            r += dr;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn dat2d(
    trace: &mut [f32],
    nt: i32,
    ft: f32,
    dt: f32,
    sx: f32,
    gx: f32,
    dat: &mut [f32],
    aperx: f32,
    nx: i32,
    fx: f32,
    dx: f32,
    nz: f32,
    fz: f32,
    dz: f32,
    mtmax: i32,
    xm: f32,
    fmax: f32,
    nxi: i32,
    fxi: f32,
    dxi: f32,
    angmax: f32,
    tb: &Grid,
    pb: &Grid,
    angb: &Grid,
    nr: i32,
    tsum: &Grid,
    nzt: i32,
    fzt: f32,
    dzt: f32,
    nxt: i32,
    fxt: f32,
    dxt: f32,
    antiali: i32,
    sgn: i32,
    szif: &Grid,
    nangl: f32,
    sigb: &Grid,
    filter: &mut Filter,
    trace_left: f32,
) {
    let mut tmt = Grid::new(nzt, nxt);
    let mut ampt = Grid::new(nzt, nxt);
    let mut ampti = Grid::new(nzt, nxt);
    let mut tm = vec![0.0f32; nzt as usize];
    let mut tzt = vec![0.0f32; nzt as usize];
    let mut amp = vec![0.0f32; nzt as usize];
    let mut ampi = vec![0.0f32; nzt as usize];
    let mut zpt = vec![0.0f32; nxt as usize];
    let mut trf = vec![0.0f32; (nt + 2 * mtmax) as usize];
    let odt = (1.0 / dt as f64) as f32;
    let z0 = (((fz - fzt) / dzt) as f64 + 0.0 * sx as f64 + 0.0 * gx as f64) as f32;
    let rdz = dz / dzt;
    let pmin = (1.0 / (2.0 * dx as f64 * fmax as f64)) as f32;
    let fxf = fxi + (nxi - 2) as f32 * dxi;
    filter.filt(trace, nt, dt, fmax, mtmax, &mut trf, sgn);
    let rxz = if angmax == 90.0 {
        0.0
    } else {
        (1.0 / (angmax as f64 * PI / 180.0).tan()) as f32
    };
    let mut nxtf = ((xm - aperx - fxt) / dxt) as i32;
    if nxtf < 0 {
        nxtf = 0;
    }
    let mut nxte = ((xm + aperx - fxt) / dxt + 1.0) as i32;
    if nxte >= nxt {
        nxte = nxt - 1;
    }
    // Compute amplitudes and filter length.
    for ix in nxtf..=nxte {
        let x = fxt + ix as f32 * dxt;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = (((dis - dxt) * rxz - fzt) / dzt - 1.0) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let ar = if xm >= x {
            (xm - x) / dx
        } else {
            (x - xm) / dx
        };
        let mut jr = ar as i32;
        if jr > nr - 2 {
            jr = nr - 2;
        }
        let sr = ar - jr as f32;
        let sr0 = (1.0 - sr as f64) as f32;
        let sig = if xm - x < 0.0 { 1.0f32 } else { -1.0f32 };
        zpt[ix as usize] = fzt + (nzt - 1) as f32 * dzt;
        for iz in izt0..nzt {
            let sigma = sr0 * sigb[(jr, iz)] + sr * sigb[(jr + 1, iz)];
            let ang = sr0 * angb[(jr, iz)] + sr * angb[(jr + 1, iz)];
            let ang = cabs(sig * ang - nangl);
            let mut ampd =
                ((ang as f64).cos() / ((nangl as f64).cos() * (sigma as f64).sqrt())) as f32;
            // Peels (1988): filter 90 degrees to the operator.
            if cabs(ang) as f64 >= PI / 2.0 {
                ampd = 0.0;
            }
            if ampd < 0.0 {
                ampd = -ampd;
            }
            ampt[(ix, iz)] = ampd;
            let mut pd = sr0 * pb[(jr, iz)] + sr * pb[(jr + 1, iz)];
            if pd < 0.0 {
                pd = -pd;
            }
            let mut temp = pd * dx * odt;
            if temp < 1.0 {
                temp = 1.0;
            }
            if temp > mtmax as f32 {
                temp = mtmax as f32;
            }
            ampti[(ix, iz)] = ampd / (temp * temp);
            tmt[(ix, iz)] = temp;
            if pd < pmin && zpt[ix as usize] as f64 > fzt as f64 + (nzt as f64 - 1.1) * dzt as f64 {
                zpt[ix as usize] = fzt + iz as f32 * dzt;
            }
        }
    }
    let mut nxf = (((xm - aperx - fx) / dx) as f64 + 0.5) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    let mut nxe = if xm + aperx >= fxf {
        (((fxf - fx) / dx) as f64 + 0.5) as i32
    } else {
        (((xm + aperx - fx) / dx) as f64 + 0.5) as i32
    };
    if nxe >= nx {
        nxe = nx - 1;
    }
    let mut am = (fx - fxi) / dxi;
    let mut mr = am as i32;
    am -= mr as f32;
    if am as f64 <= 0.01 {
        am = 0.0;
    }
    if am as f64 >= 0.99 {
        am = 1.0;
    }
    let am0 = (1.0 - am as f64) as f32;
    if mr < 0 {
        mr = 0;
    }
    if mr + nxe >= nxi - 1 {
        err!("Topography definition is out of range!\n");
    }
    // Interpolate amplitudes and filter length along lateral.
    for ix in nxf..=nxe {
        let x = fx + ix as f32 * dx;
        let dis = if xm >= x { xm - x } else { x - xm };
        let mut izt0 = ((dis * rxz - fzt) / dzt) as i32;
        if izt0 < 0 {
            izt0 = 0;
        }
        if izt0 >= nzt {
            izt0 = nzt - 1;
        }
        let mut iz0 = ((dis * rxz - fz) / dz) as i32;
        if iz0 < 0 {
            iz0 = 0;
        }
        if iz0 as f32 >= nz {
            iz0 = (nz - 1.0) as i32;
        }
        let mut ax = (x - fxt) / dxt;
        let mut jx = ax as i32;
        ax -= jx as f32;
        if ax as f64 <= 0.01 {
            ax = 0.0;
        }
        if ax as f64 >= 0.99 {
            ax = 1.0;
        }
        let ax0 = (1.0 - ax as f64) as f32;
        if jx > nxte - 1 {
            jx = nxte - 1;
        }
        if jx < nxtf {
            jx = nxtf;
        }
        let ar = if xm >= x {
            (xm - x) / dx
        } else {
            (x - xm) / dx
        };
        let mut jr = ar as i32;
        if jr > nr - 2 {
            jr = nr - 2;
        }
        let sr = ar - jr as f32;
        let sr0 = (1.0 - sr as f64) as f32;
        for iz in izt0..nzt {
            tzt[iz as usize] = ax0 * tsum[(jx, iz)]
                + ax * tsum[(jx + 1, iz)]
                + sr0 * tb[(jr, iz)]
                + sr * tb[(jr + 1, iz)];
            amp[iz as usize] = ax0 * ampt[(jx, iz)] + ax * ampt[(jx + 1, iz)];
            ampi[iz as usize] = ax0 * ampti[(jx, iz)] + ax * ampti[(jx + 1, iz)];
            tm[iz as usize] = ax0 * tmt[(jx, iz)] + ax * tmt[(jx + 1, iz)];
        }
        let mut nzp =
            (((ax0 * zpt[jx as usize] + ax * zpt[jx as usize + 1] - fz) / dz) as f64 + 0.5) as i32;
        if nzp < iz0 {
            nzp = iz0;
        }
        if nzp as f32 > nz {
            nzp = nz as i32;
        }
        let iz = ((cabs(am0 * szif[(ix + mr, 1)] + am * szif[(ix + mr + 1, 1)]) - fz) / dz) as i32;
        if iz as f32 >= nz {
            err!("Datuming surface is out of output range!\n");
        }
        let az = z0 + iz as f32 * rdz;
        let mut jz = az as i32;
        if jz >= nzt - 1 {
            jz = if nzt > 1 { nzt - 2 } else { 0 };
        }
        // A singleton depth table uses its sole sample at both endpoints.
        // C lines 1026-1052 instead read index -1 with weight zero; this
        // bounds fix leaves the round 1 case's output byte-identical.
        let jz1 = cmin(jz + 1, nzt - 1) as usize;
        let sz = az - jz as f32;
        let sz0 = (1.0 - sz as f64) as f32;
        let td = sz0 * tzt[jz as usize] + sz * tzt[jz1];
        let mut at = (sgn as f32 * td - ft) * odt;
        if iz < nzp && antiali != 0 {
            at += mtmax as f32;
            let jt = at as i32;
            let ampd = sz0 * ampi[jz as usize] + sz * ampi[jz1];
            let mt = (0.5 + (sz0 * tm[jz as usize]) as f64 + (sz * tm[jz1]) as f64) as i32;
            let res = cabs(at - jt as f32);
            let res0 = (1.0 - res as f64) as f32;
            for it in 0..nt {
                if it + jt >= mtmax && jt < nt - it + mtmax - sgn {
                    // With upward continuation C can read trf[-1]. The
                    // preceding allocator word is zero in the reference
                    // build; represent that one padding sample explicitly.
                    let left = if it + jt - mt + sgn == -1 {
                        0.0
                    } else {
                        trf[(it + jt - mt + sgn) as usize]
                    };
                    let temp = (res0 as f64
                        * (-(trf[(it + jt - mt) as usize] as f64)
                            + 2.0 * trf[(it + jt) as usize] as f64
                            - trf[(it + jt + mt) as usize] as f64)
                        + res as f64
                            * (-(left as f64) + 2.0 * trf[(it + jt + sgn) as usize] as f64
                                - trf[(it + jt + mt + sgn) as usize] as f64))
                        * ampd as f64;
                    dat[(ix * nt + it) as usize] += temp as f32;
                }
            }
        } else {
            let jt = at as i32;
            let ampd = sz0 * amp[jz as usize] + sz * amp[jz1];
            let res = cabs(at - jt as f32);
            let res0 = (1.0 - res as f64) as f32;
            for it in 0..nt {
                if it + jt >= 0 && jt < nt - it - sgn {
                    let next = if it + jt + sgn == -1 {
                        trace_left
                    } else {
                        trace[(it + jt + sgn) as usize]
                    };
                    let temp = (res0 * trace[(it + jt) as usize] + res * next) * ampd;
                    dat[(ix * nt + it) as usize] += temp;
                }
            }
        }
    }
}

#[derive(Default)]
struct Filter {
    nfft: i32,
    itaper: i32,
    nw: i32,
    nwf: i32,
    taper: Vec<f32>,
    amp: Vec<f32>,
    ampi: Vec<f32>,
    dw: f32,
}
impl Filter {
    #[allow(clippy::too_many_arguments)]
    fn filt(
        &mut self,
        trace: &mut [f32],
        nt: i32,
        dt: f32,
        mut fmax: f32,
        m: i32,
        trf: &mut [f32],
        sgn: i32,
    ) {
        fmax = (fmax as f64 * (2.0 * PI)) as f32;
        let ftaper = (0.1 * fmax as f64) as f32;
        let const2 = (0.5 * 2.0f64.sqrt()) as f32;
        if self.nfft == 0 {
            self.nfft = npfaro(2 * (nt + 2 * m), 4 * (nt + 2 * m));
            if self.nfft >= SU_NFLTS as i32 || self.nfft >= 720720 {
                err!("Padded nt=%d -- too big", self.nfft);
            }
            self.nw = self.nfft / 2 + 1;
            self.dw = (2.0 * PI / (self.nfft as f32 * dt) as f64) as f32;
            self.itaper = (0.5 + (ftaper / self.dw) as f64) as i32;
            self.taper = vec![0.0; (2 * self.itaper + 1) as usize];
            for iw in -self.itaper..=self.itaper {
                let temp = (iw as f64 / (1.0 + self.itaper as f64)) as f32;
                self.taper[(iw + self.itaper) as usize] =
                    (1.0 - temp) * (1.0 - temp) * (temp + 2.0) / 4.0;
            }
            self.nwf = (0.5 + (fmax / self.dw) as f64) as i32;
            if self.nwf > self.nw - self.itaper - 1 {
                self.nwf = self.nw - self.itaper - 1;
            }
            self.amp = vec![0.0; (self.nwf + self.itaper + 1) as usize];
            self.ampi = self.amp.clone();
            for iw in 1..=self.nwf + self.itaper {
                self.amp[iw as usize] =
                    (((self.dw * iw as f32) as f64).sqrt() / self.nfft as f64) as f32;
                self.ampi[iw as usize] =
                    (0.5 / (1.0 - ((iw as f32 * self.dw * dt) as f64).cos())) as f32;
            }
        }
        let mut rt = vec![0.0; self.nfft as usize];
        let mut ct = vec![Complex::default(); self.nw as usize];
        rt[..nt as usize].copy_from_slice(&trace[..nt as usize]);
        pfarc(1, &rt, &mut ct);
        for iw in self.nwf - self.itaper..=self.nwf + self.itaper {
            let itemp = iw - (self.nwf - self.itaper);
            ct[iw as usize].r = self.taper[itemp as usize] * ct[iw as usize].r;
            ct[iw as usize].i = self.taper[itemp as usize] * ct[iw as usize].i;
        }
        for iw in self.nwf + self.itaper + 1..self.nw {
            ct[iw as usize].r = 0.0;
            ct[iw as usize].i = 0.0;
        }
        for iw in 0..=self.nwf + self.itaper {
            let c = &mut ct[iw as usize];
            let temp = (c.r - sgn as f32 * c.i) * self.amp[iw as usize] * const2;
            c.i = (sgn as f32 * c.r + c.i) * self.amp[iw as usize] * const2;
            c.r = temp;
        }
        pfacr(-1, &ct, &mut rt);
        for it in 0..nt {
            trace[it as usize] = rt[it as usize];
        }
        for iw in 0..=self.nwf + self.itaper {
            ct[iw as usize].i *= self.ampi[iw as usize];
            ct[iw as usize].r *= self.ampi[iw as usize];
        }
        pfacr(-1, &ct, &mut rt);
        for it in 0..m {
            trf[it as usize] = rt[(self.nfft - m + it) as usize];
        }
        for it in 0..nt + m {
            trf[(it + m) as usize] = rt[it as usize];
        }
    }
}
