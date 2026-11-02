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

//! SUKDMDCR - 2.5D datuming of receivers for common-source prestack data.
//! Port of su/main/datuming/sukdmdcr.c (Revision 1.3).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: Steven D. Sheaffer (CWP), 11/97.
//! Adapts Bleistein & Jaramillo's 2.5D data mapping formula for receiver
//! extrapolation. Reference: Sheaffer, S., 1999, "2.5D Downward Continuation
//! of the Seismic Wavefield Using Kirchhoff Data Mapping," M.Sc. Thesis,
//! Dept. of Mathematical & Computer Sciences, Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{abs, max, min, PI};
use oxwave_pfa_fft::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    "  SUKDMDCR - 2.5D datuming of receivers for prestack, common source    ",
    "            data using constant-background data mapping formula.       ",
    "            (See selfdoc for specific survey requirements.)            ",
    " \t\t\t\t\t\t\t\t\t",
    "    sukdmdcr  infile=  outfile=  [parameters] \t         \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required file parameters:\t\t\t\t\t\t",
    " infile=stdin\t\tfile for input seismic traces\t\t\t",
    " outfile=stdout\tfile for output          \t\t\t",
    " ttfile\t\tfile for input traveltime tables\t\t",
    "                                                                       ",
    " Required parameters describing the traveltime tables:\t\t        ",
    " fzt \t\t\tfirst depth sample in traveltime table\t\t",
    " nzt \t\t\tnumber of depth samples in traveltime table\t",
    " dzt\t\t\tdepth interval in traveltime table\t\t",
    " fxt\t\t\tfirst lateral sample in traveltime table\t",
    " nxt\t\t\tnumber of lateral samples in traveltime table\t",
    " dxt\t\t\tlateral interval in traveltime table\t\t",
    " fs \t\t\tx-coordinate of first source in table\t\t",
    " ns \t\t\tnumber of sources in table\t\t\t",
    " ds \t\t\tx-coordinate increment of sources in table\t",
    "\t\t\t\t\t\t\t\t\t",
    " Parameters describing the input data:                                 ",
    " nxso                  number of shots                                 ",
    " dxso                  shot interval                                   ",
    " fxso=0                x-coordinate of first shot                      ",
    " nxgo                  number of receiver offsets per shot             ",
    " dxgo                  receiver offset interval                        ",
    " fxgo=0                first receiver offset                           ",
    " dt= or from header (dt)       time sampling interval of input data    ",
    " ft= or from header (ft)       first time sample of input data         ",
    "                                                                       ",
    " Parameters describing the domain of the problem:             \t\t",
    " dzo=0.2*dzt\t\tvertical spacing in surface determination       ",
    " offmax=99999\t\tmaximum absolute offset allowed          \t",
    "                                                                       ",
    " Parameters describing the recording and datuming surfaces:            ",
    " recsurf=0             recording surface (horizontal=0, topographic=1) ",
    " zrec                  defines recording surface when recsurf=0        ",
    " recfile=              defines recording surface when recsurf=1        ",
    " datsurf=0             datuming surface (horizontal=0, irregular=1)    ",
    " zdat                  defines datuming surface when datsurf=0         ",
    " datfile=              defines datuming surface when datsurf=1         ",
    "                                                                       ",
    " Optional parameters describing the extrapolation:                     ",
    " aperx=nxt*dxt/2  \tlateral half-aperture \t\t\t\t",
    " v0=1500(m/s)\t\treference wavespeed               \t\t",
    " freq=50               dominant frequency in data, used to determine   ",
    "                       the minimum distance below the datum that       ",
    "                       the stationary phase calculation is valid.      ",
    " scale=1.0             user defined scale factor for output            ",
    " jpfile=stderr\t\tjob print file name \t\t\t\t",
    " mtr=100  \t\tprint verbal information at every mtr traces\t",
    " ntr=100000\t\tmaximum number of input traces to be datumed\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    // Explicitly define downward continuation.
    let sgn = 1;
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut datain = String::from("stdin");
    let infp = if getparstring("datain", &mut datain) == 0 {
        stdin()
    } else {
        fopen(&datain, "r").unwrap_or_else(|| err!("cannot open datain=%s\n", datain))
    };
    let mut dataout = String::from("stdout");
    let outfp = if getparstring("dataout", &mut dataout) == 0 {
        stdout()
    } else {
        unchecked_file(fopen(&dataout, "w"))
    };
    fseeko(&infp, 0, SEEK_CUR);
    fseeko(&outfp, 0, SEEK_CUR);
    let mut ttfile = String::new();
    if getparstring("ttfile", &mut ttfile) == 0 {
        err!("must specify ttfile!\n");
    }
    let ttfp = fopen(&ttfile, "r").unwrap_or_else(|| err!("cannot open ttfile=%s\n", ttfile));
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };
    if fgettr(&infp, &mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if (dt as f64) < 0.0000001 {
        err!("dt must be positive!\n");
    }
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = (tr.delrt as f32 as f64 / 1000.0) as f32;
    }
    let mut nxt = 0i32;
    if getparint("nxt", &mut nxt) == 0 {
        err!("must specify nxt!\n");
    }
    let mut fxt = 0.0f32;
    if getparfloat("fxt", &mut fxt) == 0 {
        err!("must specify fxt!\n");
    }
    let mut dxt = 0.0f32;
    if getparfloat("dxt", &mut dxt) == 0 {
        err!("must specify dxt!\n");
    }
    let mut nzt = 0i32;
    if getparint("nzt", &mut nzt) == 0 {
        err!("must specify nzt!\n");
    }
    let mut fzt = 0.0f32;
    if getparfloat("fzt", &mut fzt) == 0 {
        err!("must specify fzt!\n");
    }
    let mut dzt = 0.0f32;
    if getparfloat("dzt", &mut dzt) == 0 {
        err!("must specify dzt!\n");
    }
    let mut ns = 0i32;
    if getparint("ns", &mut ns) == 0 {
        err!("must specify ns!\n");
    }
    let mut fs = 0.0f32;
    if getparfloat("fs", &mut fs) == 0 {
        err!("must specify fs!\n");
    }
    let mut ds = 0.0f32;
    if getparfloat("ds", &mut ds) == 0 {
        err!("must specify ds!\n");
    }
    let ext = fxt + (nxt - 1) as f32 * dxt;
    let ezt = fzt + (nzt - 1) as f32 * dzt;
    let es = fs + (ns - 1) as f32 * ds;
    let mut nxso = 0;
    if getparint("nxso", &mut nxso) == 0 {
        err!("must specify nxso!\n");
    }
    let mut fxso = 0.0f32;
    if getparfloat("fxso", &mut fxso) == 0 {
        fxso = 0.0;
    }
    let mut dxso = 0.0f32;
    if getparfloat("dxso", &mut dxso) == 0 {
        err!("must specify dxso!\n");
    }
    let mut nxgo = 0;
    if getparint("nxgo", &mut nxgo) == 0 {
        err!("must specify nxgo!\n");
    }
    let mut fxgo = 0.0f32;
    if getparfloat("fxgo", &mut fxgo) == 0 {
        fxgo = 0.0;
    }
    let mut dxgo = 0.0f32;
    if getparfloat("dxgo", &mut dxgo) == 0 {
        err!("must specify dxgo!\n");
    }
    if nxso > 1 && dxgo != dxso {
        err!("in this implementation, dxso must = dxgo!\n");
    }
    if dxgo < 0.0 {
        err!("in this implementation, dxgo must be >0!\n");
    }
    if dxso < 0.0 {
        err!("in this implementation, dxso must be >0!\n");
    }
    let dxi = dxgo;
    let fxi = fxgo;
    let nxi = ((fxgo + (nxgo + nxso - 2) as f32 * dxgo) / dxi) as i32 + 1;
    let exso = fxso + (nxso - 1) as f32 * dxso;
    let exgo = fxgo + (nxi - 1) as f32 * dxgo;
    let mut dzo = 0.0f32;
    if getparfloat("dzo", &mut dzo) == 0 {
        dzo = (dzt as f64 * 0.2) as f32;
    }
    let fzo = fzt * sgn as f32;
    let nzo = (1.0f32 + ((nzt - 1) as f32 * dzt) / dzo) as i32;
    dzo *= sgn as f32;
    let ezo = fzo + (nzo - 1) as f32 * dzo;
    if fxt > fxso || fxt > fxgo || ext < exso || ext < exgo || fzt > fzo || ezt < ezo {
        err!("output range is out of traveltime table!\n");
    }
    let mut recsurf = 0;
    if getparint("recsurf", &mut recsurf) == 0 {
        recsurf = 0;
    }
    let mut datsurf = 0;
    if getparint("datsurf", &mut datsurf) == 0 {
        datsurf = 0;
    }
    let mut zrec = 0.0f32;
    let mut zdat = 0.0f32;
    if recsurf == 0 && getparfloat("zrec", &mut zrec) == 0 {
        err!("must specify zrec when recsurf=0!\n");
    }
    if datsurf == 0 && getparfloat("zdat", &mut zdat) == 0 {
        err!("must specify zdat when datsurf=0!\n");
    }
    let mut v0 = 1500.0f32;
    if getparfloat("v0", &mut v0) == 0 {
        v0 = 1500.0;
    }
    let mut aperx = 0.0f32;
    if getparfloat("aperx", &mut aperx) == 0 {
        aperx = (0.5 * nxt as f64 * dxt as f64) as f32;
    }
    let mut offmax = 99999.0f32;
    if getparfloat("offmax", &mut offmax) == 0 {
        offmax = 99999.0;
    }
    let mut freq = 50.0f32;
    if getparfloat("freq", &mut freq) == 0 {
        freq = 50.0;
    }
    if freq <= 0.0 {
        err!("freq must be positive and non-zero!\n");
    }
    if 1.0f32 / (2.0f32 * freq) < dt {
        warn!("freq set too high - you may have numerical singularities !\n");
    }
    let mut scale = 1.0f32;
    if getparfloat("scale", &mut scale) == 0 {
        scale = 1.0;
    }
    let mut ntr = 100000;
    if getparint("ntr", &mut ntr) == 0 {
        ntr = 100000;
    }
    let mut mtr = 100;
    if getparint("mtr", &mut mtr) == 0 {
        mtr = 100;
    }
    let jpfp = unchecked_file(jpfp);
    fprintf!(&jpfp, "\n");
    fprintf!(&jpfp, " Datuming parameters\n");
    fprintf!(&jpfp, " ================\n");
    fprintf!(&jpfp, " datain=%s \n", datain);
    fprintf!(&jpfp, " dataout=%s \n", dataout);
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
    fprintf!(&jpfp, " nt=%d ft=%g dt=%g \n", nt, ft, dt);
    fprintf!(&jpfp, " freq=%g v0=%g\n", freq, v0);
    fprintf!(&jpfp, " aperx=%g offmax=%g \n", aperx, offmax);
    fprintf!(&jpfp, " ntr=%d mtr=%d \n", ntr, mtr);
    fprintf!(&jpfp, " ================\n");
    fflush(&jpfp);

    let mut recfile = String::new();
    let mut datfile = String::new();
    let mut recfp = None;
    let mut datfp = None;
    if recsurf == 1 {
        if getparstring("recfile", &mut recfile) == 0 {
            err!("you chose recsurf=1, so you must specify a recfile!\n");
        }
        recfp = fopen(&recfile, "r");
    }
    if datsurf == 1 {
        if getparstring("datfile", &mut datfile) == 0 {
            err!("you chose datsurf=1, so you must specify a datfile!\n");
        }
        datfp = fopen(&datfile, "r");
    }
    checkpars();
    let mut szif = vec![[0.0f32; 2]; nxi as usize];
    if recsurf == 0 {
        for i in 0..nxi as usize {
            szif[i][0] = zrec;
        }
    } else {
        let recfp = unchecked_file(recfp);
        // C tests fscanf < 0, retaining r1 on a matching failure. A failure
        // before any assignment leaves r1 indeterminate in the original.
        let mut r1 = 0.0f32;
        for i in 0..nxi as usize {
            if let Some(value) = fscanf_f32(&recfp) {
                r1 = value;
            } else if feof(&recfp) || ferror(&recfp) {
                err!("not enough points in recfile!\n");
            }
            szif[i][0] = r1;
        }
        fclose(&recfp);
    }
    if datsurf == 0 {
        for i in 0..nxi as usize {
            szif[i][1] = zdat;
        }
    } else {
        let datfp = unchecked_file(datfp);
        let mut r2 = 0.0f32;
        for i in 0..nxi as usize {
            if let Some(value) = fscanf_f32(&datfp) {
                r2 = value;
            } else if feof(&datfp) || ferror(&datfp) {
                err!("not enough points in datfile!\n");
            }
            szif[i][1] = r2;
        }
        fclose(&datfp);
    }
    // ealloc3float storage is contiguous; each source has nxt*nzt floats.
    let table_size = nxt as usize * nzt as usize;
    let mut ttab = vec![0.0f32; ns as usize * table_size];
    let mut tsum = vec![vec![0.0f32; nzt as usize]; nxt as usize];
    let mut ng = vec![0i32; nxso as usize];
    let gather_size = nxgo as usize * nt as usize;
    let mut dats = vec![0.0f32; nxso as usize * gather_size];
    fprintf!(&jpfp, " input traveltime tables \n");
    for is in 0..ns as usize {
        let nseek = nxt.wrapping_mul(nzt).wrapping_mul(is as i32) as i64;
        fseeko(&ttfp, nseek * 4, SEEK_SET);
        fread(&mut ttab[is * table_size..(is + 1) * table_size], &ttfp);
    }
    fprintf!(&jpfp, " start receiver datuming ... \n");
    fprintf!(&jpfp, " \n");
    fflush(&jpfp);
    let mut jtr = 1i32;
    let mut ktr = 0i32;
    let mut fxin = fxgo;
    let mut s = fxso;
    let hdrfp = etmpfile();
    let mut filter = Filter::default();
    loop {
        let sx = tr.sx as f32;
        let gx = tr.gx as f32;
        if sx != s {
            s = sx;
            if gx > fxgo {
                fxin = gx;
            }
        }
        let io = ((sx - fxso) / dxso) as i32;
        if min(sx, gx) >= fs && max(sx, gx) <= es && max(gx - sx, sx - gx) <= offmax {
            ng[io as usize] += 1;
            efwrite(tr.hdr.as_bytes(), &hdrfp);
            let as_ = (gx - fs) / ds;
            let mut is = as_ as i32;
            if is == ns - 1 {
                is = ns - 2;
            }
            let mut res = as_ - is as f32;
            if res as f64 <= 0.01 {
                res = 0.0;
            }
            if res as f64 >= 0.99 {
                res = 1.0;
            }
            sum2(
                nxt,
                nzt,
                1.0 - res,
                res,
                &ttab[is as usize * table_size..(is as usize + 1) * table_size],
                &ttab[(is + 1) as usize * table_size..(is as usize + 2) * table_size],
                &mut tsum,
            );
            dat2d(
                &mut tr.data,
                nt,
                ft,
                dt,
                sx,
                gx,
                &mut dats[io as usize * gather_size..(io as usize + 1) * gather_size],
                aperx,
                nxgo,
                fxin,
                dxgo,
                nzo as f32,
                fzo,
                dzo,
                nxi,
                fxi,
                dxi,
                &tsum,
                sgn,
                nzt,
                fzt,
                dzt,
                nxt,
                fxt,
                dxt,
                &szif,
                v0,
                freq,
                &mut filter,
            );
            ktr += 1;
            // ARM64 C integer remainder with a zero divisor returns the dividend.
            if (jtr - 1).checked_rem(mtr).unwrap_or(jtr - 1) == 0 {
                fprintf!(&jpfp, " Datumed receiver in trace %d\n", jtr);
                fflush(&jpfp);
            }
        }
        jtr += 1;
        if fgettr(&infp, &mut tr) == 0 || jtr > ntr {
            break;
        }
    }
    fprintf!(&jpfp, " Datumed receivers in %d total traces \n", ktr);
    erewind(&hdrfp);
    let mut scal = (dxgo as f64 * (1.0 / (2.0 * PI * v0 as f64).sqrt())) as f32;
    scal *= scale;
    for io in 0..nxso as usize {
        for ixo in 0..ng[io] as usize {
            efread(tro.hdr.as_bytes_mut(), &hdrfp);
            let start = io * gather_size + ixo * nt as usize;
            tro.data[..nt as usize].copy_from_slice(&dats[start..start + nt as usize]);
            for it in 0..nt as usize {
                tro.data[it] *= scal;
            }
            fputtr(&outfp, &mut tro);
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

// The original passes unchecked fopen results to stdio. Preserve its SIGSEGV
// on NULL instead of introducing a Rust panic or a new diagnostic.
fn unchecked_file(fp: Option<CFile>) -> CFile {
    fp.unwrap_or_else(|| {
        unsafe extern "C" {
            fn raise(sig: i32) -> i32;
        }
        unsafe {
            raise(11);
        }
        exit(139)
    })
}

// Sum of two tables.
fn sum2(nx: i32, nz: i32, a1: f32, a2: f32, t1: &[f32], t2: &[f32], t: &mut [Vec<f32>]) {
    for ix in 0..nx as usize {
        for iz in 0..nz as usize {
            t[ix][iz] = a1 * t1[ix * nz as usize + iz] + a2 * t2[ix * nz as usize + iz];
        }
    }
}

#[allow(clippy::too_many_arguments, unused_assignments)]
fn dat2d(
    trace: &mut [f32],
    nt: i32,
    ft: f32,
    dt: f32,
    sx: f32,
    gx: f32,
    dat: &mut [f32],
    aperx: f32,
    nxgo: i32,
    fxin: f32,
    dxgo: f32,
    nzo: f32,
    fzo: f32,
    dzo: f32,
    nxi: i32,
    fxi: f32,
    dxi: f32,
    tsum: &[Vec<f32>],
    sgn: i32,
    nzt: i32,
    _fzt: f32,
    dzt: f32,
    _nxt: i32,
    fxt: f32,
    dxt: f32,
    szif: &[[f32; 2]],
    v0: f32,
    freq: f32,
    filter: &mut Filter,
) {
    // Tolerance for numerical problems.
    let atol = 1e-10f32;
    let odt = (1.0 / dt as f64) as f32;
    let mut tzt = vec![0.0f32; nzt as usize];
    let h = (if gx >= sx {
        (gx - sx) as f64 / 2.0
    } else {
        (sx - gx) as f64 / 2.0
    }) as f32;
    filter.filt(trace, nt, dt, sgn);
    let mut nxf = ((gx - aperx - fxin) / dxgo) as i32;
    if nxf < 0 {
        nxf = 0;
    }
    let mut nxe = ((gx + aperx - fxin) / dxgo) as i32;
    if nxe >= nxgo {
        nxe = nxgo - 1;
    }

    let mut am = (fxin - fxi) / dxi;
    let mut mr = am as i32;
    am -= mr as f32;
    if am as f64 <= 0.01 {
        am = 0.0;
    }
    if am as f64 >= 0.99 {
        am = 1.0;
    }
    let mut am0 = (1.0 - am as f64) as f32;
    if mr < 0 {
        mr = 0;
    }
    // A single surface point makes the C endpoint adjustment set mr=-1.
    // Match the minimal C bounds fix by rejecting it before interpolation.
    if nxi < 2 || mr + nxe > nxi - 1 {
        err!("Topography definition is out of range!\n");
    }

    // Depth of input receiver on the recording surface.
    am = (gx - fxi) / dxi;
    mr = am as i32;
    am -= mr as f32;
    if am as f64 <= 0.01 {
        am = 0.0;
    }
    if am as f64 >= 0.99 {
        am = 1.0;
    }
    am0 = (1.0 - am as f64) as f32;
    if mr < 0 {
        mr = 0;
    }
    if mr >= nxi - 1 {
        mr -= 1;
        am0 = 0.0;
        am = 1.0;
    }
    let zi = abs(am0 * szif[mr as usize][0] + am * szif[(mr + 1) as usize][0]);
    if zi >= fzo + nzo * dzo {
        err!("Recording surface is out of range!\n");
    }
    // Approximate the surface derivative near this input location.
    let zp;
    let dzde;
    if mr == 0 {
        zp = szif[(mr + 1) as usize][0] - szif[mr as usize][0];
        dzde = zp / dxi;
    } else {
        zp = (szif[mr as usize][0] - szif[(mr - 1) as usize][0])
            + (szif[(mr + 1) as usize][0] - szif[mr as usize][0]);
        dzde = zp / (2.0f32 * dxi);
    }
    // Depth of source location.
    am = (sx - fxi) / dxi;
    mr = am as i32;
    am -= mr as f32;
    if am as f64 <= 0.01 {
        am = 0.0;
    }
    if am as f64 >= 0.99 {
        am = 1.0;
    }
    am0 = (1.0 - am as f64) as f32;
    if mr < 0 {
        mr = 0;
    }
    if mr >= nxi - 1 {
        mr -= 1;
        am0 = 0.0;
        am = 1.0;
    }
    let zsx = abs(am0 * szif[mr as usize][0] + am * szif[(mr + 1) as usize][0]);
    if zsx >= fzo + nzo * dzo {
        err!("Recording surface is out of range!\n");
    }

    // Loop over the lateral aperture in the output section.
    for ix in nxf..=nxe {
        let x = fxin + ix as f32 * dxgo;
        let mut ax = (x - fxt) / dxt;
        let jx = ax as i32;
        ax -= jx as f32;
        if ax as f64 <= 0.01 {
            ax = 0.0;
        }
        if ax as f64 >= 0.99 {
            ax = 1.0;
        }
        let ax0 = (1.0 - ax as f64) as f32;
        for iz in 0..nzt as usize {
            // At the last table column ax=0; the unused neighbor is outside
            // the allocation. Match the C fix by skipping that zero term.
            let right = if ax == 0.0 {
                0.0
            } else {
                ax * tsum[(jx + 1) as usize][iz]
            };
            tzt[iz] = ax0 * tsum[jx as usize][iz] + right;
        }
        // Datuming depth at output location.
        am = (x - fxi) / dxi;
        mr = am as i32;
        am -= mr as f32;
        if am as f64 <= 0.01 {
            am = 0.0;
        }
        if am as f64 >= 0.99 {
            am = 1.0;
        }
        am0 = (1.0 - am as f64) as f32;
        if mr < 0 {
            mr = 0;
        }
        if mr >= nxi - 1 {
            mr -= 1;
            am = 0.0;
            am0 = 1.0;
        }
        let z = abs(am0 * szif[mr as usize][1] + am * szif[(mr + 1) as usize][1]);
        if z <= fzo || z >= fzo + nzo * dzo {
            err!("Datuming surface is out of travel time range!\n");
        }
        let az = (z - fzo) / dzt;
        let mut jz = az as i32;
        if jz >= nzt - 1 {
            jz = nzt - 2;
        }
        let sz = az - jz as f32;
        let sz0 = (1.0 - sz as f64) as f32;
        let tio = sz0 * tzt[jz as usize] + sz * tzt[(jz + 1) as usize];
        let at = (sgn as f32 * tio - ft) * odt;
        let jt = at as i32;
        let res = abs(at - jt as f32);
        let res0 = (1.0 - res as f64) as f32;
        // Receivers on the greater side of the source form this gather.
        if x >= sx {
            let sing =
                ((((sx - x) * (sx - x) + (zsx - z) * (zsx - z)) as f64).sqrt() / v0 as f64) as f32;
            let phi = if x == sx {
                (PI / 2.0) as f32
            } else {
                (((z - zsx) / (x - sx)) as f64).atan() as f32
            };
            let xp = ((x - sx) as f64 * (phi as f64).cos() + (z - zsx) as f64 * (phi as f64).sin())
                as f32;
            let hp = (0.5 * xp as f64) as f32;
            let gxp = ((gx - sx) as f64 * (phi as f64).cos()
                + (zi - zsx) as f64 * (phi as f64).sin()) as f32;
            let zip = (-(gx - sx) as f64 * (phi as f64).sin()
                + (zi - zsx) as f64 * (phi as f64).cos()) as f32;
            if hp == 0.0 {
                err!("Datum and recording surface overlap!\n");
            }
            for it in 0..nt {
                if it + jt >= 0 && jt < nt - it - sgn {
                    let t0 = ft + it as f32 * dt;
                    let mut ampd;
                    if t0 < sing + atol {
                        ampd = 0.0f32;
                    } else {
                        let xst;
                        let zst;
                        if h == 0.0 {
                            xst = (hp as f64 + (v0 * t0) as f64 / 2.0) as f32;
                            zst = 0.0f32;
                        } else {
                            let ctau = (v0 * t0) * (v0 * t0);
                            let q = ((ctau as f64 - 4.0 * hp as f64 * hp as f64)
                                / (4.0 * ctau as f64)) as f32;
                            let mut p =
                                (2.0 * q as f64 * (gxp - xp) as f64 * xp as f64 * zip as f64)
                                    as f32;
                            p = (p as f64
                                / ((zip * zip) as f64
                                    + 4.0 * q as f64 * (gxp - xp) as f64 * (gxp - xp) as f64))
                                as f32;
                            let mut s = q * zip * zip * (ctau - (xp * xp));
                            s = (s as f64
                                / ((zip * zip) as f64
                                    + 4.0 * q as f64 * (gxp - xp) as f64 * (gxp - xp) as f64))
                                as f32;
                            zst = (((s + p * p) as f64).sqrt() - p as f64) as f32;
                            xst = ((gxp - xp) / zip) * zst + xp;
                        }
                        // Rotate back to unprimed coordinates.
                        let xs = (xst as f64 * (phi as f64).cos() - zst as f64 * (phi as f64).sin()
                            + sx as f64) as f32;
                        let zs = (xst as f64 * (phi as f64).sin()
                            + zst as f64 * (phi as f64).cos()
                            + zsx as f64) as f32;
                        let rig = (((xst - gxp) * (xst - gxp) + (zst - zip) * (zst - zip)) as f64)
                            .sqrt() as f32;
                        let rog = (((xst - xp) * (xst - xp) + zst * zst) as f64).sqrt() as f32;
                        let rs = ((xst * xst + zst * zst) as f64).sqrt() as f32;
                        let g = (zs - zi) - dzde * (xs - gx);
                        ampd = (g as f64 * ((rs + rig) as f64).sqrt() / rig as f64) as f32;
                        ampd = (ampd as f64
                            / (((rs + rog) as f64).sqrt() * (abs(rig - rog) as f64).sqrt()))
                            as f32;
                        // Check validity of stationary point.
                        if rig - rog < v0 / (4.0f32 * freq) {
                            ampd = 0.0;
                        }
                        if zs < 0.0 {
                            ampd = 0.0;
                        }
                        if rig <= rog {
                            ampd = 0.0;
                        }
                        if rig < 0.0 || rog < 0.0 || rs < 0.0 {
                            ampd = 0.0;
                        }
                    }
                    let temp = (res0 * trace[(it + jt) as usize]
                        + res * trace[(it + jt + sgn) as usize])
                        * ampd;
                    dat[ix as usize * nt as usize + it as usize] += temp;
                }
            }
        }
    }
}

// The C function's static FFT parameters, owned by this invocation of main.
#[derive(Default)]
struct Filter {
    nfft: i32,
    n: i32,
    dw: f32,
    fw: f32,
}
impl Filter {
    fn filt(&mut self, trace: &mut [f32], nt: i32, dt: f32, _sgn: i32) {
        let const2 = (0.5 * 2.0f64.sqrt()) as f32;
        if self.nfft == 0 {
            self.nfft = npfao(2 * nt, 4 * nt);
            if self.nfft >= SU_NFLTS as i32 || self.nfft >= 720720 {
                err!("Padded nt=%d -- too big", self.nfft);
            }
            self.n = self.nfft / 2;
            self.dw = (2.0 * PI / (self.nfft as f32 * dt) as f64) as f32;
            self.fw = (-PI / dt as f64) as f32;
        }
        let mut ct = vec![Complex::default(); self.nfft as usize];
        for it in 0..nt as usize {
            ct[it].r = trace[it];
        }
        pfacc(1, &mut ct);
        for iw in 0..=self.n as usize {
            let omega = iw as f32 * self.dw;
            let amp = (abs(omega) as f64).sqrt() as f32;
            let sgnw = 1.0f32;
            let temp = (ct[iw].r - sgnw * ct[iw].i) * amp * const2;
            ct[iw].i = (sgnw * ct[iw].r + ct[iw].i) * amp * const2;
            ct[iw].r = temp;
        }
        for iw in (self.n + 1) as usize..self.nfft as usize {
            let omega = self.fw + (iw as i32 - self.n) as f32 * self.dw;
            let amp = (abs(omega) as f64).sqrt() as f32;
            let sgnw = -1.0f32;
            let temp = (ct[iw].r - sgnw * ct[iw].i) * amp * const2;
            ct[iw].i = (sgnw * ct[iw].r + ct[iw].i) * amp * const2;
            ct[iw].r = temp;
        }
        pfacc(-1, &mut ct);
        for it in 0..nt as usize {
            trace[it] = ct[it].r / self.nfft as f32;
        }
    }
}
