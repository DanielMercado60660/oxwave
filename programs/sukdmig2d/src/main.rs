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

//! SUKDMIG2D - Kirchhoff Depth Migration of 2D poststack/prestack data.
//! Port of su/main/migration_inversion/sukdmig2d.c (Revision 1.26).
//! Credits: Zhenyue Liu, 03/01/95, Colorado School of Mines.
//! Modifications: Gary Billings, Talisman Energy, Sept 2005: absoff, limoff.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{max, min, nint, PI};
mod migration;
use migration::{mig2d, resit, sum2, timeb, Filter};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    "SUKDMIG2D - Kirchhoff Depth Migration of 2D poststack/prestack data\t",
    " \t\t\t\t\t\t\t\t\t",
    "    sukdmig2d  infile=  outfile=  ttfile=   [parameters] \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " infile=stdin\t\tfile for input seismic traces\t\t\t",
    " outfile=stdout\tfile for common offset migration output  \t",
    " ttfile=\t\tfile for input traveltime tables\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " ...  The following 9 parameters describe traveltime tables:\t\t",
    " fzt= \t\t\tfirst depth sample in traveltime table\t\t",
    " nzt= \t\t\tnumber of depth samples in traveltime table\t",
    " dzt=\t\t\tdepth interval in traveltime table\t\t",
    " fxt=\t\t\tfirst lateral sample in traveltime table\t",
    " nxt=\t\t\tnumber of lateral samples in traveltime table\t",
    " dxt=\t\t\tlateral interval in traveltime table\t\t",
    " fs= \t\t\tx-coordinate of first source\t\t\t",
    " ns= \t\t\tnumber of sources\t\t\t\t",
    " ds= \t\t\tx-coordinate increment of sources\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " dt= or from header (dt) \ttime sampling interval of input data\t",
    " ft= or from header (ft) \tfirst time sample of input data\t\t",
    " dxm= or from header (d2) \tsampling interval of midpoints \t\t",
    " fzo=fzt\t\t    z-coordinate of first point in output trace\t",
    " dzo=0.2*dzt\t\tvertical spacing of output trace \t\t",
    " nzo=5*(nzt-1)+1 \tnumber of points in output trace\t\t",
    " fxo=fxt\t\t    x-coordinate of first output trace \t\t",
    " dxo=0.5*dxt\t\thorizontal spacing of output trace \t\t",
    " nxo=2*(nxt-1)+1  \tnumber of output traces \t\t\t",
    " off0=0\t\t   \tfirst offest in output \t\t\t",
    " doff=99999\t\toffset increment in output \t\t\t",
    " noff=1\t \tnumber of offsets in output \t\t\t",
    " absoff=0      flag for using absolute offsets of input traces\t\t",
    "               =0 means use offset=gx-sx                \t\t",
    "               =1 means use abs(gx-sx)                  \t\t",
    " limoff=0      flag for only using input traces that fall within the range",
    "               of defined output offset bins (off0,doff,noff) \t\t",
    "               =0 means use all input traces                 \t\t",
    "               =1 means limit traces used by offset           \t\t",
    " fmax=0.25/dt\t\tfrequency-highcut for input traces\t\t",
    " offmax=99999\t\tmaximum absolute offset allowed in migration \t",
    " aperx=nxt*dxt/2  \tmigration lateral aperature \t\t\t",
    " angmax=60\t\tmigration angle aperature from vertical \t",
    " v0=1500(m/s)\t\treference velocity value at surface\t\t",
    " dvz=0.0  \t\treference velocity vertical gradient\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " ls=1\t\t\tflag for line source\t\t\t\t",
    " jpfile=stderr\t\tjob print file name \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " mtr=100  \t\tprint verbal information at every mtr traces\t",
    " ntr=100000\t\tmaximum number of input traces to be migrated\t",
    " npv=0\t\t\tflag of computing quantities for velocity analysis",
    " rscale=1000.0 \tscaling for roundoff error suppression\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   ...if npv>0 specify the following three files:\t\t\t",
    " tvfile=tvfile\t\tinput file of traveltime variation tables\t",
    "\t\t\ttv[ns][nxt][nzt]\t\t\t\t",
    " csfile=csfile\t\tinput file of cosine tables cs[ns][nxt][nzt]\t",
    " outfile1=outfile1\tfile containning additional migration output   \t",
    "\t\t\twith extra amplitude\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. Traveltime tables were generated by program rayt2d (or other ones)\t",
    "    on relatively coarse grids, with dimension ns*nxt*nzt. In the\t",
    "    migration process, traveltimes are interpolated into shot/gephone \t",
    "    positions and output grids.\t\t\t\t\t",
    " 2. Input seismic traces must be SU format and can be any type of \t",
    "    gathers (common shot, common offset, common CDP, and so on).\t",
    " 3. Migrated traces are output in CDP gathers if velocity analysis\t",
    "    is required, with dimension nxo*noff*nzo.  \t\t\t",
    " 4. If the offset value of an input trace is not in the offset array \t",
    "    of output, the nearest one in the array is chosen. \t\t",
    " 5. Memory requirement for this program is about\t\t\t",
    "    \t[ns*nxt*nzt+noff*nxo*nzo+4*nr*nzt+5*nxt*nzt+npa*(2*ns*nxt*nzt   ",
    "\t+noff*nxo*nzo+4*nxt*nzt)]*4 bytes\t\t\t\t",
    "    where nr = 1+min(nxt*dxt,0.5*offmax+aperx)/dxo. \t\t\t",
    " 6. Amplitudes are computed using the reference velocity profile, v(z),",
    "    specified by the parameters v0= and dvz=.\t\t\t\t",
    " 7. Input traces must specify source and receiver positions via the header",
    "    fields tr.sx and tr.gx. Offset is computed automatically.\t\t",
    " 8. if limoff=0, input traces from outside the range defined by off0, doff, ",
    "    noff, will get migrated into the extremal offset bins/planes.  E.g. if ",
    "    absoff=0 and limoff=0, all traces with gx<sx will get migrated into the ",
    "    off0 bin.",
    "\t\t\t\t\t\t\t\t\t",
];

// Parse parameters in the same order as the C entry point.
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    let mut infile = String::from("stdin");
    let mut outfile = String::from("stdout");
    let infp = if getparstring("infile", &mut infile) == 0 {
        stdin()
    } else {
        fopen(&infile, "r").unwrap_or_else(|| err!("cannot open infile=%s\n", infile))
    };
    let outfp = if getparstring("outfile", &mut outfile) == 0 {
        stdout()
    } else {
        fopen(&outfile, "w").unwrap_or_else(null_file)
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
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };
    if fgettr(&infp, &mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if (dt as f64) < 0.0000001 {
        err!("dt must be positive!\n");
    }
    let mut ft = 0.0f32;
    if getparfloat("ft", &mut ft) == 0 {
        ft = (tr.delrt as f64 / 1000.0) as f32;
    }
    let mut dxm = 0.0f32;
    if getparfloat("dxm", &mut dxm) == 0 {
        dxm = tr.d2;
    }
    if (dxm as f64) < 0.0000001 {
        err!("dxm must be positive!\n");
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
    let mut rscale = 1000.0f32;
    getparfloat("rscale", &mut rscale);
    let ext = nint((rscale * (fxt + (nxt - 1) as f32 * dxt)) as f64) as f32 / rscale;
    let ezt = nint((rscale * (fzt + (nzt - 1) as f32 * dzt)) as f64) as f32 / rscale;
    let es = nint((rscale * (fs + (ns - 1) as f32 * ds)) as f64) as f32 / rscale;
    let mut nxo = (nxt - 1) * 2 + 1;
    getparint("nxo", &mut nxo);
    let mut fxo = fxt;
    getparfloat("fxo", &mut fxo);
    let mut dxo = (dxt as f64 * 0.5) as f32;
    getparfloat("dxo", &mut dxo);
    let mut nzo = (nzt - 1) * 5 + 1;
    getparint("nzo", &mut nzo);
    let mut fzo = fzt;
    getparfloat("fzo", &mut fzo);
    let mut dzo = (dzt as f64 * 0.2) as f32;
    getparfloat("dzo", &mut dzo);
    let exo = nint((rscale * (fxo + (nxo - 1) as f32 * dxo)) as f64) as f32 / rscale;
    let ezo = nint((rscale * (fzo + (nzo - 1) as f32 * dzo)) as f64) as f32 / rscale;
    let jpfp = jpfp.unwrap_or_else(null_file);
    fprintf!(&jpfp, " fxt=%f fxo=%f \n", fxt, fxo);
    fprintf!(&jpfp, " ext=%f exo=%f \n", ext, exo);
    fprintf!(&jpfp, " fzt=%f fzo=%f \n", fzt, fzo);
    fprintf!(&jpfp, " ezt=%f ezo=%f \n", ezt, ezo);
    fprintf!(&jpfp, " es=%f \n", es);
    fprintf!(&jpfp, " \n");
    if fxt > fxo || ext < exo || fzt > fzo || ezt < ezo {
        warn!("This condition must NOT be satisfied: fxt>fxo || ext<exo || fzt>fzo || ezt<ezo");
        warn!(
            "fxt=%f fxo=%f ext=%f exo=%f fzt=%f fzo=%f ezt=%f ezo=%f",
            fxt, fxo, ext, exo, fzt, fzo, ezt, ezo
        );
        err!(" migration output range is out of traveltime table!\n");
    }
    let mut v0 = 1500.0f32;
    getparfloat("v0", &mut v0);
    let mut dvz = 0.0f32;
    getparfloat("dvz", &mut dvz);
    let mut angmax = 60.0f32;
    getparfloat("angmax", &mut angmax);
    if (angmax as f64) < 0.00001 {
        err!("angmax must be positive!\n");
    }
    let mut mtmax =
        ((2.0f32 * dxm) as f64 * (angmax as f64 * PI / 180.0).sin() / (v0 * dt) as f64) as i32;
    if mtmax < 1 {
        mtmax = 1;
    }
    let mut aperx = (0.5 * nxt as f64 * dxt as f64) as f32;
    getparfloat("aperx", &mut aperx);
    let mut offmax = 3000.0f32;
    getparfloat("offmax", &mut offmax);
    let mut fmax = (0.25 / dt as f64) as f32;
    getparfloat("fmax", &mut fmax);
    let mut noff = 1i32;
    getparint("noff", &mut noff);
    let mut off0 = 0.0f32;
    getparfloat("off0", &mut off0);
    let mut doff = 99999.0f32;
    getparfloat("doff", &mut doff);
    let mut ls = 1i32;
    getparint("ls", &mut ls);
    let mut absoff = 0i32;
    getparint("absoff", &mut absoff);
    let mut limoff = 0i32;
    getparint("limoff", &mut limoff);
    let mut ntr = 100000i32;
    getparint("ntr", &mut ntr);
    let mut mtr = 100i32;
    getparint("mtr", &mut mtr);
    let mut npv = 0i32;
    getparint("npv", &mut npv);
    let mut tvfile = String::new();
    let mut csfile = String::new();
    let mut outfile1 = String::from("outfile1");
    let (mut tvfp, mut csfp, mut out1fp) = (None, None, None);
    if npv != 0 {
        if getparstring("tvfile", &mut tvfile) == 0 {
            err!("must specify tvfile!\n");
        }
        tvfp = fopen(&tvfile, "r");
        if getparstring("csfile", &mut csfile) == 0 {
            err!("must specify csfile!\n");
        }
        csfp = fopen(&csfile, "r");
        getparstring("outfile1", &mut outfile1);
        out1fp = fopen(&outfile1, "w");
    }
    checkpars();
    fprintf!(&jpfp, "\n");
    fprintf!(&jpfp, " Migration parameters\n");
    fprintf!(&jpfp, " ================\n");
    fprintf!(&jpfp, " infile=%s \n", infile);
    fprintf!(&jpfp, " outfile=%s \n", outfile);
    fprintf!(&jpfp, " ttfile=%s \n", ttfile);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzt=%d fzt=%g dzt=%g\n", nzt, fzt, dzt);
    fprintf!(&jpfp, " nxt=%d fxt=%g dxt=%g\n", nxt, fxt, dxt);
    fprintf!(&jpfp, " ns=%d fs=%g ds=%g\n", ns, fs, ds);
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " nzo=%d fzo=%g dzo=%g\n", nzo, fzo, dzo);
    fprintf!(&jpfp, " nxo=%d fxo=%g dxo=%g\n", nxo, fxo, dxo);
    fprintf!(&jpfp, " \n");

    let mut rmax = max(es - fxt, ext - fs);
    rmax = min(rmax as f64, 0.5 * offmax as f64 + aperx as f64) as f32;
    let nr = 2 + (rmax / dxo) as i32;
    let mut tb = vec![vec![0.0f32; nzt as usize]; nr as usize];
    let mut pb = tb.clone();
    let mut cs0b = tb.clone();
    let mut angb = tb.clone();
    timeb(
        nr as usize,
        nzt as usize,
        dxo,
        dzt,
        fzt,
        dvz,
        v0,
        &mut tb,
        &mut pb,
        &mut cs0b,
        &mut angb,
    );
    fprintf!(&jpfp, " nt=%d ft=%g dt=%g \n", nt, ft, dt);
    fprintf!(&jpfp, " dxm=%g fmax=%g\n", dxm, fmax);
    fprintf!(&jpfp, " noff=%d off0=%g doff=%g\n", noff, off0, doff);
    fprintf!(&jpfp, " v0=%g dvz=%g \n", v0, dvz);
    fprintf!(
        &jpfp,
        " aperx=%g offmax=%g angmax=%g\n",
        aperx,
        offmax,
        angmax
    );
    fprintf!(&jpfp, " ntr=%d mtr=%d ls=%d npv=%d\n", ntr, mtr, ls, npv);
    fprintf!(&jpfp, " absoff=%d limoff=%d\n", absoff, limoff);
    if npv != 0 {
        fprintf!(
            &jpfp,
            " tvfile=%s csfile=%s outfile1=%s\n",
            tvfile,
            csfile,
            outfile1
        );
    }
    fprintf!(&jpfp, " ================\n");
    fflush(&jpfp);

    let mut mig = vec![vec![vec![0.0f32; nzo as usize]; nxo as usize]; noff as usize];
    let mut ttab = vec![vec![0.0f32; (nxt * nzt) as usize]; ns as usize];
    let mut tt = vec![0.0f32; (nxt * nzt) as usize];
    let mut tsum = tt.clone();
    let mut tv = if npv != 0 { ttab.clone() } else { Vec::new() };
    let mut cs = if npv != 0 { ttab.clone() } else { Vec::new() };
    let mut tvsum = if npv != 0 { tt.clone() } else { Vec::new() };
    let mut cssum = if npv != 0 { tt.clone() } else { Vec::new() };
    let mut mig1 = if npv != 0 {
        mig.clone()
    } else {
        vec![vec![vec![0.0f32]]; noff as usize]
    };
    fprintf!(&jpfp, " input traveltime tables \n");
    for is in 0..ns as usize {
        let nseek = nxt as i64 * nzt as i64 * is as i64;
        efseeko(&ttfp, nseek * 4, SEEK_SET);
        fread(&mut ttab[is], &ttfp);
        let s = fs + is as f32 * ds;
        resit(
            nxt as usize,
            fxt,
            dxt,
            nzt as usize,
            nr as usize,
            dxo,
            &tb,
            &mut ttab[is],
            s,
        );
        if npv != 0 {
            let tvfp = tvfp.as_ref().unwrap_or_else(|| null_file_ref());
            let csfp = csfp.as_ref().unwrap_or_else(|| null_file_ref());
            efseeko(tvfp, nseek * 4, SEEK_SET);
            fread(&mut tv[is], tvfp);
            efseeko(csfp, nseek * 4, SEEK_SET);
            fread(&mut cs[is], csfp);
        }
    }
    fprintf!(&jpfp, " start migration ... \n");
    fprintf!(&jpfp, " \n");
    fflush(&jpfp);
    let mut jtr = 1i32;
    let mut ktr = 0i32;
    fprintf!(&jpfp, " fs=%g es=%g offmax=%g\n", fs, es, offmax);
    let mut filter = Filter::default();
    loop {
        let (sx, gx) = if tr.scalco != 0 {
            if tr.scalco > 0 {
                // C multiplies as i32 before converting to float; preserve the
                // reference binary's wraparound even with overflow checks on.
                (
                    tr.sx.wrapping_mul(tr.scalco as i32) as f32,
                    tr.gx.wrapping_mul(tr.scalco as i32) as f32,
                )
            } else {
                // Both operands are integers in C, including this division.
                (
                    (tr.sx / (-(tr.scalco as i32))) as f32,
                    (tr.gx / (-(tr.scalco as i32))) as f32,
                )
            }
        } else {
            (tr.sx as f32, tr.gx as f32)
        };
        let mut offset = (gx - sx) as i32;
        if absoff != 0 && offset < 0 {
            offset = offset.wrapping_neg();
        }
        let mut io = (((offset as f32 - off0) / doff) as f64 + 0.5) as i32;
        // A C do/while continue evaluates its condition and does not increment jtr.
        if !(limoff != 0 && (io < 0 || io >= noff)) {
            if io < 0 {
                io = 0;
            }
            if io >= noff {
                io = noff - 1;
            }
            if min(sx, gx) >= fs && max(sx, gx) <= es && max(gx - sx, sx - gx) <= offmax {
                let as_ = (sx - fs) / ds;
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
                    (1.0 - res as f64) as f32,
                    res,
                    &ttab[is as usize],
                    &ttab[(is + 1) as usize],
                    &mut tsum,
                );
                if npv != 0 {
                    sum2(
                        (1.0 - res as f64) as f32,
                        res,
                        &tv[is as usize],
                        &tv[(is + 1) as usize],
                        &mut tvsum,
                    );
                    sum2(
                        (1.0 - res as f64) as f32,
                        res,
                        &cs[is as usize],
                        &cs[(is + 1) as usize],
                        &mut cssum,
                    );
                }
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
                    (1.0 - res as f64) as f32,
                    res,
                    &ttab[is as usize],
                    &ttab[(is + 1) as usize],
                    &mut tt,
                );
                for i in 0..tsum.len() {
                    tsum[i] = 1.0 * tt[i] + 1.0 * tsum[i];
                }
                if npv != 0 {
                    sum2(
                        (1.0 - res as f64) as f32,
                        res,
                        &tv[is as usize],
                        &tv[(is + 1) as usize],
                        &mut tt,
                    );
                    for i in 0..tvsum.len() {
                        tvsum[i] = 1.0 * tt[i] + 1.0 * tvsum[i];
                    }
                    sum2(
                        (1.0 - res as f64) as f32,
                        res,
                        &cs[is as usize],
                        &cs[(is + 1) as usize],
                        &mut tt,
                    );
                    for i in 0..cssum.len() {
                        cssum[i] = 1.0 * tt[i] + 1.0 * cssum[i];
                    }
                }
                mig2d(
                    &mut tr.data,
                    nt,
                    ft,
                    dt,
                    sx,
                    gx,
                    &mut mig[io as usize],
                    aperx,
                    nxo,
                    fxo,
                    dxo,
                    nzo,
                    fzo,
                    dzo,
                    ls,
                    mtmax,
                    dxm,
                    fmax,
                    angmax,
                    &tb,
                    &pb,
                    &cs0b,
                    &angb,
                    nr,
                    &tsum,
                    nzt,
                    fzt,
                    dzt,
                    nxt,
                    fxt,
                    dxt,
                    npv,
                    &cssum,
                    &tvsum,
                    &mut mig1[io as usize],
                    &mut filter,
                );
                ktr += 1;
                // arm64 C remainder with a zero divisor returns the dividend.
                if (jtr - 1).checked_rem(mtr).unwrap_or(jtr - 1) == 0 {
                    fprintf!(&jpfp, " migrated trace %d\n", jtr);
                    fflush(&jpfp);
                }
            }
            jtr += 1;
        }
        if fgettr(&infp, &mut tr) == 0 || jtr >= ntr {
            break;
        }
    }
    fprintf!(&jpfp, " migrated %d traces in total\n", ktr);
    let mut tro = Trace::new();
    tro.ns = nzo as u16;
    tro.d1 = dzo;
    tro.dt = (1000 * ((1000.0f32 * dt) as f64 + 0.5) as i32) as u16;
    tro.delrt = 0;
    tro.f1 = fzo;
    tro.f2 = fxo;
    tro.d2 = dxo;
    tro.trid = 200;
    let scal = (4.0 / PI.sqrt() * dxm as f64 / v0 as f64) as f32;
    for ixo in 0..nxo as usize {
        for io in 0..noff as usize {
            tro.data[..nzo as usize].copy_from_slice(&mig[io][ixo]);
            tro.offset = (off0 + io as f32 * doff) as i32;
            tro.tracr = 1 + ixo as i32;
            tro.tracl = 1 + io as i32 + noff * ixo as i32;
            tro.cdp = (fxo + ixo as f32 * dxo) as i32;
            tro.cdpt = 1 + io as i32;
            for izo in 0..nzo as usize {
                tro.data[izo] *= scal;
            }
            fputtr(&outfp, &mut tro);
            if npv != 0 {
                tro.data[..nzo as usize].copy_from_slice(&mig1[io][ixo]);
                for izo in 0..nzo as usize {
                    tro.data[izo] *= scal;
                }
                fputtr(out1fp.as_ref().unwrap_or_else(|| null_file_ref()), &mut tro);
            }
        }
    }
    fprintf!(&jpfp, " \n");
    fprintf!(&jpfp, " output done\n");
    fflush(&jpfp);
    efclose(&jpfp);
    efclose(&outfp);
    exit(CWP_Exit());
}

// The original dereferences unchecked fopen results on these error paths.
// Preserve its SIGSEGV termination rather than introducing a Rust panic message.
fn null_file() -> CFile {
    null_file_ref();
    unreachable!()
}
fn null_file_ref<'a>() -> &'a CFile {
    unsafe extern "C" {
        fn raise(sig: i32) -> i32;
    }
    unsafe {
        raise(11);
    }
    exit(139);
}
