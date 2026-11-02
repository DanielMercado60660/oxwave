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

//! SUKDMIG3D - Kirchhoff Depth Migration of 3D poststack/prestack data.
//! Port of su/main/migration_inversion/sukdmig3d.c (Revision 1.7).
//! Author: Zhaobo Meng, 01/10/97, Colorado School of Mines.

#![allow(unused_parens, unused_mut, unused_variables, unused_assignments)]
use oxwave_core::cio::{fclose, fflush, fopen, fread};
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::PI;
use oxwave_pfa_fft::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    "SUKDMIG3D - Kirchhoff Depth Migration of 3D poststack/prestack data\t",
    " \t\t\t\t\t\t\t\t\t",
    " \tsukdmig3d datain= dataout= [parameters] \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " ttfile\t file for input tttables\t\t\t       \t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " datain=stdin\t file for input seismic traces\t\t\t\t",
    " dataout=stdout file for common offset migration output\t\t",
    " crfile=NULL    file for cos theta and ray paths                       ",
    "\t\t\t\t\t\t\t\t\t",
    "   The following 17 parameters describe tttables: (from ttfile header)\t",
    " fxgd= or from header (f1)\t\tfirst x-sample in tttable\t",
    " nxt= or from header (ns)\t\tnumber of x-samples in tttable\t",
    " dxgd= or from header (d1)\t\tx-interval in tttable \t\t",
    " fygd= or from header (f2)\t\tfirst y-sample in tttable \t",
    " nyt= or from header (ntr)\t\tnumber of y-samples in tttable\t",
    " dygd= or from header (d2)\t\ty-interval in tttable\t\t",
    " ixsf= or from header (sdel) \t        x in dxgd of first source\t",
    " nxs= or from header (nhs) \t\tnumber of sources in x\t\t",
    " ixsr= or from header (swdep) \t        ratio of source & gd spacing    ",
    " iysf= or from header (gdel)\t        y in dygd of first source \t",
    " nys= or from header (nvs)\t\tnumber of sources in y          ",
    " iysr= or from header (gwdep)\t        ratio of source & gd spacing    ",
    " fzs= or from header (sdepth/1000)\tfirst depth sample in tttable\t",
    " nzs= or from header (duse) \t\tnumber of depth samples in tttable",
    " dzs= or from header (ep/1000)\t\tdepth interval in tttable\t",
    " nxgd= or from header (selev)\t\tx size of the traveltime region ",
    " nygd= or from header (gelev)\t\ty size of the traveltime region ",
    " multit= or from header (scalel)       number of multivalued traveltime",
    "\t\t\t\t\t\t\t\t\t",
    "   The following two parameters are from data header\t\t\t",
    " dt= or from header (dt)\ttime sampling interval of input data\t",
    " ft= or from header (ft)\tfirst time sample of input data\t\t",
    " dxm= or from header (d2)      mid point spacing of input data         ",
    "\t\t\t\t\t\t\t\t\t",
    " Default: output is 5 times finer in depth and 2 times finer in lateral",
    " fzo=fzs\tz-coordinate of first point in output trace\t\t",
    " dzo=0.2*dzs\tvertical spacing of output trace (5 times finer)\t",
    " nzo=5*(nzs-1)+1 number of points in output trace (5 times finer)\t",
    " fxo=fxgd\tx-coordinate of first output trace\t \t\t",
    " dxo=0.5*dxgd\thorizontal spacing of output trace (2 times finer)\t",
    " nxo=2*(nxgd-1)+1 number of output traces (2 times finer)\t\t",
    " fyo=fygd\ty-coordinate of first output trace\t\t\t",
    " dyo=0.5*dygd\thorizontal spacing of output trace (2 times finer)\t",
    " nyo=2*(nygd-1)+1 number of output traces (2 times finer)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Default: poststack migration\t\t\t\t\t\t",
    " fxoffset=0\t\tfirst offest in output in x\t\t\t",
    " fyoffset=0\t\tfirst offest in output in y\t\t\t",
    " dxoffset=99999\toffset increment in output in x\t \t\t",
    " dyoffset=99999\toffset increment in output in y\t\t\t",
    " nxoffset=1\t\tnumber of offsets in output in x\t\t",
    " nyoffset=1\t\tnumber of offsets in output in y\t\t",
    " xoffsetmax=99999\tx-maximum absolute offset allowed in migration \t",
    " yoffsetmax=99999\ty-maximum absolute offset allowed in migration  ",
    " xaper=nxt*dxgd/2.5\tmigration lateral aperature in x\t\t",
    " yaper=nyt*dygd/2.5\tmigration lateral aperature in y\t\t",
    " angmax=60             max angle to handle                             ",
    " fmax=0.25/dt          max frequency in the data                       ",
    " jpfile=stderr\t\tjob print file name \t\t\t\t",
    " pptr=100\t\tprint verbal information at every pptr traces\t",
    " ntrmax=100000\t\tmaximum number of input traces to be migrated\t",
    " ls=0                  point =0 line source =1                         ",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " 1. Traveltime tables were generated by program SUTETRARAY (or other\t",
    "    ones) on very sparse tetrahedral model, with dimension nys*nxs*nzs ",
    "    *nyt*nxt.                                                          ",
    " 2. Input seismic traces must be SU format and can be any type of \t",
    "    gathers (common shot, common offset, common CDP, and so on).\t",
    " 3. Migrated traces are output in CDP gathers if velocity analysis\t",
    "    is required, with dimension nyoffset*nxoffset*nyo*nxo*nzo.\t\t",
    " 4. If the offset value of an input trace is not in the offset array \t",
    "    of output, the nearest one in the array is chosen. \t\t",
    " 5. Memory requirement for this program is about\t\t\t",
    "    [nys*nxs*nzs*nyt*nxt+nyoffset*nxoffset*nxo*nyo*nzo+\t       \t",
    "    nys*nxo*nzo*nyt*nxt]                                               ",
    " 6. Input traces must specify source and receiver positions via header\t",
    "    fields tr.sx and tr.gx, as well as tr.sy and tr.gy. Offset is \t",
    "    computed automatically.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Disclaimer:\t\t\t\t\t\t\t\t",
    " This is a research code that will take considerable work to get into\t",
    " the form of a a production level 3D migration code. The code is\t",
    " offered as is, along with tetramod and sutetraray, to provide a starting",
    " point for researchers who wish to write their own 3D migration codes.",
    "\t\t\t\t\t\t\t\t\t",
];

struct Gd {
    fxgd: f32,
    dxgd: f32,
    nxgd: i32,
    fygd: f32,
    dygd: f32,
    nygd: i32,
    nxt: i32,
    nyt: i32,
    nxo: i32,
    fxo: f32,
    dxo: f32,
    nyo: i32,
    fyo: f32,
    dyo: f32,
    nzo: i32,
    fzo: f32,
    dzo: f32,
    fxs: f32,
    dxs: f32,
    nxs: i32,
    fys: f32,
    dys: f32,
    nys: i32,
    fzs: f32,
    dzs: f32,
    nzs: i32,
    xaper: f32,
    yaper: f32,
    multit: i32,
    ixsf: i32,
    ixsr: i32,
    iysf: i32,
    iysr: i32,
    sy: f32,
    gy: f32,
    sx: f32,
    gx: f32,
    nt: i32,
    ft: f32,
    dt: f32,
    v0: f32,
    fmax: f32,
    ls: i32,
    ntrimax: i32,
    crfp: Option<CFile>,
    jpfp: CFile,
}
impl Gd {
    fn new() -> Self {
        Self {
            fxgd: 0.0,
            dxgd: 0.0,
            nxgd: 0,
            fygd: 0.0,
            dygd: 0.0,
            nygd: 0,
            nxt: 0,
            nyt: 0,
            nxo: 0,
            fxo: 0.0,
            dxo: 0.0,
            nyo: 0,
            fyo: 0.0,
            dyo: 0.0,
            nzo: 0,
            fzo: 0.0,
            dzo: 0.0,
            fxs: 0.0,
            dxs: 0.0,
            nxs: 0,
            fys: 0.0,
            dys: 0.0,
            nys: 0,
            fzs: 0.0,
            dzs: 0.0,
            nzs: 0,
            xaper: 0.0,
            yaper: 0.0,
            multit: 0,
            ixsf: 0,
            ixsr: 0,
            iysf: 0,
            iysr: 0,
            sy: 0.0,
            gy: 0.0,
            sx: 0.0,
            gx: 0.0,
            nt: 0,
            ft: 0.0,
            dt: 0.0,
            v0: 0.0,
            fmax: 0.0,
            ls: 0,
            ntrimax: 0,
            crfp: None,
            jpfp: stderr(),
        }
    }
}

// mig3d - migrate a single trace.
fn mig3d(
    gd: &Gd,
    trf: &[f32],
    mig: &mut [f32],
    tbuf: &Grid<u16, 6>,
    cbuf: &Grid<u8, 5>,
    rbuf: &Grid<u8, 5>,
) {
    let mut nyof: i32 = 0;
    let mut nyoe: i32 = 0;
    let mut nxof: i32 = 0;
    let mut nxoe: i32 = 0;
    let mut iyo: i32 = 0;
    let mut ixo: i32 = 0;
    let mut ym: f32 = 0.0;
    let mut xm: f32 = 0.0;
    let mut yo: f32 = 0.0;
    let mut xo: f32 = 0.0;
    let mut xs: f32 = 0.0;
    let mut ys: f32 = 0.0;
    let mut ybeta: f32 = 0.0;
    let mut ybeta0: f32 = 0.0;
    let mut xalpha: f32 = 0.0;
    let mut xalpha0: f32 = 0.0;
    let mut ixs: i32 = 0;
    let mut iys: i32 = 0;
    let mut izo: i32 = 0;
    let mut xt: f32 = 0.0;
    let mut yt: f32 = 0.0;
    let mut xs0a: f32 = 0.0;
    let mut xs0a0: f32 = 0.0;
    let mut xg0a: f32 = 0.0;
    let mut xg0a0: f32 = 0.0;
    let mut xs1a: f32 = 0.0;
    let mut xs1a0: f32 = 0.0;
    let mut xg1a: f32 = 0.0;
    let mut xg1a0: f32 = 0.0;
    let mut ys0b: f32 = 0.0;
    let mut ys0b0: f32 = 0.0;
    let mut yg0b: f32 = 0.0;
    let mut yg0b0: f32 = 0.0;
    let mut ys1b: f32 = 0.0;
    let mut ys1b0: f32 = 0.0;
    let mut yg1b: f32 = 0.0;
    let mut yg1b0: f32 = 0.0;
    let mut tiys0: f32 = 0.0;
    let mut tiys1: f32 = 0.0;
    let mut ttotal: f32 = (0f64 as f32);
    let mut tiyg0: f32 = 0.0;
    let mut tiyg1: f32 = 0.0;
    let mut ciys0: f32 = 0.0;
    let mut ciys1: f32 = 0.0;
    let mut riys0: f32 = 0.0;
    let mut riys1: f32 = 0.0;
    let mut riyg0: f32 = 0.0;
    let mut riyg1: f32 = 0.0;
    let mut rtotal: f32 = 0.0;
    let mut cosine: f32 = 0.0;
    let mut itintr: i32 = 0;
    let mut tintr: f32 = 0.0;
    let mut wtintr: f32 = 0.0;
    let mut wtintr0: f32 = 0.0;
    let mut amp: f32 = (1f64 as f32);
    let mut ixts0: i32 = 0;
    let mut ixtg0: i32 = 0;
    let mut ixts1: i32 = 0;
    let mut ixtg1: i32 = 0;
    let mut iyts0: i32 = 0;
    let mut iytg0: i32 = 0;
    let mut iyts1: i32 = 0;
    let mut iytg1: i32 = 0;
    let mut rs: f32 = 0.0;
    let mut rg: f32 = 0.0;
    let mut zo: f32 = 0.0;
    let mut im: i32 = 0;
    let mut rhos: Grid<f32, 2>;
    let mut rhog: Grid<f32, 2>;
    let mut ntri: i32 = 0;
    let mut drho: f32 = 0.0;
    let mut ddxg: f32 = 0.0;
    let mut ddyg: f32 = 0.0;
    let mut ddxs: f32 = 0.0;
    let mut ddys: f32 = 0.0;
    let mut temp: f32 = 0.0;
    let mut ts: f32 = (0f64 as f32);
    let mut tg: f32 = (0f64 as f32);
    rhos = Grid::new([gd.nyo, gd.nxo]);
    rhog = Grid::new([gd.nyo, gd.nxo]);
    ym = ((0.5f64 * ((gd.sy + gd.gy) as f64)) as f32);
    xm = ((0.5f64 * ((gd.sx + gd.gx) as f64)) as f32);
    fprintf!(&gd.jpfp, "xm,ym=%g %g\n", (xm as f64), (ym as f64));
    nyof = ((((((ym - gd.yaper) - gd.fyo) / gd.dyo) as f64) + 0.5f64) as i32);
    if (nyof < 0i32) {
        nyof = 0i32;
    }
    nyoe = ((((((ym + gd.yaper) - gd.fyo) / gd.dyo) as f64) + 0.5f64) as i32);
    if (nyoe >= gd.nyo) {
        nyoe = (gd.nyo - 1i32);
    }
    fprintf!(
        &gd.jpfp,
        "This trace: fyo,dyo,nyof,nyoe=%g %g %d %d\n",
        (gd.fyo as f64),
        (gd.dyo as f64),
        nyof,
        nyoe
    );
    nxof = ((((((xm - gd.xaper) - gd.fxo) / gd.dxo) as f64) + 0.5f64) as i32);
    if (nxof < 0i32) {
        nxof = 0i32;
    }
    nxoe = ((((((xm + gd.xaper) - gd.fxo) / gd.dxo) as f64) + 0.5f64) as i32);
    if (nxoe >= gd.nxo) {
        nxoe = (gd.nxo - 1i32);
    }
    fprintf!(
        &gd.jpfp,
        "This trace: fxo,dxo,nxof,nxoe=%g %g %d %d\n",
        (gd.fxo as f64),
        (gd.dxo as f64),
        nxof,
        nxoe
    );
    for iyo in nyof..=nyoe {
        yo = (gd.fyo + ((iyo as f32) * gd.dyo));
        ys = ((yo - gd.fys) / gd.dys);
        iys = (ys as i32);
        iys = (if ((iys) < (gd.nys - 2i32)) {
            (iys)
        } else {
            (gd.nys - 2i32)
        });
        iys = (if ((0i32) > (iys)) { (0i32) } else { (iys) });
        ybeta = (ys - (iys as f32));
        ybeta0 = ((1f64 - (ybeta as f64)) as f32);
        yt = (((gd.sy - gd.fygd) - (ybeta * gd.dys)) / gd.dygd);
        iyts0 = (yt as i32);
        iyts0 = (if ((iyts0) < (gd.nygd - 2i32)) {
            (iyts0)
        } else {
            (gd.nygd - 2i32)
        });
        iyts0 = (if ((iyts0) > (0i32)) { (iyts0) } else { (0i32) });
        ys0b = (yt - (iyts0 as f32));
        ys0b0 = ((1f64 - (ys0b as f64)) as f32);
        iyts0 = ((iyts0 as i32) - ((gd.iysf + (gd.iysr * iys)) - (gd.nyt / 2i32))) as i32;
        iyts0 = (if ((iyts0) < (gd.nyt - 2i32)) {
            (iyts0)
        } else {
            (gd.nyt - 2i32)
        });
        iyts0 = (if ((iyts0) > (0i32)) { (iyts0) } else { (0i32) });
        yt = (((gd.gy - gd.fygd) - (ybeta * gd.dys)) / gd.dygd);
        iytg0 = (yt as i32);
        iytg0 = (if ((iytg0) < (gd.nygd - 2i32)) {
            (iytg0)
        } else {
            (gd.nygd - 2i32)
        });
        iytg0 = (if ((iytg0) > (0i32)) { (iytg0) } else { (0i32) });
        yg0b = (yt - (iytg0 as f32));
        yg0b0 = ((1f64 - (yg0b as f64)) as f32);
        iytg0 = ((iytg0 as i32) - ((gd.iysf + (gd.iysr * iys)) - (gd.nyt / 2i32))) as i32;
        iytg0 = (if ((iytg0) < (gd.nyt - 2i32)) {
            (iytg0)
        } else {
            (gd.nyt - 2i32)
        });
        iytg0 = (if ((iytg0) > (0i32)) { (iytg0) } else { (0i32) });
        yt = (((gd.sy - gd.fygd) + (ybeta0 * gd.dys)) / gd.dygd);
        iyts1 = (yt as i32);
        iyts1 = (if ((iyts1) < (gd.nygd - 2i32)) {
            (iyts1)
        } else {
            (gd.nygd - 2i32)
        });
        iyts1 = (if ((iyts1) > (0i32)) { (iyts1) } else { (0i32) });
        ys1b = (yt - (iyts1 as f32));
        ys1b0 = ((1f64 - (ys1b as f64)) as f32);
        iyts1 = ((iyts1 as i32) - ((gd.iysf + (gd.iysr * (iys + 1i32))) - (gd.nyt / 2i32))) as i32;
        iyts1 = (if ((iyts1) < (gd.nyt - 2i32)) {
            (iyts1)
        } else {
            (gd.nyt - 2i32)
        });
        iyts1 = (if ((iyts1) > (0i32)) { (iyts1) } else { (0i32) });
        yt = (((gd.gy - gd.fygd) + (ybeta0 * gd.dys)) / gd.dygd);
        iytg1 = (yt as i32);
        iytg1 = (if ((iytg1) < (gd.nygd - 2i32)) {
            (iytg1)
        } else {
            (gd.nygd - 2i32)
        });
        iytg1 = (if ((iytg1) > (0i32)) { (iytg1) } else { (0i32) });
        yg1b = (yt - (iytg1 as f32));
        yg1b0 = ((1f64 - (yg1b as f64)) as f32);
        iytg1 = ((iytg1 as i32) - ((gd.iysf + (gd.iysr * (iys + 1i32))) - (gd.nyt / 2i32))) as i32;
        iytg1 = (if ((iytg1) < (gd.nyt - 2i32)) {
            (iytg1)
        } else {
            (gd.nyt - 2i32)
        });
        iytg1 = (if ((iytg1) > (0i32)) { (iytg1) } else { (0i32) });
        for ixo in nxof..=nxoe {
            xo = (gd.fxo + ((ixo as f32) * gd.dxo));
            xs = ((xo - gd.fxs) / gd.dxs);
            ixs = (xs as i32);
            ixs = (if ((ixs) < (gd.nxs - 2i32)) {
                (ixs)
            } else {
                (gd.nxs - 2i32)
            });
            ixs = (if ((0i32) > (ixs)) { (0i32) } else { (ixs) });
            xalpha = (xs - (ixs as f32));
            xalpha0 = ((1f64 - (xalpha as f64)) as f32);
            rhos[[iyo, ixo]] = (((((xo - gd.sx) * (xo - gd.sx)) + ((yo - gd.sy) * (yo - gd.sy)))
                as f64)
                .sqrt() as f32);
            rhog[[iyo, ixo]] = (((((xo - gd.gx) * (xo - gd.gx)) + ((yo - gd.gy) * (yo - gd.gy)))
                as f64)
                .sqrt() as f32);
            ddxs = (((((xo - gd.sx) as f64).abs() * (gd.dxo as f64)) / (rhos[[iyo, ixo]] as f64))
                as f32);
            ddys = (((((yo - gd.sy) as f64).abs() * (gd.dyo as f64)) / (rhog[[iyo, ixo]] as f64))
                as f32);
            ddxg = (((((xo - gd.gx) as f64).abs() * (gd.dxo as f64)) / (rhos[[iyo, ixo]] as f64))
                as f32);
            ddyg = (((((yo - gd.gy) as f64).abs() * (gd.dyo as f64)) / (rhog[[iyo, ixo]] as f64))
                as f32);
            drho = ((0.5f64
                * (((((ddxs * ddxs) + (ddys * ddys)) + (ddxg * ddxg)) + (ddyg * ddyg)) as f64))
                .sqrt() as f32);
            temp = (((((4f64 / (gd.v0 as f64)) / (gd.v0 as f64)) * (drho as f64)) / (gd.dt as f64))
                as f32);
            rhos[[iyo, ixo]] = ((rhos[[iyo, ixo]] as f32) * temp) as f32;
            rhog[[iyo, ixo]] = ((rhog[[iyo, ixo]] as f32) * temp) as f32;
            xt = (((gd.sx - gd.fxgd) - (xalpha * gd.dxs)) / gd.dxgd);
            ixts0 = (xt as i32);
            ixts0 = (if ((ixts0) < (gd.nxgd - 2i32)) {
                (ixts0)
            } else {
                (gd.nxgd - 2i32)
            });
            ixts0 = (if ((ixts0) > (0i32)) { (ixts0) } else { (0i32) });
            xs0a = (xt - (ixts0 as f32));
            xs0a0 = ((1f64 - (xs0a as f64)) as f32);
            ixts0 = ((ixts0 as i32) - ((gd.ixsf + (gd.ixsr * ixs)) - (gd.nxt / 2i32))) as i32;
            ixts0 = (if ((ixts0) < (gd.nxt - 2i32)) {
                (ixts0)
            } else {
                (gd.nxt - 2i32)
            });
            ixts0 = (if ((ixts0) > (0i32)) { (ixts0) } else { (0i32) });
            xt = (((gd.gx - gd.fxgd) - (xalpha * gd.dxs)) / gd.dxgd);
            ixtg0 = (xt as i32);
            ixtg0 = (if ((ixtg0) < (gd.nxgd - 2i32)) {
                (ixtg0)
            } else {
                (gd.nxgd - 2i32)
            });
            ixtg0 = (if ((ixtg0) > (0i32)) { (ixtg0) } else { (0i32) });
            xg0a = (xt - (ixtg0 as f32));
            xg0a0 = ((1f64 - (xg0a as f64)) as f32);
            ixtg0 = ((ixtg0 as i32) - ((gd.ixsf + (gd.ixsr * ixs)) - (gd.nxt / 2i32))) as i32;
            ixtg0 = (if ((ixtg0) < (gd.nxt - 2i32)) {
                (ixtg0)
            } else {
                (gd.nxt - 2i32)
            });
            ixtg0 = (if ((ixtg0) > (0i32)) { (ixtg0) } else { (0i32) });
            xt = (((gd.sx - gd.fxgd) + (xalpha0 * gd.dxs)) / gd.dxgd);
            ixts1 = (xt as i32);
            ixts1 = (if ((ixts1) < (gd.nxgd - 2i32)) {
                (ixts1)
            } else {
                (gd.nxgd - 2i32)
            });
            ixts1 = (if ((ixts1) > (0i32)) { (ixts1) } else { (0i32) });
            xs1a = (xt - (ixts1 as f32));
            xs1a0 = ((1f64 - (xs1a as f64)) as f32);
            ixts1 =
                ((ixts1 as i32) - ((gd.ixsf + (gd.ixsr * (ixs + 1i32))) - (gd.nxt / 2i32))) as i32;
            ixts1 = (if ((ixts1) < (gd.nxt - 2i32)) {
                (ixts1)
            } else {
                (gd.nxt - 2i32)
            });
            ixts1 = (if ((ixts1) > (0i32)) { (ixts1) } else { (0i32) });
            xt = (((gd.gx - gd.fxgd) + (xalpha0 * gd.dxs)) / gd.dxgd);
            ixtg1 = (xt as i32);
            ixtg1 = (if ((ixtg1) < (gd.nxgd - 2i32)) {
                (ixtg1)
            } else {
                (gd.nxgd - 2i32)
            });
            ixtg1 = (if ((ixtg1) > (0i32)) { (ixtg1) } else { (0i32) });
            xg1a = (xt - (ixtg1 as f32));
            xg1a0 = ((1f64 - (xg1a as f64)) as f32);
            ixtg1 =
                ((ixtg1 as i32) - ((gd.ixsf + (gd.ixsr * (ixs + 1i32))) - (gd.nxt / 2i32))) as i32;
            ixtg1 = (if ((ixtg1) < (gd.nxt - 2i32)) {
                (ixtg1)
            } else {
                (gd.nxt - 2i32)
            });
            ixtg1 = (if ((ixtg1) > (0i32)) { (ixtg1) } else { (0i32) });
            for izo in 0i32..gd.nzo {
                // ts and tg start at zero and persist between image points,
                // as in the original triangle antialiasing operator.
                ntri = ((if (((1i32) as f32)
                    > (((rhos[[iyo, ixo]] / ts) + (rhog[[iyo, ixo]] / tg)) - (1i32 as f32)))
                {
                    ((1i32) as f32)
                } else {
                    (((rhos[[iyo, ixo]] / ts) + (rhog[[iyo, ixo]] / tg)) - (1i32 as f32))
                }) as i32);
                zo = (gd.fzo + ((izo as f32) * gd.dzo));
                if (gd.crfp != (None)) {
                    ciys0 = ((yg0b0
                        * ((xg0a0 * ((cbuf[[iys, ixo, izo, iytg0, ixtg0]] as i32) as f32))
                            + (xg0a
                                * ((cbuf[[iys, ixo, izo, iytg0, (ixtg0 + 1i32)]] as i32)
                                    as f32))))
                        + (yg0b
                            * ((xg0a0
                                * ((cbuf[[iys, ixo, izo, (iytg0 + 1i32), ixtg0]] as i32)
                                    as f32))
                                + (xg0a
                                    * ((cbuf[[iys, ixo, izo, (iytg0 + 1i32), (ixtg0 + 1i32)]]
                                        as i32) as f32)))));
                    ciys1 = ((yg1b0
                        * ((xg1a0
                            * ((cbuf[[(iys + 1i32), ixo, izo, iytg1, ixtg1]] as i32) as f32))
                            + (xg1a
                                * ((cbuf[[(iys + 1i32), ixo, izo, iytg1, (ixtg1 + 1i32)]] as i32)
                                    as f32))))
                        + (yg1b
                            * ((xg1a0
                                * ((cbuf[[(iys + 1i32), ixo, izo, (iytg1 + 1i32), ixtg1]] as i32)
                                    as f32))
                                + (xg1a
                                    * ((cbuf
                                        [[(iys + 1i32), ixo, izo, (iytg1 + 1i32), (ixtg1 + 1i32)]]
                                        as i32) as f32)))));
                    cosine = (((ybeta0 * ciys0) + (ybeta * ciys1)) / (255i32 as f32));
                    riys0 = ((ys0b0
                        * ((xs0a0 * ((rbuf[[iys, ixo, izo, iyts0, ixts0]] as i32) as f32))
                            + (xs0a
                                * ((rbuf[[iys, ixo, izo, iyts0, (ixts0 + 1i32)]] as i32)
                                    as f32))))
                        + (ys0b
                            * ((xs0a0
                                * ((rbuf[[iys, ixo, izo, (iyts0 + 1i32), ixts0]] as i32)
                                    as f32))
                                + (xs0a
                                    * ((rbuf[[iys, ixo, izo, (iyts0 + 1i32), (ixts0 + 1i32)]]
                                        as i32) as f32)))));
                    riyg0 = ((yg0b0
                        * ((xg0a0 * ((rbuf[[iys, ixo, izo, iytg0, ixtg0]] as i32) as f32))
                            + (xg0a
                                * ((rbuf[[iys, ixo, izo, iytg0, (ixtg0 + 1i32)]] as i32)
                                    as f32))))
                        + (yg0b
                            * ((xg0a0
                                * ((rbuf[[iys, ixo, izo, (iytg0 + 1i32), ixtg0]] as i32)
                                    as f32))
                                + (xg0a
                                    * ((rbuf[[iys, ixo, izo, (iytg0 + 1i32), (ixtg0 + 1i32)]]
                                        as i32) as f32)))));
                    riys1 = ((ys1b0
                        * ((xs1a0
                            * ((rbuf[[(iys + 1i32), ixo, izo, iyts1, ixts1]] as i32) as f32))
                            + (xs1a
                                * ((rbuf[[(iys + 1i32), ixo, izo, iyts1, (ixts1 + 1i32)]] as i32)
                                    as f32))))
                        + (ys1b
                            * ((xs1a0
                                * ((rbuf[[(iys + 1i32), ixo, izo, (iyts1 + 1i32), ixts1]] as i32)
                                    as f32))
                                + (xs1a
                                    * ((rbuf
                                        [[(iys + 1i32), ixo, izo, (iyts1 + 1i32), (ixts1 + 1i32)]]
                                        as i32) as f32)))));
                    riyg1 = ((yg1b0
                        * ((xg1a0
                            * ((rbuf[[(iys + 1i32), ixo, izo, iytg1, ixtg1]] as i32) as f32))
                            + (xg1a
                                * ((rbuf[[(iys + 1i32), ixo, izo, iytg1, (ixtg1 + 1i32)]] as i32)
                                    as f32))))
                        + (yg1b
                            * ((xg1a0
                                * ((rbuf[[(iys + 1i32), ixo, izo, (iytg1 + 1i32), ixtg1]] as i32)
                                    as f32))
                                + (xg1a
                                    * ((rbuf
                                        [[(iys + 1i32), ixo, izo, (iytg1 + 1i32), (ixtg1 + 1i32)]]
                                        as i32) as f32)))));
                    rs =
                        ((((ybeta0 * riys0) + (ybeta * riys1)) / (255i32 as f32)) * (10i32 as f32));
                    rg =
                        ((((ybeta0 * riyg0) + (ybeta * riyg1)) / (255i32 as f32)) * (10i32 as f32));
                    rtotal = (rs + rg);
                } else {
                    cosine = (1f64 as f32);
                    rs = ((((((xo - gd.sx) * (xo - gd.sx)) + ((yo - gd.sy) * (yo - gd.sy)))
                        + (zo * zo)) as f64)
                        .sqrt() as f32);
                    rg = ((((((xo - gd.gx) * (xo - gd.gx)) + ((yo - gd.gy) * (yo - gd.gy)))
                        + (zo * zo)) as f64)
                        .sqrt() as f32);
                    rtotal = (rs + rg);
                    if (((((((rs + rg) as f64) * 1000f64) - (ttotal as f64)).abs()
                        / (((rs + rg) as f64) * 1000f64))
                        > 0.050000000000000003f64)
                        && (izo == 0i32))
                    {
                        fprintf!(
                            &gd.jpfp,
                            "at: ixo=%d,iyo=%d,izo=%d: est. ttot=%f, err=%f\n",
                            ixo,
                            iyo,
                            izo,
                            (ttotal as f64),
                            (((((rs + rg) as f64) * 1000f64) - (ttotal as f64))
                                / (((rs + rg) as f64) * 1000f64))
                        );
                    }
                }
                amp = (((((cosine * zo) * rtotal) as f64)
                    * (((1f64 / (rs as f64)) / (rs as f64)) + ((1f64 / (rg as f64)) / (rg as f64))))
                    as f32);
                for im in 0i32..gd.multit {
                    if (((((((((((((((((tbuf[[im, iys, ixo, izo, iyts0, ixts0]]
                        as i32)
                        >= 65535i32)
                        || ((tbuf[[im, iys, ixo, izo, (iyts0 + 1i32), ixts0]]
                            as i32)
                            >= 65535i32))
                        || ((tbuf[[im, iys, ixo, izo, iyts0, (ixts0 + 1i32)]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, iys, ixo, izo, (iyts0 + 1i32), (ixts0 + 1i32)]]
                            as i32)
                            >= 65535i32))
                        || ((tbuf[[im, iys, ixo, izo, iytg0, ixtg0]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, iys, ixo, izo, (iytg0 + 1i32), ixtg0]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, iys, ixo, izo, iytg0, (ixtg0 + 1i32)]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, iys, ixo, izo, (iytg0 + 1i32), (ixtg0 + 1i32)]]
                            as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, iyts1, ixts1]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, (iyts1 + 1i32), ixts1]]
                            as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, iyts1, (ixts1 + 1i32)]]
                            as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, (iyts1 + 1i32), (ixts1 + 1i32)]]
                            as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, iytg1, ixtg1]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, (iytg1 + 1i32), ixtg1]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, iytg1, (ixtg1 + 1i32)]] as i32)
                            >= 65535i32))
                        || ((tbuf[[im, (iys + 1i32), ixo, izo, (iytg1 + 1i32), (ixtg1 + 1i32)]]
                            as i32)
                            >= 65535i32))
                    {
                        fprintf!(&gd.jpfp, " skipped ixo=%d ", ixo);
                        continue;
                    }
                    tiys0 = ((ys0b0
                        * ((xs0a0 * ((tbuf[[im, iys, ixo, izo, iyts0, ixts0]] as i32) as f32))
                            + (xs0a
                                * ((tbuf[[im, iys, ixo, izo, iyts0, (ixts0 + 1i32)]] as i32)
                                    as f32))))
                        + (ys0b
                            * ((xs0a0
                                * ((tbuf[[im, iys, ixo, izo, (iyts0 + 1i32), ixts0]] as i32)
                                    as f32))
                                + (xs0a
                                    * ((tbuf[[im, iys, ixo, izo, (iyts0 + 1i32), (ixts0 + 1i32)]]
                                        as i32) as f32)))));
                    tiyg0 = ((yg0b0
                        * ((xg0a0 * ((tbuf[[im, iys, ixo, izo, iytg0, ixtg0]] as i32) as f32))
                            + (xg0a
                                * ((tbuf[[im, iys, ixo, izo, iytg0, (ixtg0 + 1i32)]] as i32)
                                    as f32))))
                        + (yg0b
                            * ((xg0a0
                                * ((tbuf[[im, iys, ixo, izo, (iytg0 + 1i32), ixtg0]] as i32)
                                    as f32))
                                + (xg0a
                                    * ((tbuf[[im, iys, ixo, izo, (iytg0 + 1i32), (ixtg0 + 1i32)]]
                                        as i32) as f32)))));
                    tiys1 = ((ys1b0
                        * ((xs1a0
                            * ((tbuf[[im, (iys + 1i32), ixo, izo, iyts1, ixts1]] as i32)
                                as f32))
                            + (xs1a
                                * ((tbuf[[im, (iys + 1i32), ixo, izo, iyts1, (ixts1 + 1i32)]]
                                    as i32) as f32))))
                        + (ys1b
                            * ((xs1a0
                                * ((tbuf[[im, (iys + 1i32), ixo, izo, (iyts1 + 1i32), ixts1]]
                                    as i32) as f32))
                                + (xs1a
                                    * ((tbuf[[
                                        im,
                                        (iys + 1i32),
                                        ixo,
                                        izo,
                                        (iyts1 + 1i32),
                                        (ixts1 + 1i32),
                                    ]] as i32) as f32)))));
                    tiyg1 = ((yg1b0
                        * ((xg1a0
                            * ((tbuf[[im, (iys + 1i32), ixo, izo, iytg1, ixtg1]] as i32)
                                as f32))
                            + (xg1a
                                * ((tbuf[[im, (iys + 1i32), ixo, izo, iytg1, (ixtg1 + 1i32)]]
                                    as i32) as f32))))
                        + (yg1b
                            * ((xg1a0
                                * ((tbuf[[im, (iys + 1i32), ixo, izo, (iytg1 + 1i32), ixtg1]]
                                    as i32) as f32))
                                + (xg1a
                                    * ((tbuf[[
                                        im,
                                        (iys + 1i32),
                                        ixo,
                                        izo,
                                        (iytg1 + 1i32),
                                        (ixtg1 + 1i32),
                                    ]] as i32) as f32)))));
                    ts = ((((ybeta0 * tiys0) + (ybeta * tiys1)) / (65535i32 as f32))
                        * (7000i32 as f32));
                    tg = ((((ybeta0 * tiyg0) + (ybeta * tiyg1)) / (65535i32 as f32))
                        * (7000i32 as f32));
                    ts = ((if (((ts) as f64) > (0.001f64)) {
                        ((ts) as f64)
                    } else {
                        (0.001f64)
                    }) as f32);
                    tg = ((if (((tg) as f64) > (0.001f64)) {
                        ((tg) as f64)
                    } else {
                        (0.001f64)
                    }) as f32);
                    ttotal = (ts + tg);
                    if ((ttotal >= (7000i32 as f32)) || ((ttotal as f64) <= 0f64)) {
                        continue;
                    }
                    tintr = ((ttotal - gd.ft) / gd.dt);
                    itintr = (tintr as i32);
                    if (((itintr - ntri) > 0i32) && (((itintr + ntri) + 1i32) < gd.nt)) {
                        wtintr = (tintr - (itintr as f32));
                        wtintr0 = ((1f64 - (wtintr as f64)) as f32);
                        // The selected C mig[iyoffset][ixoffset] volume is
                        // contiguous, with depth as its innermost dimension.
                        let p = (iyo as usize * gd.nxo as usize + ixo as usize) * gd.nzo as usize
                            + izo as usize;
                        mig[p] = ((mig[p] as f64)
                            + ((amp as f64)
                                * (((wtintr0 as f64)
                                    * ((((-trf[(itintr - ntri) as usize]) as f64)
                                        + (2f64 * (trf[(itintr) as usize] as f64)))
                                        - (trf[(itintr + ntri) as usize] as f64)))
                                    + ((wtintr as f64)
                                        * ((((-trf[((itintr - ntri) + 1i32) as usize]) as f64)
                                            + (2f64 * (trf[(itintr + 1i32) as usize] as f64)))
                                            - (trf[((itintr + ntri) + 1i32) as usize] as f64))))))
                            as f32;
                    }
                }
            }
        }
    }
}

// get_bufs - interpolate and quantize the traveltime buffer.
// Narrow float-to-byte/short casts go through u32, matching the original
// ARM64 fcvtzu instruction followed by truncation of the low bits.
fn get_bufs(
    gd: &Gd,
    ttab: &Grid<f32, 6>,
    ctab: &Grid<f32, 5>,
    rtab: &Grid<f32, 5>,
    tbuf: &mut Grid<u16, 6>,
    cbuf: &mut Grid<u8, 5>,
    rbuf: &mut Grid<u8, 5>,
) {
    let mut iyt: i32 = 0;
    let mut ixt: i32 = 0;
    let mut iys: i32 = 0;
    let mut ixs: i32 = 0;
    let mut izs: i32 = 0;
    let mut izo: i32 = 0;
    let mut ixo: i32 = 0;
    let mut im: i32 = 0;
    let mut xalpha: f32 = 0.0;
    let mut xalpha0: f32 = 0.0;
    let mut zgamma: f32 = 0.0;
    let mut zgamma0: f32 = 0.0;
    let mut xo: f32 = 0.0;
    let mut xs: f32 = 0.0;
    let mut zo: f32 = 0.0;
    let mut zs: f32 = 0.0;
    let mut t: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut r: f32 = 0.0;
    for ixo in 0i32..gd.nxo {
        xo = (gd.fxo + ((ixo as f32) * gd.dxo));
        xs = ((xo - gd.fxs) / gd.dxs);
        ixs = (xs as i32);
        ixs = (if ((ixs) < (gd.nxs - 2i32)) {
            (ixs)
        } else {
            (gd.nxs - 2i32)
        });
        ixs = (if ((ixs) > (0i32)) { (ixs) } else { (0i32) });
        xalpha = (xs - (ixs as f32));
        xalpha0 = ((1f64 - (xalpha as f64)) as f32);
        for izo in 0i32..gd.nzo {
            zo = (gd.fzo + ((izo as f32) * gd.dzo));
            zs = ((zo - gd.fzs) / gd.dzs);
            izs = (zs as i32);
            izs = (if ((izs) < (gd.nzs - 2i32)) {
                (izs)
            } else {
                (gd.nzs - 2i32)
            });
            izs = (if ((izs) > (0i32)) { (izs) } else { (0i32) });
            zgamma = (zs - (izs as f32));
            zgamma0 = ((1f64 - (zgamma as f64)) as f32);
            for iys in 0i32..gd.nys {
                for iyt in 0i32..gd.nyt {
                    for ixt in 0i32..gd.nxt {
                        if ((((ttab[[0i32, iys, ixs, izs, iyt, ixt]] >= (7000i32 as f32))
                            || (ttab[[0i32, iys, (ixs + 1i32), izs, iyt, ixt]]
                                >= (7000i32 as f32)))
                            || (ttab[[0i32, iys, ixs, (izs + 1i32), iyt, ixt]]
                                >= (7000i32 as f32)))
                            || (ttab[[0i32, iys, (ixs + 1i32), (izs + 1i32), iyt, ixt]]
                                >= (7000i32 as f32)))
                        {
                            tbuf[[0i32, iys, ixo, izo, iyt, ixt]] = (65535i32 as u16);
                            cbuf[[iys, ixo, izo, iyt, ixt]] = (0i32 as u8);
                            rbuf[[iys, ixo, izo, iyt, ixt]] = (255i32 as u8);
                        } else {
                            t = ((zgamma0
                                * ((xalpha0 * ttab[[0i32, iys, ixs, izs, iyt, ixt]])
                                    + (xalpha * ttab[[0i32, iys, (ixs + 1i32), izs, iyt, ixt]])))
                                + (zgamma
                                    * ((xalpha0
                                        * ttab[[0i32, iys, ixs, (izs + 1i32), iyt, ixt]])
                                        + (xalpha
                                            * ttab[[
                                                0i32,
                                                iys,
                                                (ixs + 1i32),
                                                (izs + 1i32),
                                                iyt,
                                                ixt,
                                            ]]))));
                            tbuf[[0i32, iys, ixo, izo, iyt, ixt]] =
                                (((t / (7000i32 as f32)) * (65535i32 as f32)) as u32 as u16);
                            if (gd.crfp == (None)) {
                                continue;
                            }
                            c = ((zgamma0
                                * ((xalpha0 * ctab[[iys, ixs, izs, iyt, ixt]])
                                    + (xalpha * ctab[[iys, (ixs + 1i32), izs, iyt, ixt]])))
                                + (zgamma
                                    * ((xalpha0 * ctab[[iys, ixs, (izs + 1i32), iyt, ixt]])
                                        + (xalpha
                                            * ctab[[iys, (ixs + 1i32), (izs + 1i32), iyt, ixt]]))));
                            c = ((if (((c) as f64) < (1f64)) {
                                ((c) as f64)
                            } else {
                                (1f64)
                            }) as f32);
                            cbuf[[iys, ixo, izo, iyt, ixt]] = ((c * (255i32 as f32)) as u32 as u8);
                            r = ((zgamma0
                                * ((xalpha0 * rtab[[iys, ixs, izs, iyt, ixt]])
                                    + (xalpha * rtab[[iys, (ixs + 1i32), izs, iyt, ixt]])))
                                + (zgamma
                                    * ((xalpha0 * rtab[[iys, ixs, (izs + 1i32), iyt, ixt]])
                                        + (xalpha
                                            * rtab[[iys, (ixs + 1i32), (izs + 1i32), iyt, ixt]]))));
                            r = (if ((r) < ((10i32) as f32)) {
                                (r)
                            } else {
                                ((10i32) as f32)
                            });
                            rbuf[[iys, ixo, izo, iyt, ixt]] =
                                (((r / (10i32 as f32)) * (255i32 as f32)) as u32 as u8);
                        }
                        for im in 1i32..gd.multit {
                            if ((((ttab[[im, iys, ixs, izs, iyt, ixt]] >= (7000i32 as f32))
                                || (ttab[[im, iys, (ixs + 1i32), izs, iyt, ixt]]
                                    >= (7000i32 as f32)))
                                || (ttab[[im, iys, ixs, (izs + 1i32), iyt, ixt]]
                                    >= (7000i32 as f32)))
                                || (ttab[[im, iys, (ixs + 1i32), (izs + 1i32), iyt, ixt]]
                                    >= (7000i32 as f32)))
                            {
                                tbuf[[im, iys, ixo, izo, iyt, ixt]] = (65535i32 as u16);
                            } else {
                                t = ((zgamma0
                                    * ((xalpha0 * ttab[[im, iys, ixs, izs, iyt, ixt]])
                                        + (xalpha
                                            * ttab[[im, iys, (ixs + 1i32), izs, iyt, ixt]])))
                                    + (zgamma
                                        * ((xalpha0
                                            * ttab[[im, iys, ixs, (izs + 1i32), iyt, ixt]])
                                            + (xalpha
                                                * ttab[[
                                                    im,
                                                    iys,
                                                    (ixs + 1i32),
                                                    (izs + 1i32),
                                                    iyt,
                                                    ixt,
                                                ]]))));
                                tbuf[[im, iys, ixo, izo, iyt, ixt]] =
                                    (((t / (7000i32 as f32)) * (65535i32 as f32)) as u32 as u16);
                            }
                        }
                    }
                }
            }
        }
    }
}

// Flat storage follows the C allocations. In particular, tbuf retains the
// original ealloc6ushort(multit,nxt,nyt,nzo,nxo,nys) order. C pointer rows
// can overlap when indexed beyond an inner dimension; preserve those offsets.
struct Grid<T, const N: usize> {
    dims: [i32; N],
    data: Vec<T>,
}
impl<T: Default + Clone, const N: usize> Grid<T, N> {
    fn new(dims: [i32; N]) -> Self {
        let mut size = 1usize;
        for &n in &dims {
            // Negative sizes become enormous size_t arguments in C and fail
            // in ealloc. Check before converting or multiplying in Rust.
            if n < 0 {
                err!("ealloc.c: malloc failed (Cannot allocate memory)");
            }
            // Empty output grids are valid: repairing C's zero-size pointer
            // tables does not change their output (UB audit, November 2026).
            size = size
                .checked_mul(n as usize)
                .unwrap_or_else(|| err!("ealloc.c: malloc failed (Cannot allocate memory)"));
        }
        let mut data = Vec::new();
        // Also checks the element byte count and Rust's isize::MAX limit.
        if data.try_reserve_exact(size).is_err() {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        data.resize(size, T::default());
        Self { dims, data }
    }
}
impl<T, const N: usize> Grid<T, N> {
    fn offset(&self, indices: [i32; N]) -> usize {
        let mut offset = 0usize;
        let mut count = 1usize;
        for k in 0..N {
            count *= self.dims[k] as usize;
            offset = offset
                .wrapping_mul(self.dims[k] as usize)
                .wrapping_add(indices[k] as usize);
            if offset >= count {
                err!("invalid grid index");
            }
        }
        offset
    }
}
impl<T, const N: usize> std::ops::Index<[i32; N]> for Grid<T, N> {
    type Output = T;
    fn index(&self, i: [i32; N]) -> &T {
        &self.data[self.offset(i)]
    }
}
impl<T, const N: usize> std::ops::IndexMut<[i32; N]> for Grid<T, N> {
    fn index_mut(&mut self, i: [i32; N]) -> &mut T {
        let p = self.offset(i);
        &mut self.data[p]
    }
}
// filt - Low-pass filter, integration and phase shift for input data.
// Static C filter state is owned by main and shared across traces.
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
    fn filt(&mut self, trace: &mut [f32], gd: &mut Gd, trf: &mut [f32]) {
        // The original multiplies fmax on every trace, even after setup.
        gd.fmax = ((gd.fmax as f64) * (2.0 * PI)) as f32;
        let ftaper = (0.1 * gd.fmax as f64) as f32;
        let const2 = 2.0f64.sqrt() as f32;
        if self.nfft == 0 {
            let nmin = gd
                .nt
                .checked_add(gd.ntrimax)
                .unwrap_or_else(|| err!("invalid filter dimensions"));
            let nmax = nmin
                .checked_mul(2)
                .unwrap_or_else(|| err!("invalid filter dimensions"));
            self.nfft = npfaro(nmin, nmax);
            if self.nfft >= SU_NFLTS as i32 || self.nfft >= 720720 {
                err!("Padded nt=%d -- too big", self.nfft);
            }
            self.nw = self.nfft / 2 + 1;
            self.dw = (2.0 * PI / (self.nfft as f32 * gd.dt) as f64) as f32;
            self.itaper = (0.5 + (ftaper / self.dw) as f64) as i32;
            let ntaper = self
                .itaper
                .checked_mul(2)
                .and_then(|n| n.checked_add(1))
                .unwrap_or_else(|| err!("ealloc.c: malloc failed (Cannot allocate memory)"));
            self.taper = Grid::<f32, 1>::new([ntaper]).data;
            for iw in -self.itaper..=self.itaper {
                let temp = (iw as f32 as f64 / (1.0 + self.itaper as f64)) as f32;
                self.taper[(iw + self.itaper) as usize] =
                    (((1.0f32 - temp) * (1.0f32 - temp) * (temp + 2.0f32)) as f64 / 4.0) as f32;
            }
            self.nwf = (0.5 + (gd.fmax / self.dw) as f64) as i32;
            if self.nwf > self.nw - self.itaper - 1 {
                self.nwf = self.nw - self.itaper - 1;
            }
            let namp = self
                .nwf
                .checked_add(self.itaper)
                .and_then(|n| n.checked_add(1))
                .unwrap_or_else(|| err!("ealloc.c: malloc failed (Cannot allocate memory)"));
            self.amp = Grid::<f32, 1>::new([namp]).data;
            self.ampi = Grid::<f32, 1>::new([namp]).data;
            for iw in 1..=self.nwf + self.itaper {
                // sqrt returns double; division precedes narrowing to float.
                self.amp[iw as usize] =
                    (((self.dw * iw as f32) as f64).sqrt() / self.nfft as f64) as f32;
                self.ampi[iw as usize] =
                    (0.5 / (1.0 - ((iw as f32 * self.dw * gd.dt) as f64).cos())) as f32;
            }
        }
        let mut rt = Grid::<f32, 1>::new([self.nfft]).data;
        let mut ct = Grid::<Complex, 1>::new([self.nw]).data;
        rt[..gd.nt as usize].copy_from_slice(&trace[..gd.nt as usize]);
        pfarc(1, &rt, &mut ct);
        for iw in self.nwf - self.itaper..=self.nwf + self.itaper {
            if iw < 0 || iw >= self.nw {
                err!("invalid filter index");
            }
            let itemp = iw - (self.nwf - self.itaper);
            ct[iw as usize].r = self.taper[itemp as usize] * ct[iw as usize].r;
            ct[iw as usize].i = self.taper[itemp as usize] * ct[iw as usize].i;
        }
        for iw in self.nwf + self.itaper + 1..self.nw {
            ct[iw as usize] = Complex::default();
        }
        if gd.ls == 0 {
            for iw in 0..=self.nwf + self.itaper {
                let z = &mut ct[iw as usize];
                let temp = (z.r - z.i) * self.amp[iw as usize] * const2;
                z.i = (z.r + z.i) * self.amp[iw as usize] * const2;
                z.r = temp;
            }
        } else {
            for iw in 0..=self.nwf + self.itaper {
                ct[iw as usize].i *= self.amp[iw as usize];
                ct[iw as usize].r *= self.amp[iw as usize];
            }
        }
        pfacr(-1, &ct, &mut rt);
        trace[..gd.nt as usize].copy_from_slice(&rt[..gd.nt as usize]);
        for iw in 0..=self.nwf + self.itaper {
            ct[iw as usize].i *= self.ampi[iw as usize];
            ct[iw as usize].r *= self.ampi[iw as usize];
        }
        pfacr(-1, &ct, &mut rt);
        for it in 0..gd.ntrimax {
            trf[it as usize] = rt[(self.nfft - gd.ntrimax + it) as usize];
        }
        for it in 0..gd.nt + gd.ntrimax {
            trf[(it + gd.ntrimax) as usize] = rt[it as usize];
        }
    }
}

fn main() {
    let argv = args();
    let mut gd = Gd::new();
    let mut tr = Trace::new();
    let mut ttr = Trace::new();
    let mut nxoffset: i32 = 0;
    let mut nyoffset: i32 = 0;
    let mut ixoffset: i32 = 0;
    let mut iyoffset: i32 = 0;
    let mut ixo: i32 = 0;
    let mut iyo: i32 = 0;
    let mut ixt: i32 = 0;
    let mut iyt: i32 = 0;
    let mut izs: i32 = 0;
    let mut ntrmax: i32 = 0;
    let mut imigtr: i32 = 0;
    let mut itotmigtr: i32 = 0;
    let mut pptr: i32 = 0;
    let mut fxoffset: f32 = 0.0;
    let mut fyoffset: f32 = 0.0;
    let mut dxm: f32 = 0.0;
    let mut angmax: f32 = 0.0;
    let mut dxoffset: f32 = 0.0;
    let mut dyoffset: f32 = 0.0;
    let mut ext: f32 = 0.0;
    let mut eyt: f32 = 0.0;
    let mut ezs: f32 = 0.0;
    let mut ezo: f32 = 0.0;
    let mut exo: f32 = 0.0;
    let mut eyo: f32 = 0.0;
    let mut exs: f32 = 0.0;
    let mut eys: f32 = 0.0;
    let mut xoffsetmax: f32 = 0.0;
    let mut yoffsetmax: f32 = 0.0;
    let mut iys: i32 = 0;
    let mut ixs: i32 = 0;
    let mut imul: i32 = 0;
    let mut xm: f32 = 0.0;
    let mut ym: f32 = 0.0;
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut datain = String::new();
    let infp = if getparstring("datain", &mut datain) == 0 {
        stdin()
    } else {
        fopen(&datain, "r").unwrap_or_else(|| err!("cannot open datain=%s", datain))
    };
    let mut dataout = String::new();
    let outfp = if getparstring("dataout", &mut dataout) == 0 {
        Some(stdout())
    } else {
        fopen(&dataout, "w")
    };
    efseeko(&infp, 0, SEEK_CUR);
    let outfp = outfp.unwrap_or_else(|| err!("cannot open dataout=%s", dataout));
    efseeko(&outfp, 0, SEEK_CUR);
    let mut ttfile = String::new();
    if getparstring("ttfile", &mut ttfile) == 0 {
        err!("must specify ttfile!");
    }
    let mut ttfp = fopen(&ttfile, "r").unwrap_or_else(|| err!("Can not open ttfile=%s", ttfile));
    if fgettr(&ttfp, &mut ttr) == 0 {
        err!("Can't get first trace from tttable");
    }
    let mut crfile = String::new();
    if getparstring("crfile", &mut crfile) != 0 {
        gd.crfp =
            Some(fopen(&crfile, "r").unwrap_or_else(|| err!("Can not open crfile=%s", crfile)));
    }
    if (getparfloat("fxgd", &mut gd.fxgd) == 0) {
        gd.fxgd = ttr.f1;
    }
    if (getparfloat("dxgd", &mut gd.dxgd) == 0) {
        gd.dxgd = ttr.d1;
    }
    if (getparint("nxt", &mut gd.nxt) == 0) {
        gd.nxt = (ttr.ns as i32);
    }
    if ((gd.dxgd as f64) < 0.10000000000000001f64) {
        gd.v0 = (2.5f64 as f32);
    } else {
        gd.v0 = (2500f64 as f32);
    }
    if (getparint("multit", &mut gd.multit) == 0) {
        gd.multit = (ttr.scalel as i32);
    }
    if (gd.multit == 0i32) {
        gd.multit = 1i32;
    }
    gd.nxt = gd
        .nxt
        .checked_div(gd.multit)
        .unwrap_or_else(|| err!("invalid table dimensions"));
    if (getparint("ls", &mut gd.ls) == 0) {
        gd.ls = 0i32;
    }
    if (getparfloat("fygd", &mut gd.fygd) == 0) {
        gd.fygd = ttr.f2;
    }
    if (getparfloat("dygd", &mut gd.dygd) == 0) {
        gd.dygd = ttr.d2;
    }
    if (getparint("nyt", &mut gd.nyt) == 0) {
        gd.nyt = ttr.ntr;
    }
    gd.nxgd = ttr.selev;
    gd.nygd = ttr.gelev;
    if (getparint("nxs", &mut gd.nxs) == 0) {
        gd.nxs = (ttr.nhs as i32);
    }
    if (getparint("nys", &mut gd.nys) == 0) {
        gd.nys = (ttr.nvs as i32);
    }
    if (getparint("nzs", &mut gd.nzs) == 0) {
        gd.nzs = (ttr.duse as i32);
    }
    if (getparint("ixsf", &mut gd.ixsf) == 0) {
        gd.ixsf = ttr.sdel;
    }
    if (getparint("iysf", &mut gd.iysf) == 0) {
        gd.iysf = ttr.gdel;
    }
    if (getparfloat("fzs", &mut gd.fzs) == 0) {
        gd.fzs = (((ttr.sdepth as f64) / 1000f64) as f32);
    }
    if (getparint("ixsr", &mut gd.ixsr) == 0) {
        gd.ixsr = ttr.swdep;
    }
    if (getparint("iysr", &mut gd.iysr) == 0) {
        gd.iysr = ttr.gwdep;
    }
    if (getparfloat("dzs", &mut gd.dzs) == 0) {
        gd.dzs = (((ttr.ep as f64) / 1000f64) as f32);
    }
    fclose(&ttfp);
    ttfp = fopen(&ttfile, "r").unwrap_or_else(|| err!("Can not open ttfile=%s", ttfile));
    efseeko(&ttfp, 0, SEEK_CUR);
    if (getparfloat("dyoffset", &mut dyoffset) == 0) {
        dyoffset = (99999i32 as f32);
    }
    if (getparint("ntrmax", &mut ntrmax) == 0) {
        ntrmax = 100000i32;
    }
    if (getparint("pptr", &mut pptr) == 0) {
        pptr = 100i32;
    }
    gd.fxs = (gd.fxgd + ((gd.ixsf as f32) * gd.dxgd));
    gd.dxs = ((gd.ixsr as f32) * gd.dxgd);
    gd.fys = (gd.fygd + ((gd.iysf as f32) * gd.dygd));
    gd.dys = ((gd.iysr as f32) * gd.dygd);
    exs = (gd.fxs + (((gd.nxs - 1i32) as f32) * gd.dxs));
    eys = (gd.fys + (((gd.nys - 1i32) as f32) * gd.dys));
    ezs = (gd.fzs + (((gd.nzs - 1i32) as f32) * gd.dzs));
    ext = (gd.fxgd + (((gd.nxgd - 1i32) as f32) * gd.dxgd));
    eyt = (gd.fygd + (((gd.nygd - 1i32) as f32) * gd.dygd));
    if (getparint("nxo", &mut gd.nxo) == 0) {
        gd.nxo = (((gd.nxgd - 1i32) * 2i32) + 1i32);
    }
    if (getparfloat("fxo", &mut gd.fxo) == 0) {
        gd.fxo = gd.fxgd;
    }
    if (getparfloat("dxo", &mut gd.dxo) == 0) {
        gd.dxo = (((gd.dxgd as f64) * 0.5f64) as f32);
    }
    if (getparint("nyo", &mut gd.nyo) == 0) {
        gd.nyo = (((gd.nygd - 1i32) * 2i32) + 1i32);
    }
    if (getparfloat("fyo", &mut gd.fyo) == 0) {
        gd.fyo = gd.fygd;
    }
    if (getparfloat("dyo", &mut gd.dyo) == 0) {
        gd.dyo = (((gd.dygd as f64) * 0.5f64) as f32);
    }
    if (getparint("nzo", &mut gd.nzo) == 0) {
        gd.nzo = (((gd.nzs - 1i32) * 5i32) + 1i32);
    }
    if (getparfloat("fzo", &mut gd.fzo) == 0) {
        gd.fzo = gd.fzs;
    }
    if (getparfloat("dzo", &mut gd.dzo) == 0) {
        gd.dzo = (((gd.dzs as f64) * 0.20000000000000001f64) as f32);
    }
    exo = (gd.fxo + (((gd.nxo - 1i32) as f32) * gd.dxo));
    eyo = (gd.fyo + (((gd.nyo - 1i32) as f32) * gd.dyo));
    ezo = (gd.fzo + (((gd.nzo - 1i32) as f32) * gd.dzo));
    if (getparfloat("xaper", &mut gd.xaper) == 0) {
        gd.xaper = (((((gd.nxt as f32) * gd.dxgd) as f64) / 2.5f64) as f32);
    }
    if (getparfloat("yaper", &mut gd.yaper) == 0) {
        gd.yaper = (((((gd.nyt as f32) * gd.dygd) as f64) / 2.5f64) as f32);
    }
    if (getparfloat("xoffsetmax", &mut xoffsetmax) == 0) {
        xoffsetmax = (3000i32 as f32);
    }
    if (getparfloat("yoffsetmax", &mut yoffsetmax) == 0) {
        yoffsetmax = (3000i32 as f32);
    }
    if (getparint("nxoffset", &mut nxoffset) == 0) {
        nxoffset = 1i32;
    }
    if (getparint("nyoffset", &mut nyoffset) == 0) {
        nyoffset = 1i32;
    }
    if (getparfloat("fxoffset", &mut fxoffset) == 0) {
        fxoffset = (0f64 as f32);
    }
    if (getparfloat("fyoffset", &mut fyoffset) == 0) {
        fyoffset = (0f64 as f32);
    }
    if (getparfloat("dxoffset", &mut dxoffset) == 0) {
        dxoffset = (99999i32 as f32);
    }
    if fgettr(&infp, &mut tr) == 0 {
        err!("Can't get first trace");
    }
    gd.nt = (tr.ns as i32);
    if (getparfloat("dt", &mut gd.dt) == 0) {
        gd.dt = (((tr.dt as f64) / 1000f64) as f32);
    }
    if ((gd.dt as f64) < 9.9999999999999995E-8f64) {
        err!("dt must be positive!");
    }
    if (getparfloat("ft", &mut gd.ft) == 0) {
        gd.ft = ((((tr.delrt as i32) as f64) / 1000f64) as f32);
    }
    if (getparfloat("dxm", &mut dxm) == 0) {
        dxm = tr.d2;
    }
    if ((dxm as f64) < 9.9999999999999995E-8f64) {
        err!("dxm must be positive");
    }
    if (getparfloat("angmax", &mut angmax) == 0) {
        angmax = (((60f64 * (3.1415926535897931f64)) / 180f64) as f32);
    }
    if ((angmax as f64) < 9.9999999999999995E-8f64) {
        err!("angmax must be positive");
    }
    gd.ntrimax = ((((((2i32 as f32) * dxm) as f64) * (angmax as f64).sin())
        / ((gd.v0 * gd.dt) as f64)) as i32);
    if (gd.ntrimax < 1i32) {
        gd.ntrimax = 1i32;
    }
    if (getparfloat("fmax", &mut gd.fmax) == 0) {
        gd.fmax = ((0.25f64 / (gd.dt as f64)) as f32);
    }
    let mut jpfile = String::new();
    let jpfp = if getparstring("jpfile", &mut jpfile) == 0 {
        Some(stderr())
    } else {
        fopen(&jpfile, "w")
    };
    checkpars();
    gd.jpfp = jpfp.unwrap_or_else(|| err!("cannot open jpfile=%s", jpfile));
    fprintf!(&gd.jpfp, "\n");
    fprintf!(&gd.jpfp, " Migration parameters\n");
    fprintf!(&gd.jpfp, " ================\n");
    fprintf!(&gd.jpfp, " ttfile=%s \n", ttfile);
    fprintf!(&gd.jpfp, " \n");
    fprintf!(
        &gd.jpfp,
        " nxt=%d fxgd=%g dxgd=%g\n",
        gd.nxt,
        (gd.fxgd as f64),
        (gd.dxgd as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nyt=%d fygd=%g dygd=%g\n",
        gd.nyt,
        (gd.fygd as f64),
        (gd.dygd as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nxs=%d fxs=%g dxs=%g\n",
        gd.nxs,
        (gd.fxs as f64),
        (gd.dxs as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nys=%d fys=%g dys=%g\n",
        gd.nys,
        (gd.fys as f64),
        (gd.dys as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nzs=%d gd->fzs=%g dzs=%g\n",
        gd.nzs,
        (gd.fzs as f64),
        (gd.dzs as f64)
    );
    fprintf!(&gd.jpfp, " nxgd=%d nygd=%d\n\n", gd.nxgd, gd.nygd);
    fprintf!(
        &gd.jpfp,
        " nzo=%d fzo=%g dzo=%g\n",
        gd.nzo,
        (gd.fzo as f64),
        (gd.dzo as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nxo=%d fxo=%g dxo=%g\n",
        gd.nxo,
        (gd.fxo as f64),
        (gd.dxo as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nyo=%d fyo=%g dyo=%g\n\n",
        gd.nyo,
        (gd.fyo as f64),
        (gd.dyo as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nt=%d ft=%g dt=%g dxm=%g\n\n",
        gd.nt,
        (gd.ft as f64),
        (gd.dt as f64),
        (dxm as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nxoffset=%d fxoffset=%g dxoffset=%g\n",
        nxoffset,
        (fxoffset as f64),
        (dxoffset as f64)
    );
    fprintf!(
        &gd.jpfp,
        " nyoffset=%d fyoffset=%g dyoffset=%g\n",
        nyoffset,
        (fyoffset as f64),
        (dyoffset as f64)
    );
    fprintf!(
        &gd.jpfp,
        " xaper=%g xoffsetmax=%g\n",
        (gd.xaper as f64),
        (xoffsetmax as f64)
    );
    fprintf!(
        &gd.jpfp,
        " yaper=%g yoffsetmax=%g\n",
        (gd.yaper as f64),
        (yoffsetmax as f64)
    );
    fprintf!(&gd.jpfp, " ntrmax=%d pptr=%d\n", ntrmax, pptr);
    fprintf!(&gd.jpfp, " ext=%g exo=%g\n", (ext as f64), (exo as f64));
    fprintf!(&gd.jpfp, " eyt=%g eyo=%g\n", (eyt as f64), (eyo as f64));
    fprintf!(&gd.jpfp, " ezs=%g ezo=%g\n", (ezs as f64), (ezo as f64));
    fprintf!(&gd.jpfp, " ====================================\n");
    fflush(&gd.jpfp);
    if ((((gd.fxgd > gd.fxo) || (ext < exo)) || (gd.fygd > gd.fyo)) || (eyt < eyo)) {
        fprintf!(
            &gd.jpfp,
            "fxgd=%g > fxo=%g ?\n",
            (gd.fxgd as f64),
            (gd.fxo as f64)
        );
        fprintf!(&gd.jpfp, "ext=%g < exo=%g ?\n", (ext as f64), (exo as f64));
        fprintf!(
            &gd.jpfp,
            "fygd=%g > fyo=%g ?\n",
            (gd.fygd as f64),
            (gd.fyo as f64)
        );
        fprintf!(&gd.jpfp, "eyt=%g < eyo=%g ?\n", (eyt as f64), (eyo as f64));
        err!(" migration output range is out of tttable");
    }
    fprintf!(
        &gd.jpfp,
        "allocating mig of size %g Mbytes\n",
        4.0e-6 * gd.nzo as f64 * gd.nxo as f64 * gd.nyo as f64 * nxoffset as f64 * nyoffset as f64
    );
    let mut mig = Grid::<f32, 5>::new([nyoffset, nxoffset, gd.nyo, gd.nxo, gd.nzo]);
    fprintf!(
        &gd.jpfp,
        "allocating ttab, ctab and rtab of size %g Mbytes\n",
        12.0e-6 * gd.nzs as f64 * gd.nxs as f64 * gd.nys as f64 * gd.nxt as f64 * gd.nyt as f64
    );
    let mut ttab = Grid::<f32, 6>::new([gd.multit, gd.nys, gd.nxs, gd.nzs, gd.nyt, gd.nxt]);
    let mut ctab = Grid::<f32, 5>::new([gd.nys, gd.nxs, gd.nzs, gd.nyt, gd.nxt]);
    let mut rtab = Grid::<f32, 5>::new([gd.nys, gd.nxs, gd.nzs, gd.nyt, gd.nxt]);
    fprintf!(
        &gd.jpfp,
        "allocating tbuf, cbuf and rbuf of size %g Mbytes\n",
        4.0e-6 * gd.nys as f64 * gd.nxo as f64 * gd.nzo as f64 * gd.nyt as f64 * gd.nxt as f64
    );
    let mut tbuf = Grid::<u16, 6>::new([gd.nys, gd.nxo, gd.nzo, gd.nyt, gd.nxt, gd.multit]);
    let mut cbuf = Grid::<u8, 5>::new([gd.nys, gd.nxo, gd.nzo, gd.nyt, gd.nxt]);
    let mut rbuf = Grid::<u8, 5>::new([gd.nys, gd.nxo, gd.nzo, gd.nyt, gd.nxt]);
    let ntrf = gd
        .ntrimax
        .checked_mul(2)
        .and_then(|n| gd.nt.checked_add(n))
        .unwrap_or_else(|| err!("ealloc.c: malloc failed (Cannot allocate memory)"));
    let mut trf = Grid::<f32, 1>::new([ntrf]).data;
    fprintf!(&gd.jpfp, "\nReading ttables ...");
    fprintf!(&gd.jpfp, "\n");
    for iys in 0..gd.nys {
        for ixs in 0..gd.nxs {
            for izs in 0..gd.nzs {
                fprintf!(&gd.jpfp, "Reading source No. (%d,%d,%d)\n", ixs, iys, izs);
                for iyt in 0..gd.nyt {
                    if fgettr(&ttfp, &mut ttr) == 0 {
                        err!("Can't read ttfile");
                    }
                    for imul in 0..gd.multit {
                        for ixt in 0..gd.nxt {
                            ttab[[imul, iys, ixs, izs, iyt, ixt]] =
                                ttr.data[(ixt + imul * gd.nxt) as usize];
                        }
                    }
                    if let Some(crfp) = &gd.crfp {
                        // C reads both rows, then overwrites each sample with 1.0.
                        let mut row = Grid::<f32, 1>::new([gd.nxt]).data;
                        if fread(&mut row, crfp) != gd.nxt as usize {
                            err!("Can not read in crfile");
                        }
                        for ixt in 0..gd.nxt {
                            ctab[[iys, ixs, izs, iyt, ixt]] = 1.0;
                        }
                        if fread(&mut row, crfp) != gd.nxt as usize {
                            err!("Can not read in crfile");
                        }
                        for ixt in 0..gd.nxt {
                            rtab[[iys, ixs, izs, iyt, ixt]] = 1.0;
                        }
                    }
                }
            }
        }
    }
    fclose(&ttfp);
    fprintf!(&gd.jpfp, "\nTable reading done\n");
    fprintf!(&gd.jpfp, "\nGet traveltime buffer\n");
    get_bufs(&gd, &ttab, &ctab, &rtab, &mut tbuf, &mut cbuf, &mut rbuf);
    fprintf!(&gd.jpfp, "Start migration ...\n");
    fflush(&gd.jpfp);
    let mut filter = Filter::default();
    imigtr = 1;
    itotmigtr = 0;
    loop {
        gd.sx = tr.sx as f32;
        gd.gx = tr.gx as f32;
        gd.sy = tr.sy as f32;
        gd.gy = tr.gy as f32;
        // DEBUG is enabled in the original file.
        gd.sy = 3.0;
        gd.gy = 3.0;
        fprintf!(
            &gd.jpfp,
            "trace %d: sx,gx,sy,gy=%g %g %g %g\n",
            imigtr,
            gd.sx,
            gd.gx,
            gd.sy,
            gd.gy
        );
        fprintf!(
            &gd.jpfp,
            "trace in line: %d; trace in reel: %d\n",
            tr.tracl,
            tr.tracr
        );
        xm = (0.5 * (gd.sx + gd.gx) as f64) as f32;
        ym = (0.5 * (gd.sy + gd.gy) as f64) as f32;
        let xdiff = if gd.gx - gd.sx > gd.sx - gd.gx {
            gd.gx - gd.sx
        } else {
            gd.sx - gd.gx
        };
        let ydiff = if gd.gy - gd.sy > gd.sy - gd.gy {
            gd.gy - gd.sy
        } else {
            gd.sy - gd.gy
        };
        // C continue in a do/while still reads the next trace, and does not
        // increment imigtr when the geometry rejects a trace.
        if !(xm < gd.fxs
            || xm > exs
            || xdiff > xoffsetmax
            || ym < gd.fys
            || ym > eys
            || ydiff > yoffsetmax)
        {
            ixoffset = (((gd.gx - gd.sx - fxoffset) / dxoffset) as f64 + 0.5) as i32;
            if ixoffset < 0 {
                ixoffset = 0;
            }
            if ixoffset >= nxoffset {
                ixoffset = nxoffset - 1;
            }
            iyoffset = (((gd.gy - gd.sy - fyoffset) / dyoffset) as f64 + 0.5) as i32;
            if iyoffset < 0 {
                iyoffset = 0;
            }
            if iyoffset >= nyoffset {
                iyoffset = nyoffset - 1;
            }
            filter.filt(&mut tr.data, &mut gd, &mut trf);
            // Grid::new already checked the full product. Compute this
            // subvolume length in usize so it cannot overflow a C int.
            let len = gd.nyo as usize * gd.nxo as usize * gd.nzo as usize;
            // C passes a pointer to an empty subvolume when any output
            // dimension is zero. No element is read or written in that case.
            let start = if len == 0 {
                0
            } else {
                mig.offset([iyoffset, ixoffset, 0, 0, 0])
            };
            mig3d(
                &gd,
                &trf,
                &mut mig.data[start..start + len],
                &tbuf,
                &cbuf,
                &rbuf,
            );
            itotmigtr += 1;
            // ARM64 integer division by zero returns zero.
            if (imigtr - 1).checked_rem(pptr).unwrap_or(imigtr - 1) == 0 {
                fprintf!(&gd.jpfp, "Migrated trace %d\n", imigtr);
                fflush(&gd.jpfp);
            }
            imigtr += 1;
        }
        if fgettr(&infp, &mut tr) == 0 || imigtr >= ntrmax {
            break;
        }
    }
    fprintf!(&gd.jpfp, "Migrated %d traces in total\n", itotmigtr);
    let mut tro = Trace::new();
    tro.ns = gd.nzo as u16;
    tro.dt = 1000i32.wrapping_mul(((1000.0f32 * gd.dt) as f64 + 0.5) as i32) as u16;
    tro.delrt = 0;
    tro.f1 = gd.fzo;
    tro.d1 = gd.dzo;
    tro.f2 = gd.fxo;
    tro.d2 = gd.dxo;
    for iyo in 0..gd.nyo {
        for ixo in 0..gd.nxo {
            for iyoffset in 0..nyoffset {
                for ixoffset in 0..nxoffset {
                    // A zero-sample trace still emits its SU header.
                    let start = if gd.nzo == 0 {
                        0
                    } else {
                        mig.offset([iyoffset, ixoffset, iyo, ixo, 0])
                    };
                    tro.data[..gd.nzo as usize]
                        .copy_from_slice(&mig.data[start..start + gd.nzo as usize]);
                    tro.offset = (fxoffset + ixoffset as f32 * dxoffset) as i32;
                    tro.tracr = 1 + ixo;
                    tro.tracl = 1 + ixoffset + nxoffset * ixo;
                    tro.cdp = (gd.fxo + ixo as f32 * gd.dxo) as i32;
                    tro.cdpt = 1 + ixoffset;
                    fputtr(&outfp, &mut tro);
                }
            }
        }
    }
    fprintf!(&gd.jpfp, "\n");
    fprintf!(&gd.jpfp, "Output done\n");
    fflush(&gd.jpfp);
    efclose(&gd.jpfp);
    efclose(&outfp);
    exit(CWP_Exit());
}
