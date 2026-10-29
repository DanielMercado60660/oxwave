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
//! SUBACKUSH: port of well_logs/subackush.c, revision 1.6.
//! Chris Liner acknowledges Saudi Aramco's permission to release this code
//! developed at EXPEC-ARC. References: Thomsen's 2002 SEG DISC Notes;
//! Berryman, Grechka and Berge (1997), and Wang (2002).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUBACKUSH - calculate Thomsen anisotropy parameters from \t",
    " \t     well log (vp,vs,rho) data and optionally include\t",
    " \t     intrinsic VTI shale layers based on gramma ray log\t",
    " \t     via BACKUS averaging\t\t\t\t",
    " subackush < vp_vs_rho.su >stdout [options]\t\t\t",
    " subackush < vp_vs_rho_gr.su  gr=1 >stdout [options]\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    " navg=101\tnumber of depth samples in Backus avg window \t",
    "\t\t\t\t\t\t\t\t",
    " \tIntrinsic anisotropy of shale layers can be included ...",
    " gr=0\t\tno gamma ray log input for shale \t\t",
    "\t\t=1 input is vp_vs_rho_gr\t\t\t",
    " grs=100\tpure shale gamma ray value (API units)\t\t",
    " grc=10\t0% shale gamma ray value (API units)\t\t",
    " smode=1\tinclude shale anis params prop to shale volume \t",
    "\t\t=0 include shale anis for pure shale only\t",
    " se=0.209\tshale epsilon (Thomsen parameter)\t\t",
    " sd=0.033\tshale delta (Thomsen parameter)\t\t\t",
    " sg=0.203\tshale gamma (Thomsen parameter)\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " 1. Input are (vp,vs,rho) traces in metric units\t\t",
    " 2. Output are  \t\t\t\t\t\t",
    "    tracl\t=(1,2,3,4,5,6)\t\t\t\t\t",
    "    quantity\t=(vp0,vs0,<rho>,epsilon,delta,gamma) \t\t",
    "    units\t=(m/s,m/s,kg/m^3,nd,nd,nd) nd=dimensionless\t",
    "    tracl\t=(7,8)\t\t\t\t\t\t",
    "    quantity\t=(Vsh,shaleEps) Vsh=shale volume fraction\t",
    "    units\t=(nd,nd) \t\t\t\t\t",
    " 3. (epsilon,delta,etc.) can be isolated by tracl header field ",
    " 4. (vp0,vs0) are backus averaged vertical wavespeeds\t\t",
    " 5. <rho> is backus averaged density, etc.\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t",
    " las2su < logs.las nskip=34 nlog=4 > logs.su \t\t\t",
    " suwind < logs.su  key=tracl min=2 max=3 | suop op=s2vm > v.su\t",
    " suwind < logs.su  key=tracl min=4 max=4 | suop op=d2m > d.su\t",
    " fcat v.su d.su > vp_vs_rho.su\t\t\t\t\t",
    " subackus < vp_vs_rho.su > vp0_vs0_rho_eps_delta_gamma.su\t",
    " In this example we start with a well las file containing \t",
    " 34 header lines and 4 log tracks (depth,p_son,s_son,den).\t",
    " This is converted to su format by las2su.  Then we pull off\t",
    " the sonic logs and convert them to velocity in metric units.\t",
    " Then the density log is pulled off and converted to metric.\t",
    " All three metric curves are bundled into one su file which \t",
    " is the input to subackus. \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Related programs: subackus, sulprime\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn dobackus(navg: i32, nz: usize, p: &[f32], avg: &mut [f32]) {
    let m = (navg - 1) / 2;
    for i in 0..nz as i32 {
        let mut val = 0.0f32;
        if i - m > 0 && i + m < nz as i32 {
            for j in i - m..=i + m {
                val += p[j as usize];
            }
            avg[i as usize] = val / navg as f32;
        }
        if i <= m {
            avg[i as usize] = p[i as usize];
        }
        if i >= nz as i32 - m {
            avg[i as usize] = p[i as usize];
        }
    }
}

fn handl_ends(navg: i32, nz: usize, p: &mut [f32]) {
    let m = (navg - 1) / 2 + 1;
    // The original reads outside its allocation if the window is too large.
    if m < 0 || m as usize > nz {
        return;
    }
    let gt = p[m as usize];
    let gb = p[nz - m as usize];
    for i in 0..m as usize {
        p[i] = gt;
    }
    for i in nz - m as usize..nz {
        p[i] = gb;
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut navg = 101;
    let mut gr = 0;
    let mut grs = 100.0f32;
    let mut grc = 10.0f32;
    let mut se = 0.209f32;
    let mut sd = 0.033f32;
    let mut sg = 0.203f32;
    let mut smode = 1;
    getparint("navg", &mut navg);
    getparint("gr", &mut gr);
    getparfloat("grs", &mut grs);
    getparfloat("grc", &mut grc);
    getparfloat("se", &mut se);
    getparfloat("sd", &mut sd);
    getparfloat("sg", &mut sg);
    getparint("smode", &mut smode);
    checkpars();
    if navg != 201 {
        navg = navg / 2 * 2 + 1;
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nz = tr.ns as usize;
    let d1 = tr.d1;
    let f1 = tr.f1;
    let dt = tr.dt;
    let vp = tr.data[..nz].to_vec();
    gettr(&mut tr);
    let vs = tr.data[..nz].to_vec();
    gettr(&mut tr);
    let rhoo = tr.data[..nz].to_vec();
    let mut vsh = vec![0.0f32; nz];
    if gr == 1 {
        gettr(&mut tr);
        for i in 0..nz {
            let ftmp = (tr.data[i] - grc) / (grs - grc);
            vsh[i] = (0.33 * (2.0f64.powf((2.0f32 * ftmp) as f64) - 1.0)) as f32;
            if vsh[i] > 1.0 {
                vsh[i] = 1.0;
            }
            if vsh[i] < 0.0 {
                vsh[i] = 0.0;
            }
        }
    }
    let mut ao = vec![0.0f32; nz];
    let mut fo = ao.clone();
    let mut co = ao.clone();
    let mut lo = ao.clone();
    let mut mo = ao.clone();
    for i in 0..nz {
        let (mut epso, mut deltao, mut gammao) = (0.0f32, 0.0f32, 0.0f32);
        if gr == 1 {
            if smode == 0 && vsh[i] >= 1.0 {
                epso = se;
                deltao = sd;
                gammao = sg;
            }
            if smode == 1 {
                epso = se * vsh[i];
                deltao = sd * vsh[i];
                gammao = sg * vsh[i];
            }
        }
        ao[i] = ((1.0 + 2.0 * epso as f64) * rhoo[i] as f64 * vp[i] as f64 * vp[i] as f64) as f32;
        let mut ftmp = vp[i] * vp[i] - vs[i] * vs[i];
        ftmp = (ftmp as f64
            * ((1.0 + 2.0 * deltao as f64) * vp[i] as f64 * vp[i] as f64 - (vs[i] * vs[i]) as f64))
            as f32;
        fo[i] = (rhoo[i] as f64 * ((ftmp as f64).sqrt() - (vs[i] * vs[i]) as f64)) as f32;
        co[i] = rhoo[i] * vp[i] * vp[i];
        lo[i] = rhoo[i] * vs[i] * vs[i];
        mo[i] = ((1.0 + 2.0 * gammao as f64) * rhoo[i] as f64 * vs[i] as f64 * vs[i] as f64) as f32;
    }
    let mut tmp = vec![0.0f32; nz];
    let mut btmp = tmp.clone();
    let mut tmp3 = tmp.clone();
    let mut ab = tmp.clone();
    let mut fb = tmp.clone();
    let mut cb = tmp.clone();
    let mut lb = tmp.clone();
    let mut mb = tmp.clone();
    let mut rhob = tmp.clone();
    for i in 0..nz {
        tmp[i] = (1.0 / co[i] as f64) as f32;
    }
    dobackus(navg, nz, &tmp, &mut cb);
    for i in 0..nz {
        cb[i] = (1.0 / cb[i] as f64) as f32;
    }
    for i in 0..nz {
        tmp[i] = fo[i] / co[i];
    }
    dobackus(navg, nz, &tmp, &mut btmp);
    for i in 0..nz {
        fb[i] = cb[i] * btmp[i];
        tmp3[i] = btmp[i] * btmp[i];
    }
    for i in 0..nz {
        tmp[i] = ao[i] - fo[i] * fo[i] / co[i];
    }
    dobackus(navg, nz, &tmp, &mut btmp);
    for i in 0..nz {
        ab[i] = btmp[i] + cb[i] * tmp3[i];
    }
    for i in 0..nz {
        tmp[i] = (1.0 / lo[i] as f64) as f32;
    }
    dobackus(navg, nz, &tmp, &mut lb);
    for i in 0..nz {
        lb[i] = (1.0 / lb[i] as f64) as f32;
    }
    dobackus(navg, nz, &mo, &mut mb);
    dobackus(navg, nz, &rhoo, &mut rhob);
    let mut vp0 = tmp.clone();
    let mut vs0 = tmp.clone();
    let mut epsb = tmp.clone();
    let mut deltab = tmp.clone();
    let mut gammab = tmp.clone();
    for i in 0..nz {
        vp0[i] = ((cb[i] / rhob[i]) as f64).sqrt() as f32;
        vs0[i] = ((lb[i] / rhob[i]) as f64).sqrt() as f32;
        epsb[i] = ((ab[i] - cb[i]) as f64 / (2.0 * cb[i] as f64)) as f32;
        let ftmp1 = (fb[i] + lb[i]) * (fb[i] + lb[i]);
        let ftmp2 = (cb[i] - lb[i]) * (cb[i] - lb[i]);
        let ftmp3 = (2.0 * cb[i] as f64 * (cb[i] - lb[i]) as f64) as f32;
        deltab[i] = (ftmp1 - ftmp2) / ftmp3;
        gammab[i] = ((mb[i] - lb[i]) as f64 / (2.0 * lb[i] as f64)) as f32;
    }
    tr.trid = 1;
    tr.ns = nz as u16;
    tr.dt = dt;
    tr.d1 = d1;
    tr.f1 = f1;
    for (k, values) in [&vp0, &vs0, &rhob, &epsb, &deltab, &gammab]
        .iter()
        .enumerate()
    {
        tr.tracl = (k + 1) as i32;
        tr.data[..nz].copy_from_slice(values);
        handl_ends(navg, nz, &mut tr.data);
        puttr(&mut tr);
    }
    tr.tracl = 7;
    tr.data[..nz].copy_from_slice(&vsh);
    puttr(&mut tr);
    tr.tracl = 8;
    for i in 0..nz {
        tr.data[i] = if smode == 1 && vsh[i] < 1.0 {
            se * vsh[i]
        } else {
            se
        };
    }
    puttr(&mut tr);
    exit(CWP_Exit());
}
