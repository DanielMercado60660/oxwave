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
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SUBACKUS: well_logs/subackus.c, revision 1.6.
//! Chris Liner acknowledges Saudi Aramco permission to release this code
//! developed at EXPEC-ARC. References: Thomsen, 2002, SEG DISC Notes;
//! Berryman, Grechka, and Berge, 1997, SEP94.
use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUBACKUS - calculate Thomsen anisotropy parameters from \t",
    " \t     well log (vp,vs,rho) data via BACKUS averaging\t",
    " \t\t\t\t\t\t\t\t",
    " subackus < vp_vs_rho.su >stdout [options]\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameter:\t\t\t\t\t\t",
    " navg=201\tnumber of depth samples in Backus avg window \t",
    " all=0\t\t=1 to output extra parameters \t\t\t",
    "\t\t(<vp0>,<vs0>,<rho>,eta,vang,a,f,c,l,A,B,<lam>,<mu>)",
    " ang=30\tangle (deg) for use in vang\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " 1. Input are (vp,vs,rho) traces in metric units\t\t",
    " 2. Output are (epsilon,delta,gamma) dimensionless traces\t",
    "    tracl=(1,2,3)=(epsilon,delta,gamma) \t\t\t",
    "\t(all=1 output optional traces listed below)\t\t",
    "    tracl=(4,5,6,7,8)=(vp0,vs0,<rho>,eta,vang)\t\t\t",
    "    tracl=(9,10,11,12)=(a,f,c,l)=(c11,c13,c33,c44) backus avg'd",
    "    tracl=(13,14)=(<lam/lamp2mu>^2,4<mu*lampmu/lamp2mu>)=(A,B)\t",
    "       used to analyze eps=(a-c)/2c; a=c11=A*x+B;  c=c33=x\t",
    "    tracl=(15,16)=(<lambda>,<mu>)\t\t\t\t",
    "       for fluid analysis (lambda affected by fluid, mu not)   ",
    "    tracl=(17,18,19)=(vp,vs,rho)  orig log values\t\t",
    "    tracl=(20)=(m)=(c66) Backus avg'd \t\t\t\t",
    "    tracl=(21,22,23,24,25)=(a,f,c,l,m)=(c11,c13,c33,c44,c66) orig",
    " 3. (epsilon,delta,etc.) can be isolated by tracl header field ",
    " 4. (vp0,vs0) are backus averaged vertical wavespeeds\t\t",
    " 5. <rho> is backus averaged density, etc.\t\t\t",
    " 6. eta = (eps - delta) / (1 + 2 delta)\t\t\t",
    " 7. vang=vp(ang_deg)/vp0  phase velocity ratio\t\t\t",
    "    The idea being that if vang~vp0 there are small time effects",
    "    (30 deg comes from ~ max ray angle preserved in processing)",
    "\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t",
    " las2su < logs.las nskip=34 nlog=4 > logs.su \t\t\t",
    " suwind < logs.su  key=tracl min=2 max=3 | suop op=s2vm > v.su\t",
    " suwind < logs.su  key=tracl min=4 max=4 | suop op=d2m > d.su\t",
    " fcat v.su d.su > vp_vs_rho.su\t\t\t\t\t",
    " subackus < vp_vs_rho.su > eps_delta_gamma.su\t\t\t",
    " In this example we start with a well las file containing \t",
    " 34 header lines and 4 log tracks (depth,p_son,s_son,den).\t",
    " This is converted to su format by las2su.  Then we pull off\t",
    " the sonic logs and convert them to velocity in metric units.\t",
    " Then the density log is pulled off and converted to metric.\t",
    " All three metric curves are bundled into one su file which \t",
    " is the input to subackus. \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Related codes: sulprime subackush\t\t\t\t",
];

fn dobackus(navg: i32, p: &[f32]) -> Vec<f32> {
    let nz = p.len() as i32;
    let m = (navg - 1) / 2;
    let mut avg = vec![0.0; p.len()];
    for i in 0..nz {
        let mut val = 0.0f32;
        if i - m > 0 && i + m < nz {
            for j in i - m..=i + m { val += p[j as usize]; }
            avg[i as usize] = val / navg as f32;
        }
        if i <= m { avg[i as usize] = p[i as usize]; }
        if i >= nz - m { avg[i as usize] = p[i as usize]; }
    }
    avg
}

fn handl_ends(navg: i32, nz: usize, p: &mut [f32]) {
    let m = (navg - 1) / 2 + 1;
    // Oversized windows in the original access outside the trace data.
    let gt = p.get(m as usize).copied().unwrap_or(0.0);
    let gb = p.get((nz as i32 - m) as usize).copied().unwrap_or(0.0);
    for i in 0..m { if let Some(x) = p.get_mut(i as usize) { *x = gt; } }
    for i in nz as i32 - m..nz as i32 { if let Some(x) = p.get_mut(i as usize) { *x = gb; } }
}

fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    let mut navg = 201; let mut all = 0; let mut ang = 30.0f32;
    getparint("navg", &mut navg); getparint("all", &mut all); getparfloat("ang", &mut ang);
    checkpars();
    if navg != 201 { navg = navg / 2 * 2 + 1; }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nz = tr.ns as usize; let d1 = tr.d1; let f1 = tr.f1; let dt = tr.dt as f32;
    let vp = tr.data[..nz].to_vec(); gettr(&mut tr);
    let vs = tr.data[..nz].to_vec(); gettr(&mut tr);
    let rho = tr.data[..nz].to_vec();
    let mut mu = vec![0.0f32; nz]; let mut lp2mu = mu.clone(); let mut lam = mu.clone(); let mut lpmu = mu.clone();
    for i in 0..nz {
        mu[i] = rho[i]*vs[i]*vs[i]; lp2mu[i] = rho[i]*vp[i]*vp[i];
        lam[i] = rho[i]*(vp[i]*vp[i] - 2.0f32*vs[i]*vs[i]);
        lpmu[i] = rho[i]*(vp[i]*vp[i] - vs[i]*vs[i]);
    }
    let blam = dobackus(navg, &lam); let bmu = dobackus(navg, &mu);
    let mut t1 = mu.clone(); let mut t2 = mu.clone(); let mut t3 = mu.clone();
    for i in 0..nz {
        t1[i] = lam[i]/lp2mu[i]; t2[i] = (1.0 / lp2mu[i] as f64) as f32;
        t3[i] = mu[i]*lpmu[i]/lp2mu[i];
    }
    let mut t4 = dobackus(navg, &t1); let mut t5 = dobackus(navg, &t2); let t6 = dobackus(navg, &t3);
    let mut c11 = mu.clone(); let mut aa = mu.clone(); let mut bb = mu.clone();
    for i in 0..nz {
        t4[i] = t4[i]*t4[i]; t5[i] = (1.0/t5[i] as f64) as f32;
        c11[i] = ((t4[i]*t5[i]) as f64 + 4.0*t6[i] as f64) as f32;
        aa[i] = t4[i]; bb[i] = 4.0f32*t6[i];
    }
    for i in 0..nz { t1[i] = lam[i]/lp2mu[i]; t2[i] = (1.0/lp2mu[i] as f64) as f32; }
    t3 = dobackus(navg, &t1); t4 = dobackus(navg, &t2);
    let mut c13 = mu.clone();
    for i in 0..nz { t4[i] = (1.0/t4[i] as f64) as f32; c13[i] = t3[i]*t4[i]; }
    for i in 0..nz { t1[i] = (1.0/lp2mu[i] as f64) as f32; }
    let mut c33 = dobackus(navg, &t1); for x in &mut c33 { *x = (1.0/(*x as f64)) as f32; }
    for i in 0..nz { t1[i] = (1.0/mu[i] as f64) as f32; }
    let mut c44 = dobackus(navg, &t1); for x in &mut c44 { *x = (1.0/(*x as f64)) as f32; }
    let c66 = dobackus(navg, &mu); let rhoa = dobackus(navg, &rho);
    let sina = ((ang as f64)*3.1415927/180.0).sin() as f32;
    let cosa = ((ang as f64)*3.1415927/180.0).cos() as f32;
    let sina2 = sina*sina; let cosa2 = cosa*cosa; let sina4 = sina2*sina2;
    let mut eps=mu.clone(); let mut delta=mu.clone(); let mut gamma=mu.clone();
    let mut vp0=mu.clone(); let mut vs0=mu.clone(); let mut eta=mu.clone(); let mut vang=mu.clone();
    for i in 0..nz {
        eps[i] = ((c11[i]-c33[i]) as f64/(2.0*c33[i] as f64)) as f32;
        t1[i] = (c13[i]+c44[i])*(c13[i]+c44[i]); t2[i] = (c33[i]-c44[i])*(c33[i]-c44[i]);
        t3[i] = (2.0*c33[i] as f64*(c33[i]-c44[i]) as f64) as f32;
        delta[i] = (t1[i]-t2[i])/t3[i]; gamma[i] = ((c66[i]-c44[i]) as f64/(2.0*c44[i] as f64)) as f32;
        vp0[i] = ((c33[i]/rhoa[i]) as f64).sqrt() as f32; vs0[i] = ((c44[i]/rhoa[i]) as f64).sqrt() as f32;
        eta[i] = ((eps[i]-delta[i]) as f64/(1.0+2.0*delta[i] as f64)) as f32;
        let u = 1.0f32 - (vs0[i]*vs0[i])/(vp0[i]*vp0[i]); let ds = (2.0f32*delta[i]-eps[i])*u;
        let dds = (0.5*(((u*u + 4.0f32*ds*sina2*cosa2 + 4.0f32*(u+eps[i])*eps[i]*sina4) as f64).sqrt()-u as f64)) as f32;
        let vpa = (vp0[i] as f64*((1.0f32+eps[i]*sina2+dds) as f64).sqrt()) as f32; vang[i] = vpa/vp0[i];
    }
    tr.trid=1; tr.ns=nz as u16; tr.dt=dt as u16; tr.d1=d1; tr.f1=f1;
    let outputs: [&[f32];25] = [&eps,&delta,&gamma,&vp0,&vs0,&rhoa,&eta,&vang,&c11,&c13,&c33,&c44,&aa,&bb,&blam,&bmu,&vp,&vs,&rho,&c66,&lp2mu,&lam,&lp2mu,&mu,&mu];
    let count = if all == 1 {25} else {3};
    for (k, data) in outputs.iter().enumerate().take(count) {
        tr.tracl = (k+1) as i32; tr.data[..nz].copy_from_slice(data);
        if k < 16 || k == 19 { handl_ends(navg, nz, &mut tr.data); }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
