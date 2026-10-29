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
//! SUTXTAPER - taper a data panel in time and distance.
//! Port of su/main/tapering/sutxtaper.c (Revision 1.4).
//! Credits: based on sutaper; CWP: Chris Liner, Jack K. Cohen.
//! Rewrite: Tagir Galikeev, October 2002; Gerald Klein, IFM-GEOMAR, April 2004.
use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                               ",
    " SUTXTAPER - TAPER in (X,T) the edges of a data panel to zero.\t",
    "                                                               ",
    " sutxtaper <stdin >stdout [optional parameters]\t\t",
    "                                                               ",
    " Optional Parameters:                                          ",
    " low=0.    \tminimum amplitude factor of taper\t\t",
    " tbeg=0    \tlength of taper (ms) at trace start\t\t",
    " tend=0     \tlength of taper (ms) at trace end\t\t",
    " taper=1       taper type                                      ",
    "                 =1 linear (default)                           ",
    "                 =2 sine                                       ",
    "                 =3 cosine                                     ",
    "                 =4 gaussian (+/-3.8)                          ",
    "                 =5 gaussian (+/-2.0)                          ",
    " key=tr\tset key to compute x-domain taper weights\t",
    "               default is using internal tracecount (tr)       ",
    " tr1=0         number of traces to be tapered at beg (key=tr)\t",
    " tr2=tr1       number of traces to be tapered at end (key=tr)\t",
    " \t\t\t\t\t\t\t\t",
    " min=0.\tminimum value of key where taper starts (amp=1.)",
    " max=0.\tmaximum value of key where taper starts (amp=1.)",
    " dx=1. \tlength of taper (in key units)\t\t\t",
    "\t\tif key=tr (unset) length is tr1 and (ntr-tr2)\t",
    "                                                               ",
    " Notes:                                                        ",
    "   Taper type is used for trace (x-domain) tapering as well \t",
    "   as for time domain tapering.\t\t\t\t",
    "   The taper is applied to all traces <tr1 (or key<min) and    ",
    "   >tr2 (or key >max) and all time samples <tbeg and >tend. \t",
    "   Taper weights are amp*1 for traces n tr1<n<tr2 and computed\t",
    "   for all other traces corresponding to the taper typ.\t",
    "   If key is given the taper length is defined by dx, otherwise",
    "   the length of taper is tr1 and (ntr-tr2) respectively.\t",
    "   To eliminate the taper, choose tbeg=0. and tend=0. and tr1=0",
    "   If key is set, min,max values take precedence over tr1,tr2.\t",
    "                                                               ",
];

fn envelope(f: f32, kind: i32, low: f32, linear_low: bool) -> f32 {
    match kind as i8 {
        1 => if linear_low { low + (1.0f32-low)*f } else { f },
        2 => (std::f64::consts::PI*f as f64/2.0).sin() as f32,
        3 => (0.5*(1.0-(std::f64::consts::PI*f as f64).cos())) as f32,
        4 | 5 => { let x = ((if kind as i8 == 4 { 3.8090232 } else { 2.0 }) * (1.0f32-f) as f64) as f32; (-(x*x) as f64).exp() as f32 },
        _ => err!(" taper ?!"),
    }
}
fn taper(t1: f32, t2: f32, kind: i32, t: f32, dt: f32, data: &mut [f32]) {
    let nt = ((t/dt) as f64+1.0) as i32;
    let nt1 = ((t1/dt) as f64+1.0) as i32;
    let nt2 = ((t2/dt) as f64+1.0) as i32;
    if nt1 > 1 { for i in 0..nt1 { let f=i as f32/nt1 as f32; *data.get_mut(i as usize).unwrap_or_else(|| err!("invalid taper length"))*=envelope(f,kind,0.0,false); } }
    if nt2 > 1 { for i in 0..nt2 { let f=i as f32/nt2 as f32; *data.get_mut((nt-i) as usize).unwrap_or_else(|| err!("invalid taper length"))*=envelope(f,kind,0.0,false); } }
}
fn main() {
    let argv=args(); let mut tr=Trace::new(); initargs(&argv); requestdoc(0,SDOC);
    let mut verbose: i16=0; getparshort("verbose",&mut verbose);
    let mut low=0.0f32; getparfloat("low",&mut low); if low>1.0 { err!("low must be less than 1"); }
    let mut t1=0.0f32; getparfloat("tbeg",&mut t1);
    let mut t2=0.0f32; getparfloat("tend",&mut t2);
    let mut ttaper=1; getparint("taper",&mut ttaper);
    let mut key=String::from("tr"); getparstring("key",&mut key);
    let mut min=0.0f32; getparfloat("min",&mut min);
    let mut max=0.0f32; getparfloat("max",&mut max);
    let mut dkey=1.0f32; getparfloat("dx",&mut dkey);
    let mut tr1=0; getparint("tr1",&mut tr1);
    let mut tr2=tr1; getparint("tr2",&mut tr2);
    if gettr(&mut tr)==0 { err!("can't get first trace"); }
    let nt=tr.ns as usize; let dt=(tr.dt as f64/1000.0) as f32;
    let mut ntr=tr.ntr; let mut have_ntr=ntr!=0;
    if getparint("n2",&mut ntr)!=0 || getparint("ntr",&mut ntr)!=0 { have_ntr=true; }
    if !have_ntr { warn!("ntr neither set nor getparred, \n                  traces will be counted "); }
    let mut index=0; let mut typ="";
    if key!="tr" { typ=hdtype(&key); index=getindex(&key); }
    let tlen=(nt as i32-1) as f32*dt;
    if t1+t2>tlen { err!("sum of tapers tbeg=%f, tend=%f exceeds trace length(%f ms)",t1,t2,tlen); }
    checkpars(); let headerfp=etmpfile(); let tracefp=etmpfile(); let x2fp=etmpfile();
    let mut nx=0i32;
    loop {
        nx=nx.wrapping_add(1);
        let x=if index!=0 { let mut val=Value::default(); gethval(&tr,index,&mut val); vtof(typ,val) } else { nx as f32 };
        efwrite(&[x],&x2fp);
        if t1!=0.0 || t2!=0.0 { taper(t1,t2,ttaper,tlen,dt,&mut tr.data); }
        efwrite(tr.hdr.as_bytes(),&headerfp);
        efwrite(&tr.data[..nt],&tracefp);
        if gettr(&mut tr)==0 { break; }
    }
    if !have_ntr { ntr=nx; }
    if ntr<0 { err!("invalid trace count"); }
    let mut x2=vec![0.0f32;ntr as usize]; erewind(&x2fp); efread(&mut x2,&x2fp); efclose(&x2fp);
    let ntap: i32; let mut w: Vec<f32>;
    if index!=0 {
        ntap=ntr; w=vec![0.0;ntap as usize];
        if verbose!=0 { fprintf!(stderr(),"Using x values of key %s \n",key); }
        for i in 0..ntap as usize {
            let y=x2[i];
            let f=if y<min { if y>min-dkey { (1.0-((min-y)/dkey) as f64) as f32 } else { 0.0 } } else if y>max { if y<max+dkey { (max+dkey-y)/dkey } else { 0.0 } } else { 1.0 };
            w[i]=envelope(f,ttaper,low,true);
        }
        tr1=ntr.wrapping_add(1); tr2= -ntr;
    } else {
        ntap=if tr1.wrapping_sub(tr2)!=0 { tr1.wrapping_add(tr2).wrapping_add(2) } else { tr1.wrapping_add(1) }; if ntap<0 { err!("invalid taper length"); } w=vec![0.0;ntap as usize];
        if verbose!=0 { fprintf!(stderr(),"Using trace count \n"); }
        if tr2!=0 && tr1!=tr2 { for i in 0..=tr2 { *w.get_mut(1i32.wrapping_add(tr1).wrapping_add(i) as usize).unwrap_or_else(|| err!("invalid taper length"))=envelope(i as f32/tr2 as f32,ttaper,low,true); } }
        if tr1!=0 { for i in 0..=tr1 { *w.get_mut(i as usize).unwrap_or_else(|| err!("invalid taper length"))=envelope(i as f32/tr1 as f32,ttaper,low,true); } }
    }
    if verbose>1 { fprintf!(stderr(),"print taper weights: \n"); for i in 0..ntap as usize { fprintf!(stderr(),"Taper %i = %g (x=%g) \n",i as i32,w[i],x2.get(i).copied().unwrap_or(0.0)); } }
    erewind(&headerfp); erewind(&tracefp);
    let mut trbuf=vec![0.0f32;nt];
    let mut wp=0i32;
    for itr in 0..ntr {
        // At EOF the temporary-file reads leave the previous header and buffer intact.
        efread(tr.hdr.as_bytes_mut(),&headerfp);
        efread(&mut trbuf,&tracefp);
        let mut fac=1.0f32;
        if itr<tr1 { fac*=*w.get(wp as usize).unwrap_or_else(|| err!("invalid taper length")); wp+=1; }
        if itr==tr1 { wp=wp.wrapping_add(ntap.wrapping_sub(tr1).wrapping_sub(2)); }
        if itr>ntr.wrapping_sub(tr2).wrapping_sub(1) { fac*=*w.get(wp as usize).unwrap_or_else(|| err!("invalid taper length")); wp=wp.wrapping_sub(1); }
        if verbose!=0 && tr1.wrapping_add(tr2)!=0 && fac<1.0 { fprintf!(stderr(),"trace %i factor %g \n",itr,fac); }
        for i in 0..nt { tr.data[i]=fac*trbuf[i]; }
        tr.unscale=(1.0/fac as f64) as f32; puttr(&mut tr);
    }
    efclose(&tracefp); efclose(&headerfp);
    exit(CWP_Exit());
}
