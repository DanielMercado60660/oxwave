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
//! SUPICKAMP - pick amplitudes in a resampled window.
//! Port of su/main/picking/supickamp.c (Revision 1.12).
//! Credits: CWP: Andreas Rueger; MTU: David Forel (verbose=2).
use oxwave_core::prelude::*;
use oxwave_linear_interp::intlin;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t   \t\t\t\t\t\t\t\t",
    " SUPICKAMP - pick amplitudes within user defined and resampled window\t",
    " \t   \t\t\t\t\t\t\t\t",
    "   supickamp <stdin >stdout d2=  [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " d2=\t\t   sampling interval for slow dimension\t\t\t",
    "\t\t\t(required if key-parameter not specified)\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " key= \t\t\tKey header word specifying trace offset \t",
    " \t\t\t(alternatively, specify d2,x2beg)\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " x_above=\t\tarray of lateral position values   \t\t",
    " \t\t\t(upper window corner)\t\t\t\t",
    " t_above=\t\tarray of time values   \t\t\t\t",
    " \t\t\t(upper window corner)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  ... or input via files:\t\t\t\t\t\t",
    " t_xabove=\t\tfile containing time and lateral position values",
    " \t\t\t(upper window corner)\t\t\t\t",
    " t_xbelow=\t\tfile containing time and lateral position values",
    " \t\t\t(lower window corner)\t\t\t\t",
    " wl=\t\t   \twindow width if t_xbelow is not specified\t",
    "\t\t\t(No windowing if not specified)\t\t \t",
    " \t\t\t\t\t\t\t\t\t",
    " dt_resamp=dt\t  resampling interval within pick window\t  \t",
    "\t\t\t(dt has to come from trace headers)\t\t",
    " tmin=0\t\tminimum time in input trace\t\t\t",
    " x2beg=0\t\tfirst lateral position\t\t\t\t",
    " format=ascii \t\twrite ascii data to stdout\t\t\t",
    "\t\t\t =binary for binary floats to stdout\t\t",
    " verbose=1 \t\twrites complete  pick information into outpar   ",
    "\t\t\t=2 writes complete pick information into outpar\t",
    "\t\t\t   in tab-delimited column format\t\t",
    " outpar=/dev/tty\toutput parameter file; contains output\t\t",
    "\t\t\t\t\tfrom verbose\t\t\t",
    " arg1=max\t\toutput (first dimension) to stdout\t\t",
    " arg2=i2\t\toutput (second dimension) to stdout\t\t",
    "\t\t\t(see notes for other options)\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Window can be defined using\t\t\t\t\t\t",
    " (1)   vectors x_above, t_above, [wl]\t\t\t\t\t",
    " (2)   file  t_xabove, [wl]\tor\t\t\t\t\t",
    " (3)   files t_xabove, t_xbelow\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " files t_xabove, t_xbelow can be generated using xwigb's picking\t",
    " algorithm. The lateral positions have to be monotonically increasing  ",
    " or decreasing for both vector and file input.\t\t\t\t",
    " verbose=1 or 2 writes min, max, abs[max], energy and associated times\t",
    " tmin,tmax,tabs to outpar, together with global values. verbose=0\t",
    " only outputs global values.\t\t\t\t\t\t",
    " Acceptable arg-parameters for lateral positions are\t\t\t",
    " (1) x2   (2) i2 = trace number\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " If key=keyword is set, then the values of x2 are taken from the header",
    " field represented by the keyword (for example key=offset)\t\t",
    " Type\tsukeyword -o   to see the complete list of SU keywords.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

#[derive(Clone, Copy)]
struct Pick { max:f32, min:f32, abs:f32, tmax:f32, tmin:f32, tabs:f32, energy:f32, imax:i32, imin:i32, iabs:i32, ienergy:i32 }
fn initpick() -> Pick { Pick { min:f32::MAX, max:f32::MIN_POSITIVE, abs:-0.1, tmin:0., tmax:0., tabs:0., energy:-0.1, imax:0, imin:0, iabs:0, ienergy:0 } }
fn tracepick(p: &mut Pick, t: &[f32], val: &mut [f32], nt: usize) {
    let (mut min, mut max, mut abs) = (f32::MAX, -f32::MAX, 0.0f32);
    let (mut tmin, mut tmax, mut tabs) = (t[0], t[0], t[0]); let mut energy = 0.0f32;
    for it in 0..nt {
        if val[it] > max { max=val[it]; tmax=t[it]; } if val[it] < min { min=val[it]; tmin=t[it]; }
        let a = if val[it] < 0.0 { -val[it] } else { val[it] }; if a > abs { abs=a; tabs=t[it]; }
        energy += val[it]*val[it]; val[it]=0.0;
    }
    p.min=min; p.max=max; p.abs=abs; p.tmin=tmin; p.tmax=tmax; p.tabs=tabs; p.energy=energy;
}
fn globalpick(p: &Pick, g: &mut Pick) {
    if p.max > g.max { g.max=p.max; g.tmax=p.tmax; g.imax=p.iabs; } if p.min < g.min { g.min=p.min; g.tmin=p.tmin; g.imin=p.iabs; }
    if p.abs > g.abs { g.abs=p.abs; g.tabs=p.tabs; g.iabs=p.iabs; } if p.energy > g.energy { g.energy=p.energy; g.ienergy=p.iabs; }
}
fn getarg(mode: &str, p: &Pick, x:f32, i:i32) -> f32 { match mode { "max"=>p.max, "tmax"=>p.tmax, "min"=>p.min, "tmin"=>p.tmin, "abs"=>p.abs, "tabs"=>p.tabs, "energy"=>p.energy, "x2"=>x, "i2"=>i as f32, _=>0. } }
fn valid_mode(s:&str) -> bool { matches!(s,"max"|"tmax"|"min"|"tmin"|"abs"|"tabs"|"energy"|"x2"|"i2") }
fn read_picks(file:&str) -> (Vec<f32>,Vec<f32>) {
    let fp=efopen(file,"r"); let (mut x,mut t)=(Vec::new(),Vec::new());
    while let Some(a)=fscanf_f32(&fp) { if let Some(b)=fscanf_f32(&fp) { t.push(a); x.push(b); } else { break; } }
    efclose(&fp); (x,t)
}
fn interp(x:&[f32], t:&[f32], f:f32) -> f32 { let mut result=[0.]; intlin(x,t,t[0],t[t.len()-1],&[f],&mut result); result[0] }
fn main() {
    let argv=args(); let mut tr=Trace::new(); initargs(&argv); requestdoc(1,SDOC);
    let mut d2=0.0f32; getparfloat("d2",&mut d2); let mut key=String::new(); let has_key=getparstring("key",&mut key)!=0;
    if has_key { d2=0.; } else if d2==0. { err!(" Specify d2"); }
    let mut verbose=1; getparint("verbose",&mut verbose); let mut arg1="max".to_string(); getparstring("arg1",&mut arg1);
    let mut arg2="x2".to_string(); getparstring("arg2",&mut arg2); let mut wl=0.0f32; getparfloat("wl",&mut wl);
    let mut dt_resamp=0.0f32; getparfloat("dt_resamp",&mut dt_resamp); let mut x2beg=0.0f32; getparfloat("x2beg",&mut x2beg);
    let mut tmin_in=0.0f32; getparfloat("tmin",&mut tmin_in); let mut format="ascii".to_string(); getparstring("format",&mut format);
    let mut outpar="/dev/tty".to_string(); getparstring("outpar",&mut outpar); let out=efopen(&outpar,"w");
    if verbose==2 { fprintf!(out,"trace\tx\tt_above\tt_below\ttrace #\tx2\tmax\t tmax\tmin\t tmin\tabs\t tabs\tenergy in pick window\n"); }
    let (typ,index)=if d2==0. { (hdtype(&key),getindex(&key)) } else { ("",0) };
    let mut afile=String::new(); let has_a=getparstring("t_xabove",&mut afile)!=0;
    let (xa,ta)=if has_a { read_picks(&afile) } else {
        let nx=countparval("x_above"); let nt=countparval("t_above"); if nx!=nt { err!("lengths of x_above, t_above must be the same"); }
        let mut x=vec![0.0f32;nx as usize]; let mut t=vec![0.0f32;nx as usize]; getparfloat("x_above",&mut x); getparfloat("t_above",&mut t); (x,t)
    };
    let mut bfile=String::new(); let has_b=getparstring("t_xbelow",&mut bfile)!=0;
    let (xb,tb)=if has_b { if !has_a { err!(" Specify t_xabove -file"); } read_picks(&bfile) } else { (Vec::new(),Vec::new()) };
    if has_b && xb.is_empty() { err!("t_xbelow file contains no picks"); }
    if has_b && xa.is_empty() { err!("t_xabove file contains no picks"); }
    if has_b && xb.len()>xa.len() { err!("t_xbelow has more picks than t_xabove"); }
    let nx=if has_b { xb.len() } else { xa.len() }; checkpars();
    if has_b && has_a { warn!("window determined by above_ and below_ files "); } else if has_a && wl!=0. { warn!("window determined by above_file and wl"); }
    else if has_a { warn!("window determined by above_file and trace end"); } else if wl!=0. && nx!=0 { warn!("window determined by above_vector and wl"); }
    else if wl==0. && nx!=0 { warn!("window determined by above_vector and trace length"); } else if wl==0. { warn!("window determined by full trace length"); } else { err!("window not properly specified"); }
    if !valid_mode(&arg1) { err!(" Unknown arg1 "); } if !valid_mode(&arg2) { err!(" Unknown arg2"); } if format!="ascii" && format!="binary" { err!(" Unknown format"); }
    if gettr(&mut tr)==0 { err!("can't read first trace"); } if tr.dt==0 { err!("dt header field must be set"); }
    if tr.delrt!=0 { tmin_in=(tr.delrt as f64/1000.0) as f32; }
    let nt=tr.ns as usize; let dt=(tr.dt as f64/1000000.0) as f32; if dt_resamp==0. { dt_resamp=dt; }
    let nts=(((nt as i32-1) as f32*dt/dt_resamp) as f64+1.0) as usize;
    // C interpolates one additional, initially zero time sample.
    let mut ts=vec![0.0f32;nts.saturating_add(2).max(2)]; let mut vs=ts.clone(); let mut p=initpick(); let mut g=initpick(); let mut trnum=0i32;
    loop {
        p.iabs=trnum; trnum+=1; let fval=if d2==0. { let mut val=Value::default(); gethval(&tr,index,&mut val); vtof(typ,val) } else { x2beg+(trnum-1) as f32*d2 };
        // C uses the lower count for both arrays; longer lower files are
        // rejected above because C would read past the upper allocation.
        let mut tmin=if nx==0 { tmin_in } else { interp(&xa[..nx],&ta[..nx],fval) }; let end=tmin_in+(nt as i32-1) as f32*dt;
        let mut tmax=if has_b { interp(&xb,&tb,fval) } else if wl!=0. { tmin+wl } else { end }; if tmin<tmin_in { tmin=tmin_in; } if tmax>end { tmax=end; }
        if verbose==1 { fprintf!(out,"trace %i x=%f t_above=%f t_below=%f \n",trnum-1,fval,tmin,tmax); }
        let mut itime=0usize; let mut tvalue=tmin; while tvalue<=tmax { if itime>=ts.len()-1 { ts.resize(itime+2,0.); vs.resize(itime+2,0.); } ts[itime]=tvalue; itime+=1; tvalue+=dt_resamp; }
        ints8r(dt,tmin_in,&tr.data[..nt],0.,0.,&ts[..itime+1],&mut vs[..itime+1]); tracepick(&mut p,&ts,&mut vs,itime); p.energy*=dt_resamp; globalpick(&p,&mut g);
        if verbose==1 { fprintf!(out,"trace #=%i \t   x2=%e \n",p.iabs,fval); fprintf!(out,"max=%e \t tmax=%e \n",p.max,p.tmax); fprintf!(out,"min=%e \t tmin=%e \n",p.min,p.tmin); fprintf!(out,"abs=%e \t tabs=%e \n",p.abs,p.tabs); fprintf!(out,"energy in pick window =%e \n \n",p.energy); }
        if verbose==2 { fprintf!(out,"%i\t %f\t %f\t %f\t%i\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",trnum-1,fval,tmin,tmax,p.iabs,fval,p.max,p.tmax,p.min,p.tmin,p.abs,p.tabs,p.energy); }
        let a=getarg(&arg1,&p,fval,trnum-1); let b=getarg(&arg2,&p,fval,trnum-1); if format=="ascii" { fprintf!(stdout(),"%e \t %e\n",a,b); } else { efwrite(&[b],&stdout()); efwrite(&[a],&stdout()); }
        if gettr(&mut tr)==0 { break; }
    }
    fprintf!(out,"\nglobal max=%e \t #=%i \t tmax=%e \n",g.max,g.imax,g.tmax); fprintf!(out,"global min=%e \t #=%i \t tmin=%e \n",g.min,g.imin,g.tmin);
    fprintf!(out,"global abs=%e \t #=%i \t tabs=%e \n",g.abs,g.iabs,g.tabs); fprintf!(out,"max energy in pick window =%e \t #=%i \n \n",g.energy,g.ienergy);
    if format=="binary" { fprintf!(out,"\n npairs = %d\n",trnum); } exit(CWP_Exit());
}
