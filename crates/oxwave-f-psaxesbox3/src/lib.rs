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

//! Port of psplot/lib/psAxesBox3.c.
//! Credits: Dave Hale and Zhiming Li, Colorado School of Mines.
#![allow(non_snake_case, clippy::too_many_arguments)]
use oxwave_f_basic::{fontheight, fontwidth};
use oxwave_f_basic::*;
use oxwave_axis_scale::scaxis;
use oxwave_core::sprintf;

/// C: `void psAxesBBox3(float x, float y, float width, float height, char *labelFont, float labelSize, char *titleFont, float titleSize, int style, int bbox[]);`
pub fn psAxesBBox3(x: f32, y: f32, width: f32, height: f32, labelFont: &str, labelSize: f32, titleFont: &str, titleSize: f32, style: i32, bbox: &mut [i32]) {
    let labelCH = fontheight(labelFont, labelSize);
    let labelCW = fontwidth(labelFont, labelSize);
    let titleCH = fontheight(titleFont, titleSize);
    bbox[0] = (x as f64 - 6.0 * labelCH as f64) as i32;
    bbox[1] = if style == 0 { (y as f64 - 3.0 * labelCH as f64 - 2.0 * titleCH as f64) as i32 } else { (y as f64 - 2.0 * titleCH as f64) as i32 };
    bbox[2] = (x + width + labelCW) as i32;
    bbox[3] = if style == 0 { ((y + height) as f64 + 2.0 * titleCH as f64) as i32 } else { ((y + height) as f64 + 3.0 * labelCH as f64 + 2.0 * titleCH as f64) as i32 };
}

/// C: `void psAxesBox3(float x, float y, float width, float height, float x1Beg, float x1End, float p1Beg, float p1End, float d1Num, float f1Num, int n1Tic, int grid1, char *label1, float x2Beg, float x2End, float p2Beg, float p2End, float d2Num, float f2Num, int n2Tic, int grid2, char *label2, char *labelFont, float labelSize, char *title, char *titleFont, float titleSize, int style, char *title2);`
pub fn psAxesBox3(x:f32,y:f32,width:f32,height:f32,x1Beg:f32,x1End:f32,p1Beg:f32,p1End:f32,mut d1Num:f32,mut f1Num:f32,n1Tic:i32,grid1:i32,label1:&str,x2Beg:f32,x2End:f32,p2Beg:f32,p2End:f32,mut d2Num:f32,mut f2Num:f32,n2Tic:i32,grid2:i32,label2:&str,labelFont:&str,labelSize:f32,title:&str,titleFont:&str,titleSize:f32,style:i32,title2:&str) {
    let lh=fontheight(labelFont,labelSize); let lw=fontwidth(labelFont,labelSize); let la=fontascender(labelFont,labelSize); let ld=fontdescender(labelFont,labelSize);
    let th=fontheight(titleFont,titleSize); let ta=fontascender(titleFont,titleSize); let td=fontdescender(titleFont,titleSize);
    let (s1,s2)=if style==0 {(width,height)} else {(height,width)};
    if d1Num==0.0 {let mut n=(s1/(4.0*lw)) as i32; scaxis(x1Beg,x1End,&mut n,&mut d1Num,&mut f1Num);}
    if d2Num==0.0 {let mut n=(s2/(4.0*lw)) as i32; scaxis(x2Beg,x2End,&mut n,&mut d2Num,&mut f2Num);}
    gsave(); translate(x,y); if style!=0 {rotate(-90.0); translate(-height,0.0);} newpath(); setfont(labelFont,labelSize);
    let tic= -(0.3*lh as f64) as f32;
    let amin=if x1Beg<x1End{x1Beg}else{x1End}; let amax=if x1Beg>x1End{x1Beg}else{x1End}; let zero=(0.0001*(amax-amin) as f64) as f32;
    let scale=s1/(x1End+p1End-x1Beg-p1Beg); let base= -scale*(x1Beg+p1Beg); let numb=tic-la;
    let mut count=0; let mut a=f1Num;
    while a<=amax {if a>=amin {let v=base+scale*a; moveto(v,0.0); lineto(v,tic); let text=sprintf!("%.4g",if a> -zero && a<zero {0.0} else {a as f64}); moveto((v as f64+0.3*lh as f64) as f32,numb); rotate(90.0); count=count.max(text.len()); justshow(-0.5,&text); rotate(-90.0);} a+=d1Num;} stroke();
    ticks(f1Num,d1Num,n1Tic,amin,amax,base,scale,tic,true); grid(grid1,f1Num,d1Num,amin,amax,base,scale,s2,true);
    count=count.min(4); let labelb=(numb as f64-(count as f32*lw) as f64/2.0) as f32; moveto((s1 as f64/2.0) as f32,labelb); rotate(180.0); justshow(-0.5,label1); rotate(-180.0);
    let amin=if x2Beg<x2End{x2Beg}else{x2End}; let amax=if x2Beg>x2End{x2Beg}else{x2End}; let zero=(0.0001*(amax-amin) as f64) as f32;
    let scale=s2/(x2End+p2End-x2Beg-p2Beg); let base= -scale*(x2Beg+p2Beg); let numb=tic+ld; let labelb=numb-lh;
    let mut a=f2Num; while a<=amax {if a>=amin {let v=base+scale*a; moveto(0.0,v); lineto(tic,v); let text=sprintf!("%.4g",if a> -zero && a<zero {0.0} else {a as f64}); moveto(numb,v); rotate(90.0); justshow(-0.5,&text); rotate(-90.0);} a+=d2Num;} stroke();
    ticks(f2Num,d2Num,n2Tic,amin,amax,base,scale,tic,false); grid(grid2,f2Num,d2Num,amin,amax,base,scale,s1,false);
    moveto(labelb,(s2 as f64/2.0) as f32); rotate(90.0); justshow(-0.5,label2); rotate(-90.0); setfont(titleFont,titleSize);
    if style==0 {moveto((s1 as f64/2.0) as f32,(s2 as f64+0.5*th as f64-td as f64) as f32); justshow(-0.5,title); moveto((s1 as f64/2.0) as f32,(-0.5*th as f64+td as f64+labelb as f64) as f32); justshow(-0.5,title2);} else {moveto((s1 as f64+0.5*th as f64+ta as f64) as f32,(s2 as f64/2.0) as f32); rotate(90.0); justshow(-0.5,title); rotate(-90.0); moveto((labelb as f64-0.5*th as f64-ta as f64) as f32,(s2 as f64/2.0) as f32); rotate(90.0); justshow(-0.5,title2); rotate(-90.0);}
    moveto(0.0,0.0); lineto(s1,0.0); lineto(s1,s2); lineto(0.0,s2); lineto(0.0,0.0); stroke(); grestore();
}
fn ticks(first:f32,step:f32,n:i32,amin:f32,amax:f32,base:f32,scale:f32,tic:f32,axis:bool) {let dt=step/n as f32; let mut a=first-n as f32*dt-dt; while a<=amax {if a>=amin {let v=base+scale*a; if axis {moveto(v,0.0);lineto(v,tic/2.0);} else {moveto(0.0,v);lineto(tic/2.0,v);}} a+=dt;} stroke();}
fn grid(code:i32,first:f32,step:f32,amin:f32,amax:f32,base:f32,scale:f32,size:f32,axis:bool) {let dash:&[f32]=match code {3=>&[],2=>&[10.0],1=>&[1.0,5.0],_=>return}; let mut a=first; while a<=amax {if a>=amin {let v=base+scale*a; if axis {moveto(v,0.0);lineto(v,size);} else {moveto(0.0,v);lineto(size,v);}} a+=step;} setdash(dash,0.0);stroke();setdash(&[],0.0);}
