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

//! PSCONTOUR - PostScript contours. Port of psplot/main/pscontour.c (Revision 1.18).
//! Author: Dave Hale; modified by Craig Artley, Zhenyue Liu, and Morten Wendell Pedersen.

use oxwave_core::prelude::*;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;
use oxwave_f_pscontour::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " PSCONTOUR - PostScript CONTOURing of a two-dimensional function f(x1,x2)",
    "\t\t\t\t\t\t\t\t\t",
    " pscontour n1= [optional parameters] <binaryfile >postscriptfile\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1                     number of samples in 1st (fast) dimension\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " d1=1.0                 sampling interval in 1st dimension\t\t",
    " f1=d1                  first sample in 1st dimension\t\t\t",
    " x1=f1,f1+d1,...        array of monotonic sampled values in 1st dimension",
    " n2=all                 number of samples in 2nd (slow) dimension\t",
    " d2=1.0                 sampling interval in 2nd dimension\t\t",
    " f2=d2                  first sample in 2nd dimension\t\t\t",
    " x2=f2,f2+d2,...        array of monotonic sampled values in 2nd dimension",
    " nc=5                   number of contour values\t\t\t",
    " dc=(zmax-zmin)/nc      contour interval\t\t\t\t",
    " fc=min+dc              first contour\t\t\t\t\t",
    " c=fc,fc+dc,...         array of contour values\t\t\t",
    " cwidth=1.0,...         array of contour line widths\t\t\t",
    " cgray=0.0,...          array of contour grays (0.0=black to 1.0=white)",
    " ccolor=none,...        array of contour colors; none means use cgray\t",
    " cdash=0.0,...          array of dash spacings (0.0 for solid)\t\t",
    " labelcf=1              first labeled contour (1,2,3,...)\t\t",
    " labelcper=1            label every labelcper-th contour\t\t",
    " nlabelc=nc             number of labeled contours (0 no contour label)",
    " nplaces=6              number of decimal places in contour label      ",
    " xbox=1.5               offset in inches of left side of axes box\t",
    " ybox=1.5               offset in inches of bottom side of axes box\t",
    " wbox=6.0               width in inches of axes box\t\t\t",
    " hbox=8.0               height in inches of axes box\t\t\t",
    " x1beg=x1min            value at which axis 1 begins\t\t\t",
    " x1end=x1max            value at which axis 1 ends\t\t\t",
    " d1num=0.0              numbered tic interval on axis 1 (0.0 for automatic)",
    " f1num=x1min            first numbered tic on axis 1 (used if d1num not 0.0)",
    " n1tic=1                number of tics per numbered tic on axis 1\t",
    " grid1=none             grid lines on axis 1 - none, dot, dash, or solid",
    " label1=                label on axis 1\t\t\t\t",
    " x2beg=x2min            value at which axis 2 begins\t\t\t",
    " x2end=x2max            value at which axis 2 ends\t\t\t",
    " d2num=0.0              numbered tic interval on axis 2 (0.0 for automatic)",
    " f2num=x2min            first numbered tic on axis 2 (used if d2num not 0.0)",
    " n2tic=1                number of tics per numbered tic on axis 2\t",
    " grid2=none             grid lines on axis 2 - none, dot, dash, or solid",
    " label2=                label on axis 2\t\t\t\t",
    " labelfont=Helvetica    font name for axes labels\t\t\t",
    " labelsize=18           font size for axes labels\t\t\t",
    " title=                 title of plot\t\t\t\t\t",
    " titlefont=Helvetica-Bold font name for title\t\t\t\t",
    " titlesize=24           font size for title\t\t\t\t",
    " labelcfont=Helvetica-Bold font name for contour labels\t\t",
    " labelcsize=6           font size of contour labels   \t\t\t",
    " labelccolor=black      color of contour labels   \t\t\t",
    " titlecolor=black       color of title\t\t\t\t\t",
    " axescolor=black        color of axes\t\t\t\t\t",
    " gridcolor=black        color of grid\t\t\t\t\t",
    " axeswidth=1            width (in points) of axes\t\t\t",
    " ticwidth=axeswidth     width (in points) of tic marks\t\t",
    " gridwidth=axeswidth    width (in points) of grid lines\t\t",
    " style=seismic          normal (axis 1 horizontal, axis 2 vertical) or\t",
    "                        seismic (axis 1 vertical, axis 2 horizontal)\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note.\t\t\t\t\t\t\t\t\t",
    " The line width of unlabeled contours is designed as a quarter of that\t",
    " of labeled contours. \t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " All color specifications may also be made in X Window style Hex format",
    " example:   axescolor=#255\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Legal font names are:\t\t\t\t\t\t\t",
    " AvantGarde-Book AvantGarde-BookOblique AvantGarde-Demi AvantGarde-DemiOblique Bookman-Demi Bookman-DemiItalic Bookman-Light Bookman-LightItalic ",
    " Courier Courier-Bold Courier-BoldOblique Courier-Oblique ",
    " Helvetica Helvetica-Bold Helvetica-BoldOblique Helvetica-Oblique ",
    " Helvetica-Narrow Helvetica-Narrow-Bold Helvetica-Narrow-BoldOblique ",
    " Helvetica-Narrow-Oblique NewCentrySchlbk-Bold NewCenturySchlbk-BoldItalic NewCenturySchlbk-Roman Palatino-Bold  ",
    " Palatino-BoldItalic Palatino-Italics Palatino-Roman ",
    " SanSerif-Bold SanSerif-BoldItalic SanSerif-Roman ",
    " Symbol Times-Bold Times-BoldItalic ",
    " Times-Roman Times-Italic ZapfChancery-MediumItalic ",
    " type:   sudoc pscontour    for more information\t\t\t",
];

fn fp(n: &str, d: f32) -> f32 {
    let mut v = d;
    getparfloat(n, &mut v);
    v
}
fn ip(n: &str, d: i32) -> i32 {
    let mut v = d;
    getparint(n, &mut v);
    v
}
fn sp(n: &str, d: &str) -> String {
    let mut v = d.to_string();
    getparstring(n, &mut v);
    v
}
fn grid(s: &str) -> i32 {
    match s {
        "dot" => 1,
        "dash" => 2,
        "solid" => 3,
        _ => 0,
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut n1 = countparval("x1");
    if n1 == 0 && getparint("n1", &mut n1) == 0 {
        err!("Must specify n1!\n");
    }
    let mut x1 = vec![0.0f32; n1 as usize];
    if getparfloat("x1", &mut x1) == 0 {
        let d = fp("d1", 1.0);
        let f = fp("f1", d);
        for i in 0..n1 {
            x1[i as usize] = f + i as f32 * d;
        }
    }
    let (mut x1min, mut x1max) = (x1[0], x1[0]);
    for &v in &x1[1..] {
        x1min = if x1min < v { x1min } else { v };
        x1max = if x1max > v { x1max } else { v };
    }
    let infp = stdin();
    let mut n2 = countparval("x2");
    if n2 == 0 && getparint("n2", &mut n2) == 0 {
        if efseeko(&infp, 0, SEEK_END) != 0 {
            err!("must specify n2 if in a pipe!\n");
        }
        let nf = (eftello(&infp) / 4) as i32;
        efseeko(&infp, 0, SEEK_SET);
        n2 = nf.checked_div(n1).unwrap_or(0);
    }
    let mut x2 = vec![0.0f32; n2 as usize];
    if getparfloat("x2", &mut x2) == 0 {
        let d = fp("d2", 1.0);
        let f = fp("f2", d);
        for i in 0..n2 {
            x2[i as usize] = f + i as f32 * d;
        }
    }
    let (mut x2min, mut x2max) = (x2[0], x2[0]);
    for &v in &x2[1..] {
        x2min = if x2min < v { x2min } else { v };
        x2max = if x2max > v { x2max } else { v };
    }
    let mut z = vec![0.0f32; (n1 * n2) as usize];
    if fread(&mut z, &infp) != z.len() {
        err!("error reading input file!\n");
    }
    let mut w = vec![0.0f32; z.len()];
    let (mut zmin, mut zmax) = (z[0], z[0]);
    for &v in &z {
        zmin = if zmin < v { zmin } else { v };
        zmax = if zmax > v { zmax } else { v };
    }
    let mut c = vec![0.0f32; 200];
    let mut nc = getparfloat("c", &mut c);
    if nc == 0 {
        nc = ip("nc", 5);
        let dc = fp("dc", (zmax - zmin) / nc as f32);
        let fc = fp("fc", zmin + dc);
        for i in 0..nc {
            c[i as usize] = fc + i as f32 * dc;
        }
    }
    let mut cw = vec![1.0f32; 200];
    let mut cg = vec![0.0f32; 200];
    let mut cd = vec![0.0f32; 200];
    for (name, a) in [("cwidth", &mut cw), ("cgray", &mut cg), ("cdash", &mut cd)] {
        let n = getparfloat(name, &mut *a);
        if n != 0 {
            for i in n..nc {
                a[i as usize] = a[(n - 1) as usize];
            }
        }
    }
    let mut colors = vec!["none".to_string(); 200];
    let mut scolor = String::new();
    if getparstring("ccolor", &mut scolor) != 0 {
        let p: Vec<&str> = scolor.split(',').filter(|s| !s.is_empty()).collect();
        for i in 0..nc as usize {
            colors[i] = p[i.min(p.len() - 1)].to_string();
        }
    }
    let labelcf = ip("labelcf", 1);
    let labelcper = ip("labelcper", 1);
    let nlabelc = ip("nlabelc", nc);
    let labelcsize = fp("labelcsize", 6.0);
    let labelcfont = sp("labelcfont", "Helvetica-Bold");
    let labelccolor = sp("labelccolor", "black");
    let nplaces = ip("nplaces", 6);
    let mut xbox = fp("xbox", 1.5);
    let mut ybox = fp("ybox", 1.5);
    let mut wbox = fp("wbox", 6.0);
    let mut hbox = fp("hbox", 8.0);
    let x1beg = fp("x1beg", x1min);
    let x1end = fp("x1end", x1max);
    let d1num = fp("d1num", 0.0);
    let f1num = fp("f1num", x1min);
    let n1tic = ip("n1tic", 1);
    let grid1 = grid(&sp("grid1", "none"));
    let label1 = sp("label1", "");
    let x2beg = fp("x2beg", x2min);
    let x2end = fp("x2end", x2max);
    let d2num = fp("d2num", 0.0);
    let f2num = fp("f2num", 0.0);
    let n2tic = ip("n2tic", 1);
    let grid2 = grid(&sp("grid2", "none"));
    let label2 = sp("label2", "");
    let labelfont = sp("labelfont", "Helvetica");
    let labelsize = fp("labelsize", 18.0);
    let title = sp("title", "");
    let titlefont = sp("titlefont", "Helvetica-Bold");
    let titlesize = fp("titlesize", 24.0);
    let titlecolor = sp("titlecolor", "black");
    let axescolor = sp("axescolor", "black");
    let gridcolor = sp("gridcolor", "black");
    let axeswidth = fp("axeswidth", 1.0);
    let ticwidth = fp("ticwidth", axeswidth);
    let gridwidth = fp("gridwidth", axeswidth);
    let styles = sp("style", "seismic");
    checkpars();
    let style = if styles == "normal" { 0 } else { 1 };
    xbox *= 72.0;
    ybox *= 72.0;
    wbox *= 72.0;
    hbox *= 72.0;
    let mut bbox = [0; 4];
    psAxesBBox(
        xbox, ybox, wbox, hbox, &labelfont, labelsize, &titlefont, titlesize, style, &mut bbox,
    );
    boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
    begineps();
    gsave();
    rectclip(xbox, ybox, wbox, hbox);
    let x1size = if style == 0 { wbox } else { hbox };
    let x2size = if style == 0 { hbox } else { wbox };
    translate(xbox, ybox);
    if style != 0 {
        rotate(-90.0);
        translate(-hbox, 0.0);
    }
    let x1scale = x1size / (x1end - x1beg);
    let x2scale = x2size / (x2end - x2beg);
    translate(-x1beg * x1scale, -x2beg * x2scale);
    for v in &mut x1 {
        *v *= x1scale;
    }
    for v in &mut x2 {
        *v *= x2scale;
    }
    for ic in 0..nc {
        let i = ic as usize;
        setlinewidth(cw[i]);
        if colors[i] != "none" {
            setcolor(&colors[i]);
        } else {
            setgray(cg[i]);
        }
        if cd[i] != 0.0 {
            setdash(&cd[i..i + 1], 0.0);
        } else {
            setdash(&[], 0.0);
        }
        let mut lcsize = 0.0;
        if nlabelc > 0 {
            if (ic - labelcf + 1).checked_rem(labelcper).unwrap_or(0) == 0
                && ic >= labelcf - 1
                && ic < labelcf - 1 + labelcper * nlabelc
            {
                setlinewidth(cw[i]);
                lcsize = labelcsize;
            } else {
                setlinewidth((0.25 * cw[i] as f64) as f32);
            }
        }
        psContour(
            c[i],
            &x1,
            &x2,
            &mut z,
            lcsize,
            &labelcfont,
            &labelccolor,
            &mut w,
            nplaces,
        );
    }
    grestore();
    psAxesBox(
        xbox,
        ybox,
        wbox,
        hbox,
        x1beg,
        x1end,
        0.0,
        0.0,
        d1num,
        f1num,
        n1tic,
        grid1,
        &label1,
        x2beg,
        x2end,
        0.0,
        0.0,
        d2num,
        f2num,
        n2tic,
        grid2,
        &label2,
        &labelfont,
        labelsize,
        &title,
        &titlefont,
        titlesize,
        &titlecolor,
        &axescolor,
        &gridcolor,
        ticwidth,
        axeswidth,
        gridwidth,
        style,
    );
    showpage();
    endeps();
    exit(CWP_Exit());
}
