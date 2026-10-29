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
//! PSWIGP, port of psplot/main/pswigp.c, revision 1.29.
//! Credits: Dave Hale and Craig Artley, Colorado School of Mines; Dave Hale,
//! Advance Geophysical; Morten Wendell Pedersen, Aarhus University; Brian K.
//! Macy, Phillips Petroleum; Paul Michaels, Boise State University; Ekkehart
//! Tessmer, University of Hamburg; Chris Lienr, University of Arkansas.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSWIGP - PostScript WIGgle-trace plot of f(x1,x2) via Polygons\t",
    " Best for few traces.  Use PSWIGB (Bitmap version) for many traces.\t",
    " \t\t\t\t\t\t\t\t\t",
    " pswigp n1= [optional parameters] <binaryfile >postscriptfile\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1                     number of samples in 1st (fast) dimension\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " d1=1.0                 sampling interval in 1st dimension\t\t",
    " f1=0.0                 first sample in 1st dimension\t\t\t",
    " n2=all                 number of samples in 2nd (slow) dimension\t",
    " d2=1.0                 sampling interval in 2nd dimension\t\t",
    " f2=0.0                 first sample in 2nd dimension\t\t\t",
    " x2=f2,f2+d2,...        array of sampled values in 2nd dimension\t",
    " bias=0.0               data value corresponding to location along axis 2",
    " perc=100.0             percentile for determining clip\t\t",
    " clip=(perc percentile) data values < bias+clip and > bias-clip are clipped",
    " xcur=1.0               wiggle excursion in traces corresponding to clip",
    " fill=1\t\t\t=0 for no fill;\t\t\t\t",
    "\t\t\t\t>0 for pos. fill;\t\t\t",
    "\t\t\t\t<0 for neg. fill\t\t\t",
    "                               =2 for pos. fill solid, neg. fill grey  ",
    "                               =-2for neg. fill solid, pos. fill grey  ",
    "                       SHADING: 2<=abs(fill)<=5  2=lightgrey 5=black   ",
    " linewidth=1.0         linewidth in points (0.0 for thinest visible line)",
    " tracecolor=black       color of traces; should contrast with background",
    " backcolor=none         color of background; none means no background\t",
    " verbose=1              =1 for info printed on stderr (0 for no info)\t",
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
    " titlecolor=black       color of title\t\t\t\t\t",
    " axescolor=black        color of axes\t\t\t\t\t",
    " gridcolor=black        color of grid\t\t\t\t\t",
    " axeswidth=1            width (in points) of axes\t\t\t",
    " ticwidth=axeswidth     width (in points) of tic marks\t\t",
    " gridwidth=axeswidth    width (in points) of grid lines\t\t",
    " frame=1\t\t =0 wiggle traces only, no frame\t\t",
    " style=seismic          normal (axis 1 horizontal, axis 2 vertical) or\t",
    "                        seismic (axis 1 vertical, axis 2 horizontal)\t",
    "\t\t\t\t\t\t\t\t\t",
    " curve=curve1,curve2,...  file(s) containing points to draw curve(s)   ",
    " npair=n1,n2,n2,...            number(s) of pairs in each file         ",
    " curvecolor=black,..    color of curve(s)                              ",
    " curvewidth=axeswidth   width (in points) of curve(s)                  ",
    " curvedash=0            solid curve(s), dash indices 1,...,11 produce  ",
    "                        curve(s) with various dash styles              ",
    "\t\t\t\t\t\t\t\t\t",
    " Note:  linewidth=0.0 produces the thinest possible line on the output.\t",
    " device.  Thus the result is device-dependent, put generally looks the\t",
    " best for seismic traces.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The curve file is an ascii file with the points specified as x1 x2 pairs,",
    " one pair to a line.  A \"vector\" of curve files and curve colors may ",
    " be specified as curvefile=file1,file2,etc. and similarly\t\t",
    " curvecolor=color1,color2,etc, and the number of pairs of values\t",
    " in each file as npair=npair1,npair2,... .                             ",
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
    "\t\t\t\t\t\t\t\t\t",
];
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;
use oxwave_f_psdrawcurve::*;
use oxwave_f_pswiggle::*;
use oxwave_sort::qkfind;

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
    let mut n1 = 0;
    if getparint("n1", &mut n1) == 0 {
        err!("Must specify number of samples in 1st dimension!\n");
    }
    let d1 = fp("d1", 1.0);
    let f1 = fp("f1", 0.0);
    let x1min = if d1 > 0.0 {
        f1
    } else {
        f1 + (n1 - 1) as f32 * d1
    };
    let x1max = if d1 < 0.0 {
        f1
    } else {
        f1 + (n1 - 1) as f32 * d1
    };
    let infp = stdin();
    let mut n2 = countparval("x2");
    if n2 == 0 && getparint("n2", &mut n2) == 0 {
        if efseeko(&infp, 0, SEEK_END) != 0 {
            err!("must specify n2 if in a pipe!\n");
        }
        let nfloats = (eftello(&infp) / 4) as i32;
        efseek(&infp, 0, SEEK_SET);
        n2 = nfloats.checked_div(n1).unwrap_or(0);
    }
    let mut x2 = vec![0.0f32; n2 as usize];
    if getparfloat("x2", &mut x2) == 0 {
        let d2 = fp("d2", 1.0);
        let f2 = fp("f2", 0.0);
        for i in 0..n2 as usize {
            x2[i] = f2 + i as f32 * d2;
        }
    }
    let mut x2min = x2[0];
    let mut x2max = x2[0];
    for &v in &x2[1..] {
        if v < x2min {
            x2min = v;
        }
        if v > x2max {
            x2max = v;
        }
    }
    let curve = countparval("curve") as usize;
    let mut curvefile = Vec::new();
    let mut npair = vec![0i32; curve];
    let mut curvecolor = Vec::new();
    let mut x1curve = Vec::new();
    let mut x2curve = Vec::new();
    if curve != 0 {
        getparstringarray("curve", &mut curvefile);
        getparint("npair", &mut npair);
        let nc = countparval("curvecolor");
        getparstringarray("curvecolor", &mut curvecolor);
        if nc == 0 {
            curvecolor.push("black".to_string());
        }
        while curvecolor.len() < curve {
            curvecolor.push(curvecolor.last().unwrap().clone());
        }
        for j in 0..curve {
            let file = efopen(&curvefile[j], "r");
            let mut a = vec![0.0; npair[j] as usize];
            let mut b = a.clone();
            for i in 0..a.len() {
                if let Some(v) = fscanf_f32(&file) {
                    a[i] = v;
                }
                if let Some(v) = fscanf_f32(&file) {
                    b[i] = v;
                }
            }
            efclose(&file);
            x1curve.push(a);
            x2curve.push(b);
        }
    }
    let nz = n1 * n2;
    let mut z = vec![0.0f32; nz as usize];
    if fread(&mut z, &infp) != nz as usize {
        err!("Error reading input file!\n");
    }
    let mut bias = 0.0;
    if getparfloat("bias", &mut bias) != 0 && bias != 0.0 {
        for v in &mut z {
            *v -= bias;
        }
    }
    let mut clip = 0.0;
    if getparfloat("clip", &mut clip) == 0 {
        let perc = fp("perc", 100.0);
        let mut temp: Vec<f32> = z.iter().map(|v| v.abs()).collect();
        let mut iz = ((nz as f32 * perc) as f64 / 100.0) as i32;
        if iz < 0 {
            iz = 0;
        }
        if iz > nz - 1 {
            iz = nz - 1;
        }
        qkfind(iz as usize, &mut temp);
        clip = temp[iz as usize];
    }
    let verbose = ip("verbose", 1);
    if verbose != 0 {
        warn!("clip=%g", clip);
    }
    let fill = ip("fill", 1);
    let linewidth = fp("linewidth", 1.0);
    let tracecolor = sp("tracecolor", "black");
    let backcolor = sp("backcolor", "none");
    let mut xbox = fp("xbox", 1.5);
    let mut ybox = fp("ybox", 1.5);
    let mut wbox = fp("wbox", 6.0);
    let mut hbox = fp("hbox", 8.0);
    let mut x1beg = fp("x1beg", x1min);
    let mut x1end = fp("x1end", x1max);
    let d1num = fp("d1num", 0.0);
    let f1num = fp("f1num", x1min);
    let n1tic = ip("n1tic", 1);
    let grid1 = grid(&sp("grid1", "none"));
    let x2beg = fp("x2beg", x2min);
    let x2end = fp("x2end", x2max);
    let d2num = fp("d2num", 0.0);
    let f2num = fp("f2num", 0.0);
    let n2tic = ip("n2tic", 1);
    let grid2 = grid(&sp("grid2", "none"));
    let label1 = sp("label1", "");
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
    let frame = ip("frame", 1);
    let mut curvewidth = Vec::new();
    let mut curvedash = Vec::new();
    if curve != 0 {
        let nw = countparval("curvewidth");
        getparfloat("curvewidth", &mut curvewidth);
        if nw == 0 {
            curvewidth.push(axeswidth);
        }
        while curvewidth.len() < curve {
            curvewidth.push(*curvewidth.last().unwrap());
        }
        let nd = countparval("curvedash");
        getparint("curvedash", &mut curvedash);
        if nd == 0 {
            curvedash.push(0);
        }
        while curvedash.len() < curve {
            curvedash.push(*curvedash.last().unwrap());
        }
    }
    let style = if sp("style", "seismic") == "normal" {
        0
    } else {
        1
    };
    x2min = if x2beg < x2end { x2beg } else { x2end };
    x2max = if x2beg > x2end { x2beg } else { x2end };
    let mut n2in = 0;
    for &v in &x2 {
        if v >= x2min && v <= x2max {
            n2in += 1;
        }
    }
    let mut xcur = fp("xcur", 1.0).abs();
    if n2in > 1 {
        xcur *= (x2max - x2min) / (n2in - 1) as f32;
    }
    let p2beg = if x2end >= x2beg { -xcur } else { xcur };
    let p2end = if x2end >= x2beg { xcur } else { -xcur };
    checkpars();
    let i1beg = nint(((x1beg - f1) / d1) as f64).max(0).min(n1 - 1);
    x1beg = f1 + i1beg as f32 * d1;
    let i1end = nint(((x1end - f1) / d1) as f64).max(0).min(n1 - 1);
    x1end = f1 + i1end as f32 * d1;
    let if1p = i1beg.min(i1end) as usize;
    let n1p = (i1beg.max(i1end) - if1p as i32 + 1) as usize;
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
    let x1size = if style == 0 { wbox } else { hbox };
    let x2size = if style == 0 { hbox } else { wbox };
    let pscale = x2size / (x2end + p2end - x2beg - p2beg);
    let poffset = -(x2beg + p2beg) * pscale;
    let pxcur = xcur * pscale;
    let p1fz = if x1end > x1beg { 0.0 } else { x1size };
    let p1lz = if x1end > x1beg { x1size } else { 0.0 };
    gsave();
    translate(xbox, ybox);
    if backcolor != "none" {
        gsave();
        setcolor(&backcolor);
        newpath();
        moveto(0.0, 0.0);
        lineto(wbox, 0.0);
        lineto(wbox, hbox);
        lineto(0.0, hbox);
        lineto(0.0, 0.0);
        printf!("fill\n");
        grestore();
    }
    if style != 0 {
        rotate(-90.0);
        translate(-hbox, 0.0);
    }
    setcolor(&tracecolor);
    setlinewidth(linewidth);
    for i in 0..n2 as usize {
        if x2[i] < x2min || x2[i] > x2max {
            continue;
        }
        let px2 = poffset + x2[i] * pscale;
        let start = i * n1 as usize + if1p;
        psWiggle(
            &z[start..start + n1p],
            -clip,
            clip,
            0.0,
            px2 - pxcur,
            px2 + pxcur,
            p1fz,
            p1lz,
            fill,
            &tracecolor,
        );
    }
    grestore();
    for i in 0..curve {
        gsave();
        psDrawCurve(
            xbox,
            ybox,
            wbox,
            hbox,
            x1beg,
            x1end,
            0.0,
            0.0,
            x2beg,
            x2end,
            p2beg,
            p2end,
            &x1curve[i],
            &x2curve[i],
            &curvecolor[i],
            curvewidth[i],
            curvedash[i],
            style,
        );
        grestore();
    }
    if frame != 0 {
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
            p2beg,
            p2end,
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
    }
    showpage();
    endeps();
    exit(CWP_Exit());
}
