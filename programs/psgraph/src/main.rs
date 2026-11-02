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
//! PSGRAPH, port of psplot/main/psgraph.c, revision 1.40.
//! Dave Hale; modifications by Jack K. Cohen, Lydia Deng, Craig Artley,
//! John Stockwell, Morten Wendell Pedersen (Aarhus University), and James Gunning (CSIRO).

use oxwave_core::prelude::*;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " PSGRAPH - PostScript GRAPHer\t\t\t\t\t\t",
    " Graphs n[i] pairs of (x,y) coordinates, for i = 1 to nplot.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " psgraph n= [optional parameters] <binaryfile >postscriptfile\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n                      array containing number of points per plot\t",
    "\t\t\t\t\t\t\t\t\t",
    " Data formats supported:\t\t\t\t\t\t",
    "\t1.a. x1,y1,x2,y2,...,xn,yn\t\t\t\t\t",
    "\t  b. x1,x2,...,xn,y1,y2,...,yn (must set pairs=0)\t\t",
    "\t2.   y1,y2,...,yn (must give non-zero d1[]=)\t\t\t",
    "\t3.   x1,x2,...,xn (must give non-zero d2[]=)\t\t\t",
    "\t4.   nil (must give non-zero d1[]= and non-zero d2[]=)\t\t",
    "  The formats may be repeated and mixed in any order, but if\t\t",
    "  formats 2-4 are used, the d1 and d2 arrays must be specified including",
    "  d1[]=0.0 d2[]=0.0 entries for any internal occurences of format 1.\t",
    "  Similarly, the pairs array must contain place-keeping entries for\t",
    "  plots of formats 2-4 if they are mixed with both formats 1.a and 1.b.",
    "  Also, if formats 2-4 are used with non-zero f1[] or f2[] entries, then",
    "  the corresponding array(s) must be fully specified including f1[]=0.0",
    "  and/or f2[]=0.0 entries for any internal occurences of format 1 or\t",
    "  formats 2-4 where the zero entries are desired.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  Available colors are all the common ones and many more. The complete\t",
    "  list of 68 colors is in the file $CWPROOT/src/psplot/basic.c.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nplot=number of n's    number of plots\t\t\t\t",
    " d1=0.0,...             x sampling intervals (0.0 if x coordinates input)",
    " f1=0.0,...             first x values (not used if x coordinates input)",
    " d2=0.0,...             y sampling intervals (0.0 if y coordinates input)",
    " f2=0.0,...             first y values (not used if y coordinates input)",
    " pairs=1,...            =1 for data pairs in format 1.a, =0 for format 1.b",
    " linewidth=1.0,...      line widths (in points) (0.0 for no lines)\t",
    " linegray=0.0,...       line gray levels (black=0.0 to white=1.0)\t",
    " linecolor=none,...     line colors; none means use linegray\t\t",
    "                        Typical use: linecolor=red,yellow,blue,...\t",
    " lineon=1.0,...         length of line segments for dashed lines (in points)",
    " lineoff=0.0,...        spacing between dashes (0.0 for solid line)\t",
    " mark=0,1,2,3,...       indices of marks used to represent plotted points",
    " marksize=0.0,0.0,...   size of marks (0.0 for no marks)\t\t",
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
    " ticwidth=axeswidth     width (in points) of tic marks\t\t\t",
    " gridwidth=axeswidth    width (in points) of grid lines\t\t",
    " style=normal           normal (axis 1 horizontal, axis 2 vertical) or\t",
    "                        seismic (axis 1 vertical, axis 2 horizontal)\t",
    " reverse=0              =1 to reverse sequence of plotting curves      ",
    " Note:\tn1 and n2 are acceptable aliases for n and nplot, respectively.\t",
    "\t\t\t\t\t\t\t\t\t",
    " mark index:                                                           ",
    " 1. asterisk                                                           ",
    " 2. x-cross                                                            ",
    " 3. open triangle                                                      ",
    " 4. open square                                                        ",
    " 5. open circle                                                        ",
    " 6. solid triangle                                                     ",
    " 7. solid square                                                       ",
    " 8. solid circle                                                       ",
    "\t\t\t\t\t\t\t\t\t",
    " All color specifications may also be made in X Window style Hex format",
    " example:   axescolor=#255\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " psgraph n=50,100,20 d1=2.5,1,0.33 <datafile >psfile\t\t\t",
    "  plots three curves with equally spaced x values in one plot frame\t",
    "  x1-coordinates are x1(i) = f1+i*d1 for i = 1 to n (f1=0 by default)\t",
    "  number of x2's and then x2-coordinates for each curve are read\t",
    "  sequentially from datafile.\t\t\t\t\t\t",
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
];

fn floats(name: &str, default: f32, n: usize) -> Vec<f32> {
    let mut v = vec![default; 5000];
    let k = getparfloat(name, &mut v) as usize;
    if k > 0 {
        for i in k..n {
            v[i] = v[k - 1];
        }
    }
    v
}
fn ints(name: &str, default: i32, n: usize) -> Vec<i32> {
    let mut v = vec![default; 5000];
    let k = getparint(name, &mut v) as usize;
    if k > 0 {
        for i in k..n {
            v[i] = v[k - 1];
        }
    }
    v
}
fn string(name: &str, default: &str) -> String {
    let mut v = default.to_owned();
    getparstring(name, &mut v);
    v
}
fn float(name: &str, default: f32) -> f32 {
    let mut v = default;
    getparfloat(name, &mut v);
    v
}
fn int(name: &str, default: i32) -> i32 {
    let mut v = default;
    getparint(name, &mut v);
    v
}
fn read_float() -> f32 {
    let mut v = [0.0f32];
    if efread(&mut v, &stdin()) != 1 {
        err!("Error reading input!\n");
    }
    v[0]
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
    let mut pairs = ints("pairs", 1, 0);
    let mut d1 = floats("d1", 0.0, 0);
    let mut f1 = floats("f1", 0.0, 0);
    let mut d2 = floats("d2", 0.0, 0);
    let mut f2 = floats("f2", 0.0, 0);
    let reverse = int("reverse", 0);
    let mut n = vec![0i32; 5000];
    let mut nn = getparint("n", &mut n);
    if nn == 0 {
        nn = getparint("n1", &mut n);
    }
    if nn == 0 {
        err!("Must specify n, the number of points per plot!");
    }
    let mut nplot = nn;
    getparint("n2", &mut nplot);
    getparint("nplot", &mut nplot);
    if nplot > 5000 {
        err!("too many plots");
    }
    // ealloc1 receives nplot as an unsigned allocation size in C.
    if nplot < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let np = nplot.max(0) as usize;
    for i in nn as usize..np {
        n[i] = n[nn as usize - 1];
    }
    let k = countparval("pairs") as usize;
    if k > 0 {
        for i in k..np {
            pairs[i] = pairs[k - 1];
        }
    }
    for (name, v) in [
        ("d1", &mut d1),
        ("f1", &mut f1),
        ("d2", &mut d2),
        ("f2", &mut f2),
    ] {
        let k = countparval(name) as usize;
        if k > 0 {
            for i in k..np {
                v[i] = v[k - 1];
            }
        }
    }
    let (mut x1min, mut x1max, mut x2min, mut x2max) = (f32::MAX, -f32::MAX, f32::MAX, -f32::MAX);
    let mut x1data = Vec::new();
    let mut x2data = Vec::new();
    for p in 0..np {
        // ealloc1float converts a negative npoint to an allocation size.
        if n[p] < 0 {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        let len = n[p] as usize;
        let mut x = vec![0.0; len];
        let mut y = vec![0.0; len];
        if d1[p] == 0.0 && d2[p] == 0.0 && pairs[p] == 0 {
            for a in &mut x {
                *a = read_float();
            }
            for a in &mut y {
                *a = read_float();
            }
        } else {
            for i in 0..len {
                x[i] = if d1[p] != 0.0 {
                    f1[p] + i as f32 * d1[p]
                } else {
                    read_float()
                };
                y[i] = if d2[p] != 0.0 {
                    f2[p] + i as f32 * d2[p]
                } else {
                    read_float()
                };
            }
        }
        for i in 0..len {
            if x[i] > x1max {
                x1max = x[i];
            }
            if x[i] < x1min {
                x1min = x[i];
            }
            if y[i] > x2max {
                x2max = y[i];
            }
            if y[i] < x2min {
                x2min = y[i];
            }
        }
        x1data.push(x);
        x2data.push(y);
    }
    if x1min == f32::MAX {
        x1min = 0.0;
        x1max = 0.0;
    }
    if x2min == f32::MAX {
        x2min = 0.0;
        x2max = 0.0;
    }
    if x1min == x1max {
        x1min -= 1.0;
        x1max += 1.0;
    }
    if x2min == x2max {
        x2min -= 1.0;
        x2max += 1.0;
    }
    let label1 = string("label1", "");
    let label2 = string("label2", "");
    let title = string("title", "");
    let style = if string("style", "normal") == "normal" {
        0
    } else {
        1
    };
    let labelfont = string("labelfont", "Helvetica");
    let titlefont = string("titlefont", "Helvetica-Bold");
    let labelsize = float("labelsize", 18.0);
    let titlesize = float("titlesize", 24.0);
    let titlecolor = string("titlecolor", "black");
    let axescolor = string("axescolor", "black");
    let gridcolor = string("gridcolor", "black");
    let mut xbox = float("xbox", 1.5);
    let mut ybox = float("ybox", 1.5);
    let mut wbox = float("wbox", 6.0);
    let mut hbox = float("hbox", 8.0);
    let x1beg = float("x1beg", x1min);
    let x1end = float("x1end", x1max);
    let x2beg = float("x2beg", x2min);
    let x2end = float("x2end", x2max);
    let d1num = float("d1num", 0.0);
    let f1num = float("f1num", x1min);
    let n1tic = int("n1tic", 1);
    let grid1 = grid(&string("grid1", "none"));
    let d2num = float("d2num", 0.0);
    let f2num = float("f2num", 0.0);
    let n2tic = int("n2tic", 1);
    let grid2 = grid(&string("grid2", "none"));
    let mut mark: Vec<i32> = (0..5000).map(|i| i % 9).collect();
    let k = getparint("mark", &mut mark) as usize;
    if k > 0 {
        for i in k..np {
            mark[i] = mark[k - 1];
        }
    }
    let marksize = floats("marksize", 0.0, np);
    let linewidth = floats("linewidth", 1.0, np);
    let lineon = floats("lineon", 1.0, np);
    let lineoff = floats("lineoff", 0.0, np);
    let linegray = floats("linegray", 0.0, np);
    let colors = string("linecolor", "none");
    let parts: Vec<&str> = colors.split(',').filter(|s| !s.is_empty()).collect();
    let axeswidth = float("axeswidth", 1.0);
    let ticwidth = float("ticwidth", axeswidth);
    let gridwidth = float("gridwidth", axeswidth);
    checkpars();
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
    let (xsize, ysize) = if style == 0 {
        (wbox, hbox)
    } else {
        (hbox, wbox)
    };
    translate(xbox, ybox);
    if style != 0 {
        rotate(-90.0);
        translate(-hbox, 0.0);
    }
    setlinejoin(2);
    let xscale = xsize / (x1end - x1beg);
    let yscale = ysize / (x2end - x2beg);
    for j in 0..np {
        let p = if reverse == 1 { np - 1 - j } else { j };
        let x: Vec<f32> = x1data[p].iter().map(|v| (*v - x1beg) * xscale).collect();
        let y: Vec<f32> = x2data[p].iter().map(|v| (*v - x2beg) * yscale).collect();
        let color = parts
            .get(p)
            .or_else(|| parts.last())
            .copied()
            .unwrap_or("none");
        gsave();
        if linewidth[p] != 0.0 {
            setlinewidth(linewidth[p]);
            if lineoff[p] != 0.0 {
                setdash(&[lineon[p], lineoff[p]], 0.0);
            } else {
                setdash(&[], 0.0);
            }
            if color != "none" {
                setcolor(color)
            } else {
                setgray(linegray[p])
            }
            if x.is_empty() {
                polyline(&[0.0], &[0.0]);
            } else {
                polyline(&x, &y);
            }
        }
        if marksize[p] != 0.0 {
            if color != "none" {
                setcolor(color)
            } else {
                setgray(linegray[p])
            }
            for i in 0..x.len() {
                markto(x[i], y[i], mark[p], marksize[p]);
            }
        }
        grestore();
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
