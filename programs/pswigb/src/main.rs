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

//! PSWIGB - PostScript wiggle-trace bitmap plot, port of psplot/main/pswigb.c.
//! Author: Dave Hale; modifications by Craig Artley, Morten Wendell Pedersen, Brian K. Macy, P. Michaels, and Ekkehart Tessmer.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSWIGB - PostScript WIGgle-trace plot of f(x1,x2) via Bitmap\t\t",
    " Best for many traces.  Use PSWIGP (Polygon version) for few traces.\t",
    " \t\t\t\t\t\t\t\t\t",
    " pswigb n1= [optional parameters] <binaryfile >postscriptfile\t\t",
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
    " wt=1                   =0 for no wiggle-trace; =1 for wiggle-trace\t",
    " va=1                   =0 for no variable-area; =1 for variable-area fill",
    "                        =2 for variable area, solid/grey fill          ",
    "                        SHADING: 2<= va <=5  va=2 lightgrey, va=5 black",
    " nbpi=72                number of bits per inch at which to rasterize\t",
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
    " style=seismic          normal (axis 1 horizontal, axis 2 vertical) or\t",
    "                        seismic (axis 1 vertical, axis 2 horizontal)\t",
    " interp=0\t\t no display interpolation\t\t\t",
    "\t\t\t =1 use 8 point sinc interpolation\t\t",
    " curve=curve1,curve2,...  file(s) containing points to draw curve(s)   ",
    " npair=n1,n2,n2,...            number(s) of pairs in each file         ",
    " curvecolor=black,..    color of curve(s)                              ",
    " curvewidth=axeswidth   width (in points) of curve(s)                  ",
    " curvedash=0            solid curve(s), dash indices 1,...,11 produce  ",
    "                        curve(s) with various dash styles              ",
    " \t\t\t\t\t\t\t\t\t",
    " Notes: \t\t\t\t\t\t\t\t",
    " The interp option may be useful for high nbpi values, however, it\t",
    " tacitly assumes that the data are purely oscillatory.\tNon-oscillatory\t",
    " data will not be represented correctly when this option is set.\t",
    " \t\t\t\t\t\t\t\t\t",
    " The curve file is an ascii file with the points specified as x1 x2\t",
    " pairs, one pair to a line.  A \"vector\" of curve files and curve\t",
    " colors may be specified as curvefile=file1,file2,etc. and \t\t",
    " curvecolor=color1,color2,etc, and the number of pairs of values in each",
    " file as npair=npair1,npair2,... .\t\t\t\t\t",
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
use oxwave_sort::qkfind;
use oxwave_wiggle_raster::*;
fn pf(n: &str, d: f32) -> f32 {
    let mut v = d;
    getparfloat(n, &mut v);
    v
}
fn pi(n: &str, d: i32) -> i32 {
    let mut v = d;
    getparint(n, &mut v);
    v
}
fn ps(n: &str, d: &str) -> String {
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
    let d1 = pf("d1", 1.0);
    let f1 = pf("f1", 0.0);
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
    let mut n2 = countparval("x2");
    let infp = stdin();
    if n2 == 0 && getparint("n2", &mut n2) == 0 {
        if efseeko(&infp, 0, SEEK_END) != 0 {
            err!("must specify n2 if in a pipe!\n");
        }
        let nf = (eftello(&infp) / 4) as i32;
        efseeko(&infp, 0, SEEK_SET);
        n2 = nf.checked_div(n1).unwrap_or(0);
    }
    let mut x2 = vec![0.0; n2 as usize];
    if getparfloat("x2", &mut x2) == 0 {
        let d2 = pf("d2", 1.0);
        let f2 = pf("f2", 0.0);
        for i in 0..n2 as usize {
            x2[i] = f2 + i as f32 * d2;
        }
    }
    let mut x2min = x2[0];
    let mut x2max = x2[0];
    for &v in &x2[1..] {
        if x2min >= v {
            x2min = v;
        }
        if x2max <= v {
            x2max = v;
        }
    }
    let curve = countparval("curve") as usize;
    let mut files = Vec::new();
    let mut npair = vec![0; curve];
    let mut colors = Vec::new();
    let mut curves = Vec::new();
    if curve > 0 {
        getparstringarray("curve", &mut files);
        getparint("npair", &mut npair);
        let nc = countparval("curvecolor") as usize;
        getparstringarray("curvecolor", &mut colors);
        if nc == 0 {
            colors.push("black".to_string());
        }
        while colors.len() < curve {
            colors.push(colors.last().unwrap().clone());
        }
        for j in 0..curve {
            let fp = efopen(&files[j], "r");
            let mut a = vec![0.0; npair[j] as usize];
            let mut b = a.clone();
            for i in 0..a.len() {
                a[i] = fscanf_f32(&fp).unwrap_or(0.0);
                b[i] = fscanf_f32(&fp).unwrap_or(0.0);
            }
            efclose(&fp);
            curves.push((a, b));
        }
    }
    let nz = n1 * n2;
    let mut z = vec![0.0f32; nz as usize];
    if fread(&mut z, &infp) != nz as usize {
        err!("error reading input file!\n");
    }
    let mut bias = 0.0;
    if getparfloat("bias", &mut bias) != 0 && bias != 0.0 {
        for v in &mut z {
            *v -= bias;
        }
    }
    let mut clip = 0.0;
    if getparfloat("clip", &mut clip) == 0 {
        let perc = pf("perc", 100.0);
        let mut temp: Vec<f32> = z.iter().map(|v| v.abs()).collect();
        let iz = (((nz as f32 * perc) as f64 / 100.0) as i32)
            .max(0)
            .min(nz - 1) as usize;
        qkfind(iz, &mut temp);
        clip = temp[iz];
    }
    let verbose = pi("verbose", 1);
    if verbose != 0 {
        warn!("clip=%g\n", clip);
    }
    let mut wt = pi("wt", 1);
    let mut va = pi("va", 1);
    if va >= 2 {
        wt = va;
        va = 1;
    }
    let mut xbox = pf("xbox", 1.5);
    let mut ybox = pf("ybox", 1.5);
    let mut wbox = pf("wbox", 6.0);
    let mut hbox = pf("hbox", 8.0);
    let mut x1beg = pf("x1beg", x1min);
    let mut x1end = pf("x1end", x1max);
    let d1num = pf("d1num", 0.0);
    let f1num = pf("f1num", x1min);
    let n1tic = pi("n1tic", 1);
    let grid1 = grid(&ps("grid1", "none"));
    let label1 = ps("label1", "");
    let x2beg = pf("x2beg", x2min);
    let x2end = pf("x2end", x2max);
    let d2num = pf("d2num", 0.0);
    let f2num = pf("f2num", 0.0);
    let n2tic = pi("n2tic", 1);
    let grid2 = grid(&ps("grid2", "none"));
    let label2 = ps("label2", "");
    let labelfont = ps("labelfont", "Helvetica");
    let labelsize = pf("labelsize", 18.0);
    let title = ps("title", "");
    let titlefont = ps("titlefont", "Helvetica-Bold");
    let titlesize = pf("titlesize", 24.0);
    let titlecolor = ps("titlecolor", "black");
    let axescolor = ps("axescolor", "black");
    let gridcolor = ps("gridcolor", "black");
    let axeswidth = pf("axeswidth", 1.0);
    let ticwidth = pf("ticwidth", axeswidth);
    let gridwidth = pf("gridwidth", axeswidth);
    let mut widths = vec![0.0; curve];
    let mut dashes = vec![0; curve];
    if curve > 0 {
        let nw = countparval("curvewidth") as usize;
        if getparfloat("curvewidth", &mut widths) == 0 {
            widths[0] = axeswidth;
        }
        for i in nw.max(1)..curve {
            widths[i] = widths[nw.max(1) - 1];
        }
        let nd = countparval("curvedash") as usize;
        getparint("curvedash", &mut dashes);
        for i in nd.max(1)..curve {
            dashes[i] = dashes[nd.max(1) - 1];
        }
    }
    let style = if ps("style", "seismic") == "normal" {
        0
    } else {
        1
    };
    let nbpi = pi("nbpi", 72);
    let n1bits = (nbpi as f32 * if style == 0 { wbox } else { hbox }) as i32;
    let n2bits = (nbpi as f32 * if style == 0 { hbox } else { wbox }) as i32;
    let nbpr = 1 + (n2bits - 1) / 8;
    let bitmap_len = (nbpr * n1bits) as usize;
    // Grey-fill stepping can advance one byte beyond the final raster row.
    // Keep that byte out of the emitted image while allowing the C behavior.
    let mut bits = vec![0u8; bitmap_len + 1];
    x2min = if x2beg < x2end { x2beg } else { x2end };
    x2max = if x2beg > x2end { x2beg } else { x2end };
    let n2in = x2.iter().filter(|&&x| x >= x2min && x <= x2max).count();
    let mut xcur = pf("xcur", 1.0).abs();
    if n2in > 1 {
        xcur *= (x2max - x2min) / (n2in - 1) as f32;
    }
    let p2beg = if x2end >= x2beg { -xcur } else { xcur };
    let p2end = if x2end >= x2beg { xcur } else { -xcur };
    let interp = pi("interp", 0);
    checkpars();
    let bscale = (n2bits - 1) as f32 / (x2end + p2end - x2beg - p2beg);
    let boffset = -(x2beg + p2beg) * bscale;
    let bxcur = xcur * bscale;
    let i1beg = nint(((x1beg - f1) / d1) as f64).max(0).min(n1 - 1);
    x1beg = f1 + i1beg as f32 * d1;
    let i1end = nint(((x1end - f1) / d1) as f64).max(0).min(n1 - 1);
    x1end = f1 + i1end as f32 * d1;
    let if1r = i1beg.min(i1end);
    let n1r = i1beg.max(i1end) - if1r + 1;
    let b1fz = if x1end > x1beg { 0 } else { n1bits - 1 };
    let b1lz = if x1end > x1beg { n1bits - 1 } else { 0 };
    for i in 0..n2 as usize {
        if x2[i] < x2min || x2[i] > x2max {
            continue;
        }
        let bx2 = boffset + x2[i] * bscale;
        let samples =
            &z[(i * n1 as usize + if1r as usize)..(i * n1 as usize + (if1r + n1r) as usize)];
        let fun = if interp == 0 { rfwtva } else { rfwtvaint };
        fun(
            samples,
            -clip,
            clip,
            if va != 0 { 0.0 } else { clip },
            (bx2 - bxcur) as i32,
            (bx2 + bxcur) as i32,
            b1fz,
            b1lz,
            wt,
            nbpr,
            &mut bits,
            0,
        );
    }
    bits.truncate(bitmap_len);
    for b in &mut bits {
        *b = !b.reverse_bits();
    }
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
    translate(xbox, ybox);
    let matrix = if style == 0 {
        [0.0, n1bits as f32, n2bits as f32, 0.0, 0.0, 0.0]
    } else {
        [
            n2bits as f32,
            0.0,
            0.0,
            -(n1bits as f32),
            0.0,
            n1bits as f32,
        ]
    };
    scale(wbox, hbox);
    image(n2bits, n1bits, 1, &matrix, &bits);
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
            &curves[i].0,
            &curves[i].1,
            &colors[i],
            widths[i],
            dashes[i],
            style,
        );
        grestore();
    }
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
    showpage();
    endeps();
    exit(CWP_Exit());
}
