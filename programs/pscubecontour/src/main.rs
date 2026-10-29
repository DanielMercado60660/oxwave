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

//! Port of psplot/main/pscubecontour.c. Original pscube and pscontour by
//! Craig Artley, Zhiming Li and Dave Hale; merged by Claudia Vanelle,
//! Institute of Geophysics, University of Hamburg, 2000.

use oxwave_core::prelude::*;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;
use oxwave_f_pscaxesbox::*;
use oxwave_f_pscontour::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSCCONTOUR - PostScript Contour plot of a data CUBE\t\t        ",
    " \t\t\t\t\t\t\t\t\t",
    " pscubecontour n1= n2= n3= [optional parameters] <binaryfile >postscriptfile\t",
    "    or\t\t\t\t\t\t\t\t\t",
    " pscubecontour n1= n2= n3= front= side= top= [optional parameters] >postscriptfile",
    " \t\t\t\t\t\t\t\t\t",
    " Data formats supported:\t\t\t\t\t\t",
    "\t1. Entire cube read from stdin (n1*n2*n3 floats) [default format]",
    "\t2. Faces read from stdin (n1*n2 floats for front, followed by n1*n3",
    "\t   floats for side, and n2*n3 floats for top) [specify faces=1]\t",
    "\t3. Faces read from separate data files [specify filenames]\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1                     number of samples in 1st (fastest) dimension\t",
    " n2                     number of samples in 2nd dimension\t\t",
    " n3                     number of samples in 3rd (slowest) dimension\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " front                  name of file containing front panel\t\t",
    " side                   name of file containing side panel\t\t",
    " top                    name of file containing top panel\t\t",
    " faces=0                =1 to read faces from stdin (data format 2)\t",
    " d1=1.0                 sampling interval in 1st dimension\t\t",
    " f1=0.0                 first sample in 1st dimension\t\t\t",
    " d2=1.0                 sampling interval in 2nd dimension\t\t",
    " f2=0.0                 first sample in 2nd dimension\t\t\t",
    " d3=1.0                 sampling interval in 3rd dimension\t\t",
    " f3=0.0                 first sample in 3rd dimension\t\t\t",
    " d1s=1.0                factor by which to scale d1 before imaging\t",
    " d2s=1.0                factor by which to scale d2 before imaging\t",
    " d3s=1.0                factor by which to scale d3 before imaging\t",
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
    " size1=4.0              size in inches of 1st axes (vertical)\t\t",
    " size2=4.0              size in inches of 2nd axes (horizontal)\t",
    " size3=3.0              size in inches of 3rd axes (projected)\t\t",
    " angle=45               projection angle of cube in degrees (0<angle<90)",
    "                        (angle between 2nd axis and projected 3rd axis)",
    " x1end=x1max            value at which axis 1 ends\t\t\t",
    " d1num=0.0              numbered tic interval on axis 1 (0.0 for automatic)",
    " f1num=x1min            first numbered tic on axis 1 (used if d1num not 0.0)",
    " n1tic=1                number of tics per numbered tic on axis 1\t",
    " grid1=none             grid lines on axis 1 - none, dot, dash, or solid",
    " label1=                label on axis 1\t\t\t\t",
    " x2beg=x2min            value at which axis 2 begins\t\t\t",
    " d2num=0.0              numbered tic interval on axis 2 (0.0 for automatic)",
    " f2num=x2min            first numbered tic on axis 2 (used if d2num not 0.0)",
    " n2tic=1                number of tics per numbered tic on axis 2\t",
    " grid2=none             grid lines on axis 2 - none, dot, dash, or solid",
    " label2=                label on axis 2\t\t\t\t",
    " x3end=x3max            value at which axis 3 ends\t\t\t",
    " d3num=0.0              numbered tic interval on axis 3 (0.0 for automatic)",
    " f3num=x3min            first numbered tic on axis 3 (used if d3num not 0.0)",
    " n3tic=1                number of tics per numbered tic on axis 3\t",
    " grid3=none             grid lines on axis 3 - none, dot, dash, or solid",
    " label3=                label on axis 3\t\t\t\t",
    " labelfont=Helvetica    font name for axes labels\t\t\t",
    " labelsize=18           font size for axes labels\t\t\t",
    " title=                 title of plot\t\t\t\t\t",
    " titlefont=Helvetica-Bold font name for title\t\t\t\t",
    " titlesize=24           font size for title\t\t\t\t",
    " titlecolor=black       color of title\t\t\t\t\t",
    " labelcfont=Helvetica-Bold font name for contour labels\t\t",
    " labelcsize=6           font size of contour labels   \t\t\t",
    " labelccolor=black      color of contour labels   \t\t\t",
    " axescolor=black        color of axes\t\t\t\t\t",
    " gridcolor=black        color of grid\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " All color specifications may also be made in X Window style Hex format",
    " example:   axescolor=#255\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: The values of x1beg=x1min, x2end=x2max and x3beg=x3min cannot   ",
    " be changed.\t\t\t\t\t\t\t\t",
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
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut n1 = 0;
    if getparint("n1", &mut n1) == 0 {
        err!("must specify n1!\n");
    }
    let mut d1 = 1.0;
    let mut f1 = 0.0;
    getparfloat("d1", &mut d1);
    getparfloat("f1", &mut f1);
    let mut n2 = 0;
    if getparint("n2", &mut n2) == 0 {
        err!("must specify n2!\n");
    }
    let mut d2 = 1.0;
    let mut f2 = 0.0;
    getparfloat("d2", &mut d2);
    getparfloat("f2", &mut f2);
    let mut n3 = 0;
    if getparint("n3", &mut n3) == 0 {
        err!("must specify n3!\n");
    }
    if n3 < 2 {
        err!("must have n3>=2!");
    }
    let mut d3 = 1.0;
    let mut f3 = 0.0;
    getparfloat("d3", &mut d3);
    getparfloat("f3", &mut f3);
    let mut faces = 0;
    getparint("faces", &mut faces);
    let (a, b, c) = (
        n1 as usize * n2 as usize,
        n1 as usize * n3 as usize,
        n2 as usize * n3 as usize,
    );
    let mut z = vec![0.0f32; a + b + c];
    let (mut front, mut side, mut top) = (String::new(), String::new(), String::new());
    let all = getparstring("front", &mut front) != 0
        && getparstring("side", &mut side) != 0
        && getparstring("top", &mut top) != 0;
    if all {
        for (path, start, len, kind) in [(&front, 0, a, 0), (&side, a, b, 1), (&top, a + b, c, 2)] {
            let fp = match fopen(path, "r") {
                Some(fp) => fp,
                None => match kind {
                    0 => err!("error opening front file!\n"),
                    1 => err!("error opening side file!\n"),
                    _ => err!("error opening top file!\n"),
                },
            };
            if fread(&mut z[start..start + len], &fp) != len {
                match kind {
                    0 => err!("error reading front file!\n"),
                    1 => err!("error reading side file!\n"),
                    _ => err!("error reading top file!\n"),
                }
            }
        }
    } else if getparstring("front", &mut front) != 0
        || getparstring("side", &mut side) != 0
        || getparstring("top", &mut top) != 0
    {
        err!("must specify all or none of face, side, and top!\n");
    } else if faces != 0 {
        let fp = stdin();
        if fread(&mut z[..a], &fp) != a {
            err!("error reading front from input!\n");
        }
        if fread(&mut z[a..a + b], &fp) != b {
            err!("error reading side from input!\n");
        }
        if fread(&mut z[a + b..], &fp) != c {
            err!("error reading top from input!\n");
        }
    } else {
        let fp = stdin();
        let mut row = vec![0.0f32; n1 as usize];
        for i3 in 0..n3 as usize {
            for i2 in 0..n2 as usize {
                if fread(&mut row, &fp) != n1 as usize {
                    err!("error reading cube from input!\n");
                }
                for i1 in 0..n1 as usize {
                    if i3 == 0 {
                        z[i1 + i2 * n1 as usize] = row[i1];
                    }
                    if i2 == n2 as usize - 1 {
                        z[a + i1 + i3 * n1 as usize] = row[i1];
                    }
                }
                z[a + b + i2 + i3 * n2 as usize] = row[0];
            }
        }
    }
    let axis = |n: i32, d: f32, f: f32| {
        let x: Vec<f32> = (0..n).map(|i| f + i as f32 * d).collect();
        let mut lo = x[0];
        let mut hi = x[0];
        for &v in &x {
            lo = if lo < v { lo } else { v };
            hi = if hi > v { hi } else { v };
        }
        (x, lo, hi)
    };
    let (mut x1, x1min, x1max) = axis(n1, d1, f1);
    let (mut x2, x2min, x2max) = axis(n2, d2, f2);
    let (mut x3, x3min, x3max) = axis(n3, d3, f3);
    let mut w = vec![0.0f32; z.len()];
    let (mut zmin, mut zmax) = (z[0], z[0]);
    for &v in &z {
        zmin = if zmin < v { zmin } else { v };
        zmax = if zmax > v { zmax } else { v };
    }
    let mut contours = vec![0.0f32; 200];
    let mut nc = getparfloat("c", &mut contours) as usize;
    if nc == 0 {
        nc = pi("nc", 5) as usize;
        let dc = pf("dc", (zmax - zmin) / nc as f32);
        let fc = pf("fc", zmin + dc);
        for i in 0..nc {
            contours[i] = fc + i as f32 * dc;
        }
    }
    let mut widths = vec![1.0f32; 200];
    let mut grays = vec![0.0f32; 200];
    let mut dashes = vec![0.0f32; 200];
    for (name, arr) in [
        ("cwidth", &mut widths),
        ("cgray", &mut grays),
        ("cdash", &mut dashes),
    ] {
        let n = getparfloat(name, arr) as usize;
        if n > 0 {
            for i in n..nc {
                arr[i] = arr[n - 1];
            }
        }
    }
    let color = ps("ccolor", "none");
    let mut colors: Vec<&str> = color.split(',').filter(|s| !s.is_empty()).collect();
    if colors.is_empty() {
        colors.push("none");
    }
    while colors.len() < nc {
        colors.push(*colors.last().unwrap());
    }
    let labelcf = pi("labelcf", 1);
    let labelcper = pi("labelcper", 1);
    let nlabelc = pi("nlabelc", nc as i32);
    let labelcsize = pf("labelcsize", 6.);
    let labelcfont = ps("labelcfont", "Helvetica-Bold");
    let labelccolor = ps("labelccolor", "black");
    let nplaces = pi("nplaces", 6);
    let size1 = pf("size1", 4.) * 72.;
    let size2 = pf("size2", 4.) * 72.;
    let size3 = pf("size3", 3.) * 72.;
    let xbox = pf("xbox", 1.5) * 72.;
    let ybox = pf("ybox", 1.5) * 72.;
    let angle = pf("angle", 45.).max(0.00001).min(90.);
    let angle = (angle as f64 * (std::f64::consts::PI / 180.)) as f32;
    let x1beg = x1min;
    let x1end = pf("x1end", x1max);
    let d1num = pf("d1num", 0.);
    let f1num = pf("f1num", x1min);
    let n1tic = pi("n1tic", 1);
    let grid1 = grid(&ps("grid1", "none"));
    let label1 = ps("label1", "");
    let x2beg = pf("x2beg", x2min);
    let x2end = x2max;
    let d2num = pf("d2num", 0.);
    let f2num = pf("f2num", x2min);
    let n2tic = pi("n2tic", 1);
    let grid2 = grid(&ps("grid2", "none"));
    let label2 = ps("label2", "");
    let x3beg = x3min;
    let x3end = pf("x3end", x3max);
    let d3num = pf("d3num", 0.);
    let f3num = pf("f3num", x3min);
    let n3tic = pi("n3tic", 1);
    let grid3 = grid(&ps("grid3", "none"));
    let label3 = ps("label3", "");
    let labelfont = ps("labelfont", "Helvetica");
    let labelsize = pf("labelsize", 18.);
    let title = ps("title", "");
    let titlefont = ps("titlefont", "Helvetica-Bold");
    let titlesize = pf("titlesize", 24.);
    let titlecolor = ps("titlecolor", "black");
    let axescolor = ps("axescolor", "black");
    let gridcolor = ps("gridcolor", "black");
    checkpars();
    let pads = |beg: f32, end: f32, d: f32| {
        let p = (d.abs() as f64 / 2.) as f32;
        if end > beg {
            (-p, p)
        } else {
            (p, -p)
        }
    };
    let (p1beg, p1end) = pads(x1beg, x1end, d1);
    let (p2beg, p2end) = pads(x2beg, x2end, d2);
    let (p3beg, p3end) = pads(x3beg, x3end, d3);
    let sin = (angle as f64).sin();
    let cos = (angle as f64).cos();
    let tan = (angle as f64).tan();
    let mut bbox = [0; 4];
    psAxesBBox(
        xbox,
        ybox,
        (size2 as f64 + cos * size3 as f64) as f32,
        (size1 as f64 + sin * size3 as f64) as f32,
        &labelfont,
        labelsize,
        &titlefont,
        titlesize,
        1,
        &mut bbox,
    );
    boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
    begineps();
    gsave();
    translate(xbox, ybox);
    let draw = |x: &[f32], y: &[f32], z: &mut [f32], w: &mut [f32], labels: bool| {
        for ic in 0..nc {
            setlinewidth(widths[ic]);
            if colors[ic] != "none" {
                setcolor(colors[ic]);
            } else {
                setgray(grays[ic]);
            }
            if dashes[ic] != 0. {
                setdash(&dashes[ic..ic + 1], 0.);
            } else {
                setdash(&[], 0.);
            }
            let mut lcsize = 0.;
            if nlabelc > 0 {
                if (ic as i32 - labelcf + 1) % labelcper == 0
                    && ic as i32 >= labelcf - 1
                    && (ic as i32) < labelcf - 1 + labelcper * nlabelc
                {
                    setlinewidth(widths[ic]);
                    lcsize = labelcsize;
                } else {
                    setlinewidth((0.25f64 * widths[ic] as f64) as f32);
                }
            }
            psContour(
                contours[ic],
                x,
                y,
                z,
                if labels { lcsize } else { 0. },
                &labelcfont,
                &labelccolor,
                w,
                nplaces,
            );
        }
    };
    gsave();
    let x2scale = size2 / (x2end - x2beg);
    let x3scale = (size3 as f64 * sin / (x3end - x3beg) as f64) as f32;
    for v in &mut x2 {
        *v *= x2scale;
    }
    for v in &mut x3 {
        *v *= x3scale;
    }
    let matrix = [1., 0., (1. / tan) as f32, 1., 0., 0.];
    translate(
        ((-(f2 - x2min + x2beg) * x2scale) as f64 - (f3 * x3scale) as f64 / tan) as f32,
        size1 - f3 * x3scale,
    );
    concat(&matrix);
    rectclip(
        x2beg * x2scale,
        f3 * x3scale,
        size2,
        (size3 as f64 * sin) as f32,
    );
    draw(&x2, &x3, &mut z[a + b..], &mut w[a + b..], false);
    for v in &mut x2 {
        *v /= x2scale;
    }
    for v in &mut x3 {
        *v /= x3scale;
    }
    translate(
        (((f2 - x2min + x2beg) * x2scale) as f64 + (f3 * x3scale) as f64 / tan) as f32,
        -size1 + f3 * x3scale,
    );
    grestore();
    gsave();
    let x1scale = size1 / (x1end - x1beg);
    let x2scale = size2 / (x2end - x2beg);
    for v in &mut x1 {
        *v *= x1scale;
    }
    for v in &mut x2 {
        *v *= x2scale;
    }
    rotate(-90.);
    translate(-size1 - f1 * x1scale, -(f2 - x2min + x2beg) * x2scale);
    rectclip(
        size1 + f1 * x1scale,
        (f2 - x2min + x2beg) * x2scale,
        -size1,
        size2,
    );
    draw(&x1, &x2, &mut z[..a], &mut w[..a], true);
    for v in &mut x1 {
        *v /= x1scale;
    }
    for v in &mut x2 {
        *v /= x2scale;
    }
    translate(size1 + f1 * x1scale, (f2 - x2min + x2beg) * x2scale);
    rotate(90.);
    grestore();
    gsave();
    let x1scale = size1 / (x1end - x1beg);
    let x3scale = (size3 as f64 * cos / (x3end - x3beg) as f64) as f32;
    for v in &mut x1 {
        *v *= x1scale;
    }
    for v in &mut x3 {
        *v *= x3scale;
    }
    let matrix = [1., tan as f32, 0., 1., 0., 0.];
    translate(
        size2 - f3 * x3scale,
        ((size1 + f1 * x1scale) as f64 - (f3 * x3scale) as f64 * tan) as f32,
    );
    concat(&matrix);
    rectclip(
        f3 * x3scale,
        -x1end * x1scale,
        (size3 as f64 * cos) as f32,
        size1,
    );
    rotate(-90.);
    draw(&x1, &x3, &mut z[a..a + b], &mut w[a..a + b], false);
    for v in &mut x1 {
        *v /= x1scale;
    }
    for v in &mut x3 {
        *v /= x3scale;
    }
    rotate(90.);
    translate(
        -size2 + f3 * x3scale,
        ((-size1 - f1 * x1scale) as f64 + (f3 * x3scale) as f64 * tan) as f32,
    );
    grestore();
    grestore();
    psCubeAxesBox(
        xbox,
        ybox,
        size1,
        size2,
        size3,
        angle,
        x1beg,
        x1end,
        p1beg,
        p1end,
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
        x3beg,
        x3end,
        p3beg,
        p3end,
        d3num,
        f3num,
        n3tic,
        grid3,
        &label3,
        &labelfont,
        labelsize,
        &title,
        &titlefont,
        titlesize,
        &titlecolor,
        &axescolor,
        &gridcolor,
    );
    showpage();
    endeps();
    exit(CWP_Exit());
}

fn pf(name: &str, mut value: f32) -> f32 {
    getparfloat(name, &mut value);
    value
}
fn pi(name: &str, mut value: i32) -> i32 {
    getparint(name, &mut value);
    value
}
fn ps(name: &str, value: &str) -> String {
    let mut s = value.to_string();
    getparstring(name, &mut s);
    s
}
fn grid(s: &str) -> i32 {
    match s {
        "dot" => 1,
        "dash" => 2,
        "solid" => 3,
        _ => 0,
    }
}
