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

//! PSMOVIE - PostScript movie plot. Port of psplot/main/psmovie.c, revision 1.16.
//! Authors: Dave Hale and Zhiming Li; modified by Craig Artley.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " PSMOVIE - PostScript MOVIE plot of a uniformly-sampled function f(x1,x2,x3)",
    " \t\t\t\t\t\t\t\t\t",
    " psmovie n1= [optional parameters] <binaryfile >postscriptfile\t\t",
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
    " perc=100.0             percentile used to determine clip\t\t",
    " clip=(perc percentile) clip used to determine bclip and wclip\t\t",
    " bperc=perc             percentile for determining black clip value\t",
    " wperc=100.0-perc       percentile for determining white clip value\t",
    " bclip=clip             data values outside of [bclip,wclip] are clipped",
    " wclip=-clip            data values outside of [bclip,wclip] are clipped",
    " d1s=1.0                factor by which to scale d1 before imaging\t",
    " d2s=1.0                factor by which to scale d2 before imaging\t",
    " verbose=1              =1 for info printed on stderr (0 for no info)\t",
    " xbox=1.0               offset in inches of left side of axes box\t",
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
    " style=seismic          normal (axis 1 horizontal, axis 2 vertical) or\t",
    "                        seismic (axis 1 vertical, axis 2 horizontal)\t",
    " n3=1                   number of samples in third dimension\t\t",
    " title2=                second title to annotate different frames\t",
    " loopdsp=3              display loop type (1=loop over n1; 2=loop over n2;",
    "                                           3=loop over n3)\t\t",
    " d3=1.0                 sampling interval in 3rd dimension\t\t",
    " f3=d3                  first sample in 3rd dimension\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " NeXT: view movie via:   psmovie < infile n1= [optional params...] | open",
    " Note: currently only the Preview Application can handle the multipage  ",
    "       PostScript output by this program.\t\t\t\t",
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
use oxwave_byte_interp::intl2b;
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox3::*;
use oxwave_sort::qkfind;
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
    let infp = stdin();
    let mut n1 = 0;
    if getparint("n1", &mut n1) == 0 {
        err!("must specify n1!\n");
    }
    let mut d1 = 1.0f32;
    getparfloat("d1", &mut d1);
    let mut f1 = 0.0f32;
    getparfloat("f1", &mut f1);
    let mut n2 = 0;
    if getparint("n2", &mut n2) == 0 {
        efseeko(&infp, 0, SEEK_END);
        let nfloats = (eftello(&infp) / 4) as i32;
        efseeko(&infp, 0, SEEK_SET);
        n2 = nfloats.checked_div(n1).unwrap_or(0);
    }
    let mut d2 = 1.0f32;
    getparfloat("d2", &mut d2);
    let mut f2 = 0.0f32;
    getparfloat("f2", &mut f2);
    let mut n3 = 0;
    if getparint("n3", &mut n3) == 0 {
        efseeko(&infp, 0, SEEK_END);
        let nfloats = (eftello(&infp) / 4) as i32;
        efseeko(&infp, 0, SEEK_SET);
        n3 = nfloats
            .checked_div(n1)
            .unwrap_or(0)
            .checked_div(n2)
            .unwrap_or(0);
    }
    let mut d3 = 1.0f32;
    getparfloat("d3", &mut d3);
    let mut f3 = d3;
    getparfloat("f3", &mut f3);
    let mut loopdsp = 3;
    getparint("loopdsp", &mut loopdsp);
    let (n1dsp, n2dsp, n3dsp) = match loopdsp {
        3 => (n1, n2, n3),
        1 => (n2, n3, n1),
        2 => (n1, n3, n2),
        _ => (0, 0, 0),
    };
    let x1min = if d1 > 0.0 {
        f1
    } else {
        f1 + (n1dsp - 1) as f32 * d1
    };
    let x1max = if d1 < 0.0 {
        f1
    } else {
        f1 + (n1dsp - 1) as f32 * d1
    };
    let x2min = if d2 > 0.0 {
        f2
    } else {
        f2 + (n2dsp - 1) as f32 * d2
    };
    let x2max = if d2 < 0.0 {
        f2
    } else {
        f2 + (n2dsp - 1) as f32 * d2
    };
    let nz = n1 * n2 * n3;
    let mut z = vec![0.0f32; nz as usize];
    if fread(&mut z, &infp) != nz as usize {
        err!("error reading input file!\n");
    }
    let mut clip = 0.0f32;
    let mut bclip = 0.0f32;
    let mut wclip = 0.0f32;
    if getparfloat("clip", &mut clip) != 0 {
        bclip = clip;
        wclip = -clip;
    }
    if (getparfloat("bclip", &mut bclip) == 0 || getparfloat("wclip", &mut wclip) == 0)
        && getparfloat("clip", &mut clip) == 0
    {
        let mut perc = 100.0f32;
        getparfloat("perc", &mut perc);
        let mut temp = z.clone();
        if getparfloat("bclip", &mut bclip) == 0 {
            let mut bperc = perc;
            getparfloat("bperc", &mut bperc);
            let iz = (((nz as f32 * bperc) as f64 / 100.0) as i32)
                .max(0)
                .min(nz - 1) as usize;
            qkfind(iz, &mut temp);
            bclip = temp[iz];
        }
        if getparfloat("wclip", &mut wclip) == 0 {
            let mut wperc = (100.0 - perc as f64) as f32;
            getparfloat("wperc", &mut wperc);
            let iz = (((nz as f32 * wperc) as f64 / 100.0) as i32)
                .max(0)
                .min(nz - 1) as usize;
            qkfind(iz, &mut temp);
            wclip = temp[iz];
        }
    }
    if loopdsp == 1 || loopdsp == 2 {
        let temp = z.clone();
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let dst = if loopdsp == 1 {
                        i2 + i3 * n2 + i1 * n2 * n3
                    } else {
                        i1 + i3 * n1 + i2 * n1 * n3
                    };
                    z[dst as usize] = temp[(i1 + i2 * n1 + i3 * n1 * n2) as usize];
                }
            }
        }
    }
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    if verbose != 0 {
        warn!("bclip=%g wclip=%g", bclip, wclip);
    }
    let mut d1s = 1.0f32;
    getparfloat("d1s", &mut d1s);
    d1s = d1s.abs() * d1;
    let mut d2s = 1.0f32;
    getparfloat("d2s", &mut d2s);
    d2s = d2s.abs() * d2;
    let mut xbox = 1.0f32;
    getparfloat("xbox", &mut xbox);
    let mut ybox = 1.5f32;
    getparfloat("ybox", &mut ybox);
    let mut wbox = 6.0f32;
    getparfloat("wbox", &mut wbox);
    let mut hbox = 8.0f32;
    getparfloat("hbox", &mut hbox);
    let mut x1beg = x1min;
    getparfloat("x1beg", &mut x1beg);
    let mut x1end = x1max;
    getparfloat("x1end", &mut x1end);
    let mut d1num = 0.0f32;
    getparfloat("d1num", &mut d1num);
    let mut f1num = x1min;
    getparfloat("f1num", &mut f1num);
    let mut n1tic = 1;
    getparint("n1tic", &mut n1tic);
    let mut grid1s = "none".to_string();
    getparstring("grid1", &mut grid1s);
    let grid1 = grid(&grid1s);
    let mut label1 = String::new();
    getparstring("label1", &mut label1);
    let mut x2beg = x2min;
    getparfloat("x2beg", &mut x2beg);
    let mut x2end = x2max;
    getparfloat("x2end", &mut x2end);
    let mut d2num = 0.0f32;
    getparfloat("d2num", &mut d2num);
    let mut f2num = 0.0f32;
    getparfloat("f2num", &mut f2num);
    let mut n2tic = 1;
    getparint("n2tic", &mut n2tic);
    let mut grid2s = "none".to_string();
    getparstring("grid2", &mut grid2s);
    let grid2 = grid(&grid2s);
    let mut label2 = String::new();
    getparstring("label2", &mut label2);
    let mut labelfont = "Helvetica".to_string();
    getparstring("labelfont", &mut labelfont);
    let mut labelsize = 18.0f32;
    getparfloat("labelsize", &mut labelsize);
    let mut title = String::new();
    getparstring("title", &mut title);
    let mut title2 = String::new();
    let notitle2 = getparstring("title2", &mut title2) == 0;
    let mut titlefont = "Helvetica-Bold".to_string();
    getparstring("titlefont", &mut titlefont);
    let mut titlesize = 24.0f32;
    getparfloat("titlesize", &mut titlesize);
    let mut styles = "seismic".to_string();
    getparstring("style", &mut styles);
    let style = if styles == "normal" { 0 } else { 1 };
    checkpars();
    let i1beg = nint(((x1beg - f1) / d1) as f64).min(n1dsp - 1).max(0);
    x1beg = f1 + i1beg as f32 * d1;
    let i1end = nint(((x1end - f1) / d1) as f64).min(n1dsp - 1).max(0);
    x1end = f1 + i1end as f32 * d1;
    let i2beg = nint(((x2beg - f2) / d2) as f64).min(n2dsp - 1).max(0);
    x2beg = f2 + i2beg as f32 * d2;
    let i2end = nint(((x2end - f2) / d2) as f64).min(n2dsp - 1).max(0);
    x2end = f2 + i2end as f32 * d2;
    let n1c = 1 + (i1end - i1beg).abs();
    let n2c = 1 + (i2end - i2beg).abs();
    let zscale = if wclip != bclip {
        (255.0f64 / (wclip - bclip) as f64) as f32
    } else {
        1.0e10f32
    };
    let zoffset = -bclip * zscale;
    let i1step = if i1end > i1beg { 1 } else { -1 };
    let i2step = if i2end > i2beg { 1 } else { -1 };
    let n1s = nint((1.0f32 + (n1c - 1) as f32 * d1 / d1s) as f64).max(1);
    d1s = if n1s > 1 {
        d1 * (n1c - 1) as f32 / (n1s - 1) as f32
    } else {
        d1
    };
    let n2s = nint((1.0f32 + (n2c - 1) as f32 * d2 / d2s) as f64).max(1);
    d2s = if n2s > 1 {
        d2 * (n2c - 1) as f32 / (n2s - 1) as f32
    } else {
        d2
    };
    xbox *= 72.0;
    ybox *= 72.0;
    wbox *= 72.0;
    hbox *= 72.0;
    let mut bbox = [0; 4];
    psAxesBBox3(
        xbox, ybox, wbox, hbox, &labelfont, labelsize, &titlefont, titlesize, style, &mut bbox,
    );
    boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
    beginps();
    for i3 in 0..n3dsp {
        let mut cz = Vec::with_capacity((n1c * n2c) as usize);
        for i1c in 0..n1c {
            let i1 = i1beg + i1c * i1step;
            for i2c in 0..n2c {
                let i2 = i2beg + i2c * i2step;
                let mut zi = zoffset + z[(i1 + i2 * n1dsp + i3 * n1dsp * n2dsp) as usize] * zscale;
                if zi < 0.0 {
                    zi = 0.0;
                }
                if zi > 255.0 {
                    zi = 255.0;
                }
                cz.push(zi as u8);
            }
        }
        let sz = if n1s != n1c || n2s != n2c {
            let mut sz = vec![0u8; (n1s * n2s) as usize];
            intl2b(
                n2c, d2, 0.0, n1c, d1, 0.0, &cz, n2s, d2s, 0.0, n1s, d1s, 0.0, &mut sz,
            );
            sz
        } else {
            cz
        };
        let p1beg = if x1end > x1beg {
            -(d1s.abs() as f64) / 2.0
        } else {
            (d1s.abs() as f64) / 2.0
        } as f32;
        let p1end = -p1beg;
        let p2beg = if x2end > x2beg {
            -(d2s.abs() as f64) / 2.0
        } else {
            (d2s.abs() as f64) / 2.0
        } as f32;
        let p2end = -p2beg;
        newpage("1", i3 + 1);
        gsave();
        translate(xbox, ybox);
        let matrix = if style == 0 {
            [0.0, n1s as f32, n2s as f32, 0.0, 0.0, 0.0]
        } else {
            [n2s as f32, 0.0, 0.0, -(n1s as f32), 0.0, n1s as f32]
        };
        scale(wbox, hbox);
        image(n2s, n1s, 8, &matrix, &sz);
        grestore();
        let mut title2s = title2.clone();
        if !notitle2 {
            let f3s = f3 + i3 as f32 * d3;
            title2s.push_str(&sprintf!("= %.4g %s", f3s, ""));
        }
        psAxesBox3(
            xbox, ybox, wbox, hbox, x1beg, x1end, p1beg, p1end, d1num, f1num, n1tic, grid1,
            &label1, x2beg, x2end, p2beg, p2end, d2num, f2num, n2tic, grid2, &label2, &labelfont,
            labelsize, &title, &titlefont, titlesize, style, &title2s,
        );
        showpage();
    }
    endps();
    exit(CWP_Exit());
}
