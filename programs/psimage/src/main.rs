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

//! PSIMAGE - PostScript image plot of a uniformly sampled function.
//! Port of psplot/main/psimage.c (Revision 1.32).
//! AUTHOR: Dave Hale, Colorado School of Mines, 05/29/90.
//! MODIFIED: Craig Artley, Colorado School of Mines, 08/30/91: bounding box;
//! 12/16/93: color options (courtesy Dave Hale, Advance Geophysical).
//! Morten Wendell Pedersen, Aarhus University, 23/3-97: line widths.
//! Torsten Schoenfelder, Koeln, Germany, 006/07/97: legend (as in ximage,
//! by Berend Scheffers, Delft); 02/10/99: three-color interpolation.
//! Brian K. Macy, Phillips Petroleum, 01/14/99: curve plotting.
//! Ekkehart Tessmer, University of Hamburg, Germany, 08/22/2007: curve dashes.

#![allow(unused_assignments)]
use oxwave_byte_interp::*;
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{max, min, nint};
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::*;
use oxwave_f_psdrawcurve::*;
use oxwave_f_pslegendbox::*;
use oxwave_sort::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " PSIMAGE - PostScript IMAGE plot of a uniformly-sampled function f(x1,x2)",
    "            with the option to display a second attribute\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " psimage n1= [optional parameters] <binaryfile >postscriptfile\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1\t\t\t number of samples in 1st (fast) dimension\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " d1=1.0\t\t sampling interval in 1st dimension\t\t",
    " f1=0.0\t\t first sample in 1st dimension\t\t\t",
    " n2=all\t\t number of samples in 2nd (slow) dimension\t",
    " d2=1.0\t\t sampling interval in 2nd dimension\t\t",
    " f2=0.0\t\t first sample in 2nd dimension\t\t\t",
    " perc=100.0\t\t percentile used to determine clip\t\t",
    " clip=(perc percentile) clip used to determine bclip and wclip\t\t",
    " bperc=perc\t\t percentile for determining black clip value\t",
    " wperc=100.0-perc\t percentile for determining white clip value\t",
    " bclip=clip\t\t data values outside of [bclip,wclip] are clipped",
    " wclip=-clip\t\t data values outside of [bclip,wclip] are clipped",
    "                        bclip and wclip will be set to be inside       ",
    "                        [lbeg,lend] if lbeg and/or lend are supplied   ",
    " threecolor=1\t\t supply 3 color values instead of only two,\t",
    "                        using not only black and white, but f.e. red,\t",
    "                        green and blue\t\t\t\t\t",
    " brgb=0.0,0.0,0.0\t red, green, blue values corresponding to black\t",
    " grgb=1.0,1.0,1.0\t red, green, blue values corresponding to grey\t",
    " wrgb=1.0,1.0,1.0\t red, green, blue values corresponding to white\t",
    " bhls=0.0,0.0,0.0\t hue, lightness, saturation corresponding to black",
    " ghls=0.0,1.0,0.0\t hue, lightness, saturation corresponding to grey",
    " whls=0.0,1.0,0.0\t hue, lightness, saturation corresponding to white",
    " bps=12\t\t bits per sample for color plots, either 12 or 24",
    " d1s=1.0\t\t factor by which to scale d1 before imaging\t",
    " d2s=1.0\t\t factor by which to scale d2 before imaging\t",
    " verbose=1\t\t =1 for info printed on stderr (0 for no info)\t",
    " xbox=1.5\t\t offset in inches of left side of axes box\t",
    " ybox=1.5\t\t offset in inches of bottom side of axes box\t",
    " width=6.0\t\t width in inches of axes box\t\t\t",
    " height=8.0\t\t height in inches of axes box\t\t\t",
    " x1beg=x1min\t\t value at which axis 1 begins\t\t\t",
    " x1end=x1max\t\t value at which axis 1 ends\t\t\t",
    " d1num=0.0\t\t numbered tic interval on axis 1 (0.0 for automatic)",
    " f1num=x1min\t\t first numbered tic on axis 1 (used if d1num not 0.0)",
    " n1tic=1\t\t number of tics per numbered tic on axis 1\t",
    " grid1=none\t\t grid lines on axis 1 - none, dot, dash, or solid",
    " label1=\t\t label on axis 1\t\t\t\t",
    " x2beg=x2min\t\t value at which axis 2 begins\t\t\t",
    " x2end=x2max\t\t value at which axis 2 ends\t\t\t",
    " d2num=0.0\t\t numbered tic interval on axis 2 (0.0 for automatic)",
    " f2num=x2min\t\t first numbered tic on axis 2 (used if d2num not 0.0)",
    " n2tic=1\t\t number of tics per numbered tic on axis 2\t",
    " grid2=none\t\t grid lines on axis 2 - none, dot, dash, or solid",
    " label2=\t\t label on axis 2\t\t\t\t",
    " labelfont=Helvetica\t font name for axes labels\t\t\t",
    " labelsize=18\t\t font size for axes labels\t\t\t",
    " title=\t\t title of plot\t\t\t\t\t",
    " titlefont=Helvetica-Bold font name for title\t\t\t\t",
    " titlesize=24\t\t  font size for title\t\t\t\t",
    " titlecolor=black\t color of title\t\t\t\t\t",
    " axescolor=black\t color of axes\t\t\t\t\t",
    " gridcolor=black\t color of grid\t\t\t\t\t",
    " axeswidth=1            width (in points) of axes                      ",
    " ticwidth=axeswidth     width (in points) of tic marks\t\t\t",
    " gridwidth=axeswidth    width (in points) of grid lines\t\t",
    " style=seismic\t\t normal (axis 1 horizontal, axis 2 vertical) or\t",
    "\t\t\t seismic (axis 1 vertical, axis 2 horizontal)\t",
    " legend=0\t         =1 display the color scale\t\t\t",
    " lnice=0                =1 nice legend arrangement                     ",
    "                        (overrides ybox,lx,width,height parameters)    ",
    " lstyle=vertleft \tVertical, axis label on left side   \t\t",
    "\t\t\t =vertright (Vertical, axis label on right side)",
    "\t\t\t =horibottom (Horizontal, axis label on bottom)\t",
    " units=\t\t unit label for legend\t\t\t\t",
    " legendfont=times_roman10    font name for title\t\t\t",
    " following are defaults for lstyle=0. They are changed for other lstyles",
    " lwidth=1.2\t\t colorscale (legend) width in inches \t\t",
    " lheight=height/3     \t colorscale (legend) height in inches\t\t",
    " lx=1.0\t\t colorscale (legend) x-position in inches\t",
    " ly=(height-lheight)/2+xybox    colorscale (legend) y-position in pixels",
    " lbeg= lmin or wclip-5*perc    value at which legend axis begins\t",
    " lend= lmax or bclip+5*perc    value at which legend axis ends        \t",
    " ldnum=0.0\t numbered tic interval on legend axis (0.0 for automatic)",
    " lfnum=lmin\t first numbered tic on legend axis (used if d1num not 0.0)",
    " lntic=1\t number of tics per numbered tic on legend axis ",
    " lgrid=none\t grid lines on legend axis - none, dot, dash, or solid",
    "\t\t\t\t\t\t\t\t\t",
    " curve=curve1,curve2,...  file(s) containing points to draw curve(s)   ",
    " npair=n1,n2,n2,...            number(s) of pairs in each file         ",
    " curvecolor=black,..\t color of curve(s)\t\t\t\t",
    " curvewidth=axeswidth\t width (in points) of curve(s)\t\t\t",
    " curvedash=0            solid curve(s), dash indices 1,...,11 produce  ",
    "                        curve(s) with various dash styles              ",
    "\t\t\t\t\t\t\t\t\t",
    " infile=none            filename of second attribute n1xn2 array       ",
    "                        values must be from range 0.0 - 1.0            ",
    "                        (plain unformatted C-style file)               ",
    " bckgr=0.5              background gray value\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " NOTES:\t\t\t\t\t\t\t\t",
    " The curve file is an ascii file with the points specified as x1 x2 \t",
    " pairs, one pair to a line.  A \"vector\" of curve files and curve\t",
    " colors may be specified as curvefile=file1,file2,etc.\t\t\t",
    " and curvecolor=color1,color2,etc, and the number of pairs of values   ",
    " in each file as npair=npair1,npair2,... .\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " You may eliminate the blocky appearance of psimages by adjusting the  ",
    " d1s= and d2s= parameters, so that psimages appear similar to ximages.\t",
    "\t\t\t\t\t\t\t\t\t",
    " All color specifications may also be made in X Window style Hex format",
    " example:   axescolor=#255\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Some example colormap settings:\t\t\t\t\t",
    " red white blue: wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 \t\t",
    " white red blue: wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0 \t",
    " blue red white: wrgb=0.0,0.0,1.0 grgb=1.0,0.0,0.0 brgb=1.0,1.0,1.0 \t",
    " red green blue: wrgb=1.0,0,0 grgb=0,1.0,0 brgb=0,0,1.0\t\t",
    " orange light-blue green: wrgb=1.0,.5,0 grgb=0,.7,1.0 brgb=0,1.0,0\t",
    " red light-blue dark blue: wrgb=0.0,0,1.0 grgb=0,1.0,1.0 brgb=0,0,1.0 \t",
    " \t\t\t\t\t\t\t\t\t",
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

macro_rules! floatpar {
    ($name:literal, $default:expr) => {{
        let mut value: f32 = $default;
        getparfloat($name, &mut value);
        value
    }};
}
macro_rules! intpar {
    ($name:literal, $default:expr) => {{
        let mut value: i32 = $default;
        getparint($name, &mut value);
        value
    }};
}
macro_rules! stringpar {
    ($name:literal, $default:expr) => {{
        let mut value = String::from($default);
        getparstring($name, &mut value);
        value
    }};
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let infp = stdin();

    // Parameters describing the first dimension.
    let mut n1 = 0i32;
    if getparint("n1", &mut n1) == 0 {
        err!("must specify n1!\n");
    }
    let d1 = floatpar!("d1", 1.0);
    let f1 = floatpar!("f1", 0.0);
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

    // Parameters describing the second dimension.
    let mut n2 = 0i32;
    if getparint("n2", &mut n2) == 0 {
        if efseeko(&infp, 0, SEEK_END) != 0 {
            err!("must specify n2 if in a pipe!\n");
        }
        let nfloats = (eftello(&infp) / 4) as i32;
        efseeko(&infp, 0, SEEK_SET);
        // arm64 C integer division by zero yields zero.
        n2 = nfloats.checked_div(n1).unwrap_or(0);
    }
    let d2 = floatpar!("d2", 1.0);
    let f2 = floatpar!("f2", 0.0);
    let x2min = if d2 > 0.0 {
        f2
    } else {
        f2 + (n2 - 1) as f32 * d2
    };
    let x2max = if d2 < 0.0 {
        f2
    } else {
        f2 + (n2 - 1) as f32 * d2
    };

    // Color[][0] is black, color[][2] is white.
    let threecolor = intpar!("threecolor", 1);
    let mut bps = 8;
    let mut hls = false;
    let mut colors = [[0.0f32, 0.5, 1.0]; 3];
    if countparval("brgb") != 0 || countparval("wrgb") != 0 {
        let mut brgb = [0.0f32; 3];
        let mut wrgb = [1.0f32; 3];
        getparfloat("brgb", &mut brgb[..]);
        getparfloat("wrgb", &mut wrgb[..]);
        let mut grgb = [0.0f32; 3];
        for i in 0..3 {
            grgb[i] = ((brgb[i] + wrgb[i]) as f64 / 2.0) as f32;
        }
        if threecolor == 1 {
            getparfloat("grgb", &mut grgb[..]);
        }
        for i in 0..3 {
            brgb[i] = max(0.0, min(1.0, brgb[i]));
            grgb[i] = max(0.0, min(1.0, grgb[i]));
            wrgb[i] = max(0.0, min(1.0, wrgb[i]));
            colors[i] = [brgb[i], grgb[i], wrgb[i]];
        }
        bps = intpar!("bps", 12);
        if bps != 12 && bps != 24 {
            err!("bps must equal 12 or 24 for color plots!\n");
        }
    } else if countparval("bhls") != 0 || countparval("whls") != 0 {
        hls = true;
        let mut bhls = [0.0f32, 0.0, 0.0];
        let mut whls = [0.0f32, 1.0, 0.0];
        getparfloat("bhls", &mut bhls[..]);
        getparfloat("whls", &mut whls[..]);
        let mut ghls = [0.0f32; 3];
        for i in 0..3 {
            ghls[i] = ((bhls[i] + whls[i]) as f64 / 2.0) as f32;
        }
        if threecolor == 1 {
            getparfloat("ghls", &mut ghls[..]);
        }
        for i in 1..3 {
            bhls[i] = max(0.0, min(1.0, bhls[i]));
            ghls[i] = max(0.0, min(1.0, ghls[i]));
            whls[i] = max(0.0, min(1.0, whls[i]));
        }
        for i in 0..3 {
            colors[i] = [bhls[i], ghls[i], whls[i]];
        }
        bps = intpar!("bps", 12);
        if bps != 12 && bps != 24 {
            err!("bps must equal 12 or 24 for color plots!\n");
        }
    }

    let legend = intpar!("legend", 0);
    let units = stringpar!("units", "");
    let _legendfont = stringpar!("legendfont", "times_roman10");

    // Set up curve plotting.
    let curve = countparval("curve") as usize;
    let mut npair = vec![0i32; curve];
    let mut curvefile = Vec::<String>::new();
    let mut curvecolor = Vec::<String>::new();
    let mut x1curve = Vec::<Vec<f32>>::new();
    let mut x2curve = Vec::<Vec<f32>>::new();
    if curve != 0 {
        getparstringarray("curve", &mut curvefile);
        getparint("npair", &mut npair[..]);
        let mut ncurvecolor = countparval("curvecolor") as usize;
        if ncurvecolor < curve {
            if getparstringarray("curvecolor", &mut curvecolor) == 0 {
                curvecolor = vec!["black".into()];
                ncurvecolor = 1;
            }
            curvecolor.resize(curve, curvecolor[ncurvecolor - 1].clone());
        } else if ncurvecolor != 0 {
            getparstringarray("curvecolor", &mut curvecolor);
        }
        for j in 0..curve {
            let curvefp = efopen(&curvefile[j], "r");
            let mut x1 = vec![0f32; npair[j] as usize];
            let mut x2 = vec![0f32; npair[j] as usize];
            for i in 0..npair[j] as usize {
                if let Some(x) = fscanf_f32(&curvefp) {
                    x1[i] = x;
                }
                if let Some(x) = fscanf_f32(&curvefp) {
                    x2[i] = x;
                }
            }
            efclose(&curvefp);
            x1curve.push(x1);
            x2curve.push(x2);
        }
    }

    let mut infile = String::new();
    let has_infile = getparstring("infile", &mut infile) != 0;
    let in2fp = if has_infile {
        match fopen(&infile, "r") {
            Some(fp) => Some(fp),
            None => {
                fprintf!(stderr(), "error opening infile: %s\n", infile);
                err!("stopped.\n");
            }
        }
    } else {
        None
    };

    // Read the binary data and optional second attribute.
    let nz = n1.wrapping_mul(n2);
    let mut z = Vec::<f32>::new();
    if z.try_reserve_exact(nz as usize).is_err() {
        // Vec rejects the same impossible allocation as ealloc1float.
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    z.resize(nz as usize, 0.0);
    if fread(&mut z, &infp) != nz as usize {
        err!("error reading input file!\n");
    }
    let mut y = Vec::<f32>::new();
    let mut bckgr = 0.5f32;
    if let Some(in2fp) = &in2fp {
        bckgr = floatpar!("bckgr", 0.5);
        y = vec![0f32; nz as usize];
        if fread(&mut y, in2fp) != nz as usize {
            err!("error reading modulating input file!\n");
        }
        let mut infile_min = 1.0f32;
        let mut infile_max = 0.0f32;
        for iz in 0..nz as usize {
            infile_min = min(infile_min, y[iz]);
            infile_max = max(infile_max, y[iz]);
        }
        if infile_min < 0.0 || infile_max > 1.0 {
            err!("error: value of infile out of range 0..1!\n");
        }
    }
    // For a zero dimension C reads past malloc(0) while selecting clips and
    // rasterizing. The reference build consistently observes zero-filled
    // allocator storage. Provide those guard samples explicitly.
    if nz == 0 {
        z.resize(max(1, n1) as usize, 0.0);
        if has_infile {
            y.resize(max(1, n1) as usize, 0.0);
        }
    }

    // Determine clips from percentiles, keeping getpar's short circuit order.
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
        let perc = floatpar!("perc", 100.0);
        let mut temp = z.clone();
        if getparfloat("bclip", &mut bclip) == 0 {
            let bperc = floatpar!("bperc", perc);
            let mut iz = ((nz as f32 * bperc) as f64 / 100.0) as i32;
            if iz < 0 {
                iz = 0;
            }
            if iz > nz - 1 {
                iz = nz - 1;
            }
            if nz == 0 {
                bclip = 0.0;
            } else {
                qkfind(iz as usize, &mut temp);
                bclip = temp[iz as usize];
            }
        }
        if getparfloat("wclip", &mut wclip) == 0 {
            let wperc = floatpar!("wperc", (100.0 - perc as f64) as f32);
            let mut iz = ((nz as f32 * wperc) as f64 / 100.0) as i32;
            if iz < 0 {
                iz = 0;
            }
            if iz > nz - 1 {
                iz = nz - 1;
            }
            if nz == 0 {
                wclip = 0.0;
            } else {
                qkfind(iz as usize, &mut temp);
                wclip = temp[iz as usize];
            }
        }
    }
    let verbose = intpar!("verbose", 1);
    if verbose != 0 {
        warn!("bclip=%g wclip=%g", bclip, wclip);
    }

    let mut d1s = floatpar!("d1s", 1.0);
    let mut d2s = floatpar!("d2s", 1.0);
    d1s = (d1s as f64).abs() as f32;
    d1s *= d1;
    d2s = (d2s as f64).abs() as f32;
    d2s *= d2;

    // Axes parameters.
    let mut xbox = floatpar!("xbox", 1.5);
    let mut ybox = floatpar!("ybox", 1.5);
    let mut width = floatpar!("wbox", 6.0);
    getparfloat("width", &mut width);
    let mut height = floatpar!("hbox", 8.0);
    getparfloat("height", &mut height);
    let lnice = intpar!("lnice", 0);
    if lnice == 1 {
        ybox = 2.2;
        width = 5.4;
        height = 7.2;
    }
    let mut x1beg = floatpar!("x1beg", x1min);
    let mut x1end = floatpar!("x1end", x1max);
    let d1num = floatpar!("d1num", 0.0);
    let f1num = floatpar!("f1num", x1min);
    let n1tic = intpar!("n1tic", 1);
    let grid1s = stringpar!("grid1", "none");
    let grid1 = match grid1s.as_str() {
        "dot" => 1,
        "dash" => 2,
        "solid" => 3,
        _ => 0,
    };
    let label1 = stringpar!("label1", "");
    let mut x2beg = floatpar!("x2beg", x2min);
    let mut x2end = floatpar!("x2end", x2max);
    let d2num = floatpar!("d2num", 0.0);
    let f2num = floatpar!("f2num", 0.0);
    let n2tic = intpar!("n2tic", 1);
    let grid2s = stringpar!("grid2", "none");
    let grid2 = match grid2s.as_str() {
        "dot" => 1,
        "dash" => 2,
        "solid" => 3,
        _ => 0,
    };
    let label2 = stringpar!("label2", "");
    let labelfont = stringpar!("labelfont", "Helvetica");
    let labelsize = floatpar!("labelsize", 18.0);
    let title = stringpar!("title", "");
    let titlefont = stringpar!("titlefont", "Helvetica-Bold");
    let titlesize = floatpar!("titlesize", 24.0);
    let titlecolor = stringpar!("titlecolor", "black");
    let axescolor = stringpar!("axescolor", "black");
    let gridcolor = stringpar!("gridcolor", "black");
    let axeswidth = floatpar!("axeswidth", 1.0);
    let ticwidth = floatpar!("ticwidth", axeswidth);
    let gridwidth = floatpar!("gridwidth", axeswidth);

    let mut curvewidth = vec![0f32; curve];
    let mut curvedash = vec![0i32; curve];
    if curve != 0 {
        let mut ncurvewidth = countparval("curvewidth") as usize;
        if ncurvewidth < curve {
            if getparfloat("curvewidth", &mut curvewidth[..]) == 0 {
                curvewidth[0] = axeswidth;
                ncurvewidth = 1;
            }
            for i in ncurvewidth..curve {
                curvewidth[i] = curvewidth[ncurvewidth - 1];
            }
        } else {
            curvewidth.resize(ncurvewidth, 0.0);
            getparfloat("curvewidth", &mut curvewidth[..]);
        }
        let mut ncurvedash = countparval("curvedash") as usize;
        if ncurvedash < curve {
            if getparint("curvedash", &mut curvedash[..]) == 0 {
                curvedash[0] = 0;
                ncurvedash = 1;
            }
            for i in ncurvedash..curve {
                curvedash[i] = curvedash[ncurvedash - 1];
            }
        } else {
            curvedash.resize(ncurvedash, 0);
            getparint("curvedash", &mut curvedash[..]);
        }
    }
    let styles = stringpar!("style", "seismic");
    let style = if styles == "normal" { 0 } else { 1 };

    let mut lmin = f32::MAX;
    let mut lmax = -f32::MAX;
    if legend != 0 {
        for lz in 0..nz as usize {
            lmin = min(lmin, z[lz]);
            lmax = max(lmax, z[lz]);
        }
        if verbose == 2 {
            warn!("lmin=%g lmax=%g", lmin, lmax);
        }
    }
    let mut lbeg = 0f32;
    let mut lend = 0f32;
    let mut lbegsup = false;
    let mut lendsup = false;
    let mut lstyle = 0;
    let mut lwidth = 0f32;
    let mut lheight = 0f32;
    let mut lx = 0f32;
    let mut ly = 0f32;
    let mut ugrid = 3;
    if legend != 0 {
        lbeg = lmin;
        lbegsup = getparfloat("lbeg", &mut lbeg) != 0;
        lend = lmax;
        lendsup = getparfloat("lend", &mut lend) != 0;
        wclip = max(lbeg, wclip);
        bclip = min(lend, bclip);
        if !lbegsup {
            let rangeperc = ((bclip - wclip) as f64 / 20.0) as f32;
            lbeg = wclip - rangeperc;
        }
        if !lendsup {
            let rangeperc = ((bclip - wclip) as f64 / 20.0) as f32;
            lend = bclip + rangeperc;
        }
        let _lfnum = floatpar!("lfnum", lmin);
        let lstyles = stringpar!("lstyle", "vertleft");
        if lstyles == "vertright" {
            lstyle = 1;
        } else if lstyles == "horibottom" {
            lstyle = 2;
        }
        lwidth = 0.1;
        lheight = height / 2.0f32;
        if lstyle == 2 {
            lwidth = (width as f64 / 1.2) as f32;
            lheight = 0.24;
        }
        getparfloat("lwidth", &mut lwidth);
        getparfloat("lheight", &mut lheight);
        lx = 0.8;
        ly = ybox + (height - lheight) / 2.0f32;
        if lstyle == 1 {
            lx = ((xbox + width) as f64 + 0.1) as f32;
        } else if lstyle == 2 {
            lx = (xbox as f64 + (width - lwidth) as f64 / 2.0) as f32;
            ly = 1.0;
        }
        getparfloat("lx", &mut lx);
        if lnice == 1 {
            lx = 8.0;
        }
        getparfloat("ly", &mut ly);
        let lgrids = stringpar!("lgrid", "none");
        ugrid = match lgrids.as_str() {
            "dot" => 1,
            "dash" => 2,
            "solid" => 3,
            _ => 0,
        };
    }

    // Snap each extent to the sampled coordinates.
    let i1beg = max(0, min(n1 - 1, nint(((x1beg - f1) / d1) as f64)));
    x1beg = f1 + i1beg as f32 * d1;
    let i1end = max(0, min(n1 - 1, nint(((x1end - f1) / d1) as f64)));
    x1end = f1 + i1end as f32 * d1;
    let i2beg = max(0, min(n2 - 1, nint(((x2beg - f2) / d2) as f64)));
    x2beg = f2 + i2beg as f32 * d2;
    let i2end = max(0, min(n2 - 1, nint(((x2end - f2) / d2) as f64)));
    x2end = f2 + i2end as f32 * d2;

    let mut data_legend = Vec::<u8>::new();
    let mut ln = 256i32;
    let mut lf = 0f32;
    let mut ldnum = 0f32;
    if legend != 0 {
        if !lbegsup {
            ln += 13;
        }
        if !lendsup {
            ln += 13;
        }
        if !lbegsup {
            for _ in 0..13 {
                data_legend.push(255);
            }
        }
        for lcount in (0..=255).rev() {
            data_legend.push(lcount as u8);
        }
        data_legend.resize(ln as usize, 0);
        lf = lbeg;
        let ld = (lend - lbeg) / (ln - 1) as f32;
        ldnum = floatpar!("ldnum", 0.0);
        checkpars();
        let ilbeg = max(0, min(ln - 1, nint(((lbeg - lf) / ld) as f64)));
        lbeg = lf + ilbeg as f32 * ld;
        let ilend = max(0, min(ln - 1, nint(((lend - lf) / ld) as f64)));
        lend = lf + ilend as f32 * ld;
    }

    let n1c = 1 + (i1end - i1beg).abs();
    let n2c = 1 + (i2end - i2beg).abs();
    let mut cz = vec![0u8; (n1c * n2c) as usize];
    let mut cy = if has_infile {
        vec![0u8; (n1c * n2c) as usize]
    } else {
        Vec::new()
    };
    let zscale = if wclip != bclip {
        (255.0 / (wclip - bclip) as f64) as f32
    } else {
        1.0e10f32
    };
    let zoffset = -bclip * zscale;
    let i1step = if i1end > i1beg { 1 } else { -1 };
    let i2step = if i2end > i2beg { 1 } else { -1 };
    let mut i1 = i1beg;
    let mut k = 0usize;
    for _i1c in 0..n1c {
        let mut i2 = i2beg;
        for _i2c in 0..n2c {
            let mut zi = zoffset + z[(i1 + i2 * n1) as usize] * zscale;
            if zi < 0.0 {
                zi = 0.0;
            }
            if zi > 255.0 {
                zi = 255.0;
            }
            cz[k] = zi as u8;
            if has_infile {
                cy[k] = ((y[(i1 + i2 * n1) as usize] as f64 * 255.0) as i32) as u8;
            }
            k += 1;
            i2 += i2step;
        }
        i1 += i1step;
    }
    drop(z);
    drop(y);

    let n1s = max(1, nint((1.0f32 + (n1c - 1) as f32 * d1 / d1s) as f64));
    d1s = if n1s > 1 {
        d1 * (n1c - 1) as f32 / (n1s - 1) as f32
    } else {
        d1
    };
    let n2s = max(1, nint((1.0f32 + (n2c - 1) as f32 * d2 / d2s) as f64));
    d2s = if n2s > 1 {
        d2 * (n2c - 1) as f32 / (n2s - 1) as f32
    } else {
        d2
    };
    let sz = if n1s != n1c || n2s != n2c {
        let mut sz = vec![0u8; (n1s * n2s) as usize];
        intl2b(
            n2c, d2, 0.0, n1c, d1, 0.0, &cz, n2s, d2s, 0.0, n1s, d1s, 0.0, &mut sz,
        );
        sz
    } else {
        cz
    };
    let sy = if has_infile && (n1s != n1c || n2s != n2c) {
        let mut sy = vec![0u8; (n1s * n2s) as usize];
        intl2b(
            n2c, d2, 0.0, n1c, d1, 0.0, &cy, n2s, d2s, 0.0, n1s, d1s, 0.0, &mut sy,
        );
        sy
    } else {
        cy
    };
    let p1beg = (if x1end > x1beg {
        -(d1s as f64).abs() / 2.0
    } else {
        (d1s as f64).abs() / 2.0
    }) as f32;
    let p1end = (if x1end > x1beg {
        (d1s as f64).abs() / 2.0
    } else {
        -(d1s as f64).abs() / 2.0
    }) as f32;
    let p2beg = (if x2end > x2beg {
        -(d2s as f64).abs() / 2.0
    } else {
        (d2s as f64).abs() / 2.0
    }) as f32;
    let p2end = (if x2end > x2beg {
        (d2s as f64).abs() / 2.0
    } else {
        -(d2s as f64).abs() / 2.0
    }) as f32;

    xbox = (xbox as f64 * 72.0) as f32;
    ybox = (ybox as f64 * 72.0) as f32;
    width = (width as f64 * 72.0) as f32;
    height = (height as f64 * 72.0) as f32;
    if legend != 0 {
        lx = (lx as f64 * 72.0) as f32;
        ly = (ly as f64 * 72.0) as f32;
        lwidth = (lwidth as f64 * 72.0) as f32;
        lheight = (lheight as f64 * 72.0) as f32;
    }
    let mut bbox = [0i32; 4];
    psAxesBBox(
        xbox, ybox, width, height, &labelfont, labelsize, &titlefont, titlesize, style, &mut bbox,
    );
    if legend != 0 {
        let mut lbbox = [0i32; 4];
        psLegendBBox(
            lx, ly, lwidth, lheight, &labelfont, labelsize, lstyle, &mut lbbox,
        );
        bbox[0] = min(bbox[0], lbbox[0]);
        bbox[1] = min(bbox[1], lbbox[1]);
        bbox[2] = max(bbox[2], lbbox[2]);
        bbox[3] = max(bbox[3], lbbox[3]);
    }
    boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
    begineps();
    gsave();
    translate(xbox, ybox);
    let matrix = if style == 0 {
        [0.0, n1s as f32, n2s as f32, 0.0, 0.0, 0.0]
    } else {
        [n2s as f32, 0.0, 0.0, -(n1s as f32), 0.0, n1s as f32]
    };
    scale(width, height);
    if !has_infile {
        drawimage(hls, &colors, n2s, n1s, bps, &matrix, &sz);
    } else {
        drawimage2(hls, &colors, n2s, n1s, bps, &matrix, &sz, &sy, bckgr);
    }
    grestore();

    if legend != 0 {
        gsave();
        translate(lx, ly);
        scale(lwidth, lheight);
        if lstyle == 0 || lstyle == 1 {
            let labmatrix = [1.0, 0.0, 0.0, ln as f32, 0.0, 0.0];
            drawimage(hls, &colors, 1, ln, bps, &labmatrix, &data_legend);
        } else {
            let labmatrix = [-1.0, 0.0, 0.0, ln as f32, 0.0, 0.0];
            rotate(-90.0);
            drawimage(hls, &colors, 1, ln, bps, &labmatrix, &data_legend);
            rotate(90.0);
        }
        grestore();
    }
    for i in 0..curve {
        gsave();
        psDrawCurve(
            xbox,
            ybox,
            width,
            height,
            x1beg,
            x1end,
            p1beg,
            p1end,
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
    gsave();
    psAxesBox(
        xbox,
        ybox,
        width,
        height,
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
    grestore();
    if legend != 0 {
        gsave();
        psLegendBox(
            lx, ly, lwidth, lheight, lbeg, lend, 0.0, 0.0, ldnum, lf, 1, ugrid, &units, &labelfont,
            labelsize, &axescolor, &gridcolor, lstyle,
        );
        grestore();
    }
    showpage();
    endeps();
    exit(CWP_Exit());
}

// Image helpers below are private routines from psimage.c.
fn grayscale(hls: bool, colors: &[[f32; 3]; 3]) -> bool {
    (hls && colors[1] == [0.0, 0.5, 1.0] && colors[2] == [0.0, 0.0, 0.0])
        || (colors[0] == [0.0, 0.5, 1.0]
            && colors[1] == [0.0, 0.5, 1.0]
            && colors[2] == [0.0, 0.5, 1.0])
}

fn drawimage(
    hls: bool,
    colors: &[[f32; 3]; 3],
    width: i32,
    height: i32,
    bps: i32,
    matrix: &[f32],
    z: &[u8],
) {
    if grayscale(hls, colors) {
        image(width, height, 8, matrix, z);
        return;
    }
    let mut rtab = [0u8; 256];
    let mut gtab = [0u8; 256];
    let mut btab = [0u8; 256];
    if hls {
        let hbase_low = colors[0][0];
        let hscale_low = ((colors[0][1] - colors[0][0]) as f64 / 127.0) as f32;
        let lbase_low = colors[1][0];
        let lscale_low = ((colors[1][1] - colors[1][0]) as f64 / 127.0) as f32;
        let sbase_low = colors[2][0];
        let sscale_low = ((colors[2][1] - colors[2][0]) as f64 / 127.0) as f32;
        let hbase_hi = colors[0][1];
        let hscale_hi = ((colors[0][2] - colors[0][1]) as f64 / 127.0) as f32;
        let lbase_hi = colors[1][1];
        let lscale_hi = ((colors[1][2] - colors[1][1]) as f64 / 127.0) as f32;
        let sbase_hi = colors[2][1];
        let sscale_hi = ((colors[2][2] - colors[2][1]) as f64 / 127.0) as f32;
        for zi in 0..128 {
            let h = hbase_low + zi as f32 * hscale_low;
            let l = lbase_low + zi as f32 * lscale_low;
            let s = sbase_low + zi as f32 * sscale_low;
            let (r, g, b) = hlsrgb((h as f64 * 360.0) as f32, l, s);
            rtab[zi] = (r as f64 * 255.0) as i32 as u8;
            gtab[zi] = (g as f64 * 255.0) as i32 as u8;
            btab[zi] = (b as f64 * 255.0) as i32 as u8;
        }
        for zi in 128..256 {
            let h = hbase_hi + (zi - 128) as f32 * hscale_hi;
            let l = lbase_hi + (zi - 128) as f32 * lscale_hi;
            let s = sbase_hi + (zi - 128) as f32 * sscale_hi;
            let (r, g, b) = hlsrgb((h as f64 * 360.0) as f32, l, s);
            rtab[zi] = (r as f64 * 255.0) as i32 as u8;
            gtab[zi] = (g as f64 * 255.0) as i32 as u8;
            btab[zi] = (b as f64 * 255.0) as i32 as u8;
        }
    } else {
        let rbase_low = (colors[0][0] as f64 * 255.0) as f32;
        let rscale_low = (2.0 * (colors[0][1] - colors[0][0]) as f64) as f32;
        let gbase_low = (colors[1][0] as f64 * 255.0) as f32;
        let gscale_low = (2.0 * (colors[1][1] - colors[1][0]) as f64) as f32;
        let bbase_low = (colors[2][0] as f64 * 255.0) as f32;
        let bscale_low = (2.0 * (colors[2][1] - colors[2][0]) as f64) as f32;
        let rbase_hi = (colors[0][1] as f64 * 255.0) as f32;
        let rscale_hi = (2.0 * (colors[0][2] - colors[0][1]) as f64) as f32;
        let gbase_hi = (colors[1][1] as f64 * 255.0) as f32;
        let gscale_hi = (2.0 * (colors[1][2] - colors[1][1]) as f64) as f32;
        let bbase_hi = (colors[2][1] as f64 * 255.0) as f32;
        let bscale_hi = (2.0 * (colors[2][2] - colors[2][1]) as f64) as f32;
        for zi in 0..128 {
            rtab[zi] = (rbase_low + zi as f32 * rscale_low) as i32 as u8;
            gtab[zi] = (gbase_low + zi as f32 * gscale_low) as i32 as u8;
            btab[zi] = (bbase_low + zi as f32 * bscale_low) as i32 as u8;
        }
        for zi in 128..256 {
            rtab[zi] = (rbase_hi + (zi - 128) as f32 * rscale_hi) as i32 as u8;
            gtab[zi] = (gbase_hi + (zi - 128) as f32 * gscale_hi) as i32 as u8;
            btab[zi] = (bbase_hi + (zi - 128) as f32 * bscale_hi) as i32 as u8;
        }
    }
    let mut rgb = vec![0u8; (width * height * 3) as usize];
    let mut zp = 0usize;
    let mut rgbp = 0usize;
    for _row in 0..height {
        for _col in 0..width {
            let zi = z[zp] as usize;
            zp += 1;
            rgb[rgbp] = rtab[zi];
            rgbp += 1;
            rgb[rgbp] = gtab[zi];
            rgbp += 1;
            rgb[rgbp] = btab[zi];
            rgbp += 1;
        }
    }
    if bps == 12 {
        compress12(width, height, &mut rgb);
    }
    rgbimage(width, height, bps / 3, matrix, &rgb);
}

#[allow(clippy::too_many_arguments)]
fn drawimage2(
    hls: bool,
    colors: &[[f32; 3]; 3],
    width: i32,
    height: i32,
    bps: i32,
    matrix: &[f32],
    z: &[u8],
    y: &[u8],
    bckgr: f32,
) {
    if grayscale(hls, colors) {
        image(width, height, 8, matrix, z);
        return;
    }
    let mut htab = [0f32; 256];
    let mut ltab = [0f32; 256];
    let mut stab = [0f32; 256];
    if hls {
        let hbase_low = colors[0][0];
        let hscale_low = ((colors[0][1] - colors[0][0]) as f64 / 127.0) as f32;
        let lbase_low = colors[1][0];
        let lscale_low = ((colors[1][1] - colors[1][0]) as f64 / 127.0) as f32;
        let sbase_low = colors[2][0];
        let sscale_low = ((colors[2][1] - colors[2][0]) as f64 / 127.0) as f32;
        let hbase_hi = colors[0][1];
        let hscale_hi = ((colors[0][2] - colors[0][1]) as f64 / 127.0) as f32;
        let lbase_hi = colors[1][1];
        let lscale_hi = ((colors[1][2] - colors[1][1]) as f64 / 127.0) as f32;
        let sbase_hi = colors[2][1];
        let sscale_hi = ((colors[2][2] - colors[2][1]) as f64 / 127.0) as f32;
        for zi in 0..128 {
            htab[zi] = hbase_low + zi as f32 * hscale_low;
            ltab[zi] = lbase_low + zi as f32 * lscale_low;
            stab[zi] = sbase_low + zi as f32 * sscale_low;
        }
        for zi in 128..256 {
            htab[zi] = hbase_hi + (zi - 128) as f32 * hscale_hi;
            ltab[zi] = lbase_hi + (zi - 128) as f32 * lscale_hi;
            stab[zi] = sbase_hi + (zi - 128) as f32 * sscale_hi;
        }
    } else {
        let rbase_low = (colors[0][0] as f64 * 255.0) as f32;
        let rscale_low = (2.0 * (colors[0][1] - colors[0][0]) as f64) as f32;
        let gbase_low = (colors[1][0] as f64 * 255.0) as f32;
        let gscale_low = (2.0 * (colors[1][1] - colors[1][0]) as f64) as f32;
        let bbase_low = (colors[2][0] as f64 * 255.0) as f32;
        let bscale_low = (2.0 * (colors[2][1] - colors[2][0]) as f64) as f32;
        let rbase_hi = (colors[0][1] as f64 * 255.0) as f32;
        let rscale_hi = (2.0 * (colors[0][2] - colors[0][1]) as f64) as f32;
        let gbase_hi = (colors[1][1] as f64 * 255.0) as f32;
        let gscale_hi = (2.0 * (colors[1][2] - colors[1][1]) as f64) as f32;
        let bbase_hi = (colors[2][1] as f64 * 255.0) as f32;
        let bscale_hi = (2.0 * (colors[2][2] - colors[2][1]) as f64) as f32;
        for zi in 0..128 {
            let rval = ((rbase_low + zi as f32 * rscale_low) as f64 / 255.0) as f32;
            let gval = ((gbase_low + zi as f32 * gscale_low) as f64 / 255.0) as f32;
            let bval = ((bbase_low + zi as f32 * bscale_low) as f64 / 255.0) as f32;
            (htab[zi], ltab[zi], stab[zi]) = rgbhls(rval, gval, bval);
            htab[zi] = (htab[zi] as f64 / 360.0) as f32;
        }
        for zi in 128..256 {
            let rval = ((rbase_hi + (zi - 128) as f32 * rscale_hi) as f64 / 255.0) as f32;
            let gval = ((gbase_hi + (zi - 128) as f32 * gscale_hi) as f64 / 255.0) as f32;
            let bval = ((bbase_hi + (zi - 128) as f32 * bscale_hi) as f64 / 255.0) as f32;
            (htab[zi], ltab[zi], stab[zi]) = rgbhls(rval, gval, bval);
            htab[zi] = (htab[zi] as f64 / 360.0) as f32;
        }
    }
    let mut rgb = vec![0u8; (width * height * 3) as usize];
    let mut zp = 0usize;
    let mut rgbp = 0usize;
    let mut yi = 0usize;
    for _row in 0..height {
        for _col in 0..width {
            let zi = z[zp] as usize;
            zp += 1;
            let hval = htab[zi];
            let lval = (bckgr as f64 + ((ltab[zi] - bckgr) * y[yi] as f32) as f64 / 255.0) as f32;
            let sval = ((stab[zi] * y[yi] as f32) as f64 / 255.0) as f32;
            yi += 1;
            let (rval, gval, bval) = hlsrgb((hval as f64 * 360.0) as f32, lval, sval);
            rgb[rgbp] = (rval as f64 * 255.0) as i32 as u8;
            rgbp += 1;
            rgb[rgbp] = (gval as f64 * 255.0) as i32 as u8;
            rgbp += 1;
            rgb[rgbp] = (bval as f64 * 255.0) as i32 as u8;
            rgbp += 1;
        }
    }
    if bps == 12 {
        compress12(width, height, &mut rgb);
    }
    rgbimage(width, height, bps / 3, matrix, &rgb);
}

// Identical in-place packing loops from drawimage and drawimage2.
fn compress12(width: i32, height: i32, rgb: &mut [u8]) {
    let rowbytes = 1 + (width * 12 - 1) / 8;
    let mut crgbp = 0usize;
    let mut rgbp = 0usize;
    for _row in 0..height {
        for _byte in 0..rowbytes - 1 {
            let hi = rgb[rgbp] >> 4;
            rgbp += 1;
            let lo = rgb[rgbp] >> 4;
            rgbp += 1;
            rgb[crgbp] = (hi << 4) | lo;
            crgbp += 1;
        }
        if width % 2 != 0 {
            let hi = rgb[rgbp] >> 4;
            rgbp += 1;
            rgb[crgbp] = hi << 4;
            crgbp += 1;
        } else {
            let hi = rgb[rgbp] >> 4;
            rgbp += 1;
            let lo = rgb[rgbp] >> 4;
            rgbp += 1;
            rgb[crgbp] = (hi << 4) | lo;
            crgbp += 1;
        }
    }
}

fn rgbvalue(n1: f32, n2: f32, mut hue: f32) -> f32 {
    while hue > 360.0 {
        hue = (hue as f64 - 360.0) as f32;
    }
    while hue < 0.0 {
        hue = (hue as f64 + 360.0) as f32;
    }
    if hue < 60.0 {
        (n1 as f64 + ((n2 - n1) * hue) as f64 / 60.0) as f32
    } else if hue < 180.0 {
        n2
    } else if hue < 240.0 {
        (n1 as f64 + (n2 - n1) as f64 * (240.0 - hue as f64) / 60.0) as f32
    } else {
        n1
    }
}

fn hlsrgb(h: f32, l: f32, s: f32) -> (f32, f32, f32) {
    let m2 = if l <= 0.5 {
        (l as f64 * (1.0 + s as f64)) as f32
    } else {
        l + s - l * s
    };
    let m1 = 2.0f32 * l - m2;
    if s == 0.0 {
        (l, l, l)
    } else {
        (
            rgbvalue(m1, m2, (h as f64 + 120.0) as f32),
            rgbvalue(m1, m2, h),
            rgbvalue(m1, m2, (h as f64 - 120.0) as f32),
        )
    }
}

fn rgbhls(r: f32, g: f32, b: f32) -> (f32, f32, f32) {
    let mut s = 0f32;
    let mut h = 0f32;
    let mut m1 = r;
    if g > m1 {
        m1 = g;
    }
    if b > m1 {
        m1 = b;
    }
    let mut m2 = r;
    if g < m2 {
        m2 = g;
    }
    if b < m2 {
        m2 = b;
    }
    let sss = m1 + m2;
    let ddd = m1 - m2;
    let l = (sss as f64 * 0.5) as f32;
    if ddd > 0.0 {
        if l <= 0.5 {
            s = ddd / sss;
        } else {
            s = (ddd as f64 / (2.0 - sss as f64)) as f32;
        }
        if r == m1 {
            h = (((g - b) / ddd) as f64 * 60.0) as f32;
        } else if g == m1 {
            h = (120.0 + ((b - r) / ddd) as f64 * 60.0) as f32;
        } else {
            h = (240.0 + ((r - g) / ddd) as f64 * 60.0) as f32;
            if h < 0.0 {
                h = (h as f64 + 360.0) as f32;
            }
        }
    }
    (h, l, s)
}
