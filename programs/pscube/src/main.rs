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
// Original source: Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.
//! PSCUBE, port of psplot/main/pscube.c, revision 1.15.
//! Authors: Craig Artley; original by Zhiming Li and Dave Hale.
//! Legend additions by Claudia Vanelle, University of Hamburg.
//! HLS-to-RGB routines adapted from Foley & Van Dam.
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSCUBE - PostScript image plot with Legend of a data CUBE       ",
    " \t\t\t\t\t\t\t\t\t",
    " pscube n1= n2= n3= [optional parameters] <binaryfile >postscriptfile\t",
    "    or\t\t\t\t\t\t\t\t\t",
    " pscube n1= n2= n3= front= side= top= [optional parameters] >postscriptfile",
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
    " perc=100.0             percentile used to determine clip\t\t",
    " clip=(perc percentile) clip used to determine bclip and wclip\t\t",
    " bperc=perc             percentile for determining black clip value\t",
    " wperc=100.0-perc       percentile for determining white clip value\t",
    " bclip=clip             data values outside of [bclip,wclip] are clipped",
    " wclip=-clip            data values outside of [bclip,wclip] are clipped",
    " brgb=0.0,0.0,0.0       red, green, blue values corresponding to black\t",
    " wrgb=1.0,1.0,1.0       red, green, blue values corresponding to white\t",
    " bhls=0.0,0.0,0.0       hue, lightness, saturation corresponding to black",
    " whls=0.0,1.0,0.0       hue, lightness, saturation corresponding to white",
    " bps=12                 bits per sample for color plots, either 12 or 24",
    " d1s=1.0                factor by which to scale d1 before imaging\t",
    " d2s=1.0                factor by which to scale d2 before imaging\t",
    " d3s=1.0                factor by which to scale d3 before imaging\t",
    " verbose=1              =1 for info printed on stderr (0 for no info)\t",
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
    " axescolor=black        color of axes\t\t\t\t\t",
    " gridcolor=black        color of grid\t\t\t\t\t",
    " legend=0               =1 display the color scale                     ",
    "                        if ==1, resize xbox,ybox,width,height          ",
    " lstyle=vertleft       Vertical, axis label on left side               ",
    "                        =vertright (Vertical, axis label on right side)",
    "                        =horibottom (Horizontal, axis label on bottom) ",
    " units=                 unit label for legend                          ",
    " legendfont=times_roman10    font name for title                       ",
    " following are defaults for lstyle=0. They are changed for other lstyles",
    " lwidth=1.2             colorscale (legend) width in inches            ",
    " lheight=height/3       colorscale (legend) height in inches           ",
    " lx=1.0                 colorscale (legend) x-position in inches       ",
    " ly=(height-lheight)/2+xybox    colorscale (legend) y-position in pixels",
    " lbeg= lmin or wclip-5*perc    value at which legend axis begins       ",
    " lend= lmax or bclip+5*perc    value at which legend axis ends         ",
    " ldnum=0.0      numbered tic interval on legend axis (0.0 for automatic)",
    " lfnum=lmin     first numbered tic on legend axis (used if d1num not 0.0)",
    " lntic=1        number of tics per numbered tic on legend axis ",
    " lgrid=none     grid lines on legend axis - none, dot, dash, or solid",
    " \t\t\t\t\t\t\t\t\t",
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
    " Note: The values of x1beg=x1min, x2end=x2max and x3beg=x3min cannot   ",
    " be changed.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

use oxwave_byte_interp::intl2b;
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_f_basic::*;
use oxwave_f_psaxesbox::psAxesBBox;
use oxwave_f_pscaxesbox::psCubeAxesBox;
use oxwave_f_pslegendbox::{psLegendBBox, psLegendBox};
use oxwave_sort::qkfind;
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut n1 = 0;
    let mut n2 = 0;
    let mut n3 = 0;
    let (mut d1, mut d2, mut d3) = (1.0f32, 1.0f32, 1.0f32);
    let (mut f1, mut f2, mut f3) = (0.0f32, 0.0f32, 0.0f32);
    if getparint("n1", &mut n1) == 0 {
        err!("must specify n1!\n");
    }
    getparfloat("d1", &mut d1);
    getparfloat("f1", &mut f1);
    if getparint("n2", &mut n2) == 0 {
        err!("must specify n2!\n");
    }
    getparfloat("d2", &mut d2);
    getparfloat("f2", &mut f2);
    if getparint("n3", &mut n3) == 0 {
        err!("must specify n3!\n");
    }
    getparfloat("d3", &mut d3);
    getparfloat("f3", &mut f3);
    let mut faces = 0;
    getparint("faces", &mut faces);
    let mut colors = [[0.0, 1.0]; 3];
    let mut hls = false;
    let mut bits = 8;
    if countparval("brgb") != 0 || countparval("wrgb") != 0 {
        let mut b = [0.0f32; 3];
        let mut w = [1.0f32; 3];
        getparfloat("brgb", &mut b[..]);
        getparfloat("wrgb", &mut w[..]);
        let mut bps = 12;
        getparint("bps", &mut bps);
        if bps != 12 && bps != 24 {
            err!("bps must equal 12 or 24 for color plots!\n");
        }
        bits = bps;
        for i in 0..3 {
            colors[i] = [b[i].clamp(0.0, 1.0), w[i].clamp(0.0, 1.0)];
        }
    } else if countparval("bhls") != 0 || countparval("whls") != 0 {
        let mut b = [0.0f32; 3];
        let mut w = [0.0f32, 1.0, 0.0];
        getparfloat("bhls", &mut b[..]);
        getparfloat("whls", &mut w[..]);
        let mut bps = 12;
        getparint("bps", &mut bps);
        if bps != 12 && bps != 24 {
            err!("bps must equal 12 or 24 for color plots!\n");
        }
        bits = bps;
        hls = true;
        for i in 0..3 {
            colors[i] = [b[i], w[i]];
        }
        for i in 1..3 {
            colors[i][0] = colors[i][0].clamp(0.0, 1.0);
            colors[i][1] = colors[i][1].clamp(0.0, 1.0);
        }
    }
    /* Allocate one contiguous array for the three faces. */
    let nz = n1
        .wrapping_mul(n2)
        .wrapping_add(n1.wrapping_mul(n3))
        .wrapping_add(n2.wrapping_mul(n3));
    let mut z = Vec::new();
    if z.try_reserve_exact(nz as usize).is_err() {
        // Vec rejects an impossible allocation before calling malloc. Preserve
        // ealloc1float's diagnostic, including syserr's ENOMEM suffix.
        err!(
            "%s: malloc failed (%s)",
            "ealloc.c",
            oxwave_core::err::strerror(12)
        );
    }
    z.resize(nz as usize, 0.0f32);
    let (zfront, remaining) = z.split_at_mut((n1 * n2) as usize);
    let (zside, ztop) = remaining.split_at_mut((n1 * n3) as usize);
    let mut front = String::new();
    let mut side = String::new();
    let mut top = String::new();
    if getparstring("front", &mut front) != 0
        && getparstring("side", &mut side) != 0
        && getparstring("top", &mut top) != 0
    {
        let frontfp = fopen(&front, "r").unwrap_or_else(|| err!("error opening front file!\n"));
        if fread(zfront, &frontfp) != (n1 * n2) as usize {
            err!("error reading front file!\n");
        }
        let sidefp = fopen(&side, "r").unwrap_or_else(|| err!("error opening side file!\n"));
        if fread(zside, &sidefp) != (n1 * n3) as usize {
            err!("error reading side file!\n");
        }
        let topfp = fopen(&top, "r").unwrap_or_else(|| err!("error opening top file!\n"));
        if fread(ztop, &topfp) != (n2 * n3) as usize {
            err!("error reading top file!\n");
        }
    } else if getparstring("front", &mut front) != 0
        || getparstring("side", &mut side) != 0
        || getparstring("top", &mut top) != 0
    {
        err!("must specify all or none of face, side, and top!\n");
    } else if faces != 0 {
        if fread(zfront, &stdin()) != (n1 * n2) as usize {
            err!("error reading front from input!\n");
        }
        if fread(zside, &stdin()) != (n1 * n3) as usize {
            err!("error reading side from input!\n");
        }
        if fread(ztop, &stdin()) != (n2 * n3) as usize {
            err!("error reading top from input!\n");
        }
    } else {
        // A zero-width cube still reads temp[0] in C. This build's fresh
        // zero-size allocation supplies a zero float; keep a guard for it.
        let mut temp = vec![0.0f32; (n1 as usize).max(1)];
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                if fread(&mut temp[..n1 as usize], &stdin()) != n1 as usize {
                    err!("error reading cube from input!\n");
                }
                if i3 == 0 {
                    for i1 in 0..n1 {
                        zfront[(i1 + i2 * n1) as usize] = temp[i1 as usize];
                    }
                }
                if i2 == n2 - 1 {
                    for i1 in 0..n1 {
                        zside[(i1 + i3 * n1) as usize] = temp[i1 as usize];
                    }
                }
                ztop[(i2 + i3 * n2) as usize] = temp[0];
            }
        }
    }
    render_faces(
        &z,
        [n1, n2, n3],
        [d1, d2, d3],
        [f1, f2, f3],
        hls,
        colors,
        bits,
    );
    exit(CWP_Exit());
}

// The arrays below group the three axes; operations retain the C statement order
// and single-precision intermediates.
fn render_faces(
    z: &[f32],
    n: [i32; 3],
    d: [f32; 3],
    f: [f32; 3],
    hls: bool,
    colors: [[f32; 2]; 3],
    bits: i32,
) {
    let nz = z.len();
    let (mut clip, mut bc, mut wc) = (0.0f32, 0.0f32, 0.0f32);
    if getparfloat("clip", &mut clip) != 0 {
        bc = clip;
        wc = -clip;
    }
    if (getparfloat("bclip", &mut bc) == 0 || getparfloat("wclip", &mut wc) == 0)
        && getparfloat("clip", &mut clip) == 0
    {
        let mut perc = 100.0f32;
        getparfloat("perc", &mut perc);
        let mut temp = z.to_vec();
        if getparfloat("bclip", &mut bc) == 0 {
            let mut bp = perc;
            getparfloat("bperc", &mut bp);
            let iz = ((nz as f32 * bp) as f64 / 100.0) as i32;
            let iz = iz.max(0).min(nz as i32 - 1) as usize;
            qkfind(iz, &mut temp);
            bc = temp[iz];
        }
        if getparfloat("wclip", &mut wc) == 0 {
            let mut wp = (100.0 - perc as f64) as f32;
            getparfloat("wperc", &mut wp);
            let iz = ((nz as f32 * wp) as f64 / 100.0) as i32;
            let iz = iz.max(0).min(nz as i32 - 1) as usize;
            qkfind(iz, &mut temp);
            wc = temp[iz];
        }
    }
    let mut verbose = 1;
    getparint("verbose", &mut verbose);
    if verbose != 0 {
        warn!("bclip=%g wclip=%g", bc, wc);
    }
    let mut ds = [1.0f32; 3];
    for i in 0..3 {
        getparfloat(&format!("d{}s", i + 1), &mut ds[i]);
        ds[i] = ds[i].abs();
        ds[i] *= d[i];
    }
    let mut angle = 45.0f32;
    getparfloat("angle", &mut angle);
    angle = angle.max(0.00001).min(90.0);
    angle = (angle as f64 * (std::f64::consts::PI / 180.0)) as f32;
    let mut size = [4.0f32, 4.0, 3.0];
    let (mut xbox, mut ybox) = (1.5f32, 1.5f32);
    for i in 0..3 {
        getparfloat(&format!("size{}", i + 1), &mut size[i]);
    }
    getparfloat("xbox", &mut xbox);
    getparfloat("ybox", &mut ybox);
    let mut beg = [0.0f32; 3];
    let mut end = [0.0f32; 3];
    let mut first = [0.0f32; 3];
    for i in 0..3 {
        let v = f[i] + (n[i] - 1) as f32 * d[i];
        beg[i] = if d[i] > 0.0 { f[i] } else { v };
        end[i] = if d[i] < 0.0 { f[i] } else { v };
        first[i] = beg[i];
    }
    let mut labels = [String::new(), String::new(), String::new()];
    let mut grids = [0; 3];
    let mut nums = [0.0f32; 3];
    let mut tics = [1; 3];
    for i in 0..3 {
        if i == 1 {
            getparfloat("x2beg", &mut beg[i]);
        } else {
            getparfloat(&format!("x{}end", i + 1), &mut end[i]);
        }
        getparfloat(&format!("d{}num", i + 1), &mut nums[i]);
        getparfloat(&format!("f{}num", i + 1), &mut first[i]);
        getparint(&format!("n{}tic", i + 1), &mut tics[i]);
        let mut g = "none".to_string();
        getparstring(&format!("grid{}", i + 1), &mut g);
        grids[i] = match g.as_str() {
            "dot" => 1,
            "dash" => 2,
            "solid" => 3,
            _ => 0,
        };
        getparstring(&format!("label{}", i + 1), &mut labels[i]);
    }
    let mut labelfont = "Helvetica".to_string();
    let mut titlefont = "Helvetica-Bold".to_string();
    let mut title = String::new();
    let mut titlecolor = "black".to_string();
    let mut axescolor = titlecolor.clone();
    let mut gridcolor = titlecolor.clone();
    let mut labelsize = 18.0f32;
    let mut titlesize = 24.0f32;
    getparstring("labelfont", &mut labelfont);
    getparfloat("labelsize", &mut labelsize);
    getparstring("title", &mut title);
    getparstring("titlefont", &mut titlefont);
    getparfloat("titlesize", &mut titlesize);
    getparstring("titlecolor", &mut titlecolor);
    getparstring("axescolor", &mut axescolor);
    getparstring("gridcolor", &mut gridcolor);
    let mut ib = [0; 3];
    let mut ie = [0; 3];
    let mut nc = [0; 3];
    let mut step = [0; 3];
    for i in 0..3 {
        ib[i] = nint(((beg[i] - f[i]) / d[i]) as f64)
            .min(if i == 0 { n[i] } else { n[i] - 1 })
            .max(0);
        beg[i] = f[i] + ib[i] as f32 * d[i];
        ie[i] = nint(((end[i] - f[i]) / d[i]) as f64).min(n[i] - 1).max(0);
        end[i] = f[i] + ie[i] as f32 * d[i];
        nc[i] = 1 + (ie[i] - ib[i]).abs();
        step[i] = if ie[i] > ib[i] { 1 } else { -1 };
    }
    let zs = if wc != bc {
        (255.0 / (wc - bc) as f64) as f32
    } else {
        1.0e10
    };
    let zo = -bc * zs;
    let offsets = [
        0,
        (n[0] * n[1]) as usize,
        (n[0] * n[1] + n[0] * n[2]) as usize,
    ];
    // With an empty face, C's snapped image is still one sample wide. Its
    // reads past the face allocation return zero in the original build.
    // Preserve that observed result with explicit guard storage.
    let guarded_z;
    let image_z = if n.contains(&0) {
        let mut v = z.to_vec();
        v.resize(nz + n[0].max(n[1]).max(1) as usize, 0.0);
        guarded_z = v;
        &guarded_z[..]
    } else {
        z
    };
    let axes = [(0, 1), (0, 2), (1, 2)];
    let mut cz = [Vec::new(), Vec::new(), Vec::new()];
    for (face, &(a, b)) in axes.iter().enumerate() {
        let mut j = ib[b];
        for _ in 0..nc[b] {
            let mut k = ib[a];
            for _ in 0..nc[a] {
                let mut zi = zo + image_z[offsets[face] + (k + j * n[a]) as usize] * zs;
                if zi < 0.0 {
                    zi = 0.0;
                }
                if zi > 255.0 {
                    zi = 255.0;
                }
                cz[face].push(zi as u8);
                k += step[a];
            }
            j += step[b];
        }
    }
    let mut ns = [0; 3];
    for i in 0..3 {
        ns[i] = nint((1.0 + (nc[i] - 1) as f32 * d[i] / ds[i]) as f64).max(1);
        ds[i] = if ns[i] > 1 {
            d[i] * (nc[i] - 1) as f32 / (ns[i] - 1) as f32
        } else {
            d[i]
        };
    }
    for (face, &(a, b)) in axes.iter().enumerate() {
        if ns[a] != nc[a] || ns[b] != nc[b] {
            let mut scaled = vec![0; (ns[a] * ns[b]) as usize];
            intl2b(
                nc[a],
                d[a],
                0.0,
                nc[b],
                d[b],
                0.0,
                &cz[face],
                ns[a],
                ds[a],
                0.0,
                ns[b],
                ds[b],
                0.0,
                &mut scaled,
            );
            cz[face] = scaled;
        }
    }
    let mut pb = [0.0f32; 3];
    let mut pe = [0.0f32; 3];
    for i in 0..3 {
        pb[i] = if end[i] > beg[i] {
            -ds[i].abs() / 2.0
        } else {
            ds[i].abs() / 2.0
        };
        pe[i] = -pb[i];
    }
    let mut legend = 0;
    getparint("legend", &mut legend);
    let mut units = String::new();
    let mut legendfont = "times_roman10".to_string();
    getparstring("units", &mut units);
    getparstring("legendfont", &mut legendfont);
    let (ca, sa, ta) = (
        (angle as f64).cos(),
        (angle as f64).sin(),
        (angle as f64).tan(),
    );
    let (mut lstyle, mut ugrid) = (0, 3);
    let (mut lx, mut ly, mut lw, mut lh, mut lbeg, mut lend, mut lf, mut ldnum) = (
        0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32,
    );
    let mut data_legend = Vec::new();
    let mut ln = 256;
    if legend != 0 {
        let (mut lmin, mut lmax) = (f32::MAX, -f32::MAX);
        for &v in z {
            lmin = if lmin > v { v } else { lmin };
            lmax = if lmax > v { lmax } else { v };
        }
        if verbose == 2 {
            warn!("lmin=%g lmax=%g", lmin, lmax);
        }
        lbeg = lmin;
        let lbegsup = getparfloat("lbeg", &mut lbeg) != 0;
        lend = lmax;
        let lendsup = getparfloat("lend", &mut lend) != 0;
        wc = if lbeg > wc { lbeg } else { wc };
        bc = if lend > bc { bc } else { lend };
        if !lbegsup {
            let rangeperc = ((bc - wc) as f64 / 20.0) as f32;
            lbeg = wc - rangeperc;
        }
        if !lendsup {
            let rangeperc = ((bc - wc) as f64 / 20.0) as f32;
            lend = bc + rangeperc;
        }
        let mut lfnum = lmin;
        getparfloat("lfnum", &mut lfnum);
        let mut lstyles = "vertleft".to_string();
        getparstring("lstyle", &mut lstyles);
        if lstyles == "vertright" {
            lstyle = 1;
        } else if lstyles == "horibottom" {
            lstyle = 2;
        }
        lw = 0.1;
        lh = (size[0] as f64 + sa * size[2] as f64 / 2.0) as f32;
        if lstyle == 2 {
            lw = (size[1] as f64 + ca * size[2] as f64 / 1.2) as f32;
            lh = 0.24;
        }
        getparfloat("lwidth", &mut lw);
        getparfloat("lheight", &mut lh);
        lx = 0.8;
        ly = (ybox as f64 + (size[0] as f64 + sa * size[2] as f64 - lh as f64) / 2.0) as f32;
        if lstyle == 1 {
            lx = ((xbox + size[1]) as f64 + ca * size[2] as f64 + 0.1) as f32;
        } else if lstyle == 2 {
            lx = (xbox as f64 + (size[1] as f64 + ca * size[2] as f64 - lw as f64) / 2.0) as f32;
            ly = 1.0;
        }
        getparfloat("lx", &mut lx);
        getparfloat("ly", &mut ly);
        let mut lgrid = "none".to_string();
        getparstring("lgrid", &mut lgrid);
        ugrid = match lgrid.as_str() {
            "dot" => 1,
            "dash" => 2,
            "solid" => 3,
            _ => 0,
        };
        if !lbegsup {
            ln += 13;
            data_legend.extend_from_slice(&[255; 13]);
        }
        if !lendsup {
            ln += 13;
        }
        for v in (0..=255).rev() {
            data_legend.push(v as u8);
        }
        data_legend.resize(ln as usize, 0);
        lf = lbeg;
        let ld = (lend - lbeg) / (ln - 1) as f32;
        getparfloat("ldnum", &mut ldnum);
        checkpars();
        let ilbeg = nint(((lbeg - lf) / ld) as f64).max(0).min(ln - 1);
        lbeg = lf + ilbeg as f32 * ld;
        let ilend = nint(((lend - lf) / ld) as f64).max(0).min(ln - 1);
        lend = lf + ilend as f32 * ld;
        lx *= 72.0;
        ly *= 72.0;
        lw *= 72.0;
        lh *= 72.0;
    }
    for s in &mut size {
        *s *= 72.0;
    }
    xbox *= 72.0;
    ybox *= 72.0;
    let mut bbox = [0; 4];
    check_metric_font(&labelfont);
    check_metric_font(&titlefont);
    psAxesBBox(
        xbox,
        ybox,
        (size[1] as f64 + ca * size[2] as f64) as f32,
        (size[0] as f64 + sa * size[2] as f64) as f32,
        &labelfont,
        labelsize,
        &titlefont,
        titlesize,
        1,
        &mut bbox,
    );
    if legend != 0 {
        let mut lbbox = [0; 4];
        psLegendBBox(lx, ly, lw, lh, &labelfont, labelsize, lstyle, &mut lbbox);
        for i in 0..2 {
            bbox[i] = bbox[i].min(lbbox[i]);
            bbox[i + 2] = bbox[i + 2].max(lbbox[i + 2]);
        }
    }
    boundingbox(bbox[0], bbox[1], bbox[2], bbox[3]);
    begineps();
    gsave();
    translate(xbox, ybox);
    gsave();
    translate(0.0, 0.0);
    scale(size[1], size[0]);
    draw_image(
        hls,
        colors,
        bits,
        ns[0],
        ns[1],
        &[0.0, ns[1] as f32, -(ns[0] as f32), 0.0, ns[0] as f32, 0.0],
        &cz[0],
    );
    grestore();
    gsave();
    translate(size[1], 0.0);
    concat(&[1.0, ta as f32, 0.0, 1.0, 0.0, 0.0]);
    scale((size[2] as f64 * ca) as f32, size[0]);
    draw_image(
        hls,
        colors,
        bits,
        ns[0],
        ns[2],
        &[0.0, ns[2] as f32, -(ns[0] as f32), 0.0, ns[0] as f32, 0.0],
        &cz[1],
    );
    grestore();
    gsave();
    translate(0.0, size[0]);
    concat(&[1.0, 0.0, (1.0 / ta) as f32, 1.0, 0.0, 0.0]);
    scale(size[1], (size[2] as f64 * sa) as f32);
    draw_image(
        hls,
        colors,
        bits,
        ns[1],
        ns[2],
        &[ns[1] as f32, 0.0, 0.0, ns[2] as f32, 0.0, 0.0],
        &cz[2],
    );
    grestore();
    if legend != 0 {
        gsave();
        translate(-xbox, -ybox);
        translate(lx, ly);
        scale(lw, lh);
        if lstyle == 0 || lstyle == 1 {
            draw_image(
                hls,
                colors,
                bits,
                1,
                ln,
                &[1.0, 0.0, 0.0, ln as f32, 0.0, 0.0],
                &data_legend,
            );
        } else {
            rotate(-90.0);
            draw_image(
                hls,
                colors,
                bits,
                1,
                ln,
                &[-1.0, 0.0, 0.0, ln as f32, 0.0, 0.0],
                &data_legend,
            );
            rotate(90.0);
        }
        grestore();
    }
    grestore();
    psCubeAxesBox(
        xbox,
        ybox,
        size[0],
        size[1],
        size[2],
        angle,
        beg[0],
        end[0],
        pb[0],
        pe[0],
        nums[0],
        first[0],
        tics[0],
        grids[0],
        &labels[0],
        beg[1],
        end[1],
        pb[1],
        pe[1],
        nums[1],
        first[1],
        tics[1],
        grids[1],
        &labels[1],
        beg[2],
        end[2],
        pb[2],
        pe[2],
        nums[2],
        first[2],
        tics[2],
        grids[2],
        &labels[2],
        &labelfont,
        labelsize,
        &title,
        &titlefont,
        titlesize,
        &titlecolor,
        &axescolor,
        &gridcolor,
    );
    if legend != 0 {
        gsave();
        psLegendBox(
            lx, ly, lw, lh, lbeg, lend, 0.0, 0.0, ldnum, lf, 1, ugrid, &units, &labelfont,
            labelsize, &axescolor, &gridcolor, lstyle,
        );
        grestore();
    }
    showpage();
    endeps();
}

// basic.c declares 38 metrics but initializes only 37. Unknown names reach
// strcmp with the last entry's null name. Reject them as the fixed C does,
// before the shared font-metrics library can panic.
fn check_metric_font(name: &str) {
    match name {
        "AvantGarde-Book"
        | "AvantGarde-BookOblique"
        | "AvantGarde-Demi"
        | "AvantGarde-DemiOblique"
        | "Bookman-Demi"
        | "Bookman-DemiItalic"
        | "Bookman-Light"
        | "Bookman-LightItalic"
        | "Courier"
        | "Courier-Bold"
        | "Courier-BoldOblique"
        | "Courier-Oblique"
        | "Helvetica"
        | "Helvetica-Bold"
        | "Helvetica-BoldOblique"
        | "Helvetica-Oblique"
        | "Helvetica-Narrow"
        | "Helvetica-Narrow-Bold"
        | "Helvetica-Narrow-BoldOblique"
        | "Helvetica-Narrow-Oblique"
        | "NewCentrySchlbk-Bold"
        | "NewCenturySchlbk-BoldItalic"
        | "NewCenturySchlbk-Roman"
        | "Palatino-Bold"
        | "Palatino-BoldItalic"
        | "Palatino-Italics"
        | "Palatino-Roman"
        | "SanSerif-Bold"
        | "SanSerif-BoldItalic"
        | "SanSerif-Roman"
        | "Symbol"
        | "Times-Bold"
        | "Times-BoldItalic"
        | "Times-Roman"
        | "Times-Italic"
        | "ZapfChancery-MediumItalic"
        | "Unknown" => (),
        _ => err!("unsupported font: %s", name),
    }
}

fn rgb_value(n1: f32, n2: f32, mut hue: f32) -> f32 {
    while hue > 360.0 {
        hue -= 360.0;
    }
    while hue < 0.0 {
        hue += 360.0;
    }
    if hue < 60.0 {
        ((n1 as f64) + ((n2 - n1) * hue) as f64 / 60.0) as f32
    } else if hue < 180.0 {
        n2
    } else if hue < 240.0 {
        (n1 as f64 + (n2 - n1) as f64 * (240.0 - hue as f64) / 60.0) as f32
    } else {
        n1
    }
}
fn draw_image(hls: bool, c: [[f32; 2]; 3], bits: i32, w: i32, h: i32, m: &[f32], z: &[u8]) {
    if (hls && c[1] == [0.0, 1.0] && c[2] == [0.0, 0.0]) || c == [[0.0, 1.0]; 3] {
        image(w, h, 8, m, z);
        return;
    }
    let mut tab = [[0u8; 256]; 3];
    for i in 0..256 {
        if hls {
            let mut x = [0.0; 3];
            for j in 0..3 {
                let scale = ((c[j][1] - c[j][0]) as f64 / 255.0) as f32;
                x[j] = c[j][0] + i as f32 * scale;
            }
            let (h, l, s) = ((x[0] as f64 * 360.0) as f32, x[1], x[2]);
            let m2 = if l <= 0.5 {
                (l as f64 * (1.0 + s as f64)) as f32
            } else {
                l + s - l * s
            };
            let m1 = 2.0 * l - m2;
            let v = if s == 0.0 {
                [l, l, l]
            } else {
                [
                    rgb_value(m1, m2, (h as f64 + 120.0) as f32),
                    rgb_value(m1, m2, h),
                    rgb_value(m1, m2, (h as f64 - 120.0) as f32),
                ]
            };
            for j in 0..3 {
                tab[j][i] = (v[j] as f64 * 255.0) as u8;
            }
        } else {
            for j in 0..3 {
                let base = (c[j][0] as f64 * 255.0) as f32;
                tab[j][i] = (base + i as f32 * (c[j][1] - c[j][0])) as u8;
            }
        }
    }
    let mut rgb = Vec::with_capacity(z.len() * 3);
    for &v in z {
        for table in &tab {
            rgb.push(table[v as usize]);
        }
    }
    if bits == 12 {
        let rowbytes = (1 + (w * 12 - 1) / 8) as usize;
        let mut packed = Vec::new();
        let mut k = 0;
        for _ in 0..h {
            for _ in 0..rowbytes - 1 {
                packed.push((rgb[k] & 0xf0) | (rgb[k + 1] >> 4));
                k += 2;
            }
            if w % 2 != 0 {
                packed.push(rgb[k] & 0xf0);
                k += 1;
            } else {
                packed.push((rgb[k] & 0xf0) | (rgb[k + 1] >> 4));
                k += 2;
            }
        }
        rgb = packed;
    }
    rgbimage(w, h, bits / 3, m, &rgb);
}
