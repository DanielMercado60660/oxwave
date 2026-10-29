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
//! PSLABEL - output a PostScript text string. Port of psplot/main/pslabel.c (Revision 1.7).
//! Author: John E. Anderson, Visiting Scientist from Mobil, 1994.

use oxwave_core::prelude::*;
use oxwave_f_basic::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " PSLABEL - output PostScript file consisting of a single TEXT string\t",
    "          on a specified background. (Use with psmerge to label plots.)",
    " \t\t\t\t\t\t\t\t\t",
    " pslabel t= [t=] [optional parameters] > epsfile\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "Required Parameters (can have multiple specifications to mix fonts):\t",
    "  t=                 text string to write to output\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "Optional Parameters:\t\t\t\t\t\t\t",
    "  f=Times-Bold       font for text string\t\t\t\t",
    "                      (multiple specifications for each t)\t\t",
    "  size=30            size of characters in points (72 points/inch)\t",
    "  tcolor=black       color of text string\t\t\t\t",
    "  bcolor=white       color of background box\t\t\t\t",
    "  nsub=0             number of characters to subtract when\t\t",
    "                     computing size of background box (not all\t\t",
    "                     characters are the same size so the\t\t",
    "                     background box may be too big at times.)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t\t",
    " pslabel t=\"(a) \" f=Times-Bold t=\"h\" f=Symbol t=\"=0.04\" nsub=3 > epsfile",
    "\t\t\t\t\t\t\t\t\t",
    " This example yields the PostScript equivalent of the string\t\t",
    "  (written here in LaTeX notation) $ (a)\\; \\eta=0.04 $\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program produces a (color if desired) PostScript text string that",
    " can be positioned and pasted on a PostScript plot using   psmerge \t",
    "     see selfdoc of   psmerge for further information.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Possible fonts:   Helvetica, Helvetica-Oblique, Helvetica-Bold,\t",
    "  Helvetica-BoldOblique,Times-Roman,Times-Italic,Times-Bold,\t\t",
    "  Times-BoldItalic,Courier,Courier-Bold,Courier-Oblique,\t\t",
    "  Courier-BoldOblique,Symbol\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Possible colors:  greenyellow,yellow,goldenrod,dandelion,apricot,\t",
    "  peach,melon,yelloworange,orange,burntorange,bittersweet,redorange,\t",
    "  mahogany,maroon,brickred,red,orangered,rubinered,wildstrawberry,\t",
    "  salmon,carnationpink,magenta,violetred,rhodamine,mulberry,redviolet,\t",
    "  fuchsia,lavender,thistle,orchid,darkorchid,purple,plum,violet,royalpurple,",
    "  blueviolet,periwinkle,cadetblue,cornflowerblue,midnightblue,naveblue,",
    "  royalblue,blue,cerulean,cyan,processblue,skyblue,turquoise,tealblue,\t",
    "  aquamarine,bluegreen,emerald,junglegreen,seagreen,green,forestgreen,\t",
    "  pinegreen,limegreen,yellowgreen,springgreen,olivegreen,rawsienna,sepia,",
    "  brown,tan,white,black,gray\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " All color specifications may also be made in X Window style Hex format",
    " example:   tcolor=#255\t\t\t\t\t\t",
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

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let mut nsub = 0;
    let mut size = 30.0f32;
    let mut textcolor = String::from("black");
    let mut boxcolor = String::from("white");
    getparint("nsub", &mut nsub);
    getparfloat("size", &mut size);
    getparstring("tcolor", &mut textcolor);
    getparstring("bcolor", &mut boxcolor);
    checkpars();
    let eps = (0.25 * size as f64) as f32;
    let eps2 = (0.1 * size as f64) as f32;
    let n = countparname("t");
    if n == 0 { err!("must enter at least one PSTEXT text stream as parameter t"); }
    if n != countparname("f") { warn!("suggest specify same number of values for t and f"); }
    let mut text = Vec::new();
    let mut font = Vec::new();
    let mut x = Vec::new();
    let mut bigx = eps;
    let mut bigy = 0.0f32;
    let mut label_cd = 0.0f32;
    let mut label_cw = 0.0f32;
    for j in 0..n {
        x.push(bigx);
        let mut t = String::new();
        let mut f = String::new();
        if getnparstring(j + 1, "t", &mut t) == 0 { t = "hello".into(); }
        if getnparstring(j + 1, "f", &mut f) == 0 { f = "Times-Bold".into(); }
        let label_ch = fontheight(&f, size);
        label_cw = fontwidth(&f, size);
        let label_ca = fontascender(&f, size);
        let descender = fontdescender(&f, size);
        label_cd = if label_cd < descender { label_cd } else { descender };
        let nchar = t.len();
        bigx = (bigx as f64 + 0.5 * (nchar as f64 * label_cw as f64)) as f32;
        let height = (label_ch + eps) as f64 + 0.0 * label_ca as f64;
        bigy = if bigy as f64 > height { bigy } else { height as f32 };
        text.push(t);
        font.push(f);
    }
    bigx += eps;
    bigx = (bigx as f64 - 0.5 * nsub as f64 * label_cw as f64) as f32;
    boundingbox((-eps2) as i32, (-eps2) as i32, (bigx + eps2) as i32, (bigy + eps2) as i32);
    begineps();
    gsave();
    rectclip(0.0, 0.0, bigx, bigy);
    newpath();
    moveto(0.0, 0.0);
    lineto(bigx, 0.0);
    lineto(bigx, bigy);
    lineto(0.0, bigy);
    closepath();
    setcolor(&boxcolor);
    fill();
    setcolor(&textcolor);
    moveto(eps, eps - label_cd);
    for j in 0..n as usize {
        setfont(&font[j], size);
        show(&text[j]);
    }
    grestore();
    showpage();
    endeps();
    exit(CWP_Exit());
}
