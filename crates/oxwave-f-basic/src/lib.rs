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

//! Port of psplot/lib/basic.c.
//! Credits: Dave Hale, Craig Artley, and Robert Krueger.

use oxwave_core::printf;
use std::sync::Mutex;

struct State {
    bbox: [i32; 4],
    set: bool,
    out: bool,
    pages: i32,
}
static STATE: Mutex<State> = Mutex::new(State {
    bbox: [0, 0, 612, 792],
    set: false,
    out: false,
    pages: 0,
});
const PROLOGUE: &str = "/M {moveto} def\n/RM {rmoveto} def\n/L {lineto} def\n/RL {rlineto} def\n/S {stroke} def\n/F {fill} def\n/GS {gsave} def\n/GR {grestore} def\n/SH {show} def\n/SW {stringwidth} def\n/NP {newpath} def\n/CP {closepath} def\n/SC {scale} def\n/RO {rotate} def\n/TR {translate} def\n/CAT {concat} def\n/CLW {currentlinewidth} def\n/SLW {setlinewidth} def\n";
fn begin(pages: bool) {
    let mut s = STATE.lock().unwrap();
    printf!("%%!PS-Adobe-2.0 EPSF-1.2\n%%%%DocumentFonts:\n");
    if s.set {
        let b = s.bbox;
        printf!("%%%%BoundingBox: %d %d %d %d\n", b[0], b[1], b[2], b[3]);
        s.out = true;
    } else {
        printf!("%%%%BoundingBox: (atend)\n");
    }
    if pages {
        printf!("%%%%Pages: (atend)\n");
    }
    printf!("%%%%EndComments\n%s\n%%%%EndProlog\nGS\n", PROLOGUE);
}
fn end(pages: bool) {
    let mut s = STATE.lock().unwrap();
    printf!("GR\n%%%%Trailer\n");
    if !s.out {
        let b = s.bbox;
        printf!("%%%%BoundingBox: %d %d %d %d\n", b[0], b[1], b[2], b[3]);
        s.out = true;
    }
    if pages {
        printf!("%%%%Pages: %d\n", s.pages);
    }
}
/// C: `void beginps(void);`
pub fn beginps() {
    begin(true);
}
/// C: `void begineps(void);`
pub fn begineps() {
    begin(false);
}
/// C: `void endps(void);`
pub fn endps() {
    end(true);
}
/// C: `void endeps(void);`
pub fn endeps() {
    end(false);
}
/// C: `void boundingbox(int llx, int lly, int urx, int ury);`
pub fn boundingbox(llx: i32, lly: i32, urx: i32, ury: i32) {
    let mut s = STATE.lock().unwrap();
    s.bbox = [llx, lly, urx, ury];
    s.set = true;
}
/// C: `void newpage(const char *label, int ordinal);`
pub fn newpage(label: &str, ordinal: i32) {
    printf!("%%%%Page: %s %d\n", label, ordinal);
    let mut s = STATE.lock().unwrap();
    s.pages = s.pages.wrapping_add(1);
}
/// C: `void showpage(void);`
pub fn showpage() {
    printf!("showpage\n");
}
/// C: `void gsave(void);`
pub fn gsave() {
    printf!("GS\n");
}
/// C: `void grestore(void);`
pub fn grestore() {
    printf!("GR\n");
}
/// C: `void newpath(void);`
pub fn newpath() {
    printf!("NP\n");
}
/// C: `void closepath(void);`
pub fn closepath() {
    printf!("CP\n");
}
/// C: `void clip(void);`
pub fn clip() {
    printf!("clip\n");
}
/// C: `void stroke(void);`
pub fn stroke() {
    printf!("S\n");
}
/// C: `void fill(void);`
pub fn fill() {
    printf!("F\n");
}
/// C: `void translate(float tx, float ty);`
pub fn translate(x: f32, y: f32) {
    printf!("%.7g %.7g TR\n", x, y);
}
/// C: `void scale(float sx, float sy);`
pub fn scale(x: f32, y: f32) {
    printf!("%.4g %.4g SC\n", x, y);
}
/// C: `void moveto(float x, float y);`
pub fn moveto(x: f32, y: f32) {
    printf!("%.4g %.4g M\n", x, y);
}
/// C: `void rmoveto(float x, float y);`
pub fn rmoveto(x: f32, y: f32) {
    printf!("%.4g %.4g RM\n", x, y);
}
/// C: `void lineto(float x, float y);`
pub fn lineto(x: f32, y: f32) {
    printf!("%.4g %.4g L\n", x, y);
}
/// C: `void rlineto(float x, float y);`
pub fn rlineto(x: f32, y: f32) {
    printf!("%.4g %.4g RL\n", x, y);
}
/// C: `void rotate(float angle);`
pub fn rotate(x: f32) {
    printf!("%.4g RO\n", x);
}
/// C: `void setgray(float gray);`
pub fn setgray(x: f32) {
    printf!("%.4g setgray\n", x);
}
/// C: `void setlinewidth(float width);`
pub fn setlinewidth(x: f32) {
    printf!("%.4g SLW\n", x);
}

/// C: `void concat(float m[]);`
pub fn concat(m: &[f32]) {
    printf!(
        "[%.4g %.4g %.4g %.4g %.4g %.4g] CAT\n",
        m[0],
        m[1],
        m[2],
        m[3],
        m[4],
        m[5]
    );
}
/// C: `void setrgbcolor(float red, float green, float blue);`
pub fn setrgbcolor(r: f32, g: f32, b: f32) {
    printf!("%.4g %.4g %.4g setrgbcolor\n", r, g, b);
}
/// C: `void setcymkcolor(float cyan, float magenta, float yellow, float black);`
pub fn setcymkcolor(c: f32, m: f32, y: f32, k: f32) {
    printf!("%.4g %.4g %.4g %.4g setcymkcolor\n", c, m, y, k);
}
/// C: `void setlinejoin(int code);`
pub fn setlinejoin(code: i32) {
    printf!("%d setlinejoin\n", code);
}
/// C: `void setdash(float dash[], int ndash, float offset);`
pub fn setdash(dash: &[f32], offset: f32) {
    printf!("[ ");
    for x in dash {
        printf!("%.4g ", x);
    }
    printf!("] %.4g setdash\n", offset);
}
/// C: `void arc(float x, float y, float r, float ang1, float ang2);`
pub fn arc(x: f32, y: f32, r: f32, a: f32, b: f32) {
    printf!("%.4g %.4g %.4g %.4g %.4g arc\n", x, y, r, a, b);
}
/// C: `void show(const char *str);`
pub fn show(s: &str) {
    printf!("(%s) SH\n", s);
}
/// C: `void justshow(float just, const char *str);`
pub fn justshow(j: f32, s: &str) {
    printf!(
        "(%s) SW exch %.4g mul\nexch %.4g mul RM (%s) SH\n",
        s,
        j,
        j,
        s
    );
}
/// C: `void setfont(const char *fontname, float fontsize);`
pub fn setfont(name: &str, size: f32) {
    printf!("/%s findfont %.4g scalefont setfont\n", name, size);
}
const FONTS: &[(&str, [i32; 4], i32, i32, i32, i32)] = &[
    (
        "AvantGarde-Book",
        [-113, -222, 1148, 955],
        740,
        547,
        -192,
        740,
    ),
    (
        "AvantGarde-BookOblique",
        [-113, -222, 1279, 955],
        740,
        547,
        -192,
        740,
    ),
    (
        "AvantGarde-Demi",
        [-123, -251, 1222, 1021],
        740,
        555,
        -185,
        740,
    ),
    (
        "AvantGarde-DemiOblique",
        [-123, -251, 1256, 1021],
        740,
        555,
        -185,
        740,
    ),
    ("Bookman-Demi", [-194, -250, 1346, 934], 681, 502, -212, 725),
    (
        "Bookman-DemiItalic",
        [-231, -250, 1333, 941],
        681,
        515,
        -213,
        732,
    ),
    (
        "Bookman-Light",
        [-188, -251, 1266, 908],
        681,
        484,
        -228,
        717,
    ),
    (
        "Bookman-LightItalic",
        [-228, -250, 1269, 883],
        681,
        494,
        -212,
        717,
    ),
    ("Courier", [-28, -250, 628, 805], 562, 426, -157, 629),
    ("Courier-Bold", [-113, -250, 749, 801], 562, 439, -142, 626),
    (
        "Courier-BoldOblique",
        [-56, -250, 868, 801],
        562,
        439,
        -142,
        626,
    ),
    (
        "Courier-Oblique",
        [-28, -250, 742, 805],
        562,
        426,
        -157,
        629,
    ),
    ("Helvetica", [-166, -225, 1000, 931], 718, 523, -207, 718),
    (
        "Helvetica-Bold",
        [-170, -228, 1003, 962],
        718,
        532,
        -207,
        718,
    ),
    (
        "Helvetica-BoldOblique",
        [-174, -228, 1114, 962],
        718,
        532,
        -207,
        718,
    ),
    (
        "Helvetica-Oblique",
        [-170, -225, 1116, 931],
        718,
        523,
        -207,
        718,
    ),
    (
        "Helvetica-Narrow",
        [-136, -225, 820, 931],
        718,
        523,
        -207,
        718,
    ),
    (
        "Helvetica-Narrow-Bold",
        [-139, -228, 822, 962],
        718,
        532,
        -207,
        718,
    ),
    (
        "Helvetica-Narrow-BoldOblique",
        [-143, -228, 913, 962],
        718,
        532,
        -207,
        718,
    ),
    (
        "Helvetica-Narrow-Oblique",
        [-139, -225, 915, 931],
        718,
        523,
        -207,
        228,
    ),
    (
        "NewCentrySchlbk-Bold",
        [-165, -250, 1000, 988],
        722,
        475,
        -205,
        737,
    ),
    (
        "NewCenturySchlbk-BoldItalic",
        [-205, -250, 1147, 991],
        722,
        477,
        -204,
        737,
    ),
    (
        "NewCenturySchlbk-Roman",
        [-195, -250, 1000, 965],
        722,
        464,
        -205,
        737,
    ),
    (
        "Palatino-Bold",
        [-152, -266, 1000, 924],
        681,
        471,
        -258,
        720,
    ),
    (
        "Palatino-BoldItalic",
        [-170, -271, 1073, 926],
        681,
        469,
        -271,
        726,
    ),
    (
        "Palatino-Italics",
        [-170, -276, 1010, 918],
        692,
        482,
        -276,
        733,
    ),
    (
        "Palatino-Roman",
        [-166, -283, 1021, 927],
        692,
        469,
        -281,
        726,
    ),
    (
        "SanSerif-Bold",
        [-165, -250, 1000, 988],
        722,
        475,
        -205,
        737,
    ),
    (
        "SanSerif-BoldItalic",
        [-205, -250, 1147, 991],
        722,
        477,
        -204,
        737,
    ),
    (
        "SanSerif-Roman",
        [-195, -250, 1000, 965],
        722,
        464,
        -205,
        737,
    ),
    ("Symbol", [-180, -293, 1090, 1010], 729, 525, -219, 729),
    ("Times-Bold", [-168, -218, 1000, 935], 676, 461, -205, 676),
    (
        "Times-BoldItalic",
        [-200, -218, 996, 921],
        669,
        462,
        -205,
        699,
    ),
    ("Times-Roman", [-168, -218, 1000, 898], 662, 450, -217, 683),
    ("Times-Italic", [-169, -217, 1010, 883], 653, 441, -205, 683),
    (
        "ZapfChancery-MediumItalic",
        [-181, -314, 1065, 831],
        708,
        438,
        -314,
        714,
    ),
    ("Unknown", [-174, -220, 1001, 944], 729, 525, -219, 729),
];
fn font(name: &str) -> &'static (&'static str, [i32; 4], i32, i32, i32, i32) {
    FONTS
        .iter()
        .find(|f| f.0 == name)
        .expect("unknown font: the C reference has undefined behavior for this input")
}
/// C: `void fontbbox(const char *fontname, float fontsize, float bbox[]);`
pub fn fontbbox(name: &str, size: f32, bbox: &mut [f32]) {
    let b = font(name).1;
    for i in 0..4 {
        bbox[i] = (b[i] as f64 * 0.001 * size as f64) as f32;
    }
}
/// C: `float fontheight(const char *fontname, float fontsize);`
pub fn fontheight(name: &str, size: f32) -> f32 {
    let b = font(name).1;
    ((b[3] - b[1]) as f64 * (0.001 * size as f64)) as f32
}
/// C: `float fontwidth(const char *fontname, float fontsize);`
pub fn fontwidth(name: &str, size: f32) -> f32 {
    let b = font(name).1;
    ((b[2] - b[0]) as f64 * (0.001 * size as f64)) as f32
}
/// C: `float fontcapheight(const char *fontname, float fontsize);`
pub fn fontcapheight(name: &str, size: f32) -> f32 {
    (font(name).2 as f64 * 0.001 * size as f64) as f32
}
/// C: `float fontxheight(const char *fontname, float fontsize);`
pub fn fontxheight(name: &str, size: f32) -> f32 {
    (font(name).3 as f64 * 0.001 * size as f64) as f32
}
/// C: `float fontdescender(const char *fontname, float fontsize);`
pub fn fontdescender(name: &str, size: f32) -> f32 {
    (font(name).4 as f64 * 0.001 * size as f64) as f32
}
/// C: `float fontascender(const char *fontname, float fontsize);`
pub fn fontascender(name: &str, size: f32) -> f32 {
    (font(name).5 as f64 * 0.001 * size as f64) as f32
}
const COLORS: &[(&str, [f32; 3])] = &[
    ("greenyellow", [0.85, 0.0, 0.31]),
    ("yellow", [1.0, 1.0, 0.0]),
    ("goldenrod", [1.0, 0.9, 0.16]),
    ("dandelion", [1.0, 0.71, 0.16]),
    ("apricot", [1.0, 0.68, 0.48]),
    ("peach", [1.0, 0.50, 0.30]),
    ("melon", [1.0, 0.54, 0.50]),
    ("yelloworange", [1.0, 0.58, 0.0]),
    ("orange", [1.0, 0.39, 0.13]),
    ("burntorange", [1.0, 0.49, 0.0]),
    ("bittersweet", [0.76, 0.01, 0.0]),
    ("redorange", [1.0, 0.23, 0.13]),
    ("mahogany", [0.65, 0.0, 0.0]),
    ("maroon", [0.68, 0.0, 0.0]),
    ("brickred", [1.0, 0.0, 0.0]),
    ("red", [1.0, 0.0, 0.0]),
    ("orangered", [1.0, 0.0, 0.50]),
    ("rubinered", [1.0, 0.0, 0.87]),
    ("wildstrawberry", [1.0, 0.04, 0.61]),
    ("salmon", [1.0, 0.47, 0.62]),
    ("carnationpink", [1.0, 0.37, 1.0]),
    ("magenta", [1.0, 0.0, 1.0]),
    ("violetred", [1.0, 0.19, 1.0]),
    ("rhodamine", [1.0, 0.18, 1.0]),
    ("mulberry", [0.64, 0.08, 0.98]),
    ("redviolet", [0.59, 0.0, 0.64]),
    ("fuchsia", [0.45, 0.01, 0.92]),
    ("lavender", [1.0, 0.52, 1.0]),
    ("thistle", [0.88, 0.41, 1.0]),
    ("orchid", [0.68, 0.36, 1.0]),
    ("darkorchid", [0.60, 0.20, 0.80]),
    ("purple", [0.55, 0.14, 1.0]),
    ("plum", [0.50, 0.0, 1.0]),
    ("violet", [0.21, 0.12, 1.0]),
    ("royalpurple", [0.25, 0.10, 1.0]),
    ("blueviolet", [0.10, 0.05, 0.96]),
    ("periwinkle", [0.43, 0.45, 1.0]),
    ("cadetblue", [0.38, 0.43, 0.77]),
    ("cornflowerblue", [0.35, 0.87, 1.0]),
    ("midnightblue", [0.0, 0.44, 0.57]),
    ("naveblue", [0.06, 0.46, 1.0]),
    ("royalblue", [0.0, 0.50, 1.0]),
    ("blue", [0.0, 0.0, 1.0]),
    ("cerulean", [0.06, 0.89, 1.0]),
    ("cyan", [0.0, 1.0, 1.0]),
    ("processblue", [0.04, 1.0, 1.0]),
    ("skyblue", [0.38, 1.0, 0.88]),
    ("turquoise", [0.15, 1.0, 0.80]),
    ("tealblue", [0.12, 0.98, 0.64]),
    ("aquamarine", [0.18, 1.0, 0.70]),
    ("bluegreen", [0.15, 1.0, 0.67]),
    ("emerald", [0.0, 1.0, 0.50]),
    ("junglegreen", [0.01, 1.0, 0.48]),
    ("seagreen", [0.31, 1.0, 0.50]),
    ("green", [0.0, 1.0, 0.0]),
    ("forestgreen", [0.0, 0.88, 0.0]),
    ("pinegreen", [0.0, 0.75, 0.16]),
    ("limegreen", [0.50, 1.0, 0.0]),
    ("yellowgreen", [0.56, 1.0, 0.26]),
    ("springgreen", [0.74, 1.0, 0.24]),
    ("olivegreen", [0.0, 0.60, 0.0]),
    ("rawsienna", [0.55, 0.0, 0.0]),
    ("sepia", [0.30, 0.0, 0.0]),
    ("brown", [0.40, 0.0, 0.0]),
    ("tan", [0.86, 0.58, 0.44]),
    ("white", [1.0, 1.0, 1.0]),
    ("black", [0.0, 0.0, 0.0]),
    ("gray", [0.5, 0.5, 0.5]),
];
/// C: `void setcolor(const char *name);`
pub fn setcolor(name: &str) {
    let name = name.to_ascii_lowercase();
    let mut c = [0.5, 0.5, 0.5];
    if let Some((_, v)) = COLORS.iter().find(|(n, _)| *n == name) {
        c = *v;
    } else if let Some(hex) = name.strip_prefix('#') {
        let n = hex.len();
        if [3, 6, 9, 12].contains(&n) {
            let w = n / 3;
            let max = ((1u32 << (w * 4)) - 1) as f32;
            if let (Ok(r), Ok(g), Ok(b)) = (
                u32::from_str_radix(&hex[..w], 16),
                u32::from_str_radix(&hex[w..2 * w], 16),
                u32::from_str_radix(&hex[2 * w..], 16),
            ) {
                c = [r as f32 / max, g as f32 / max, b as f32 / max];
            }
        }
    }
    setrgbcolor(c[0], c[1], c[2]);
}
fn image_impl(w: i32, h: i32, b: i32, m: &[f32], samples: &[u8], n: i32) {
    let row = 1 + (w * b * n - 1) / 8;
    let bytes = (row * h) as usize;
    printf!(
        "/picstr %d string def\n%d %d %d [%.4g %.4g %.4g %.4g %.4g %.4g]\n",
        row,
        w,
        h,
        b,
        m[0],
        m[1],
        m[2],
        m[3],
        m[4],
        m[5]
    );
    if n == 1 {
        printf!("{currentfile picstr readhexstring pop} image\n");
    } else {
        printf!(
            "{currentfile picstr readhexstring pop}\nfalse %d colorimage\n",
            n
        );
    }
    for chunk in samples[..bytes].chunks(39) {
        for v in chunk {
            printf!("%02x", *v as u32);
        }
        printf!("\n");
    }
}
/// C: `void image(int w, int h, int bps, float m[], unsigned char *samples);`
pub fn image(w: i32, h: i32, bps: i32, m: &[f32], samples: &[u8]) {
    image_impl(w, h, bps, m, samples, 1);
}
/// C: `void rgbimage(int w, int h, int bpc, float m[], unsigned char *samples);`
pub fn rgbimage(w: i32, h: i32, bpc: i32, m: &[f32], samples: &[u8]) {
    image_impl(w, h, bpc, m, samples, 3);
}
/// C: `void polyline(const float *x, const float *y, int n);`
pub fn polyline(x: &[f32], y: &[f32]) {
    assert_eq!(x.len(), y.len());
    gsave();
    newpath();
    moveto(x[0], y[0]);
    for i in 1..x.len() {
        lineto(x[i], y[i]);
        if i % 200 == 0 {
            stroke();
            moveto(x[i], y[i]);
        }
    }
    stroke();
    grestore();
}
/// C: `void markto(float x, float y, int index, float size);`
pub fn markto(x: f32, y: f32, index: i32, size: f32) {
    gsave();
    translate(x, y);
    scale(size, size);
    printf!("CLW %0.4g div SLW\n", size);
    newpath();
    match index % 9 {
        0 => {
            moveto(-0.5, 0.0);
            rlineto(1.0, 0.0);
            moveto(0.0, -0.5);
            rlineto(0.0, 1.0);
            stroke();
        }
        1 => {
            moveto(-0.5, 0.0);
            rlineto(1.0, 0.0);
            moveto(-0.25, -0.433);
            rlineto(0.5, 0.866);
            moveto(-0.25, 0.433);
            rlineto(0.5, -0.866);
            stroke();
        }
        2 => {
            moveto(-0.5, -0.5);
            rlineto(1.0, 1.0);
            moveto(-0.5, 0.5);
            rlineto(1.0, -1.0);
            stroke();
        }
        3 | 6 => {
            moveto(-0.5, -0.25);
            rlineto(1.0, 0.0);
            rlineto(-0.5, 0.809);
            closepath();
            if index % 9 == 3 {
                stroke()
            } else {
                fill()
            }
        }
        4 | 7 => {
            moveto(-0.5, -0.5);
            rlineto(1.0, 0.0);
            rlineto(0.0, 1.0);
            rlineto(-1.0, 0.0);
            closepath();
            if index % 9 == 4 {
                stroke()
            } else {
                fill()
            }
        }
        5 | 8 => {
            arc(0.0, 0.0, 0.5, 0.0, 360.0);
            if index % 9 == 5 {
                stroke()
            } else {
                fill()
            }
        }
        _ => {}
    }
    grestore();
}
fn rect(x: f32, y: f32, w: f32, h: f32) {
    newpath();
    moveto(x, y);
    rlineto(w, 0.0);
    rlineto(0.0, h);
    rlineto(-w, 0.0);
    closepath();
}
/// C: `void rectclip(float x, float y, float width, float height);`
pub fn rectclip(x: f32, y: f32, w: f32, h: f32) {
    rect(x, y, w, h);
    clip();
    newpath();
}
/// C: `void rectfill(float x, float y, float width, float height);`
pub fn rectfill(x: f32, y: f32, w: f32, h: f32) {
    gsave();
    rect(x, y, w, h);
    fill();
    grestore();
}
/// C: `void rectstroke(float x, float y, float width, float height);`
pub fn rectstroke(x: f32, y: f32, w: f32, h: f32) {
    gsave();
    rect(x, y, w, h);
    stroke();
    grestore();
}
