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

// Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! SUAZIMUTH - compute trace azimuth, offset, and midpoint coordinates.
//! Port of su/main/headers/suazimuth.c (Revision 1.21).
//! Credits: based on suchw and su3dchart; CWP: John Stockwell;
//! UTulsa: Chris Liner; RISSC: Nils Maercklin.
//! setval credits: Einar Kajartansson, Jack K. Cohen, John Stockwell.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{nint, PI};
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                        ",
    " SUAZIMUTH - compute trace AZIMUTH, offset, and midpoint coordinates    ",
    "             and set user-specified header fields to these values       ",
    "                                                                        ",
    "  suazimuth <stdin >stdout [optional parameters]                        ",
    "                                                                        ",
    " Required parameters:                                                   ",
    "     none                                                               ",
    "                                                                        ",
    " Optional parameters:                                                   ",
    "   key=otrav      header field to store computed azimuths in            ",
    "   scale=1.0      value(key) = scale * azimuth                          ",
    "   az=0           azimuth convention flag                               ",
    "                   0: 0-179.999 deg, reciprocity assumed                ",
    "                   1: 0-359.999 deg, points from receiver to source     ",
    "                  -1: 0-359.999 deg, points from source to receiver     ",
    "   sector=1.0     if set, defines output in sectors of size             ",
    "                  sector=degrees_per_sector, the default mode is        ",
    "                  the full range of angles specified by az              ",
    "                                                                        ",
    "   offset=0       if offset=1 then set offset header field              ",
    "   signedflag=0   is offset signed? if =1 signedflag=1\t\t\t ",
    "   offkey=offset  header field to store computed offsets in             ",
    "   offkey=offset  header field to store computed offsets in             ",
    "                                                                        ",
    "   cmp=0          if cmp=1, then compute midpoint coordinates and       ",
    "                  set header fields for (cmpx, cmpy)                    ",
    "   mxkey=ep       header field to store computed cmpx in                ",
    "   mykey=cdp      header field to store computed cmpy in                ",
    "                                                                        ",
    " Notes:                                                                 ",
    "   All values are computed from the values in the coordinate fields     ",
    "   sx,sy (source) and gx,gy (receiver).                                 ",
    "   The output field \"otrav\" for the azimuth was chosen arbitrarily as ",
    "   an example of a little-used header field, however, the user may      ",
    "   choose any field that is convenient for his or her application.      ",
    "                                                                        ",
    "   Setting the sector=number_of_degrees_per_sector sets key field to    ",
    "   sector number rather than an angle in degrees.                       ",
    "                                                                        ",
    "   For az=0, azimuths are measured from the North, however, reciprocity ",
    "   is assumed, so azimuths go from 0 to 179.9999 degrees. If sector     ",
    "   option is set, then the range is from 0 to 180/sector.               ",
    "                                                                        ",
    "   For az=1, azimuths are measured from the North, with the assumption  ",
    "   that the direction vector points from the receiver to the source.    ",
    "   For az=-1, the direction vector points from the source to the        ",
    "   receiver. No reciprocity is assumed in these cases, so the angles go ",
    "   from 0 to 359.999 degrees.                                           ",
    "   If the sector option is set, then the range is from 0 to 360/sector. ",
    "                                                                        ",
    " Caveat:                                                                ",
    "   This program honors the value of scalco in scaling the values of     ",
    "   sx,sy,gx,gy. Type \"sukeyword scalco\" for more information.         ",
    "                                                                        ",
    "   Type \"sukeyword -o\" to see the keywords and descriptions of all    ",
    "   header fields.                                                       ",
    "                                                                        ",
    "   To plot midpoints, use: su3dchart                                    ",
    "                                                                        ",
];

fn compute_azimuth(sx: f64, sy: f64, mx: f64, my: f64, scale: f64, az: i32) -> f64 {
    let mut dval_out;
    if az == 0 {
        dval_out = 90.0 - (180.0 / PI) * ((sy - my) / (sx - mx)).atan();
        if dval_out == 180.0 { dval_out = 0.0; }
    } else if az < 0 {
        dval_out = (180.0 / PI) * (mx - sx).atan2(my - sy);
        if dval_out < 0.0 { dval_out += 360.0; }
    } else {
        dval_out = (180.0 / PI) * (sx - mx).atan2(sy - my);
        if dval_out < 0.0 { dval_out += 360.0; }
    }
    dval_out * scale
}

fn setval(type_out: &str, val_out: &mut Value, dval_out: f64) {
    match type_out.as_bytes()[0] {
        b's' => err!("can't change char header word"),
        b'h' => val_out.set_h(nint(dval_out) as i16),
        b'u' => val_out.set_u(nint(dval_out) as u16),
        b'l' => val_out.set_l(nint(dval_out) as i64),
        b'v' => val_out.set_v(nint(dval_out) as u64),
        b'i' => val_out.set_i(nint(dval_out)),
        b'p' => val_out.set_p(nint(dval_out) as u32),
        b'f' => val_out.set_f(dval_out as f32),
        b'd' => val_out.set_d(dval_out),
        _ => err!("unknown type %s", type_out),
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::new();
    let mut okey = String::new();
    let mut mxkey = String::new();
    let mut mykey = String::new();
    let mut scale = 0.0f64;
    let mut az = 0i32;
    let mut sector = 0.0f64;
    let mut oflag = 0i32;
    let mut signedflag = 0i32;
    let mut cmpflag = 0i32;
    if getparstring("key", &mut key) == 0 { key = "otrav".into(); }
    if getparstring("offkey", &mut okey) == 0 { okey = "offset".into(); }
    if getparstring("mxkey", &mut mxkey) == 0 { mxkey = "ep".into(); }
    if getparstring("mykey", &mut mykey) == 0 { mykey = "cdp".into(); }
    if getpardouble("scale", &mut scale) == 0 { scale = 1.0; }
    if getparint("az", &mut az) == 0 { az = 0; }
    if getpardouble("sector", &mut sector) == 0 { sector = 1.0; }
    if getparint("offset", &mut oflag) == 0 { oflag = 0; }
    if getparint("signedflag", &mut signedflag) == 0 { signedflag = 0; }
    if getparint("cmp", &mut cmpflag) == 0 { cmpflag = 0; }
    checkpars();
    scale /= sector;
    let type_out = hdtype(&key);
    let index_out = getindex(&key);
    let (otype_out, oindex_out) = if oflag != 0 { (hdtype(&okey), getindex(&okey)) } else { ("", 0) };
    let (mxtype_out, mxindex_out, mytype_out, myindex_out) = if cmpflag != 0 {
        (hdtype(&mxkey), getindex(&mxkey), hdtype(&mykey), getindex(&mykey))
    } else { ("", 0, "", 0) };
    let mut val_out = Value::default();
    while gettr(&mut tr) != 0 {
        let mut sx = tr.sx as f64;
        let mut sy = tr.sy as f64;
        let mut gx = tr.gx as f64;
        let mut gy = tr.gy as f64;
        let mut factor = if tr.scalco == 0 { 1.0 } else { tr.scalco as f64 };
        if factor < 0.0 { factor = -1.0 / factor; }
        sx *= factor; sy *= factor; gx *= factor; gy *= factor;
        let mx = 0.5 * (sx + gx);
        let my = 0.5 * (sy + gy);
        let azimuth = compute_azimuth(sx, sy, mx, my, scale, az);
        setval(type_out, &mut val_out, azimuth);
        puthval(&mut tr, index_out, &val_out);
        if oflag != 0 {
            let mut offsgn = 1;
            if signedflag != 0 {
                if gx - sx != 0.0 { offsgn = if (gx - sx) as f32 >= 0.0 { 1 } else { -1 }; }
                else if gy - sy != 0.0 { offsgn = if (gy - sy) as f32 >= 0.0 { 1 } else { -1 }; }
            }
            let mut offset = ((gx-sx)*(gx-sx) + (gy-sy)*(gy-sy)).sqrt();
            offset *= offsgn as f64;
            setval(otype_out, &mut val_out, offset);
            puthval(&mut tr, oindex_out, &val_out);
        }
        if cmpflag != 0 {
            let cmpx = (sx + gx) / 2.0 / factor;
            let cmpy = (sy + gy) / 2.0 / factor;
            setval(mxtype_out, &mut val_out, cmpx);
            puthval(&mut tr, mxindex_out, &val_out);
            setval(mytype_out, &mut val_out, cmpy);
            puthval(&mut tr, myindex_out, &val_out);
        }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
