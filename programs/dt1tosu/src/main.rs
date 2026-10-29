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

//! DT1TOSU - Convert Sensors & Software DT1 radar data to SU.
//! Port of su/main/data_conversion/dt1tosu.c (Revision 1.17).
//! Credits: CWP John Stockwell; based on UTULSA sugpr by Chris Liner and Bill Underwood.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t",
    " DT1TOSU - Convert ground-penetrating radar data in the\t",
    "\tSensors & Software X.DT1 GPR format to SU format.\t",
    "\t\t\t\t\t\t\t\t",
    " dt1tosu < gpr_data_in_dt1_format  > stdout\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    " ns=from header\tnumber of samples per trace\t\t",
    " dt=.8\t\ttime sample interval (see below)\t\t",
    " swap=endian\tendian is auto-determined =1 (big endian) swap\t",
    "\t\t=0 don't swap bytes (little endian machines)\t",
    " verbose=0\tsilent\t\t\t\t\t\t",
    "\t\t=1 S & S header values from first trace\t\t",
    "\t\t\tsent to outpar\t\t\t\t",
    "\t\t=2 S & S header values from all traces\t\t",
    "\t\t\tsent to outpar\t\t\t\t",
    " outpar=/dev/tty\toutput parameter file\t\t\t",
    " list=0\tsilent\t\t\t\t\t\t",
    "\t\t=1 list explaining labels used in verbose\t",
    "\t\t     is printed to stderr\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Caution: An incorrect ns field will munge subsequent processing.",
    "\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t",
    " For compatiblity with SEGY header, apparent dt is set to\t",
    " .8 ms (800 microsecs).  Actual dt is .8 nanosecs.\t\t",
    " Using TRUE DISTANCES, this scales velocity\t\t\t",
    " and frequency by a factor of 1 million.\t\t\t",
    "\tExample: v_air = 9.83X10^8 ft/s\t (real)\t\t\t",
    "\t\t v_air = 983 ft/s\t(apparent for su)\t",
    "\tExample: fnyquist = 625 MHz\t(real)\t\t\t",
    "\t\tfnyquist = 625 Hz\t(apparent for su)\t",
    "\t\t\t\t\t\t\t\t",
    " IBM RS6000, NeXT, SUN are examples of big endian machines\t",
    " PC's and DEC are examples of little endian machines\t\t",
    "\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t",
    " This program has not been tested on DEC, some modification of the",
    " byte swapping routines may be required.\t\t\t",
    "\t\t\t\t\t\t\t\t",
];

fn value(b: &[u8], offset: usize) -> f32 {
    f32::from_ne_bytes(b[offset..offset + 4].try_into().unwrap())
}

fn swap_header(b: &mut [u8]) {
    for offset in [0, 4, 8, 12, 20, 24, 28, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92] {
        b[offset..offset + 4].reverse();
    }
    for offset in [32, 40, 48] { b[offset..offset + 8].reverse(); }
}

fn fprintf_sstr(fp: &CFile, b: &[u8]) {
    fprintf!(fp, "tracl = %d\n", value(b, 0) as i32);
    fprintf!(fp, "posit = %f\n", value(b, 4));
    fprintf!(fp, "ns = %d\n", value(b, 8) as i32);
    fprintf!(fp, "topog = %f\n", value(b, 12));
    fprintf!(fp, "bpp = %d\n", value(b, 20) as i32);
    fprintf!(fp, "twind = %f\n", value(b, 24));
    fprintf!(fp, "nstks = %d\n", value(b, 28) as i32);
    fprintf!(fp, "tza = %f\n", value(b, 80));
    fprintf!(fp, "zflag = %d\n", value(b, 84) as i32);
    fprintf!(fp, "nuc = %f\n", value(b, 88));
    let cflag = value(b, 96) as i32;
    fprintf!(fp, "cflag = %d\n", cflag);
    if cflag != 0 { efwrite(&b[100..128], fp); } else { fprintf!(fp, "no comm\n"); }
}

fn assign_header(tr: &mut Trace, b: &[u8], dt: f32) {
    tr.tracl = value(b, 0) as i32;
    tr.sx = value(b, 68) as i32;
    tr.sy = value(b, 72) as i32;
    tr.selev = value(b, 76) as i32;
    tr.gx = value(b, 56) as i32;
    tr.gy = value(b, 60) as i32;
    tr.gelev = value(b, 64) as i32;
    tr.ns = value(b, 8) as i32 as u16;
    tr.dt = ((dt as f64 + if dt >= 0.0 { 0.5 } else { -0.5 }) as i32) as u16;
    tr.delrt = value(b, 80) as i32 as i16;
    tr.trid = ((value(b, 84) as f64 + 1.0) as i32) as i16;
    let hour = value(b, 92) / 3600.0_f32;
    let minute = (hour - hour as i32 as f32) * 60.0_f32;
    tr.hour = hour as i32 as i16;
    tr.minute = minute as i32 as i16;
    tr.sec = ((minute - minute as i32 as f32) * 60.0_f32) as i32 as i16;
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    match filestat(STDOUT) {
        FileType::Tty => err!("stdout can't be tty"),
        FileType::Directory => err!("stdout must be a file, not a directory"),
        FileType::BadFileType => err!("stdout is illegal filetype"),
        _ => ()
    }
    let mut swap = if cfg!(target_endian = "big") { 1 } else { 0 };
    getparint("swap", &mut swap);
    // ssdt1 includes trailing padding for its double fields' alignment.
    let mut b = vec![0u8; (128 + 2 * SU_NFLTS + 7) & !7];
    efread(&mut b[..128], &stdin());
    if swap != 0 { swap_header(&mut b); }
    let mut ns = value(&b, 8) as i32;
    getparint("ns", &mut ns);
    if ns == 0 { err!("ns not set!"); }
    warn!("ns= %d", ns);
    let mut dt = 0.8f32;
    getparfloat("dt", &mut dt);
    dt = (dt as f64 * 1000.0) as f32;
    let mut list = 0;
    let mut verbose = 0;
    let mut outpar = String::from("/dev/tty");
    getparint("list", &mut list);
    getparint("verbose", &mut verbose);
    getparstring("outpar", &mut outpar);
    checkpars();
    let outparfp = efopen(&outpar, "w");
    if verbose == 1 { fprintf_sstr(&outparfp, &b); }
    if list != 0 { for line in LIST { fprintf!(stderr(), "%s\n", *line); } }
    let mut tr = Trace::new();
    assign_header(&mut tr, &b, dt);
    if ns > SU_NFLTS as i32 {
        // C would write beyond tr.data in the first sample loop.
        err!("ns exceeds SU_NFLTS: %d", ns);
    }
    for i in 0..ns.max(0) as usize {
        let mut temp = [0u8; 2];
        efread(&mut temp, &stdin());
        if swap != 0 { temp.reverse(); }
        tr.data[i] = i16::from_ne_bytes(temp) as f32;
    }
    puttr(&mut tr);
    // sizeof(short) promotes ns to size_t before the multiplication in C.
    // Thus ns=-64 requests zero bytes, while ns=-65 requests SIZE_MAX-1.
    let ssdatabytes = (ns as usize).wrapping_mul(2);
    let ssbytes = 128usize.wrapping_add(ssdatabytes);
    loop {
        // Negative ns can wrap ssbytes to a huge size_t. C still permits
        // a short read (including EOF) that fits in sstr. Read up to its
        // size, then reject only an actual write beyond the C allocation.
        let readbytes = ssbytes.min(b.len());
        let nread = efread(&mut b[..readbytes], &stdin());
        if ssbytes > b.len() && nread == b.len() {
            let mut extra = [0u8; 1];
            if efread(&mut extra, &stdin()) != 0 {
                err!("DT1 record exceeds ssdt1 buffer for ns=%d", ns);
            }
        }
        if nread == 0 { break; }
        if swap != 0 { swap_header(&mut b); }
        assign_header(&mut tr, &b, dt);
        for i in 0..ns.max(0) as usize {
            let mut temp = [b[128 + 2*i], b[129 + 2*i]];
            if swap != 0 { temp.reverse(); }
            tr.data[i] = i16::from_ne_bytes(temp) as f32;
        }
        if verbose == 2 { fprintf_sstr(&outparfp, &b); }
        puttr(&mut tr);
    }
    if verbose != 0 { efclose(&outparfp); }
    exit(CWP_Exit());
}

const LIST: &[&str] = &[
    " float\t tracl;\t\ttrace number\t\t\t\t",
    " float\t posit;\t\tposition\t\t\t\t",
    " float\t ns;\t\tnumber of points per trace\t\t",
    " float\t topog;\t\tTopographic data if any\t\t\t",
    " float\t nua;\t\tnot used \"a\"\t\t\t\t",
    " float\t bpp;\t\tbytes/point (always 2 for Rev 3 firmware)",
    " float\t twind;\t\ttime window\t\t\t\t",
    " float\t nstks;\t\tnumber of stacks\t\t\t",
    " double  gps_x;\t/* GPS X X-position */\t\t\t",
    " double  gps_y;\t/* GPS Y Y-position */ \t\t\t",
    " double  gps_z;\t/* GPS Z Z-position */\t\t\t",
    " float  rx;\t\t/* receiver x-position */\t\t",
    " float  ry;\t\t/* receiver y-position */\t\t",
    " float  rz;\t\t/* receiver z-position */\t\t",
    " float  tx;\t\t/* transmitter x-position */\t\t",
    " float  ty;\t\t/* transmitter y-position */ \t\t",
    " float  tz;\t\t/* transmitter z-position */ \t\t",
    " float\t tza;\t\ttime zero adjustment\t\t\t",
    "\t\t\twhere point(x) = point(x + adjustment)\t",
    " float\t zflag;\t\tzero flag: 0=data ok ,\t1 = zero data\t",
    " float\t nuc;\t\tnot used \"c\"\t\t\t\t",
    " float\t tod;\t\ttime of day data are collected\t\t",
    "\t\t\t\tin secs past midnight\t\t",
    " float\t cflag;\t\tcomment flag =1 if comment is attached\t",
    " char comm[28];\tcomment 28 characters\t\t\t",
    " short data[SU_NFLTS];\tdata part of the ssdt1 structure\t",
];
