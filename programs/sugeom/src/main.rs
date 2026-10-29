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

//! SUGEOM - fill geometry headers. Port of su/main/headers/sugeom.c (Revision 1.2).
//! Credits: Fernando Roxo <petro@roxo.org>.

use oxwave_core::cio::CFile;
use oxwave_core::prelude::*;
use oxwave_sps_records::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                          ",
    " SUGEOM - Fill up geometry in trace headers.                              ",
    "                                                                          ",
    " sugeom rps= sps= xps=  <stdin >stdout [optional parameters]              ",
    "                                                                          ",
    " Required parameters:                                                     ",
    "       rps=  filename for the receiver points definitions. (AKA R file)   ",
    "       sps=  filename for the source   points definitions. (AKA S file)   ",
    "       xps=  filename for the relation specification.      (AKA X file)   ",
    "                                                                          ",
    " OBS:  if 'prefix' (see bellow) is given the last three parameters become ",
    "       optional, but any one can be used as an override for the name      ",
    "       built using prefix.                                                ",
    "          Ex: sugeom prefix=blah  xps=bleh.xxx                            ",
    "              Used filenames:  blah.s blah.r and bleh.xxx                 ",
    "                                                                          ",
    " Optional parameters:                                                     ",
    "       prefix=  prefix name for the SPS files.  If given the filename will",
    "                be constructed as '<prefix>.s', '<prefix>.r' and          ",
    "                '<prefix>.x'.                                             ",
    "       rev=0             for SPS revision 0 format.                       ",
    "          =2 (Default)   for SPS revision 2.1 format.                     ",
    "   verbose=0 (Default)   silent run.                                      ",
    "          =1             verbose stats.                                   ",
    "     ibin=0  (default)   inline binsize (cdp interval)                    ",
    "                         0 -> do not compute cdp number                   ",
    "                                                                          ",
    "  For SPS format specification consult:                                   ",
    "    http://www.seg.org/resources/publications/misc/technical-standards    ",
    "                                                                          ",
    "                                                                          ",
    "  WARNING: This is not a full fledged geometry program.                   ",
    "                                                                          ",
    "  The SPS format is fully described in two documents available at the     ",
    "  above SEG address.                                                      ",
    "                                                                          ",
    "  To make short what can be a loooong discussion, the SPS file format is  ",
    "  intended to describe completely a land (or tranzition zone) survey. The ",
    "  main difference between land and marine survey is that in land all      ",
    "  points are meant to be previously located and remain fixed along        ",
    "  acquisition operation.  In a marine survey the cables and boat can be   ",
    "  dragged by current and wind, making the prosition of source and receivers",
    "  unpredictable.  The planning and survey description of a marine program ",
    "  must take into account all moving points, being necessary a full        ",
    "  positioning description of source and receivers at every single shot.   ",
    "                                                                          ",
    "  The SPS format standard is composed of three text files with 80 columns ",
    "  each.  Two of those files are Point Record Specification and are used to",
    "  describe all survey points, that is the receiver stations and source    ",
    "  points.  The remaining file is the Relation Record Specification and is ",
    "  used to describe how each source point is related to a set of           ",
    "  corresponding receiver points during the registration of each record    ",
    "  (fldr).                                                                 ",
    "                                                                          ",
    "  These files usually have a number of lines at the start describing the  ",
    "  survey and the projections used for coordinates.  This program just     ",
    "  skip them.  Those are header entries and have an 'H' in first column.   ",
    "                                                                          ",
    "  Each line (entry) of the  Point Record Specification contain the        ",
    "  information of a single point.  The only difference of a source point   ",
    "  specification entry and a receiver point specification entry is the     ",
    "  first  column of the file, it has an 'S' for a source point description ",
    "  and an 'R' for a receiver point description.  For each point entry there",
    "  is an identification with a pair of informations, the Line it belongs   ",
    "  and the Point Number, a set of coordinates (X, Y, surface Elevation,    ",
    "  depth of element), and the static correction.  All source points        ",
    "  description are in a single file (known as S-File), and all receiver    ",
    "  station informations are in another file (R-File).                      ",
    "                                                                          ",
    "  The Relation Record Specification (X-File) is a file with as many       ",
    "  entries (text lines) as necessary to completely describe each record    ",
    "  (fldr) acquired.   Each entry containing a record information starts    ",
    "  with an 'X' in first column.  The informations in the entry starts with ",
    "  the tape number in which the register was recorded (not used for this   ",
    "  program) and the Field Record Number (fldr).  Next comes the source     ",
    "  point description (the same Line and Point(ep) from the S-File). Then a ",
    "  sequence of channels (tracf) numbers (just first and last) and a channel",
    "  increment.  Next comes the receiver description using the same Line and ",
    "  Point/station from the R-File).  The Receiver Points are specified as the",
    "  first and last station used at that line.  Finally comes the recording  ",
    "  time.   If the spread has a gap, very common for older lines, it will   ",
    "  require at least two entries (text lines) to describe this record, one  ",
    "  describing the channels and receiver stations before the gap, and another",
    "  for the channels and stations after the gap.   The initial informations ",
    "  (tape, record, and source) are repeated for all entries.   In a 3D there",
    "  is one entry for each line of the patch.                                ",
    "                                                                          ",
    "  To use this program it is necessary to use all three SPS files.  The    ",
    "  S-File and the R-File can describe more points than will be used for the",
    "  processing.  For example, if one is processing just a section of a 2D   ",
    "  line the Point Record Specification files can describe all points of the",
    "  complete line, the information in excess will be disregarded.           ",
    "                                                                          ",
    "  The X-File must have the information in the same order of the input SU  ",
    "  file, if not so all files (fldr) that do not match the order will be    ",
    "  skipped until the program find the fldr corresponding to the sequence of",
    "  the X-File.                                                             ",
    "                                                                          ",
    "  Although this program read all of the SPS files, it is not ready for 3D ",
    "  geometry processing.   It is just a basic 2D straight line geometry     ",
    "  processing program that (hopefully) will fill correctly the informations",
    "  in header.  It can be used for a crooked line for the coordinates (X, Y ",
    "  Elevation, and element depth), the static correction, and offset.  The  ",
    "  cdp numbering will need an specialized program for computation.         ",
    "                                                                          ",
    "  All coordinates are expected to be an UTM projection in meters or feet. ",
    "                                                                          ",
    "  The X file order must be the same as the input file.  Upon reading      ",
    "  a trace whose trcf is not the next record to be processed in X file,    ",
    "  this record, and the following, will be skipped until a match is found  ",
    "  for the current entry in X file.   This is a way to, e.g., skip a noisy ",
    "  record, just remove its entry from X file.                              ",
    "                                                                          ",
    "  If a trace (understand fldr/tracf pair) is not represented in the X file",
    "  it will be skipped until a trace in current fldr matches the current    ",
    "  set of X entry for that fldr.  If it is desirable to process just a     ",
    "  subset of channels, just keep in the X-File information about those     ",
    "  channels.                                                               ",
    "                                                                          ",
    "  Updated trace header positions:                                         ",
    "                                                                          ",
    "  tracl - it keeps counting to overcome possible trace skipping.          ",
    "  ep    - If zeroed out in header it uses the relation (X) file info.     ",
    "                                                                          ",
    "  sx, sy, gx, gy -  will be updated with the coordinates in point files   ",
    "                    R and S. The scalco value used is the already         ",
    "                    stored in header, almost ever equal zero.             ",
    "  gelev, selev, sdepth - will be updated with values in point files       ",
    "                    R and S. The scalel value used is the already         ",
    "                    stored in header, almost ever equal zero.             ",
    "  sstat, gstatt - these values are filled with the information of the     ",
    "           Point Record Specification files.  If not available will be zero.",
    "  offset - is computed from source and receiver coordinates.  The offset  ",
    "           signal is made negative if source station number is greater    ",
    "           than the receiver station number.  It consider that source and ",
    "           receiver numbering are the same.   Fractionary stations, e.g.  ",
    "           source point halfway between stations are ok.                  ",
    "  cdp    - if parameter 'ibin' is passed it is considered the cdp spacing ",
    "           and the cdp number will be computed as follow:                 ",
    "           The cdp position is computed as the midpoint between source    ",
    "           and receiver.  The distance from this point to the first       ",
    "           station of R point file is divided by 'ibin' and added to 100. ",
    "           This 100 is absolutely artitrary.  ;)                          ",
    "                                                                          ",
];

#[derive(Default)]
struct Point {
    line: f32,
    point: f32,
    stat: i32,
    depth: f32,
    x: f64,
    y: f64,
    elev: f32,
}
#[derive(Default)]
struct Relation {
    num: i32,
    inc: i32,
    sl: f32,
    sp: f32,
    fc: i32,
    tc: i32,
    ic: i32,
    rl: f32,
    fr: f32,
    tr: f32,
    ri: i32,
}
fn count(fp: &CFile, typ: u8, verbose: i32) -> i32 {
    let (mut n, mut total) = (0, 0);
    while let Some(s) = fgets(120, fp) {
        total += 1;
        if s.as_bytes().first() == Some(&typ) {
            n += 1;
        }
    }
    if verbose != 0 {
        warn!("Read %d entries, %d of type %c\n", total, n, typ as i32);
    }
    erewind(fp);
    n
}
fn points(fp: &CFile, typ: u8, rev: i32) -> Vec<Point> {
    let mut v = Vec::new();
    while let Some(s) = fgets(120, fp) {
        let b = s.as_bytes();
        if b.first() != Some(&typ) {
            continue;
        }
        let mut p = Point::default();
        getSPSfield(b, rev, 1, SPSField::Float(&mut p.line));
        getSPSfield(b, rev, 2, SPSField::Float(&mut p.point));
        getSPSfield(b, rev, 5, SPSField::Int(&mut p.stat));
        getSPSfield(b, rev, 6, SPSField::Float(&mut p.depth));
        getSPSfield(b, rev, 10, SPSField::Double(&mut p.x));
        getSPSfield(b, rev, 11, SPSField::Double(&mut p.y));
        getSPSfield(b, rev, 12, SPSField::Float(&mut p.elev));
        v.push(p);
    }
    v
}
fn relations(fp: &CFile, rev: i32) -> Vec<Relation> {
    let mut v = Vec::new();
    while let Some(s) = fgets(120, fp) {
        let b = s.as_bytes();
        if b.first() != Some(&b'X') {
            continue;
        }
        let mut p = Relation {
            inc: 1,
            ic: 1,
            ri: 1,
            ..Default::default()
        };
        getSPSfield(b, rev, 2, SPSField::Int(&mut p.num));
        getSPSfield(b, rev, 3, SPSField::Int(&mut p.inc));
        getSPSfield(b, rev, 5, SPSField::Float(&mut p.sl));
        getSPSfield(b, rev, 6, SPSField::Float(&mut p.sp));
        getSPSfield(b, rev, 8, SPSField::Int(&mut p.fc));
        getSPSfield(b, rev, 9, SPSField::Int(&mut p.tc));
        getSPSfield(b, rev, 10, SPSField::Int(&mut p.ic));
        getSPSfield(b, rev, 11, SPSField::Float(&mut p.rl));
        getSPSfield(b, rev, 12, SPSField::Float(&mut p.fr));
        getSPSfield(b, rev, 13, SPSField::Float(&mut p.tr));
        getSPSfield(b, rev, 14, SPSField::Int(&mut p.ri));
        v.push(p);
    }
    v
}
fn index(v: &[Point], point: f32, line: f32, kind: &str) -> usize {
    v.iter()
        .position(|p| {
            ((p.line - line) as f64).abs() <= 0.01 && ((p.point - point) as f64).abs() <= 0.01
        })
        // The C lookup returns an invalid index here; fail cleanly instead of
        // accessing beyond the point array or panicking.
        .unwrap_or_else(|| err!("**** %s point not found: line=%g, point=%g", kind, line, point))
}
fn filename(key: &str, prefix: &Option<String>, suffix: &str) -> String {
    let mut s = String::new();
    if getparstring(key, &mut s) != 0 {
        return s;
    }
    if let Some(p) = prefix {
        return format!("{p}.{suffix}");
    }
    err!("**** Missing required parameter '%s=', see self-doc", key);
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let (mut verbose, mut rev, mut ibin) = (0, 2, 0f32);
    getparint("verbose", &mut verbose);
    getparint("rev", &mut rev);
    getparfloat("ibin", &mut ibin);
    let mut pref = String::new();
    let prefix = if getparstring("prefix", &mut pref) != 0 {
        Some(pref)
    } else {
        None
    };
    let name = filename("rps", &prefix, "r");
    if verbose != 0 {
        warn!("Receiver file: %s", name);
    }
    let fp = efopen(&name, "r");
    let nr = count(&fp, b'R', verbose);
    warn!("Found %d receiver points entries.", nr);
    if nr == 0 {
        err!("**** No receiver point information found!!");
    }
    let rec = points(&fp, b'R', rev);
    if verbose != 0 {
        warn!("Stored: %d Receiver Point Informations.", rec.len() as i32);
        for p in &rec {
            warn!(
                "R %10.2f,%10.2f,%4.1f,%10.2lf,%10.2lf,%5.1f",
                p.line, p.point, p.depth, p.x, p.y, p.elev
            );
        }
    }
    let name = filename("sps", &prefix, "s");
    if verbose != 0 {
        warn!("Source file: %s", name);
    }
    let fp = efopen(&name, "r");
    let ns = count(&fp, b'S', verbose);
    warn!("Found %d source points entries.", ns);
    if ns == 0 {
        err!("**** No source point information found!!");
    }
    let src = points(&fp, b'S', rev);
    if verbose != 0 {
        warn!("Stored: %d Source Point Informations.", src.len() as i32);
        for p in &src {
            warn!(
                "S %10.2f,%10.2f,%4.1f,%10.2lf,%10.2lf,%5.1f",
                p.line, p.point, p.depth, p.x, p.y, p.elev
            );
        }
    }
    let name = filename("xps", &prefix, "x");
    if verbose != 0 {
        warn!("Relation file: %s", name);
    }
    let fp = efopen(&name, "r");
    let nx = count(&fp, b'X', verbose);
    warn!("Found %d relation information entries.", nx);
    if nx == 0 {
        err!("**** No relation information found!!");
    }
    let files = relations(&fp, rev);
    if verbose != 0 {
        warn!(
            "Stored: %d Register Information entries.",
            files.len() as i32
        );
        for p in &files {
            warn!(
                "X %4d,%1d,%9.1f,%9.1f,%4d,%4d,%1d,%9.1f,%9.1f,%10.1f,%1d",
                p.num, p.inc, p.sl, p.sp, p.fc, p.tc, p.ic, p.rl, p.fr, p.tr, p.ri
            );
        }
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    if verbose != 0 {
        if (0..=17).contains(&tr.trid) || tr.trid == 109 {
            warn!("input is seismic data, trid=%d", tr.trid);
        } else {
            warn!("input is not seismic data, trid=%d", tr.trid);
        }
    }
    let mut dt = 0f32;
    if getparfloat("dt", &mut dt) == 0 {
        warn!("WARNING:  neither dt nor d1 are set, nor is dt getparred!");
        warn!("WARNING:  It is not used, just checked!!");
    }
    let (mut oldfldr, mut oldskip, mut curx, mut curtrl, mut nproct) =
        (0i32, 0i32, -1i32, 0i32, 0i32);
    'outer: loop {
        'trace: {
            if tr.fldr != oldfldr {
                if curx + 1 < nx && tr.fldr == files[(curx + 1) as usize].num {
                    curx += 1;
                    if curx == nx {
                        break 'outer;
                    }
                    if verbose != 0 {
                        warn!("** processing fldr=%d", tr.fldr);
                    }
                    oldfldr = tr.fldr;
                } else {
                    if tr.fldr != oldskip {
                        if curx + 1 == nx {
                            curx += 1;
                            break 'outer;
                        }
                        warn!(
                            "** skipping fldr=%d, looking for %d",
                            tr.fldr,
                            files[(curx + 1) as usize].num
                        );
                        oldskip = tr.fldr;
                    }
                    break 'trace;
                }
            }
            if curx < 0 {
                break 'trace;
            }
            if tr.tracf < files[curx as usize].fc {
                warn!(
                    "** skipping fldr=%d, tracf=%d, looking for tracf=%d ",
                    tr.fldr, tr.tracf, files[curx as usize].fc
                );
                break 'trace;
            } else if tr.tracf > files[curx as usize].tc {
                if curx + 1 == nx {
                    curx += 1;
                    break 'outer;
                }
                if tr.fldr == files[(curx + 1) as usize].num {
                    loop {
                        curx += 1;
                        let f = &files[curx as usize];
                        if tr.tracf >= f.fc && tr.tracf <= f.tc {
                            break;
                        }
                        if curx + 1 >= nx || tr.fldr != files[(curx + 1) as usize].num {
                            break;
                        }
                    }
                    let f = &files[curx as usize];
                    if tr.tracf < f.fc || tr.tracf > f.tc {
                        break 'trace;
                    }
                }
            }
            let f = &files[curx as usize];
            curtrl += 1;
            tr.tracl = curtrl;
            if tr.ep == 0 {
                tr.ep = f.sp as i32;
            }
            tr.cdp = 0;
            // ARM C integer division returns zero for a zero divisor.
            let channel = (tr.tracf - f.fc).checked_div(f.ic).unwrap_or(0);
            let station = f.fr + (f.tr - f.fr) * channel as f32 / (f.tc - f.fc) as f32;
            let r = &rec[index(&rec, station, f.rl, "Receiver")];
            let s = &src[index(&src, f.sp, f.sl, "Source")];
            if tr.scalel == 0 {
                tr.scalel = 1
            }
            if tr.scalco == 0 {
                tr.scalco = 1
            }
            let mut el = (1.0 / tr.scalel as f64) as f32;
            let mut co = (1.0 / tr.scalco as f64) as f32;
            if tr.scalel < 0 {
                el = -(tr.scalel as i32) as f32
            }
            if tr.scalco < 0 {
                co = -(tr.scalco as i32) as f32
            }
            tr.sx = (s.x * co as f64) as i32;
            tr.sy = (s.y * co as f64) as i32;
            tr.gx = (r.x * co as f64) as i32;
            tr.gy = (r.y * co as f64) as i32;
            if ibin != 0.0 {
                let x = (r.x + s.x) / 2.0;
                let y = (r.y + s.y) / 2.0;
                let temp = ((x - rec[0].x).powf(2.0) + (y - rec[0].y).powf(2.0)).sqrt() as f32;
                tr.cdp = (100.5 + (temp / ibin) as f64) as i32;
            }
            tr.gelev = (r.elev * el) as i32;
            tr.selev = (s.elev * el) as i32;
            tr.sdepth = (s.depth * el) as i32;
            tr.offset = ((s.x - r.x).powf(2.0) + (s.y - r.y).powf(2.0)).sqrt() as i32;
            if s.point > r.point {
                tr.offset = -tr.offset;
            }
            tr.counit = 1;
            tr.sstat = s.stat as i16;
            tr.gstat = r.stat as i16;
            puttr(&mut tr);
            nproct += 1;
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    if curx >= nx {
        warn!("***  Last record processed -> fldr %d  ***", oldfldr);
        warn!("***  Records after this were discarded ***");
    }
    warn!("***  %d traces processed...", nproct);
    exit(CWP_Exit());
}
