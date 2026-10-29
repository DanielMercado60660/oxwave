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

//! SUBINCSV - grid binning; port of su/main/headers/subincsv.c, revision 1.1.
//! Credits: Andre Latour.

use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_header_case::{fromhead, GetCase};
use oxwave_k_file::{readkfile, writekfile};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                          ",
    " SUBINCSV - 3D Grid Resolve and output and assign traces to cdps (cells). ",
    "                                                                          ",
    " subincsv  rfile=in.csv wfile=out.csv <in.su >out.su                      ",
    "                                                                          ",
    " Parameter overview:                                                      ",
    "                                                                          ",
    "       rfile= if specified, read a K-file containing parameter values.    ",
    "       wfile= if specified, write a K-file with parameter values.         ",
    "                                                                          ",
    "       in.su and out.su do not have to be specified (if you just want     ",
    "       to specify command line parameters and get the output K-file).     ",
    "                                                                          ",
    "       Parameters specified on the command line override parameters from  ",
    "       the rfile.                                                         ",
    "                                                                          ",
    "       bintype=        See below for which other parameters are required  ",
    "                       for the option numbers here.                       ",
    "               30      Compute trace midpoint coordinates (sx+gx)/2 and   ",
    "                       (gx+gy)/2 and update cdp,igi,igc keys (grid cell   ",
    "                       number and inline and crossline index numbers).    ",
    "                       Option 30 is usually used for prestack traces.     ",
    "              -30      Use input cdp (cell) number and update keys        ",
    "                       igi,igc,sx,sy,gx,gy. Where igi,igc are grid index  ",
    "                       numbers of the cdp (cell) number and sx,sy is the  ",
    "                       cell centre in raw coordinates and gx,gy is cell   ",
    "                       centre in grid coordinates (which are shifted and  ",
    "                       aligned with grid definition, but not scaled by    ",
    "                       cell widths). Note that sx,sy are only approximate ",
    "                       cell centres since scalco rounds after sin,cosine  ",
    "                       computations, but gx,gy are usually much more      ",
    "                       precise since they are multiples of cell widths.   ",
    "                       Option -30 is usually used for poststack traces.   ",
    "              -31      Use input cdp (cell) number and update igi,igc.    ",
    "              -32      Use input igi,igc and update cdp number as well as ",
    "                       sx,sy,gx,gy as described in option -30.            ",
    "               20      Use point numbers to compute a 2d cdp number.      ",
    "                       Option 20 is usually used for prestack traces.     ",
    "                                                                          ",
    "                                                                          ",
    "       offset=         By default, bintype=30 and 20 recompute offset,    ",
    "                       but other bintypes leave it as-is.                 ",
    "             =1        Recompute offset key.                              ",
    "             =0        Do not recompute offset key.                       ",
    "                                                                          ",
    "       check=0         Do not print checking details.                     ",
    "             1         For grid bintypes, after grid defintion is set,    ",
    "                       run some grid functions on the 4 corner points     ",
    "                       and print the results. The intention here is to    ",
    "                       exercise many functions in case of issues created  ",
    "                       by coding or compiler or optimizer errors/changes. ",
    "                       But the output print may also be useful for users. ",
    "                       For instance, you will see slight differences in   ",
    "                       the coordinates of those 4 corners when produced   ",
    "                       by different functions, and when run on different  ",
    "                       hardware or with different compilers/optimizers.   ",
    "                                                                          ",
    " Grid parameters (required, either on command line or in rfile).          ",
    "                                                                          ",
    "    grid_xa=  X coordinate of corner A.                                   ",
    "    grid_ya=  Y coordinate of corner A.                                   ",
    "    grid_xb=  X coordinate for corner B.                                  ",
    "    grid_yb=  Y coordinate for corner B.                                  ",
    "    grid_xc=  X coordinate for corner C.                                  ",
    "    grid_yc=  Y coordinate for corner C.                                  ",
    "    grid_wb=  width of cells in A-->B direction.                          ",
    "    grid_wc=  width of cells in A-->C direction.                          ",
    "                                                                          ",
    " Note that corner A coordinates are used exactly, but corner B is reset   ",
    " to an exact multiple distance of the cell width in A-->B direction.      ",
    " And corner C is reset onto a line at right angle to A-->B direction      ",
    " through A and also to an exact multiple distance of A-->C cell width.    ",
    "                                                                          ",
    " -----------------------------------------------------------------------  ",
    "                                                                          ",
    " Cdps from Points parameters (either on command line or in rfile).        ",
    "                                                                          ",
    "    point_crz= cdp number of receiver zero                                ",
    "    point_cru= cdps per one receiver point unit                           ",
    "    point_csz= cdp number of source zero                                  ",
    "    point_csu= cdps per one source point unit                             ",
    "                                                                          ",
    " Cdps from Points parameters (only on command line).                      ",
    "                                                                          ",
    "    rpkey=     key containing receiver point numbers (default is gaps)    ",
    "    spkey=     key containing source point numbers (default is grnlof).   ",
    "               Note these defaults are same as sutoolcsv and sugeomcsv.   ",
    "                                                                          ",
    " The computed cdp value is simply:                                        ",
    "   cdp = (rpkey*point_cru + point_crz + spkey*point_csu + point_csz) / 2. ",
    " with value rounded to nearest integer.                                   ",
    "                                                                          ",
    " Note: The survey does not need an actual receiver or source point zero.  ",
    "                                                                          ",
    " Example: Assuming a 2D situation where you want 2 cdps per receiver.     ",
    "          If receiver points are numbered from 100 to 1000 by 1, then     ",
    "          point_crz=0 and point_cru=2 means receiver point 100 will       ",
    "          have a cdp number of 200 and receiver 101 will have cdp 202.    ",
    "          But if receivers are numbered 100 to 10000 by 10, then you can  ",
    "          specify point_cru=0.2 and receiver point 100 will have cdp 20   ",
    "          and receiver point 110 will have cdp 22. If you want, you can   ",
    "          set point_crz=180 to get the same cdp numbers as the first      ",
    "          situation. Both point_cru and point_csu can be negative. This   ",
    "          means your cdps will number backwards relative to the point     ",
    "          numbering, which is sometimes desirable/required. But you       ",
    "          should set point_crz and point_csz to avoid negative cdps.      ",
    "          Exactly the same computations are used with point_csz and       ",
    "          point_csu to produce the cdp number of the sources (shots).     ",
    "          Often for 2d surveys, the same point numbering values have been ",
    "          used for rpkey and spkey, so the parameters should be the same. ",
    "          An exception to this is when the sources are skidded to halfway ",
    "          between receivers. Even though the source is given the same     ",
    "          point number, you should adjust point_csz slightly to get the   ",
    "          most accurate cdp number at sources. Sometimes, of course, the  ",
    "          source point numbering will be completely unrelated to the      ",
    "          receiver point numbering, so you will have to specify different ",
    "          parameter values to get corresponding cdp numbering.            ",
    "          And there are situations where you want to produce 4 cdps per   ",
    "          receiver and source (or 3 or 1 or whatever). That just requires ",
    "          setting point_cru and point_csu properly. In any case, once     ",
    "          computed the trace receiver cdp and the trace source cdp are    ",
    "          summed and divided by 2 for the output trace cdp key value.     ",
    "                                                                          ",
    " Advice:  OK, why not just use the grid options with a very wide cell     ",
    "          width in the A-->C direction and put corner C on corner A?      ",
    "          That creates one cell in the A-->C direction and will function. ",
    "          If your input coordinates are just simple increments along a    ",
    "          straight line, then the wide-cell will make cdps having a nice  ",
    "          fold. But with actual surveyed coordinates, the cdp fold is     ",
    "          likely to vary a lot (with some cdps having roughtly twice the  ",
    "          average fold and some cdps having low fold (and not right next  ",
    "          to each other either). This occurs because surveyed coordinates ",
    "          tend to produce trace midpoints in clumps, and those clumps are ",
    "          not exactly the same interval apart. Slowly the cell boundaries ",
    "          in the A-->B direction get misaligned with those clumps. So,    ",
    "          at some location, two clumps get into one cell, and at other    ",
    "          locations, no clump gets into a cell. (Even when the cdp fold   ",
    "          is well-behaved, the offset distribution for wide-cell cdps     ",
    "          tends to vary a lot).                                           ",
    "          On the other hand, points tend to be numbered exactly as the    ",
    "          survey was designed. So using them means you get cdps that have ",
    "          the designed fold and offset distribution. But it also means    ",
    "          that cmps may have traces which spatially overlap other cmps.   ",
    "          In other words, wide-cell grids make spatially exact boundaries ",
    "          at the expense of nice fold and offset distribution, whereas    ",
    "          points can make nice fold and offset distribution at the        ",
    "          expense of exact spatial boundaries.                            ",
    "                                                                          ",
    " ***********************************************************              ",
    "   To output this documentation:  subincsv 2> bindoc.txt                  ",
    " ***********************************************************              ",
    "                                                                          ",
    " Seismic Unix has a 240 byte header which already has defined key names.  ",
    " In other seismic processing systems the ability to expand trace headers  ",
    " and insert intricate grid transform values is both a blessing and a curse",
    " (trust me on that). For SU, intricate grid-related values have no keys   ",
    " to be stored. So only 3 grid values are updated (to keys cdp, igi, igc). ",
    " Where igi is set to the cell index in the direction from corner A to     ",
    " corner B (a mnemonic for igi is index-grid-inline). And where igc is set ",
    " to cell index in the direction from corner A to corner C (a mnemonic for ",
    " igc is index-grid-crossline).                                            ",
    "                                                                          ",
    " The input grid definition command line parameters are processed and      ",
    " written to an external file. That file follows conventions established   ",
    " by SUTOOLCSV and SUGEOMCSV. I call this the K-file (K for Konstants).    ",
    " Reading this K-file should allow other SU programs to perform intricate  ",
    " grid transforms on-the-fly using sx,sy,gx,gy coordinates as well as      ",
    " backwards transforms from cdp,igi,igc to cell centre XYs. The grid in    ",
    " the K-file will also allow transforms of XYs values in S and R tables in ",
    " spreadsheets (eventually). And, once output, the K-file can be re-input  ",
    " to this program instead of using command line parameters.                ",
    "                                                                          ",
    " The grid will be an exact rectangle with 4 corner points A,B,C, and D.   ",
    " But you can only specify XYs for corners A,B,C (D is computed herein).   ",
    " Corner A coordinates are used exactly as input. Then the direction from  ",
    " corner A to input corner B is determined exactly. After that, corner B   ",
    " coordinates are adjusted to an exact multiple of the cell width you      ",
    " specify for the A-->B direction. Then your input coordinates of corner C ",
    " are used to compute the distance from corner A to corner C. The right    ",
    " angle to A-->B gives direction for output corner C (along line thru A).  ",
    " Corner C is then adjusted to an exact multiple of the cell width you     ",
    " specify for the A-->C direction. Note this means input corner C only     ",
    " determines how wide the rectangle is, and which side of A-->B the        ",
    " processed corner C is on. Corner D is computed from the other corners.   ",
    "                                                                          ",
    " The first cell is centred on corner A and has igi=1 and igc=1.           ",
    " Cell centres then increment by their corresponding widths in the A-->B   ",
    " and A-->C directions. The igi and igc values increment by 1. cdp numbers ",
    " start at 1 in the first cell and increment by 1 in the A-->B direction   ",
    " until reaching corner B, then moves 1 cell in the A-->C direction        ",
    " (near corner A) and continues to increment by 1 in the A-->B direction.  ",
    "                                                                          ",
    " Cells only contain one-half of their boundaries. This ensures that a     ",
    " trace midpoint that is exactly between 2 or 4 cell centres is assigned   ",
    " to a specific cell. Note: This is why you cannot use proximity to cell   ",
    " centres to assign traces to cells. You need to actually compute the cdp  ",
    " and igi,igc numbers the way that it is done herein.                      ",
    "                                                                          ",
    " For bintypes that use the grid, grid_xa,grid_ya,grid_xb,grid_yb,         ",
    " grid_xc,grid_yc,grid_wb,grid_wc must be specified either on the          ",
    " command line or in the input K-file. These values are the XYs of         ",
    " corners A,B,C as well as the cell widths. The output K-file contains     ",
    " the following processed grid definition values:                          ",
    "                                                                          ",
    "              Name      Definition                                        ",
    "              ----      ----------                                        ",
    "                                                                          ",
    "             bintype = bintype number                                     ",
    "             grid_lf = which side of A-->B is C on? 1=left, -1=right      ",
    "             grid_xa = raw, real world X coordinate of corner A           ",
    "             grid_ya = raw, real world Y coordinate of corner A           ",
    "             grid_xb = raw, real world X coordinate of corner B           ",
    "             grid_yb = raw, real world Y coordinate of corner B           ",
    "             grid_xc = raw, real world X coordinate of corner C           ",
    "             grid_yc = raw, real world Y coordinate of corner C           ",
    "             grid_xd = raw, real world X coordinate of corner D           ",
    "             grid_yd = raw, real world Y coordinate of corner D           ",
    "             grid_wb = width of cells in A-->B direction                  ",
    "             grid_wc = width of cells in A-->C direction                  ",
    "             grid_nb = number of cells in A-->B direction                 ",
    "             grid_nc = number of cells in A-->C direction                 ",
    "             grid_fp = first cdp (cell) number                            ",
    "             grid_lp = last  cdp (cell) number                            ",
    "             grid_sb = sine   of A-->B to X-axis.                         ",
    "             grid_cb = cosine of A-->B to X-axis.                         ",
    "                                                                          ",
    "  Note that corners A,B,C,D are at the centres of cells. Thus coordinates ",
    "  can be half a cell width outside A,B,C,D and still be inside the grid.  ",
    "                                                                          ",
    " Warning and advice:                                                      ",
    "  Cell boundaries and other grid computations use double precision values ",
    "  and are therefore extremely precise. This very precision causes issues. ",
    "  When a trace midpoint is very near a cell boundary, it only takes a     ",
    "  slight difference in hardware/compiler/optimizer computations for the   ",
    "  boundaries to move a bit, and therefore assign some traces to different ",
    "  cells. You should expect that. Similarly, trying to reverse or invert a ",
    "  grid by exchanging corners A and B and so on, is also not likely to     ",
    "  result in exactly the same distribution of traces in the cells.         ",
    "                                                                          ",
    "  -------------------------------------------------------------------     ",
    "                                                                          ",
    "  For bintype=20, point_crz,point_cru,point_csz,point_csu must be         ",
    "  specified either on the command line or in the input K-file.            ",
    "  The output K-file contains the same values:                             ",
    "                                                                          ",
    "              Name      Definition                                        ",
    "              ----      ----------                                        ",
    "                                                                          ",
    "             bintype   = bintype number                                   ",
    "             point_crz = cdp number of receiver zero                      ",
    "             point_cru = cdps per one receiver point unit                 ",
    "             point_csz = cdp number of source zero                        ",
    "             point_csu = cdps per one source point unit                   ",
    "                                                                          ",
    " Note that K-files can contain these values in any order, and             ",
    " K-files can contain other values that this program does not use.         ",
    "                                                                          ",
    " ----------------------------------------------------------------------   ",
    "                                                                          ",
];
fn main() {
    let a = args();
    initargs(&a);
    requestdoc(1, SDOC);
    let (mut r, mut w) = (String::new(), String::new());
    let hr = getparstring("rfile", &mut r) != 0;
    let hw = getparstring("wfile", &mut w) != 0;
    if hr && hw && r == w {
        err!("**** Error: wfile= output K-file must be different than rfile= input K-file.");
    }
    let traces = !isatty(0);
    if !traces {
        if !hw {
            err!("**** Error: wfile= output K-file name must be specified when no input traces.");
        }
        if !isatty(1) {
            err!("**** Error: Cannot specify output trace file with no input trace file.");
        }
    } else if isatty(1) {
        err!("**** Error: Must have output trace file when input trace file is specified.");
    }
    let (mut bt, mut off, mut check) = (-1, -1, 0);
    getparint("bintype", &mut bt);
    getparint("offset", &mut off);
    getparint("check", &mut check);
    let mut names = vec![String::new(); 999];
    let mut forms = names.clone();
    let mut vals = vec![0.; 999];
    let (mut count, mut ew) = (0, 0);
    if hr {
        let f = fopen(&r, "r")
            .unwrap_or_else(|| err!("rfile error: input K-file did not open correctly."));
        readkfile(&f, &mut names, &mut forms, &mut vals, &mut count, &mut ew);
        match ew {
            1 => err!("K-file read error: more than one C_SU_NAMES parameter record."),
            2 => err!("K-file read error: no C_SU_NAMES parameter record."),
            3 => err!("K-file read error: more than one C_SU_FORMS parameter record."),
            4 => err!("K-file read error: no C_SU_FORMS parameter record."),
            5 => err!("K-file read error: more than one C_SU_SETID record."),
            6 => err!("K-file read error: no C_SU_SETID record."),
            7 => {
                err!("K-file read error: different number of values on C_SU_NAMES and C_SU_FORMS.")
            }
            8 => err!("K-file read error: unable to allocate memory."),
            9 => err!("K-file read error: name exists at least twice in C_SU_NAMES list."),
            10 => err!("K-file read error: at least 1 field-unreadable as a number."),
            11 => err!("K-file read error: at least 1 field containing 2 numbers."),
            12 => {
                err!("K-file read error: not-enough-commas to get all values for C_SU_NAMES list.")
            }
            -1 => warn!("K-file read warning: at least 1 all-blank field, assumed zero for all."),
            n if n > 0 => err!("K-file read error: returned with some unrecognized error code."),
            _ => (),
        }
        for x in &mut names[..count as usize] {
            x.make_ascii_lowercase();
        }
    }
    let mut g = vec![-1.1e308; 18];
    if bt == -1 {
        g[0] = -1.;
        for j in 0..count as usize {
            if names[j] == "bintype" {
                g[0] = vals[j];
            }
        }
        bt = if g[0] > 0. {
            (g[0] + 0.1) as i32
        } else {
            (g[0] - 0.1) as i32
        };
    } else {
        g[0] = bt as f64;
    }
    if ![30, -30, -31, -32, 20].contains(&bt) {
        err!("**** Error: bintype number %d is not recognized.", bt);
    }
    if off == -1 {
        off = if bt == 30 || bt == 20 { 1 } else { 0 };
    }
    let gn: Vec<&str>;
    let (mut rp, mut sp) = (0, 0);
    if bt != 20 {
        gn = vec![
            "bintype", "grid_lf", "grid_xa", "grid_ya", "grid_xb", "grid_yb", "grid_xc", "grid_yc",
            "grid_xd", "grid_yd", "grid_wb", "grid_wc", "grid_nb", "grid_nc", "grid_fp", "grid_lp",
            "grid_sb", "grid_cb",
        ];
    } else {
        let (mut rk, mut sk) = (String::new(), String::new());
        if getparstring("rpkey", &mut rk) == 0 {
            rk = "gaps".into();
        }
        rp = GetCase(&rk);
        if rp < 1 {
            err!("Error: rpkey= %s is not recognized.", rk);
        }
        if getparstring("spkey", &mut sk) == 0 {
            sk = "grnlof".into();
        }
        sp = GetCase(&sk);
        if sp < 1 {
            err!("Error: spkey= %s is not recognized.", sk);
        }
        gn = vec![
            "bintype",
            "point_crz",
            "point_cru",
            "point_csz",
            "point_csu",
        ];
    }
    let range = if bt == 20 { 1..5 } else { 2..12 };
    for i in range {
        if bt != 20 && (i == 8 || i == 9) {
            continue;
        }
        if getpardouble(gn[i], &mut g[i]) == 0 {
            for j in 0..count as usize {
                if names[j] == gn[i] {
                    g[i] = vals[j];
                }
            }
            if g[i] < -1.0e308 {
                if bt == 20 {
                    err!(
                        "**** Error: bintype=%d and parameter %s not found.",
                        bt,
                        gn[i]
                    );
                } else {
                    err!("Error: bintype=%d and parameter %s not found.", bt, gn[i]);
                }
            }
        }
    }
    if bt != 20 {
        gridset(&mut g, &mut ew);
        match ew {
            1 => err!("gridset error: grid_wb cell width must be positive."),
            2 => err!("gridset error: grid_wc cell width must be positive."),
            3 => err!("gridset error: corner B is within grid_wb cell width of corner A."),
            -1 => warn!("gridset warning: corner C is near A and is reset to A."),
            n if n > 0 => err!("gridset error: returned with some unrecognized error code."),
            _ => (),
        }
        gridcheck(&g, check, &mut ew);
        if ew > 0 {
            err!("gridcheck error: returned with some unrecognized error code.");
        }
    }
    if hw {
        let f = fopen(&w, "w")
            .unwrap_or_else(|| err!("wfile error: output K-file did not open correctly."));
        let mut n = count as usize;
        for i in 0..gn.len() {
            let mut found = false;
            for j in 0..count as usize {
                if names[j] == gn[i] {
                    vals[j] = g[i];
                    found = true;
                    break;
                }
            }
            if !found {
                vals[n] = g[i];
                names[n] = gn[i].into();
                forms[n] = "%.20g".into();
                n += 1;
            }
        }
        writekfile(&f, &names[..n], &forms[..n], &vals[..n], &mut ew);
        if ew > 0 {
            err!("K-file write error: returned with an unrecognized error code.");
        }
        efclose(&f);
    }
    if !traces {
        exit(0);
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Error: cannot get first trace");
    }
    let mut n = 0;
    loop {
        let (mut dx, mut dy, mut tx, mut ty) = (0., 0., 0., 0.);
        let (mut c, mut i, mut j) = (0, 0, 0);
        if off == 1 {
            dx = tr.sx as f64;
            dy = tr.sy as f64;
            tx = tr.gx as f64;
            ty = tr.gy as f64;
            if tr.scalco > 1 {
                let s = tr.scalco as f64;
                dx *= s;
                dy *= s;
                tx *= s;
                ty *= s;
            } else if tr.scalco < 0 {
                let s = -(tr.scalco as f64);
                dx /= s;
                dy /= s;
                tx /= s;
                ty /= s;
            }
            tr.offset = ((dx - tx) * (dx - tx) + (dy - ty) * (dy - ty))
                .sqrt()
                .round_ties_even() as i32;
        }
        if bt == 30 {
            dx = 0.5 * tr.sx.wrapping_add(tr.gx) as f64;
            dy = 0.5 * tr.sy.wrapping_add(tr.gy) as f64;
            if tr.scalco > 1 {
                dx *= tr.scalco as f64;
                dy *= tr.scalco as f64;
            } else if tr.scalco < 0 {
                dx /= -(tr.scalco as f64);
                dy /= -(tr.scalco as f64);
            }
            gridrawxycdpic(&g, dx, dy, &mut c, &mut i, &mut j);
            if c < -2147483644 {
                err!(
                    "Error: input midpoint XYs not in grid (cannot compute cdp number). Trace= %d",
                    n + 1
                );
            }
            tr.cdp = c;
            tr.igi = i as i16;
            tr.igc = j as i16;
        } else if bt == -30 || bt == -31 {
            c = tr.cdp;
            gridcdpic(&g, c, &mut i, &mut j);
            if i < -2147483644 {
                err!(
                    "Error: input cdp number not in grid (cannot compute igi,igc). Trace= %d",
                    n + 1
                );
            }
            tr.igi = i as i16;
            tr.igc = j as i16;
        } else if bt == -32 {
            i = tr.igi as i32;
            j = tr.igc as i32;
            gridiccdp(&g, i, j, &mut c);
            if c < -2147483644 {
                err!(
                    "Error: input igi,igc numbers not in grid (cannot compute cdp). Trace= %d",
                    n + 1
                );
            }
            tr.cdp = c;
        } else {
            tx = fromhead(&tr, rp) * g[2] + g[1];
            ty = fromhead(&tr, sp) * g[4] + g[3];
            tr.cdp = ((tx + ty) / 2.).round_ties_even() as i32;
        }
        if bt == -30 || bt == -32 {
            gridicrawxy(&g, i, j, &mut dx, &mut dy);
            gridicgridxy(&g, i, j, &mut tx, &mut ty);
            if tr.scalco > 1 {
                let s = tr.scalco as f64;
                dx /= s;
                dy /= s;
                tx /= s;
                ty /= s;
            } else if tr.scalco < 0 {
                let s = -(tr.scalco as f64);
                dx *= s;
                dy *= s;
                tx *= s;
                ty *= s;
            }
            tr.sx = dx as i32;
            tr.sy = dy as i32;
            tr.gx = tx as i32;
            tr.gy = ty as i32;
        }
        puttr(&mut tr);
        n += 1;
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    warn!("Number of traces %d ", n);
    exit(CWP_Exit());
}
