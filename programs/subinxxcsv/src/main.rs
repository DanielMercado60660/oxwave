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

//! SUBINXXCSV - expand and exchange corners of a 3D grid.
//! Port of su/main/headers/subinxxcsv.c (Revision 1.01).
//! Credits: Andre Latour, December 2021.

use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;
use oxwave_k_file::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUBINXXCSV - Expand Size and Exchange Corners of 3D Grid Definition.       ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  subinxxcsv <stdin >stdout [parameters]                                    ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t\t\t          \t\t\t\t     ",
    "                                                                            ",
    " rfile=   If specified, read a K-file containing a 3D grid definition.      ",
    "          Alternately you can supply the grid definition via command line   ",
    "          parameters. The required command line parameters are:             ",
    "           grid_xa,grid_ya,grid_xb,grid_yb,grid_xc,grid_yc,grid_wb,grid_wc. ",
    "          Note that input grid definitions always run a routine which       ",
    "          adjusts the input values into compliance with grid standards.     ",
    "          (See program SUBINCSV for 3D grid details).                       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " wfile=   File to write the new grid definition after changes caused by     ",
    "          the remaining parameters. Must be specified.                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " igilow=0  Add cells at low igi end (default is not to add cells).          ",
    "       =n  Add n cells at low igi end (negative makes grid smaller).        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " igihigh=0 Add cells at high igi end (default is not to add cells).         ",
    "        =n Add n cells at high igi end (negative makes grid smaller).       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " igclow=0  Add cells at low igc end (default is not to add cells).          ",
    "       =n  Add n cells at low igc end (negative makes grid smaller).        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " igchigh=0 Add cells at high igc end (default is not to add cells).         ",
    "        =n Add n cells at high igc end (negative makes grid smaller).       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " iflip=0   Defaults to not exchanging any grid corners.                     ",
    "           This option is done AFTER any add-cell options listed above.     ",
    "      =1   Exchange corner B with corner C and also exchange cell widths.   ",
    "           This results in cdp numbers incrementing by 1 at right-angle     ",
    "           to input grid and also exchanges igi,igc directions.             ",
    "      =2   Exchange corner A with corner D and also exchange cell widths.   ",
    "           This results in cdp numbers incrementing by 1 at right-angle     ",
    "           to input grid and also exchanges igi,igc directions. But cdp     ",
    "           numbering will now start where input corner D was.               ",
    "      =-1  Exchange corner A with corner B and corner C with corner D.      ",
    "      =-2  Exchange corner A with corner C and corner B with corner D.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " check=0   Do not print grid checking details.                              ",
    "     Note: The options of this program cannot repair a corrupt input grid.  ",
    "           Grid checking is always done on input grid definition and will   ",
    "           error-halt if corrupt.                                           ",
    "      =1   print grid checking details. This runs some grid functions on    ",
    "           the 4 corner points of input grid, expanded grid (if any), and   ",
    "           flipped grid (if any). This information can be written to a file ",
    "           by putting 2>yourfile on command line.                           ",
    "\t\t\t\t\t\t\t\t\t     ",
    " CAUTION: Any grid change is likely to result in some pre-stack traces      ",
    "          being assigned to different cells (cdps) than pure mathematics    ",
    "          would predict. That is, NEVER EXPECT the group of traces assigned ",
    "          to each cdp by the input grid definition to be exactly the same   ",
    "          as the group of traces assigned to any corresponding cdps by the  ",
    "          output grid definition. This is true despite computations herein  ",
    "          being performed with 8 byte floating point (double precision).    ",
    "          This issue occurs because some trace midpoints lie exactly on     ",
    "          cell boundaries - and even the slightest change in the grid       ",
    "          definition results in slight adjustments of all cell boundaries.  ",
    "            So, you should make grid adjustments before you start ANY trace ",
    "            processing, or re-run all trace processes, including, but not   ",
    "            limited to, cdp sorting.                                        ",
    "          Notwithstanding this caution, you have a reasonable chance of     ",
    "          using an altered grid post-stack because the post-stack cdps are  ",
    "          essentially located at the cell centre XYs, and cell centres are  ",
    "          usually far away from your output grid cell boundaries.           ",
    "\t\t\t\t\t\t\t\t\t     ",
];

fn set_and_check(g: &mut [f64], check: i32) {
    let mut e = 0;
    gridset(g, &mut e);
    match e {
        1 => err!("gridset error: grid_wb cell width must be positive."),
        2 => err!("gridset error: grid_wc cell width must be positive."),
        3 => err!("gridset error: corner B is within grid_wb cell width of corner A."),
        x if x > 0 => err!("gridset error: returned with some unrecognized error code."),
        -1 => warn!("gridset warning: corner C is near A and is reset to A."),
        _ => ()
    }
    gridcheck(g, check, &mut e);
    if e > 0 { err!("gridcheck error: returned with some unrecognized error code."); }
}
fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    if !isatty(0) || !isatty(1) { err!("**** Error: this program does not input or output traces."); }
    let (mut r, mut w) = (String::new(), String::new());
    let has_r = getparstring("rfile", &mut r) != 0;
    let has_w = getparstring("wfile", &mut w) != 0;
    if has_r && has_w && r == w { err!("**** Error: wfile= output K-file must be different than rfile= input K-file."); }
    if !has_r && has_w { err!("**** Error: wfile= cannot output K-file with no rfile= input K-file."); }
    let mut check = 0; getparint("check", &mut check);
    if !(0..=1).contains(&check) { err!("**** Error: check= value out of range."); }
    let mut il = 0; getparint("igilow", &mut il);
    let mut ih = il; getparint("igihigh", &mut ih);
    let mut cl = 0; getparint("igclow", &mut cl);
    let mut ch = cl; getparint("igchigh", &mut ch);
    let mut flip = 0; getparint("iflip", &mut flip);
    if !(-2..=2).contains(&flip) { err!("**** Error: iflip= value out of range."); }
    let mut maygrid = 0; gridcommand(&mut maygrid);
    if maygrid == 1 && has_r { err!("error: input k-file not allowed when full grid on command line."); }
    if maygrid == -1 && !has_r { err!("error: input k-file required when partial grid on command line."); }
    let mut names = vec![String::new(); 999]; let mut forms = names.clone();
    let mut fields = vec![0.0; 999]; let mut count = 0; let mut e = 0;
    if maygrid != 1 {
        let fp = fopen(&r, "r").unwrap_or_else(|| err!("error: input K-file did not open correctly."));
        readkfile(&fp, &mut names, &mut forms, &mut fields, &mut count, &mut e);
        match e {
            1 => err!("K-file read error: more than one C_SU_NAMES parameter record."),
            2 => err!("K-file read error: no C_SU_NAMES parameter record."),
            3 => err!("K-file read error: more than one C_SU_FORMS parameter record."),
            4 => err!("K-file read error: no C_SU_FORMS parameter record."),
            5 => err!("K-file read error: more than one C_SU_SETID record."),
            6 => err!("K-file read error: no C_SU_SETID record."),
            7 => err!("K-file read error: different number of values on C_SU_NAMES and C_SU_FORMS."),
            8 => err!("K-file read error: unable to allocate memory."),
            9 => err!("K-file read error: name exists at least twice in C_SU_NAMES list."),
            10 => err!("K-file read error: at least 1 field-unreadable as a number."),
            11 => err!("K-file read error: at least 1 field containing 2 numbers."),
            12 => err!("K-file read error: not-enough-commas to get all values for C_SU_NAMES list."),
            x if x > 0 => err!("K-file read error: returned with some unrecognized error code."),
            -1 => warn!("K-file read warning: at least 1 all-blank field, assumed zero for all."), _ => ()
        }
        for name in &mut names[..count as usize] { name.make_ascii_lowercase(); }
    }
    let gn = ["bintype", "grid_lf", "grid_xa", "grid_ya", "grid_xb", "grid_yb", "grid_xc", "grid_yc", "grid_xd", "grid_yd", "grid_wb", "grid_wc", "grid_nb", "grid_nc", "grid_fp", "grid_lp", "grid_sb", "grid_cb"];
    let mut g = [0.0; 18];
    // The C leaves required entries uninitialized; make its missing-value
    // check deterministic while retaining the zero default for bintype.
    for i in [2, 3, 4, 5, 6, 7, 10, 11] { g[i] = -1.1e308; }
    for i in 0..12 {
        if i == 1 || i == 8 || i == 9 { continue; }
        if getpardouble(gn[i], &mut g[i]) == 0 {
            for j in 0..count as usize { if names[j] == gn[i] { g[i] = fields[j]; } }
            if g[i] < -1.0e308 { g[i] = (i + 100) as f64; if i > 0 { err!("Error: parameter %s not found.", gn[i]); } }
        }
    }
    checkpars(); set_and_check(&mut g, check);
    if il != 0 || ih != 0 || cl != 0 || ch != 0 {
        let mut b = g;
        let (lo, hi) = (1i32.wrapping_sub(il), ih.wrapping_add((g[12] + 0.1) as i32));
        let (lc, hc) = (1i32.wrapping_sub(cl), ch.wrapping_add((g[13] + 0.1) as i32));
        let (mut x, mut y) = (0.0, 0.0);
        gridicrawxy(&g, lo, lc, &mut x, &mut y); b[2]=x; b[3]=y;
        gridicrawxy(&g, hi, lc, &mut x, &mut y); b[4]=x; b[5]=y;
        gridicrawxy(&g, lo, hc, &mut x, &mut y); b[6]=x; b[7]=y;
        g = b; set_and_check(&mut g, check);
    }
    if flip != 0 {
        let mut b = g;
        match flip {
            1 => { b[4]=g[6]; b[5]=g[7]; b[6]=g[4]; b[7]=g[5]; b[10]=g[11]; b[11]=g[10]; },
            2 => { b[2]=g[8]; b[3]=g[9]; b[10]=g[11]; b[11]=g[10]; },
            -1 => { b[2]=g[4]; b[3]=g[5]; b[4]=g[2]; b[5]=g[3]; b[6]=g[8]; b[7]=g[9]; },
            -2 => { b[2]=g[6]; b[3]=g[7]; b[6]=g[2]; b[7]=g[3]; b[4]=g[8]; b[5]=g[9]; }, _ => ()
        }
        g = b; set_and_check(&mut g, check);
    }
    let fp = fopen(&w, "w").unwrap_or_else(|| err!("wfile error: output K-file did not open correctly."));
    let mut out = count as usize;
    for i in 0..18 {
        let mut found = false;
        for j in 0..count as usize { if names[j] == gn[i] { fields[j] = g[i]; found = true; break; } }
        if !found { fields[out]=g[i]; names[out]=gn[i].into(); forms[out]="%.20g".into(); out += 1; }
    }
    writekfile(&fp, &names[..out], &forms[..out], &fields[..out], &mut e);
    if e > 0 { err!("K-file write error: returned with an unrecognized error code."); }
    exit(CWP_Exit());
}
