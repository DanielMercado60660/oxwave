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
//! SUBINBIGCSV - Make some CDPs larger.
//! Port of su/main/headers/subinbigcsv.c (Revision 1.00).
//! Credits: Andre Latour.

use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                            ",
    " SUBINBIGCSV - Make Some CDPs Larger and Delete All Traces Outside Them.    ",
    "                                                                            ",
    " subinbigcsv  rfile=in.csv <in.su >out.su  (other parameters)               ",
    "                                                                            ",
    " Parameter overview:                                                        ",
    "                                                                            ",
    "  rfile=    If specified, read a K-file containing 3D Grid definition.      ",
    "            See subincsv for 3D Grid documentation.                         ",
    "            If a 3D Grid is not specified, assume input has 2D cdp numbers. ",
    "                                                                            ",
    " igiout=f,l,i  Output igi range. If specified, all 3 integers required.     ",
    "               Start at f, increment by i, end at l or nearest igi that     ",
    "               is not greater than maximum igi of Grid. See EXAMPLE below.  ",
    "           *** Or, specify the next igilist paramter.                       ",
    "                                                                            ",
    " igilist=      Output igi list. Or specify previous igiout parame3ter.      ",
    "               This parameter still produces locations which satisfy the    ",
    "               enclosing-rectangle requirements for bilinear interpolation. ",
    "                                                                            ",
    " igcout=f,l,i  Output igc range. If specified, all 3 integers required.     ",
    "               Start at f, increment by i, end at l or nearest igc that     ",
    "               is not greater than maximum igc of Grid. See EXAMPLE below.  ",
    "           *** Or, specify the next igclist paramter.                       ",
    "                                                                            ",
    " igclist=      Output igc list. Or specify previous igcout parame3ter.      ",
    "               This parameter still produces locations which satisfy the    ",
    "               enclosing-rectangle requirements for bilinear interpolation. ",
    "                                                                            ",
    " igiext=0   Output igi extent (default 0 means extent is only igiout itself)",
    "            List with one or two non-negative integers which form an        ",
    "            extent of igi values around the central igiout values.          ",
    "       =l     with only 1 value, h defaults to l (symmetric extent)         ",
    "       =l,h   start at l lower than each central igiout value, continue     ",
    "              until reaching h higher than each central igiout value.       ",
    "      Note: Must be less than the igiout increment.                         ",
    "                                                                            ",
    " igcext=0   Output igc extent (default 0 means extent is only igcout itself)",
    "            List with one or two non-negative integers which form an        ",
    "            extent of igc values around the central igcout values.          ",
    "       =l     with only 1 value, h defaults to l (symmetric extent)         ",
    "       =l,h   start at l lower than each central igcout value, continue     ",
    "              until reaching h higher than each central igcout value.       ",
    "      Note: Must be less than the igcout increment.                         ",
    "                                                                            ",
    " intype=0   Default 3 if a 3D Grid is input, otherwise -2.                  ",
    "       =3   Compute trace midpoint coordinates (sx+gx)/2 and (gx+gy)/2 and  ",
    "            determine cdp,igi,igc values from them (using the 3D Grid).     ",
    "            For this option traces do not need cdp,igi,igc numbers on input ",
    "            (they do not need to have been through grid program subincsv).  ",
    "      NOTE: For 3D STACKed traces, options 2 or 1 are required (usually).   ",
    "       =2   Use input cdp (cell) key and 3D Grid to determine igi,igc.      ",
    "       =1   Use input igi,igc keys and 3D Grid to determine cdp.            ",
    "       =-2  Use input cdp (cell) key. A 3D Grid cannot be input.            ",
    "      NOTE: For option -2 parameters igiout and igiext refer to cdp key     ",
    "            numbers and igcout, igcext, reigi, reigc cannot be specified    ",
    "            (reigi and reigc both use option 0).                            ",
    "                                                                            ",
    "  qout=     Output file containing the locations specified by the           ",
    "            igiout (or igilist) and igcout (or igclist) parameters.         ",
    "                                                                            ",
    "                                                                            ",
    " recdp=1    Reset the trace cdp key to its central cdp number.              ",
    "      =0    Do not reset cdp key numbers.                                   ",
    "                                                                            ",
    " reigi=1    Reset the trace igi key to its central igi number.              ",
    "      =0    Do not reset igi key numbers.                                   ",
    "                                                                            ",
    " reigc=1    Reset the trace igc key to its central igc number.              ",
    "      =0    Do not reset igc key numbers.                                   ",
    "                                                                            ",
    "   check=0  Do not print checking details.                                  ",
    "         1  Run functions on the 4 grid corner points and print results.    ",
    "            This output print may be useful for users.                      ",
    "                                                                            ",
    " EXAMPLE ****************************************************************** ",
    "   igiout=40,9999,40                                                        ",
    "   igcout=5,9999,5                                                          ",
    "   igiext=2                                                                 ",
    "   igcext=1                                                                 ",
    " This defines a central cdp every 40 cells in igi direction by 5 cells in   ",
    " igc direction. Traces will ONLY BE OUTPUT if they are within 2 cells in    ",
    " igi direction and also within 1 cell in igc direction of the central cdps. ",
    " By default, the cdp,igi,igc key values in the output traces will be reset  ",
    " to the values of the central cdp of those traces. One purpose of this is   ",
    " to create super-cdps for Velocity analysis.                                ",
    " ************************************************************************** ",
    " This program does not sort, it deletes traces, and renumbers cdp key.      ",
    " Since cdp is renumbered the output is NOT cdp ordered even if you input a  ",
    " cdp ordered dataset.                                                       ",
    " ************************************************************************** ",
    "                                                                            ",
    " ***********************************************************                ",
    "   To output this documentation:  subinbigcsv 2> binbigdoc.txt              ",
    " ***********************************************************                ",
    "                                                                            ",
    " ------------------------------------------------------------------------   ",
    "                                                                            ",
];

// The corresponding C arithmetic has undefined behavior when it overflows.
fn checked_int(value: i64) -> i32 {
    i32::try_from(value).unwrap_or_else(|_| err!("error: integer overflow in selection range."))
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut r = String::new();
    let hr = getparstring("rfile", &mut r) != 0;
    let mut mg = 0;
    gridcommand(&mut mg);
    if mg == 1 && hr {
        err!("error: input k-file not allowed when full grid on command line.");
    }
    if mg == -1 && !hr {
        err!("error: input k-file required when partial grid on command line.");
    }
    let is3d = mg != 0 || hr;
    let mut check = 0;
    getparint("check", &mut check);
    let mut g = [0.0; 999];
    if is3d {
        let fp = if mg != 1 {
            Some(
                fopen(&r, "r")
                    .unwrap_or_else(|| err!("error: input K-file did not open correctly.")),
            )
        } else {
            None
        };
        let mut e = 1;
        gridread(fp.as_ref(), &mut g, &mut e);
        if e > 0 {
            err!("error reading grid (from K-file or command line parameters)");
        }
        gridset(&mut g, &mut e);
        match e {
            1 => err!("gridset error: grid_wb cell width must be positive."),
            2 => err!("gridset error: grid_wc cell width must be positive."),
            3 => err!("gridset error: corner B is within grid_wb cell width of corner A."),
            n if n > 0 => err!("gridset error: returned with some unrecognized error code."),
            -1 => warn!("gridset warning: corner C is near A and is reset to A."),
            _ => (),
        }
        gridcheck(&g, check, &mut e);
        if e > 0 {
            err!("gridcheck error: returned with some unrecognized error code.");
        }
    }
    let mut typ = 0;
    getparint("intype", &mut typ);
    if typ == 0 {
        typ = if is3d { 3 } else { -2 };
    }
    if ![-2, 1, 2, 3].contains(&typ) {
        err!("error: intype= option not in range.");
    }
    let (mut ri, mut rc) = (1, 1);
    if typ == -2 {
        if is3d {
            err!("error: intype -2 not allowed when a 3D Grid is input");
        }
        if countparval("igcout") > 0 {
            err!("error: igcout= cannot be specified when intype is -2 (no 3D Grid).");
        }
        if countparval("igcext") > 0 {
            err!("error: igcext= cannot be specified when intype is -2 (no 3D Grid).");
        }
        if countparval("igclist") > 0 {
            err!("error: igclist= cannot be specified when intype is -2 (no 3D Grid).");
        }
        if countparval("reigi") > 0 {
            err!("error: reigi= cannot be specified when intype is -2 (no 3D Grid).");
        }
        if countparval("reigc") > 0 {
            err!("error: reigc= cannot be specified when intype is -2 (no 3D Grid).");
        }
    } else {
        if !is3d {
            err!("error: only intype -2 allowed when no 3D Grid is input");
        }
        getparint("reigi", &mut ri);
        if !(0..=1).contains(&ri) {
            err!("error: reigi= option not in range.");
        }
        getparint("reigc", &mut rc);
        if !(0..=1).contains(&rc) {
            err!("error: reigc= option not in range.");
        }
        if countparval("igclist") > 0 && countparval("igcout") > 0 {
            err!("error: igcout= and igclist= cannot both be specified.");
        }
        if countparval("igcext") > 2 {
            err!("error: igcext= has too long a list.");
        }
    }
    let mut rd = 1;
    getparint("recdp", &mut rd);
    if !(0..=1).contains(&rd) {
        err!("error: recdp= option not in range.");
    }
    if countparval("igilist") > 0 && countparval("igiout") > 0 {
        err!("error: igiout= and igilist= cannot both be specified.");
    }
    if countparval("igiext") > 2 {
        err!("error: igiext= has too long a list.");
    }
    let (mut li, mut lc) = (Vec::<i32>::new(), Vec::<i32>::new());
    if countparval("igilist") > 0 {
        getparint("igilist", &mut li);
    } else if countparval("igiout") != 3 {
        err!("error: igiout= must have 3 integers in list.");
    }
    if is3d {
        if countparval("igclist") > 0 {
            getparint("igclist", &mut lc);
        } else if countparval("igcout") != 3 {
            err!("error: igcout= must have 3 integers in list.");
        }
    }
    let (mut ei, mut ec) = (vec![0, 0], vec![0, 0]);
    if countparval("igiext") > 0 {
        getparint("igiext", &mut ei);
    }
    if countparval("igiext") == 1 {
        ei[1] = ei[0];
    }
    if ei[0] < 0 || ei[1] < 0 {
        err!("error: igiext= values cannot be less than 0");
    }
    if countparval("igcext") > 0 {
        getparint("igcext", &mut ec);
    }
    if countparval("igcext") == 1 {
        ec[1] = ec[0];
    }
    if ec[0] < 0 || ec[1] < 0 {
        err!("error: igcext= values cannot be less than 0");
    }
    checked_int(ei[0] as i64 + ei[1] as i64 + 1);
    checked_int(ec[0] as i64 + ec[1] as i64 + 1);
    if is3d {
        warn!("igi extent total cells=%d, total width=%g    igc extent total cells=%d, total width=%g",ei[0]+ei[1]+1,(ei[0]+ei[1]+1) as f64*g[10],ec[0]+ec[1]+1,(ec[0]+ec[1]+1) as f64*g[11]);
    } else {
        warn!("igi extent total cells=%d.", ei[0] + ei[1] + 1);
    }
    let (mut mi, mut ma, mut mc, mut mx) = if is3d {
        (1, (0.1 + g[12]) as i32, 1, (0.1 + g[13]) as i32)
    } else {
        (-2000000000, 2000000000, 1, 1)
    };
    let (mut oi, mut oc) = (vec![mi, ma, 1], vec![mc, mx, 1]);
    if li.is_empty() {
        getparint("igiout", &mut oi);
        if oi[2] < 1 {
            err!("error: igiout= increment is less than 1");
        }
        if is3d {
            if oi[0] < mi {
                err!(
                    "error: igiout= first is less than minimum grid igi index %d",
                    mi
                );
            }
            if oi[1] > ma {
                oi[1] = ma;
            }
        }
        if oi[0] > oi[1] {
            err!("error: igiout= first is greater than last. ");
        }
        if ei[0] + ei[1] >= oi[2] {
            err!("error: igiext= total igi extent must be less than igiout increment.");
        }
        mi = checked_int(oi[0] as i64 - ei[0] as i64);
        let span = checked_int(oi[1] as i64 - oi[0] as i64);
        let steps = ((span as f64 + 0.1) / oi[2] as f64) as i32;
        let offset = checked_int(oi[2] as i64 * steps as i64);
        let last = checked_int(oi[0] as i64 + offset as i64);
        ma = checked_int(last as i64 + ei[1] as i64);
        ei[0] = oi[2] - ei[0] - 1;
        ei[1] += 1;
    } else {
        let (mut lo, mut hi) = (2000000000, -2000000000);
        for i in 0..li.len() {
            lo = lo.min(li[i]);
            hi = hi.max(li[i]);
            for j in i + 1..li.len() {
                if (li[i] as i64 - li[j] as i64).abs() < ei[0] as i64 + ei[1] as i64 + 1 {
                    err!("error: igilist= has values that are too close to each other (extents overlap)");
                }
            }
        }
        if is3d {
            if lo < mi {
                err!(
                    "error: igilist= has value less than minimum grid igi index %d",
                    mi
                );
            }
            if hi > ma {
                err!(
                    "error: igilist= has value greater than maximum grid igi index %d",
                    ma
                );
            }
        }
        mi = checked_int(lo as i64 - ei[0] as i64);
        ma = checked_int(hi as i64 + ei[1] as i64);
        ei[0] += 1;
        ei[1] += 1;
    }
    if is3d {
        if lc.is_empty() {
            getparint("igcout", &mut oc);
            if oc[2] < 1 {
                err!("error: igcout= increment is less than 1");
            }
            if oc[0] < mc {
                err!(
                    "error: igcout= first is less than minimum grid igc index %d",
                    mc
                );
            }
            if oc[1] > mx {
                oc[1] = mx;
            }
            if oc[0] > oc[1] {
                err!("error: igcout= first is greater than last. ");
            }
            if ec[0] + ec[1] >= oc[2] {
                err!("error: igcext= total igc extent must be less than igcout increment.");
            }
            mc = checked_int(oc[0] as i64 - ec[0] as i64);
            let span = checked_int(oc[1] as i64 - oc[0] as i64);
            let steps = ((span as f64 + 0.1) / oc[2] as f64) as i32;
            let offset = checked_int(oc[2] as i64 * steps as i64);
            let last = checked_int(oc[0] as i64 + offset as i64);
            mx = checked_int(last as i64 + ec[1] as i64);
            ec[0] = oc[2] - ec[0] - 1;
            ec[1] += 1;
        } else {
            let (mut lo, mut hi) = (2000000000, -2000000000);
            for i in 0..lc.len() {
                lo = lo.min(lc[i]);
                hi = hi.max(lc[i]);
                for j in i + 1..lc.len() {
                    if (lc[i] as i64 - lc[j] as i64) < ec[0] as i64 && (lc[j] as i64 - lc[i] as i64) < ec[1] as i64 {
                        err!("error: igclist= has values that are too close to each other (extents overlap)");
                    }
                }
            }
            if lo < mc {
                err!(
                    "error: igclist= has value less than minimum grid igc index %d",
                    mc
                );
            }
            if hi > mx {
                err!(
                    "error: igclist= has value greater than maximum grid igc index %d",
                    mx
                );
            }
            mc = checked_int(lo as i64 - ec[0] as i64);
            mx = checked_int(hi as i64 + ec[1] as i64);
            ec[0] += 1;
            ec[1] += 1;
        }
    }
    let mut q = String::new();
    let hq = getparstring("qout", &mut q) != 0;
    checkpars();
    if hq {
        let fp = fopen(&q, "w").unwrap_or_else(|| err!("qout error: file did not open correctly."));
        fputs("C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID", &fp);
        if is3d {
            fputs(",%.20g,%.20g,%.20g,%.20g,%.2f,%.2f,%.2f,%.2f", &fp);
        } else {
            fputs(",%.20g,%.20g", &fp);
        }
        fputs("\nC_SU_NAMES\nC_SU_ID", &fp);
        if is3d {
            fputs(",cdp,cdpt,igi,igc,gx,gy,sx,sy\n", &fp);
        } else {
            fputs(",cdp,cdpt\n", &fp);
        }
        let ni = if li.is_empty() {
            checked_int((oi[1] as i64 - oi[0] as i64) / oi[2] as i64 + 1)
        } else {
            li.len() as i32
        };
        let nc = if lc.is_empty() {
            checked_int((oc[1] as i64 - oc[0] as i64) / oc[2] as i64 + 1)
        } else {
            lc.len() as i32
        };
        for i in 0..ni {
            let a = if li.is_empty() {
                checked_int(oi[0] as i64 + i as i64 * oi[2] as i64)
            } else {
                li[i as usize]
            };
            if is3d {
                for j in 0..nc {
                    let b = if lc.is_empty() {
                        checked_int(oc[0] as i64 + j as i64 * oc[2] as i64)
                    } else {
                        lc[j as usize]
                    };
                    let (mut d, mut t, mut x, mut y, mut u, mut v) = (0, 0, 0.0, 0.0, 0.0, 0.0);
                    gridiccdp(&g, a, b, &mut d);
                    gridiccdp90(&g, a, b, &mut t);
                    gridicgridxy(&g, a, b, &mut x, &mut y);
                    gridicrawxy(&g, a, b, &mut u, &mut v);
                    fprintf!(
                        fp,
                        "Q,%d,%d,%d,%d,%.2f,%.2f,%.2f,%.2f\n",
                        d,
                        t,
                        a,
                        b,
                        x,
                        y,
                        u,
                        v
                    );
                }
            } else {
                fprintf!(fp, "Q,%d,1\n", a);
            }
        }
        fclose(&fp);
    }
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Error: cannot get first trace");
    }
    let (mut n, mut m) = (0, 0);
    loop {
        n += 1;
        let (mut i, mut c, mut d) = (1, 1, 1);
        if typ == 3 {
            let mut x = 0.5 * tr.sx.wrapping_add(tr.gx) as f64;
            let mut y = 0.5 * tr.sy.wrapping_add(tr.gy) as f64;
            if tr.scalco > 1 {
                x *= tr.scalco as f64;
                y *= tr.scalco as f64;
            } else if tr.scalco < 0 {
                x /= -(tr.scalco as f64);
                y /= -(tr.scalco as f64);
            }
            gridrawxycdpic(&g, x, y, &mut d, &mut i, &mut c);
            if d < -2147483644 {
                err!(
                    "Error: input midpoint XYs not in grid (cannot compute cdp number). Trace= %d",
                    n
                );
            }
        } else if typ == 2 {
            d = tr.cdp;
            gridcdpic(&g, d, &mut i, &mut c);
            if i < -2147483644 {
                err!(
                    "Error: input cdp number not in grid (cannot compute igi,igc). Trace= %d",
                    n
                );
            }
        } else if typ == 1 {
            i = tr.igi as i32;
            c = tr.igc as i32;
            gridiccdp(&g, i, c, &mut d);
            if d < -2147483644 {
                err!(
                    "Error: input igi,igc numbers not in grid (cannot compute cdp).      Trace= %d",
                    n
                );
            }
        } else {
            i = tr.cdp;
        }
        if i >= mi && i <= ma && c >= mc && c <= mx {
            let (mut ki, mut kc) = (false, true);
            if li.is_empty() {
                let z = checked_int(i as i64 - oi[0] as i64) % oi[2];
                if z > ei[0] {
                    ki = true;
                    i = checked_int(checked_int(i as i64 + oi[2] as i64) as i64 - z as i64);
                } else if z < ei[1] {
                    ki = true;
                    i = checked_int(i as i64 - z as i64);
                }
            } else {
                for &v in &li {
                    if (v as i64 - i as i64) < ei[0] as i64 && (i as i64 - v as i64) < ei[1] as i64 {
                        ki = true;
                        i = v;
                        break;
                    }
                }
            }
            if typ != -2 {
                kc = false;
                if lc.is_empty() {
                    let z = checked_int(c as i64 - oc[0] as i64) % oc[2];
                    if z > ec[0] {
                        kc = true;
                        c = checked_int(checked_int(c as i64 + oc[2] as i64) as i64 - z as i64);
                    } else if z < ec[1] {
                        kc = true;
                        c = checked_int(c as i64 - z as i64);
                    }
                } else {
                    for &v in &lc {
                        if (v as i64 - c as i64) < ec[0] as i64 && (c as i64 - v as i64) < ec[1] as i64 {
                            kc = true;
                            c = v;
                            break;
                        }
                    }
                }
            }
            if ki && kc {
                if typ == -2 {
                    if rd > 0 {
                        tr.cdp = i;
                    }
                } else {
                    if rd > 0 {
                        gridiccdp(&g, i, c, &mut d);
                        if d < -2147483644 {
                            err!("Error: input igi,igc numbers not in grid (cannot compute cdp). Trace= %d   %d",n,d);
                        }
                        tr.cdp = d;
                    }
                    if ri > 0 {
                        tr.igi = i as i16;
                    }
                    if rc > 0 {
                        tr.igc = c as i16;
                    }
                }
                puttr(&mut tr);
                m += 1;
            }
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    warn!("Number of traces input=%d.  Number output=%d.", n, m);
    exit(CWP_Exit());
}
