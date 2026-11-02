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

//! SUPROFCSV - port of su/main/headers/suprofcsv.c (Revision 1.02).
//! Copyright (c) Colorado School of Mines, 2022.
//! All rights reserved.
//! Created June 2022 by Andre Latour; started from subinqcsv.
//! Modified November 2023 by Andre Latour.

use oxwave_core::prelude::*;
use oxwave_q_file::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUPROFCSV - Input Q-file or parameters, make profile, output other Q-file. ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  suprofcsv [parameters].   (No traces in or out).                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Note: You should start by looking at the 2d cdp numbering computation that ",
    "       exists in the subincsv program. For 2d lines with consistent shot    ",
    "       and receiver coordinates and bends of less than 10 degrees or so,    ",
    "       subincsv can usually give nicer results for people who are not very  ",
    "       familiar with crooked profiles (and even for people who are).        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Note: All defaults here are roughtly set for a simple situation consisting ",
    "       of an input profile of receiver XY locations with a 25 metre spacing.",
    "       Those receivers might have come from an SPS2 R-file converted to     ",
    "       a q-file using program sutoolcsv.                                    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " qin=       Input file containing q-records. This parameter is optional,    ",
    "            but if you do not specify it, you must specify a set of         ",
    "            parameters which you can see lower down in this document.       ",
    "            Those parameters allow specification similar to sunmo.          ",
    "      Note: The keyp,keyx,keyy parameters recognize root key names even if  ",
    "            they are surrounded by underscores and digits in the qin file.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keyp=gaps Point Key Name. The default (gaps) is the receiver point name    ",
    "           used by sugeomcsv and sutoolcsv.                                 ",
    "           Name must be in the input Q-file (or parameters listed below).   ",
    "           This is used to order (sort) the input locations based on point  ",
    "           numbers, station numbers, receiver numbers, or shot numbers.     ",
    "     =asis Means use the points in the order they are input.                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keyx=gx   X coordinate name. The default (gx) is the receiver x name.      ",
    "           Name must be in the input Q-file (or parameters listed below).   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keyy=gy   Y coordinate name. The default (gy) is the receiver y name.      ",
    "           Name must be in the input Q-file (or parameters listed below).   ",
    "     =none Means use only keyx as coordinate of points. Usually this means  ",
    "           that keyx is a value such as a station number, or a cdp number,  ",
    "           or the accumulated distance between points.                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " chordi=12.5   Initial Chord Distance. Start at lowest point and linearly   ",
    "           interplate new points at increments of this distance using the   ",
    "           surrounding two input points. Original points are not preserved. ",
    "           For typical situations, this should be your desired cpd spacing. ",
    "           Note: This chording occurs before profile averaging or smoothing.",
    "       =0.0 is a flag that means to not perform initial chording. So the    ",
    "            input is the initial profile. Generally this means the input    ",
    "            points should already have a constant interval and profile ends ",
    "            should have more than max(nmaxa,nmaxs) spare points.            ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  nmina=0  Minimum Averaging Number. Must be less or equal to nmaxa.        ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  nmaxa=10   Maximum Averaging Number.                                      ",
    "           A value of 0 means no averaging is done.                         ",
    "           Each location in the current profile is averaged symmetrically   ",
    "           with its neighbours (from nmina to nmaxa on each side).          ",
    "           For an input profile with a simple bend, this results in a       ",
    "           new profile towards the inside of the bend. The smaller this     ",
    "           value, the nearer the resulting profile follows input profile.   ",
    "           For typical surveys with shots located approximately along the   ",
    "           input (receiver) profile, the maximum reasonable value here is   ",
    "           nominal maximum source-receiver offset dist divided by chordi.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " dextra=nmaxa*chordi  Extra distance to extrapolate before Initial Chording.",
    "           This linear extrapolation uses the first 2 input points and the  ",
    "           last 2 input points. Cannot be specified if chordi=0.0           ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  nmins=0  Minimum Centred Smoothing Number. Must be less or equal to nmaxs.",
    "           Note: Usually should be 0 even if nmina is non-zero.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  nmaxs=nmaxa   Maximum Centred Smoothing Number.                           ",
    "           A value of 0 means no smoothing is done.                         ",
    "           Each location in the current profile is smoothed symmetrically   ",
    "           with its neighbours (from nmins to nmaxs on each side).          ",
    "           This smoothing uses a technique which leaves the profile         ",
    "           approximately centred along the same positions.                  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " chordf=chordi   Final Chord Distance. Start at lowest point and linearly   ",
    "           interpolate new points at increments of this distance using the  ",
    "           surrounding two points. For typical situations, this should be   ",
    "           your desired cpd spacing.                                        ",
    "       =0.0 is a flag that means to not perform final chording.             ",
    "           Note: This chording is done after initial chording, averaging,   ",
    "                 and smoothing.                                             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " ename=     Names to eliminate. These names will be ignored in the input    ",
    "            q-records and therefore not be in the output q-records.         ",
    "            Recommend not eliminating any (unless filesize is an issue).    ",
    "      Note: Must be the full names, including any underscores and digits.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " formxy=%.20g  The C format code for printing all values to the q-records.  ",
    "              Note that the default format prints up to 20 digits           ",
    "              (but not trailing zeroes to the right of the decimal point).  ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "            The following two parameters default to doing nothing. To use   ",
    "            them, view/plot the cdp locations after this program and compare",
    "            them to trace midpoint locations. Sometimes it is beneficial to ",
    "            align the cdp locations to the centre of clumps of midpoints    ",
    "            (or align them to midpoints with small source-receiver offsets).",
    "              For more information, consult src/demos/Geom3D/Suprofcsv      ",
    "                                    or my YouTube video about profiles      ",
    " shiftrec=  List of Q-record numbers to apply small shifts along profile.   ",
    "            The first record starting with Q in the qin file is record 1.   ",
    "            If the Q-file was previously output by this program,            ",
    "            and Q-records were not deleted or rearranged, the cdp numbers   ",
    "            in that Q-file match the record numbers to specify here.        ",
    " shiftdist= List of distances to shift the Q-records in shiftrec. These     ",
    "            shifts are along the profile. Positive shifts are in direction  ",
    "            of increasing Q-record numbers. The amount of shift for the     ",
    "            records between your specified points is partitioned according  ",
    "            to the distances between records. Note that, if records between ",
    "            your specified points have a constant interval on input, they   ",
    "            will have a different but still constant interval on output.    ",
    "            The shifts specified for the lowest and highest records are     ",
    "            held constant outward to their ends of the profile.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  ***  Note that the following two parameters also exist in subincsv.       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " point_crz=1.0  cdp number of point zero.                                   ",
    "                Note: NOT cdp of first output point unless point_cru=0.0    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " point_cru=0.0  cdps per one point unit.     Must be 0.0 if keyp=asis.      ",
    "            The first output cdp is set to point_crz plus first output      ",
    "            value of your specified keyp name multiplied by point_cru.      ",
    "              (first cdp) = (first keyp) * point_cru + point_crz            ",
    "       ***  Subsequent output cdps just increment by 1.                     ",
    "            Subsequent cdps are NOT computed using output keyp values and   ",
    "            will deviate from cdp numbers computed from output keyp values  ",
    "            unless you use the setup described next.                        ",
    "            If you have specified these parameters in subincsv then you     ",
    "            can output similar cdp numbers here - using some specific setup.",
    "             - keyp here should be rpkey name in subincsv (typically gaps). ",
    "             - keyx here should be the same as keyp (typically gaps).       ",
    "             - keyy here should be none.                                    ",
    "             - chordi here should be 1/point_cru (typically 0.5)            ",
    "            You can still average and smooth other values using this setup. ",
    "            (But for this setup the keyp values are linear and so they will ",
    "            be the same values because averaging/smoothing are symmetrical).",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qout=      Output file for q-records. Must be specified.                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    "      Note: A sequential cdp number is always generated in the qout file.   ",
    "      Note: The name cdp can be used from the input qin file for parameters ",
    "            keyp or keyx or keyy, but its name will be null in qout file    ",
    "            because another cdp number is created by this program.          ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "   ------------------------------------------------------------------       ",
    "   The following parameters can only be specified if you                    ",
    "   do not specify an input q-file using parameter qin=.                     ",
    " *** For consistency, these are the same parameters as for subinqcsv. ***   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cdp=             The parameters to the left are all optional but at least  ",
    " cdpt=            one of them must be specified. These are the names which  ",
    " fldr=            contain single values at each location (also known as     ",
    " grnors=    \t   the non-tuple names). Typically, gaps is the most likely  ",
    " grnofr=    \t   name to specify. Each name specified here must have the   ",
    " grnlof=    \t   same amount of values (separated by commas).              ",
    " gaps=      \t   For example: gaps=2,7,23,44 means that you must also      ",
    " igi=       \t   list 4 values if you specify numa=.                       ",
    " igc=       \t                                                             ",
    " gx=                                                                        ",
    " gy=                                                                        ",
    " sx=                                                                        ",
    " sy=                                                                        ",
    " numa=                                                                      ",
    " numb=                                                                      ",
    " vuma=                                                                      ",
    " vumb=                                                                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " tupa=            The parameters to the left are all optional. These are    ",
    " offs=            the names which have multiple values at each location     ",
    " tims=            (also known as the tuple names). Typically, tims= and     ",
    " tnmo=            vnmo= are the most likely names to specify. Each name     ",
    " dpth=      \t   specified here must be repeated same number of times as   ",
    " vels=            amount of values in non-tuple names above (except Note 1).",
    " vnmo=                                                                      ",
    " tupb=    Note 1: Usually, any of the previous parameter names must be      ",
    " vupa=            repeated for each value in the non-tuple lists. But the   ",
    "                  independent dimension name can also be specified just one ",
    "                  time if all the other names here have the same amount of  ",
    "                  values as that single set of independent dimension values.",
    "          Note 2: The independent dimension defaults to first name that is  ",
    "            \t   found in order here (NOT their order on the command line).",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "\t\t\t\t\t\t\t\t\t     ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    if !isatty(0) || !isatty(1) {
        err!("**** Error: this program does not input or output traces.");
    }
    let mut keyp = "gaps".to_owned();
    let mut keyx = "gx".to_owned();
    let mut keyy = "gy".to_owned();
    if countparval("keyp") > 0 {
        getparstring("keyp", &mut keyp);
    }
    if countparval("keyx") > 0 {
        getparstring("keyx", &mut keyx);
    }
    if countparval("keyy") > 0 {
        getparstring("keyy", &mut keyy);
    }
    let mut chordi = 12.5;
    getpardouble("chordi", &mut chordi);
    if chordi < 0.0 {
        err!("**** Error: chordi= must be greater than or equal to 0.0 ");
    }
    let mut nmina = 0;
    let mut nmaxa = 10;
    getparint("nmina", &mut nmina);
    getparint("nmaxa", &mut nmaxa);
    if nmaxa < 0 {
        err!("**** Error: nmaxa= must be greater or equal to 0 ");
    }
    if nmina < 0 || nmina > nmaxa {
        err!("**** Error: nmaxa= must be greater or equal to nmina= ");
    }
    let mut dextra = nmaxa as f64 * chordi;
    getpardouble("dextra", &mut dextra);
    if chordi == 0.0 && dextra != 0.0 {
        err!("**** Error: dextra= must be 0.0 if chordi=0.0 ");
    }
    if dextra < 0.0 {
        err!("**** Error: dextra= cannot be negative. ");
    }
    let mut nmins = 0;
    let mut nmaxs = nmaxa;
    getparint("nmins", &mut nmins);
    getparint("nmaxs", &mut nmaxs);
    if nmaxs < 0 {
        err!("**** Error: nmaxs= must be greater or equal to 0 ");
    }
    if nmins < 0 || nmins > nmaxs {
        err!("**** Error: nmaxs= must be greater or equal to nmins= ");
    }
    let mut chordf = chordi;
    getpardouble("chordf", &mut chordf);
    if chordf < 0.0 {
        err!("**** Error: chordf= must be greater than or equal to 0.0 ");
    }
    let mut pcrz = 1.0;
    let mut pcru = 0.0;
    getpardouble("point_crz", &mut pcrz);
    getpardouble("point_cru", &mut pcru);
    if keyp == "asis" && pcru != 0.0 {
        err!("**** Error: point_cru is not 0.0 but keyp=asis ");
    }
    let mut qin = String::new();
    let has_qin = getparstring("qin", &mut qin) != 0;
    let mut qout = String::new();
    let has_qout = getparstring("qout", &mut qout) != 0;
    // C compares names only when both filename pointers are present.
    if has_qin && has_qout && qin == qout {
        err!("**** Error: qout= output Q-file name must be different than qin= input Q-file name.");
    }
    let fp = match fopen(&qout, "w") {
        Some(fp) => fp,
        None => err!("qfile error: output Q-file did not open correctly."),
    };

    let mut formxy = "%.20g".to_owned();
    getparstring("formxy", &mut formxy);
    let formxylong = format!(",{}", formxy);
    let mut nshift = countparval("shiftrec") as usize;
    if nshift != countparval("shiftdist") as usize {
        err!("**** Error: shiftrec and shiftdist different lengths.");
    }
    let mut shiftrec = Vec::<i32>::new();
    let mut shiftdist = Vec::<f64>::new();
    if nshift > 0 {
        getparint("shiftrec", &mut shiftrec);
        getpardouble("shiftdist", &mut shiftdist);
        shiftrec.insert(0, 0);
        shiftrec.push(0);
        shiftdist.insert(0, 0.0);
        shiftdist.push(0.0);
        nshift += 2;
    }
    let mut q = if !has_qin {
        let zname = [
            "cdp", "cdpt", "fldr", "grnors", "grnofr", "grnlof", "gaps", "igi", "igc", "gx", "gy",
            "sx", "sy", "numa", "numb", "vuma", "vumb", "tupa", "offs", "tims", "tnmo", "dpth",
            "vels", "vnmo", "tupb", "vupa",
        ];
        let mut pname = Vec::new();
        let mut iztuple = 0;
        for j in 0..26 {
            if countparname(zname[j]) > 0 {
                pname.push(zname[j].to_owned());
                if j < 17 {
                    iztuple += 1;
                }
            }
        }
        match getviacommand(&mut pname, iztuple, 0) {
            Ok(q) => q,
            Err(errwarn) => match errwarn {
                1 => err!("getviacommand error: no non-tuple name passed in."),
                2 => err!("getviacommand error: non-tuple names have different amounts of values."),
                3 => err!("getviacommand error: independent dimension parameter is empty."),
                4 => err!("getviacommand error: an independent dimension parameter is empty."),
                5 => err!("getviacommand error: members of tuple have different amounts at same location."),
                _ => err!("getviacommand error: returned unrecognized error number = %d", errwarn),
            }
        }
    } else {
        let fpp = match fopen(&qin, "r") {
            Some(f) => f,
            None => err!("error: input Q-file did not open correctly."),
        };
        let mut ename = Vec::<String>::new();
        if countparval("ename") > 0 {
            getparstringarray("ename", &mut ename);
        }
        match getviaqfile(&fpp, &ename, -(ename.len() as i32), 0) {
            Ok(q) => q,
            Err(errwarn) => match errwarn {
                1 => err!("getqinfo error: extra C_SU_NAMES record in q-file"),
                2 => err!("getqinfo error: extra C_SU_NDIMS record in q-file"),
                3 => err!("getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),
                11 => err!("readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"),
                12 => err!("readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."),
                22 => err!("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),
                23 => err!("getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
                24 => err!("getviaqfile error: C_SU_NDIMS blank where valid number expected"),
                25 => err!("getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
                26 => err!("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
                27 => err!("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),
                n if n > 100 => err!("getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)", n-99),
                _ => err!("getviaqfile error: unrecognized error code %d", errwarn),
            }
        }
    };
    checkpars();
    if q.ifixd == 0 {
        err!("error: input with varying number of tuples is not allowed.");
    }
    let iztuple = q.iztuple as usize;
    let ktuple = q.ktuple as usize;
    let klast = iztuple + (ktuple - 1) * q.records[0].nto as usize;
    let mut ncdp = q.records.len();

    let mut locx = None;
    for i in 0..iztuple {
        if q.names[i] == keyx {
            locx = Some(i);
        }
    }
    if locx.is_none() {
        let ukey = format!("_{}_", keyx);
        for i in 0..iztuple {
            if q.names[i].contains(&ukey) {
                locx = Some(i);
            }
        }
    }
    let locx = locx.unwrap_or_else(|| {
        err!(
            "error: input q-file must have your keyx=%s (among non-tuple names).",
            keyx
        )
    });
    let locy;
    if keyy != "none" {
        let mut found = None;
        for i in 0..iztuple {
            if q.names[i] == keyy {
                found = Some(i);
            }
        }
        let ukey = format!("_{}_", keyy);
        for i in 0..iztuple {
            if q.names[i].contains(&ukey) {
                found = Some(i);
            }
        }
        locy = found.unwrap_or_else(|| {
            err!(
                "error: input q-file must have your keyy=%s (among non-tuple names).",
                keyy
            )
        });
    } else {
        locy = locx;
        chordi *= 2.0f64.sqrt();
        chordf *= 2.0f64.sqrt();
        dextra *= 2.0f64.sqrt();
        for i in 0..nshift {
            shiftdist[i] *= 2.0f64.sqrt();
        }
    }
    let mut locp = None;
    if keyp != "asis" {
        for i in 0..iztuple {
            if q.names[i] == keyp {
                locp = Some(i);
            }
        }
        let ukey = format!("_{}_", keyp);
        for i in 0..iztuple {
            if q.names[i].contains(&ukey) {
                locp = Some(i);
            }
        }
        let lp = locp.unwrap_or_else(|| {
            err!(
                "error: input q-file must have your keyp=%s (among non-tuple names).",
                keyp
            )
        });
        q.records.sort_by(|a, b| {
            if a.dlots[lp] < b.dlots[lp] {
                std::cmp::Ordering::Less
            } else if a.dlots[lp] > b.dlots[lp] {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        for jcdp in 1..ncdp {
            if q.records[jcdp - 1].dlots[lp] == q.records[jcdp].dlots[lp] {
                err!(
                    "error: Two sets of values input for location = %f",
                    q.records[jcdp].dlots[lp]
                );
            }
        }
    }

    // Initial chording, with extra points extrapolated at both ends.
    let mut dalla;
    let mut msize;
    if chordi > 0.0 {
        let mut dpatha = vec![0.0; ncdp];
        for jcdp in 1..ncdp {
            let pint = (q.records[jcdp].dlots[locx] - q.records[jcdp - 1].dlots[locx])
                * (q.records[jcdp].dlots[locx] - q.records[jcdp - 1].dlots[locx])
                + (q.records[jcdp].dlots[locy] - q.records[jcdp - 1].dlots[locy])
                    * (q.records[jcdp].dlots[locy] - q.records[jcdp - 1].dlots[locy]);
            dpatha[jcdp] = dpatha[jcdp - 1] + pint.sqrt();
        }
        let mut chordhere = 0.0 - chordi * nmaxa as f64 - dextra;
        let chordimax = dpatha[ncdp - 1] + chordi * nmaxa as f64 + dextra;
        let mout = (9.0 + (chordimax - chordhere) / chordi) as usize;
        dalla = vec![vec![0.0; klast]; mout];
        let mut nhere = 1;
        msize = 0;
        while chordhere <= chordimax {
            while chordhere > dpatha[nhere] && nhere < ncdp - 1 {
                nhere += 1;
            }
            let wi = if (dpatha[nhere] - dpatha[nhere - 1]).abs() < 1.0e-15 {
                0.5
            } else {
                (chordhere - dpatha[nhere - 1]) / (dpatha[nhere] - dpatha[nhere - 1])
            };
            for i in 0..klast {
                dalla[msize][i] =
                    q.records[nhere - 1].dlots[i] * (1.0 - wi) + q.records[nhere].dlots[i] * wi;
            }
            chordhere += chordi;
            msize += 1;
        }
    } else {
        dalla = vec![vec![0.0; klast]; ncdp + 2];
        for jcdp in 0..ncdp {
            for i in 0..klast {
                dalla[jcdp + 1][i] = q.records[jcdp].dlots[i];
            }
        }
        for i in 0..klast {
            dalla[0][i] = dalla[1][i] * 2.0 - dalla[2][i];
        }
        for i in 0..klast {
            dalla[ncdp + 1][i] = dalla[ncdp][i] * 2.0 - dalla[ncdp - 1][i];
        }
        msize = ncdp + 2;
    }

    if nmaxa > 0 {
        let ierr = runav(&mut dalla, msize, klast, nmina, nmaxa);
        if ierr == 1 {
            err!("**** Error in runav routine for averaging. Profile too short for averaging length.");
        } else if ierr > 0 {
            err!("**** Error in runav routine for averaging. Parameter out-of-range.");
        }
        let na = nmaxa as usize;
        for jcdp in 0..na {
            let wi = (na - jcdp + 1) as f64;
            for i in 0..klast {
                dalla[jcdp][i] = dalla[na][i] * wi + dalla[na + 1][i] * (1.0 - wi);
            }
        }
        for jcdp in msize - na - 1..msize {
            let wi = (nmaxa + jcdp as i32 - msize as i32 + 3) as f64;
            for i in 0..klast {
                dalla[jcdp][i] =
                    dalla[msize - na - 3][i] * (1.0 - wi) + dalla[msize - na - 2][i] * wi;
            }
        }
    }
    if nmaxs > 0 {
        let ierr = runav(&mut dalla, msize, klast, nmins, -nmaxs);
        if ierr == 1 {
            err!("**** Error in runav routine for pre-smoothing. Profile too short for smoothing length.");
        } else if ierr > 0 {
            err!("**** Error in runav routine for pre-smoothing. Parameter out-of-range.");
        }
        let ns = nmaxs as usize;
        for jcdp in 0..ns {
            let wi = (ns - jcdp + 1) as f64;
            for i in 0..klast {
                dalla[jcdp][i] = dalla[ns][i] * wi + dalla[ns + 1][i] * (1.0 - wi);
            }
        }
        for jcdp in msize - ns - 1..msize {
            let wi = (nmaxs + jcdp as i32 - msize as i32 + 3) as f64;
            for i in 0..klast {
                dalla[jcdp][i] =
                    dalla[msize - ns - 3][i] * (1.0 - wi) + dalla[msize - ns - 2][i] * wi;
            }
        }
        let ierr = runav(&mut dalla, msize, klast, nmins, nmaxs);
        if ierr == 1 {
            err!("**** Error in runav routine for smoothing. Profile too short for smoothing length.");
        } else if ierr > 0 {
            err!("**** Error in runav routine for smoothing. Parameter out-of-range.");
        }
        for jcdp in 0..ns {
            let wi = (ns - jcdp + 1) as f64;
            for i in 0..klast {
                dalla[jcdp][i] = dalla[ns][i] * wi + dalla[ns + 1][i] * (1.0 - wi);
            }
        }
        for jcdp in msize - ns - 2..msize {
            let wi = (nmaxs + jcdp as i32 - msize as i32 + 3) as f64;
            for i in 0..klast {
                dalla[jcdp][i] =
                    dalla[msize - ns - 3][i] * (1.0 - wi) + dalla[msize - ns - 2][i] * wi;
            }
        }
    }

    let mut dalls;
    if chordf > 0.0 {
        ncdp = msize;
        // Initial chording can leave one point. Provide an initialized boundary
        // distance for nhere=1, as in the minimally fixed C. The spare dalla
        // row is already zeroed; the output loop omits this single endpoint.
        let mut dpaths = vec![0.0; ncdp + 1];
        for jcdp in 1..ncdp {
            let pint = (dalla[jcdp][locx] - dalla[jcdp - 1][locx])
                * (dalla[jcdp][locx] - dalla[jcdp - 1][locx])
                + (dalla[jcdp][locy] - dalla[jcdp - 1][locy])
                    * (dalla[jcdp][locy] - dalla[jcdp - 1][locy]);
            dpaths[jcdp] = dpaths[jcdp - 1] + pint.sqrt();
        }
        dpaths[ncdp] = dpaths[ncdp - 1];
        let mut chordhere = 0.0 - chordf * nmaxs as f64;
        let chordfmax = dpaths[ncdp - 1] + chordf * nmaxs as f64;
        let mout = (9.0 + (chordfmax - chordhere) / chordf) as usize;
        dalls = vec![vec![0.0; klast]; mout];
        let mut nhere = 1;
        msize = 0;
        while chordhere <= chordfmax {
            while chordhere > dpaths[nhere] && nhere < ncdp - 1 {
                nhere += 1;
            }
            let wi = if (dpaths[nhere] - dpaths[nhere - 1]).abs() < 1.0e-15 {
                0.5
            } else {
                (chordhere - dpaths[nhere - 1]) / (dpaths[nhere] - dpaths[nhere - 1])
            };
            for i in 0..klast {
                dalls[msize][i] = dalla[nhere - 1][i] * (1.0 - wi) + dalla[nhere][i] * wi;
            }
            chordhere += chordf;
            msize += 1;
        }
    } else {
        dalls = dalla;
    }

    if nshift > 0 {
        shiftrec[0] = 0;
        shiftrec[nshift - 1] = msize as i32 - 1;
        shiftdist[0] = shiftdist[1];
        shiftdist[nshift - 1] = shiftdist[nshift - 2];
        dalla = dalls;
        ncdp = msize;
        let mut dpaths = vec![0.0; ncdp + 1];
        let mut pint = 0.0;
        for jcdp in 1..ncdp {
            pint = (dalla[jcdp][locx] - dalla[jcdp - 1][locx])
                * (dalla[jcdp][locx] - dalla[jcdp - 1][locx])
                + (dalla[jcdp][locy] - dalla[jcdp - 1][locy])
                    * (dalla[jcdp][locy] - dalla[jcdp - 1][locy]);
            dpaths[jcdp] = dpaths[jcdp - 1] + pint.sqrt();
        }
        dpaths[ncdp] = dpaths[ncdp - 1] + pint.sqrt();
        dalls = vec![vec![0.0; klast]; 9 + ncdp];
        let mut nhere = 1;
        pint = 0.0;
        let mut pint2 = shiftdist[0];
        for jcdp in 1..ncdp {
            if jcdp as i32 > shiftrec[nhere] {
                nhere += 1;
                pint = (shiftdist[nhere] - shiftdist[nhere - 1])
                    / (dpaths[shiftrec[nhere] as usize] - dpaths[shiftrec[nhere - 1] as usize]);
                pint2 = shiftdist[nhere - 1];
            }
            let wi = if (dpaths[jcdp] - dpaths[jcdp - 1]).abs() < 1.0e-15 {
                0.5
            } else {
                pint2 / (dpaths[jcdp] - dpaths[jcdp - 1])
            };
            for i in 0..klast {
                dalls[jcdp - 1][i] = dalla[jcdp - 1][i] * (1.0 - wi) + dalla[jcdp][i] * wi;
            }
            pint2 += pint * (dpaths[jcdp] - dpaths[jcdp - 1]);
        }
    }

    fprintf!(fp, "C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
    for _i in 0..klast + 1 {
        fprintf!(fp, ",%s", formxy);
    }
    if q.ifixd != 2 {
        fprintf!(fp, "\nC_SU_NDIMS,%s", q.ndims[0]);
        for _i in 1..iztuple + 1 {
            fprintf!(fp, ",");
        }
        for i in 0..q.records[0].nto as usize {
            for _k in 0..ktuple - 1 {
                // An empty numeric parameter leaves C's allocated independent
                // value uninitialized. Fresh allocations in the reference
                // binary yield zero; the reader returns an empty Vec.
                let independent = q.pindepa.get(i).copied().unwrap_or(0.0);
                fprintf!(fp, ",%.15g", independent);
            }
        }
    }
    fprintf!(fp, "\nC_SU_NAMES\nC_SU_ID");
    fprintf!(fp, ",cdp");
    for i in 0..iztuple {
        if q.names[i] == "cdp" {
            fprintf!(fp, ",null");
        } else {
            fprintf!(fp, ",%s", q.names[i]);
        }
    }
    for _i in 0..q.records[0].nto {
        for k in 0..ktuple - 1 {
            if q.names[k + iztuple] == "cdp" {
                fprintf!(fp, ",null");
            } else {
                fprintf!(fp, ",%s", q.names[k + iztuple]);
            }
        }
    }
    fprintf!(fp, "\n");
    // C's lrint returns a long: subtract in i64 before narrowing to int.
    // Converting directly to i32 would saturate large profile origins.
    let jcrf = if let Some(lp) = locp {
        (dalls[1][lp] * pcru + pcrz).round_ties_even() as i64
    } else {
        pcrz.round_ties_even() as i64
    }
    .wrapping_sub(1) as i32;
    for jcdp in 1..msize - 1 {
        fprintf!(fp, "Q");
        fprintf!(fp, &formxylong, (jcdp as i32 + jcrf) as f64);
        for i in 0..klast {
            fprintf!(fp, &formxylong, dalls[jcdp][i]);
        }
        fprintf!(fp, "\n");
    }
    exit(CWP_Exit());
}

// Running average from suprofcsv.c. Negative nmax performs pre-smoothing.
// End values are scratch space, and the first total is deliberately advanced
// before being stored, exactly as in the original routine.
fn runav(dall: &mut [Vec<f64>], ncdp: usize, klast: usize, nmin: i32, mut nmax: i32) -> i32 {
    let mut iavr = 1;
    if nmax < 0 {
        nmax = -nmax;
        iavr = 0;
    }
    if ncdp < 1 || klast < 1 || nmin < 0 || nmax < nmin {
        return 2;
    }
    if 2 * (nmax + 3) as usize > ncdp {
        return 1;
    }
    let nmax = nmax as usize;
    let nmin = nmin as usize;
    let mut dtog = vec![0.0; klast];
    let mut n = 0;
    for i in 0..nmax * 2 + 1 {
        n += 1;
        for k in 0..klast {
            dtog[k] += dall[i][k];
        }
    }
    for i in nmax - nmin + 1..nmax + nmin {
        n -= 1;
        for k in 0..klast {
            dtog[k] -= dall[i][k];
        }
    }
    for i in nmax..ncdp - nmax - 1 {
        for k in 0..klast {
            dtog[k] += dall[i + 1 + nmax][k] - dall[i - nmax][k];
        }
        if nmin > 0 {
            for k in 0..klast {
                dtog[k] += dall[i + 1 - nmin][k] - dall[i + nmin][k];
            }
        }
        if iavr == 1 {
            for k in 0..klast {
                dall[i - nmax][k] = dtog[k] / n as f64;
            }
        } else {
            for k in 0..klast {
                dall[i - nmax][k] = dall[i + 1][k] + dall[i + 1][k] - dtog[k] / n as f64;
            }
        }
    }
    for i in (nmax..=ncdp - nmax - 1).rev() {
        for k in 0..klast {
            dall[i][k] = dall[i - nmax][k];
        }
    }
    0
}
