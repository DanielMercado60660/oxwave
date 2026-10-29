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

//! SUSMOQCSV - Input Q-records, alter and spatially smooth, output Q-records.
//! Port of su/main/headers/susmoqcsv.c (Revision 1.00, February 2024).
//! Created February 2024: Andre Latour.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! The running-average routine was copied from sufarldcsv.
//! Copyright (c) Colorado School of Mines, 2023. All rights reserved.

use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;
use oxwave_linear_interp::*;
use oxwave_q_file::*;
use std::ffi::c_void;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUSMOQCSV - Input Q-records, alter and spatially smooth, output Q-records. ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  susmoqcsv [parameters].   (No traces in or out).                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " qin=       Input file containing q-records. Must be specified.             ",
    "            For 3D there must be the same number of q-records in the        ",
    "            qin file as cdps (cells) defined by the 3D Grid.                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qout=      Output file for q-records. Must be specified.                   ",
    "            Output records will be in cdp order.                            ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rfile=     If set, read a K-record containing 3D grid definition.          ",
    "            Assumed to be 2D if no file is specified here.                  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " check=0   Do not print checking details. The options of this program       ",
    "           cannot repair a corrupt input grid. Grid checking is always      ",
    "           done on input grid definition and will error-halt if corrupt.    ",
    "      =1   print checking details.                                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " inloc=1   Location options. The option here determines which input names   ",
    "           contain the locations of all other input values. The option here ",
    "           must find the corresponding names within the non-tuple names in  ",
    "           the input q-file.                                                ",
    "      =1   cdp (default). Option 1 is required if no grid is input.         ",
    "      =2   cdpt (grid transposed cdp numbers).                              ",
    "      =3   igi,igc (grid index numbers).                                    ",
    "      =4   gx,gy (grid coordinates).  Note option 6.                        ",
    "      =5   sx,sy (world coordinates). Note option 7.                        ",
    "      =6   gx,gy (world coordinates).                                       ",
    "      =7   sx,sy (grid coordinates).                                        ",
    "     Note: The q-records in qin file do not have to be in order but there   ",
    "           must be exactly one record for each cdp defined by the 3D Grid.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keyloc=   Location key. Can only be specified if inloc=1 and 2D (no rfile).",
    "           Default is cdp. Name must be in the input Q-file.                ",
    "           This allows smoothing based on values such as point numbers,     ",
    "           station numbers, reciever numbers, or shot numbers.              ",
    "           Note that q-records are internally sorted by the values of       ",
    "           this name. Only cdps are required to be contiguous with an       ",
    "           increment of 1 (no gaps). For other names navrg, nsmth, nback    ",
    "           refer to the sorted record order, not to their actual values.    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " extrapt=0        do not extrapolate at ends in independent values.         ",
    "                  (values outside range of functions are held constant).    ",
    "               =1 extrapolate outside range of functions                    ",
    "               =2 extrapolate when less than range of functions             ",
    "               =3 extrapolate when greater than range of functions          ",
    "                                                                            ",
    " outind=0   Output independent dimension range. One, two, or three numbers. ",
    "            No default, but outlist can be specified instead.               ",
    "            Cannot be specified if input has no tuples.                     ",
    "       =f,l,i start at f, increment by i, end at last value <= l.           ",
    "       =f,l   start at f, increment by 1000, end at last value <= l.        ",
    "       =f     start at f, increment by 1000, end at last value <= 3000.     ",
    "            The i value cannot be 0, but can be negative, which requires    ",
    "            an l value less than or equal to f.                             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " outlist=   Output independent dimension value list. Cannot be specified    ",
    "            if outind is used. Cannot be specified if input has no tuples.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " ename=     Names to eliminate. These names will be ignored in the input    ",
    "            q-records and therefore not be in the output q-records (note    ",
    "            that standard 2d and 3d names will still be output even if you  ",
    "            specify them here). A primary purpose of this option is to get  ",
    "            rid of names in input tuples in order to reduce the complexity  ",
    "            of output q-records (in case certain programs cannot handle it).",
    "            But you can also get rid of non-tuple names here.               ",
    "\t\t\t\t\t\t\t\t\t     ",
    " sname=     Names to smooth. By default all names with 1 element per record ",
    "            are NOT averaged and smoothed. Listing names here overrides     ",
    "            that and averages and smooths their values.                     ",
    "                                                                            ",
    " cname=     Names to just copy. By default all names with multiple values   ",
    "            per record (tuples) are averaged and smoothed. Listing them here",
    "            overrides that and does NOT average and smooth their values.    ",
    "                                                                            ",
    "                                                                            ",
    " navrg=10       CDP Averaging Operator Size.                                ",
    "                The averaging is symmetrical and extends navrg on each side ",
    "                of the location being averaged.                             ",
    "                For 2D surveys, this is number of cin Q-records to average  ",
    "                symmetrically with its neighbours (navrg on each side).     ",
    "                For 3D surveys, two values can be listed. The first value   ",
    "                is for the inline direction and the second is for the       ",
    "                crossline direction. If no second value is specified, the   ",
    "                same value is used for both directions.                     ",
    "            =0  means no averaging is done.                                 ",
    "            =-1 means average the values of all locations in this direction.",
    "                You also MUST set nsmth to 0 for this direction.            ",
    "                Typically this option is used for 3D surveys where the      ",
    "                number of inline or crossline cells is less than 5 or one   ",
    "                of the directions has little change (or you want to produce ",
    "                a constant value in that direction).                        ",
    "           ***  You must use 0 or -1 if number of locations is less than 5. ",
    "\t\t\t\t\t\t\t\t\t     ",
    " nsmth=navrg    CDP Smoothing Operator Size.                                ",
    "                This smoothing uses a technique which leaves the values     ",
    "                approximately centred with their initial values.            ",
    "                The smoothing affect is symmetrical and extends 2*nsmth on  ",
    "                each side of the location being smoothed.                   ",
    "                For 3D surveys, two values can be listed. The first value   ",
    "                is for the inline direction and the second is for the       ",
    "                crossline direction. If no second value is specified, the   ",
    "                same value is used for both directions.                     ",
    "             =0 means no smoothing is done.                                 ",
    "           ***  You must use 0 if number of locations is <5 or navrg=-1     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " nback=3        Linear Extrapolation Location.                              ",
    "                For 2D surveys, this is the number of cin Q-records inwards ",
    "                from the ends of the survey used to linearly extrapolate    ",
    "                some padding outside the ends of the survey. These          ",
    "                extrapolated values affect the averaging and smoothing of   ",
    "                the navrg and nsmth parameters (at both ends of the survey).",
    "                This value determines which cin records inside the survey   ",
    "                are used to determine extrapolation points.                 ",
    "                Must be less than the number of Q-records in the cin file.  ",
    "            =1  means use the first and first+1 records to extrapolate into ",
    "                the padding at the beginning of the survey -and- use the    ",
    "                last and last-1 records to extrapolate into the end padding.",
    "            =n  means use the first and first+n records to extrapolate into ",
    "                the padding at the beginning of the survey -and- use the    ",
    "                last and last-n records to extrapolate into the end padding.",
    "            =0  means duplicate the values of first record into the padding ",
    "                at the beginning of the survey -and- duplicate the values   ",
    "                of the last record into the padding at the end of survey.   ",
    "           ***  For 3D surveys, two values can be listed. The first value   ",
    "                indicates how many cells inwards in the inline direction    ",
    "                and the second value indicates how many cells inwards in    ",
    "                the crossline direction. If no second value is specified,   ",
    "                the same value is used for both directions.                 ",
    "                Must be less than the number of 3D Grid cells defined       ",
    "                in the same direction.                                      ",
    "         Note:  The AMOUNT of added padding is NOT set by this parameter,   ",
    "                it is actually set as needed for navrg and nsmth parameters.",
    "\t\t\t\t\t\t\t\t\t     ",
    " cdpt=0      Output in transposed cdp order.                                ",
    "          =0 No. Output in cdp order.                                       ",
    "          =1 Yes. Output in transposed cdp order. Can only specify for 3D.  ",
    "             This has no effect on values within the q-records.             ",
    "             This order may be convenient (notably 2D by 2D migrations).    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " formtv=%.2f The C format code for printing all values in the q-records     ",
    "             except the standard names cdp,cdpt,igi,igc.                    ",
    "             These names are specially handled within this program.         ",
    "             For 3D, they are generated for output using grid definition.   ",
    "             For 2D, igi,igc are not output.                                ",
    "             The default is 2 digits below the decimal point.               ",
    "\t\t\t\t\t\t\t\t\t     ",
    " formxy=%.2f The C format code for printing gx,gy,sx,sy in the q-records.   ",
    "             These names are specially handled within this program.         ",
    "             For 3D, they are generated for output using grid definition.   ",
    "             For 2D, they are not output.                                   ",
    "             The default is 2 digits below the decimal point.               ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "   ------------------------------------------------------------------       ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  Notes:                                                       \t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   All input values are stored in double-precision (8 byte float).          ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   This program does not care if the values are seconds or milliseconds     ",
    "   or feet or metres. To this program they are all just numbers.            ",
    "\t\t\t\t\t\t\t\t\t     ",
];

// libc qsort is also used by the C program. Its ordering of equal
// non-cdp keyloc values is observable in the copied and smoothed values.
unsafe extern "C" {
    fn qsort(
        base: *mut c_void,
        count: usize,
        size: usize,
        compare: unsafe extern "C" fn(*const c_void, *const c_void) -> i32,
    );
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    if !isatty(0) || !isatty(1) {
        err!("**** Error: this program does not input or output traces.");
    }

    let mut i = 0;
    gridcommand(&mut i);
    if i != 0 {
        err!("error: grid definition parameters on command line (input grid only via K-file).");
    }
    let mut rname = String::new();
    let is3d = getparstring("rfile", &mut rname) != 0;
    let mut mgtextr = 0;
    getparint("extrapt", &mut mgtextr);
    if mgtextr < 0 || mgtextr > 3 {
        err!("error: extrapt= option not in range ");
    }
    let mut inloc = 1;
    getparint("inloc", &mut inloc);
    if inloc < 1 || inloc > 7 {
        err!("error: inloc= option not in range ");
    }
    if !is3d && inloc > 1 {
        err!("**** Error: cannot specify inloc>1 with no input 3d grid.");
    }
    let mut keyloc = "cdp".to_owned();
    if countparval("keyloc") > 0 {
        if inloc != 1 || is3d {
            err!("Error: to use keyloc= must be inloc=1 and also not a 3d.");
        }
        getparstring("keyloc", &mut keyloc);
    }
    let mut icheck = 0;
    getparint("check", &mut icheck);
    if icheck < 0 || icheck > 1 {
        err!("**** Error: check= value out of range.");
    }
    if !is3d && icheck != 0 {
        err!("**** Error: cannot specify icheck= with no input 3d grid.");
    }
    let mut pname_file = String::new();
    let mut qname_file = String::new();
    let has_p = getparstring("qin", &mut pname_file) != 0;
    let has_q = getparstring("qout", &mut qname_file) != 0;
    if has_p && has_q && pname_file == qname_file {
        err!("**** Error: qout= output Q-file name must be different than qin= input Q-file name.");
    }
    let fpq = fopen(&qname_file, "w").unwrap_or_else(|| {
        err!("qfile error: output Q-file did not open correctly.");
    });

    let mut numdind = 0;
    let mut dind = Vec::<f64>::new();
    if countparval("outind") > 0 {
        if countparval("outlist") > 0 {
            err!("error: outlist= and outind= cannot both be specified.");
        }
        if countparval("outind") > 3 {
            err!("error: outind= too long (maximum 3 values: first,last,inc).");
        }
        let mut outind = [0.0, 3000.0, 1000.0];
        getpardouble("outind", &mut outind[..]);
        if outind[2] == 0.0 || (outind[1] - outind[0]) / outind[2] < -0.000005 {
            err!("error: outind= values are wrong.");
        }
        numdind = (1.00001 + (outind[1] - outind[0]) / outind[2]) as i32;
        dind = vec![0.0; numdind as usize];
        for i in 0..numdind as usize {
            dind[i] = outind[0] + i as f64 * outind[2];
        }
    } else if countparval("outlist") > 0 {
        numdind = countparval("outlist");
        dind = vec![0.0; numdind as usize];
        getpardouble("outlist", &mut dind);
    }

    let mut navrg = [10i32; 2];
    if countparval("navrg") > 2 {
        err!("**** Error: navrg= list can only have 1 or 2 values ");
    }
    if countparval("navrg") >= 1 {
        getparint("navrg", &mut navrg[..]);
        if countparval("navrg") == 1 {
            navrg[1] = navrg[0];
        }
    }
    if navrg[0] < -1 {
        err!("**** Error in first navrg value. It must be greater or equal to -1 ");
    }
    if navrg[1] < -1 {
        err!("**** Error in second navrg value. It must be greater or equal to -1 ");
    }

    let mut nsmth = navrg;
    if countparval("nsmth") > 2 {
        err!("**** Error: nsmth= list can only have 1 or 2 values ");
    }
    if countparval("nsmth") < 1 {
        if nsmth[0] == -1 {
            nsmth[0] = 0;
        }
        if nsmth[1] == -1 {
            nsmth[1] = 0;
        }
    } else {
        getparint("nsmth", &mut nsmth[..]);
        if countparval("nsmth") == 1 {
            nsmth[1] = nsmth[0];
        }
    }
    if nsmth[0] < 0 {
        err!("**** Error in first nsmth value. It must be greater or equal to 0 ");
    }
    if nsmth[1] < 0 {
        err!("**** Error in second nsmth value. It must be greater or equal to 0 ");
    }

    let mut nback = [3i32; 2];
    if countparval("nback") > 2 {
        err!("**** Error: nback= list can only have 1 or 2 values ");
    }
    if countparval("nback") >= 1 {
        getparint("nback", &mut nback[..]);
        if countparval("nback") == 1 {
            nback[1] = nback[0];
        }
    }
    if nback[0] < 0 {
        err!("**** Error: first nback= must be greater than or equal to 0 ");
    }
    if nback[1] < 0 {
        err!("**** Error: second nback= must be greater than or equal to 0 ");
    }

    let mut cdpt = 0;
    getparint("cdpt", &mut cdpt);
    if cdpt < 0 || cdpt > 1 {
        err!("**** Error: cdpt parameter out of range.");
    }
    if !is3d && cdpt == 1 {
        err!("**** Error: cannot specify cdpt=1 with no input 3d grid.");
    }

    let mut formxy = "%.2f".to_owned();
    if getparstring("formxy", &mut formxy) != 0 && !is3d {
        err!("**** Error: cannot specify formxy= with no input 3d grid.");
    }
    let formxylong = format!("Q,%d,%d,%d,%d,{formxy},{formxy},{formxy},{formxy}");
    // Standard names are generated explicitly in the output, rather than
    // copied from input. Their formats also go into the Q-file header.
    let mut qname = vec![keyloc.clone(), "cdpt".to_owned()];
    let mut qform = vec!["%.20g".to_owned(); 2];
    if is3d {
        for name in ["igi", "igc", "gx", "gy", "sx", "sy"] {
            qname.push(name.to_owned());
        }
        qform.push("%.20g".to_owned());
        qform.push("%.20g".to_owned());
        for _ in 0..4 {
            qform.push(formxy.clone());
        }
    }
    let numstandard = qname.len();
    let mut formtv = "%.2f".to_owned();
    getparstring("formtv", &mut formtv);
    let formtvlong = format!(",{formtv}");

    let mut gvals = [0.0; 99];
    if is3d {
        let fpr = fopen(&rname, "r").unwrap_or_else(|| {
            err!("error: input K-file did not open correctly.");
        });
        let mut errwarn = 1;
        gridread(Some(&fpr), &mut gvals, &mut errwarn);
        if errwarn > 0 {
            err!("error reading grid from K-file.");
        }
        gridset(&mut gvals, &mut errwarn);
        if errwarn == 1 {
            err!("gridset error: grid_wb cell width must be positive.");
        } else if errwarn == 2 {
            err!("gridset error: grid_wc cell width must be positive.");
        } else if errwarn == 3 {
            err!("gridset error: corner B is within grid_wb cell width of corner A.");
        } else if errwarn > 0 {
            err!("gridset error: returned with some unrecognized error code.");
        } else if errwarn == -1 {
            warn!("gridset warning: corner C is near A and is reset to A.");
        }
        gridcheck(&gvals, icheck, &mut errwarn);
        if errwarn > 0 {
            err!("gridcheck error: returned with some unrecognized error code.");
        }
    }

    let fpp = fopen(&pname_file, "r").unwrap_or_else(|| {
        err!("error: input Q-file did not open correctly.");
    });
    let mut pname = Vec::<String>::new();
    if countparval("ename") > 0 {
        getparstringarray("ename", &mut pname);
    }
    let mut sname = Vec::<String>::new();
    if countparval("sname") > 0 {
        getparstringarray("sname", &mut sname);
    }
    let mut cname = Vec::<String>::new();
    if countparval("cname") > 0 {
        getparstringarray("cname", &mut cname);
    }
    if is3d {
        if inloc != 1 {
            pname.push("cdp".to_owned());
        }
        if inloc != 2 {
            pname.push("cdpt".to_owned());
        }
        if inloc != 3 {
            pname.push("igi".to_owned());
            pname.push("igc".to_owned());
        }
        if inloc != 4 && inloc != 6 {
            pname.push("gx".to_owned());
            pname.push("gy".to_owned());
        }
        if inloc != 5 && inloc != 7 {
            pname.push("sx".to_owned());
            pname.push("sy".to_owned());
        }
    }
    // Nonpositive numpname means exclude names on the list.
    let qfile = match getviaqfile(&fpp, &pname, -(pname.len() as i32), numdind) {
        Ok(q) => q,
        Err(e) => match e {
            1 => err!("getqinfo error: extra C_SU_NAMES record in q-file"),
            2 => err!("getqinfo error: extra C_SU_NDIMS record in q-file"),
            3 => err!(
                "getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."
            ),
            11 => err!(
                "readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"
            ),
            12 => err!(
                "readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."
            ),
            22 => {
                err!("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.")
            }
            23 => err!("getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
            24 => err!("getviaqfile error: C_SU_NDIMS blank where valid number expected"),
            25 => err!("getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
            26 => err!("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
            27 => {
                err!("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.")
            }
            n if n > 100 => err!(
                "getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",
                n - 99
            ),
            _ => err!("getviaqfile error: unrecognized error code %d", e),
        },
    };
    pname = qfile.names;
    let iztuple = qfile.iztuple as usize;
    let ktuple = qfile.ktuple as usize;
    let ifixd = qfile.ifixd;
    let mut recinfo = qfile.records;
    let pindepa = qfile.pindepa;
    let ndims = qfile.ndims;
    let ncdp = recinfo.len();
    if ifixd == 0 {
        err!("error: input q-records have varying number of values.");
    }
    if numdind > 0 && ifixd == 2 {
        err!("error: input does not have tuples so outind= and outlist= cannot be specified.");
    }
    let mut lastab = 0usize;
    let mut lastac = 0usize;
    if is3d {
        lastab = (0.1 + gvals[12]) as usize;
        lastac = (0.1 + gvals[13]) as usize;
        if lastab * lastac != ncdp {
            err!("error: input q-file does not have same amount of q-records as 3D Grid cdps (cells).");
        }
    }
    if ifixd == 1 {
        for i in 1..recinfo[0].nto as usize {
            if pindepa[i - 1] >= pindepa[i] {
                err!("error: independent dimension values not increasing.");
            }
        }
        if numdind == 0 {
            numdind = recinfo[0].nto;
            dind = vec![0.0; numdind as usize];
            for i in 0..numdind as usize {
                dind[i] = pindepa[i];
            }
        }
    }

    for icdp in 0..ncdp {
        recinfo[icdp].kinf = vec![0; 3];
    }
    let mut jnloc1 = None;
    let mut jnloc2 = None;
    if inloc == 1 {
        for i in 0..iztuple {
            if pname[i] == keyloc {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!(
                "error: for inloc=1, input must have your keyloc=%s (among non-tuple names).",
                keyloc
            );
        }
        for icdp in 0..ncdp {
            let k = recinfo[icdp].dlots[jnloc1.unwrap()].round_ties_even() as i64 as i32;
            recinfo[icdp].kinf[0] = k;
            if is3d {
                let (mut igi, mut igc) = (0, 0);
                gridcdpic(&gvals, k, &mut igi, &mut igc);
                recinfo[icdp].kinf[1] = igi;
                recinfo[icdp].kinf[2] = igc;
                if igi < -2147483644 {
                    err!("error: for inloc=1, input cdp %d is not in grid", k);
                }
            } else {
                recinfo[icdp].kinf[1] = k;
                recinfo[icdp].kinf[2] = 1;
            }
        }
    } else if inloc == 2 {
        for i in 0..iztuple {
            if pname[i] == "cdpt" {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!("error: for inloc=2, input must have cdpt (among non-tuple names).");
        }
        for icdp in 0..ncdp {
            let k = recinfo[icdp].dlots[jnloc1.unwrap()].round_ties_even() as i64 as i32;
            let (mut igi, mut igc, mut cdp) = (0, 0, 0);
            gridcdpic90(&gvals, k, &mut igi, &mut igc);
            if igi < -2147483644 {
                err!(
                    "error: for inloc=2, input 90 degree cdpt %d is not in grid",
                    k
                );
            }
            gridiccdp(&gvals, igi, igc, &mut cdp);
            recinfo[icdp].kinf = vec![cdp, igi, igc];
        }
    } else if inloc == 3 {
        for i in 0..iztuple {
            if pname[i] == "igi" {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!("error: for inloc=3, input must have igi (among non-tuple names).");
        }
        for i in 0..iztuple {
            if pname[i] == "igc" {
                jnloc2 = Some(i);
                break;
            }
        }
        if jnloc2.is_none() {
            err!("error: for inloc=3, input must have igc (among non-tuple names).");
        }
        for icdp in 0..ncdp {
            let mut cdp = 0;
            let igi = recinfo[icdp].dlots[jnloc1.unwrap()].round_ties_even() as i64 as i32;
            let igc = recinfo[icdp].dlots[jnloc2.unwrap()].round_ties_even() as i64 as i32;
            gridiccdp(&gvals, igi, igc, &mut cdp);
            if cdp < -2147483644 {
                err!(
                    "error: for inloc=3 input igi,igc = %d,%d are not in grid",
                    igi,
                    igc
                );
            }
            recinfo[icdp].kinf = vec![cdp, igi, igc];
        }
    } else if inloc == 4 {
        for i in 0..iztuple {
            if pname[i] == "gx" {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!("error: for inloc=4, input must have gx (among non-tuple names).");
        }
        for i in 0..iztuple {
            if pname[i] == "gy" {
                jnloc2 = Some(i);
                break;
            }
        }
        if jnloc2.is_none() {
            err!("error: for inloc=4, input must have gy (among non-tuple names).");
        }
        for icdp in 0..ncdp {
            let (mut igi, mut igc, mut cdp) = (0, 0, 0);
            let xg = recinfo[icdp].dlots[jnloc1.unwrap()];
            let yg = recinfo[icdp].dlots[jnloc2.unwrap()];
            let (mut xw, mut yw) = (0.0, 0.0);
            gridgridxyrawxy(&gvals, xg, yg, &mut xw, &mut yw);
            gridrawxycdpic(&gvals, xw, yw, &mut cdp, &mut igi, &mut igc);
            if cdp < -2147483644 {
                err!(
                    "error: for inloc=4, input gx,gy = %g,%g are not in grid",
                    xg,
                    yg
                );
            }
            recinfo[icdp].kinf = vec![cdp, igi, igc];
        }
    } else if inloc == 5 {
        for i in 0..iztuple {
            if pname[i] == "sx" {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!("error: for inloc=5, input must have sx (among non-tuple names).");
        }
        for i in 0..iztuple {
            if pname[i] == "sy" {
                jnloc2 = Some(i);
                break;
            }
        }
        if jnloc2.is_none() {
            err!("error: for inloc=5, input must have sy (among non-tuple names).");
        }
        for icdp in 0..ncdp {
            let (mut igi, mut igc, mut cdp) = (0, 0, 0);
            let xw = recinfo[icdp].dlots[jnloc1.unwrap()];
            let yw = recinfo[icdp].dlots[jnloc2.unwrap()];
            gridrawxycdpic(&gvals, xw, yw, &mut cdp, &mut igi, &mut igc);
            if cdp < -2147483644 {
                err!(
                    "error: for inloc=5, input sx,sy = %g,%g are not in grid",
                    xw,
                    yw
                );
            }
            recinfo[icdp].kinf = vec![cdp, igi, igc];
        }
    } else if inloc == 6 {
        for i in 0..iztuple {
            if pname[i] == "gx" {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!("error: for inloc=6, input must have gx (among non-tuple names).");
        }
        for i in 0..iztuple {
            if pname[i] == "gy" {
                jnloc2 = Some(i);
                break;
            }
        }
        if jnloc2.is_none() {
            err!("error: for inloc=6, input must have gy (among non-tuple names).");
        }
        for icdp in 0..ncdp {
            let (mut igi, mut igc, mut cdp) = (0, 0, 0);
            let xw = recinfo[icdp].dlots[jnloc1.unwrap()];
            let yw = recinfo[icdp].dlots[jnloc2.unwrap()];
            gridrawxycdpic(&gvals, xw, yw, &mut cdp, &mut igi, &mut igc);
            if cdp < -2147483644 {
                err!(
                    "error: for inloc=6, input gx,gy = %g,%g are not in grid",
                    xw,
                    yw
                );
            }
            recinfo[icdp].kinf = vec![cdp, igi, igc];
        }
    } else if inloc == 7 {
        for i in 0..iztuple {
            if pname[i] == "sx" {
                jnloc1 = Some(i);
                break;
            }
        }
        if jnloc1.is_none() {
            err!("error: for inloc=7, input must have sx (among non-tuple names).");
        }
        for i in 0..iztuple {
            if pname[i] == "sy" {
                jnloc2 = Some(i);
                break;
            }
        }
        if jnloc2.is_none() {
            err!("error: for inloc=7, input must have sy (among non-tuple names).");
        }
        for icdp in 0..ncdp {
            let (mut igi, mut igc, mut cdp) = (0, 0, 0);
            let xg = recinfo[icdp].dlots[jnloc1.unwrap()];
            let yg = recinfo[icdp].dlots[jnloc2.unwrap()];
            let (mut xw, mut yw) = (0.0, 0.0);
            gridgridxyrawxy(&gvals, xg, yg, &mut xw, &mut yw);
            gridrawxycdpic(&gvals, xw, yw, &mut cdp, &mut igi, &mut igc);
            if cdp < -2147483644 {
                err!(
                    "error: for inloc=7, input sx,sy = %g,%g are not in grid",
                    xg,
                    yg
                );
            }
            recinfo[icdp].kinf = vec![cdp, igi, igc];
        }
    }

    // SAFETY: all records are initialized and have three kinf entries.
    // qsort only relocates whole QInfo values; the callback cannot unwind.
    unsafe {
        qsort(
            recinfo.as_mut_ptr().cast(),
            ncdp,
            std::mem::size_of::<QInfo>(),
            comp_sort,
        );
    }
    if is3d {
        for jcdp in 1..ncdp {
            if recinfo[jcdp - 1].kinf[0] == recinfo[jcdp].kinf[0] {
                err!(
                    "error: Two q-records input for cdp,igi,igc = %d %d %d",
                    recinfo[jcdp].kinf[0],
                    recinfo[jcdp].kinf[1],
                    recinfo[jcdp].kinf[2]
                );
            }
        }
    } else if keyloc == "cdp" {
        for jcdp in 1..ncdp {
            if recinfo[jcdp - 1].kinf[0].wrapping_add(1) != recinfo[jcdp].kinf[0] {
                err!(
                    "error: input q-records have missing cdp(s) after cdp= %d",
                    recinfo[jcdp - 1].kinf[0]
                );
            }
        }
    }

    checkpars();
    let mut kindx = vec![0; iztuple + ktuple - 1];
    let mut ksize = vec![0; iztuple + ktuple - 1];
    qelementout(
        iztuple as i32,
        ktuple as i32,
        numdind,
        &mut kindx,
        &mut ksize,
    );
    let klast = (kindx[iztuple + ktuple - 2] + ksize[iztuple + ktuple - 2]) as usize;
    let mut dswap = vec![0.0; klast];
    let mut dtog = vec![0.0; klast];
    let mut ksmth = vec![0; klast];
    for jcdp in 0..ncdp {
        for j in 0..iztuple {
            dswap[kindx[j] as usize] = recinfo[jcdp].dlots[j];
        }
        for i in 0..numdind as usize {
            for k in 0..ktuple - 1 {
                let nto = recinfo[jcdp].nto as usize;
                linterpd(
                    dind[i],
                    &pindepa[..nto],
                    &recinfo[jcdp].dlots[iztuple + k * nto..iztuple + (k + 1) * nto],
                    mgtextr,
                    &mut dswap[kindx[k + iztuple] as usize + i],
                );
            }
        }
        for j in 0..klast {
            recinfo[jcdp].dlots[j] = dswap[j];
        }
    }

    // Negative ksize flags standard names so they are not output twice.
    for j in 0..iztuple + ktuple - 1 {
        for i in 0..numstandard {
            if qname[i] == pname[j] {
                ksize[j] = -ksize[j];
                break;
            }
        }
    }
    // Smooth all columns, but only copy selected columns back to RecInfo.
    // Single values default to copying; tuple values default to smoothing.
    for k in 0..iztuple {
        ksmth[k] = 0;
    }
    for i in 0..iztuple {
        if ksize[i] > 0 {
            for j in 0..sname.len() {
                if sname[j] == pname[i] {
                    ksmth[i] = 1;
                }
            }
        }
    }
    for k in iztuple..klast {
        ksmth[k] = 1;
    }
    for k in 0..ktuple - 1 {
        let i = ktuple - 2 - k;
        if ksize[k + iztuple] > 0 {
            for j in 0..cname.len() {
                if cname[j] == pname[i + iztuple] {
                    for n in 0..ksize[iztuple + k] as usize {
                        ksmth[kindx[iztuple + k] as usize + n] = 0;
                    }
                }
            }
        }
    }
    let mut isize = 2 * navrg[0] + 3;
    if is3d && 2 * navrg[1] + 3 > isize {
        isize = 2 * navrg[1] + 3;
    }
    if 2 * nsmth[0] + 1 > isize {
        isize = 2 * nsmth[0] + 1;
    }
    if is3d && 2 * nsmth[1] + 1 > isize {
        isize = 2 * nsmth[1] + 1;
    }
    let mut dalla = vec![vec![0.0; klast]; ncdp + isize as usize];

    fprintf!(fpq, "C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
    for i in 0..numstandard {
        fprintf!(fpq, ",%s", qform[i]);
    }
    for i in 0..iztuple {
        if ksize[i] > 0 {
            fprintf!(fpq, ",%s", formtv);
        }
    }
    for _i in 0..numdind {
        for k in 0..ktuple - 1 {
            if ksize[k + iztuple] > 0 {
                fprintf!(fpq, ",%s", formtv);
            }
        }
    }
    if ifixd != 2 {
        fprintf!(fpq, "\nC_SU_NDIMS,%s", ndims[0]);
        for _ in 1..numstandard {
            fprintf!(fpq, ",");
        }
        for i in 0..iztuple {
            if ksize[i] > 0 {
                fprintf!(fpq, ",");
            }
        }
        for i in 0..numdind as usize {
            for k in 0..ktuple - 1 {
                if ksize[k + iztuple] > 0 {
                    fprintf!(fpq, ",%.15g", dind[i]);
                }
            }
        }
    }
    fprintf!(fpq, "\nC_SU_NAMES\nC_SU_ID");
    for i in 0..numstandard {
        fprintf!(fpq, ",%s", qname[i]);
    }
    for i in 0..iztuple {
        if ksize[i] > 0 {
            fprintf!(fpq, ",%s", pname[i]);
        }
    }
    for _i in 0..numdind {
        for k in 0..ktuple - 1 {
            if ksize[k + iztuple] > 0 {
                fprintf!(fpq, ",%s", pname[k + iztuple]);
            }
        }
    }
    fprintf!(fpq, "\n");

    let mut ierr = 0;
    if !is3d {
        for j in 0..ncdp {
            for i in 0..klast {
                dalla[j][i] = recinfo[j].dlots[i];
            }
        }
        if navrg[0] > 0 {
            runsmo(
                &mut dalla,
                ncdp as i32 + 2 * navrg[0] + 1,
                klast,
                nback[0],
                navrg[0],
                &mut dtog,
                &mut ierr,
            );
            if ierr == 1 {
                err!("**** Error in averaging. Less than 5 Q-records in cin file and navrg>0");
            } else if ierr > 0 {
                err!("**** Error in navrg averaging for cin file. Parameter out-of-range.");
            }
        } else if navrg[0] == -1 {
            for i in 0..klast {
                dtog[i] = 0.0;
            }
            for j in 0..ncdp {
                for i in 0..klast {
                    dtog[i] += dalla[j][i];
                }
            }
            for j in 0..ncdp {
                for i in 0..klast {
                    dalla[j][i] = dtog[i] / ncdp as f64;
                }
            }
        }
        if nsmth[0] > 0 {
            runsmo(
                &mut dalla,
                ncdp as i32 + 2 * nsmth[0] + 1,
                klast,
                nback[0],
                -nsmth[0],
                &mut dtog,
                &mut ierr,
            );
            if ierr == 1 {
                err!("**** Error in smoothing. Less than 5 Q-records in cin file and nsmth>0");
            } else if ierr > 0 {
                err!("**** Error in nsmth smoothing for cin file. Parameter out-of-range.");
            }
        }
        for j in 0..ncdp {
            for i in 0..klast {
                if ksmth[i] > 0 {
                    recinfo[j].dlots[i] = dalla[j][i];
                }
            }
        }
    } else {
        if navrg[0] != 0 || nsmth[0] > 0 {
            for n in 0..lastac {
                for j in 0..lastab {
                    for i in 0..klast {
                        dalla[j][i] = recinfo[j + n * lastab].dlots[i];
                    }
                }
                if navrg[0] > 0 {
                    runsmo(
                        &mut dalla,
                        lastab as i32 + 2 * navrg[0] + 1,
                        klast,
                        nback[0],
                        navrg[0],
                        &mut dtog,
                        &mut ierr,
                    );
                    if ierr == 1 {
                        err!("**** Error in averaging. Less than 5 inline cells defined in 3D Grid and navrg>0");
                    } else if ierr > 0 {
                        err!("**** Error in navrg averaging for 3D Grid inlines. Parameter out-of-range.");
                    }
                } else if navrg[0] == -1 {
                    for i in 0..klast {
                        dtog[i] = 0.0;
                    }
                    for j in 0..lastab {
                        for i in 0..klast {
                            dtog[i] += dalla[j][i];
                        }
                    }
                    for j in 0..lastab {
                        for i in 0..klast {
                            dalla[j][i] = dtog[i] / lastab as f64;
                        }
                    }
                }
                if nsmth[0] > 0 {
                    runsmo(
                        &mut dalla,
                        lastab as i32 + 2 * nsmth[0] + 1,
                        klast,
                        nback[0],
                        -nsmth[0],
                        &mut dtog,
                        &mut ierr,
                    );
                    if ierr == 1 {
                        err!("**** Error in smoothing. Less than 5 inline cells defined in 3D Grid and nsmth>0.");
                    } else if ierr > 0 {
                        err!("**** Error in nsmth smooth for 3D Grid inlines. Parameter out-of-range.");
                    }
                }
                for j in 0..lastab {
                    for i in 0..klast {
                        if ksmth[i] > 0 {
                            recinfo[j + n * lastab].dlots[i] = dalla[j][i];
                        }
                    }
                }
            }
        }
        if navrg[1] != 0 || nsmth[1] > 0 {
            for n in 0..lastab {
                for j in 0..lastac {
                    for i in 0..klast {
                        dalla[j][i] = recinfo[n + j * lastab].dlots[i];
                    }
                }
                if navrg[1] > 0 {
                    runsmo(
                        &mut dalla,
                        lastac as i32 + 2 * navrg[1] + 1,
                        klast,
                        nback[1],
                        navrg[1],
                        &mut dtog,
                        &mut ierr,
                    );
                    if ierr == 1 {
                        err!("**** Error averaging. Less than 5 crossline cells defined in 3D Grid and navrg>0");
                    } else if ierr > 0 {
                        err!("**** Error in navrg average for 3D Grid crosslines. Parameter out-of-range.");
                    }
                }
                if nsmth[1] > 0 {
                    runsmo(
                        &mut dalla,
                        lastac as i32 + 2 * nsmth[1] + 1,
                        klast,
                        nback[1],
                        -nsmth[1],
                        &mut dtog,
                        &mut ierr,
                    );
                    if ierr == 1 {
                        err!("**** Error smoothing. Less than 5 crossline cells defined in 3D Grid and nsmth>0");
                    } else if ierr > 0 {
                        err!("**** Error in nsmth smooth for 3D Grid crosslines. Parameter out-of-range.");
                    }
                } else if navrg[1] == -1 {
                    for i in 0..klast {
                        dtog[i] = 0.0;
                    }
                    for j in 0..lastac {
                        for i in 0..klast {
                            dtog[i] += dalla[j][i];
                        }
                    }
                    for j in 0..lastac {
                        for i in 0..klast {
                            dalla[j][i] = dtog[i] / lastac as f64;
                        }
                    }
                }
                for j in 0..lastac {
                    for i in 0..klast {
                        if ksmth[i] > 0 {
                            recinfo[n + j * lastab].dlots[i] = dalla[j][i];
                        }
                    }
                }
            }
        }
    }
    if cdpt == 1 {
        for jcdp in 0..ncdp {
            let kinf = &mut recinfo[jcdp].kinf;
            let mut cdp = 0;
            gridiccdp90(&gvals, kinf[1], kinf[2], &mut cdp);
            kinf[0] = cdp;
        }
        // SAFETY: all records are initialized and have three kinf entries.
        // qsort only relocates whole QInfo values; the callback cannot unwind.
        unsafe {
            qsort(
                recinfo.as_mut_ptr().cast(),
                ncdp,
                std::mem::size_of::<QInfo>(),
                comp_sort,
            );
        }
        for jcdp in 0..ncdp {
            let kinf = &mut recinfo[jcdp].kinf;
            let mut cdp = 0;
            gridiccdp(&gvals, kinf[1], kinf[2], &mut cdp);
            kinf[0] = cdp;
        }
    }

    for jcdp in 0..ncdp {
        if is3d {
            let kinf = &recinfo[jcdp].kinf;
            let mut icdpt = 0;
            let (mut xg, mut yg, mut xw, mut yw) = (0.0, 0.0, 0.0, 0.0);
            gridiccdp90(&gvals, kinf[1], kinf[2], &mut icdpt);
            gridicgridxy(&gvals, kinf[1], kinf[2], &mut xg, &mut yg);
            gridicrawxy(&gvals, kinf[1], kinf[2], &mut xw, &mut yw);
            fprintf!(
                fpq,
                &formxylong,
                kinf[0],
                icdpt,
                kinf[1],
                kinf[2],
                xg,
                yg,
                xw,
                yw
            );
        } else {
            fprintf!(fpq, "Q,%d,1", recinfo[jcdp].kinf[0]);
        }
        for j in 0..klast {
            dswap[j] = recinfo[jcdp].dlots[j];
        }
        for j in 0..iztuple {
            if ksize[j] > 0 {
                fprintf!(fpq, &formtvlong, dswap[kindx[j] as usize]);
            }
        }
        for i in 0..numdind as usize {
            for k in (0..ktuple - 1).rev() {
                if ksize[k + iztuple] > 0 {
                    fprintf!(fpq, &formtvlong, dswap[kindx[k + iztuple] as usize + i]);
                }
            }
        }
        fprintf!(fpq, "\n");
    }
    exit(CWP_Exit());
}

// Running average, copied from sufarldcsv February 2024.
// The running total is stored behind the active window, then shifted back.
fn runav(
    dall: &mut [Vec<f64>],
    ncdp: usize,
    klast: usize,
    nmin: usize,
    mut nmax: i32,
    dtog: &mut [f64],
    ierr: &mut i32,
) {
    let mut iavr = 1;
    if nmax < 0 {
        nmax = -nmax;
        iavr = 0;
    }
    let nmax = nmax as usize;
    if ncdp < 1 || klast < 1 || nmax < nmin {
        *ierr = 2;
        return;
    }
    if 2 * (nmax + 3) > ncdp {
        *ierr = 1;
        return;
    }
    *ierr = 0;
    let mut n = 0i32;
    for k in 0..klast {
        dtog[k] = 0.0;
    }
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
}

// Running smoother (or averager). nmaxi < 0 runs the signed average
// followed by the ordinary average, with endpoint extrapolation each time.
fn runsmo(
    dalla: &mut [Vec<f64>],
    msize: i32,
    klast: usize,
    nback: i32,
    nmaxi: i32,
    dtog: &mut [f64],
    ierr: &mut i32,
) {
    let nmax = if nmaxi < 0 { -nmaxi } else { nmaxi };
    if msize < 1 || klast < 1 || nmax < 1 || nback < 0 {
        *ierr = 2;
        return;
    }
    if 2 * nmax + 6 > msize {
        *ierr = 1;
        return;
    }
    let lastcdp = (msize - 2 * nmax - 1) as usize;
    // Both extrapolation passes need real locations at first+nback and
    // last-nback. C can access invalid memory when this documented bound
    // is exceeded; report the caller's parameter-range error instead.
    if nback as usize >= lastcdp {
        *ierr = 2;
        return;
    }
    let mut dback = 0.0;
    if nback > 0 {
        dback = 1.0 / nback as f64;
    }
    *ierr = 0;
    let nmax = nmax as usize;
    let nback = nback as usize;
    for j in (0..lastcdp).rev() {
        for i in 0..klast {
            dalla[j + 1 + nmax][i] = dalla[j][i];
        }
    }
    for j in 0..nmax + 1 {
        for i in 0..klast {
            let dt = (dalla[nmax + 1][i] - dalla[nmax + 1 + nback][i]) * dback;
            dalla[j][i] = dalla[nmax + 1][i] + (nmax - j + 1) as f64 * dt;
        }
    }
    for j in lastcdp + nmax + 1..lastcdp + 2 * nmax + 1 {
        for i in 0..klast {
            let dt = (dalla[lastcdp + nmax][i] - dalla[lastcdp + nmax - nback][i]) * dback;
            dalla[j][i] = dalla[lastcdp + nmax][i] + (j - lastcdp - nmax) as f64 * dt;
        }
    }
    runav(dalla, lastcdp + 2 * nmax + 1, klast, 0, nmaxi, dtog, ierr);
    if *ierr > 0 {
        return;
    }
    for j in 0..lastcdp {
        for i in 0..klast {
            dalla[j][i] = dalla[j + nmax][i];
        }
    }
    if nmaxi > 0 {
        return;
    }

    for j in (0..lastcdp).rev() {
        for i in 0..klast {
            dalla[j + 1 + nmax][i] = dalla[j][i];
        }
    }
    for j in 0..nmax + 1 {
        for i in 0..klast {
            let dt = (dalla[nmax + 1][i] - dalla[nmax + 1 + nback][i]) * dback;
            dalla[j][i] = dalla[nmax + 1][i] + (nmax - j + 1) as f64 * dt;
        }
    }
    for j in lastcdp + nmax + 1..lastcdp + 2 * nmax + 1 {
        for i in 0..klast {
            let dt = (dalla[lastcdp + nmax][i] - dalla[lastcdp + nmax - nback][i]) * dback;
            dalla[j][i] = dalla[lastcdp + nmax][i] + (j - lastcdp - nmax) as f64 * dt;
        }
    }
    runav(dalla, lastcdp + 2 * nmax + 1, klast, 0, -nmaxi, dtog, ierr);
    if *ierr > 0 {
        return;
    }
    for j in 0..lastcdp {
        for i in 0..klast {
            dalla[j][i] = dalla[j + nmax][i];
        }
    }
}

// C compSort: sort by cdp (or by the temporary transposed cdp).
unsafe extern "C" fn comp_sort(q1: *const c_void, q2: *const c_void) -> i32 {
    // SAFETY: qsort supplies pointers to initialized QInfo elements.
    let p1 = unsafe { &*q1.cast::<QInfo>() };
    let p2 = unsafe { &*q2.cast::<QInfo>() };
    if p1.kinf[0] < p2.kinf[0] {
        return -1;
    }
    if p1.kinf[0] > p2.kinf[0] {
        return 1;
    }
    0
}
