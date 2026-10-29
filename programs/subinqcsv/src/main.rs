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

//! SUBINQCSV - Input Q-records or parameters, interpolate, output Q-records.
//! Port of su/main/headers/subinqcsv.c (Revision 1.02).
//! Credits: Created Nov 2021, Andre Latour.

use oxwave_bilinear_interp::*;
use oxwave_core::cio::fopen;
use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;
use oxwave_linear_interp::*;
use oxwave_q_file::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUBINQCSV - Input Q-records or parameters. interpolate, output Q-records.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  subinqcsv [parameters].   (No traces in or out).                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " qin=       Input file containing q-records. This parameter is optional,    ",
    "            but if you do not specify it, you must specify a set of         ",
    "            parameters which you can see lower down in this document.       ",
    "            Those parameters allow specification similar to sunmo.          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qout=      Output file for q-records. Must be specified.                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rfile=     If set, read a K-file containing 3D grid definition. Assumed 2D ",
    "            if no file is specified here. A 3D forces restrictions on the   ",
    "            input locations (they have to form aligned rectangles).         ",
    "            A 2D means that parameters for the igi component of 3Ds will be ",
    "            honored, but parameters for the igc component are not allowed.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " check=0   Do not print checking details. The options of this program       ",
    "           cannot repair a corrupt input grid. Grid checking is always      ",
    "           done on input grid definition and will error-halt if corrupt.    ",
    "      =1   print checking details.                                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " inloc=1   Location options. The option here determines which input names   ",
    "           contain the locations of all other input values. The option here ",
    "           must find the corresponding names within the non-tuple names in  ",
    "           the input q-file or in the parameter names you specify (later).  ",
    "      =1   cdp (default). Option 1 is required if no grid is input.         ",
    "      =2   cdpt (grid transposed cdp numbers).                              ",
    "      =3   igi,igc (grid index numbers).                                    ",
    "      =4   gx,gy (grid coordinates).  Note option 6.                        ",
    "      =5   sx,sy (world coordinates). Note option 7.                        ",
    "      =6   gx,gy (world coordinates).                                       ",
    "      =7   sx,sy (grid coordinates).                                        ",
    "           You may want to use iecho=1 or 2 initially with these options    ",
    "           since the output q-records will contain the cell centre values   ",
    "           corresponding to your choice here (giving you a chance to        ",
    "           adjust your input values before going through the error-check    ",
    "           that forces the input locations to be in aligned rectangles).    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keyloc=   Location key. Can only be specified if inloc=1 and 2D (no rfile).",
    "           Name must be in the input Q-file (or parameters listed below).   ",
    "           This allows linear interpolation based on values such as point   ",
    "           numbers, station numbers, reciever numbers, or shot numbers.     ",
    "           Parameters extrapi and igiout refer to these values.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " extrapi=0        do not extrapolate at ends in igi direction.              ",
    "                  (values beyond outer functions are held constant).        ",
    "               =1 extrapolate beyond outer functions at both ends           ",
    "               =2 extrapolate beyond outer functions only at lower end      ",
    "               =3 extrapolate beyond outer functions only at higher end     ",
    "               For 2D, this parameter controls extrapolation at line ends.  ",
    "                                                                            ",
    " extrapc=0        do not extrapolate at ends in igc direction.              ",
    "                  (values beyond outer functions are held constant).        ",
    "               =1 extrapolate beyond outer functions at both ends           ",
    "               =2 extrapolate beyond outer functions only at lower end      ",
    "               =3 extrapolate beyond outer functions only at higher end     ",
    "                                                                            ",
    " extrapt=0        do not extrapolate at ends in independent values.         ",
    "                  (values outside range of functions are held constant).    ",
    "               =1 extrapolate outside range of functions                    ",
    "               =2 extrapolate when less than range of functions             ",
    "               =3 extrapolate when greater than range of functions          ",
    "                                                                            ",
    " igiout=1   Output igi range (default is entire igi range of grid).         ",
    "            List with one, two, or three integers.                          ",
    "       =f     start at f, increment by 1, end at maximum igi of grid.       ",
    "       =f,l   start at f, increment by 1, end at l or maximum igi of grid.  ",
    "       =f,l,i start at f, increment by i, end at l or nearest igi that      ",
    "              is not greater than maximum igi of grid.                      ",
    "             For 2D, f and l default to minimum and maximum numbers within  ",
    "             the cdp= list or in the q-records input via qin=. For 2D, they ",
    "             are allowed to be set outside that range (you can extrapolate  ",
    "             to the true first and last cdp of line).                       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " igcout=1   Output igc range (default is entire igc range of grid).         ",
    "            Same options as the previous igiout parameter.                  ",
    "\t\t\t\t\t\t\t\t\t     ",
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
    " invp2=v    Square and inverse values before applying bilinear weighting.   ",
    "            (This is how velocities are treated in sunmo and sunmocsv).     ",
    "            The sequence is: linear interpolate to outind or outlist values,",
    "            square and inverse, perform weighted-sum based on location,     ",
    "            then inverse square-root that sum. The default is to do this if ",
    "            the value name starts with v. You can also specify any other    ",
    "            characters such as vq, tzm, qtal and so on.                     ",
    "       =all use square and inverse for all values.                          ",
    "      =none use square and inverse for no values.                           ",
    "     =tuple use square and inverse for all values in the tuples.            ",
    "  =nontuple use square and inverse for all values not in the tuples.        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " ename=     Names to eliminate. These names will be ignored in the input    ",
    "            q-records and therefore not be in the output q-records (note    ",
    "            that standard 2d and 3d names will still be output even if you  ",
    "            specify them here). A primary purpose of this option is to get  ",
    "            rid of names in input tuples in order to reduce the complexity  ",
    "            of output q-records (in case certain programs cannot handle it).",
    "            But you can also get rid of non-tuple names here.               ",
    "       ***  Note this option can be used to eliminate values even if you    ",
    "            choose to input values via the parameters below rather than     ",
    "            via a q-file in qin=. This reduces manual editing of input      ",
    "            values in certain situations.                                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " iecho=0    Full spatial output (do not just output at the input locations).",
    "      =1    Just output at the same locations as input.                     ",
    "            This still performs interpolation in the independent dimension  ",
    "            as specified by outind or outlist. (This option bypasses the    ",
    "            error check that input locations must form aligned rectangles). ",
    "              Parameters igiout= and igcout= cannot be specified.           ",
    "      =2    Just output at the same locations as input. This option         ",
    "            does NOT perform interpolation in the independent dimension,    ",
    "            BUT IT DOES output q-files with varying numbers of tuples.      ",
    "            For instance if your input parameters have varying amounts of   ",
    "            velocity picks, it outputs them as a varying tuples q-file.     ",
    "       ***  Note this option allows you to input values via parameters      ",
    "            below and then output a varying tuples q-file. A q-file is much ",
    "            easier to input into spreadsheet programs in order to perform   ",
    "            computations (such as Dix equation velocity transforms).        ",
    "              Parameters outind and outlist cannot be specified.            ",
    "              Parameters igiout= and igcout= cannot be specified.           ",
    "                                                                            ",
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
    "   The following parameters can only be specified if you                    ",
    "   do not specify an input q-file using parameter qin=.                     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cdp=             cdp is a likely parameter name to specify here but        ",
    "                  all seismic unix key names are allowed as parameter names.",
    "                  So, for instance, fldr= is also allowed to be specified.  ",
    "                  These parameters are all optional but at least one of     ",
    "                  them must be specified. These are the names which         ",
    "                  contain single values at each location (also known as     ",
    "            \t   the non-tuple names). Each name specified here must have  ",
    "            \t   the same amount of values (separated by commas).          ",
    "            \t   For example: cdp=2,7,23,44 means that you must also       ",
    "            \t   list 4 values if you specify tstat=.                      ",
    "              *** Note the names specified must include the names needed    ",
    "                  for your inloc= option (by default inloc= requires cdp).  ",
    " numa=        *** In addition to all seismic unix key names, the four names ",
    " numb=            to the left are ALSO allowed.                             ",
    " vuma=        *** See invp2= for why some names here start with v.          ",
    " vumb=                                                                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " tupa=            The parameters to the left are all optional. These are    ",
    " offs=            the names which have multiple values at each location     ",
    " tims=            (also known as the tuple names). Typically, tims= and     ",
    " tnmo=            vnmo= are the most likely names to specify. Each name     ",
    " dpth=      \t   specified here must be repeated same number of times as   ",
    " vels=            amount of values in non-tuple names above (except Note 1).",
    " vnmo=            (See invp2= for why some names here start with v).        ",
    " tupb=    Note 1: Usually, any of the previous parameter names must be      ",
    " vupa=            repeated for each value in the non-tuple lists. But the   ",
    "                  independent dimension name can also be specified just one ",
    "                  time if all the other names here have the same amount of  ",
    "                  values as that single set of independent dimension values.",
    "          Note 2: The independent dimension defaults to first name that is  ",
    "            \t   found in order here (NOT their order on the command line).",
    "\t\t\t\t\t\t\t\t\t     ",
    "\tExample using 3 cpds with velocity semblance scan picks:             ",
    "\t\t\t\t\t\t\t\t\t     ",
    "         cdp=7,27,57                                                        ",
    "         tims=500,780,1240,2100,4000                                        ",
    "         vels=1500,1930,2260,3100,4400                                      ",
    "         tims=400,680,2400,3800                                             ",
    "         vels=1300,1730,3700,4300                                           ",
    "         tims=400,680,880,1140,2400,3800                                    ",
    "         vels=1300,1500,1730,2300,3700,4300                                 ",
    "\t\t\t\t\t\t\t\t\t     ",
    "       And, if you have a floating datum static at these 3 cdps, then also: ",
    "         tstat=5.30,17.44,12.63                                             ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  Notes:                                                       \t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   All input values are stored in double-precision (8 byte float).          ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   This program does not care if the values are seconds or milliseconds     ",
    "   or feet or metres. To this program they are all just numbers. However,   ",
    "   default of parameter invp2= treats names starting with v differently     ",
    "   for weighting.                                                           ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " For 3D, user needs to input locations which form aligned rectangles.       ",
    " That is, how-ever-many grid inlines the user chooses to put locations on,  ",
    " there must always be the same number of locations on each inline and those ",
    " functions must be located at the same grid crosslines. For instance, if    ",
    " the user inputs locations for inline 7 at crosslines 15,25,40 then the     ",
    " the user must alos input locations at crosslines 15,25,40 for any other    ",
    " inlines that the user wants to supply locations for. (If user is lucky     ",
    " enough that the grid has 100 inlines, then the input locations could be at ",
    " CDPs 715,725,740 and 1115,1125,1140 and 2015,2025,2040 and 2515,2525,2540. ",
    " Note that the inloc= option allows different values to be used to specify  ",
    " locations but they are all translated to inline and crossline numbers      ",
    " of the cell centre using the input 3D grid definition.                     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Bilinear interpolation is done if the output location is surrounded by 4   ",
    " input locations. If the output location is not surrounded by 4 input       ",
    " locations, the result depends on the extrapi and extrapc options.          ",
    " The result can be any combination of linear interpolation, linear          ",
    " extrapolation, and constant extrapolation. If input locations only exist   ",
    " along 1 inline or 1 crossline the result is linear interpolation in that   ",
    " direction (and linear or constant extrapolation at locations beyond ends). ",
    "\t\t\t\t\t\t\t\t\t     ",
];

const ZNAME: &[&str] = &[
    "cdp", "cdpt", "fldr", "grnors", "grnofr", "grnlof", "gaps", "igi", "igc", "gx", "gy", "sx",
    "sy", "numa", "numb", "vuma", "vumb", "tracl", "tracr", "tracf", "ep", "trid", "nvs", "nhs",
    "duse", "offset", "gelev", "selev", "sdepth", "gdel", "sdel", "swdep", "gwdep", "scalel",
    "scalco", "counit", "wevel", "swevel", "sut", "gut", "sstat", "gstat", "tstat", "laga", "lagb",
    "delrt", "muts", "mute", "ns", "dt", "gain", "corr", "sfs", "sfe", "slen", "styp", "stas",
    "stae", "tatyp", "afilf", "afils", "nofilf", "nofils", "lcf", "hcf", "lcs", "hcs", "year",
    "day", "hour", "minute", "sec", "timbas", "trwf", "otrav", "d1", "f1", "d2", "f2", "ungpow",
    "unscale", "ntr", "mark", "shortpad", "tupa", "offs", "tims", "tnmo", "dpth", "vels", "vnmo",
    "tupb", "vupa",
];

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
    let mut mgiextr = 0;
    getparint("extrapi", &mut mgiextr);
    if mgiextr < 0 || mgiextr > 3 {
        err!("error: extrapi= option not in range ");
    }
    let mut mgcextr = 0;
    getparint("extrapc", &mut mgcextr);
    if !is3d && mgcextr != 0 {
        err!("**** Error: cannot specify mgcextr= with no input 3d grid.");
    }
    if mgcextr < 0 || mgcextr > 3 {
        err!("error: extrapc= option not in range ");
    }
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
    let has_qin = getparstring("qin", &mut pname_file) != 0;
    let mut qname_file = String::new();
    let has_qout = getparstring("qout", &mut qname_file) != 0;
    if has_qin && has_qout && pname_file == qname_file {
        err!("**** Error: qout= output Q-file name must be different than qin= input Q-file name.");
    }
    let fpq = fopen(&qname_file, "w").unwrap_or_else(|| {
        err!("qfile error: output Q-file did not open correctly.");
    });
    let mut iecho = 0;
    getparint("iecho", &mut iecho);
    if iecho < 0 || iecho > 2 {
        err!("error: iecho= option not in range ");
    }
    if iecho > 1 && (countparval("outind") > 0 || countparval("outlist") > 0) {
        err!("error: iecho=2 cannot be specified with outind= or outlist=");
    }
    if iecho > 0 && (countparval("igiout") > 0 || countparval("igcout") > 0) {
        err!("error: iecho>0 cannot be specified with igiout= or igcout=");
    }
    if countparval("igiout") > 3 {
        err!("error: igiout= has too long a list.");
    }
    if countparval("igcout") > 3 {
        err!("error: igcout= has too long a list.");
    }
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
        let numdind = (1.00001 + (outind[1] - outind[0]) / outind[2]) as i32;
        for i in 0..numdind {
            dind.push(outind[0] + i as f64 * outind[2]);
        }
    } else if countparval("outlist") > 0 {
        getpardouble("outlist", &mut dind);
    }
    let numdind = dind.len() as i32;
    let mut nvp2 = "v".to_owned();
    getparstring("invp2", &mut nvp2);
    let invp2 = match nvp2.as_str() {
        "none" => 0,
        "all" => 2,
        "tuple" => 3,
        "nontuple" => 4,
        _ => 1,
    };
    let mut formxy = "%.2f".to_owned();
    if getparstring("formxy", &mut formxy) != 0 && !is3d {
        err!("**** Error: cannot specify formxy= with no input 3d grid.");
    }
    let formxylong = format!("Q,%d,%d,%d,%d,{0},{0},{0},{0}", formxy);
    let mut qname = vec![keyloc.clone(), "cdpt".to_owned()];
    let mut qform = vec!["%.20g".to_owned(); 2];
    if is3d {
        for s in ["igi", "igc", "gx", "gy", "sx", "sy"] {
            qname.push(s.to_owned());
        }
        for _ in 0..2 {
            qform.push("%.20g".to_owned());
        }
        for _ in 0..4 {
            qform.push(formxy.clone());
        }
    }
    let numstandard = qname.len();
    let mut formtv = "%.2f".to_owned();
    getparstring("formtv", &mut formtv);
    let formtvlong = format!(",{}", formtv);
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
        match errwarn {
            1 => err!("gridset error: grid_wb cell width must be positive."),
            2 => err!("gridset error: grid_wc cell width must be positive."),
            3 => err!("gridset error: corner B is within grid_wb cell width of corner A."),
            n if n > 0 => err!("gridset error: returned with some unrecognized error code."),
            -1 => warn!("gridset warning: corner C is near A and is reset to A."),
            _ => {}
        }
        gridcheck(&gvals, icheck, &mut errwarn);
        if errwarn > 0 {
            err!("gridcheck error: returned with some unrecognized error code.");
        }
    }

    // Read scalar names first, followed by tuple names in their C order.
    let mut ename = Vec::<String>::new();
    let mut q = if !has_qin {
        let mut pname = Vec::new();
        let mut iztuple = 0;
        for (j, name) in ZNAME.iter().enumerate() {
            if countparname(name) > 0 {
                pname.push((*name).to_owned());
                if j < 84 {
                    iztuple += 1;
                }
            }
        }
        getviacommand(&mut pname, iztuple, numdind).unwrap_or_else(|n| command_error(n))
    } else {
        let fpp = fopen(&pname_file, "r").unwrap_or_else(|| {
            err!("error: input Q-file did not open correctly.");
        });
        if countparval("ename") > 0 {
            getparstringarray("ename", &mut ename);
        }
        let mut exclude = ename.clone();
        if is3d {
            if inloc != 1 {
                exclude.push("cdp".to_owned());
            }
            if inloc != 2 {
                exclude.push("cdpt".to_owned());
            }
            if inloc != 3 {
                exclude.extend(["igi".to_owned(), "igc".to_owned()]);
            }
            if inloc != 4 && inloc != 6 {
                exclude.extend(["gx".to_owned(), "gy".to_owned()]);
            }
            if inloc != 5 && inloc != 7 {
                exclude.extend(["sx".to_owned(), "sy".to_owned()]);
            }
        }
        getviaqfile(&fpp, &exclude, -(exclude.len() as i32), numdind)
            .unwrap_or_else(|n| qfile_error(n))
    };
    let ifixd = q.ifixd;
    let iz = q.iztuple as usize;
    let kt = q.ktuple as usize;
    if iecho == 2 {
        if ifixd == 1 {
            err!("error: for iecho=2 input must have varying number of tuples.");
        }
    } else if numdind == 0 && ifixd != 2 {
        err!("error: input has tuples so outind= or outlist= must be specified.");
    } else if numdind > 0 && ifixd == 2 {
        err!("error: input does not have tuples so outind= and outlist= cannot be specified.");
    }
    if ifixd == 0 && iecho < 2 {
        // The C reader allocates dimension-name slots even without NDIMS.
        if q.ndims.is_empty() {
            q.ndims.push(String::new());
        }
        q.ndims[0] = q.names[iz].clone();
        q.names.remove(iz);
    }
    // getnpar's empty numeric parameter has a count but writes no values.
    // Keep the allocated independent-dimension slot for the one-point case;
    // its value is immaterial to linterpd with a single input point.
    if ifixd == 1 && q.records.first().is_some_and(|r| r.nto == 1) && q.pindepa.is_empty() {
        q.pindepa.push(0.0);
    }
    if ifixd == 0 {
        for (j, rec) in q.records.iter().enumerate() {
            let pindepa = &rec.dlots[iz + (kt - 1) * rec.nto as usize..];
            for i in 1..rec.nto as usize {
                if pindepa[i - 1] >= pindepa[i] {
                    err!(
                        "error: independent dimension values not increasing (input location = %d)",
                        j + 1
                    );
                }
            }
        }
    } else if ifixd == 1 {
        for i in 1..q.records[0].nto as usize {
            if q.pindepa[i - 1] >= q.pindepa[i] {
                err!("error: independent dimension values not increasing.");
            }
        }
    }
    set_locations(&mut q, is3d, inloc, &keyloc, &gvals);
    let (mut minigi, mut maxigi, minigc, maxigc) = if is3d {
        (1, (0.1 + gvals[12]) as i32, 1, (0.1 + gvals[13]) as i32)
    } else {
        (2147483645, -2147483645, 1, 1)
    };
    if !is3d {
        for rec in &q.records {
            if rec.kinf[0] < minigi {
                minigi = rec.kinf[0];
            }
            if rec.kinf[0] > maxigi {
                maxigi = rec.kinf[0];
            }
        }
    }
    let mut igiout = [minigi, maxigi, 1];
    if countparval("igiout") > 0 {
        getparint("igiout", &mut igiout[..]);
    }
    if igiout[2] < 1 {
        err!("error: igiout= increment is less than 1");
    }
    if is3d {
        if igiout[0] < minigi {
            err!(
                "error: igiout= first is less than minimum grid igi index %d",
                minigi
            );
        }
        if igiout[1] > maxigi {
            igiout[1] = maxigi;
        }
    }
    if igiout[0] > igiout[1] {
        err!("error: igiout= first is greater than last. ");
    }
    let mut igcout = [minigc, maxigc, 1];
    if countparval("igcout") > 0 {
        getparint("igcout", &mut igcout[..]);
    }
    if igcout[2] < 1 {
        err!("error: igcout= increment is less than 1");
    }
    if is3d {
        if igcout[0] < minigc {
            err!(
                "error: igcout= first is less than minimum grid igc index %d",
                minigc
            );
        }
        if igcout[1] > maxigc {
            igcout[1] = maxigc;
        }
    }
    if igcout[0] > igcout[1] {
        err!("error: igcout= first is greater than last. ");
    }

    if iecho == 2 {
        if countparval("ename") > 0 {
            getparstringarray("ename", &mut ename);
        }
        checkpars();
        // Scalar input has an unused tuple-name slot in the C allocation.
        q.names.resize(iz + kt, "null".to_owned());
        let mut ksize = vec![1; iz + kt];
        for j in 0..iz + kt {
            if ename.contains(&q.names[j]) || qname.contains(&q.names[j]) {
                ksize[j] = -1;
            }
        }
        fprintf!(fpq, "C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
        for s in &qform {
            fprintf!(fpq, ",%s", s);
        }
        for j in 0..iz {
            if ksize[j] > 0 {
                fprintf!(fpq, ",%s", formtv);
            }
        }
        if ifixd != 2 {
            for _ in 0..2 {
                for k in 0..kt {
                    if ksize[k + iz] > 0 {
                        fprintf!(fpq, ",%s", formtv);
                    }
                }
            }
        }
        fprintf!(fpq, "\nC_SU_NAMES\nC_SU_ID");
        for s in &qname {
            fprintf!(fpq, ",%s", s);
        }
        for j in 0..iz {
            if ksize[j] > 0 {
                fprintf!(fpq, ",%s", q.names[j]);
            }
        }
        if ifixd != 2 {
            for _ in 0..2 {
                for k in 0..kt {
                    if ksize[k + iz] > 0 {
                        fprintf!(fpq, ",%s", q.names[k + iz]);
                    }
                }
            }
        }
        fprintf!(fpq, "\n");
        for rec in &q.records {
            write_location(&fpq, is3d, &gvals, &formxylong, &rec.kinf);
            for j in 0..iz {
                if ksize[j] > 0 {
                    fprintf!(fpq, &formtvlong, rec.dlots[j]);
                }
            }
            for i in 0..rec.nto as usize {
                for k in (0..kt).rev() {
                    if ksize[k + iz] > 0 {
                        fprintf!(fpq, &formtvlong, rec.dlots[iz + k * rec.nto as usize + i]);
                    }
                }
            }
            fprintf!(fpq, "\n");
        }
        exit(CWP_Exit());
    }
    checkpars();
    let (mut mgi, mut mgc) = (Vec::new(), Vec::new());
    if iecho == 0 {
        q.records
            .sort_by(|a, b| a.kinf[2].cmp(&b.kinf[2]).then(a.kinf[1].cmp(&b.kinf[1])));
        for j in 1..q.records.len() {
            if q.records[j - 1].kinf[1] == q.records[j].kinf[1]
                && q.records[j - 1].kinf[2] == q.records[j].kinf[2]
            {
                err!(
                    "error: Two sets of values input for cdp,igi,igc = %d %d %d",
                    q.records[j].kinf[0],
                    q.records[j].kinf[1],
                    q.records[j].kinf[2]
                );
            }
        }
        (mgi, mgc) = qsplit(&q.records).unwrap_or_else(|n| {
            err!(
                "qsplit error: locations do not form aligned rectangles (detected at cdp= %d)",
                q.records[n as usize - 1].kinf[0]
            );
        });
    }
    let mut kindx = vec![0; iz + kt - 1];
    let mut ksize = vec![0; iz + kt - 1];
    qelementout(iz as i32, kt as i32, numdind, &mut kindx, &mut ksize);
    let klast = (kindx[iz + kt - 2] + ksize[iz + kt - 2]) as usize;
    let mut dswap = vec![0.0; klast];
    for rec in &mut q.records {
        for j in 0..iz {
            dswap[kindx[j] as usize] = rec.dlots[j];
        }
        let pindepa = if ifixd == 0 {
            &rec.dlots[iz + (kt - 1) * rec.nto as usize..]
        } else {
            &q.pindepa[..]
        };
        for i in 0..numdind as usize {
            for k in 0..kt - 1 {
                let start = iz + k * rec.nto as usize;
                linterpd(
                    dind[i],
                    &pindepa[..rec.nto as usize],
                    &rec.dlots[start..start + rec.nto as usize],
                    mgtextr,
                    &mut dswap[kindx[k + iz] as usize + i],
                );
            }
        }
        rec.dlots[..klast].copy_from_slice(&dswap);
    }
    let mut knvp2 = vec![1; iz + kt - 1];
    // Preserve C's (invp2==1 && strncmp(...)) == 0 expression.
    for k in 0..iz {
        if invp2 != 0
            && (invp2 == 2 || invp2 == 4 || !(invp2 == 1 && !q.names[k].starts_with(&nvp2)))
        {
            knvp2[k] = -2;
        }
    }
    for k in 0..kt - 1 {
        let i = kt - 2 - k;
        if invp2 != 0
            && (invp2 == 2 || invp2 == 3 || !(invp2 == 1 && !q.names[i + iz].starts_with(&nvp2)))
        {
            knvp2[k + iz] = -2;
        }
    }
    let mut kflag = false;
    for j in 0..iz + kt - 1 {
        if qname.contains(&q.names[j]) {
            ksize[j] = -ksize[j];
        } else if knvp2[j] == -2 {
            kflag = true;
        }
    }
    if kflag {
        for rec in &mut q.records {
            for k in 0..iz + kt - 1 {
                if knvp2[k] == -2 {
                    for i in kindx[k]..kindx[k] + ksize[k] {
                        let x = &mut rec.dlots[i as usize];
                        if *x > 1.0e-98 {
                            *x = 1.0 / (*x * *x);
                        } else if *x < -1.0e-98 {
                            *x = -1.0 / (*x * *x);
                        }
                    }
                }
            }
        }
    }
    fprintf!(fpq, "C_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
    for s in &qform {
        fprintf!(fpq, ",%s", s);
    }
    for j in 0..iz {
        if ksize[j] > 0 {
            fprintf!(fpq, ",%s", formtv);
        }
    }
    for _ in 0..numdind {
        for k in 0..kt - 1 {
            if ksize[k + iz] > 0 {
                fprintf!(fpq, ",%s", formtv);
            }
        }
    }
    if ifixd != 2 {
        fprintf!(fpq, "\nC_SU_NDIMS,%s", q.ndims[0]);
        for _ in 1..numstandard {
            fprintf!(fpq, ",");
        }
        for j in 0..iz {
            if ksize[j] > 0 {
                fprintf!(fpq, ",");
            }
        }
        for i in 0..numdind as usize {
            for k in 0..kt - 1 {
                if ksize[k + iz] > 0 {
                    fprintf!(fpq, ",%.15g", dind[i]);
                }
            }
        }
    }
    fprintf!(fpq, "\nC_SU_NAMES\nC_SU_ID");
    for s in &qname {
        fprintf!(fpq, ",%s", s);
    }
    for j in 0..iz {
        if ksize[j] > 0 {
            fprintf!(fpq, ",%s", q.names[j]);
        }
    }
    for _ in 0..numdind {
        for k in 0..kt - 1 {
            if ksize[k + iz] > 0 {
                fprintf!(fpq, ",%s", q.names[k + iz]);
            }
        }
    }
    fprintf!(fpq, "\n");
    if iecho == 1 {
        for rec in &q.records {
            write_location(&fpq, is3d, &gvals, &formxylong, &rec.kinf);
            dswap.copy_from_slice(&rec.dlots[..klast]);
            if kflag {
                inverse(&mut dswap, &kindx, &ksize, &knvp2);
            }
            write_values(&fpq, &formtvlong, &dswap, iz, kt, numdind, &kindx, &ksize);
        }
        exit(CWP_Exit());
    }
    let mgi_tot = mgi.len() as i32;
    let mgc_tot = mgc.len() as i32;
    let mgi_totdeg = if mgi_tot == 1 || mgc_tot == 1 {
        0
    } else {
        mgi_tot
    };
    let (mut ndxi, mut ndxc, mut mdxi, mut mdxc) = (0, 0, 0, 0);
    let (mut wi, mut wc) = (0.0, 0.0);
    let mut kigc = igcout[0];
    while kigc <= igcout[1] {
        let mut kigi = igiout[0];
        while kigi <= igiout[1] {
            let mut icdp = kigi;
            if is3d {
                gridiccdp(&gvals, kigi, kigc, &mut icdp);
            }
            write_location(&fpq, is3d, &gvals, &formxylong, &[icdp, kigi, kigc]);
            if q.records.len() > 1 {
                let (mut mgix, mut mgcx) = (0, 0);
                binterpfind(
                    kigi, &mgi, mgiextr, kigc, &mgc, mgcextr, &mut mgix, &mut mgcx, &mut wi,
                    &mut wc,
                );
                ndxi = mgix + mgi_tot * (mgcx - 1);
                ndxc = ndxi + mgi_totdeg;
                mdxi = ndxi - 1;
                mdxc = ndxc - 1;
            }
            binterpapply(
                &q.records[mdxi as usize].dlots,
                &q.records[ndxi as usize].dlots,
                mgi_tot,
                wi,
                &q.records[mdxc as usize].dlots,
                &q.records[ndxc as usize].dlots,
                mgc_tot,
                wc,
                &mut dswap,
            );
            if kflag {
                inverse(&mut dswap, &kindx, &ksize, &knvp2);
            }
            write_values(&fpq, &formtvlong, &dswap, iz, kt, numdind, &kindx, &ksize);
            kigi = kigi.wrapping_add(igiout[2]);
        }
        kigc = kigc.wrapping_add(igcout[2]);
    }
    exit(CWP_Exit());
}

fn write_location(fp: &CFile, is3d: bool, gvals: &[f64], form: &str, kinf: &[i32]) {
    if is3d {
        let mut icdpt = 0;
        let (mut xg, mut yg, mut xw, mut yw) = (0.0, 0.0, 0.0, 0.0);
        gridiccdp90(gvals, kinf[1], kinf[2], &mut icdpt);
        gridicgridxy(gvals, kinf[1], kinf[2], &mut xg, &mut yg);
        gridicrawxy(gvals, kinf[1], kinf[2], &mut xw, &mut yw);
        fprintf!(fp, form, kinf[0], icdpt, kinf[1], kinf[2], xg, yg, xw, yw);
    } else {
        fprintf!(fp, "Q,%d,1", kinf[0]);
    }
}
fn inverse(dswap: &mut [f64], kindx: &[i32], ksize: &[i32], knvp2: &[i32]) {
    for k in 0..kindx.len() {
        if knvp2[k] == -2 {
            for i in kindx[k]..kindx[k] + ksize[k] {
                let x = &mut dswap[i as usize];
                if *x > 1.0e-98 {
                    *x = 1.0 / x.sqrt();
                } else if *x < -1.0e-98 {
                    *x = -1.0 / (0.0 - *x).sqrt();
                }
            }
        }
    }
}
fn write_values(
    fp: &CFile,
    form: &str,
    vals: &[f64],
    iz: usize,
    kt: usize,
    numdind: i32,
    kindx: &[i32],
    ksize: &[i32],
) {
    for j in 0..iz {
        if ksize[j] > 0 {
            fprintf!(fp, form, vals[kindx[j] as usize]);
        }
    }
    for i in 0..numdind as usize {
        for k in (0..kt - 1).rev() {
            if ksize[k + iz] > 0 {
                fprintf!(fp, form, vals[kindx[k + iz] as usize + i]);
            }
        }
    }
    fprintf!(fp, "\n");
}
fn set_locations(q: &mut QFile, is3d: bool, inloc: i32, keyloc: &str, gvals: &[f64]) {
    let (name1, name2) = match inloc {
        1 => (keyloc, ""),
        2 => ("cdpt", ""),
        3 => ("igi", "igc"),
        4 | 6 => ("gx", "gy"),
        _ => ("sx", "sy"),
    };
    let find = |name: &str| -> usize {
        for i in 0..q.iztuple as usize {
            if q.names[i] == name {
                return i;
            }
        }
        if inloc == 1 {
            err!(
                "error: for inloc=1, input must have your keyloc=%s (among non-tuple names).",
                keyloc
            );
        }
        match (inloc, name) {
            (2, "cdpt") => {
                err!("error: for inloc=2, input must have cdpt (among non-tuple names).")
            }
            (3, "igi") => err!("error: for inloc=3, input must have igi (among non-tuple names)."),
            (3, "igc") => err!("error: for inloc=3, input must have igc (among non-tuple names)."),
            (4, "gx") => err!("error: for inloc=4, input must have gx (among non-tuple names)."),
            (4, "gy") => err!("error: for inloc=4, input must have gy (among non-tuple names)."),
            (5, "sx") => err!("error: for inloc=5, input must have sx (among non-tuple names)."),
            (5, "sy") => err!("error: for inloc=5, input must have sy (among non-tuple names)."),
            (6, "gx") => err!("error: for inloc=6, input must have gx (among non-tuple names)."),
            (6, "gy") => err!("error: for inloc=6, input must have gy (among non-tuple names)."),
            (7, "sx") => err!("error: for inloc=7, input must have sx (among non-tuple names)."),
            (7, "sy") => err!("error: for inloc=7, input must have sy (among non-tuple names)."),
            _ => unreachable!(),
        }
    };
    let jnloc1 = find(name1);
    let jnloc2 = if name2.is_empty() { 0 } else { find(name2) };
    for rec in &mut q.records {
        let (mut icdp, mut igi, mut igc) = (0, 0, 0);
        if inloc == 1 {
            icdp = rec.dlots[jnloc1].round_ties_even() as i64 as i32;
            if is3d {
                gridcdpic(gvals, icdp, &mut igi, &mut igc);
                if igi < -2147483644 {
                    err!("error: for inloc=1, input cdp %d is not in grid", icdp);
                }
            } else {
                igi = icdp;
                igc = 1;
            }
        } else if inloc == 2 {
            let k = rec.dlots[jnloc1].round_ties_even() as i64 as i32;
            gridcdpic90(gvals, k, &mut igi, &mut igc);
            if igi < -2147483644 {
                err!(
                    "error: for inloc=2, input 90 degree cdpt %d is not in grid",
                    k
                );
            }
            gridiccdp(gvals, igi, igc, &mut icdp);
        } else if inloc == 3 {
            igi = rec.dlots[jnloc1].round_ties_even() as i64 as i32;
            igc = rec.dlots[jnloc2].round_ties_even() as i64 as i32;
            gridiccdp(gvals, igi, igc, &mut icdp);
            if icdp < -2147483644 {
                err!(
                    "error: for inloc=3 input igi,igc = %d,%d are not in grid",
                    igi,
                    igc
                );
            }
        } else {
            let x = rec.dlots[jnloc1];
            let y = rec.dlots[jnloc2];
            let (mut xw, mut yw) = (x, y);
            if inloc == 4 || inloc == 7 {
                gridgridxyrawxy(gvals, x, y, &mut xw, &mut yw);
            }
            gridrawxycdpic(gvals, xw, yw, &mut icdp, &mut igi, &mut igc);
            if icdp < -2147483644 {
                match inloc {
                    4 => err!(
                        "error: for inloc=4, input gx,gy = %g,%g are not in grid",
                        x,
                        y
                    ),
                    5 => err!(
                        "error: for inloc=5, input sx,sy = %g,%g are not in grid",
                        x,
                        y
                    ),
                    6 => err!(
                        "error: for inloc=6, input gx,gy = %g,%g are not in grid",
                        x,
                        y
                    ),
                    _ => err!(
                        "error: for inloc=7, input sx,sy = %g,%g are not in grid",
                        x,
                        y
                    ),
                }
            }
        }
        rec.kinf = vec![icdp, igi, igc];
    }
}

fn command_error(n: i32) -> QFile {
    match n {
        1 => err!("getviacommand error: no non-tuple name passed in."),
        2 => err!("getviacommand error: non-tuple names have different amounts of values."),
        3 => err!("getviacommand error: independent dimension parameter is empty."),
        4 => err!("getviacommand error: an independent dimension parameter is empty."),
        5 => err!("getviacommand error: members of tuple have different amounts at same location."),
        n if n > 0 => err!(
            "getviacommand error: returned unrecognized error number = %d",
            n
        ),
        _ => unreachable!(),
    }
}

fn qfile_error(n: i32) -> QFile {
    match n {
        1 => err!("getqinfo error: extra C_SU_NAMES record in q-file"),
        2 => err!("getqinfo error: extra C_SU_NDIMS record in q-file"),
        3 => err!("getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),
        11 => {
            err!("readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES")
        }
        12 => {
            err!("readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.")
        }
        22 => err!("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),
        23 => err!("getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
        24 => err!("getviaqfile error: C_SU_NDIMS blank where valid number expected"),
        25 => err!("getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
        26 => err!("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
        27 => err!("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),
        n if n > 100 => err!(
            "getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",
            n - 99
        ),
        n if n > 0 => err!("getviaqfile error: unrecognized error code %d", n),
        _ => unreachable!(),
    }
}
