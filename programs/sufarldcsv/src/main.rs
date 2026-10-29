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

//! SUFARLDCSV - Floating Approximate Robust Local Datum From Q-files.
//! Port of su/main/statics/sufarldcsv.c (Revision 1.01).
//! Credits: Created Oct 2023: Andre Latour (Andre Roger Latour).

use oxwave_core::prelude::*;
use oxwave_grid_geometry::*;
use oxwave_grid_parameters::*;
use oxwave_q_file::{getviaqfile, QInfo};

mod local;
use local::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUFARLDCSV - Floating Approximate Robust Local Datum From Q-files.         ",
    "\t\t\t\t\t\t\t\t\t     ",
    "              The Approximate Robust Local (ARL) algorithm is intended      ",
    "              to produce a reasonable smooth floating datum for the         ",
    "              typical spatial distributions of shots and receivers that     ",
    "              exist in many 2D and 3D surveys.                              ",
    "\t\t\t\t\t\t\t\t\t     ",
    "              This documentation uses the familiar words shot, receiver,    ",
    "              cdp, coordinates, and statics. And this program defaults to   ",
    "              Seismic Unix key names for various parameters. But nothing    ",
    "              in this program behaves differently if other words or names   ",
    "              are used. To this program all input values are just numbers.  ",
    "              For example, the same values are output whether the input     ",
    "              values are named statics, or elevations.                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    "     Caution: When using the output of this program you should use          ",
    "              cdp or igi,igc in the SUSHIFTCSV program to find the          ",
    "              cdp of the traces. Using nearest COORDINATES to find the      ",
    "              cdp does not produce cdp-consistent results because cdps      ",
    "              only contain one-half of their coordinate boundaries.         ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sufarldcsv [parameters]. (No traces in or out).                           ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " sin=sstat.csv  Input csv Q-file (usually containing Source information).   ",
    "    =none       No sin file. Value 0 is used for all sin statics.           ",
    "\t\t\t\t\t\t\t\t\t     ",
    " sloc=sx,sy     List of sin file location Names.                            ",
    "                For 2D surveys, 1 or 2 names are allowed.                   ",
    "                For 3D surveys, there must be 2 names and the first name    ",
    "                must have the x value and the second name must have the y.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " sstat=sstat    This is the default source static name.                     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " smult=1.0      Sin file Static Multiplier. Statics in Q-files should be    ",
    "                milliseconds, if not, you can scale them before use here.   ",
    "          Note: Also, to reverse static sign, specify negative here.        ",
    "                                                                            ",
    " rin=rstat.csv  Input csv Q-file (usually containing Receiver information). ",
    "    =none       No rin file. Value 0 is used for all receiver statics.      ",
    "          Note: Both sin and rin can be none since it might be usefull to   ",
    "                get a dout file formatted with 2D values or 3D Grid values. ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rloc=gx,gy     List of rin file location Names.                            ",
    "                For 2D surveys, 1 or 2 names are allowed.                   ",
    "                For 3D surveys, there must be 2 names and the first name    ",
    "                must have the x value and the second name must have the y.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rstat=gstat    This is the default rin file static name.                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rmult=1.0      Rin File Static Multiplier. Statics in Q-files should be    ",
    "                milliseconds, if not, you can scale them before use here.   ",
    "          Note: Also, to reverse static sign, specify negative here.        ",
    "                                                                            ",
    "  rfile=        If specified, read the 3D Grid defined in the K-record.     ",
    "           ***  MUST be specified for 3D and must NOT be specified for 2D.  ",
    "           ***  3D Grid must have at least 5 cells in the corresponding     ",
    "                direction if navrg>0,0 or nsmth>0,0                         ",
    "                See subincsv for 3D Grid documentation.                     ",
    "                If not a 3D survey, specify the next parameter (cin file).  ",
    "                                                                            ",
    " cin=none       Input csv Q-file containing CDP Information.                ",
    "    =none       Default. No cin file. An rfile with a 3D Grid must be input.",
    "           ***  MUST be specified for 2D and must NOT be specified for 3D.  ",
    "           ***  At least 5 Q-records must exist if nsmth>0 or navrg>0.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cloc=gx,gy     List of cin file Location Names. Only 1 or 2 names allowed. ",
    "           ***  Can only be specified for 2D.                               ",
    "          Note: The default is receiver coordinate key names because the    ",
    "                cin file is often created from receiver locations.          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cseq=cdp       This is the cin Q-file record sequencing name.              ",
    "           ***  Can only be specified for 2D.                               ",
    "     =asis      Means use cin Q-records in the order they are input.        ",
    "          Note: The input cin Q-record are sorted by these values and that  ",
    "                sequence is used by navrg, nsmth, and nback parameters.     ",
    "                The averaging and smoothing controlled by those parameters  ",
    "                just uses these sequential numbers. That is, those options  ",
    "                do not care if cseq actually contains a cdp number, and     ",
    "                they do not care whether gaps exist in input cseq numbers.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cstat=tstat    This is the default output static name. Typically it is the ",
    "                sum of the (averaged, smoothed) output shot static and the  ",
    "                (averaged, smoothed) output receiver static divided by 2,   ",
    "                but if either sin or rin files are not input then it is     ",
    "                just a copy of the other (averaged, smoothed) static.       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " avrad=1000.0   Radius to Average Shot statics and Average Receiver statics.",
    "                Statics are averaged separately for shots and receivers.    ",
    "          =-1.  Do not average any statics. Just use nearest shot static    ",
    "                and nearest receiver static. You must also set nedge=0.     ",
    "          Note: Specify this parameter in units of sloc, rloc, and cloc.    ",
    "                By default those are coordinate names so this parameter     ",
    "                is a true distance. But for 2D surveys sloc, rloc, and cloc ",
    "                can be single values such as point numbers in which case    ",
    "                this value is the difference in point numbers.              ",
    "\t\t\t\t\t\t\t\t\t     ",
    " nedge=2        Minimum Number of Points To Perform Edge Compensation.      ",
    "                Edge compensation is done by a relatively complicated       ",
    "                method based on the difference between each shot location   ",
    "                and the average location of all shots within avrad distance ",
    "                of that shot (the centroid location). The larger the        ",
    "                distance to the centroid, the larger the compensation.      ",
    "                The same method is used independently for the receivers.    ",
    "                Must be greater than 1.                                     ",
    "           ***  Cannot specify if avrad is -1.                              ",
    "          Note: For 2D surveys this compensation mostly affects points near ",
    "                line ends, but will also be larger for gaps within the line.",
    "                For 3D surveys this compensation mostly affects points near ",
    "                the edges of the areas of coverage, but compensation will   ",
    "                also be larger for gaps within the areas of coverage.       ",
    "      =0        Do not perform Edge Compensation. Typically, for surveys    ",
    "                with slants in the statics near edges or gaps, this results ",
    "                in values curving in the opposite direction as the slant.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cedge=avrad/10 Minimum Distance to Centroid to Perform Edge Compensation.  ",
    "                (See previous parameter). If the centroid location is       ",
    "                within this distance of the shot or receiver location,      ",
    "                do not perform edge compensation for this shot or receiver. ",
    "           ***  Cannot specify if avrad is -1. Cannot be greater than avrad.",
    "          Note: Specify this parameter in units of sloc, rloc, and cloc.    ",
    "                By default those are coordinate names so this parameter     ",
    "                is a true distance. But for 2D surveys sloc, rloc, and cloc ",
    "                can be single values such as point numbers in which case    ",
    "                this value is the difference in point numbers.              ",
    "\t\t\t\t\t\t\t\t\t     ",
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
    " --------------                                                             ",
    "                                                                            ",
    " formxy=%.20g  The C format code for printing all values to the q-records.  ",
    "              Note that the default format prints up to 20 digits           ",
    "              (but not trailing zeroes to the right of the decimal point).  ",
    "         Note: This default is used because all values in input q-files     ",
    "               can have full double precision and this program computes in  ",
    "               double precision. In particular, 2D cdp profile coordinates  ",
    "               and 3D cdp cell centre coordinates often need to be output   ",
    "               in full double precision (which is approximately 16 digits). ",
    "               In some situations this makes it hard to easily check the    ",
    "               output values so you might use %.10g here during testing.    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " dout=tstat.csv  Output csv Q-file name.                                    ",
    "                 For 2D this file contains:                                 ",
    "                    cseq,cloc1,cloc2,sstat,gstat,tstat                      ",
    "                     - where cseq,cloc1,cloc2 are the names specified in    ",
    "                       those parameters (with names set to numb1 and numb2  ",
    "                       and values of 0.0 if they are not being used).       ",
    "                     - and sstat,rstat,cstat are the statics associated     ",
    "                       with the names specified in those parameters.        ",
    "                 For 3D this file contains:                                 ",
    "                    cdp,cdpt,igi,igc,sx,sy,gx,gy,sstat,rstat,cstat          ",
    "                     - where cdp,cdpt,igi,igc are cell and transposed cell  ",
    "                       number and inline and crossline index computed from  ",
    "                       the 3D Grid defined in K-record of the input rfile.  ",
    "                     - and sx,sy,gx,gy are real-world and grid-relative     ",
    "                       cell centre coordinates computed from the 3D Grid.   ",
    "                     - and sstat,rstat,cstat are the statics associated     ",
    "                       with the names specified in those parameters.        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " --------------                                                             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The following 2 parameters affect cpu time, but not results. The search    ",
    " for the nearest shot and/or the nearest receiver is done using kdtrees.    ",
    "                                                                            ",
    " kdist=5   Initial search distance. If positive, this value is added to     ",
    "           the distance between the previous 2D cin file sorted Q-record    ",
    "           or 3D cell and its nearest shot and nearest receiver and used as ",
    "           the initial search distance to find the nearest shot and nearest ",
    "           receiver for the next 2D cin file sorted Q-record or 3D cell.    ",
    "           If negative, the initial search distance for all locations       ",
    "           is set to this absolute value.                                   ",
    " kmult=2   Search multiplier.                                               ",
    "           If the nearest shot or receiver is not found after searching with",
    "           the initial search distance, the search distance is multiplied   ",
    "           by this value, and the searches are performed again.             ",
    "           This repeats until finding the nearest shot and receiver.        ",
    "\t\t\t\t\t\t\t\t\t     ",
];

#[allow(unused_variables, unused_assignments, non_snake_case)]
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    if !isatty(0) || !isatty(1) {
        err!("**** Error: this program does not input or output traces.");
    }

    let mut avrad = 1000.0f64;
    let mut nedge = 2i32;
    let mut cedge = 100.0f64;
    let mut navrg = [0i32; 2];
    let mut nsmth = [0i32; 2];
    let mut nback = [0i32; 2];
    let mut kdist = 5.0f64;
    let mut kmult = 2.0f64;
    let mut kdadd = 1i32;

    if getpardouble("avrad", &mut avrad) != 0 {
        if avrad < 0.0 && avrad != -1.0 {
            err!("**** Error: avrad must be positive or -1.");
        }
    }
    if getparint("nedge", &mut nedge) != 0 {
        if nedge == 1 || nedge < 0 {
            err!("**** Error: nedge must greater than or equal to 2 (or 0).");
        }
        if avrad == -1.0 && nedge != 0 {
            err!("**** Error: nedge cannot be specified if avrad is -1");
        }
    }
    if getpardouble("cedge", &mut cedge) != 0 {
        if avrad == -1.0 {
            err!("**** Error: cedge cannot be specified if avrad is -1");
        }
        if cedge > avrad {
            err!("**** Error: cedge cannot be greater than avrad");
        }
        if cedge < 0.0 {
            err!("**** Error: cedge must be greater than or equal to 0.0");
        }
    } else {
        cedge = avrad / 10.0;
    }

    if countparval("navrg") > 2 {
        err!("**** Error: navrg= list can only have 1 or 2 values ");
    }
    if countparval("navrg") < 1 {
        navrg[0] = 10;
        navrg[1] = 10;
    } else {
        // The C array is uninitialized when navrg= contains no numbers.
        // This vendor binary consistently takes its first-value error path;
        // reproduce that diagnostic without reading uninitialized Rust memory.
        if getparint("navrg", &mut navrg[..]) == 0 {
            err!("**** Error in first navrg value. It must be greater or equal to -1 ");
        }
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

    if countparval("nsmth") > 2 {
        err!("**** Error: nsmth= list can only have 1 or 2 values ");
    }
    if countparval("nsmth") < 1 {
        nsmth[0] = navrg[0];
        nsmth[1] = navrg[1];
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

    if countparval("nback") > 2 {
        err!("**** Error: nback= list can only have 1 or 2 values ");
    }
    if countparval("nback") < 1 {
        nback[0] = 3;
        nback[1] = 3;
    } else {
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

    if getpardouble("kdist", &mut kdist) == 0 {
        kdist = 5.0;
    }
    if kdist == 0.0 {
        err!("**** Error: kdist cannot be 0.");
    }
    if getpardouble("kmult", &mut kmult) == 0 {
        kmult = 2.0;
    }
    if kmult < 0.0 {
        err!("**** Error: kmult must be non-negative.");
    }
    if kdist < 0.0 {
        kdist = 0.0 - kdist;
        kdadd = 0;
    }
    let sdist = kdist;
    let sdist2 = kdist;
    let rdist = kdist;
    let rdist2 = kdist;

    // Infinite extents, and all dimensions participate in Pythagorean nearest.
    let mut extent_dmin = [0.0f64; 9];
    let mut extent_dmax = [0.0f64; 9];
    let mut tree_nt = [0i32; 9];
    for i in 0..9 {
        extent_dmin[i] = -f64::MAX;
        extent_dmax[i] = f64::MAX;
        tree_nt[i] = 1;
    }

    let mut Kname = String::new();
    getparstring("rfile", &mut Kname);

    let mut maygrid = 0;
    gridcommand(&mut maygrid);
    let mut is3d = 1;
    let has_kname = countparval("rfile") > 0;
    if maygrid == 1 && has_kname {
        err!("error: input rfile not allowed when full grid on command line.");
    }
    if maygrid == -1 && !has_kname {
        err!("error: input rfile required when partial grid on command line.");
    }
    if maygrid == 0 && !has_kname {
        is3d = 0;
    }
    let mut icheck = 0;
    let mut iprint = 0;
    if getparint("check", &mut icheck) == 0 {
        icheck = 0;
    }
    if getparint("print", &mut iprint) == 0 {
        iprint = 0;
    }
    let mut gvals = [0.0; 999];
    let mut errwarn = 0;
    let mut j3dcdp = 1;
    let mut igi = 1;
    let mut igc = 1;
    let mut lastab = 0;
    let mut lastac = 0;
    let mut lastcdp;
    let mut fpK = None;
    if is3d == 1 {
        if countparval("cin") > 0 {
            err!("error: There is a 3D Grid in rfile. So cin= file cannot be input.");
        }
        if maygrid != 1 {
            fpK = fopen(&Kname, "r");
            if fpK.is_none() {
                err!("error: input rfile did not open correctly.");
            }
        }
        errwarn = 1;
        gridread(fpK.as_ref(), &mut gvals, &mut errwarn);
        if errwarn > 0 {
            err!("error reading grid (from rfile or command line parameters)");
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
        gridiccdp(&gvals, 1, 1, &mut j3dcdp);
        lastab = (gvals[12] + 0.001) as i32;
        lastac = (gvals[13] + 0.001) as i32;
        lastcdp = lastab * lastac;
        if nback[0] >= lastab {
            err!("error: 3D Grid inlines do not have enough cells for the nback parameter.");
        }
        if nback[1] >= lastac {
            err!("error: 3D Grid crosslines do not have enough cells for the nback parameter.");
        }
    } else {
        lastcdp = 0;
        if countparval("cin") < 1 {
            err!("error: No 3D Grid in rfile. So cin= file MUST be input.");
        }
    }
    let mut formxy = String::new();
    if getparstring("formxy", &mut formxy) == 0 {
        formxy = "%.20g".to_owned();
    }
    let formxylong = format!(",{}", formxy);
    let mut dout = String::new();
    let fpD = if countparval("dout") > 0 {
        getparstring("dout", &mut dout);
        fopen(&dout, "w")
    } else {
        fopen("tstat.csv", "w")
    };
    let Some(fpD) = fpD else {
        err!("error: output floating datum file did not open correctly.");
    };

    // Read parameters and non-tuple values for SInfo.
    let mut sstat = String::new();
    if getparstring("sstat", &mut sstat) == 0 {
        sstat = "sstat".to_owned();
    }
    let mut stree_numd = 1usize;
    let mut Sname = String::new();
    if countparval("sin") > 0 {
        getparstring("sin", &mut Sname);
        if Sname == "none" {
            stree_numd = 0;
        }
    } else {
        Sname = "sstat.csv".to_owned();
    }
    let mut skey = Vec::new();
    let mut slocn = [0usize; 2];
    let mut osloc = 0usize;
    let mut smult = 1.0;
    let mut SInfo: Vec<QInfo> = Vec::new();
    let mut stree_dl = Vec::new();
    let mut stree_nodes = Vec::new();
    if stree_numd > 0 {
        let numd = countparval("sloc");
        if numd < 1 {
            stree_numd = 2;
            skey = vec!["sx".to_owned(), "sy".to_owned()];
        } else {
            if numd > 2 {
                err!("**** Error for sin parameters: Only 1 or 2 sloc names can be specified.");
            }
            stree_numd = numd as usize;
            getparstringarray("sloc", &mut skey);
        }
        if getpardouble("smult", &mut smult) == 0 {
            smult = 1.0;
        }
        if smult == 0.0 {
            err!("**** Error for sin parameters: smult cannot be 0.0");
        }
        let Some(fpS) = fopen(&Sname, "r") else {
            err!("error sin file did not open correctly.");
        };
        let mut pname = skey.clone();
        pname.push(sstat.clone());
        let result = getviaqfile(&fpS, &pname, pname.len() as i32, 0);
        let q = match result {
            Ok(q) => q,
            Err(errwarn) => {
                if errwarn == 1 {
                    err!("sin file getqinfo error: extra C_SU_NAMES record in q-file");
                } else if errwarn == 2 {
                    err!("sin file getqinfo error: extra C_SU_NDIMS record in q-file");
                } else if errwarn == 3 {
                    err!("sin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
                } else if errwarn == 11 {
                    err!("sin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
                } else if errwarn == 12 {
                    err!("sin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
                } else if errwarn == 22 {
                    err!("sin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
                } else if errwarn == 23 {
                    err!(
                        "sin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"
                    );
                } else if errwarn == 24 {
                    err!(
                        "sin file getviaqfile error: C_SU_NDIMS blank where valid number expected"
                    );
                } else if errwarn == 25 {
                    err!("sin file getviaqfile error: C_SU_NDIMS non-number where valid number expected");
                } else if errwarn == 26 {
                    err!("sin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
                } else if errwarn == 27 {
                    err!("sin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
                } else if errwarn > 100 {
                    err!("sin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
                } else if errwarn > 0 {
                    err!(
                        "sin file getviaqfile error: unrecognized error code %d",
                        errwarn
                    );
                }

                unreachable!()
            }
        };
        if q.ifixd == 0 {
            err!("sin file error: input with varying number of tuples is not allowed.");
        }
        for i in 0..stree_numd {
            let mut found = None;
            for j in 0..q.iztuple as usize {
                if q.names[j] == skey[i] {
                    found = Some(j);
                }
            }
            let Some(j) = found else {
                err!(
                    "sin file error: sloc name %s not found in non-tuple part of sin Q-file.",
                    skey[i].as_str()
                );
            };
            slocn[i] = j;
        }
        let mut found = None;
        for j in 0..q.iztuple as usize {
            if q.names[j] == sstat {
                found = Some(j);
            }
        }
        let Some(j) = found else {
            err!("sin file error: Q-file must have sstat parameter value %s (among non-tuple names).", sstat.as_str());
        };
        osloc = j;
        SInfo = q.records;
        if SInfo.is_empty() {
            err!("sin file error: no records starting with Q found in file.");
        }
        stree_nodes = vec![Node::default(); SInfo.len()];
        for i in 0..stree_numd {
            let mut d = vec![0.0; SInfo.len()];
            for j in 0..SInfo.len() {
                d[j] = SInfo[j].dlots[slocn[i]];
            }
            stree_dl.push(d);
        }
        connect_nodes(&mut stree_nodes, SInfo.len(), &stree_dl, stree_numd, 1);
    }

    // Read parameters and non-tuple values for RInfo.
    let mut rstat = String::new();
    if getparstring("rstat", &mut rstat) == 0 {
        rstat = "gstat".to_owned();
    }
    let mut rtree_numd = 1usize;
    let mut Rname = String::new();
    if countparval("rin") > 0 {
        getparstring("rin", &mut Rname);
        if Rname == "none" {
            rtree_numd = 0;
        }
    } else {
        Rname = "rstat.csv".to_owned();
    }
    let mut rkey = Vec::new();
    let mut rlocn = [0usize; 2];
    let mut orloc = 0usize;
    let mut rmult = 1.0;
    let mut RInfo: Vec<QInfo> = Vec::new();
    let mut rtree_dl = Vec::new();
    let mut rtree_nodes = Vec::new();
    if rtree_numd > 0 {
        let numd = countparval("rloc");
        if numd < 1 {
            rtree_numd = 2;
            rkey = vec!["gx".to_owned(), "gy".to_owned()];
        } else {
            if numd > 2 {
                err!("**** Error for rin parameters: Only 1 or 2 rloc names can be specified.");
            }
            rtree_numd = numd as usize;
            getparstringarray("rloc", &mut rkey);
        }
        if getpardouble("rmult", &mut rmult) == 0 {
            rmult = 1.0;
        }
        if rmult == 0.0 {
            err!("**** Error for rin parameters: rmult cannot be 0.0");
        }
        let Some(fpR) = fopen(&Rname, "r") else {
            err!("error: rin file did not open correctly.");
        };
        let mut pname = rkey.clone();
        pname.push(rstat.clone());
        let result = getviaqfile(&fpR, &pname, pname.len() as i32, 0);
        let q = match result {
            Ok(q) => q,
            Err(errwarn) => {
                if errwarn == 1 {
                    err!("rin file getqinfo error: extra C_SU_NAMES record in q-file");
                } else if errwarn == 2 {
                    err!("rin file getqinfo error: extra C_SU_NDIMS record in q-file");
                } else if errwarn == 3 {
                    err!("rin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
                } else if errwarn == 11 {
                    err!("rin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
                } else if errwarn == 12 {
                    err!("rin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
                } else if errwarn == 22 {
                    err!("rin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
                } else if errwarn == 23 {
                    err!(
                        "rin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"
                    );
                } else if errwarn == 24 {
                    err!(
                        "rin file getviaqfile error: C_SU_NDIMS blank where valid number expected"
                    );
                } else if errwarn == 25 {
                    err!("rin file getviaqfile error: C_SU_NDIMS non-number where valid number expected");
                } else if errwarn == 26 {
                    err!("rin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
                } else if errwarn == 27 {
                    err!("rin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
                } else if errwarn > 100 {
                    err!("rin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
                } else if errwarn > 0 {
                    err!(
                        "rin file getviaqfile error: unrecognized error code %d",
                        errwarn
                    );
                }

                unreachable!()
            }
        };
        if q.ifixd == 0 {
            err!("rin file error: input with varying number of tuples is not allowed.");
        }
        for i in 0..rtree_numd {
            let mut found = None;
            for j in 0..q.iztuple as usize {
                if q.names[j] == rkey[i] {
                    found = Some(j);
                }
            }
            let Some(j) = found else {
                err!(
                    "rin file error: rloc %s not found in non-tuple part of rin Q-file.",
                    rkey[i].as_str()
                );
            };
            rlocn[i] = j;
        }
        let mut found = None;
        for j in 0..q.iztuple as usize {
            if q.names[j] == rstat {
                found = Some(j);
            }
        }
        let Some(j) = found else {
            err!(
                "error: rin Q-file must have rstat parameter value %s (among non-tuple names).",
                rstat.as_str()
            );
        };
        orloc = j;
        RInfo = q.records;
        if RInfo.is_empty() {
            err!("rin file error: no records starting with Q found in file.");
        }
        rtree_nodes = vec![Node::default(); RInfo.len()];
        for i in 0..rtree_numd {
            let mut d = vec![0.0; RInfo.len()];
            for j in 0..RInfo.len() {
                d[j] = RInfo[j].dlots[rlocn[i]];
            }
            rtree_dl.push(d);
        }
        connect_nodes(&mut rtree_nodes, RInfo.len(), &rtree_dl, rtree_numd, 1);
    }
    // Read CDP locations and sequencing values.
    let mut ctree_numd = countparval("cloc") as usize;
    let mut ckey = Vec::new();
    if ctree_numd < 1 {
        ctree_numd = 2;
        ckey = vec!["gx".to_owned(), "gy".to_owned()];
    } else {
        if ctree_numd > 2 {
            err!("**** Error for cin parameters: Only 1 or 2 cloc names can be specified.");
        }
        getparstringarray("cloc", &mut ckey);
        if ctree_numd == 1 {
            ckey.push("numb2".to_owned());
        }
    }
    let mut cstat = String::new();
    if getparstring("cstat", &mut cstat) == 0 {
        cstat = "tstat".to_owned();
    }
    let mut cseq = String::new();
    let mut CInfo: Vec<QInfo> = Vec::new();
    let mut clocn = [0usize; 2];
    let mut locp = None;
    if is3d != 1 {
        let mut Cname = String::new();
        getparstring("cin", &mut Cname);
        let Some(fpC) = fopen(&Cname, "r") else {
            err!("error: cin file did not open correctly.");
        };
        if countparval("cseq") > 0 {
            getparstring("cseq", &mut cseq);
        } else {
            cseq = "cdp".to_owned();
        }
        let mut pname = Vec::new();
        for i in 0..ctree_numd {
            pname.push(ckey[i].clone());
        }
        pname.push(cstat.clone());
        if cseq != "asis" {
            pname.push(cseq.clone());
        }
        let result = getviaqfile(&fpC, &pname, pname.len() as i32, 0);
        let q = match result {
            Ok(q) => q,
            Err(errwarn) => {
                if errwarn == 1 {
                    err!("cin file getqinfo error: extra C_SU_NAMES record in q-file");
                } else if errwarn == 2 {
                    err!("cin file getqinfo error: extra C_SU_NDIMS record in q-file");
                } else if errwarn == 3 {
                    err!("cin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record.");
                } else if errwarn == 11 {
                    err!("cin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES");
                } else if errwarn == 12 {
                    err!("cin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record.");
                } else if errwarn == 22 {
                    err!("cin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record.");
                } else if errwarn == 23 {
                    err!(
                        "cin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"
                    );
                } else if errwarn == 24 {
                    err!(
                        "cin file getviaqfile error: C_SU_NDIMS blank where valid number expected"
                    );
                } else if errwarn == 25 {
                    err!("cin file getviaqfile error: C_SU_NDIMS non-number where valid number expected");
                } else if errwarn == 26 {
                    err!("cin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple");
                } else if errwarn == 27 {
                    err!("cin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found.");
                } else if errwarn > 100 {
                    err!("cin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)",errwarn-99);
                } else if errwarn > 0 {
                    err!(
                        "cin file getviaqfile error: unrecognized error code %d",
                        errwarn
                    );
                }

                unreachable!()
            }
        };
        if q.ifixd == 0 {
            err!("cin file error: input with varying number of tuples is not allowed.");
        }
        for i in 0..ctree_numd {
            let mut found = None;
            for j in 0..q.iztuple as usize {
                if q.names[j] == ckey[i] {
                    found = Some(j);
                }
            }
            let Some(j) = found else {
                err!(
                    "cin file error: cloc %s not found in non-tuple part of cin Q-file.",
                    ckey[i].as_str()
                );
            };
            clocn[i] = j;
        }
        CInfo = q.records;
        if CInfo.is_empty() {
            err!("cin file error: no records starting with Q found in file.");
        }
        lastcdp = CInfo.len() as i32;
        if nback[0] >= lastcdp {
            err!("error: cin file does not have enough Q-records for the nback parameter.");
        }
        if cseq == "asis" {
            cseq = "numb1".to_owned();
        } else {
            for i in 0..q.iztuple as usize {
                if q.names[i] == cseq {
                    locp = Some(i);
                }
            }
            let ukey = format!("_{}_", cseq);
            for i in 0..q.iztuple as usize {
                if q.names[i].contains(&ukey) {
                    locp = Some(i);
                }
            }
            let Some(p) = locp else {
                err!(
                    "error: input cin Q-file must have your cseq=%s (among non-tuple names).",
                    cseq.as_str()
                );
            };
            CInfo.sort_by(|a, b| comp_sort1(&a.dlots, &b.dlots, p));
            for j in 1..lastcdp as usize {
                if CInfo[j - 1].dlots[p] == CInfo[j].dlots[p] {
                    err!(
                        "error: Two sets of values input for cin cseq = %f",
                        CInfo[j].dlots[p]
                    );
                }
            }
        }
    }
    if stree_numd != 0 && stree_numd != ctree_numd {
        err!("**** Error: Must have the same number of sloc as cloc names.");
    }
    if rtree_numd != 0 && rtree_numd != ctree_numd {
        err!("**** Error: Must have the same number of rloc as cloc names.");
    }
    // Duplicate the single coordinate for the circle/centroid calculations only.
    if ctree_numd == 1 {
        clocn[1] = clocn[0];
        slocn[1] = slocn[0];
        rlocn[1] = rlocn[0];
        avrad *= 2.0f64.sqrt();
    }
    let mut out_elem = vec![
        0usize;
        if SInfo.len() > RInfo.len() {
            SInfo.len()
        } else {
            RInfo.len()
        }
    ];
    let mut isize = 2 * navrg[0] + 3;
    if is3d == 1 && 2 * navrg[1] + 3 > isize {
        isize = 2 * navrg[1] + 3;
    }
    if 2 * nsmth[0] + 1 > isize {
        isize = 2 * nsmth[0] + 1;
    }
    if is3d == 1 && 2 * nsmth[1] + 1 > isize {
        isize = 2 * nsmth[1] + 1;
    }
    let klast = 2;
    let mut dalla = vec![vec![0.0; klast as usize]; (lastcdp + isize) as usize];
    let mut dtog = vec![0.0; klast as usize];
    let mut dall3 = Vec::new();
    if is3d == 1 {
        let n = if lastab > lastac { lastab } else { lastac };
        dall3 = vec![vec![0.0; klast as usize]; (n + isize) as usize];
    }
    checkpars();
    let mut target = [0.0; 9];
    let mut extent_min = [0.0; 9];
    let mut extent_max = [0.0; 9];
    let mut snear_dist: f64 = 0.0;
    let mut rnear_dist: f64 = 0.0;
    let mut sdist2 = sdist2;
    let mut rdist2 = rdist2;
    let mut near_elem = 0;
    let mut num_out = 0;
    let mut snum_found = 0;
    let mut rnum_found = 0;
    let mut kcycles = 0;
    for j in 0..lastcdp as usize {
        if is3d == 1 {
            gridcdpic(&gvals, j3dcdp + j as i32, &mut igi, &mut igc);
            let (mut x, mut y) = (0.0, 0.0);
            gridicrawxy(&gvals, igi, igc, &mut x, &mut y);
            target[0] = x;
            target[1] = y;
        } else {
            target[0] = CInfo[j].dlots[clocn[0]];
            target[1] = CInfo[j].dlots[clocn[1]];
        }

        if stree_numd > 0 {
            if kdadd == 1 {
                sdist2 = sdist + snear_dist;
            }
            cycle_for_near(
                &stree_nodes,
                &stree_dl,
                &tree_nt,
                stree_numd,
                &extent_dmin,
                &extent_dmax,
                &target,
                sdist2,
                kmult,
                &mut near_elem,
                &mut snear_dist,
                &mut snum_found,
                &mut kcycles,
            );
            snear_dist = snear_dist.sqrt();
            if avrad < 0.0 {
                dalla[j][0] = SInfo[near_elem].dlots[osloc];
            } else {
                extent_min[0] = target[0] - snear_dist - avrad;
                extent_max[0] = target[0] + snear_dist + avrad;
                extent_min[1] = target[1] - snear_dist - avrad;
                extent_max[1] = target[1] + snear_dist + avrad;
                find_in(
                    &stree_nodes,
                    &stree_dl,
                    stree_numd,
                    &extent_min,
                    &extent_max,
                    &mut out_elem,
                    &mut num_out,
                );
                dalla[j][0] = 0.0;
                let mut xc = 0.0;
                let mut yc = 0.0;
                let mut k = 0usize;
                for i in 0..num_out {
                    if (target[0] - SInfo[out_elem[i]].dlots[slocn[0]])
                        * (target[0] - SInfo[out_elem[i]].dlots[slocn[0]])
                        + (target[1] - SInfo[out_elem[i]].dlots[slocn[1]])
                            * (target[1] - SInfo[out_elem[i]].dlots[slocn[1]])
                        <= (snear_dist + avrad) * (snear_dist + avrad)
                    {
                        dalla[j][0] += SInfo[out_elem[i]].dlots[osloc];
                        xc += SInfo[out_elem[i]].dlots[slocn[0]];
                        yc += SInfo[out_elem[i]].dlots[slocn[1]];
                        out_elem[k] = out_elem[i];
                        k += 1;
                    }
                }
                if k > 1 {
                    dalla[j][0] /= k as f64;
                    xc /= k as f64;
                    yc /= k as f64;
                }
                let mut dx = target[0] - xc;
                let mut dy = target[1] - yc;
                let mut dr = (dx * dx + dy * dy).sqrt();
                if k as i32 >= nedge && dr >= cedge {
                    let mut ds = dy / dr;
                    let mut dc = dx / dr;
                    let (mut xl, mut yl, mut sl, mut kl) = (0.0, 0.0, 0.0, 0);
                    let (mut xh, mut yh, mut sh, mut kh) = (0.0, 0.0, 0.0, 0);
                    let mut tx;
                    for i in 0..k {
                        tx = (xc - SInfo[out_elem[i]].dlots[slocn[0]]) * dc
                            + (yc - SInfo[out_elem[i]].dlots[slocn[1]]) * ds;
                        if tx < 0.0 {
                            xl += SInfo[out_elem[i]].dlots[slocn[0]];
                            yl += SInfo[out_elem[i]].dlots[slocn[1]];
                            sl += SInfo[out_elem[i]].dlots[osloc];
                            kl += 1;
                        } else {
                            xh += SInfo[out_elem[i]].dlots[slocn[0]];
                            yh += SInfo[out_elem[i]].dlots[slocn[1]];
                            sh += SInfo[out_elem[i]].dlots[osloc];
                            kh += 1;
                        }
                    }
                    if kl > 0 && kh > 0 && (xh != xl || yh != yl) {
                        xl /= kl as f64;
                        yl /= kl as f64;
                        sl /= kl as f64;
                        xh /= kh as f64;
                        yh /= kh as f64;
                        sh /= kh as f64;
                        dx = xh - xl;
                        dy = yh - yl;
                        dr = (dx * dx + dy * dy).sqrt();
                        ds = dy / dr;
                        dc = dx / dr;
                        tx = (xh - target[0]) * dc + (yh - target[1]) * ds;
                        dalla[j][0] = tx / dr * sl - (tx / dr - 1.0) * sh;
                    }
                }
            }
        }

        if rtree_numd > 0 {
            if kdadd == 1 {
                rdist2 = rdist + rnear_dist;
            }
            cycle_for_near(
                &rtree_nodes,
                &rtree_dl,
                &tree_nt,
                rtree_numd,
                &extent_dmin,
                &extent_dmax,
                &target,
                rdist2,
                kmult,
                &mut near_elem,
                &mut rnear_dist,
                &mut rnum_found,
                &mut kcycles,
            );
            rnear_dist = rnear_dist.sqrt();
            if avrad < 0.0 {
                dalla[j][1] = RInfo[near_elem].dlots[orloc];
            } else {
                extent_min[0] = target[0] - rnear_dist - avrad;
                extent_max[0] = target[0] + rnear_dist + avrad;
                extent_min[1] = target[1] - rnear_dist - avrad;
                extent_max[1] = target[1] + rnear_dist + avrad;
                find_in(
                    &rtree_nodes,
                    &rtree_dl,
                    rtree_numd,
                    &extent_min,
                    &extent_max,
                    &mut out_elem,
                    &mut num_out,
                );
                dalla[j][1] = 0.0;
                let mut xc = 0.0;
                let mut yc = 0.0;
                let mut k = 0usize;
                for i in 0..num_out {
                    if (target[0] - RInfo[out_elem[i]].dlots[rlocn[0]])
                        * (target[0] - RInfo[out_elem[i]].dlots[rlocn[0]])
                        + (target[1] - RInfo[out_elem[i]].dlots[rlocn[1]])
                            * (target[1] - RInfo[out_elem[i]].dlots[rlocn[1]])
                        <= (rnear_dist + avrad) * (rnear_dist + avrad)
                    {
                        dalla[j][1] += RInfo[out_elem[i]].dlots[orloc];
                        xc += RInfo[out_elem[i]].dlots[rlocn[0]];
                        yc += RInfo[out_elem[i]].dlots[rlocn[1]];
                        out_elem[k] = out_elem[i];
                        k += 1;
                    }
                }
                if k > 1 {
                    dalla[j][1] /= k as f64;
                    xc /= k as f64;
                    yc /= k as f64;
                }
                let mut dx = target[0] - xc;
                let mut dy = target[1] - yc;
                let mut dr = (dx * dx + dy * dy).sqrt();
                if k as i32 >= nedge && dr >= cedge {
                    let mut ds = dy / dr;
                    let mut dc = dx / dr;
                    let (mut xl, mut yl, mut sl, mut kl) = (0.0, 0.0, 0.0, 0);
                    let (mut xh, mut yh, mut sh, mut kh) = (0.0, 0.0, 0.0, 0);
                    let mut tx;
                    for i in 0..k {
                        tx = (xc - RInfo[out_elem[i]].dlots[rlocn[0]]) * dc
                            + (yc - RInfo[out_elem[i]].dlots[rlocn[1]]) * ds;
                        if tx < 0.0 {
                            xl += RInfo[out_elem[i]].dlots[rlocn[0]];
                            yl += RInfo[out_elem[i]].dlots[rlocn[1]];
                            sl += RInfo[out_elem[i]].dlots[orloc];
                            kl += 1;
                        } else {
                            xh += RInfo[out_elem[i]].dlots[rlocn[0]];
                            yh += RInfo[out_elem[i]].dlots[rlocn[1]];
                            sh += RInfo[out_elem[i]].dlots[orloc];
                            kh += 1;
                        }
                    }
                    if kl > 0 && kh > 0 && (xh != xl || yh != yl) {
                        xl /= kl as f64;
                        yl /= kl as f64;
                        sl /= kl as f64;
                        xh /= kh as f64;
                        yh /= kh as f64;
                        sh /= kh as f64;
                        dx = xh - xl;
                        dy = yh - yl;
                        dr = (dx * dx + dy * dy).sqrt();
                        ds = dy / dr;
                        dc = dx / dr;
                        tx = (xh - target[0]) * dc + (yh - target[1]) * ds;
                        dalla[j][1] = tx / dr * sl - (tx / dr - 1.0) * sh;
                    }
                }
            }
        }
    } // end loop over CDPs
      // Perform averaging and smoothing, retaining the C pass order.
    let mut ierr = 0;
    if is3d != 1 {
        if navrg[0] > 0 {
            runsmo(
                &mut dalla,
                lastcdp + 2 * navrg[0] + 1,
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
            for i in 0..klast as usize {
                dtog[i] = 0.0;
            }
            for j in 0..lastcdp as usize {
                for i in 0..klast as usize {
                    dtog[i] += dalla[j][i];
                }
            }
            for j in 0..lastcdp as usize {
                for i in 0..klast as usize {
                    dalla[j][i] = dtog[i] / lastcdp as f64;
                }
            }
        }
        if nsmth[0] > 0 {
            runsmo(
                &mut dalla,
                lastcdp + 2 * nsmth[0] + 1,
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
    } else {
        if navrg[0] != 0 || nsmth[0] > 0 {
            for n in 0..lastac as usize {
                for j in 0..lastab as usize {
                    for i in 0..klast as usize {
                        dall3[j][i] = dalla[j + n * lastab as usize][i];
                    }
                }
                if navrg[0] > 0 {
                    runsmo(
                        &mut dall3,
                        lastab + 2 * navrg[0] + 1,
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
                    for i in 0..klast as usize {
                        dtog[i] = 0.0;
                    }
                    for j in 0..lastab as usize {
                        for i in 0..klast as usize {
                            dtog[i] += dall3[j][i];
                        }
                    }
                    for j in 0..lastab as usize {
                        for i in 0..klast as usize {
                            dall3[j][i] = dtog[i] / lastab as f64;
                        }
                    }
                }
                if nsmth[0] > 0 {
                    runsmo(
                        &mut dall3,
                        lastab + 2 * nsmth[0] + 1,
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
                for j in 0..lastab as usize {
                    for i in 0..klast as usize {
                        dalla[j + n * lastab as usize][i] = dall3[j][i];
                    }
                }
            }
        }
        // Repeat in the crossline direction. The C average-all branch follows
        // smoothing's else, so it is deliberately conditional on nsmth[1]==0.
        if navrg[1] != 0 || nsmth[1] > 0 {
            for n in 0..lastab as usize {
                for j in 0..lastac as usize {
                    for i in 0..klast as usize {
                        dall3[j][i] = dalla[n + j * lastab as usize][i];
                    }
                }
                if navrg[1] > 0 {
                    runsmo(
                        &mut dall3,
                        lastac + 2 * navrg[1] + 1,
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
                        &mut dall3,
                        lastac + 2 * nsmth[1] + 1,
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
                    for i in 0..klast as usize {
                        dtog[i] = 0.0;
                    }
                    for j in 0..lastac as usize {
                        for i in 0..klast as usize {
                            dtog[i] += dall3[j][i];
                        }
                    }
                    for j in 0..lastac as usize {
                        for i in 0..klast as usize {
                            dall3[j][i] = dtog[i] / lastac as f64;
                        }
                    }
                }
                for j in 0..lastac as usize {
                    for i in 0..klast as usize {
                        dalla[n + j * lastab as usize][i] = dall3[j][i];
                    }
                }
            }
        }
    }
    if stree_numd < 1 {
        smult = 0.0;
    }
    if rtree_numd < 1 {
        rmult = 0.0;
    }
    let mut cmult = 1.0;
    if stree_numd > 0 && rtree_numd > 0 {
        cmult = 0.5;
    }
    let mut d4print = [0.0; 11];
    if is3d != 1 {
        fprintf!(
            &fpD,
            "C_SU_MATCH,%s,,,,,\nC_SU_SETID,Q,,,,,\nC_SU_FORMS\nC_SU_ID",
            cseq.as_str()
        );
        for _i in 0..6 {
            fprintf!(&fpD, ",%s", formxy.as_str());
        }
        fprintf!(
            &fpD,
            "\nC_SU_NAMES,,,,,,\nC_SU_ID,%s,%s,%s,%s,%s,%s\n",
            cseq.as_str(),
            ckey[0].as_str(),
            ckey[1].as_str(),
            sstat.as_str(),
            rstat.as_str(),
            cstat.as_str()
        );
        for j in 0..lastcdp as usize {
            if let Some(p) = locp {
                d4print[0] = CInfo[j].dlots[p];
            } else {
                d4print[0] = (j + 1) as f64;
            }
            d4print[1] = CInfo[j].dlots[clocn[0]];
            if ctree_numd == 2 {
                d4print[2] = CInfo[j].dlots[clocn[1]];
            } else {
                d4print[2] = 0.0;
            }
            d4print[3] = dalla[j][0] * smult;
            d4print[4] = dalla[j][1] * rmult;
            d4print[5] = (dalla[j][0] * smult + dalla[j][1] * rmult) * cmult;
            fprintf!(&fpD, "Q");
            for i in 0..6 {
                fprintf!(&fpD, formxylong.as_str(), d4print[i]);
            }
            fprintf!(&fpD, "\n");
        }
    } else {
        fprintf!(
            &fpD,
            "C_SU_MATCH,cdp,,,,,,\nC_SU_SETID,Q,,,,,,\nC_SU_FORMS\nC_SU_ID"
        );
        for _i in 0..10 {
            fprintf!(&fpD, ",%s", formxy.as_str());
        }
        fprintf!(
            &fpD,
            "\nC_SU_NAMES,,,,,,,,,,,\nC_SU_ID,cdp,cdpt,igi,igc,sx,sy,gx,gy,%s,%s,%s\n",
            sstat.as_str(),
            rstat.as_str(),
            cstat.as_str()
        );
        for j in 0..lastcdp as usize {
            gridcdpic(&gvals, j3dcdp + j as i32, &mut igi, &mut igc);
            let mut n = 0;
            gridiccdp90(&gvals, igi, igc, &mut n);
            let (mut xc, mut yc, mut dx, mut dy) = (0.0, 0.0, 0.0, 0.0);
            gridicrawxy(&gvals, igi, igc, &mut xc, &mut yc);
            gridicgridxy(&gvals, igi, igc, &mut dx, &mut dy);
            d4print[0] = (j3dcdp + j as i32) as f64;
            d4print[1] = n as f64;
            d4print[2] = igi as f64;
            d4print[3] = igc as f64;
            d4print[4] = xc;
            d4print[5] = yc;
            d4print[6] = dx;
            d4print[7] = dy;
            d4print[8] = dalla[j][0] * smult;
            d4print[9] = dalla[j][1] * rmult;
            d4print[10] = (dalla[j][0] * smult + dalla[j][1] * rmult) * cmult;
            fprintf!(&fpD, "Q");
            for i in 0..11 {
                fprintf!(&fpD, formxylong.as_str(), d4print[i]);
            }
            fprintf!(&fpD, "\n");
        }
    }
    exit(CWP_Exit());
}
