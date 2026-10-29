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

//! SUSHIFTCSV - Apply Static Shifts From Nearest Q-file Records.
//! Port of su/main/statics/sushiftcsv.c (Revision 1.01).
//! Created: July 2023: Andre Latour. This program started from sunearcsv.

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUSHIFTCSV - Apply Static Shifts From Nearest Q-file Records.              ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sushiftcsv [parameters].                                                  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " sin=sstat.csv  Input csv Q-file containing Source Information.             ",
    "    =none       No sin file. Value 0 is used for all source statics.        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " skeyloc=fldr   List of Source Location Key Names. At least 1 is needed.    ",
    "                Must be a key name and exist in the sin Q-file.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " sstat=sstat    This is the default source static name. This name must be   ",
    "                in the sin Q-file, but does not have to be a key name.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " smult=1.0      Source Static Multiplier. Statics in Q-files should be      ",
    "                milliseconds, if not, you can scale them before use here.   ",
    "          Note: Also, to reverse static application, specify negative here. ",
    "                                                                            ",
    " slimit=        Source Location Limit. This is the maximum Pythagorean      ",
    "                difference between skeyloc values from trace and nearest    ",
    "                sin Q-file point for trace to use source static, otherwise  ",
    "                a source static of zero is used for this trace.             ",
    "           ***  If not specified, there is no limit, so the source static   ",
    "                value comes from whichever sin Q-file point is nearest.     ",
    "           ***  If you wish to only use the static at exact locations then  ",
    "                a value such as 0.1 is recommended here because trace keys  ",
    "                are integer but Q-file values are double precision float.   ",
    "                (Float values may not be exact whole numbers).              ",
    "                                                                            ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rin=rstat.csv  Input csv Q-file containing Receiver Information.           ",
    "    =none       No rin file. Value 0 is used for all receiver statics.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rkeyloc=gaps   List of Receiver Location Key Names. At least 1 is needed.  ",
    "                Must be a key name and exist in the rin Q-file.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rstat=gstat    This is the default receiver static name. This name must be ",
    "                in the rin Q-file, but does not have to be a key name.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " rmult=1.0      Receiver Static Multiplier. Statics in Q-files should be    ",
    "                milliseconds, if not, you can scale them before use here.   ",
    "          Note: Also, to reverse static application, specify negative here. ",
    "                                                                            ",
    " rlimit=        Receiver Location Limit. This is the maximum Pythagorean    ",
    "                difference between rkeyloc values from trace and nearest    ",
    "                rin Q-file point for trace to use receiver static, otherwise",
    "                a receiver static of zero is used for this trace.           ",
    "           ***  If not specified, there is no limit, so the receiver static ",
    "                value comes from whichever rin Q-file point is nearest.     ",
    "                                                                            ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cin=none       Input csv Q-file containing CDP Information.                ",
    "    =none       Default. No cin file. Value 0 is used for all CDP statics.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " ckeyloc=cdp    List of CDP Location Key Names. At least 1 is needed.       ",
    "                Must be a key name and exist in the cin Q-file.             ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cstat=tstat    This is the default CDP static name. This name must be      ",
    "                in the cin Q-file, but does not have to be a key name.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " cmult=-2.0     CDP Static Multiplier.                                      ",
    "           ***  The default here is negative 2. ***                         ",
    "                By default, using all 3 files means the applied static is:  ",
    "                    sstat + rstat - 2*cstat                                 ",
    "                This corresponds to a common situation when CDP datum       ",
    "                statics are used (full statics in sin and rin Q-files and   ",
    "                smoothed and/or averaged statics in cin Q-file).            ",
    "\t\t\t\t\t\t\t\t\t     ",
    " climit=        CDP Location Limit. This is the maximum Pythagorean         ",
    "                difference between ckeyloc values from trace and nearest    ",
    "                cin Q-file point for trace to use CDP static, otherwise     ",
    "                a CDP static of zero is used for this trace.                ",
    "           ***  If not specified, there is no limit, so the CDP static      ",
    "                value comes from whichever cin Q-file point is nearest.     ",
    "                                                                            ",
    "                                                                            ",
    " SIGN CONVENTION:                                                           ",
    "      ***       A negative static value results in shifting seismic values  ",
    "                towards the beginning of traces.                            ",
    "                For example -18.3 will shift a peak at 1000 to 981.7        ",
    "                                                                            ",
    " OTHER INFORMATION:                                                         ",
    "                This program uses the static values FROM NEAREST LOCATIONS. ",
    "                This is the Pythagorean nearest record in the Q-files       ",
    "                considering skeyloc,rkeyloc,ckeyloc values just as numbers. ",
    "                It is deliberately coded to be the nearest location rather  ",
    "                than the exact location. This is done to simplify certain   ",
    "                static situations (beyond the scope of this document).      ",
    "         ***    In usual useage, the sin,rin,cin files contain all locations",
    "                and therefor nearest location is also the exact location.   ",
    "         ***    If you specify x,y keys for any locations the scalco flag   ",
    "                in the trace header is used to scale the trace x,y but      ",
    "                Q-file standards require they contain unscaled values.      ",
    "                Similary, trace scalel is applied to elevations and related ",
    "                values from the trace header, but not to the Q-file values. ",
    "\t\t\t\t\t\t\t\t\t     ",
    "         ***    This program does not use the input header statics at all.  ",
    "         ***    This program does not change output header statics at all.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   - - - - - - - - - -                                                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The following 2 parameters affect cpu time, but not results. The search    ",
    " for the nearest source and/or receiver and/or cdp is done using kdtrees    ",
    " built for skeyloc,rkeyloc,ckeyloc values from the sin,rin,cin records.     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " kdist=5   Initial search distance. If positive, this value is added to the ",
    "           distance between the previous trace locations and their nearest  ",
    "           kdtree locations and used as the initial search distance to find ",
    "           the kdtree locations of the current trace locations.             ",
    "           If negative, the initial search distance for all trace locations ",
    "           is set to this absolute value.                                   ",
    " kmult=2   Search multiplier.                                               ",
    "           If the nearest kdtree location is not found after searching with ",
    "           the initial search distance, the search distance is multiplied   ",
    "           by this value, and the kdtree searches are performed again.      ",
    "           This repeats until finding the nearest locations in the kdtrees. ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "   ------------------------------------------------------------------       ",
    "\t\t\t\t\t\t\t\t\t     ",
];

use oxwave_core::prelude::*;
use oxwave_header_case::{fromhead, GetCase};
use oxwave_q_file::{getviaqfile, QFile};
use oxwave_sinc_interp::ints8r;

struct Statics {
    kase: Vec<i32>,
    info: QFile,
    oloc: usize,
    mult: f64,
    limit: f64,
    tree_dl: Vec<Vec<f64>>,
    nodes: Vec<Node>,
    near_dist: f64,
}

fn setup_s() -> Option<Statics> {
    let mut filename = String::new();
    if countparval("sin") > 0 {
        getparstring("sin", &mut filename);
    } else {
        filename = "sstat.csv".to_string();
    }
    if filename == "none" {
        return None;
    }
    let mut key = Vec::new();
    let numd = countparval("skeyloc");
    if numd < 1 {
        key.push("fldr".to_string());
    } else {
        getparstringarray("skeyloc", &mut key);
        if numd > 9 {
            err!("**** Error for sin parameters: Maximum of 9 location names can be specified.");
        }
    }
    let mut kase = Vec::new();
    for name in &key {
        let k = match name.as_str() {
            "sx" => 101,
            "sy" => 102,
            "gx" => 103,
            "gy" => 104,
            "selev" => 111,
            "gelev" => 112,
            "sdepth" => 113,
            "sdel" => 114,
            "gdel" => 115,
            "swdep" => 116,
            "gwdep" => 117,
            _ => {
                let k = GetCase(name);
                if k < 1 {
                    err!("**** Error for sin parameters: Specified location name %s is not recognized.", name.as_str());
                }
                k
            }
        };
        kase.push(k);
    }
    let mut stat = String::new();
    if getparstring("sstat", &mut stat) == 0 {
        stat = "sstat".to_string();
    }
    let mut mult: f64 = 1.0;
    if getpardouble("smult", &mut mult) == 0 {
        mult = 1.0;
    }
    if mult == 0.0 {
        err!("**** Error for sin parameters: smult cannot be 0.0");
    }
    let mut limit = f64::MAX;
    if getpardouble("slimit", &mut limit) != 0 {
        if limit < 0.0 {
            err!("**** Error for sin parameters: slimit cannot be less than 0.0");
        }
        limit *= limit;
    }
    let fp = match fopen(&filename, "r") {
        Some(fp) => fp,
        None => {
            err!("error sin file did not open correctly.");
        }
    };
    // Only store the requested locations and static, in Q-file order.
    let mut pname = key.clone();
    pname.push(stat.clone());
    let info = match getviaqfile(&fp, &pname, pname.len() as i32, 0) {
        Ok(info) => info,
        Err(errwarn) => match errwarn {
            1 => err!("sin file getqinfo error: extra C_SU_NAMES record in q-file"),
            2 => err!("sin file getqinfo error: extra C_SU_NDIMS record in q-file"),
            3 => err!("sin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            11 => err!("sin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"),
            12 => err!("sin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            22 => err!("sin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),
            23 => err!("sin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
            24 => err!("sin file getviaqfile error: C_SU_NDIMS blank where valid number expected"),
            25 => err!("sin file getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
            26 => err!("sin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
            27 => err!("sin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),
            n if n > 100 => err!("sin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)", n - 99),
            n if n > 0 => err!("sin file getviaqfile error: unrecognized error code %d",n),
            _ => unreachable!(),
        }
    };
    if info.ifixd == 0 {
        err!("sin file error: input with varying number of tuples is not allowed.");
    }
    let mut locn = Vec::new();
    for name in &key {
        let mut found = None;
        for j in 0..info.iztuple as usize {
            if info.names[j] == *name {
                found = Some(j);
            }
        }
        match found {
            Some(j) => locn.push(j),
            None => err!(
                "sin file error: skeyloc %s not found in non-tuple part of sin Q-file.",
                name.as_str()
            ),
        }
    }
    let mut oloc = None;
    for j in 0..info.iztuple as usize {
        if info.names[j] == stat {
            oloc = Some(j);
        }
    }
    let oloc = match oloc {
        Some(j) => j,
        None => err!(
            "sin file error: Q-file must have sstat parameter value %s (among non-tuple names).",
            stat.as_str()
        ),
    };
    let npoints = info.records.len();
    if npoints < 1 {
        err!("sin file error: no records starting with Q found in file.");
    }
    let mut tree_dl = vec![vec![0.0; npoints]; key.len()];
    for i in 0..key.len() {
        for j in 0..npoints {
            tree_dl[i][j] = info.records[j].dlots[locn[i]];
        }
    }
    let nodes = connect_nodes(npoints, &tree_dl, key.len(), 1);
    Some(Statics {
        kase,
        info,
        oloc,
        mult,
        limit,
        tree_dl,
        nodes,
        near_dist: 0.0,
    })
}

fn setup_r() -> Option<Statics> {
    let mut filename = String::new();
    if countparval("rin") > 0 {
        getparstring("rin", &mut filename);
    } else {
        filename = "rstat.csv".to_string();
    }
    if filename == "none" {
        return None;
    }
    let mut key = Vec::new();
    let numd = countparval("rkeyloc");
    if numd < 1 {
        key.push("gaps".to_string());
    } else {
        getparstringarray("rkeyloc", &mut key);
        if numd > 9 {
            err!("**** Error for rin parameters: Maximum of 9 location names can be specified.");
        }
    }
    let mut kase = Vec::new();
    for name in &key {
        let k = match name.as_str() {
            "sx" => 101,
            "sy" => 102,
            "gx" => 103,
            "gy" => 104,
            "selev" => 111,
            "gelev" => 112,
            "sdepth" => 113,
            "sdel" => 114,
            "gdel" => 115,
            "swdep" => 116,
            "gwdep" => 117,
            _ => {
                let k = GetCase(name);
                if k < 1 {
                    err!("**** Error for rin parameters: Specified location name %s is not recognized.", name.as_str());
                }
                k
            }
        };
        kase.push(k);
    }
    let mut stat = String::new();
    if getparstring("rstat", &mut stat) == 0 {
        stat = "gstat".to_string();
    }
    let mut mult: f64 = 1.0;
    if getpardouble("rmult", &mut mult) == 0 {
        mult = 1.0;
    }
    if mult == 0.0 {
        err!("**** Error for rin parameters: rmult cannot be 0.0");
    }
    let mut limit = f64::MAX;
    if getpardouble("rlimit", &mut limit) != 0 {
        if limit < 0.0 {
            err!("**** Error for rin parameters: rlimit cannot be less than 0.0");
        }
        limit *= limit;
    }
    let fp = match fopen(&filename, "r") {
        Some(fp) => fp,
        None => {
            err!("error: rin file did not open correctly.");
        }
    };
    // Only store the requested locations and static, in Q-file order.
    let mut pname = key.clone();
    pname.push(stat.clone());
    let info = match getviaqfile(&fp, &pname, pname.len() as i32, 0) {
        Ok(info) => info,
        Err(errwarn) => match errwarn {
            1 => err!("rin file getqinfo error: extra C_SU_NAMES record in q-file"),
            2 => err!("rin file getqinfo error: extra C_SU_NDIMS record in q-file"),
            3 => err!("rin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            11 => err!("rin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"),
            12 => err!("rin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            22 => err!("rin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),
            23 => err!("rin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
            24 => err!("rin file getviaqfile error: C_SU_NDIMS blank where valid number expected"),
            25 => err!("rin file getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
            26 => err!("rin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
            27 => err!("rin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),
            n if n > 100 => err!("rin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)", n - 99),
            n if n > 0 => err!("rin file getviaqfile error: unrecognized error code %d",n),
            _ => unreachable!(),
        }
    };
    if info.ifixd == 0 {
        err!("rin file error: input with varying number of tuples is not allowed.");
    }
    let mut locn = Vec::new();
    for name in &key {
        let mut found = None;
        for j in 0..info.iztuple as usize {
            if info.names[j] == *name {
                found = Some(j);
            }
        }
        match found {
            Some(j) => locn.push(j),
            None => err!(
                "rin file error: rkeyloc %s not found in non-tuple part of rin Q-file.",
                name.as_str()
            ),
        }
    }
    let mut oloc = None;
    for j in 0..info.iztuple as usize {
        if info.names[j] == stat {
            oloc = Some(j);
        }
    }
    let oloc = match oloc {
        Some(j) => j,
        None => err!(
            "error: rin Q-file must have rstat parameter value %s (among non-tuple names).",
            stat.as_str()
        ),
    };
    let npoints = info.records.len();
    if npoints < 1 {
        err!("rin file error: no records starting with Q found in file.");
    }
    let mut tree_dl = vec![vec![0.0; npoints]; key.len()];
    for i in 0..key.len() {
        for j in 0..npoints {
            tree_dl[i][j] = info.records[j].dlots[locn[i]];
        }
    }
    let nodes = connect_nodes(npoints, &tree_dl, key.len(), 1);
    Some(Statics {
        kase,
        info,
        oloc,
        mult,
        limit,
        tree_dl,
        nodes,
        near_dist: 0.0,
    })
}

fn setup_c() -> Option<Statics> {
    let mut filename = String::new();
    if countparval("cin") > 0 {
        getparstring("cin", &mut filename);
    } else {
        filename = "none".to_string();
    }
    if filename == "none" {
        return None;
    }
    let mut key = Vec::new();
    let numd = countparval("ckeyloc");
    if numd < 1 {
        key.push("cdp".to_string());
    } else {
        getparstringarray("ckeyloc", &mut key);
        if numd > 9 {
            err!("**** Error for cin parameters: Maximum of 9 location names can be specified.");
        }
    }
    let mut kase = Vec::new();
    for name in &key {
        let k = match name.as_str() {
            "sx" => 101,
            "sy" => 102,
            "gx" => 103,
            "gy" => 104,
            "selev" => 111,
            "gelev" => 112,
            "sdepth" => 113,
            "sdel" => 114,
            "gdel" => 115,
            "swdep" => 116,
            "gwdep" => 117,
            _ => {
                let k = GetCase(name);
                if k < 1 {
                    err!("**** Error for cin parameters: Specified location name %s is not recognized.", name.as_str());
                }
                k
            }
        };
        kase.push(k);
    }
    let mut stat = String::new();
    if getparstring("cstat", &mut stat) == 0 {
        stat = "tstat".to_string();
    }
    let mut mult: f64 = -2.0;
    if getpardouble("cmult", &mut mult) == 0 {
        mult = -2.0;
    }
    if mult == 0.0 {
        err!("**** Error for cin parameters: cmult cannot be 0.0");
    }
    let mut limit = f64::MAX;
    if getpardouble("climit", &mut limit) != 0 {
        if limit < 0.0 {
            err!("**** Error for cin parameters: climit cannot be less than 0.0");
        }
        limit *= limit;
    }
    let fp = match fopen(&filename, "r") {
        Some(fp) => fp,
        None => {
            err!("error: cin file did not open correctly.");
        }
    };
    // Only store the requested locations and static, in Q-file order.
    let mut pname = key.clone();
    pname.push(stat.clone());
    let info = match getviaqfile(&fp, &pname, pname.len() as i32, 0) {
        Ok(info) => info,
        Err(errwarn) => match errwarn {
            1 => err!("cin file getqinfo error: extra C_SU_NAMES record in q-file"),
            2 => err!("cin file getqinfo error: extra C_SU_NDIMS record in q-file"),
            3 => err!("cin file getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            11 => err!("cin file readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"),
            12 => err!("cin file readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            22 => err!("cin file getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),
            23 => err!("cin file getviaqfile error: C_SU_NAMES tupled names out-of-order, changed"),
            24 => err!("cin file getviaqfile error: C_SU_NDIMS blank where valid number expected"),
            25 => err!("cin file getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
            26 => err!("cin file getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
            27 => err!("cin file getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),
            n if n > 100 => err!("cin file getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)", n - 99),
            n if n > 0 => err!("cin file getviaqfile error: unrecognized error code %d",n),
            _ => unreachable!(),
        }
    };
    if info.ifixd == 0 {
        err!("cin file error: input with varying number of tuples is not allowed.");
    }
    let mut locn = Vec::new();
    for name in &key {
        let mut found = None;
        for j in 0..info.iztuple as usize {
            if info.names[j] == *name {
                found = Some(j);
            }
        }
        match found {
            Some(j) => locn.push(j),
            None => err!(
                "cin file error: ckeyloc %s not found in non-tuple part of cin Q-file.",
                name.as_str()
            ),
        }
    }
    let mut oloc = None;
    for j in 0..info.iztuple as usize {
        if info.names[j] == stat {
            oloc = Some(j);
        }
    }
    let oloc = match oloc {
        Some(j) => j,
        None => err!(
            "error: cin Q-file must have cstat parameter value %s (among non-tuple names).",
            stat.as_str()
        ),
    };
    let npoints = info.records.len();
    if npoints < 1 {
        err!("cin file error: no records starting with Q found in file.");
    }
    let mut tree_dl = vec![vec![0.0; npoints]; key.len()];
    for i in 0..key.len() {
        for j in 0..npoints {
            tree_dl[i][j] = info.records[j].dlots[locn[i]];
        }
    }
    let nodes = connect_nodes(npoints, &tree_dl, key.len(), 1);
    Some(Statics {
        kase,
        info,
        oloc,
        mult,
        limit,
        tree_dl,
        nodes,
        near_dist: 0.0,
    })
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    // Infinite extents; all dimensions contribute to Pythagorean distance.
    let extent_min = [-f64::MAX; 9];
    let extent_max = [f64::MAX; 9];
    let tree_nt = [1; 9];
    let mut kdist = 5.0;
    let mut kmult = 2.0;
    let mut kdadd = 1;
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

    let mut s = setup_s();
    let mut r = setup_r();
    let mut c = setup_c();
    checkpars();

    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Error: cannot get first trace");
    }
    let lentrc = tr.ns as usize;
    let mut onetrc = vec![0.0f32; lentrc];
    let mut oneshift = vec![0.0f32; lentrc];
    let sampi = (tr.dt as f32 as f64 / 1000.0) as f32;
    let mut nproct: i32 = 0;
    loop {
        let mut dt = 0.0f64;
        if let Some(s) = &mut s {
            if let Some(value) = static_for_trace(
                s,
                &tr,
                kdadd,
                kdist,
                kmult,
                &tree_nt,
                &extent_min,
                &extent_max,
            ) {
                dt += value;
            }
        }
        if let Some(r) = &mut r {
            if let Some(value) = static_for_trace(
                r,
                &tr,
                kdadd,
                kdist,
                kmult,
                &tree_nt,
                &extent_min,
                &extent_max,
            ) {
                dt += value;
            }
        }
        if let Some(c) = &mut c {
            if let Some(value) = static_for_trace(
                c,
                &tr,
                kdadd,
                kdist,
                kmult,
                &tree_nt,
                &extent_min,
                &extent_max,
            ) {
                dt += value;
            }
        }
        for i in 0..lentrc {
            onetrc[i] = tr.data[i];
            oneshift[i] = (i as f32 as f64 - dt / sampi as f64) as f32;
        }
        ints8r(
            1.0,
            0.0,
            &onetrc,
            0.0,
            0.0,
            &oneshift,
            &mut tr.data[..lentrc],
        );
        puttr(&mut tr);
        nproct = nproct.wrapping_add(1);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    warn!("Number of traces=%d", nproct);
    exit(CWP_Exit());
}

#[allow(clippy::too_many_arguments)]
fn static_for_trace(
    s: &mut Statics,
    tr: &Trace,
    kdadd: i32,
    kdist: f64,
    kmult: f64,
    tree_nt: &[i32],
    extent_min: &[f64],
    extent_max: &[f64],
) -> Option<f64> {
    let numd = s.kase.len();
    let mut target = [0.0; 9];
    for i in 0..numd {
        target[i] = scaledfromhead(tr, s.kase[i]);
    }
    let mut dist2 = kdist;
    if kdadd == 1 {
        dist2 = kdist + s.near_dist.sqrt();
    }
    let near = cycle_for_near(
        &s.nodes, &s.tree_dl, tree_nt, numd, extent_min, extent_max, &target, dist2, kmult,
    );
    s.near_dist = near.dist;
    if s.near_dist <= s.limit {
        Some(s.info.records[near.elem].dlots[s.oloc] * s.mult)
    } else {
        None
    }
}

#[derive(Clone, Default)]
struct Node {
    elem: usize,
    l: Option<usize>,
    r: Option<usize>,
}

// These kd-tree routines are local functions in sushiftcsv.c.
fn connect_nodes(tree_numc: usize, tree_dl: &[Vec<f64>], tree_numd: usize, ihop: i32) -> Vec<Node> {
    let mut tree_nodes = vec![Node::default(); tree_numc];
    for m in 0..tree_numc {
        tree_nodes[m].elem = m;
    }
    if ihop == 1 {
        let mut m = 0;
        let mut nrat = tree_numc;
        while nrat > 0 {
            let nd = if nrat > 6 {
                (nrat as f64).sqrt() as usize
            } else {
                0
            };
            let mut n = nd + nrat / 2;
            while n < tree_numc {
                if tree_nodes[n].l.is_none() {
                    tree_nodes[n].l = Some(0); // Temporarily flags a loaded point.
                    tree_nodes[m].elem = n;
                    m += 1;
                }
                n += nd + nrat;
            }
            nrat = (nrat as f64 * 0.6) as usize;
        }
        for m in 0..tree_numc {
            tree_nodes[m].l = None;
        }
    }
    connect_all(&mut tree_nodes, tree_dl, tree_numd);
    tree_nodes
}

fn connect_all(tree_nodes: &mut [Node], tree_dl: &[Vec<f64>], tree_numd: usize) {
    for m in 1..tree_nodes.len() {
        let mut now_node = 0;
        let mut next_node = Some(0);
        let mut naxe = 0;
        while let Some(next) = next_node {
            now_node = next;
            naxe += 1;
            if naxe == tree_numd {
                naxe = 0;
            }
            if tree_dl[naxe][tree_nodes[m].elem] < tree_dl[naxe][tree_nodes[now_node].elem] {
                next_node = tree_nodes[now_node].l;
            } else {
                next_node = tree_nodes[now_node].r;
            }
        }
        if tree_dl[naxe][tree_nodes[m].elem] < tree_dl[naxe][tree_nodes[now_node].elem] {
            tree_nodes[now_node].l = Some(m);
        } else {
            tree_nodes[now_node].r = Some(m);
        }
    }
}

struct Near {
    elem: usize,
    dist: f64,
    num_found: u64,
    ncycles: i32,
}

#[allow(clippy::too_many_arguments)]
fn find_near_rcur(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_nt: &[i32],
    tree_numd: usize,
    now_node: Option<usize>,
    mut naxe: usize,
    extent_min: &[f64],
    extent_max: &[f64],
    target: &[f64],
    near: &mut Near,
) {
    let Some(now) = now_node else {
        return;
    };
    let node = &tree_nodes[now];
    let mut inside = true;
    for i in 0..tree_numd {
        if tree_dl[i][node.elem] < extent_min[i] || tree_dl[i][node.elem] >= extent_max[i] {
            inside = false;
            break;
        }
    }
    if inside {
        let mut rad = 0.0;
        for i in 0..tree_numd {
            if tree_nt[i] != 0 {
                rad += (target[i] - tree_dl[i][node.elem]) * (target[i] - tree_dl[i][node.elem]);
            }
        }
        if rad <= near.dist {
            if rad < near.dist {
                near.num_found = 1;
                near.elem = node.elem;
                near.dist = rad;
            } else {
                near.num_found += 1;
                if node.elem > near.elem {
                    near.elem = node.elem;
                }
            }
        }
    }
    if naxe == tree_numd {
        naxe = 0;
    }
    if tree_dl[naxe][node.elem] >= extent_min[naxe] && node.l.is_some() {
        find_near_rcur(
            tree_nodes,
            tree_dl,
            tree_nt,
            tree_numd,
            node.l,
            naxe + 1,
            extent_min,
            extent_max,
            target,
            near,
        );
    }
    if tree_dl[naxe][node.elem] < extent_max[naxe] && node.r.is_some() {
        find_near_rcur(
            tree_nodes,
            tree_dl,
            tree_nt,
            tree_numd,
            node.r,
            naxe + 1,
            extent_min,
            extent_max,
            target,
            near,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn cycle_for_near(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_nt: &[i32],
    tree_numd: usize,
    extent_min: &[f64],
    extent_max: &[f64],
    target: &[f64],
    init_rad: f64,
    rad_scal: f64,
) -> Near {
    let mut now_rad = init_rad;
    let mut loc_min = [0.0; 9];
    let mut loc_max = [0.0; 9];
    let mut nset = 0;
    for i in 0..tree_numd {
        if tree_nt[i] == 0 {
            loc_min[i] = extent_min[i];
            loc_max[i] = extent_max[i];
            nset += 2;
        }
    }
    let mut near = Near {
        elem: 0,
        dist: f64::MAX,
        num_found: 0,
        ncycles: 0,
    };
    loop {
        near.ncycles = near.ncycles.wrapping_add(1);
        let mut nbig = nset;
        for i in 0..tree_numd {
            if tree_nt[i] != 0 {
                loc_min[i] = target[i] - now_rad;
                loc_max[i] = target[i] + now_rad;
                if loc_min[i] <= extent_min[i] {
                    loc_min[i] = extent_min[i];
                    nbig += 1;
                }
                if loc_max[i] >= extent_max[i] {
                    loc_max[i] = extent_max[i];
                    nbig += 1;
                }
            }
        }
        near.num_found = 0;
        near.dist = f64::MAX;
        near.elem = 0;
        find_near_rcur(
            tree_nodes,
            tree_dl,
            tree_nt,
            tree_numd,
            Some(0),
            1,
            &loc_min,
            &loc_max,
            target,
            &mut near,
        );
        if near.num_found < 1 {
            if nbig == 2 * tree_numd {
                return near;
            }
            now_rad *= rad_scal;
            continue;
        }
        if near.dist > now_rad * now_rad && nbig != 2 * tree_numd {
            now_rad = near.dist.sqrt() * 1.001;
            continue;
        }
        return near;
    }
}

fn scaledfromhead(tr: &Trace, k: i32) -> f64 {
    let (mut dval, scale) = match k {
        101 => (tr.sx as f64, tr.scalco),
        102 => (tr.sy as f64, tr.scalco),
        103 => (tr.gx as f64, tr.scalco),
        104 => (tr.gy as f64, tr.scalco),
        111 => (tr.selev as f64, tr.scalel),
        112 => (tr.gelev as f64, tr.scalel),
        113 => (tr.sdepth as f64, tr.scalel),
        114 => (tr.sdel as f64, tr.scalel),
        115 => (tr.gdel as f64, tr.scalel),
        116 => (tr.swdep as f64, tr.scalel),
        117 => (tr.gwdep as f64, tr.scalel),
        _ => return fromhead(tr, k),
    };
    if scale > 1 {
        dval *= scale as f64;
    } else if scale < 0 {
        dval /= -(scale as i32) as f64;
    }
    dval
}
