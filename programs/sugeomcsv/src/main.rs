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

//! SUGEOMCSV - land/OBotm geometry from CSV or fixed text.
//! Port of su/main/headers/sugeomcsv.c (Revision 1.2).
//! Credits: Andre Latour. Started from sugeom by Fernando Roxo <petro@roxo.org>.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;

#[derive(Clone)]
struct PointInfo {
    dfield: Vec<f64>,
    lfield: Vec<i64>,
}

fn longt(dvalue: f64, dtolh: f64, dtol: f64) -> i64 {
    if dvalue >= 0.0 {
        ((dvalue + dtolh) * dtol) as i64
    } else {
        ((dvalue - dtolh) * dtol) as i64
    }
}

fn comp_sort(p1: &PointInfo, p2: &PointInfo, num: usize) -> i32 {
    for n in 0..num {
        if p1.lfield[n] < p2.lfield[n] {
            return -1;
        }
        if p1.lfield[n] > p2.lfield[n] {
            return 1;
        }
    }
    0
}

fn bhigh(all: &[PointInfo], guy: &PointInfo, num: usize) -> usize {
    let mut low = 0;
    let mut high = all.len();
    while low < high {
        let mid = low + (high - low) / 2;
        if comp_sort(guy, &all[mid], num) >= 0 {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low
}

fn tparse(text: &str, d: char) -> Vec<String> {
    text.split(d)
        .map(|s| {
            if s.is_empty() {
                "null".into()
            } else {
                s.into()
            }
        })
        .collect()
}

// The same record selection is used by countRec and the read-in pass.
fn data_records<'a>(
    fp: &'a CFile,
    rid: &'a str,
    lenid: usize,
    nicerecord: i32,
) -> impl Iterator<Item = (i32, String)> + 'a {
    let mut ncount = 0;
    let mut nextrow = false;
    std::iter::from_fn(move || {
        while let Some(raw) = fgets(10001, fp) {
            ncount += 1;
            if ncount < nicerecord {
                continue;
            }
            let front = raw.to_ascii_lowercase();
            if front.starts_with("c_su") || nextrow {
                nextrow = front.starts_with("c_su_names") || front.starts_with("c_su_forms");
            } else if lenid == 0 || raw.as_bytes().starts_with(&rid.as_bytes()[..lenid]) {
                return Some((ncount, raw));
            }
        }
        None
    })
}

#[derive(Default)]
struct ReadErrors {
    lenerr: i32,
    comerr: i32,
    morerr: i32,
    numerr: i32,
    lrgerr: i32,
    nblank: i32,
}

fn get_csv(
    raw: &str,
    rdel: u8,
    dfield: &mut [f64],
    nspot: &[usize],
    ncount: i32,
    e: &mut ReadErrors,
) {
    let mut nbeg: i32 = -1;
    let mut nfield = 0;
    let mut ineed = 0;
    // C scans the entire zero-filled work buffer even after a short record.
    let mut text = vec![0u8; 10001];
    let len = raw.len().min(10000);
    text[..len].copy_from_slice(&raw.as_bytes()[..len]);
    for n in 0..10001 {
        let ch = text[n];
        if ch == rdel || ch == 0 || ch == b'\n' || ch == b'\r' {
            if ineed < dfield.len() && nfield == nspot[ineed] {
                let mut dval = 1.1e308;
                let mut nb: i32 = -1;
                if n as i32 - nbeg - 1 > 0 {
                    let field = &text[(nbeg + 1) as usize..n];
                    let fieldstr = String::from_utf8_lossy(field);
                    let mut ib: i32 = -1;
                    for (m, &c) in field.iter().enumerate() {
                        if c != b' ' {
                            nb = m as i32;
                            if ib > -1 {
                                e.morerr += 1;
                                if e.morerr < 4 {
                                    warn!(
                                        "Error at record %d   two-numbers in field (%s)",
                                        ncount,
                                        fieldstr.as_ref()
                                    );
                                    if e.morerr == 3 {
                                        warn!("Have 3 two-numbers in field warnings, no more will be printed.");
                                    }
                                }
                                nb = -1;
                                break;
                            }
                        }
                        if c == b' ' && nb > -1 {
                            ib = m as i32;
                        }
                    }
                    if nb > -1 {
                        // strtod has the same number grammar as sscanf("%lf").
                        let (value, consumed) = strtod(fieldstr.as_ref());
                        if consumed == 0 {
                            e.numerr += 1;
                            if e.numerr < 4 {
                                warn!(
                                    "Error at record %d   field-unreadable as a number (%s)",
                                    ncount,
                                    fieldstr.as_ref()
                                );
                                if e.numerr == 3 {
                                    warn!("Have 3 field-unreadable warnings, no more will be printed.");
                                }
                            }
                            break;
                        }
                        dval = value;
                    }
                }
                if nb < 0 {
                    e.nblank += 1;
                    dval = 0.0;
                }
                dfield[ineed] = dval;
                ineed += 1;
                if ineed >= dfield.len() {
                    break;
                }
            }
            if ch == 0 || ch == b'\n' || ch == b'\r' {
                if ineed < dfield.len() {
                    e.comerr += 1;
                    if e.comerr < 4 {
                        warn!(
                            "Error at record %d   Not-enough-commas in record to get all values",
                            ncount
                        );
                        if e.comerr == 3 {
                            warn!("Have 3 Not-enough-comma warnings, no more will be printed.");
                        }
                    }
                }
            }
            nbeg = n as i32;
            nfield += 1;
        }
    }
}

// State used in both unrepeat passes. The text pass starts at 1; the
// trace pass starts at the user's unrepeat value.
struct Unrepeat {
    nup: i32,
    nlfield: i64,
    nrep: i64,
}
const INCINT: i64 = 100000000000000;
impl Unrepeat {
    fn new(start: i64) -> Self {
        Self {
            nup: 0,
            nlfield: -999999999999999,
            nrep: start.wrapping_mul(INCINT),
        }
    }
    fn apply(&mut self, value: &mut i64) {
        if self.nup > 0 {
            if self.nlfield > *value {
                self.nup = -1;
                self.nrep = self.nrep.wrapping_add(INCINT);
            }
        } else if self.nup < 0 {
            if self.nlfield < *value {
                self.nup = 1;
                self.nrep = self.nrep.wrapping_add(INCINT);
            }
        } else if self.nlfield > -999999999999999 {
            self.nup = 1;
            if self.nlfield > *value {
                self.nup = -1;
            }
        }
        self.nlfield = *value;
        *value = value.wrapping_add(self.nrep);
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    let mut nicerecord = 1;
    getparint("nicerecord", &mut nicerecord);
    if nicerecord < 1 {
        err!("**** Error: nicerecord= cannot be less than 1");
    }
    let mut num_r: i32 = 0;
    getparint("maxrecords", &mut num_r);
    if num_r < 0 {
        err!("**** Error: maxrecords= cannot be less than 0");
    }
    let mut ftol: f32 = 0.01;
    getparfloat("tolr", &mut ftol);
    if (ftol as f64) < 0.000000001 {
        err!("**** Error: tolr= must be larger.");
    }
    let dtolh = ftol as f64 / 2.0;
    let dtol = 1.0 / ftol as f64;
    let mut unrepeat: i32 = -2147483645;
    getparint("unrepeat", &mut unrepeat);

    let mut iscalel: i32 = 10;
    let iscaleldef = getparint("scalel", &mut iscalel) == 0;
    if iscalel == -1 {
        iscalel = 1;
    }
    if ![1, 10, 100, 1000, 10000, 100000, 1000000, 10000000].contains(&iscalel.wrapping_abs()) {
        err!("**** Error: scalel= must be signed powers of 10 (1,10,100...-10,-100,...)");
    }
    let dscalel = if iscalel > 0 {
        iscalel as f64
    } else {
        -1.0 / iscalel as f64
    };
    iscalel *= -1;
    let mut iscalco: i32 = 10;
    let iscalcodef = getparint("scalco", &mut iscalco) == 0;
    if iscalco == -1 {
        iscalco = 1;
    }
    if ![1, 10, 100, 1000, 10000, 100000, 1000000, 10000000].contains(&iscalco.wrapping_abs()) {
        err!("**** Error: scalco= must be signed powers of 10 (1,10,100...-10,-100,...)");
    }
    let dscalco = if iscalco > 0 {
        iscalco as f64
    } else {
        -1.0 / iscalco as f64
    };
    iscalco *= -1;
    let mut tmiss = String::new();
    let mut missing = -1;
    if getparstring("missing", &mut tmiss) != 0 {
        match tmiss.to_ascii_lowercase().as_str() {
            "delete" => missing = 0,
            "pass" => missing = 1,
            _ => err!("**** Error: missing= option not recognized."),
        }
    }
    let mut tact = String::new();
    let mut iaction = 0;
    if getparstring("action", &mut tact) != 0 {
        match tact.to_ascii_lowercase().as_str() {
            "add" => iaction = 1,
            "subtract" => iaction = -1,
            _ => err!("**** Error: action= option not recognized."),
        }
    }
    let mut tall = String::new();
    let mut ncreate: i32 = 0;
    if getparstring("create", &mut tall) != 0 {
        if tall.eq_ignore_ascii_case("all") {
            ncreate = 2147483645;
        }
        if ncreate == 0 {
            getparint("create", &mut ncreate);
        }
        if ncreate < 0 {
            err!("**** Error create= cannot be negative.");
        }
    }
    let mut spikes = Vec::<f32>::new();
    let mut ireplicate = 0;
    let nspikes = countparval("spikes");
    if nspikes > 0 {
        if ncreate < 1 {
            err!("**** Error: spikes= can only be specified on command line if create= is also there.");
        }
        if nspikes > 1 {
            if nspikes % 2 == 1 {
                err!("**** Error: spikes= must have an even number of values (time amplitude pairs).");
            }
            getparfloat("spikes", &mut spikes);
        } else {
            getparint("spikes", &mut ireplicate);
            if ireplicate < 1 {
                err!("**** Error: For spikes=n the n must be +integer (the trace to input).");
            }
        }
    }
    if ncreate > 0 {
        if !isatty(0) {
            if ireplicate < 1 {
                err!("**** Error: Input trace file not allowed when create= (unless spikes=n).");
            }
        } else if ireplicate > 0 {
            err!("**** Error: Input trace file must be specified when spikes=n.");
        }
        if unrepeat > -2147483645 && unrepeat != 1 {
            err!("**** Error: If specified, only unrepeat=1 allowed when creating traces.");
        }
    } else if isatty(0) {
        err!("**** Error: Input trace file must be specified when not create= ");
    }

    let mut rname = String::new();
    if getparstring("rfile", &mut rname) == 0 {
        err!("**** Error: rfile= text file name must be specified.");
    }
    let mut rdel = b',';
    let mut tdel = String::new();
    let has_del = getparstring("rdelim", &mut tdel) != 0;
    if has_del {
        let b = tdel.as_bytes();
        if b.len() == 1 {
            rdel = b[0];
        } else if b.len() == 3 && b[0] == b'"' && b[2] == b'"' {
            rdel = b[1];
        } else {
            err!("**** Error: rdelim= specification not recognized.");
        }
    }
    let mut irtype = true;
    if countparval("rtype") > 0 {
        let mut rtype = String::new();
        getparstring("rtype", &mut rtype);
        match rtype.to_ascii_lowercase().as_str() {
            "csv" => irtype = true,
            "fixed" => irtype = false,
            _ => err!("**** Error: rtype= option not recognized."),
        }
    } else if rname.len() > 2 {
        irtype = rname.ends_with("csv");
    }
    if !irtype && has_del {
        err!("**** Error: you cannot specify rdelim= if rtype=fixed.");
    }
    let fpr = fopen(&rname, "r").unwrap_or_else(|| err!("**** Error opening the rfile text file."));
    let mut uname = String::new();
    let has_u = getparstring("ufile", &mut uname) != 0;
    let iread_c_su = !(has_u && uname == "command");
    let fpu = if !has_u || uname == rname || !iread_c_su {
        fpr.clone()
    } else {
        fopen(&uname, "r").unwrap_or_else(|| err!("**** Error opening the ufile"))
    };
    let mut matches = Vec::<String>::new();
    if countparval("match") > 0 {
        getparstringarray("match", &mut matches);
    }
    let mut nerr = false;
    if matches.is_empty() && !iread_c_su {
        nerr = true;
        warn!("**** Must have match= on command line when ufile=command.");
    }
    let mut rid = String::new();
    if countparval("setid") > 0 {
        getparstring("setid", &mut rid);
    }
    let mut lenid = rid.len();
    let mut names = Vec::<String>::new();
    if countparval("names") > 0 {
        getparstringarray("names", &mut names);
        for s in &mut names {
            s.make_ascii_lowercase();
        }
    } else if !iread_c_su {
        nerr = true;
        warn!("**** Must have names= on command line when ufile=command.");
    }
    if nerr {
        if !iread_c_su {
            warn!("**** If trying to read C_SU records in text file, cannot have ufile=command.");
        }
        err!("**** A required parameter has not been specified. See previous messages.");
    }
    let mut names_more = false;
    if names.is_empty() || matches.is_empty() || lenid == 0 {
        let in_names = names.len();
        let in_match = matches.len();
        let in_lenid = lenid;
        let mut num_names = 0;
        let mut num_match = 0;
        let mut num_setid = 0;
        let mut read_names = 0;
        while let Some(raw) = fgets(10001, &fpu) {
            if (in_names > 0 || read_names == -2) && !matches.is_empty() && lenid > 0 {
                break;
            }
            let text: String = raw
                .chars()
                .take_while(|&c| c != '\0' && c != '\n' && c != '\r')
                .filter(|&c| c != ' ' && c != '\t')
                .collect();
            let lower = text.to_ascii_lowercase();
            if lenid == 0 && lower.starts_with("c_su_setid") {
                let ids = tparse(&text, ',');
                if ids.len() > 1 {
                    rid = ids[1].clone();
                    lenid = rid.len();
                }
            }
            if lower.starts_with("c_su_match") {
                num_match += 1;
            }
            if lower.starts_with("c_su_setid") {
                num_setid += 1;
            }
            if lower.starts_with("c_su_names") {
                num_names += 1;
            }
            if matches.is_empty() && lower.starts_with("c_su_matc") {
                let fields = tparse(&lower, ',');
                for s in fields.into_iter().skip(1) {
                    if s == "null" {
                        break;
                    }
                    matches.push(s);
                }
            }
            if read_names == -1 {
                if lower.starts_with("c_su_more") {
                    read_names = 1;
                    names_more = true;
                } else {
                    read_names = -2;
                }
            }
            if read_names > 0 {
                names.extend(tparse(&lower, ','));
                read_names = -1;
            }
            if in_names == 0 && lower.starts_with("c_su_names") {
                read_names = 1;
            }
        }
        fseeko(&fpu, 0, SEEK_SET);
        let mut ierr = false;
        if in_names < 1 {
            if num_names > 1 {
                ierr = true;
                warn!("**** Error: Text file has more than one C_SU_NAMES parameter record.");
                warn!("****        Remove duplicates or override with names= on command line.");
            } else if num_names == 0 {
                ierr = true;
                warn!("**** Error: Text file has no C_SU_NAMES parameter record.");
                warn!("****        Add it, or specify names= on command line.");
            }
        }
        if in_match < 1 {
            if num_match > 1 {
                ierr = true;
                warn!("**** Error: Text file has more than one C_SU_MATCH record.");
                warn!("****        Remove duplicates or override with match= on command line.");
            } else if num_match == 0 {
                ierr = true;
                warn!("**** Error: Text file has no C_SU_MATCH record.");
                warn!("****        Add it, or specify match= command line.");
            }
        }
        if in_lenid == 0 {
            if num_setid > 1 {
                ierr = true;
                warn!("**** Error: Text file has more than one C_SU_SETID record.");
                warn!("****        Remove duplicates or override with setid= on command line.");
            } else if num_setid == 0 {
                ierr = true;
                warn!("**** Error: Text file has no C_SU_SETID record.");
                warn!("****        Add it, or specify setid= on command line.");
            }
        }
        if ierr {
            err!("**** Error: Text file has duplicate or missing C_SU_NAMES or C_SU_MATCH or C_SU_SETID records.");
        }
    }
    if lenid > 3 && rid.starts_with('"') && rid.ends_with('"') {
        rid = rid[1..lenid - 1].into();
        lenid -= 2;
    } else {
        rid.make_ascii_uppercase();
    }
    if (lenid == 3 && rid.eq_ignore_ascii_case("any"))
        || (lenid == 4 && rid.eq_ignore_ascii_case("none"))
    {
        lenid = 0;
    }
    if matches[0] == "sps2" || matches[0] == "sps1" {
        let p = match rid.as_str() {
            "X" => Some(["fldr", "tracf"]),
            "S" => Some(["grnofr", "grnlof"]),
            "R" => Some(["grnors", "gaps"]),
            _ => None,
        };
        if let Some(p) = p {
            matches = p.iter().map(|s| s.to_string()).collect();
        }
    }
    if let Some(p) = sps_names(&names[0], &rid) {
        names = p;
    }
    names.retain(|s| !s.starts_with("c_su_more"));
    if names[0].starts_with("c_su_id") {
        names[0] = "null".into();
    }
    let all_names = names.len();
    let mut ilead = vec![-1i32; all_names];
    let mut itrail = vec![-1i32; all_names];
    let mut namex = vec![String::new(); all_names];
    let mut incomma = false;
    let mut extra_parts = 0;
    for n in 0..all_names {
        let parts = tparse(&names[n], '_');
        let num_parts = parts.len();
        match num_parts {
            1 => names[n] = parts[0].clone(),
            2 => {
                names[n] = parts[0].clone();
                namex[n] = parts[1].clone();
            }
            3 | 4 => {
                ilead[n] = atoi(&parts[0]);
                itrail[n] = atoi(&parts[num_parts - 1]);
                if itrail[n] < ilead[n] || ilead[n] < 1 {
                    err!(
                        "**** Error: Your names= list has an entry with incorrect integers: %s",
                        names[n]
                    );
                }
                names[n] = parts[1].clone();
                if num_parts == 4 {
                    namex[n] = parts[2].clone();
                }
            }
            _ => err!(
                "**** Error: Your names= list has an entry that parses incorrectly: %s",
                names[n]
            ),
        }
        if names[n] != "null" && names[n] != "c_su_id" && names[n] != "numb" {
            if num_parts < 3 {
                incomma = true;
            }
            if num_parts == 2 || num_parts == 4 {
                extra_parts += 1;
            }
        }
    }
    let mut lerr = false;
    if incomma && !irtype {
        lerr = true;
        warn!("The rtype=fixed but at least one non-null name has no leading and trailing range.");
    }
    if extra_parts != 0 && extra_parts != 5 {
        lerr = true;
        warn!("**** Error: Your names= list only has some of _cf _ct _ci _rf _rt. Need all 5 or none.");
    }
    if incomma && names_more {
        lerr = true;
        warn!("**** Error: C_SU_MORE not permitted after C_SU_NAMES for comma-separated files.");
    }
    let num_to_sort_by = matches.len();
    for n in 0..all_names {
        if !names[n].starts_with("null") && !names[n].starts_with("numb") {
            for i in 1..=9 {
                if names[n] == format!("match{}", i) && num_to_sort_by >= i {
                    names[n] = matches[i - 1].clone();
                }
                if names[n] == format!("matche{}", i) && num_to_sort_by >= i {
                    names[n] = matches[num_to_sort_by - i].clone();
                }
            }
            if names[n].starts_with("match") {
                lerr = true;
                warn!(
                    "**** Error: Name  %s  could not be substituted from match= list.",
                    names[n]
                );
            }
            for m in n + 1..all_names {
                if names[n] == names[m] && namex[n] == namex[m] {
                    lerr = true;
                    warn!(
                        "**** Error: Name  %s  exists at least twice in the names list.",
                        names[n]
                    );
                }
            }
        }
    }
    if lerr {
        err!("**** Error: Related to names= or C_SU_NAMES record (details above).");
    }
    let mut kcase = vec![0; num_to_sort_by];
    for i in 0..num_to_sort_by {
        kcase[i] = get_case(&matches[i]);
        if kcase[i] < 1 {
            err!("**** Error: a match name not recognized (or not allowed).");
        }
    }
    let mut ncase = Vec::<i32>::new();
    let mut nspot = Vec::<usize>::new();
    let mut valmx = Vec::<f64>::new();
    let mut jscalel = false;
    let mut jscalco = false;
    let mut joffset = false;
    let mut numcases = 0;
    for i in 0..all_names {
        let iamcase = get_case(&names[i]);
        if iamcase < 0 {
            err!(
                "**** Error: Name  %s  in the names list is not recognized.",
                names[i]
            );
        }
        if iamcase > 0 && iamcase < 1000 {
            ncase.push(iamcase);
            ilead[numcases] = ilead[i];
            itrail[numcases] = itrail[i];
            names[numcases] = names[i].clone();
            namex[numcases] = namex[i].clone();
            nspot.push(if irtype { i } else { numcases });
            let ktype = hdtype(&names[i]);
            let mut vmax = match ktype {
                "i" => 2147483645.0,
                "h" => 32765.0,
                "u" => 65533.0,
                "f" => 3.4e37,
                _ => 1.7e307,
            };
            // Keep names[i], as in C: the compaction above may have changed it.
            if ["gelev", "selev", "sdepth", "gdel", "sdel", "swdep", "gwdep"]
                .contains(&names[i].as_str())
            {
                vmax /= dscalel;
                jscalel = true;
            }
            if ["sx", "sy", "gx", "gy"].contains(&names[i].as_str()) {
                vmax /= dscalco;
                jscalco = true;
            }
            if names[i] == "offset" {
                joffset = true;
            }
            valmx.push(vmax);
            numcases += 1;
        }
    }
    names.truncate(numcases);
    namex.truncate(numcases);
    if !jscalel {
        if iscaleldef {
            iscalel = 0;
        } else {
            err!("**** Error: Cannot specify scalel= if not updating any of the 7 elevation related values.");
        }
    }
    if !jscalco {
        if iscalcodef {
            iscalco = 0;
        } else {
            err!("**** Error: Cannot specify scalco= if not updating any of the 4 coordinate values.");
        }
    }
    let ioffset = !joffset && jscalco;
    let mut mapx = [-1i32; 10];
    if extra_parts == 5 {
        let mut kerr = false;
        for i in 0..numcases {
            let pos = match namex[i].as_str() {
                "cf" => 0,
                "ct" => 1,
                "ci" => 2,
                "rf" => 3,
                "rt" => 4,
                _ => 10,
            };
            if pos < 5 {
                if mapx[pos] != -1 {
                    kerr = true;
                    match pos {
                        0 => warn!("**** Error: Only one name with _cf appended is allowed."),
                        1 => warn!("**** Error: Only one name with _ct appended is allowed."),
                        2 => warn!("**** Error: Only one name with _ci appended is allowed."),
                        3 => warn!("**** Error: Only one name with _rf appended is allowed."),
                        _ => warn!("**** Error: Only one name with _rt appended is allowed."),
                    }
                }
                mapx[pos] = i as i32;
                if pos == 1 || pos == 2 || pos == 4 {
                    ncase[i] = 0;
                }
            } else if names[i] == "grnofr" {
                mapx[7] = i as i32;
            } else if names[i] == "grnlof" {
                mapx[8] = i as i32;
            }
        }
        if mapx[..5].contains(&-1) {
            err!("**** Error: Your _cf _ct _ci _rf _rt specification is incorrect.");
        }
        if names[mapx[0] as usize] != names[mapx[1] as usize]
            || names[mapx[0] as usize] != names[mapx[2] as usize]
        {
            kerr = true;
            warn!("**** Error: _cf _ct _ci must be appended to the same name.");
        }
        if names[mapx[3] as usize] != names[mapx[4] as usize] {
            kerr = true;
            warn!("**** Error: _rf _rt must be appended to the same name.");
        }
        if kerr {
            err!("**** Error: Your _cf _ct _ci _rf _rt specification is incorrect.");
        }
        for k in 0..num_to_sort_by {
            if names[mapx[0] as usize] == matches[k] {
                mapx[5] = k as i32;
                mapx[6] = kcase[k];
            }
        }
        if mapx[5] == -1 {
            err!("**** Error: Cannot find the name appended with _cf _ct and _ci in match= list.");
        }
    }
    let mut ktol = vec![0; num_to_sort_by];
    for k in 0..num_to_sort_by {
        if let Some(n) = ncase.iter().position(|&n| n == kcase[k]) {
            ktol[k] = n;
        } else {
            err!("**** Error: match= name not found in names=");
        }
    }
    let klocn = ktol.clone();
    for n in 0..numcases {
        for k in 0..num_to_sort_by {
            if ncase[n] == kcase[k] {
                ncase[n] = 0;
            }
        }
    }

    if num_r < 1 {
        num_r = data_records(&fpr, &rid, lenid, nicerecord).count() as i32;
        fseeko(&fpr, 0, SEEK_SET);
        warn!("Counted %d data records.", num_r);
        if num_r == 0 {
            err!("**** No data records found. Wrong setid value? Wrong file?");
        }
    }
    let mut rec_info = Vec::<PointInfo>::new();
    let mut e = ReadErrors::default();
    let mut ur = Unrepeat::new(1);
    let mut count = 0;
    for (ncount, mut raw) in data_records(&fpr, &rid, lenid, nicerecord) {
        count += 1;
        if count > num_r {
            continue;
        }
        let mut rec = PointInfo {
            dfield: vec![0.0; numcases],
            lfield: vec![0; num_to_sort_by],
        };
        if !irtype {
            let bytes = raw.as_bytes();
            let mut textbeg = Vec::new();
            for ineed in 0..numcases {
                if itrail[ineed] as usize >= bytes.len() {
                    e.lenerr += 1;
                    if e.lenerr < 4 {
                        warn!(
                            "Error at record %d   Record-too-short for requested fied ranges",
                            ncount
                        );
                        if e.lenerr == 3 {
                            warn!("Have 3 Record-too-short warnings, no more will be printed.");
                        }
                    }
                    break;
                }
                textbeg
                    .extend_from_slice(&bytes[(ilead[ineed] - 1) as usize..itrail[ineed] as usize]);
                textbeg.push(rdel);
            }
            raw = String::from_utf8_lossy(&textbeg).into_owned();
        }
        get_csv(&raw, rdel, &mut rec.dfield, &nspot, ncount, &mut e);
        for j in 0..numcases {
            let d = rec.dfield[j];
            if d < 1.0e308 && (d > valmx[j] || 0.0 - d > valmx[j]) {
                e.lrgerr += 1;
                if e.lrgerr < 4 {
                    warn!(
                        "Error at record %d number-too-large for SU name (%.2f)",
                        ncount, d
                    );
                    if e.lrgerr == 3 {
                        warn!(
                            "Have 3 number-too-large-to-output warnings, no more will be printed."
                        );
                    }
                }
            }
        }
        if mapx[3] > -1 {
            let cfr = rec.dfield[mapx[0] as usize];
            let ctr = rec.dfield[mapx[1] as usize];
            let cir = rec.dfield[mapx[2] as usize];
            let rfr = rec.dfield[mapx[3] as usize];
            let rtr = rec.dfield[mapx[4] as usize];
            if cfr > ctr {
                rec.dfield[mapx[0] as usize] = ctr;
                rec.dfield[mapx[1] as usize] = cfr;
                rec.dfield[mapx[3] as usize] = rtr;
                rec.dfield[mapx[4] as usize] = rfr;
            }
            rec.dfield[mapx[2] as usize] = if (rec.dfield[mapx[0] as usize]
                - rec.dfield[mapx[1] as usize])
                .abs()
                < dtolh * 2.0
            {
                1.0
            } else {
                cir.abs()
            };
        }
        for k in 0..num_to_sort_by {
            rec.lfield[k] = longt(rec.dfield[ktol[k]], dtolh, dtol);
        }
        if unrepeat > -2147483645 {
            ur.apply(&mut rec.lfield[0]);
        }
        rec_info.push(rec);
    }
    if e.nblank > 0 {
        warn!(
            "Total all-blank fields: %d. Assumed zero for all.",
            e.nblank
        );
    }
    if e.numerr > 0 {
        warn!(
            "Total errors for Field-unreadable as a number:        %d",
            e.numerr
        );
    }
    if e.morerr > 0 {
        warn!(
            "Total errors for Two-numbers in one field:            %d",
            e.morerr
        );
    }
    if e.lenerr > 0 {
        warn!(
            "Total errors for Record-too-short to get all values:  %d",
            e.lenerr
        );
    }
    if e.comerr > 0 {
        warn!(
            "Total errors for Not-enough-commas to get all values: %d",
            e.comerr
        );
    }
    if e.lrgerr > 0 {
        warn!(
            "Total errors for Number-too-large for SU name:        %d",
            e.lrgerr
        );
    }
    if e.lenerr > 0 || e.comerr > 0 || e.morerr > 0 || e.numerr > 0 || e.lrgerr > 0 {
        err!("Errors detected while reading text file (see details above).");
    }
    if unrepeat > -2147483645 {
        warn!(
            "For unrepeat option, the text-based incrementing integer ended at: %d ",
            ur.nrep / INCINT
        );
    }
    warn!(
        "Have allocated memory to store values from %d records. Found %d records.",
        num_r, count
    );
    if count > num_r {
        err!(
            "Too many records read (%d) for your maxrecords= value (%d).",
            count,
            num_r
        );
    }
    rec_info.sort_by(|a, b| comp_sort(a, b, num_to_sort_by).cmp(&0));
    let lapover = check_records(&rec_info, &matches, &mapx, dtolh, dtol);
    let record_count = rec_info.len();
    if record_count == 0 {
        // With explicit maxrecords, C retains its calloc-initialized first
        // record even when numR becomes zero. bhigh searches zero records,
        // then the matching and creation paths still access RecInfo[0].
        rec_info.resize_with(num_r as usize, || PointInfo {
            dfield: vec![0.0; numcases],
            lfield: vec![0; num_to_sort_by],
        });
    }

    if ncreate < 1 || ireplicate > 0 {
        if gettr(&mut tr) == 0 {
            err!("can't read first trace");
        }
        for _ in 1..ireplicate {
            if gettr(&mut tr) == 0 {
                err!("can't read enough traces for spikes=n option.");
            }
        }
    } else if nspikes == 0 {
        tr.ns = 1;
        tr.dt = 4000;
        tr.data[0] = 0.0;
    } else {
        let mut tmax = 0.0f32;
        for n in (2..nspikes as usize).step_by(2) {
            if spikes[n] > tmax {
                tmax = spikes[n];
            }
        }
        let nsamps = nint((tmax / spikes[0]) as f64);
        tr.ns = nsamps as u16;
        tr.dt = nint(1000.0 * spikes[0] as f64) as u16;
        if nsamps > tr.data.len() as i32 {
            err!("Error: Spike times exceed trace sample capacity.");
        }
        for n in 0..nsamps {
            tr.data[n as usize] = spikes[1];
        }
        for n in 0..nspikes as usize / 2 {
            let m = nint((spikes[n * 2] / spikes[0]) as f64) - 1;
            if m < 0 {
                err!("Error: Spike times less than first sample are not permitted.");
            }
            if m as usize >= tr.data.len() {
                err!("Error: Spike times exceed trace sample capacity.");
            }
            tr.data[m as usize] = spikes[n * 2 + 1];
        }
    }
    ur = Unrepeat::new(unrepeat as i64);
    let mut notfound = 0;
    let mut nmade = 0;
    let mut nproct = 0;
    let mut dmade = 0.0;
    if mapx[3] > -1 {
        dmade = rec_info[0].dfield[mapx[0] as usize] - rec_info[0].dfield[mapx[2] as usize];
    }
    let mut igo = true;
    loop {
        if ncreate > 0 {
            if nproct == ncreate {
                break;
            }
            if nmade >= rec_info.len() {
                // Advancing past the allocated records is undefined in C.
                err!("**** No data records found. Wrong setid value? Wrong file?");
            }
            for n in 0..num_to_sort_by {
                tohead(&mut tr, kcase[n], rec_info[nmade].dfield[klocn[n]]);
            }
            if mapx[3] < 0 {
                nmade += 1;
                if nmade == record_count {
                    igo = false;
                }
            } else {
                dmade += rec_info[nmade].dfield[mapx[2] as usize];
                if dmade > rec_info[nmade].dfield[mapx[1] as usize] {
                    nmade += 1;
                    if nmade == record_count {
                        break;
                    }
                    if nmade >= rec_info.len() {
                        err!("**** No data records found. Wrong setid value? Wrong file?");
                    }
                    for n in 0..num_to_sort_by {
                        tohead(&mut tr, kcase[n], rec_info[nmade].dfield[klocn[n]]);
                    }
                    dmade = rec_info[nmade].dfield[mapx[0] as usize];
                }
                tohead(&mut tr, mapx[6], dmade);
            }
        }
        let mut dvals = vec![0.0; num_to_sort_by];
        let mut guy = PointInfo {
            dfield: Vec::new(),
            lfield: vec![0; num_to_sort_by],
        };
        for n in 0..num_to_sort_by {
            dvals[n] = fromhead(&tr, kcase[n]);
            guy.lfield[n] = longt(dvals[n], dtolh, dtol);
        }
        if unrepeat > -2147483645 {
            ur.apply(&mut guy.lfield[0]);
        }
        let mut idx_r = bhigh(&rec_info[..record_count], &guy, num_to_sort_by);
        if idx_r > 0 {
            idx_r -= 1;
        }
        let mut ifound = 1;
        if mapx[3] < 0 {
            if comp_sort(&rec_info[idx_r], &guy, num_to_sort_by) != 0 {
                ifound = -1;
            }
        } else if comp_sort(&rec_info[idx_r], &guy, num_to_sort_by - 1) != 0 {
            ifound = -2;
        } else if lapover == 0 {
            if guy.lfield[num_to_sort_by - 1]
                > longt(rec_info[idx_r].dfield[mapx[1] as usize], dtolh, dtol)
            {
                ifound = -3;
            }
        } else {
            ifound = -4;
            for m in (0..=idx_r).rev() {
                if comp_sort(&rec_info[m], &guy, num_to_sort_by - 1) != 0 {
                    break;
                }
                if guy.lfield[num_to_sort_by - 1]
                    > longt(rec_info[m].dfield[mapx[1] as usize], dtolh, dtol)
                {
                    continue;
                }
                let mcf = longt(rec_info[m].dfield[mapx[0] as usize], dtolh, dtol);
                let mci = longt(rec_info[m].dfield[mapx[2] as usize], dtolh, dtol);
                if guy.lfield[num_to_sort_by - 1].checked_rem(mci).unwrap_or(0)
                    == mcf.checked_rem(mci).unwrap_or(0)
                {
                    ifound = 2;
                    idx_r = m;
                    break;
                }
            }
        }
        if ifound > 0 {
            if iscalel != 0 {
                tr.scalel = iscalel as i16;
            }
            if iscalco != 0 {
                tr.scalco = iscalco as i16;
                tr.counit = 1;
            }
            for n in 0..numcases {
                let mut dvalue = rec_info[idx_r].dfield[n];
                if mapx[3] == n as i32 {
                    let cfo = rec_info[idx_r].dfield[mapx[0] as usize];
                    let cto = rec_info[idx_r].dfield[mapx[1] as usize];
                    let rfo = rec_info[idx_r].dfield[mapx[3] as usize];
                    let rto = rec_info[idx_r].dfield[mapx[4] as usize];
                    dvalue = if cto - cfo < dtolh * 2.0 {
                        rfo
                    } else {
                        rfo + (dvals[mapx[5] as usize] - cfo) * (rto - rfo) / (cto - cfo)
                    };
                }
                if ncase[n] > 12 && ncase[n] < 20 {
                    dvalue *= dscalel;
                }
                if ncase[n] > 21 && ncase[n] < 26 {
                    dvalue *= dscalco;
                }
                if iaction != 0 {
                    dvalue = iaction as f64 * dvalue + fromhead(&tr, ncase[n]);
                }
                tohead(&mut tr, ncase[n], dvalue);
            }
            if ioffset {
                let dx = tr.sx.wrapping_sub(tr.gx) as f64;
                let dy = tr.sy.wrapping_sub(tr.gy) as f64;
                tr.offset =
                    (dx / dscalco * dx / dscalco + dy / dscalco * dy / dscalco).sqrt() as i32;
            }
        }
        if ifound < 0 {
            notfound += 1;
            if missing == -1 {
                err!(" Error: Could not find a match= text record for trace and missing= option is error-halt.");
            }
        }
        if ifound > 0 || missing == 1 {
            puttr(&mut tr);
            nproct += 1;
        }
        if ncreate < 1 && gettr(&mut tr) == 0 {
            igo = false;
        }
        if !igo {
            break;
        }
    }
    if unrepeat > -2147483645 {
        warn!(
            "For unrepeat option, the trace-based incrementing integer ended at: %d ",
            ur.nrep / INCINT
        );
    }
    warn!("Number of traces output: %d ", nproct);
    if notfound > 0 {
        if missing == 0 {
            warn!(
                "Traces not output due to no match= in text file: %d (your missing= option)",
                notfound
            );
        } else {
            warn!(
                "Traces output anyway despite no match= in text file: %d (your missing= option)",
                notfound
            );
        }
    }
    exit(CWP_Exit());
}

fn check_records(
    rec: &[PointInfo],
    matches: &[String],
    mapx: &[i32; 10],
    dtolh: f64,
    dtol: f64,
) -> i32 {
    let num = matches.len();
    let mut lapover = 0;
    if mapx[3] > -1 {
        let mut verr = [0i32; 8];
        let mut ntop = 0;
        let mut npchan = 0;
        let mut npsegs = 0;
        let mut rpinc = 0.0;
        let d = |m: usize, k: usize| rec[m].dfield[mapx[k] as usize];
        let l = |m: usize, k: usize| longt(d(m, k), dtolh, dtol);
        for n in 1..=rec.len() {
            if n == rec.len() || comp_sort(&rec[ntop], &rec[n], num - 1) != 0 {
                let mut same = [false; 8];
                let mut nchan = 0i32;
                let mut nsegs = 0;
                let mut top_grnofr = 0;
                let mut top_grnlof = 0;
                for m in ntop..n {
                    let mcf = l(m, 0);
                    let mct = l(m, 1);
                    let mut mci = l(m, 2);
                    if mci == 0 {
                        mci = 1;
                    }
                    // arm64 integer division returns zero for a zero divisor.
                    nchan =
                        nchan.wrapping_add(((mct - mcf).checked_div(mci).unwrap_or(0) + 1) as i32);
                    nsegs += 1;
                    let mut rinc = 1.0;
                    if mct != mcf {
                        rinc = (d(m, 3) - d(m, 4))
                            / ((mct - mcf).checked_div(mci).unwrap_or(0) as f64);
                    }
                    let j = (rec[m].lfield[0] as f64 / dtol + 0.5) as i32;
                    if ntop > 0 && !same[4] && (rinc - rpinc).abs() > dtolh * 2.0 {
                        same[4] = true;
                        verr[4] += 1;
                        if verr[4] < 4 {
                            warn!(
                                "Warning: Receiver points per channel changed at %s= %d ",
                                matches[0], j
                            );
                            if verr[4] == 3 {
                                warn!("Have 3 Receiver-points-per-channel changed warnings, no more will be printed.");
                            }
                        }
                    }
                    rpinc = rinc;
                    if !same[0]
                        && mcf.checked_rem(mci).unwrap_or(0) != mct.checked_rem(mci).unwrap_or(0)
                    {
                        same[0] = true;
                        verr[0] += 1;
                        if verr[0] < 4 {
                            warn!(
                                "Error: Layout ends do not conform to same increment at %s= %d ",
                                matches[0], j
                            );
                            if verr[0] == 3 {
                                warn!("Have 3 Layout-ends-do-not-conform-to-same-increment errors, no more will be printed.");
                            }
                        }
                    }
                    if mapx[7] > -1 {
                        if m == ntop {
                            top_grnofr = l(m, 7);
                        } else if !same[6] && top_grnofr != l(m, 7) {
                            same[6] = true;
                            verr[6] += 1;
                            if verr[6] < 4 {
                                warn!("Error: The grnofr values are different for same shot at %s= %d ", matches[0], j);
                                if verr[6] == 3 {
                                    warn!("Have 3 grnofr-values-are-different-for-same-shot, no more will be printed.");
                                }
                            }
                        }
                    }
                    if mapx[8] > -1 {
                        if m == ntop {
                            top_grnlof = l(m, 8);
                        } else if !same[7] && top_grnlof != l(m, 8) {
                            same[7] = true;
                            verr[7] += 1;
                            if verr[7] < 4 {
                                warn!("Error: The grnlof values are different for same shot at %s= %d ", matches[0], j);
                                if verr[7] == 3 {
                                    warn!("Have 3 grnlof-values-are-different-for-same-shot, no more will be printed.");
                                }
                            }
                        }
                    }
                    // The C loop deliberately stops at n-1, leaving the last
                    // segment out of the pairwise overlap check.
                    for i in m + 1..n - 1 {
                        let icf = l(i, 0);
                        let ici = l(i, 2);
                        if icf > mct {
                            break;
                        }
                        let j = (rec[i].lfield[0] as f64 / dtol + 0.5) as i32;
                        if lapover < 1 {
                            warn!("Warning: Overlapping channel range at %s= %d  Unusual, but not always an error.", matches[0], j);
                        }
                        lapover += 1;
                        if !same[1] && ici != mci {
                            same[1] = true;
                            verr[1] += 1;
                            if verr[1] < 4 {
                                warn!(
                                    "Error: Different increments in overlapping layout at %s= %d ",
                                    matches[0], j
                                );
                                if verr[1] == 3 {
                                    warn!("Have 3 Different-increments-in-overlapping layout errors, no more will be printed.");
                                }
                            }
                        }
                        if !same[2]
                            && mcf.checked_rem(mci).unwrap_or(0)
                                == icf.checked_rem(ici).unwrap_or(0)
                        {
                            same[2] = true;
                            verr[2] += 1;
                            if verr[2] < 4 {
                                warn!(
                                    "Error: Overlapping layout hits same channels at %s= %d ",
                                    matches[0], j
                                );
                                if verr[2] == 3 {
                                    warn!("Have 3 Overlapping-layout-hits-same-channels errors, no more will be printed.");
                                }
                            }
                        }
                    }
                }
                let j = (rec[ntop].lfield[0] as f64 / dtol + 0.5) as i32;
                if ntop > 0 && !same[3] && nchan != npchan {
                    verr[3] += 1;
                    if verr[3] < 4 {
                        warn!(
                            "Warning: Number of channels in layout changed at %s= %d ",
                            matches[0], j
                        );
                        if verr[3] == 3 {
                            warn!("Have 3 Number-of-channels-in-layout-changed warnings, no more will be printed.");
                        }
                    }
                }
                npchan = nchan;
                if ntop > 0 && !same[5] && nsegs != npsegs {
                    verr[5] += 1;
                    if verr[5] < 4 {
                        if nsegs == npsegs * 2 || nsegs * 2 == npsegs {
                            warn!("Warning: Number of segments in layout changed by ratio 2 at %s= %d (see unrepeat=).", matches[0], j);
                        } else {
                            warn!(
                                "Warning: Number of segments in layout changed at %s= %d ",
                                matches[0], j
                            );
                        }
                        if verr[5] == 3 {
                            warn!("Have 3 Number-of-segments-in-layout-changed warnings, no more will be printed.");
                        }
                    }
                }
                npsegs = nsegs;
                ntop = n;
            }
        }
        if lapover > 0 {
            warn!(
                "Total Overlapping-channel-ranges in layouts:        %d (very unusual)",
                lapover
            );
        }
        if verr[1] > 0 {
            warn!(
                "Total Different-increments-in-overlapping layout:   %d (error-halt)",
                verr[1]
            );
        }
        if verr[2] > 0 {
            warn!(
                "Total Overlapping-layout-hits-same-channels:        %d (error-halt)",
                verr[2]
            );
        }
        if verr[0] > 0 {
            warn!(
                "Total Layout-ends-do-not-conform-to-same-increment: %d (error-halt)",
                verr[0]
            );
        }
        if verr[3] > 0 {
            warn!(
                "Total Number-of-channels-in-layout-changed:         %d (unusual)",
                verr[3]
            );
        }
        if verr[4] > 0 {
            warn!(
                "Total Receiver-points-per-channel changed:          %d (very unusual)",
                verr[4]
            );
        }
        if verr[5] > 0 {
            warn!(
                "Total Number-of-segments-in-layout-changed:         %d (unusual)",
                verr[5]
            );
        }
        if verr[6] > 0 {
            warn!(
                "Total grnofr-values-are-different-for-same-shot:    %d (error-halt)",
                verr[6]
            );
        }
        if verr[7] > 0 {
            warn!(
                "Total grnlof-values-are-different-for-same-shot:    %d (error-halt)",
                verr[7]
            );
        }
        if verr[0] > 0 || verr[1] > 0 || verr[2] > 0 || verr[6] > 0 || verr[7] > 0 {
            err!("Errors detected while reading text file (see details above).");
        }
    } else {
        let mut l7verr = 0;
        for n in 1..rec.len() {
            if comp_sort(&rec[n], &rec[n - 1], num) == 0 {
                l7verr += 1;
                if l7verr < 4 {
                    let j = (rec[n].lfield[0] as f64 / dtol + 0.5) as i32;
                    if num == 1 {
                        warn!(
                            "Error: Records have duplicate match values of %s=%d",
                            matches[0], j
                        );
                    } else {
                        let k = (rec[n].lfield[1] as f64 / dtol + 0.5) as i32;
                        if num == 2 {
                            warn!(
                                "Error: Records have duplicate match values of %s=%d   %s=%d",
                                matches[0], j, matches[1], k
                            );
                        } else {
                            let l = (rec[n].lfield[2] as f64 / dtol + 0.5) as i32;
                            if num == 3 {
                                warn!("Error: Records have duplicate match values of %s=%d   %s=%d   %s=%d",matches[0],j,matches[1],k,matches[2],l);
                            } else {
                                let m = (rec[n].lfield[3] as f64 / dtol + 0.5) as i32;
                                warn!("Error: Records have duplicate match values of %s=%d   %s=%d   %s=%d   %s=%d",matches[0],j,matches[1],k,matches[2],l,matches[3],m);
                            }
                        }
                    }
                    if l7verr == 3 {
                        warn!("Have 3 Records-have-duplicate-match-values errors, no more will be printed.");
                    }
                }
            }
        }
        if l7verr > 0 {
            err!(
                "Total  errors  for Records-have-duplicate-match-values: %d",
                l7verr
            );
        }
        if num == 1 || num == 2 {
            let mut jp = 0;
            let mut jd = 0;
            let mut jc = 1;
            let mut je = -1;
            let mut kp = 0;
            let mut kd = 0;
            let mut kc = -1;
            for n in 0..rec.len() {
                let jn = if num == 2 {
                    (rec[n].lfield[0] as f64 / dtol + 0.5) as i32
                } else {
                    0
                };
                let kn = (rec[n].lfield[num - 1] as f64 / dtol + 0.5) as i32;
                if n == 0 {
                    jp = jn;
                    kp = kn;
                }
                if jn != jp {
                    jc += 1;
                    if jd != jn - jp {
                        je += 1;
                    }
                    jd = jn - jp;
                } else if kn - kp != kd {
                    kc += 1;
                    kd = kn - kp;
                }
                jp = jn;
                kp = kn;
            }
            if num == 2 {
                warn!(
                    "Note: There are: %d sets of %s values (lines?).",
                    jc, matches[0]
                );
                if je > 1 {
                    warn!("Warning: There are: %d irregular %s increments between the sets (missing lines?).",(je+1)/2,matches[0]);
                }
                if kc > 1 {
                    warn!("Warning: There are: %d irregular %s increments within the lines (missing points?).",(kc+1)/2,matches[1]);
                }
            } else if kc > 1 {
                warn!("Warning: There are: %d irregular %s increments within the line (missing points?).",(kc+1)/2,matches[0]);
            }
        }
    }
    lapover
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                          ",
    " SUGEOMCSV - 3D/2D land/OBotm geometry from Comma Separated or fixed text.",
    "                                                                          ",
    " sugeomcsv <intraces >outtraces rfile=in.txt   (other parameters)         ",
    "                                                                          ",
    " Parameter overview:                                                      ",
    "                                                                          ",
    "       rfile= input text values file                                      ",
    "       ufile= text file containing C_SU options (default is same rfile)   ",
    "       rtype= type of records in input text file (comma separated, fixed) ",
    "       names= assign SU names to text values (with SPS2 and SPS1 options) ",
    "       setid= accept text records based on first characters (X,S,R etc.)  ",
    "       match= trace SU key values to be used to find correct text records ",
    "      rdelim= rfile delimiter when rtype is csv (default is comma).       ",
    "  nicerecord= allows skipping bad records at start of some text files     ",
    "  maxrecords= maximum number of records to allocate memory for            ",
    "    unrepeat= help when duplicate fldr values exist (default is off)      ",
    "      scalco= set scaler for sx,sy,gx,gy coordinate values (default 10)   ",
    "      scalel= set scaler for elevation and related values (default 10)    ",
    "     missing= options when trace match= values not found in text file     ",
    "      action= how to update trace header values (default is to replace)   ",
    "      create= create output traces with no input trace file               ",
    "      spikes= put spikes at specified times in the created output traces  ",
    "                                                                          ",
    "  **********************************************************              ",
    "   To output this documentation:  sugeomcsv 2> geomdoc.txt                ",
    "  **********************************************************              ",
    "                                                                          ",
    " This program updates headers from general csv or fixed format text files.",
    " (You cannot use UKOOA streamer-marine R-records as input here).          ",
    "                                                                          ",
    " This program can update full geometry using the 3 SPS fixed format files.",
    " Parameter setting for SPS2 and SPS1 is so complex that this program has  ",
    " explicit option names=sps2 and names=sps1 for those files. But the SPS   ",
    " files often have numbers that are too large to be stored in some SU keys.",
    " Or, the SPS files may be unusable for other reasons. In these cases, you ",
    " need to use the SUTOOLCSV program to repair the SPS files or to convert  ",
    " those files to comma-separated files for input to SpreadSheet programs   ",
    " in order to make other repairs before running this program.              ",
    "                                                                          ",
    " Refer to the documentation of SUTOOLCSV for SPS files since it has       ",
    " extensive explanations and examples.                                     ",
    "                                                                          ",
    "                                                                          ",
    " For SPS2 files, 3 runs of this program are needed. Initially, I suggest  ",
    " using create=all option for the X-file. This reads no input traces but   ",
    " creates traces as defined by the X-file records (see create= below).     ",
    "                                                                          ",
    " sugeomcsv rfile=X.txt setid=x create=all >dummy1.su names=sps2 match=sps2",
    "                                                                          ",
    " sugeomcsv rfile=S.txt setid=s <dummy1.su >dummy2.su names=sps2 match=sps2",
    "                                                                          ",
    " sugeomcsv rfile=R.txt setid=r <dummy2.su >dummy3.su names=sps2 match=sps2",
    "                                                                          ",
    "                                                                          ",
    " If the SPS2 files have been processed by SUTOOLCSV then these 3 runs     ",
    " can be done instead.                                                     ",
    "                                                                          ",
    " sugeomcsv rfile=X.csv create=all >dummy1.su                              ",
    " sugeomcsv rfile=S.csv <dummy1.su >dummy2.su                              ",
    " sugeomcsv rfile=R.csv <dummy2.su >dummy3.su                              ",
    "                                                                          ",
    " ----------------------------------------------------------------------   ",
    "                                                                          ",
    " The documentation herein covers simple text files, but mostly concerns   ",
    " itself with other aspects of the use of SUGEOMCSV.                       ",
    "                                                                          ",
    " In other words, the documentation for SUGEOMCSV is complex enough that   ",
    " I do not want to repeat all the documentation contained in SUTOOLCSV.    ",
    "                                                                          ",
    "                                                                          ",
    " -------------                                                            ",
    "                                                                          ",
    " Consider a comma separated text file with 3 values per record.           ",
    " The 3 values are: record id, energy source number, and elevation.        ",
    "                                                                          ",
    "       S,11,343                                                           ",
    "       S,42,342                                                           ",
    "       S,25,340                                                           ",
    "       C,some comment,                                                    ",
    "       S,45,347                                                           ",
    "                                                                          ",
    "  To update shot elevation to trace headers requires this:                ",
    "                                                                          ",
    "  sugeomcsv <intraces >outtraces rfile=values3.csv setid=S match=es       ",
    "            names=C_su_id,es,selev                                        ",
    "                                                                          ",
    "  setid=S says only read values from text records that start with S.      ",
    "  match=es says get the value of es from input trace header.              ",
    " names=C_su_id,es,selev assigns names to the 3 values in the_text_file    ",
    " Other details for the above parameters will be discussed later.          ",
    "                                                                          ",
    " So, the value of es is gotten from input trace, then the text record     ",
    " with same es value is found, and the selev is updated to output trace.   ",
    " Note that every trace with the same es value is updated with same selev. ",
    "   You can use multiple matches to find text records (match=sx,sy).       ",
    " If the text records contain a fourth value with depth then you can       ",
    " update it at the same time using (names=C_su_id,es,selev,sdepth).        ",
    " If you want to update just sdepth then (names=C_su_id,es,null,sdepth)    ",
    "                                                                          ",
    " Note that trace input order does not matter. Neither does text record    ",
    " order. The text file is read and all non-null and non-numb name values   ",
    " are stored and sorted by match=. But there is no interpolation, so the   ",
    " match= numbers must be found in text records (within a small tolerance). ",
    "                                                                          ",
    " Now consider the same 3 values but in a fixed format text file:          ",
    "                                                                          ",
    "       S       11  343                                                    ",
    "       S       42  342                                                    ",
    "       S       25  340                                                    ",
    "       C  some comment                                                    ",
    "       S       45  347                                                    ",
    "                                                                          ",
    "  To update shot elevation to trace headers requires this:                ",
    "                                                                          ",
    "  sugeomcsv <intraces >outtraces rfile=values3.txt setid=S match=es       ",
    "            names=C_su_id,2_es_10,11_selev_15                             ",
    "                                                                          ",
    "  The first difference is that the rfile now has .txt rather than .csv    ",
    "  The other difference is in names=C_su_id,2_es_10,11_selev_15 where      ",
    "  both es and selev have leading and trailing numbers indicating the      ",
    "  character range that they occupy in the fixed format records of rfile.  ",
    "  Note that C_su_id does not have the leading and trailing numbers since  ",
    "  it is special in several ways (detailed later).                         ",
    "                                                                          ",
    "  When parsing fixed format files it is not actually needed to specify    ",
    "  fields that you decide to null. But it is easier to make changes in     ",
    "  the future if you retain all known fields in the names= lists. This is  ",
    "  why SUTOOLCSV produces C_SU records that contain all officially defined ",
    "  fields for SPS2 and SPS1 formats (even though it nulls some names).     ",
    "                                                                          ",
    "  --------------                                                          ",
    "                                                                          ",
    "  Alternately, instead of putting parameters on the command line,         ",
    "  you can add similar parameter records into the rfile.                   ",
    "  Details follow later, but here is an example:                           ",
    "                                                                          ",
    "  sugeomcsv <intraces >outtraces rfile=values3a.txt                       ",
    "                                                                          ",
    "       C_SU_MATCH,es                                                      ",
    "       C_SU_SETID,S                                                       ",
    "       C_SU_NAMES                                                         ",
    "       C_su_id,2_es_10,11_selev_15                                        ",
    "       S       11  343                                                    ",
    "       S       42  342                                                    ",
    "       S       25  340                                                    ",
    "       C  some comment                                                    ",
    "       S       45  347                                                    ",
    "                                                                          ",
    "  Note that the default is to look for C_SU_ parameter records within     ",
    "  the same rfile as the values records. And match= setid= names= are      ",
    "  optional (but their C_SU_ records must be found if they are not on the  ",
    "  command line). If they are specified on the command line and also exist ",
    "  in the rfile, their command line options override rfile.                ",
    "                                                                          ",
    "  Note that commas are still specified in the C_SU_ parameter records     ",
    "  even when the rest of the rfile is fixed format.                        ",
    "                                                                          ",
    "  --------------                                                          ",
    "                                                                          ",
    "  For comma separated rfiles, it looks like this:                         ",
    "                                                                          ",
    "  sugeomcsv <intraces >outtraces rfile=values3a.csv                       ",
    "                                                                          ",
    "       C_SU_MATCH,es,                                                     ",
    "       C_SU_SETID,S,                                                      ",
    "       C_SU_NAMES,,                                                       ",
    "       C_su_id,2_es_10,11_selev_15                                        ",
    "       S,11,343                                                           ",
    "       S,42,342                                                           ",
    "       S,25,340                                                           ",
    "       C,some comment,                                                    ",
    "       S,45,347                                                           ",
    "                                                                          ",
    "  Notice that the C_SU_NAMES are allowed to have leading and trailing     ",
    "  integers even though the data records have commas. When rtype=csv they  ",
    "  are ignored in this program (but you should retain leading and trailing ",
    "  numbers because SUTOOLCSV uses them to convert csv back to fixed files).",
    "                                                                          ",
    "  Notice also that C_SU_MATCH and C_SU_SETID have extra commas. This may  ",
    "  occur for comma-separated files output by SpreadSheet programs. These   ",
    "  extra commas are also permitted for C_SU_MATCH and C_SU_SETID records   ",
    "  in fixed format files. (The C_SU_MATCH read-in stops at the first null, ",
    "  all-blank, or sequential commas. And C_SU_SETID only reads first id).   ",
    "                                                                          ",
    "  If you are wondering why the actual names are in the row following      ",
    "  C_SU_NAMES, consider how this will be output from a SpreadSheet program.",
    "  (The 3 actual names need to be in the same columns with the data that   ",
    "  they refer to).                                                         ",
    "                                                                          ",
    "  --------------                                                          ",
    "                                                                          ",
    "  For simple text files, using the command line to specify everything may ",
    "  be more convenient. But specification for SPS and similar files gets    ",
    "  quite intricate and is usually better done via C_SU_ parameter records. ",
    "  For this reason, the SUTOOLCSV program actually outputs the high-level  ",
    "  SPS2 and SPS1 options in expanded version into its output text files.   ",
    "  This give you the opportunity to see and modify them as needed.         ",
    "                                                                          ",
    "  You can also specify ufile= to read the parameter records from another  ",
    "  file. This allows you to use the same ufile for different surveys.      ",
    "                                                                          ",
    " ------------------------------------------------------------------------ ",
    "                                                                          ",
    " For full initial geometry this program usually must be run multiple times",
    " For instance, if you have 3 SPS files, this program must be run 3 times. ",
    " The SPS X-file must be used in first run, then S then R (or R then S).   ",
    "                                                                          ",
    " The first run of this program inputs the SPS X-file and then uses the    ",
    " trace header fldr and tracf numbers to find the correct X-record inside  ",
    " the X-file. The computed gaps numbers is updated to the output trace     ",
    " (gaps is the SU key name I decided to use, it is usually called          ",
    " receiver point or receiver number in SPS and other documents).           ",
    " Other values from the X-record are also updated (these are usually       ",
    " receiver line number and shot line number and shot point number.         ",
    " For 3D surveys, in the next 2 runs of SUGEOMCSV:                         ",
    "      Two of the values output from the X-file are used                   ",
    "      to find the correct records within the S-file.                      ",
    "      Two of the values output from the X-file are used                   ",
    "      to find the correct records within the R-file.                      ",
    " (For 2D surveys, in the next 2 runs of SUGEOMCSV only 1 value may be     ",
    "  needed to find the correct records for S-file and 1 for R-file).        ",
    "                                                                          ",
    " NOTE: The SPS discussions above assume that fldr contains field record   ",
    "       number, and tracf contains channel number. This is not necessarily ",
    "       the case. In particular, the channel number might not be in tracf. ",
    "       It might be in some other SU key name, or it might not exist       ",
    "       at all, in which case you are going to have to figure out how to   ",
    "       alter some input trace key name to correspond to the channel       ",
    "       ranges on the X-records.                                           ",
    "                                                                          ",
    " NOTE: Sometimes SPS X-files deliberately specify receiver line numbers   ",
    "       or point ranges that do not exist in the R-files. For instance,    ",
    "       with a 10 line by 240 point pattern layout, it is sometimes easier ",
    "       for recording instruments to simply pretend that the line numbering",
    "       and point numbering continues infinitely. The traces corresponding ",
    "       to this are simply allowed to record noise. Othertimes the traces  ",
    "       are not recorded at all. For this situaton you need to specify     ",
    "       missing=delete or pass (or pre-delete the traces).                 ",
    "                                                                          ",
    " ------------------------------------------------------------------------ ",
    "                                                                          ",
    " Required parameters:                                                     ",
    "                                                                          ",
    "       rfile=  text values file (fixed or comma separated values)         ",
    "                                                                          ",
    " Optional parameters:                                                     ",
    "                                                                          ",
    "       rtype= type of records in rfile text file. Default is csv if the   ",
    "              file name ends in csv, otherwise defaults to fixed.         ",
    "            =csv     comma separated values                               ",
    "            =fixed   This option is required for standard SPS files       ",
    "                     that have not been converted to csv by SUTOOLCSV     ",
    "                     - along with other specifications in names= list     ",
    "                     (see extensive examples below).                      ",
    "                                                                          ",
    "       ufile=  file to search for C_SU parameter records. The default is  ",
    "               to search in the rfile text file.                          ",
    "            =command   means ignore any C_SU records in rfile text file.  ",
    "            All required parameters must be specified on the command line.",
    "                                                                          ",
    " The following 3 parameters must be found on the command line or in       ",
    " their corresponding C_SU records in the rfile (or ufile) text file.      ",
    "                                                                          ",
    "       match= any number of input trace header keys needed to find exact  ",
    "             record in the rfile text file. Must also be listed in names= ",
    "             Example on command line:                                     ",
    "             match=fldr,tracf                                             ",
    "       match=SPS2 This is just a standard way to specify the match= list  ",
    "                  for SPS Revison 2 files (see examples in SUTOOLCSV).    ",
    "                  The setid= option must also be X,S, or R.               ",
    "       match=SPS1 This is just a standard way to specify the match= list  ",
    "                  for SPS Revison 1 files (see examples in SUTOOLCSV).    ",
    "                  The setid= option must also be X,S, or R.               ",
    "                                                                          ",
    "       names= is used to assign names to values in rfile text file (you   ",
    "              are telling this program the names of values in text file). ",
    "             For files with comma-separated values, a name must be listed ",
    "             sequentially for each field in the rfile text file.          ",
    "             The names must also include the match= SU keys above.        ",
    "             Note C_su_id means this is field used for record acceptance. ",
    "        ***  Read the note   c_su_id IS SPECIAL   later. ***              ",
    "             Special names: null and numb with integer (such as null1)    ",
    "                           Means do not read/output this field. (You can  ",
    "                           also put nothing between sequential commas).   ",
    "               Example on command line:                                   ",
    "               names=C_su_id,cdp,null3,cx,cy,,ce                          ",
    "       names=SPS2 This is just a standard way to specify the names= list  ",
    "                  for SPS Revison 2 files (see examples in SUTOOLCSV).    ",
    "                  The setid= option must also be X,S, or R.               ",
    "       names=SPS1 This is just a standard way to specify the names= list  ",
    "                  for SPS Revison 1 files (see examples in SUTOOLCSV).    ",
    "                  The setid= option must also be X,S, or R.               ",
    "       names=SPS2ALL  Same as SPS2. The SUTOOLCSV program has this option ",
    "                      to output fields with no SU key (using numb names), ",
    "                      but numb names are treated like null names herein.  ",
    "       names=SPS1ALL  Same as SPS1.                                       ",
    "                                                                          ",
    "      setid= is used to accept text records based on their first field.   ",
    "        setid=S     means accept text records if their first field is S   ",
    "                          (any characters allowed, such as R,X,cdp,FRED)  ",
    "                    Note: this value is automatically upper-cased unless  ",
    "                          you surround it by double-quotes.               ",
    "                          So s becomes S unless you use double-quotes.    ",
    "        setid=ANY   means read all records (except those starting C_SU)   ",
    "                    and those records have an id field at front.          ",
    "        setid=NONE  means read all records (except those starting C_SU)   ",
    "                    but those records do not have an id field at front.   ",
    "                    (For csv files this means the field before the first  ",
    "                     comma is a value, not an identifier).                ",
    "             Example on command line:                                     ",
    "             setid=S                                                      ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "       If match= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_MATCH and reads keys from it. ",
    "             Example within the text file:                                ",
    "             C_SU_MATCH,fldr,tracf                                        ",
    "                                                                          ",
    "       If names= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_NAMES and reads names from    ",
    "       the record after the C_SU_NAMES record.                            ",
    "             Example within the text file:                                ",
    "             C_SU_NAMES                                                   ",
    "             C_su_id,cdp,null,cx,cy,null,ce                               ",
    "                                                                          ",
    "       If setid= is not specified on command line, this program searches  ",
    "       for a text record starting with C_SU_SETID and reads id from it.   ",
    "             Example within the text file:                                ",
    "             C_SU_SETID,S                                                 ",
    "                                                                          ",
    " Note these C_SU_ parameter records can be in any order within text file. ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "      rdelim= If rtype=fixed you cannot specify this parameter.           ",
    "              If rtype=csv the default is comma. You can specify any      ",
    "              single character here either by itself or surrounded by     ",
    "              double-quotes (some characters such as semi-colon may have  ",
    "              trouble getting through the command line).                  ",
    "     **Note** Specifying a blank here usually will not give good results  ",
    "              unless your input rfile has exactly 1 blank between fields. ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "      scalco= multiply coordinates by this power of 10 (1,10,100...)      ",
    "              before putting them in traces. Default is 10 which means    ",
    "              that sx,sy,gx,gy from the text file are multiplied by 10.   ",
    "              The actual value of scalco in the traces is therefore       ",
    "              set to -10 (meaning these values need to be divided by 10). ",
    "          *** You cannot specify scalco= if you are not updating any      ",
    "          *** of the 4 coordinate values.                                 ",
    "         **** The SU key offset value is recomputed if any of sx,sy,gx,gy ",
    "         **** is in the list of names (and offset itself is NOT in list). ",
    "            * If you are confident your text files contain coordinates    ",
    "            * with only whole numbers, you can set this to 1. But you     ",
    "            * cannot change this just for some of these 4 coordinates.    ",
    "           ** This SEEMS like a problem about size, but it is actually a  ",
    "           ** problem about decimal digits. It you use 1 here you will    ",
    "           ** find that values like 544444.6 get rounded to 544445        ",
    "           ** when SUGEOMCSV updates them to traces.                      ",
    "          *** If scalco is in the list of names= then that value is set   ",
    "              in output header. But only the value specified or defaulted ",
    "              HERE is actually applied. So do not put it in names= unless ",
    "              you are repairing some odd situation.                       ",
    "                                                                          ",
    "      scalel= multiply elevation and other related values by this power   ",
    "              of 10 (1,10,100...) before putting them in traces. Default  ",
    "              is 10 which means gelev,selev,sdepth,gdel,sdel,swdep,gwdep  ",
    "              from the text file are multiplied by 10.                    ",
    "              The actual value of scalel in the traces is therefore       ",
    "              set to -10 (meaning these values need to be divided by 10). ",
    "          *** You cannot specify scalel= if you are not updating any      ",
    "          *** of the 7 elevation related values.                          ",
    "            * If you are confident your text files contain these values   ",
    "            * with only whole numbers, you can set this to 1. But you     ",
    "            * cannot change this just for some of these 7 values.         ",
    "           ** This SEEMS like a problem about size, but it is actually a  ",
    "           ** problem about decimal digits. It you use 1 here you will    ",
    "           ** find that values like 3333.6 get rounded to 3334            ",
    "           ** when SUGEOMCSV updates them to traces.                      ",
    "          *** If scalel is in the list of names= then that value is set   ",
    "              in output header. But only the value specified or defaulted ",
    "              HERE is actually applied. So do not put it in names= unless ",
    "              you are repairing some odd situation.                       ",
    "                                                                          ",
    "     unrepeat= The default is not to enable this option.                  ",
    "               This option is general but most likely usefull for X-files ",
    "               where the field record number (fldr) increases but then    ",
    "               re-starts at a lower number. Such as 1->7800 then 5->4000. ",
    "               Normally, the finding-logic within this program would not  ",
    "               be able to distinguish the first fldr 5 from the second 5  ",
    "               and so on. (For that situation if you do not use this      ",
    "               option you will most likely get a channel-range error as   ",
    "               that X-file is read-in and stored).                        ",
    "      unrepeat=1 Read the text file and generate an integer from 1 and    ",
    "                 increment by 1 every time the first match= reverses.     ",
    "                 Typically, the first match= is fldr for X-files so this  ",
    "                 increments +1 when fldr is increasing and then decreases,",
    "                 and also increments +1 if fldr is decreasing and then    ",
    "                 increases. The comparison is done using order of records ",
    "                 as they exist in the text file (before sorting herein).  ",
    "                 Another incrementing integer is generated the same way   ",
    "                 except using the order of the input traces. These two    ",
    "                 integers are used to match which (fldr) value in the     ",
    "                 traces belongs to which (fldr) value from X-file.        ",
    "             **  Note that this option is not fool-proof and should be    ",
    "             **  used with caution and copious checking of results. This  ",
    "             **  option is better than a pure sequential re-numbering     ",
    "             **  because it still primarily uses the field record numbers,",
    "             **  but it is not guaranteed.                                ",
    "      unrepeat=n (any integer). The input text record incrementing integer",
    "                 will still start at 1 and increment +1. But the trace    ",
    "                 incrementing integer starts at this number. This allows  ",
    "                 you to input reduced sets of fldr records without need   ",
    "                 to edit the X-file.                                      ",
    "                                                                          ",
    " -----------------                                                        ",
    " -----------------                                                        ",
    "                                                                          ",
    "     nicerecord= record number to start trying to read data values from   ",
    "                 the rfile text file (default is 1). The beginning records",
    "                 of some text files are badly composed (comments or info.)",
    "                 When the setid= option is not able to reject them,       ",
    "                 specify a record number here where setid= will work.     ",
    "                 (This program also always knows that records starting    ",
    "                 with C_SU are not data records, and will not try to      ",
    "                 read data values from them even when setid=ALL). But     ",
    "                 it will read C_SU parameter records even if they are     ",
    "                 previous to this nicerecord number.                      ",
    "                                                                          ",
    "     maxrecords= maximum number of records to allocate memory for.        ",
    "                 If not specified, this program reads through the records ",
    "                 once and allocates for the number found. Then reads them ",
    "                 again. This double reading takes more time. If you want  ",
    "                 to avoid this, specify a maximum via this parameter.     ",
    "                                                                          ",
    "    missing= specifies what to do if an input trace match= values cannot  ",
    "             be found in the rfile text file. Default is to error-halt    ",
    "             at the first trace where this occurs.                        ",
    "     missing=delete   means delete all traces for which this occurs.      ",
    "     missing=pass     means pass those traces to output without update.   ",
    "                      This allows you to, for instance, change statics    ",
    "                      values for just some shots or receivers without     ",
    "                      having every shot or receiver in the rfile text file",
    "                                                                          ",
    "    action= specifies how to update the trace header values               ",
    "            Default is to set the headers to the text record values.      ",
    "            The option here is not applied to the values of               ",
    "            names on the match= list.                                     ",
    "     action=add       Add the text values to the trace header values.     ",
    "     action=subtract  Subtract text values from the trace header values.  ",
    "                                                                          ",
    " -----------------                                                        ",
    " -----------------                                                        ",
    "                                                                          ",
    "    create=    create output traces and update from the rfile text file.  ",
    "          ***  Attempting to also input a trace file will error-halt      ",
    "          ***  unless spikes=n).                                          ",
    "               Values of all non-null names on names= list are updated to ",
    "               traces (including those on match= list). The behavior of   ",
    "               this program is the same as when traces are input. This    ",
    "               includes sorting of the text records into match= order and ",
    "               therefore created traces are output in that order. Most    ",
    "               text files result in one trace per record, but SPS X-files ",
    "               (and similar) usually create many traces for each record.  ",
    "  create=all   create traces for all data records in the rfile text file  ",
    "  create=n     create n or all (whichever is fewer)                       ",
    "                                                                          ",
    "  spikes= ***  This option will error-halt if create= is not specified.   ",
    "               By default, create= outputs only 1 zero-sample at 4 ms.    ",
    "               This parameter accepts a list of time,amplitude pairs.     ",
    "                                                                          ",
    "                 Example: spikes=4,0.001,800,1000,1500,-2000,3000,0       ",
    "                                                                          ",
    "               The maximum time sets the trace length (3000 ms. above).   ",
    "               The first pair in the list is special, it specifies the    ",
    "               sample interval (in ms.) and base-amplitude.               ",
    "               In this example, the sample interval is 4, and the         ",
    "               base-amplitude is 0.001. The base-amplitude is the value   ",
    "               samples are set to if they are not spiked by remainder of  ",
    "               the list. (For geometry checking and other QC and survey   ",
    "               design reasons, base-amplitude of 0.0 is not a nice value. ",
    "               For instance, 0.0 makes it difficult to check mutes).      ",
    "               All times should be whole multiples of the sample interval ",
    "               but will be rounded to nearest sample time if not.         ",
    "               Times do not have to increase but cannot be negative.      ",
    "     spikes=n  Where n is the sequential number of an input trace. The    ",
    "               sample values of this trace will be copied to all created  ",
    "               output traces. Values of all non-null names on names= list ",
    "               are updated to traces (including those on the match= list).",
    "               The intention here is to allow you to input a reflectivity ",
    "               or well-log type of trace, duplicate the samples and change",
    "               some header values via the text file (for testing purposes)",
    "                                                                          ",
    " -----------------                                                        ",
    " -----------------                                                        ",
    "                                                                          ",
    " ADVICE: Run create=all for your first relation file (SPS X-file) setup.  ",
    "         Then check the values in the output trace headers.               ",
    "         Then run SPS S-file setup to update those headers (not create).  ",
    "         Then check the values in the output trace headers.               ",
    "         Then run SPS R-file setup to update those headers (not create).  ",
    "         Then check the values in the output trace headers.               ",
    "         This will allow you to check the mechanics of your setups.       ",
    "                                                                          ",
    " -----------------                                                        ",
    "                                                                          ",
    "  MORE ADVICE: Use spike options to test your setups of other programs.   ",
    "               For instance, after running the sequence advised above,    ",
    "               put positive and negative spikes on the traces.            ",
    "               Then run a simple bandpass filter. At this point the       ",
    "               results will begin to look like post-NMO shot gathers.     ",
    "               Apply inverse-NMO and offset-Mute to make your traces look ",
    "               like pre-NMO shot gathers.                                 ",
    "               Then you can run some interesting tests, like:             ",
    "           1.  Apply surface-consistent residual statics to some shots and",
    "               receivers. Surface-consistant statics analysis programs    ",
    "               should easily be able to find the statics you applied.     ",
    "           2.  Put surface-consistant gains on some shots and receivers.  ",
    "               Surface-consistant amplitude analysis programs should      ",
    "               easily find the gains you applied.                         ",
    "           3.  Put surface-consistent signatures on some shots and        ",
    "               receivers. Surface-consistent deconvolution programs       ",
    "               should easily find the signatures you applied.             ",
    "                                                                          ",
    "   CAUTION: Just because you apply something and another program computes ",
    "            the correct inverse, it does not mean the program is setup    ",
    "            well for actual seismic situations. These kinds of tests are  ",
    "            too simple to be considered forward-modelling.                ",
    "                                                                          ",
    "   OF COURSE: You might do this kind of thing without having any actual   ",
    "              seismic data, just to test a few survey-design ideas.       ",
    "                                                                          ",
    "                                                                          ",
    " ------------------------------------------------------------------------ ",
    " ------------------------------------------------------------------------ ",
    "                                                                          ",
    " C_su_id IS SPECIAL ********************                                  ",
    "                                                                          ",
    "  (a) Its character range is taken from the length of the value specified ",
    "      for setid (for instance S has range 1 to 1, FRED has range 1 to 4). ",
    "  (b) The value specified for id is case-sensitive (r is not R). The id   ",
    "      value is the only thing case-sensitive in this program except for   ",
    "      the file names. So this program does not care if you use parameter  ",
    "      records starting with C_SU or c_su, but you want other programs to  ",
    "      ignore C_SU records, so use a capital C, not a lower case c.        ",
    "                                                                          ",
    " -----------------------------------------------------------------        ",
    "                                                                          ",
];

fn get_case(cbuf: &str) -> i32 {
    if cbuf.starts_with("null") {
        return 0;
    }
    if let Some(s) = cbuf.strip_prefix("numb") {
        return 1000 + atoi(s);
    }
    match cbuf {
        "tracl" => 1,
        "tracr" => 2,
        "fldr" => 3,
        "tracf" => 4,
        "ep" => 5,
        "cdp" => 6,
        "cdpt" => 7,
        "trid" => 8,
        "nvs" => 9,
        "nhs" => 10,
        "duse" => 11,
        "offset" => 12,
        "gelev" => 13,
        "selev" => 14,
        "sdepth" => 15,
        "gdel" => 16,
        "sdel" => 17,
        "swdep" => 18,
        "gwdep" => 19,
        "scalel" => 20,
        "scalco" => 21,
        "sx" => 22,
        "sy" => 23,
        "gx" => 24,
        "gy" => 25,
        "counit" => 26,
        "wevel" => 27,
        "swevel" => 28,
        "sut" => 29,
        "gut" => 30,
        "sstat" => 31,
        "gstat" => 32,
        "tstat" => 33,
        "laga" => 34,
        "lagb" => 35,
        "delrt" => 36,
        "muts" => 37,
        "mute" => 38,
        "ns" => 39,
        "dt" => 40,
        "gain" => 41,
        "igc" => 42,
        "igi" => 43,
        "corr" => 44,
        "sfs" => 45,
        "sfe" => 46,
        "slen" => 47,
        "styp" => 48,
        "stas" => 49,
        "stae" => 50,
        "tatyp" => 51,
        "afilf" => 52,
        "afils" => 53,
        "nofilf" => 54,
        "nofils" => 55,
        "lcf" => 56,
        "hcf" => 57,
        "lcs" => 58,
        "hcs" => 59,
        "year" => 60,
        "day" => 61,
        "hour" => 62,
        "minute" => 63,
        "sec" => 64,
        "timbas" => 65,
        "trwf" => 66,
        "grnors" => 67,
        "grnofr" => 68,
        "grnlof" => 69,
        "gaps" => 70,
        "otrav" => 71,
        "d1" => 72,
        "f1" => 73,
        "d2" => 74,
        "f2" => 75,
        "ungpow" => 76,
        "unscale" => 77,
        "ntr" => 78,
        "mark" => 79,
        _ => -1,
    }
}

fn fromhead(tr: &Trace, k: i32) -> f64 {
    match k {
        1 => tr.tracl as f64,
        2 => tr.tracr as f64,
        3 => tr.fldr as f64,
        4 => tr.tracf as f64,
        5 => tr.ep as f64,
        6 => tr.cdp as f64,
        7 => tr.cdpt as f64,
        8 => tr.trid as f64,
        9 => tr.nvs as f64,
        10 => tr.nhs as f64,
        11 => tr.duse as f64,
        12 => tr.offset as f64,
        13 => tr.gelev as f64,
        14 => tr.selev as f64,
        15 => tr.sdepth as f64,
        16 => tr.gdel as f64,
        17 => tr.sdel as f64,
        18 => tr.swdep as f64,
        19 => tr.gwdep as f64,
        20 => tr.scalel as f64,
        21 => tr.scalco as f64,
        22 => tr.sx as f64,
        23 => tr.sy as f64,
        24 => tr.gx as f64,
        25 => tr.gy as f64,
        26 => tr.counit as f64,
        27 => tr.wevel as f64,
        28 => tr.swevel as f64,
        29 => tr.sut as f64,
        30 => tr.gut as f64,
        31 => tr.sstat as f64,
        32 => tr.gstat as f64,
        33 => tr.tstat as f64,
        34 => tr.laga as f64,
        35 => tr.lagb as f64,
        36 => tr.delrt as f64,
        37 => tr.muts as f64,
        38 => tr.mute as f64,
        39 => tr.ns as f64,
        40 => tr.dt as f64,
        41 => tr.gain as f64,
        42 => tr.igc as f64,
        43 => tr.igi as f64,
        44 => tr.corr as f64,
        45 => tr.sfs as f64,
        46 => tr.sfe as f64,
        47 => tr.slen as f64,
        48 => tr.styp as f64,
        49 => tr.stas as f64,
        50 => tr.stae as f64,
        51 => tr.tatyp as f64,
        52 => tr.afilf as f64,
        53 => tr.afils as f64,
        54 => tr.nofilf as f64,
        55 => tr.nofils as f64,
        56 => tr.lcf as f64,
        57 => tr.hcf as f64,
        58 => tr.lcs as f64,
        59 => tr.hcs as f64,
        60 => tr.year as f64,
        61 => tr.day as f64,
        62 => tr.hour as f64,
        63 => tr.minute as f64,
        64 => tr.sec as f64,
        65 => tr.timbas as f64,
        66 => tr.trwf as f64,
        67 => tr.grnors as f64,
        68 => tr.grnofr as f64,
        69 => tr.grnlof as f64,
        70 => tr.gaps as f64,
        71 => tr.otrav as f64,
        72 => tr.d1 as f64,
        73 => tr.f1 as f64,
        74 => tr.d2 as f64,
        75 => tr.f2 as f64,
        76 => tr.ungpow as f64,
        77 => tr.unscale as f64,
        78 => tr.ntr as f64,
        79 => tr.mark as f64,
        80 => tr.shortpad as f64,
        _ => 0.0,
    }
}

fn tohead(tr: &mut Trace, k: i32, dvalue: f64) {
    match k {
        1 => tr.tracl = (dvalue.round_ties_even() as i64) as _,
        2 => tr.tracr = (dvalue.round_ties_even() as i64) as _,
        3 => tr.fldr = (dvalue.round_ties_even() as i64) as _,
        4 => tr.tracf = (dvalue.round_ties_even() as i64) as _,
        5 => tr.ep = (dvalue.round_ties_even() as i64) as _,
        6 => tr.cdp = (dvalue.round_ties_even() as i64) as _,
        7 => tr.cdpt = (dvalue.round_ties_even() as i64) as _,
        8 => tr.trid = (dvalue.round_ties_even() as i64) as _,
        9 => tr.nvs = (dvalue.round_ties_even() as i64) as _,
        10 => tr.nhs = (dvalue.round_ties_even() as i64) as _,
        11 => tr.duse = (dvalue.round_ties_even() as i64) as _,
        12 => tr.offset = (dvalue.round_ties_even() as i64) as _,
        13 => tr.gelev = (dvalue.round_ties_even() as i64) as _,
        14 => tr.selev = (dvalue.round_ties_even() as i64) as _,
        15 => tr.sdepth = (dvalue.round_ties_even() as i64) as _,
        16 => tr.gdel = (dvalue.round_ties_even() as i64) as _,
        17 => tr.sdel = (dvalue.round_ties_even() as i64) as _,
        18 => tr.swdep = (dvalue.round_ties_even() as i64) as _,
        19 => tr.gwdep = (dvalue.round_ties_even() as i64) as _,
        20 => tr.scalel = (dvalue.round_ties_even() as i64) as _,
        21 => tr.scalco = (dvalue.round_ties_even() as i64) as _,
        22 => tr.sx = (dvalue.round_ties_even() as i64) as _,
        23 => tr.sy = (dvalue.round_ties_even() as i64) as _,
        24 => tr.gx = (dvalue.round_ties_even() as i64) as _,
        25 => tr.gy = (dvalue.round_ties_even() as i64) as _,
        26 => tr.counit = (dvalue.round_ties_even() as i64) as _,
        27 => tr.wevel = (dvalue.round_ties_even() as i64) as _,
        28 => tr.swevel = (dvalue.round_ties_even() as i64) as _,
        29 => tr.sut = (dvalue.round_ties_even() as i64) as _,
        30 => tr.gut = (dvalue.round_ties_even() as i64) as _,
        31 => tr.sstat = (dvalue.round_ties_even() as i64) as _,
        32 => tr.gstat = (dvalue.round_ties_even() as i64) as _,
        33 => tr.tstat = (dvalue.round_ties_even() as i64) as _,
        34 => tr.laga = (dvalue.round_ties_even() as i64) as _,
        35 => tr.lagb = (dvalue.round_ties_even() as i64) as _,
        36 => tr.delrt = (dvalue.round_ties_even() as i64) as _,
        37 => tr.muts = (dvalue.round_ties_even() as i64) as _,
        38 => tr.mute = (dvalue.round_ties_even() as i64) as _,
        39 => tr.ns = (dvalue.round_ties_even() as i64) as _,
        40 => tr.dt = (dvalue.round_ties_even() as i64) as _,
        41 => tr.gain = (dvalue.round_ties_even() as i64) as _,
        42 => tr.igc = (dvalue.round_ties_even() as i64) as _,
        43 => tr.igi = (dvalue.round_ties_even() as i64) as _,
        44 => tr.corr = (dvalue.round_ties_even() as i64) as _,
        45 => tr.sfs = (dvalue.round_ties_even() as i64) as _,
        46 => tr.sfe = (dvalue.round_ties_even() as i64) as _,
        47 => tr.slen = (dvalue.round_ties_even() as i64) as _,
        48 => tr.styp = (dvalue.round_ties_even() as i64) as _,
        49 => tr.stas = (dvalue.round_ties_even() as i64) as _,
        50 => tr.stae = (dvalue.round_ties_even() as i64) as _,
        51 => tr.tatyp = (dvalue.round_ties_even() as i64) as _,
        52 => tr.afilf = (dvalue.round_ties_even() as i64) as _,
        53 => tr.afils = (dvalue.round_ties_even() as i64) as _,
        54 => tr.nofilf = (dvalue.round_ties_even() as i64) as _,
        55 => tr.nofils = (dvalue.round_ties_even() as i64) as _,
        56 => tr.lcf = (dvalue.round_ties_even() as i64) as _,
        57 => tr.hcf = (dvalue.round_ties_even() as i64) as _,
        58 => tr.lcs = (dvalue.round_ties_even() as i64) as _,
        59 => tr.hcs = (dvalue.round_ties_even() as i64) as _,
        60 => tr.year = (dvalue.round_ties_even() as i64) as _,
        61 => tr.day = (dvalue.round_ties_even() as i64) as _,
        62 => tr.hour = (dvalue.round_ties_even() as i64) as _,
        63 => tr.minute = (dvalue.round_ties_even() as i64) as _,
        64 => tr.sec = (dvalue.round_ties_even() as i64) as _,
        65 => tr.timbas = (dvalue.round_ties_even() as i64) as _,
        66 => tr.trwf = (dvalue.round_ties_even() as i64) as _,
        67 => tr.grnors = (dvalue.round_ties_even() as i64) as _,
        68 => tr.grnofr = (dvalue.round_ties_even() as i64) as _,
        69 => tr.grnlof = (dvalue.round_ties_even() as i64) as _,
        70 => tr.gaps = (dvalue.round_ties_even() as i64) as _,
        71 => tr.otrav = (dvalue.round_ties_even() as i64) as _,
        72 => tr.d1 = dvalue as _,
        73 => tr.f1 = dvalue as _,
        74 => tr.d2 = dvalue as _,
        75 => tr.f2 = dvalue as _,
        76 => tr.ungpow = dvalue as _,
        77 => tr.unscale = dvalue as _,
        78 => tr.ntr = (dvalue.round_ties_even() as i64) as _,
        79 => tr.mark = (dvalue.round_ties_even() as i64) as _,
        80 => tr.shortpad = (dvalue.round_ties_even() as i64) as _,
        _ => {}
    }
}

fn sps_names(kind: &str, rid: &str) -> Option<Vec<String>> {
    let mut result: Option<Vec<String>> = match (kind, rid) {
        ("sps2" | "sps2all", "X") => Some(vec![
            "c_su_id".into(),
            "2_null2_7".into(),
            "8_match1_15".into(),
            "16_null4_16".into(),
            "17_null5_17".into(),
            "18_grnofr_27".into(),
            "28_grnlof_37".into(),
            "38_null8_38".into(),
            "39_matche1_cf_43".into(),
            "44_matche1_ct_48".into(),
            "49_matche1_ci_49".into(),
            "50_grnors_59".into(),
            "60_gaps_rf_69".into(),
            "70_gaps_rt_79".into(),
            "80_null15_80".into(),
        ]),
        ("sps2" | "sps2all", "S") => Some(vec![
            "c_su_id".into(),
            "2_grnofr_11".into(),
            "12_grnlof_21".into(),
            "22_null4_23".into(),
            "24_null5_24".into(),
            "25_null6_26".into(),
            "27_sstat_30".into(),
            "31_sdepth_34".into(),
            "35_sdel_38".into(),
            "39_sut_40".into(),
            "41_swdep_46".into(),
            "47_sx_55".into(),
            "56_sy_65".into(),
            "66_selev_71".into(),
            "72_null15_74".into(),
            "75_null16_76".into(),
            "77_null17_78".into(),
            "79_null18_80".into(),
        ]),
        ("sps2" | "sps2all", "R") => Some(vec![
            "c_su_id".into(),
            "2_grnors_11".into(),
            "12_gaps_21".into(),
            "22_null4_23".into(),
            "24_null5_24".into(),
            "25_null6_26".into(),
            "27_gstat_30".into(),
            "31_null8_34".into(),
            "35_gdel_38".into(),
            "39_gut_40".into(),
            "41_gwdep_46".into(),
            "47_gx_55".into(),
            "56_gy_65".into(),
            "66_gelev_71".into(),
            "72_null15_74".into(),
            "75_null16_76".into(),
            "77_null17_78".into(),
            "79_null18_80".into(),
        ]),
        ("sps1" | "sps1all", "X") => Some(vec![
            "c_su_id".into(),
            "2_null2_7".into(),
            "8_match1_11".into(),
            "12_null4_12".into(),
            "13_null5_13".into(),
            "14_grnofr_29".into(),
            "30_grnlof_37".into(),
            "38_null8_38".into(),
            "39_matche1_cf_42".into(),
            "43_matche1_ct_46".into(),
            "47_matche1_ci_47".into(),
            "48_grnors_63".into(),
            "64_gaps_rf_71".into(),
            "72_gaps_rt_79".into(),
            "80_null15_80".into(),
        ]),
        ("sps1" | "sps1all", "S") => Some(vec![
            "c_su_id".into(),
            "2_grnofr_17".into(),
            "18_grnlof_25".into(),
            "26_null4_26".into(),
            "26_null5_26".into(),
            "27_null6_28".into(),
            "29_sstat_32".into(),
            "33_sdepth_36".into(),
            "37_sdel_40".into(),
            "41_sut_42".into(),
            "43_swdep_46".into(),
            "47_sx_55".into(),
            "56_sy_65".into(),
            "66_selev_71".into(),
            "72_null15_74".into(),
            "75_null16_76".into(),
            "77_null17_78".into(),
            "79_null18_80".into(),
        ]),
        ("sps1" | "sps1all", "R") => Some(vec![
            "c_su_id".into(),
            "2_grnors_17".into(),
            "18_gaps_25".into(),
            "26_null4_26".into(),
            "26_null5_26".into(),
            "27_null6_28".into(),
            "29_gstat_32".into(),
            "33_null8_36".into(),
            "37_gdel_40".into(),
            "41_gut_42".into(),
            "43_gwdep_46".into(),
            "47_gx_55".into(),
            "56_gy_65".into(),
            "66_gelev_71".into(),
            "72_null15_74".into(),
            "75_null16_76".into(),
            "77_null17_78".into(),
            "79_null18_80".into(),
        ]),
        _ => None,
    };
    // These four X-preset strings fill their C allocations exactly. Their
    // terminating NULs are overwritten by the next strcpy. On the reference
    // build the first name therefore contains the following three names too.
    // Preserve the observed parse error instead of silently repairing the SPS
    // presets. Explicit names= lists do not have this allocation bug.
    if rid == "X" {
        if let Some(names) = &mut result {
            names[8] = names[8..12].concat();
        }
    }
    result
}
