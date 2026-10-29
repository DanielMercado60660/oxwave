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

//! SUNEARCSV - Assign Trace Key Values From Nearest Q-file Record.
//! Port of su/main/headers/sunearcsv.c (Revision 1.01).
//! Credits: Created July 2022: Andre Latour. Started from suprofcsv.

use oxwave_core::prelude::*;
use oxwave_header_case::*;
use oxwave_q_file::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUNEARCSV - Assign Trace Key Values From Nearest Q-file Record.            ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sunearcsv [parameters].                                                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " qin=      Input file containing q-records. This parameter is required.     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " dimx=     Dimension name X. (defaults to not used). Name must be a key and ",
    "           must exist in the qin file, or be a special name below:          ",
    "             mgx = use mid X (sx+gx)/2 from trace, but gx from qin file     ",
    "             msx = use mid X (sx+gx)/2 from trace, but sx from qin file     ",
    "             mgy = use mid Y (sy+gy)/2 from trace, but gy from qin file     ",
    "             msy = use mid Y (sy+gy)/2 from trace, but sy from qin file     ",
    "             mgaps   = use (grnlof+gaps)/2 from trace, but gaps from qin    ",
    "             mgrnlof = use (grnlof+gaps)/2 from trace, but grnlof from qin  ",
    " dimy=     Dimension name Y. (defaults to not used). Same options as dimx.  ",
    " dimr=     Dimension name R. (defaults to not used). Same options as dimx.  ",
    " dima=     Dimension name A. (defaults to not used). Same options as dimx.  ",
    " dimb=     Dimension name B. (defaults to not used). Same options as dimx.  ",
    " dimc=     Dimension name C. (defaults to not used). Same options as dimx.  ",
    " dimd=     Dimension name D. (defaults to not used). Same options as dimx.  ",
    " dime=     Dimension name E. (defaults to not used). Same options as dimx.  ",
    " dimf=     Dimension name F. (defaults to not used). Same options as dimx.  ",
    "     Note: At least 1 dimension must be specified.                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The following 3 parameters can be used for any of the dimensions.\t     ",
    " Just substitute x,y,r,a,b,c,d,e,f as the ending character,                 ",
    " for instance: typer,minr,maxr for dimension R.                             ",
    " Typically these 3 are used just to specify an additional dimension that is ",
    " used to restrict the search range (see typer options -1 and -2 below).     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " typer=1 Use in Pythagorean Nearest (squared difference to trace location), ",
    "         and minr,maxr are unchanging extent ranges for this dimension.     ",
    "      =2 Use in Pythagorean Nearest and minr,maxr are the relative extent   ",
    "         range to include for this dimension (the input trace value for     ",
    "         this dimension is added to minr,maxr to form the extent range).    ",
    "     =-1 Do not use in Pythagorean Nearest, and minr,maxr are unchanging    ",
    "         extent ranges for this dimension.                                  ",
    "     =-2 Do not use in Pythagorean Nearest and minr,maxr are the relative   ",
    "         extent range to include for this dimension (the input trace value  ",
    "         for this dimension is added to minr,maxr to form extent range).    ",
    "    Note: Negative types means the difference between the profile point and ",
    "       the trace is NOT added to Pythagorean sum for specified dimension(s).",
    "       So, the nearest point is determined as if that dimension was NOT     ",
    "       specified at all. However, the min,max for that dimension are still  ",
    "       used. This results in finding the nearest point considering only the ",
    "       type>0 dimensions, but restricted by ranges of type<0 dimension(s).  ",
    "       Type=-2 can be used for situations where a profile approaches itself ",
    "       or intersects itself or overlaps itself. In those case, a third      ",
    "       dimension value (such as station number) can be used to restrict the ",
    "       search range for each trace to only the part of the profile with     ",
    "       approximately the same station numbers as the trace midpoint station.",
    " minr=     Extent Range Min. Typically negative. Default is no min limit.   ",
    "           Greater or equal to this value is in range.                      ",
    " maxr=     Extent Range Max. Typically positive. Default is no max limit.   ",
    "           Strictly less than this value is in range. The min must be less  ",
    "           than max, but the range does not have to be symmetric or centred ",
    "           (that is, ranges such as min=-1000 and max=-200 are allowed).    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " keyp=cdp  Point Order Key Name. Name must be in the input Q-file.          ",
    "           This parameter is the difference between a PROFILE and just a    ",
    "           set of points. After the nearest point to a trace is found, this ",
    "           parameter determines which point is considered the NEXT point.   ",
    "           The angle between the nearest point and the next point is used   ",
    "           in the computation of inline distance to output in the igi key   ",
    "           and the crossline distance to output in the igc key.             ",
    "           (The igi,igc distances are NOT scaled to the trace scalco flag). ",
    "     Note: To compute the angle, whatever is specified as dimx is used as   ",
    "           the X coordinate, whatever is specified as dimy is used as Y.    ",
    "     =asis The next point to determine angle is the next record in q-file.  ",
    "     =dist Set igi to 0 and igc to the distance between the trace and the   ",
    "           nearest point. The igc distance is NOT scaled by trace scalco.   ",
    "     =none Do not reset the igi and igc output values. This option must be  ",
    "           used if dimx or dimy is not specified.                           ",
    "     Note: Sometimes there are two-or-more points which are equally near    ",
    "           the trace. In those cases, the nearest point is set to the       ",
    "           lower ordered point. For options asis,dist,none the order is     ",
    "           simply the order of points in the qin file.                      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " okeys=cdp Output key names. Name must be in the input Q-file.              ",
    "           These values are copied from qin file to the output traces.      ",
    "           Note: q-files are assumed to contain raw (non-scaled) values.    ",
    "                 The scalco,scalel key values from the input traces are     ",
    "                 used to scale their related keys in this list.             ",
    "           Note: scalco,scalel key names are not allowed in this list.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " nopoint=1 If no qin point is found within the extents, error-halt (before  ",
    "           outputting the first trace that could not find a qin point).     ",
    "        =0 Continue without any update of the okeys values for traces that  ",
    "           could not find a qin point in extents (print their count at end).",
    "   Advice: If you choose to use option 0, then you should probably include  ",
    "           a key name in the qin file which has a value that you can check  ",
    "           in the output traces (start with a 0 in the input traces, have a ",
    "           value of 1 in the qin file, so you can see if trace was updated).",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " The following 3 parameters affect cpu time, but not results. The search    ",
    " is done by building a kdtree from the dimension values in the qin file.    ",
    " The code that builds the kdtree is reasonably standard and simplistic.     ",
    " But the kdtree search code is slightly unusual due to some typer options.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    " sfunc=2   Search via the cycle_for_near function. This option is usually   ",
    "           fastest. This option uses the sdist and smult parameters.        ",
    "      =1   Search via the find_near function. This option may be faster     ",
    "           if you are specifying small extent ranges.                       ",
    "           This option does not use the sdist and smult parameters.         ",
    "      =0   Search via the brute_near function. This option may be faster    ",
    "           if there are only a small number of points in the qin file.      ",
    "           This option does not use the sdist and smult parameters.         ",
    "           (To keep the code simpler, this option still allocates and       ",
    "            builds the kdtree. But the tree is not used for searching).     ",
    "      =-1  Search via the find_near function and confirm the results using  ",
    "           the brute_near function. This tests the find_near function and   ",
    "           can also be used to determine the speed difference between them. ",
    "           This option does not use the sdist and smult parameters.         ",
    "      =-2  Search via the cycle_for_near function and confirm the results   ",
    "           using brute_near function. This tests cycle_for_near function    ",
    "           and can also be used to determine speed difference between them. ",
    "           This option uses the sdist and smult parameters.                 ",
    " sdist=100 Initial search distance. This parameter is only used             ",
    "           if sfunc=2 or -2. If positive, this value is added to the        ",
    "           distance between the previous trace and its nearest qin point    ",
    "           and used as the initial search distance for the current trace.   ",
    "           If negative, the initial search distance for all traces is set   ",
    "           to this absolute value.                                          ",
    " smult=2   Search multiplier. This parameter is only used if sfunc=2 or -2. ",
    "           If the nearest qin point is not found after searching with the   ",
    "           initial search distance, the search distance is multiplied by    ",
    "           this value, and search is performed again. This repeats until    ",
    "           finding the nearest point (or all min,max ranges are exceeded).  ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "   ------------------------------------------------------------------       ",
    "\t\t\t\t\t\t\t\t\t     ",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut tr = Trace::new();
    let mut keyn = Vec::new();
    let mut kase = Vec::new();
    let mut in_rel_extent = [1i32; 9];
    let mut extent_in_min = [-f64::MAX; 9];
    let mut extent_in_max = [f64::MAX; 9];
    let mut tree_nt = [0i32; 9];
    let mut extent_min = [0.0; 9];
    let mut extent_max = [0.0; 9];
    let mut target = [0.0; 9];
    let mut locx = None;
    let mut locy = None;

    // Preserve the C parameter-check order: target[0:2] is XY.
    for (suffix, missing) in [
        (
            'x',
            "**** Error: Parameter typex, minx, or maxx specified, but not dimx.",
        ),
        (
            'y',
            "**** Error: Parameter typey, miny, or maxy specified, but not dimy.",
        ),
        (
            'r',
            "**** Error: Parameter typer, minr, or maxr specified, but not dimr.",
        ),
        (
            'a',
            "**** Error: Parameter typea, mina, or maxa specified, but not dima.",
        ),
        (
            'b',
            "**** Error: Parameter typeb, minb, or maxb specified, but not dimb.",
        ),
        (
            'c',
            "**** Error: Parameter typec, minc, or maxc specified, but not dimc.",
        ),
        (
            'd',
            "**** Error: Parameter typed, mind, or maxd specified, but not dimd.",
        ),
        (
            'e',
            "**** Error: Parameter typee, mine, or maxe specified, but not dime.",
        ),
        (
            'f',
            "**** Error: Parameter typef, minf, or maxf specified, but not dimf.",
        ),
    ] {
        let dim = format!("dim{suffix}");
        let typ = format!("type{suffix}");
        let min = format!("min{suffix}");
        let max = format!("max{suffix}");
        if countparval(&dim) > 0 {
            let i = keyn.len();
            let mut name = String::new();
            getparstring(&dim, &mut name);
            keyn.push(name);
            if suffix == 'x' {
                locx = Some(i);
            }
            if suffix == 'y' {
                locy = Some(i);
            }
            if countparval(&typ) > 0 {
                getparint(&typ, &mut in_rel_extent[i]);
            }
            if countparval(&min) > 0 {
                getpardouble(&min, &mut extent_in_min[i]);
            }
            if countparval(&max) > 0 {
                getpardouble(&max, &mut extent_in_max[i]);
            }
        } else if countparval(&typ) > 0 || countparval(&min) > 0 || countparval(&max) > 0 {
            err!(missing);
        }
    }
    let tree_numd = keyn.len();
    if tree_numd < 1 {
        err!("**** Error: At least 1 dimension name must be specified.");
    }
    for i in 0..tree_numd {
        extent_min[i] = extent_in_min[i];
        extent_max[i] = extent_in_max[i];
        if in_rel_extent[i] > 0 {
            tree_nt[i] = 1;
        } else {
            tree_nt[i] = 0;
            in_rel_extent[i] = in_rel_extent[i].wrapping_neg();
        }
        let k = match keyn[i].as_str() {
            "sx" => 101,
            "sy" => 102,
            "gx" => 103,
            "gy" => 104,
            "msx" => {
                keyn[i] = "sx".into();
                105
            }
            "msy" => {
                keyn[i] = "sy".into();
                106
            }
            "mgrnlof" => {
                keyn[i] = "grnlof".into();
                107
            }
            "mgx" => {
                keyn[i] = "gx".into();
                105
            }
            "mgy" => {
                keyn[i] = "gy".into();
                106
            }
            "mgaps" => {
                keyn[i] = "gaps".into();
                107
            }
            "selev" => 111,
            "gelev" => 112,
            "sdepth" => 113,
            "sdel" => 114,
            "gdel" => 115,
            "swdep" => 116,
            "gwdep" => 117,
            _ => {
                let k = GetCase(&keyn[i]);
                if k < 1 {
                    err!(
                        "**** Error: Specified dimension name %s is not recognized.",
                        keyn[i]
                    );
                }
                k
            }
        };
        kase.push(k);
    }
    let mut use_keyp = 1;
    let mut keyp = String::new();
    if countparval("keyp") > 0 {
        getparstring("keyp", &mut keyp);
        match keyp.as_str() {
            "asis" => use_keyp = 0,
            "dist" => use_keyp = -1,
            "none" => use_keyp = -2,
            _ => {}
        }
        if use_keyp != -2 && (locx.is_none() || locy.is_none()) {
            err!("**** Error: keyp=none must be specified if dimx or dimy is not specified.");
        }
    } else {
        keyp = "cdp".into();
    }
    let mut okeys = Vec::new();
    if countparval("okeys") > 0 {
        getparstringarray("okeys", &mut okeys);
    } else {
        okeys.push("cdp".into());
    }
    let mut okase = Vec::new();
    for key in &okeys {
        let k = GetCase(key);
        if k < 1 {
            err!(
                "**** Error: Specified okeys name %s is not recognized.",
                key
            );
        }
        if key == "scalel" || key == "scalco" {
            err!("**** Error: Specified okeys name %s is not allowed.", key);
        }
        okase.push(k);
    }
    let mut nopoint = 1i32;
    if getparint("nopoint", &mut nopoint) == 0 {
        nopoint = 1;
    }
    if nopoint > 1 || nopoint < 0 {
        err!("**** Error: nopoint must be 1 or 0.");
    }
    let mut sfunc = 2i32;
    if getparint("sfunc", &mut sfunc) == 0 {
        sfunc = 2;
    }
    if sfunc < -2 || sfunc > 2 {
        err!("**** Error: sfunc option is out-of-range.");
    }
    let scheck = sfunc < 0;
    if scheck {
        sfunc = -sfunc;
    }
    let mut sdist = 100.0f64;
    if getpardouble("sdist", &mut sdist) == 0 {
        sdist = 100.0;
    }
    if sdist == 0.0 {
        err!("**** Error: sdist cannot be 0.");
    }
    let mut smult = 2.0f64;
    if getpardouble("smult", &mut smult) == 0 {
        smult = 2.0;
    }
    if smult < 0.0 {
        err!("**** Error: smult must be non-negative.");
    }
    let mut sdadd = true;
    if sdist < 0.0 {
        sdist = 0.0 - sdist;
        sdadd = false;
    }
    let mut sdist2 = sdist;
    let mut pname_file = String::new();
    getparstring("qin", &mut pname_file);
    let fp = match fopen(&pname_file, "r") {
        Some(fp) => fp,
        None => err!("error: input Q-file did not open correctly."),
    };
    let mut pname = keyn.clone();
    if use_keyp > 0 {
        pname.push(keyp.clone());
    }
    let mut oscal = Vec::new();
    for key in &okeys {
        pname.push(key.clone());
        oscal.push(match key.as_str() {
            "gx" | "gy" | "sx" | "sy" => 1,
            "gelev" | "selev" | "sdepth" | "sdel" | "gdel" | "swdep" | "gwdep" => 2,
            _ => 0,
        });
    }
    let mut q = match getviaqfile(&fp, &pname, pname.len() as i32, 0) {
        Ok(q) => q,
        Err(e) => qerror(e),
    };
    checkpars();
    if q.ifixd == 0 {
        err!("error: input with varying number of tuples is not allowed.");
    }
    let pname = &q.names[..q.iztuple as usize];
    let mut locn = Vec::new();
    for key in &keyn {
        let loc = pname.iter().rposition(|n| n == key);
        if loc.is_none() {
            err!("error: key %s not found in non-tuple part of q-file.", key);
        }
        locn.push(loc.unwrap());
    }
    let mut oloc = Vec::new();
    for key in &okeys {
        let loc = pname.iter().rposition(|n| n == key);
        if loc.is_none() {
            err!(
                "error: input q-file must have your okeys=%s (among non-tuple names).",
                key
            );
        }
        oloc.push(loc.unwrap());
    }
    let rec = &mut q.records;
    let ncdp = rec.len();
    if use_keyp > 0 {
        let locp = match pname.iter().rposition(|n| n == &keyp) {
            Some(i) => i,
            None => err!(
                "error: input q-file must have your keyp=%s (among non-tuple names).",
                keyp
            ),
        };
        rec.sort_by(|a, b| {
            if a.dlots[locp] < b.dlots[locp] {
                std::cmp::Ordering::Less
            } else if a.dlots[locp] > b.dlots[locp] {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        for jcdp in 1..ncdp {
            if rec[jcdp - 1].dlots[locp] == rec[jcdp].dlots[locp] {
                err!(
                    "error: Two sets of values input for location = %f",
                    rec[jcdp].dlots[locp]
                );
            }
        }
    }
    let mut locxp = 0;
    let mut locyp = 0;
    let mut igikase = 0;
    let mut igckase = 0;
    if use_keyp > -2 {
        // With implicit keyp and a missing XY dimension C reads keyn[-1],
        // emitting a garbage name in this diagnostic. Use an empty name for
        // that undefined read, retaining the observed failure and exit code.
        let xkey = locx.map(|x| keyn[x].as_str()).unwrap_or("");
        let ykey = locy.map(|y| keyn[y].as_str()).unwrap_or("");
        locxp = match pname.iter().rposition(|n| n == xkey) {
            Some(i) => i,
            None => err!(
                "error: dimx %s not found in non-tuple part of q-file.",
                xkey
            ),
        };
        locyp = match pname.iter().rposition(|n| n == ykey) {
            Some(i) => i,
            None => err!(
                "error: dimy %s not found in non-tuple part of q-file.",
                ykey
            ),
        };
        igikase = GetCase("igi");
        igckase = GetCase("igc");
    }
    // C reads profsin[1]/profcos[1] outside its allocation for one point.
    // Retain a zero-filled spare slot matching the original build's allocator.
    let mut profsin = vec![0.0; ncdp.max(2)];
    let mut profcos = vec![0.0; ncdp.max(2)];
    if use_keyp > -1 {
        profsin[0] = 0.0;
        profcos[0] = 1.0;
        let mut jcdp_prev = 0;
        for jcdp in 1..ncdp {
            let dx = rec[jcdp].dlots[locxp] - rec[jcdp - 1].dlots[locxp];
            let dy = rec[jcdp].dlots[locyp] - rec[jcdp - 1].dlots[locyp];
            if dx * dx + dy * dy > 1.0e-09 {
                profsin[jcdp] = dy / (dx * dx + dy * dy).sqrt();
                profcos[jcdp] = dx / (dx * dx + dy * dy).sqrt();
                jcdp_prev = jcdp;
            } else {
                profsin[jcdp] = profsin[jcdp_prev];
                profcos[jcdp] = profcos[jcdp_prev];
            }
        }
        profsin[0] = profsin[1];
        profcos[0] = profcos[1];
    }
    let mut tree_dl = vec![vec![0.0; ncdp]; tree_numd];
    for i in 0..tree_numd {
        for j in 0..ncdp {
            tree_dl[i][j] = rec[j].dlots[locn[i]];
        }
    }
    let tree_nodes = connect_nodes(&tree_dl, ncdp, tree_numd, 1);
    if gettr(&mut tr) == 0 {
        err!("Error: cannot get first trace");
    }
    let mut near = Near::default();
    let mut tcycles = 0i32;
    let mut nproct = 0i32;
    let mut noupdate = 0i32;
    let mut ncheck = 0i32;
    loop {
        for i in 0..tree_numd {
            target[i] = match kase[i] {
                101 => tr.sx as f64,
                102 => tr.sy as f64,
                103 => tr.gx as f64,
                104 => tr.gy as f64,
                105 => 0.5 * (tr.sx as f64 + tr.gx as f64),
                106 => 0.5 * (tr.sy as f64 + tr.gy as f64),
                107 => 0.5 * (tr.grnlof as f64 + tr.gaps as f64),
                111 => tr.selev as f64,
                112 => tr.gelev as f64,
                113 => tr.sdepth as f64,
                114 => tr.sdel as f64,
                115 => tr.gdel as f64,
                116 => tr.swdep as f64,
                117 => tr.gwdep as f64,
                k => fromhead(&tr, k),
            };
            if (101..=106).contains(&kase[i]) {
                if tr.scalco > 1 {
                    target[i] *= tr.scalco as f64;
                } else if tr.scalco < 0 {
                    target[i] /= -(tr.scalco as i32) as f64;
                }
            } else if (111..=117).contains(&kase[i]) {
                if tr.scalel > 1 {
                    target[i] *= tr.scalel as f64;
                } else if tr.scalel < 0 {
                    target[i] /= -(tr.scalel as i32) as f64;
                }
            }
            if in_rel_extent[i] == 2 {
                extent_min[i] = extent_in_min[i] + target[i];
                extent_max[i] = extent_in_max[i] + target[i];
            }
        }
        if sfunc == 2 {
            if sdadd && near.num_found > 0 {
                sdist2 = sdist + near.near_dist.sqrt();
            }
            let ncycles;
            (near, ncycles) = cycle_for_near(
                &tree_nodes,
                &tree_dl,
                &tree_nt[..tree_numd],
                &extent_min[..tree_numd],
                &extent_max[..tree_numd],
                &target[..tree_numd],
                sdist2,
                smult,
            );
            tcycles = tcycles.wrapping_add(ncycles);
        } else if sfunc == 1 {
            near = find_near(
                &tree_nodes,
                &tree_dl,
                &tree_nt[..tree_numd],
                &extent_min[..tree_numd],
                &extent_max[..tree_numd],
                &target[..tree_numd],
            );
        } else {
            near = brute_near(
                &tree_dl,
                ncdp,
                &tree_nt[..tree_numd],
                &extent_min[..tree_numd],
                &extent_max[..tree_numd],
                &target[..tree_numd],
            );
        }
        if scheck {
            let c = brute_near(
                &tree_dl,
                ncdp,
                &tree_nt[..tree_numd],
                &extent_min[..tree_numd],
                &extent_max[..tree_numd],
                &target[..tree_numd],
            );
            if c.near_elem != near.near_elem
                || c.near_dist != near.near_dist
                || c.num_found != near.num_found
            {
                ncheck += 1;
                if ncheck < 10 {
                    warn!("near point: brute=%ld tree=%ld  dist diff (squared)=%g  equally near: brute=%ld tree=%ld  trace counter=%d ", c.near_elem as u64, near.near_elem as u64, c.near_dist-near.near_dist, c.num_found, near.num_found, nproct+1);
                }
            }
        }
        if near.num_found < 1 {
            if nopoint > 0 {
                err!("error: No qin point within extents for trace number %d and nopoint option is 1.", nproct+1);
            }
            noupdate += 1;
        } else {
            for i in 0..okeys.len() {
                let mut dt = rec[near.near_elem].dlots[oloc[i]];
                if oscal[i] == 1 {
                    if tr.scalco > 1 {
                        dt /= tr.scalco as f64;
                    } else if tr.scalco < 0 {
                        dt *= -(tr.scalco as i32) as f64;
                    }
                } else if oscal[i] == 2 {
                    if tr.scalel > 1 {
                        dt /= tr.scalel as f64;
                    } else if tr.scalel < 0 {
                        dt *= -(tr.scalel as i32) as f64;
                    }
                }
                tohead(&mut tr, okase[i], dt);
            }
            if use_keyp > -1 {
                let dx = target[0] - rec[near.near_elem].dlots[locxp];
                let dy = target[1] - rec[near.near_elem].dlots[locyp];
                tohead(
                    &mut tr,
                    igikase,
                    dx * profcos[near.near_elem] + dy * profsin[near.near_elem],
                );
                tohead(
                    &mut tr,
                    igckase,
                    -dx * profsin[near.near_elem] + dy * profcos[near.near_elem],
                );
            } else if use_keyp == -1 {
                let dx = target[0] - rec[near.near_elem].dlots[locxp];
                let dy = target[1] - rec[near.near_elem].dlots[locyp];
                tohead(&mut tr, igikase, 0.0);
                tohead(&mut tr, igckase, (dx * dx + dy * dy).sqrt());
            }
        }
        puttr(&mut tr);
        nproct = nproct.wrapping_add(1);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    if sfunc == 2 {
        // arm64 C integer division returns zero for a zero divisor.
        warn!(
            "Number of traces=%d  Total search cycles=%d  Cycles per trace=%d ",
            nproct,
            tcycles,
            tcycles.checked_div(nproct).unwrap_or(0)
        );
    } else {
        warn!("Number of traces=%d ", nproct);
    }
    if noupdate > 0 {
        warn!("Number of traces not updated due to no qin point within extents = %d (permitted when nopoint=0).", noupdate);
    }
    if scheck {
        warn!("There were %d total traces where brute_near disagreed with cycle_for_near or find_near.", ncheck);
    }
    exit(CWP_Exit());
}

fn qerror(errwarn: i32) -> ! {
    match errwarn {
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
        n => err!("getviaqfile error: unrecognized error code %d", n),
    }
}

#[derive(Clone, Default)]
struct Node {
    elem: usize,
    l: Option<usize>,
    r: Option<usize>,
}

fn connect_nodes(tree_dl: &[Vec<f64>], tree_numc: usize, tree_numd: usize, ihop: i32) -> Vec<Node> {
    let mut tree_nodes: Vec<Node> = (0..tree_numc)
        .map(|m| Node {
            elem: m,
            l: None,
            r: None,
        })
        .collect();
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
                    tree_nodes[n].l = Some(0); // C uses this pointer as a flag.
                    tree_nodes[m].elem = n;
                    m += 1;
                }
                n += nd + nrat;
            }
            nrat = (nrat as f64 * 0.6) as usize;
        }
        for node in &mut tree_nodes {
            node.l = None;
        }
    }
    connect_all(&mut tree_nodes, tree_dl, tree_numd);
    tree_nodes
}

fn connect_all(tree_nodes: &mut [Node], tree_dl: &[Vec<f64>], tree_numd: usize) {
    for m in 1..tree_nodes.len() {
        let mut next_node = Some(0);
        let mut now_node = 0;
        let mut naxe = 0;
        while let Some(next) = next_node {
            now_node = next;
            naxe += 1;
            if naxe == tree_numd {
                naxe = 0;
            }
            next_node =
                if tree_dl[naxe][tree_nodes[m].elem] < tree_dl[naxe][tree_nodes[now_node].elem] {
                    tree_nodes[now_node].l
                } else {
                    tree_nodes[now_node].r
                };
        }
        if tree_dl[naxe][tree_nodes[m].elem] < tree_dl[naxe][tree_nodes[now_node].elem] {
            tree_nodes[now_node].l = Some(m);
        } else {
            tree_nodes[now_node].r = Some(m);
        }
    }
}

#[derive(Default)]
struct Near {
    near_elem: usize,
    near_dist: f64,
    num_found: u64,
}

fn find_near(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_nt: &[i32],
    extent_min: &[f64],
    extent_max: &[f64],
    target: &[f64],
) -> Near {
    let mut near = Near {
        near_elem: 0,
        near_dist: f64::MAX,
        num_found: 0,
    };
    if !tree_nodes.is_empty() {
        find_near_rcur(
            tree_nodes,
            tree_dl,
            tree_nt,
            Some(0),
            1,
            extent_min,
            extent_max,
            target,
            &mut near,
        );
    }
    near
}

#[allow(clippy::too_many_arguments)]
fn find_near_rcur(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_nt: &[i32],
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
    let tree_numd = tree_nt.len();
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
        if rad <= near.near_dist {
            if rad < near.near_dist {
                near.num_found = 1;
                near.near_elem = node.elem;
                near.near_dist = rad;
            } else {
                near.num_found += 1;
                if node.elem > near.near_elem {
                    near.near_elem = node.elem;
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
    extent_min: &[f64],
    extent_max: &[f64],
    target: &[f64],
    init_rad: f64,
    rad_scal: f64,
) -> (Near, i32) {
    let tree_numd = tree_nt.len();
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
    let mut ncycles = 0i32;
    loop {
        ncycles = ncycles.wrapping_add(1);
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
        let near = find_near(
            tree_nodes,
            tree_dl,
            tree_nt,
            &loc_min[..tree_numd],
            &loc_max[..tree_numd],
            target,
        );
        if near.num_found < 1 {
            if nbig == 2 * tree_numd {
                return (near, ncycles);
            }
            now_rad *= rad_scal;
            continue;
        }
        if near.near_dist > now_rad * now_rad && nbig != 2 * tree_numd {
            now_rad = near.near_dist.sqrt() * 1.001;
            continue;
        }
        return (near, ncycles);
    }
}

fn brute_near(
    tree_dl: &[Vec<f64>],
    tree_numc: usize,
    tree_nt: &[i32],
    extent_min: &[f64],
    extent_max: &[f64],
    target: &[f64],
) -> Near {
    let mut near = Near {
        near_elem: 0,
        near_dist: f64::MAX,
        num_found: 0,
    };
    for ibrute in 0..tree_numc {
        let mut inside = true;
        for i in 0..tree_nt.len() {
            if tree_dl[i][ibrute] < extent_min[i] || tree_dl[i][ibrute] >= extent_max[i] {
                inside = false;
                break;
            }
        }
        if inside {
            let mut rad = 0.0;
            for i in 0..tree_nt.len() {
                if tree_nt[i] != 0 {
                    rad += (target[i] - tree_dl[i][ibrute]) * (target[i] - tree_dl[i][ibrute]);
                }
            }
            if rad <= near.near_dist {
                if rad < near.near_dist {
                    near.num_found = 1;
                    near.near_elem = ibrute;
                    near.near_dist = rad;
                } else {
                    near.num_found += 1;
                    if ibrute > near.near_elem {
                        near.near_elem = ibrute;
                    }
                }
            }
        }
    }
    near
}
