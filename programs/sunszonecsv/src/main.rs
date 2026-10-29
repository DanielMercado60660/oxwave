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

//! SUNSZONECSV - Near Surface Zone Extract and Shift.
//! Port of su/main/statics/sunszonecsv.c (Revision 1.01).
//! Credits: June 2024, Andre Latour; started from sunearcsv.

use oxwave_core::prelude::*;
use oxwave_q_file::*;
use oxwave_sinc_interp::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUNSZONECSV - Near Surface Zone Extract and Shift                          ",
    "\t\t\t\t\t\t\t\t\t     ",
    "  sunszonecsv [parameters].                                                 ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Parameters:\t         \t\t\t\t\t\t     ",
    "                                                                            ",
    " qin=      Input file containing q-records. This parameter is required.     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qx=sx         Name of X coordinate at points in the q-file.                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qy=sy         Name of Y coordinate at points in the q-file.                ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qvw=vw        Name of top layer velocity at points in the q-file.          ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qzw=zw_true   Name of top layer thickness at points in the q-file.         ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qvb=vb        Name of second layer velocity at points in the q-file.       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qzb=zb_true   Name of second layer thickness at points in the q-file.      ",
    "    =unknown   Second layer has unknown thickness.                          ",
    "   ***   Note: Once unknown is specified (or defaults) for a thickness,     ",
    "               parameters for deeper layers are not used.                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qvc=vc        Name of third layer velocity at points in the q-file.        ",
    "\t\t\t\t\t\t\t\t\t     ",
    " qzc=zc_true   Name of third layer thickness at points in the q-file.       ",
    "    =unknown   Third layer has unknown thickness.                           ",
    "                                                                            ",
    " qvd=vd        Name of fourth layer velocity at points in the q-file.       ",
    "\t\t\t\t\t\t\t\t\t     ",
    " locs=0        Use near surface model from shot and receiver ends of traces.",
    "     =1        Use near surface model from shot end only.                   ",
    "     =-1       Use near surface model from receiver end only.               ",
    "         Note: The model is assumed to be flat at both ends, and raypaths   ",
    "               are assumed to be straight with a linear change in layer     ",
    "               velocities from shot to receiver. The distance between the   ",
    "               shot and receiver ends of the raypaths are adjusted to take  ",
    "               differences in elevation and layer thicknesses into account  ",
    "               and the shot hole depth is used at the shot end of the paths.",
    "               The shot HOLE DEPTHS and the ELEVATION DIFFERENCES ARE USED  ",
    "               REGARDLESS of which option is chozen here. Option 1 simply   ",
    "               uses the velocities and thicknesses of the shot location     ",
    "               for both ends of the raypaths. Option -1 simply uses the     ",
    "               velocities and thicknesses of the receiver location for      ",
    "               both ends of the raypaths.                                   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " tpath=0       Travel path shift value.                                     ",
    "               The default (0) is to set the shift to the minimum time of   ",
    "               the raypaths computed for all layers and the direct arrival  ",
    "               path (also known as the first arrival or first break time).  ",
    "      =1       Set shift to top-of-first layer direct arrival time.         ",
    "      =2       Set shift to top-of-second layer refraction arrival time.    ",
    "      =3       Set shift to top-of-third layer refraction arrival time.     ",
    "      =4       Set shift to top-of-fourth layer refraction arrival time.    ",
    "\t\t\t\t\t\t\t\t\t     ",
    " tapp=1        Apply the total time shift determined by the tpath option    ",
    "               and also set the sstat,gstat,tstat,sec key values.           ",
    "               The sstat is set to the time computed for the shot end of    ",
    "               the path, the gstat is set to the time computed for the      ",
    "               receiver end of the path, and tstat is set to time computed  ",
    "               for the adjusted travelpath between shot and receiver.       ",
    "               If tpath=0 the sec key is set to the layer number that       ",
    "               produced the minimum travel path shift value, otherwise it   ",
    "               is set to the same layer number as the tpath parameter.      ",
    "      =0       Apply the time shift determined by the tpath option          ",
    "               but do not alter the sstat,gstat,tstat,sec key values.       ",
    "      =-1      Do not shift. Only set the stat,gstat,tstat,sec key values.  ",
    "         Note: The shift is done by the same routine that does shifts for   ",
    "               statics and therefor applies fractional millisecond values   ",
    "               BUT sstat,gstat,tstat keys are set to nearest millisecond.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " tbulk=200.0   Shift output traces by this additional time in milliseconds  ",
    "               (except if tapp=-1). This allows seeing more of the data     ",
    "               before the computed refraction wavelet times. This has no    ",
    "               effect on computations or sstat,gstat,tstat,sec values.      ",
    "\t\t\t\t\t\t\t\t\t     ",
    " apply=-1      Shift Positive trace times towards zero time.                ",
    "        1      Shift Positive trace times away from zero time.              ",
    "         Note: The concept here is that this parameter allows you to        ",
    "               flatten refraction arrivals, apply some data processing,     ",
    "               then inverse the flattening.                                 ",
    "         Note: In both cases the shift is computed from model (and tbulk),  ",
    "               the input sstat,gstat,tstat,sec key values are not used.     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ***  The trace header key values sx,sy,selev,sdepth,gx,gy,gelev are      ",
    "        always accessed from the trace headers. The nearest velocities      ",
    "        and thicknesses are accessed from the qin file. Elevations and      ",
    "        shot hole depths are only accessed from traces, not from qin file.  ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ***  Input trace header values are scaled using scalco and scalel keys.  ",
    "        But q-file standards means qin values should have their true        ",
    "        magnitudes and therefor do not get scaled by this program.          ",
    "\t\t\t\t\t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " sdist=100 Initial search distance. If positive, add this value to the      ",
    "           distance between the previous trace and its nearest qin point    ",
    "           and used as the initial search distance for the current trace.   ",
    "           If negative, the initial search distance for all traces is set   ",
    "           to this absolute value.                                          ",
    " smult=2   Search multiplier.                                               ",
    "           If the nearest qin point is not found after searching with the   ",
    "           initial search distance, the search distance is multiplied by    ",
    "           this value, and search is performed again. This repeats until    ",
    "           finding the nearest point (or all min,max ranges are exceeded).  ",
    "\t\t\t\t\t\t\t\t\t     ",
    "   ------------------------------------------------------------------       ",
    "   ------------------------------------------------------------------       ",
    "\t\t\t\t\t\t\t\t\t     ",
];

fn read_model(fp: &CFile, pname: &[String]) -> QFile {
    match getviaqfile(fp, pname, pname.len() as i32, 0) {
        Ok(model) => model,
        Err(errwarn) => match errwarn {
            1 => err!("getqinfo error: extra C_SU_NAMES record in q-file"),
            2 => err!("getqinfo error: extra C_SU_NDIMS record in q-file"),
            3 => err!("getqinfo error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            11 => err!("readqhead error: if C_SU_NDIMS not vary, its numbers must align with C_SU_NAMES"),
            12 => err!("readqhead error: C_SU_ID record not found immediately after C_SU_NAMES record."),
            22 => err!("getviaqfile error: C_SU_NDIMS record not same length as C_SU_NAMES record."),
            23 => err!("getviaqfile error: C_SU_NAMES tupled names out-of-order, or accidental duplicate name"),
            24 => err!("getviaqfile error: C_SU_NDIMS blank where valid number expected"),
            25 => err!("getviaqfile error: C_SU_NDIMS non-number where valid number expected"),
            26 => err!("getviaqfile error: C_SU_NDIMS value must be same for all members of tuple"),
            27 => err!("getviaqfile error: C_SU_NAMES record followed by C_SU_ID record not found."),
            n if n > 100 => err!("getviaqfile error: record %d (wrong comma count, damaged, non-numbers, ...)", n - 99),
            n => err!("getviaqfile error: unrecognized error code %d", n),
        },
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut locs = 0;
    getparint("locs", &mut locs);
    if locs < -1 || locs > 1 {
        err!("**** Error: locs parameter out-of-range.");
    }
    let mut iapply = -1;
    getparint("apply", &mut iapply);
    if iapply != -1 && iapply != 1 {
        err!("**** Error: apply parameter out-of-range.");
    }
    let mut sdist = 100.0;
    getpardouble("sdist", &mut sdist);
    if sdist == 0.0 {
        err!("**** Error: sdist cannot be 0.");
    }
    let mut smult = 2.0;
    getpardouble("smult", &mut smult);
    if smult < 0.0 {
        err!("**** Error: smult must be non-negative.");
    }
    let mut sdadd = 1;
    if sdist < 0.0 {
        sdist = 0.0 - sdist;
        sdadd = 0;
    }
    let mut sdist2 = sdist;
    let gmult = smult;
    let gdadd = sdadd;
    let gdist = sdist;
    let mut gdist2 = sdist2;

    let mut pname_file = String::new();
    getparstring("qin", &mut pname_file);
    let fp = match fopen(&pname_file, "r") {
        Some(fp) => fp,
        None => err!("error: input Q-file did not open correctly."),
    };
    let mut pname: Vec<String> = Vec::new();
    for (key, default) in [
        ("qx", "sx"),
        ("qy", "sy"),
        ("qvw", "vw"),
        ("qzw", "zw_true"),
        ("qvb", "vb"),
        ("qzb", "zb_true"),
        ("qvc", "vc"),
        ("qzc", "zc_true"),
        ("qvd", "vd"),
    ] {
        let mut name = default.to_string();
        if countparval(key) > 0 {
            getparstring(key, &mut name);
        }
        pname.push(name);
    }
    // nlayer is the last FULL layer; one more layer has a velocity only.
    let mut nlayer = 3;
    let mut numpname = 9;
    if pname[7] == "unknown" {
        nlayer = 2;
        numpname = 7;
    }
    if pname[5] == "unknown" {
        nlayer = 1;
        numpname = 5;
    }
    let mut tpath = 0;
    getparint("tpath", &mut tpath);
    if tpath < 0 || tpath > nlayer as i32 + 1 {
        err!("**** Error: tpath parameter out-of-range.");
    }
    let mut tapp = 1;
    getparint("tapp", &mut tapp);
    if tapp < -1 || tapp > 1 {
        err!("**** Error: tapp parameter out-of-range.");
    }
    let mut tbulk = 200.0;
    getpardouble("tbulk", &mut tbulk);
    let iname = &pname[..numpname];
    let model = read_model(&fp, iname);
    // C's tree search assumes at least one Q record; an empty model is undefined.
    if model.records.is_empty() {
        err!("error: input Q-file contains no Q records.");
    }

    checkpars();
    if model.ifixd == 0 {
        err!("error: input with varying number of tuples is not allowed.");
    }
    let mut locn = [0; 9];
    for i in 0..numpname {
        let mut loc = None;
        for j in 0..model.iztuple as usize {
            if model.names[j] == iname[i] {
                loc = Some(j);
            }
        }
        locn[i] = match loc {
            Some(j) => j,
            None => err!(
                "qin file error: name %s not found in non-tuple part of qin Q-file.",
                iname[i]
            ),
        };
    }
    let mut tree_dl = [
        vec![0.0; model.records.len()],
        vec![0.0; model.records.len()],
    ];
    for j in 0..model.records.len() {
        tree_dl[0][j] = model.records[j].dlots[locn[0]];
    }
    for j in 0..model.records.len() {
        tree_dl[1][j] = model.records[j].dlots[locn[1]];
    }
    let tree_nodes = connect_nodes(&tree_dl, 1);
    let mut tr = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Error: cannot get first trace");
    }
    let lentrc = tr.ns as usize;
    let mut onetrc = vec![0.0f32; lentrc];
    let mut oneshift = vec![0.0f32; lentrc];
    let sampi = (tr.dt as f32 as f64 / 1000.0) as f32;
    let mut near_dist_s: f64 = 0.0;
    let mut near_dist_g: f64 = 0.0;
    let mut vls = [0.0f64; 4];
    let mut vlg = [0.0f64; 4];
    let mut zps = [0.0f64; 4];
    let mut zpg = [0.0f64; 4];
    let mut eps = [0.0f64; 4];
    let mut epg = [0.0f64; 4];
    let mut ptm = [0.0f64; 4];
    let mut stm = [0.0f64; 4];
    let mut gtm = [0.0f64; 4];
    let mut nproct: i32 = 0;
    loop {
        let mut tars = [tr.sx as f64, tr.sy as f64];
        eps[0] = tr.selev as f64;
        let mut shole = tr.sdepth as f64;
        let mut targ = [tr.gx as f64, tr.gy as f64];
        epg[0] = tr.gelev as f64;
        if tr.scalco > 1 {
            tars[0] *= tr.scalco as f64;
            tars[1] *= tr.scalco as f64;
            targ[0] *= tr.scalco as f64;
            targ[1] *= tr.scalco as f64;
        } else if tr.scalco < 0 {
            tars[0] /= -(tr.scalco as i32) as f64;
            tars[1] /= -(tr.scalco as i32) as f64;
            targ[0] /= -(tr.scalco as i32) as f64;
            targ[1] /= -(tr.scalco as i32) as f64;
        }
        if tr.scalel > 1 {
            eps[0] *= tr.scalel as f64;
            shole *= tr.scalel as f64;
            epg[0] *= tr.scalel as f64;
        } else if tr.scalel < 0 {
            eps[0] /= -(tr.scalel as i32) as f64;
            shole /= -(tr.scalel as i32) as f64;
            epg[0] /= -(tr.scalel as i32) as f64;
        }
        if locs >= 0 {
            if sdadd == 1 {
                sdist2 = sdist + near_dist_s.sqrt();
            }
            let (near_s, dist, _, _) = cycle_for_near(&tree_nodes, &tree_dl, &tars, sdist2, smult);
            near_dist_s = dist;
            for i in 0..nlayer + 1 {
                vls[i] = model.records[near_s].dlots[locn[2 + i * 2]];
            }
            for i in 0..nlayer {
                zps[i] = model.records[near_s].dlots[locn[3 + i * 2]];
            }
            for i in 1..nlayer + 1 {
                eps[i] = eps[i - 1] - zps[i - 1];
            }
        }
        if locs <= 0 {
            if gdadd == 1 {
                gdist2 = gdist + near_dist_g.sqrt();
            }
            let (near_g, dist, _, _) = cycle_for_near(&tree_nodes, &tree_dl, &targ, gdist2, gmult);
            near_dist_g = dist;
            for i in 0..nlayer + 1 {
                vlg[i] = model.records[near_g].dlots[locn[2 + i * 2]];
            }
            for i in 0..nlayer {
                zpg[i] = model.records[near_g].dlots[locn[3 + i * 2]];
            }
            for i in 1..nlayer + 1 {
                epg[i] = epg[i - 1] - zpg[i - 1];
            }
        }
        if locs == 1 {
            for i in 0..nlayer + 1 {
                vlg[i] = vls[i];
            }
            for i in 0..nlayer {
                zpg[i] = zps[i];
            }
            for i in 1..nlayer + 1 {
                epg[i] = epg[i - 1] - zps[i - 1];
            }
        }
        if locs == -1 {
            for i in 0..nlayer + 1 {
                vls[i] = vlg[i];
            }
            for i in 0..nlayer {
                zps[i] = zpg[i];
            }
            for i in 1..nlayer + 1 {
                eps[i] = eps[i - 1] - zpg[i - 1];
            }
        }
        let poff =
            (tars[0] - targ[0]) * (tars[0] - targ[0]) + (tars[1] - targ[1]) * (tars[1] - targ[1]);
        for i in 0..nlayer + 1 {
            ptm[i] =
                (poff + (eps[i] - epg[i]) * (eps[i] - epg[i])).sqrt() / ((vls[i] + vlg[i]) / 2.0);
        }
        if shole > 0.0 {
            for i in 0..nlayer {
                if shole >= zps[i] {
                    shole -= zps[i];
                    zps[i] = 0.0;
                } else {
                    zps[i] -= shole;
                    break;
                }
            }
        }
        for i in 0..nlayer + 1 {
            stm[i] = 0.0;
            gtm[i] = 0.0;
            for j in 0..i {
                stm[i] += (zps[j] * (vls[i] * vls[i] - vls[j] * vls[j]).sqrt()) / vls[j];
                gtm[i] += (zpg[j] * (vlg[i] * vlg[i] - vlg[j] * vlg[j]).sqrt()) / vlg[j];
            }
            stm[i] /= vls[i];
            gtm[i] /= vlg[i];
        }
        let (mut tshift, mut ilayer);
        if tpath == 0 {
            tshift = stm[0] + ptm[0] + gtm[0];
            ilayer = 0;
            for i in 1..nlayer + 1 {
                if stm[i] + ptm[i] + gtm[i] < tshift {
                    tshift = stm[i] + ptm[i] + gtm[i];
                    ilayer = i;
                }
            }
        } else {
            ilayer = tpath as usize - 1;
            tshift = stm[ilayer] + ptm[ilayer] + gtm[ilayer];
        }
        if tapp > -1 {
            for i in 0..lentrc {
                onetrc[i] = tr.data[i];
                oneshift[i] = (i as f32 as f64
                    - iapply as f64 * (1000.0 * tshift - tbulk) / sampi as f64)
                    as f32;
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
        }
        if tapp != 0 {
            // C converts double through int before narrowing to signed short.
            tr.sstat = (stm[ilayer] * 1000.0 + 0.5000000001) as i32 as i16;
            tr.gstat = (gtm[ilayer] * 1000.0 + 0.5000000001) as i32 as i16;
            tr.tstat = (ptm[ilayer] * 1000.0 + 0.5000000001) as i32 as i16;
            tr.sec = (ilayer + 1) as i16;
        }
        puttr(&mut tr);
        nproct = nproct.wrapping_add(1);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    warn!("Number of traces=%d ", nproct);
    exit(CWP_Exit());
}

#[derive(Clone, Default)]
struct Node {
    elem: usize,
    l: Option<usize>,
    r: Option<usize>,
}

// These four tree functions are local to sunszonecsv.c, not library routines.
fn connect_nodes(tree_dl: &[Vec<f64>; 2], ihop: i32) -> Vec<Node> {
    let tree_numc = tree_dl[0].len();
    let mut nodes = vec![Node::default(); tree_numc];
    for m in 0..tree_numc {
        nodes[m].elem = m;
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
                if nodes[n].l.is_none() {
                    nodes[n].l = Some(0);
                    nodes[m].elem = n;
                    m += 1;
                }
                n += nd + nrat;
            }
            nrat = (nrat as f64 * 0.6) as usize;
        }
        for m in 0..tree_numc {
            nodes[m].l = None;
        }
    }
    connect_all(&mut nodes, tree_dl);
    nodes
}

fn connect_all(nodes: &mut [Node], tree_dl: &[Vec<f64>; 2]) {
    for m in 1..nodes.len() {
        let mut next_node = Some(0);
        let mut now_node = 0;
        let mut naxe = 0;
        while let Some(next) = next_node {
            now_node = next;
            naxe += 1;
            if naxe == 2 {
                naxe = 0;
            }
            next_node = if tree_dl[naxe][nodes[m].elem] < tree_dl[naxe][nodes[now_node].elem] {
                nodes[now_node].l
            } else {
                nodes[now_node].r
            };
        }
        if tree_dl[naxe][nodes[m].elem] < tree_dl[naxe][nodes[now_node].elem] {
            nodes[now_node].l = Some(m);
        } else {
            nodes[now_node].r = Some(m);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn find_near_rcur(
    nodes: &[Node],
    tree_dl: &[Vec<f64>; 2],
    now: Option<usize>,
    mut naxe: usize,
    extent_min: &[f64; 2],
    extent_max: &[f64; 2],
    target: &[f64; 2],
    near_elem: &mut usize,
    near_dist: &mut f64,
    num_found: &mut u64,
) {
    let Some(now) = now else {
        return;
    };
    let node = &nodes[now];
    let mut inside = true;
    for i in 0..2 {
        if tree_dl[i][node.elem] < extent_min[i] || tree_dl[i][node.elem] >= extent_max[i] {
            inside = false;
            break;
        }
    }
    if inside {
        let mut rad = 0.0;
        for i in 0..2 {
            rad += (target[i] - tree_dl[i][node.elem]) * (target[i] - tree_dl[i][node.elem]);
        }
        if rad <= *near_dist {
            if rad < *near_dist {
                *num_found = 1;
                *near_elem = node.elem;
                *near_dist = rad;
            } else {
                *num_found += 1;
                if node.elem > *near_elem {
                    *near_elem = node.elem;
                }
            }
        }
    }
    if naxe == 2 {
        naxe = 0;
    }
    if tree_dl[naxe][node.elem] >= extent_min[naxe] && node.l.is_some() {
        find_near_rcur(
            nodes,
            tree_dl,
            node.l,
            naxe + 1,
            extent_min,
            extent_max,
            target,
            near_elem,
            near_dist,
            num_found,
        );
    }
    if tree_dl[naxe][node.elem] < extent_max[naxe] && node.r.is_some() {
        find_near_rcur(
            nodes,
            tree_dl,
            node.r,
            naxe + 1,
            extent_min,
            extent_max,
            target,
            near_elem,
            near_dist,
            num_found,
        );
    }
}

fn cycle_for_near(
    nodes: &[Node],
    tree_dl: &[Vec<f64>; 2],
    target: &[f64; 2],
    init_rad: f64,
    rad_scal: f64,
) -> (usize, f64, u64, i32) {
    let extent_min = [-f64::MAX; 2];
    let extent_max = [f64::MAX; 2];
    let mut now_rad = init_rad;
    let mut ncycles = 0;
    loop {
        ncycles += 1;
        let mut loc_min = [0.0; 2];
        let mut loc_max = [0.0; 2];
        let mut nbig = 0;
        for i in 0..2 {
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
        let mut num_found = 0;
        let mut near_dist = f64::MAX;
        let mut near_elem = 0;
        find_near_rcur(
            nodes,
            tree_dl,
            Some(0),
            1,
            &loc_min,
            &loc_max,
            target,
            &mut near_elem,
            &mut near_dist,
            &mut num_found,
        );
        if num_found < 1 {
            if nbig == 4 {
                return (near_elem, near_dist, num_found, ncycles);
            }
            now_rad *= rad_scal;
            continue;
        }
        if near_dist > now_rad * now_rad && nbig != 4 {
            now_rad = near_dist.sqrt() * 1.001;
            continue;
        }
        return (near_elem, near_dist, num_found, ncycles);
    }
}
