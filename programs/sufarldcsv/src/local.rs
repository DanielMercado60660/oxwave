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

//! C-local helpers from su/main/statics/sufarldcsv.c (Revision 1.01).
//! Credits: Created Oct 2023: Andre Latour (Andre Roger Latour).

use std::cmp::Ordering;

/// C pointers to nodes become optional indices in the same allocation.
#[derive(Clone, Debug, Default)]
pub(super) struct Node {
    pub elem: usize,
    pub l: Option<usize>,
    pub r: Option<usize>,
}

pub(super) fn connect_nodes(
    tree_nodes: &mut [Node],
    tree_numc: usize,
    tree_dl: &[Vec<f64>],
    tree_numd: usize,
    ihop: i32,
) {
    for m in 0..tree_numc {
        tree_nodes[m].l = None;
        tree_nodes[m].r = None;
        tree_nodes[m].elem = m;
    }

    // Hop through the input before connecting the nodes. The l pointer is
    // temporarily a flag, exactly as in the C implementation.
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
                    tree_nodes[n].l = Some(0);
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
    connect_all(tree_nodes, tree_numc, tree_dl, tree_numd);
}

pub(super) fn connect_all(
    tree_nodes: &mut [Node],
    tree_numc: usize,
    tree_dl: &[Vec<f64>],
    tree_numd: usize,
) {
    for m in 1..tree_numc {
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

pub(super) fn find_in(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_numd: usize,
    extent_min: &[f64],
    extent_max: &[f64],
    out_elem: &mut [usize],
    num_out: &mut usize,
) {
    *num_out = 0;
    find_in_rcur(
        tree_nodes,
        tree_dl,
        tree_numd,
        Some(0),
        1,
        extent_min,
        extent_max,
        out_elem,
        num_out,
    );
}

#[allow(clippy::too_many_arguments)]
fn find_in_rcur(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_numd: usize,
    now_node: Option<usize>,
    mut naxe: usize,
    extent_min: &[f64],
    extent_max: &[f64],
    out_elem: &mut [usize],
    num_out: &mut usize,
) {
    let Some(now) = now_node else { return };
    let now_node = &tree_nodes[now];
    let mut inside = true;
    for i in 0..tree_numd {
        if tree_dl[i][now_node.elem] < extent_min[i] || tree_dl[i][now_node.elem] >= extent_max[i] {
            inside = false;
            break;
        }
    }
    if inside {
        out_elem[*num_out] = now_node.elem;
        *num_out += 1;
    }

    if naxe == tree_numd {
        naxe = 0;
    }
    if tree_dl[naxe][now_node.elem] >= extent_min[naxe] && now_node.l.is_some() {
        find_in_rcur(
            tree_nodes,
            tree_dl,
            tree_numd,
            now_node.l,
            naxe + 1,
            extent_min,
            extent_max,
            out_elem,
            num_out,
        );
    }
    if tree_dl[naxe][now_node.elem] < extent_max[naxe] && now_node.r.is_some() {
        find_in_rcur(
            tree_nodes,
            tree_dl,
            tree_numd,
            now_node.r,
            naxe + 1,
            extent_min,
            extent_max,
            out_elem,
            num_out,
        );
    }
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
    near_elem: &mut usize,
    near_dist: &mut f64,
    num_found: &mut usize,
) {
    let Some(now) = now_node else { return };
    let now_node = &tree_nodes[now];
    let mut inside = true;
    for i in 0..tree_numd {
        if tree_dl[i][now_node.elem] < extent_min[i] || tree_dl[i][now_node.elem] >= extent_max[i] {
            inside = false;
            break;
        }
    }
    if inside {
        let mut rad = 0.0;
        for i in 0..tree_numd {
            if tree_nt[i] != 0 {
                rad += (target[i] - tree_dl[i][now_node.elem])
                    * (target[i] - tree_dl[i][now_node.elem]);
            }
        }
        if rad <= *near_dist {
            if rad < *near_dist {
                *num_found = 1;
                *near_elem = now_node.elem;
                *near_dist = rad;
            } else {
                *num_found += 1;
                if now_node.elem > *near_elem {
                    *near_elem = now_node.elem;
                }
            }
        }
    }

    if naxe == tree_numd {
        naxe = 0;
    }
    if tree_dl[naxe][now_node.elem] >= extent_min[naxe] && now_node.l.is_some() {
        find_near_rcur(
            tree_nodes,
            tree_dl,
            tree_nt,
            tree_numd,
            now_node.l,
            naxe + 1,
            extent_min,
            extent_max,
            target,
            near_elem,
            near_dist,
            num_found,
        );
    }
    if tree_dl[naxe][now_node.elem] < extent_max[naxe] && now_node.r.is_some() {
        find_near_rcur(
            tree_nodes,
            tree_dl,
            tree_nt,
            tree_numd,
            now_node.r,
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

#[allow(clippy::too_many_arguments)]
pub(super) fn cycle_for_near(
    tree_nodes: &[Node],
    tree_dl: &[Vec<f64>],
    tree_nt: &[i32],
    tree_numd: usize,
    extent_min: &[f64],
    extent_max: &[f64],
    target: &[f64],
    init_rad: f64,
    rad_scal: f64,
    near_elem: &mut usize,
    near_dist: &mut f64,
    num_found: &mut usize,
    ncycles: &mut i32,
) {
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
    *ncycles = 0;
    loop {
        *ncycles = ncycles.wrapping_add(1);
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
        *num_found = 0;
        *near_dist = f64::MAX;
        *near_elem = 0;
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
            near_elem,
            near_dist,
            num_found,
        );
        if *num_found < 1 {
            if nbig == 2 * tree_numd {
                return;
            }
            now_rad *= rad_scal;
            continue;
        }
        // Searching a square may miss a closer point outside that square.
        if *near_dist > now_rad * now_rad && nbig != 2 * tree_numd {
            now_rad = near_dist.sqrt() * 1.001;
            continue;
        }
        return;
    }
}

/// The original comparator's global locp becomes an explicit argument.
pub(super) fn comp_sort1(q1: &[f64], q2: &[f64], locp: usize) -> Ordering {
    if q1[locp] < q2[locp] {
        return Ordering::Less;
    }
    if q1[locp] > q2[locp] {
        return Ordering::Greater;
    }
    Ordering::Equal
}

/// Running average; negative nmax applies the original smoothing step.
pub(super) fn runav(
    dall: &mut [Vec<f64>],
    ncdp: i32,
    klast: i32,
    nmin: i32,
    mut nmax: i32,
    dtog: &mut [f64],
    ierr: &mut i32,
) {
    let mut iavr = 1;
    if nmax < 0 {
        nmax = nmax.wrapping_neg();
        iavr = 0;
    }
    if ncdp < 1 || klast < 1 || nmin < 0 || nmax < nmin {
        *ierr = 2;
        return;
    }
    if 2i32.wrapping_mul(nmax.wrapping_add(3)) > ncdp {
        *ierr = 1;
        return;
    }
    *ierr = 0;
    let mut n = 0;
    for k in 0..klast as usize {
        dtog[k] = 0.0;
    }
    for i in 0..(nmax * 2 + 1) as usize {
        n += 1;
        for k in 0..klast as usize {
            dtog[k] += dall[i][k];
        }
    }
    for i in nmax - nmin + 1..nmax + nmin {
        n -= 1;
        for k in 0..klast as usize {
            dtog[k] -= dall[i as usize][k];
        }
    }
    for i in nmax..ncdp - nmax - 1 {
        for k in 0..klast as usize {
            dtog[k] += dall[(i + 1 + nmax) as usize][k] - dall[(i - nmax) as usize][k];
        }
        if nmin > 0 {
            for k in 0..klast as usize {
                dtog[k] += dall[(i + 1 - nmin) as usize][k] - dall[(i + nmin) as usize][k];
            }
        }
        if iavr == 1 {
            for k in 0..klast as usize {
                dall[(i - nmax) as usize][k] = dtog[k] / n as f64;
            }
        } else {
            for k in 0..klast as usize {
                dall[(i - nmax) as usize][k] =
                    dall[(i + 1) as usize][k] + dall[(i + 1) as usize][k] - dtog[k] / n as f64;
            }
        }
    }
    for i in (nmax..=ncdp - nmax - 1).rev() {
        for k in 0..klast as usize {
            dall[i as usize][k] = dall[(i - nmax) as usize][k];
        }
    }
}

/// Extrapolate the ends, average (or smooth twice), and unshift the output.
pub(super) fn runsmo(
    dalla: &mut [Vec<f64>],
    msize: i32,
    klast: i32,
    nback: i32,
    nmaxi: i32,
    dtog: &mut [f64],
    ierr: &mut i32,
) {
    let nmax = if nmaxi < 0 {
        nmaxi.wrapping_neg()
    } else {
        nmaxi
    };
    if msize < 1 || klast < 1 || nmax < 1 || nback < 0 {
        *ierr = 2;
        return;
    }
    if 2i32.wrapping_mul(nmax).wrapping_add(6) > msize {
        *ierr = 1;
        return;
    }
    let mut dback = 0.0;
    if nback > 0 {
        dback = 1.0 / nback as f64;
    }
    *ierr = 0;
    let lastcdp = msize - 2 * nmax - 1;
    for j in (0..lastcdp).rev() {
        for i in 0..klast as usize {
            dalla[(j + 1 + nmax) as usize][i] = dalla[j as usize][i];
        }
    }
    for j in 0..nmax + 1 {
        for i in 0..klast as usize {
            let dt =
                (dalla[(nmax + 1) as usize][i] - dalla[(nmax + 1 + nback) as usize][i]) * dback;
            dalla[j as usize][i] = dalla[(nmax + 1) as usize][i] + (nmax - j + 1) as f64 * dt;
        }
    }
    for j in lastcdp + nmax + 1..lastcdp + 2 * nmax + 1 {
        for i in 0..klast as usize {
            let dt = (dalla[(lastcdp + nmax) as usize][i]
                - dalla[(lastcdp + nmax - nback) as usize][i])
                * dback;
            dalla[j as usize][i] =
                dalla[(lastcdp + nmax) as usize][i] + (j - lastcdp - nmax) as f64 * dt;
        }
    }
    runav(dalla, lastcdp + 2 * nmax + 1, klast, 0, nmaxi, dtog, ierr);
    if *ierr > 0 {
        return;
    }
    for j in 0..lastcdp {
        for i in 0..klast as usize {
            dalla[j as usize][i] = dalla[(j + nmax) as usize][i];
        }
    }
    if nmaxi > 0 {
        return;
    }

    // The second pass uses the opposite sign of nmaxi.
    for j in (0..lastcdp).rev() {
        for i in 0..klast as usize {
            dalla[(j + 1 + nmax) as usize][i] = dalla[j as usize][i];
        }
    }
    for j in 0..nmax + 1 {
        for i in 0..klast as usize {
            let dt =
                (dalla[(nmax + 1) as usize][i] - dalla[(nmax + 1 + nback) as usize][i]) * dback;
            dalla[j as usize][i] = dalla[(nmax + 1) as usize][i] + (nmax - j + 1) as f64 * dt;
        }
    }
    for j in lastcdp + nmax + 1..lastcdp + 2 * nmax + 1 {
        for i in 0..klast as usize {
            let dt = (dalla[(lastcdp + nmax) as usize][i]
                - dalla[(lastcdp + nmax - nback) as usize][i])
                * dback;
            dalla[j as usize][i] =
                dalla[(lastcdp + nmax) as usize][i] + (j - lastcdp - nmax) as f64 * dt;
        }
    }
    runav(
        dalla,
        lastcdp + 2 * nmax + 1,
        klast,
        0,
        nmaxi.wrapping_neg(),
        dtog,
        ierr,
    );
    if *ierr > 0 {
        return;
    }
    for j in 0..lastcdp {
        for i in 0..klast as usize {
            dalla[j as usize][i] = dalla[(j + nmax) as usize][i];
        }
    }
}
