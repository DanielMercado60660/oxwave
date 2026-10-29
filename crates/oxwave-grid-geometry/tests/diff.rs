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

//! Differential tests of `su/lib/gridxy.c` against the original SU library.
//! Credits: no individual author is named in gridxy.c; gridxy.h credits
//! "Copyright (c) Colorado School of Mines, 2021. All rights reserved."

use oxwave_cref::{CArgv, XorShift};
use oxwave_grid_geometry as rust;

extern "C" {
    fn gridset(gvals: *mut f64, errwarn: *mut i32);
    fn gridrawxycdpic(
        gvals: *mut f64,
        dx: f64,
        dy: f64,
        icdp: *mut i32,
        igi: *mut i32,
        igc: *mut i32,
    );
    fn gridicrawxy(gvals: *mut f64, igi: i32, igc: i32, dx: *mut f64, dy: *mut f64);
    fn gridicgridxy(gvals: *mut f64, igi: i32, igc: i32, dx: *mut f64, dy: *mut f64);
    fn gridiccdp(gvals: *mut f64, igi: i32, igc: i32, icdp: *mut i32);
    fn gridiccdp90(gvals: *mut f64, igi: i32, igc: i32, jcdp: *mut i32);
    fn gridcdpic(gvals: *mut f64, icdp: i32, igi: *mut i32, igc: *mut i32);
    fn gridcdpic90(gvals: *mut f64, jcdp: i32, igi: *mut i32, igc: *mut i32);
    fn gridrawxygridxy(gvals: *mut f64, dx: f64, dy: f64, tx: *mut f64, ty: *mut f64);
    fn gridgridxyrawxy(gvals: *mut f64, dx: f64, dy: f64, tx: *mut f64, ty: *mut f64);
    fn gridcheck(gvals: *mut f64, icheck: i32, errwarn: *mut i32);
}

const TRIALS: usize = 4096;
const SIZES: [i32; 12] = [1, 2, 3, 7, 17, 31, 64, 127, 256, 1024, 8191, 32768];
type Grid = [f64; 20]; // Two trailing sentinels detect unintended writes.

// All floating results on this build are bit-identical, stronger than 1e-5.
// NaN payloads are not specified by C; compare NaNs by classification.
fn same(a: f64, b: f64) {
    assert!(
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
        "Rust={a:?} ({:016x}), C={b:?} ({:016x})",
        a.to_bits(),
        b.to_bits()
    );
}

fn same_grid(a: &Grid, b: &Grid) {
    for i in 0..a.len() {
        same(a[i], b[i]);
    }
}

fn raw_grid(r: &mut XorShift, trial: usize) -> Grid {
    let mut g = [0.0; 20];
    for x in &mut g {
        *x = r.uniform(-1.0e6, 1.0e6);
    }
    let scale = 10.0_f64.powf(r.range_i32(-9, 6) as f64);
    g[2] = r.uniform(-1.0e4, 1.0e4) * scale;
    g[3] = r.uniform(-1.0e4, 1.0e4) * scale;
    if trial % 9 == 0 {
        g[2] = 0.0;
        g[3] = -0.0;
    }
    g[10] = r.uniform(0.25, 100.0) * scale;
    g[11] = r.uniform(0.25, 100.0) * scale;
    let nb = SIZES[trial % SIZES.len()].max(2);
    let nc = SIZES[(trial / SIZES.len()) % SIZES.len()];
    // Exact axis directions as well as arbitrary orientations in every quadrant.
    let (sin, cos) = match trial % 8 {
        0 => (0.0, 1.0),
        1 => (1.0, 0.0),
        2 => (0.0, -1.0),
        3 => (-1.0, 0.0),
        _ => r
            .uniform(-std::f64::consts::PI, std::f64::consts::PI)
            .sin_cos(),
    };
    let b = (nb - 1) as f64 * g[10] + r.uniform(-0.4, 0.4) * g[10];
    let c = (nc - 1) as f64 * g[11] + r.uniform(0.0, 0.4) * g[11];
    let side = if trial % 2 == 0 { 1.0 } else { -1.0 };
    g[4] = g[2] + b * cos;
    g[5] = g[3] + b * sin;
    // C need not be perpendicular to AB: gridset uses its radial distance.
    let tilt = if trial % 3 == 0 {
        r.uniform(-0.7, 0.7)
    } else {
        0.0
    };
    g[6] = g[2] + c * (tilt * cos - side * sin);
    g[7] = g[3] + c * (tilt * sin + side * cos);
    g
}

// Conversions use C-prepared fixtures so errors in Rust gridset cannot mask
// conversion errors. gridset is tested independently below.
fn prepared(r: &mut XorShift, trial: usize) -> Grid {
    let mut g = raw_grid(r, trial);
    let mut status = 777;
    unsafe { gridset(g.as_mut_ptr(), &mut status) };
    assert!(status == 0 || status == -1);
    if trial % SIZES.len() == 0 {
        // gridset rejects one-column grids. Exercise converters on a consistent
        // processed one-column definition too, including the 1x1 case.
        g[12] = 1.0;
        g[4] = g[2];
        g[5] = g[3];
        g[8] = g[6];
        g[9] = g[7];
        g[15] = g[13];
    }
    g
}

fn index_cases(g: &Grid, r: &mut XorShift) -> Vec<(i32, i32)> {
    let (nb, nc) = (g[12] as i32, g[13] as i32);
    vec![
        (1, 1),
        (nb, 1),
        (1, nc),
        (nb, nc),
        (0, 0),
        (-1, -1),
        (0, 1),
        (1, 0),
        (nb + 1, nc),
        (nb, nc + 1),
        (r.range_i32(1, nb), r.range_i32(1, nc)),
        (
            r.range_i32(-100_000, 100_000),
            r.range_i32(-100_000, 100_000),
        ),
        (i32::MIN, i32::MAX),
        (i32::MAX, i32::MIN),
    ]
}

#[test]
fn gridset_matches_c() {
    let mut r = XorShift::new(0xabc123);
    let mut statuses = [0; 5];
    for trial in 0..TRIALS * 4 {
        let mut c = raw_grid(&mut r, trial);
        match trial % 16 {
            0 => c[10] = 0.0,
            1 => c[10] = -r.uniform(0.001, 1.0e6),
            2 => c[11] = -0.0,
            3 => c[11] = -r.uniform(0.001, 1.0e6),
            4 => {
                c[4] = c[2];
                c[5] = c[3];
            }
            5 => {
                c[4] = c[2] + 0.49 * c[10];
                c[5] = c[3];
            }
            6 => {
                c[6] = c[2];
                c[7] = c[3];
            }
            7..=10 => {
                // Exact and immediately adjacent half-cell rounding thresholds.
                c[2] = 0.0;
                c[3] = 0.0;
                c[10] = 1.0;
                c[11] = 1.0;
                let half = r.range_i32(0, 32767) as f64 + 0.5;
                c[4] = match trial % 16 {
                    7 => f64::from_bits(half.to_bits() - 1),
                    8 => half,
                    _ => f64::from_bits(half.to_bits() + 1),
                };
                c[5] = 0.0;
                c[6] = 0.0;
                c[7] = if trial % 2 == 0 { -c[4] } else { c[4] };
            }
            11 => {
                c[10] = -1.0;
                c[11] = -1.0;
            }
            _ => {}
        }
        let mut mine = c;
        let (mut ce, mut re) = (123, -456);
        unsafe { gridset(c.as_mut_ptr(), &mut ce) };
        rust::gridset(&mut mine, &mut re);
        assert_eq!(re, ce, "trial {trial}");
        same_grid(&mine, &c);
        statuses[(ce + 1) as usize] += 1;
    }
    assert!(
        statuses.iter().all(|&n| n > 0),
        "missing status: {statuses:?}"
    );
}

#[test]
fn gridiccdp_matches_c() {
    let mut r = XorShift::new(2);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        for (i, j) in index_cases(&g, &mut r) {
            let (mut c, mut mine) = (0, 0);
            unsafe { gridiccdp(g.as_mut_ptr(), i, j, &mut c) };
            rust::gridiccdp(&g, i, j, &mut mine);
            assert_eq!(mine, c, "trial={trial}, index=({i},{j}), grid={g:?}");
        }
    }
}

#[test]
fn gridiccdp90_matches_c() {
    let mut r = XorShift::new(3);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        for (i, j) in index_cases(&g, &mut r) {
            let (mut c, mut mine) = (0, 0);
            unsafe { gridiccdp90(g.as_mut_ptr(), i, j, &mut c) };
            rust::gridiccdp90(&g, i, j, &mut mine);
            assert_eq!(mine, c, "trial={trial}, index=({i},{j}), grid={g:?}");
        }
    }
}

#[test]
fn gridcdpic_matches_c() {
    let mut r = XorShift::new(4);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        let last = g[15] as i32;
        for cdp in [
            i32::MIN,
            -2147483645,
            -1,
            0,
            1,
            last,
            last + 1,
            i32::MAX,
            g[12] as i32,
            g[13] as i32,
            r.range_i32(1, last),
            r.range_i32(-100_000, last + 100_000),
        ] {
            let (mut ci, mut cj, mut ri, mut rj) = (0, 0, 0, 0);
            unsafe { gridcdpic(g.as_mut_ptr(), cdp, &mut ci, &mut cj) };
            rust::gridcdpic(&g, cdp, &mut ri, &mut rj);
            assert_eq!((ri, rj), (ci, cj), "trial={trial}, cdp={cdp}, grid={g:?}");
        }
    }
}

#[test]
fn gridcdpic90_matches_c() {
    let mut r = XorShift::new(5);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        let last = g[15] as i32;
        for cdp in [
            i32::MIN,
            -2147483645,
            -1,
            0,
            1,
            last,
            last + 1,
            i32::MAX,
            g[12] as i32,
            g[13] as i32,
            r.range_i32(1, last),
            r.range_i32(-100_000, last + 100_000),
        ] {
            let (mut ci, mut cj, mut ri, mut rj) = (0, 0, 0, 0);
            unsafe { gridcdpic90(g.as_mut_ptr(), cdp, &mut ci, &mut cj) };
            rust::gridcdpic90(&g, cdp, &mut ri, &mut rj);
            assert_eq!((ri, rj), (ci, cj), "trial={trial}, cdp={cdp}, grid={g:?}");
        }
    }
}

#[test]
fn gridicgridxy_matches_c() {
    let mut r = XorShift::new(6);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        for (i, j) in index_cases(&g, &mut r) {
            let (mut cx, mut cy, mut rx, mut ry) = (0.0, 0.0, 0.0, 0.0);
            unsafe { gridicgridxy(g.as_mut_ptr(), i, j, &mut cx, &mut cy) };
            rust::gridicgridxy(&g, i, j, &mut rx, &mut ry);
            same(rx, cx);
            same(ry, cy);
        }
    }
}

#[test]
fn gridicrawxy_matches_c() {
    let mut r = XorShift::new(7);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        for (i, j) in index_cases(&g, &mut r) {
            let (mut cx, mut cy, mut rx, mut ry) = (0.0, 0.0, 0.0, 0.0);
            unsafe { gridicrawxy(g.as_mut_ptr(), i, j, &mut cx, &mut cy) };
            rust::gridicrawxy(&g, i, j, &mut rx, &mut ry);
            same(rx, cx);
            same(ry, cy);
        }
    }
}

#[test]
fn gridrawxygridxy_matches_c() {
    let mut r = XorShift::new(8);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        let extent = (g[12] * g[10]).max(g[13] * g[11]);
        for (x, y) in [
            (0.0, -0.0),
            (-0.0, 0.0),
            (g[2], g[3]),
            (g[4], g[5]),
            (g[6], g[7]),
            (g[8], g[9]),
            (r.uniform(-2.0, 2.0) * extent, r.uniform(-2.0, 2.0) * extent),
            (r.uniform(-1.0e9, 1.0e9), r.uniform(-1.0e9, 1.0e9)),
            (1.0e150, -1.0e150),
            (1.0e-300, -1.0e-300),
            (f64::INFINITY, f64::NEG_INFINITY),
            (f64::NAN, 0.0),
        ] {
            let (mut cx, mut cy, mut rx, mut ry) = (0.0, 0.0, 0.0, 0.0);
            unsafe { gridrawxygridxy(g.as_mut_ptr(), x, y, &mut cx, &mut cy) };
            rust::gridrawxygridxy(&g, x, y, &mut rx, &mut ry);
            same(rx, cx);
            same(ry, cy);
        }
    }
}

#[test]
fn gridgridxyrawxy_matches_c() {
    let mut r = XorShift::new(9);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        let extent = (g[12] * g[10]).max(g[13] * g[11]);
        for (x, y) in [
            (0.0, -0.0),
            (-0.0, 0.0),
            (g[2], g[3]),
            (g[4], g[5]),
            (g[6], g[7]),
            (g[8], g[9]),
            (r.uniform(-2.0, 2.0) * extent, r.uniform(-2.0, 2.0) * extent),
            (r.uniform(-1.0e9, 1.0e9), r.uniform(-1.0e9, 1.0e9)),
            (1.0e150, -1.0e150),
            (1.0e-300, -1.0e-300),
            (f64::INFINITY, f64::NEG_INFINITY),
            (f64::NAN, 0.0),
        ] {
            let (mut cx, mut cy, mut rx, mut ry) = (0.0, 0.0, 0.0, 0.0);
            unsafe { gridgridxyrawxy(g.as_mut_ptr(), x, y, &mut cx, &mut cy) };
            rust::gridgridxyrawxy(&g, x, y, &mut rx, &mut ry);
            same(rx, cx);
            same(ry, cy);
        }
    }
}

#[test]
fn gridrawxycdpic_matches_c() {
    let mut r = XorShift::new(10);
    for trial in 0..TRIALS {
        let mut g = prepared(&mut r, trial);
        let nb = g[12];
        let nc = g[13];
        let mut raw = vec![
            (0.0, -0.0),
            (g[2], g[3]),
            (g[4], g[5]),
            (g[6], g[7]),
            (g[8], g[9]),
        ];
        for (i, j) in [
            (-1.0, -1.0),
            (-0.5, -0.5),
            (0.5, 0.5),
            (nb - 0.5, nc - 0.5),
            (nb, nc),
            (r.uniform(-2.0, nb + 2.0), r.uniform(-2.0, nc + 2.0)),
            (r.uniform(0.0, nb), r.uniform(0.0, nc)),
            (-100_000.0, 100_000.0),
        ] {
            let (mut x, mut y) = (0.0, 0.0);
            let (gx, gy) = (i * g[10], j * g[11]);
            unsafe { gridgridxyrawxy(g.as_mut_ptr(), gx, gy, &mut x, &mut y) };
            raw.push((x, y));
            // Both floating-point neighbors around cell edges, including negatives.
            raw.push((x.next_up(), y.next_up()));
            raw.push((x.next_down(), y.next_down()));
        }
        for (x, y) in raw {
            let (mut cc, mut ci, mut cj, mut rc, mut ri, mut rj) = (0, 0, 0, 0, 0, 0);
            unsafe { gridrawxycdpic(g.as_mut_ptr(), x, y, &mut cc, &mut ci, &mut cj) };
            rust::gridrawxycdpic(&g, x, y, &mut rc, &mut ri, &mut rj);
            assert_eq!(
                (rc, ri, rj),
                (cc, ci, cj),
                "trial={trial}, xy=({x},{y}), grid={g:?}"
            );
        }
    }
}

#[test]
fn gridcheck_matches_c() {
    // Run diagnostic calls in two subprocesses. This compares the complete
    // stderr byte stream without redirecting descriptors shared by other tests.
    const WORKER: &str = "OXWAVE_GRIDCHECK_DIFF_WORKER";
    if let Ok(mode) = std::env::var(WORKER) {
        assert!(mode == "c" || mode == "rust");
        let mut argv = CArgv::new(&["gridcheck-diff"]);
        unsafe { oxwave_cref::initargs(argv.argc(), argv.argv()) };
        oxwave_core::par::initargs(&["gridcheck-diff".to_owned()]);
        let mut r = XorShift::new(11);
        for trial in 0..TRIALS {
            let mut g = prepared(&mut r, trial);
            let before = g;
            let icheck = [0, 1, -1, 2, i32::MIN, i32::MAX][trial % 6];
            let mut status = r.next_u32() as i32;
            if mode == "c" {
                unsafe { gridcheck(g.as_mut_ptr(), icheck, &mut status) };
            } else {
                rust::gridcheck(&g, icheck, &mut status);
            }
            assert_eq!(status, 0);
            same_grid(&g, &before);
        }
        return;
    }
    let run = |mode| {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "gridcheck_matches_c", "--nocapture"])
            .env(WORKER, mode)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{mode} worker failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        output.stderr
    };
    let c = run("c");
    let mine = run("rust");
    assert!(
        c.len() > 1_000_000,
        "expected diagnostics from thousands of calls"
    );
    assert_eq!(mine.len(), c.len(), "diagnostic lengths differ");
    if let Some(i) = mine.iter().zip(&c).position(|(a, b)| a != b) {
        let start = i.saturating_sub(120);
        let end = (i + 120).min(c.len());
        panic!(
            "diagnostic byte {i} differs: Rust={:?}, C={:?}",
            String::from_utf8_lossy(&mine[start..end]),
            String::from_utf8_lossy(&c[start..end])
        );
    }
}
