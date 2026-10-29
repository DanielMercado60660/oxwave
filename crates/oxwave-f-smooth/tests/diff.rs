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

//! Differential tests against par/lib/smooth.c.
//! Credits: Gabriel Alvarez (1995), Zhenyue Liu (1992, 1994),
//! and Balasz Nemeth, Potash Corporation (2008).
//! Only the source header and declarations were consulted in this round.

use oxwave_cref::XorShift;
use oxwave_f_smooth as rust;
use std::sync::Mutex;

static C_LOCK: Mutex<()> = Mutex::new(());

struct CwdGuard(std::path::PathBuf);
impl CwdGuard {
    fn scratch(name: &str) -> Self {
        let previous = std::env::current_dir().unwrap();
        let directory = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../.targets/lib-f_smooth")
            .join(name);
        std::fs::create_dir_all(&directory).unwrap();
        std::env::set_current_dir(directory).unwrap();
        Self(previous)
    }
}
impl Drop for CwdGuard {
    fn drop(&mut self) {
        std::env::set_current_dir(&self.0).unwrap();
    }
}

fn error_file() -> Vec<f32> {
    let bytes = std::fs::read("error").unwrap();
    assert_eq!(bytes.len() % 4, 0);
    bytes
        .chunks_exact(4)
        .map(|b| f32::from_ne_bytes(b.try_into().unwrap()))
        .collect()
}

extern "C" {
    fn SG_smoothing_filter(np: i32, nl: i32, nr: i32, ld: i32, m: i32, filter: *mut f32);
    fn dlsq_smoothing(
        nt: i32,
        nx: i32,
        ift: i32,
        ilt: i32,
        ifx: i32,
        ilx: i32,
        r1: f32,
        r2: f32,
        rw: f32,
        traces: *mut *mut f32,
    );
    fn gaussian1d_smoothing(ns: i32, nsr: i32, data: *mut f32);
    fn gaussian2d_smoothing(nx: i32, nt: i32, nsx: i32, nst: i32, data: *mut *mut f32);
    fn rwa_smoothing_filter(flag: i32, nl: i32, nr: i32, filter: *mut f32);
    fn smooth_1(x: *mut f32, z: *mut f32, r: f32, n: i32);
    fn smooth_histogram(nintlh: i32, pdf: *mut f32);
    fn smooth_segmented_array(index: *mut f32, val: *mut f32, n: i32, sm: i32, inc: i32, m: i32);
}

fn compare(actual: &[f32], expected: &[f32]) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected
        .iter()
        .filter(|x| x.is_finite())
        .fold(0.0_f32, |v, x| v.max(x.abs()));
    for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        if b.is_nan() {
            assert!(a.is_nan(), "sample {i}: {a} != NaN");
        } else if b.is_infinite() {
            assert_eq!(a, b, "sample {i}");
        } else {
            assert!(
                a.is_finite() && (a - b).abs() <= 1e-5 * scale + 1e-30,
                "sample {i}: Rust={a}, C={b}, scale={scale}, error={}",
                (a - b).abs()
            );
        }
    }
}

fn size(rng: &mut XorShift, case: usize) -> usize {
    // Includes singleton, pair, primes, powers of two and a large odd array.
    const EDGES: &[usize] = &[1, 2, 3, 5, 7, 16, 31, 64, 127, 256, 1021];
    if case % 4 == 0 {
        EDGES[(case / 4) % EDGES.len()]
    } else {
        rng.range_i32(1, 129) as usize
    }
}

fn samples(rng: &mut XorShift, n: usize, case: usize) -> Vec<f32> {
    let scale = [1.0, 1e-10, 1e10, 0.01][(case / 7) % 4];
    (0..n)
        .map(|i| {
            scale
                * match case % 7 {
                    0 => 0.0,
                    1 => -2.5,
                    2 => {
                        if i == n / 2 {
                            -1.0
                        } else {
                            0.0
                        }
                    }
                    3 => {
                        if i % 2 == 0 {
                            1.0
                        } else {
                            -1.0
                        }
                    }
                    _ => rng.uniform_f32(-10.0, 10.0),
                }
        })
        .collect()
}

#[test]
fn sg_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(11);
    for case in 0..6000 {
        let (nl, nr) = if case % 4 == 0 {
            let n = size(&mut rng, case);
            (n / 2, n - 1 - n / 2)
        } else {
            (rng.range_i32(0, 40) as usize, rng.range_i32(0, 40) as usize)
        };
        let m = rng.range_i32(0, (nl + nr).min(4) as i32) as usize;
        let ld = rng.range_i32(0, m as i32) as usize;
        let np = nl + nr + 1 + rng.range_i32(0, 7) as usize;
        let mut c = vec![-123.0; np];
        let mut r = c.clone();
        unsafe {
            SG_smoothing_filter(
                np as i32,
                nl as i32,
                nr as i32,
                ld as i32,
                m as i32,
                c.as_mut_ptr(),
            );
        }
        rust::SG_smoothing_filter(nl, nr, ld, m, &mut r);
        compare(&r, &c);
    }
}

#[test]
fn rwa_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(10);
    for case in 0..6000 {
        let (nl, nr) = if case % 3 == 0 {
            let n = size(&mut rng, case);
            if case % 2 == 0 {
                (0, n - 1)
            } else {
                (n - 1, 0)
            }
        } else {
            (
                rng.range_i32(0, 128) as usize,
                rng.range_i32(0, 128) as usize,
            )
        };
        let flag = 1 + (case % 2) as i32;
        let initial = [0.0, 1.0, -123.0][case % 3];
        let mut c = vec![initial; nl + nr + 4];
        let mut r = c.clone();
        unsafe {
            rwa_smoothing_filter(flag, nl as i32, nr as i32, c.as_mut_ptr());
        }

        rust::rwa_smoothing_filter(flag, nl, nr, &mut r);
        compare(&r, &c);
        if flag == 1 {
            assert_eq!(r, c);
        }
    }
}

#[test]
fn gaussian1d_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(98);
    for case in 0..4000 {
        let n = size(&mut rng, case);
        let radius = if case % 5 == 0 {
            n + case % 3
        } else {
            rng.range_i32(0, 30) as usize
        };
        let mut c = samples(&mut rng, n, case);
        let mut r = c.clone();
        unsafe {
            gaussian1d_smoothing(n as i32, radius as i32, c.as_mut_ptr());
        }
        rust::gaussian1d_smoothing(radius, &mut r);
        compare(&r, &c);
    }
}

#[test]
fn gaussian2d_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(99);
    for case in 0..3000 {
        let nx = [1, 2, 3, 7, 16, 31, 64][case % 7];
        let nt = [1, 2, 5, 13, 32, 67, 257][(case / 7) % 7];
        let nsx = if case % 3 == 0 {
            nx
        } else {
            rng.range_i32(0, 15) as usize
        };
        let nst = if case % 5 == 0 {
            nt + 1
        } else {
            rng.range_i32(0, 15) as usize
        };
        let mut c: Vec<Vec<f32>> = samples(&mut rng, nx * nt, case)
            .chunks(nt)
            .map(|r| r.to_vec())
            .collect();
        let mut r = c.clone();
        let mut ptrs: Vec<_> = c.iter_mut().map(|r| r.as_mut_ptr()).collect();
        unsafe {
            gaussian2d_smoothing(
                nx as i32,
                nt as i32,
                nsx as i32,
                nst as i32,
                ptrs.as_mut_ptr(),
            );
        }
        rust::gaussian2d_smoothing(nsx, nst, &mut r);
        compare(&r.concat(), &c.concat());
    }
}

#[test]
fn histogram_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(100);
    for case in 0..3000 {
        let n = size(&mut rng, case);
        // Histograms normally have nonnegative counts; exercise negative
        // and zero counts as well, including the C zero-sum NaNs.
        let mut c: Vec<f32> = (0..n)
            .map(|_| match case % 4 {
                0 => 0.0,
                1 => rng.uniform_f32(0.1, 1000.0),
                2 => rng.uniform_f32(-1000.0, -0.1),
                _ => rng.uniform_f32(-100.0, 1000.0),
            })
            .collect();
        let mut r = c.clone();
        unsafe {
            smooth_histogram(n as i32, c.as_mut_ptr());
        }
        rust::smooth_histogram(&mut r);
        compare(&r, &c);
    }
    for input in [vec![-1.0, 0.0, 1.0], vec![0.0], vec![-2.0]] {
        let mut c = input.clone();
        let mut r = input;
        unsafe {
            smooth_histogram(c.len() as i32, c.as_mut_ptr());
        }
        rust::smooth_histogram(&mut r);
        compare(&r, &c);
    }
}

#[test]
fn smooth_1_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(328);
    for case in 0..4000 {
        let n = size(&mut rng, case);
        let mut x = vec![-100.0; n];
        let direction = if case % 2 == 0 { 1.0 } else { -1.0 };
        for i in 1..n {
            x[i] = x[i - 1] + direction * rng.uniform_f32(0.1, 3.0);
        }
        let strength = if case % 9 == 0 {
            0.0
        } else {
            rng.uniform_f32(-20.0, 20.0)
        };
        let mut c = samples(&mut rng, n, case);
        let mut r = c.clone();
        let original_x = x.clone();
        unsafe {
            smooth_1(x.as_mut_ptr(), c.as_mut_ptr(), strength, n as i32);
        }
        rust::smooth_1(&original_x, &mut r, strength);
        assert_eq!(x, original_x);
        compare(&r, &c);
    }
}

#[test]
fn dlsq_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let _cwd = CwdGuard::scratch("diff-dlsq");
    let mut rng = XorShift::new(129);
    for case in 0..4000 {
        let nx = [1, 2, 3, 7, 16, 31, 64][case % 7];
        let nt = [2, 3, 5, 13, 32, 67, 257][(case / 7) % 7];
        let (ift, ilt, ifx, ilx) = if case % 3 == 0 {
            (0, nt, 0, nx)
        } else {
            let ift = rng.range_i32(0, nt as i32 - 1) as usize;
            let ifx = rng.range_i32(0, nx as i32 - 1) as usize;
            (
                ift,
                rng.range_i32(ift as i32 + 1, nt as i32) as usize,
                ifx,
                rng.range_i32(ifx as i32 + 1, nx as i32) as usize,
            )
        };
        let r1 = if case % 5 == 0 {
            0.0
        } else {
            rng.uniform_f32(-20.0, 20.0)
        };
        let r2 = if case % 7 == 0 {
            0.0
        } else {
            rng.uniform_f32(-20.0, 20.0)
        };
        let rw = if case % 2 == 0 {
            0.0
        } else {
            rng.uniform_f32(0.0, 5.0)
        };
        let mut c: Vec<Vec<f32>> = samples(&mut rng, nx * nt, case)
            .chunks(nt)
            .map(|r| r.to_vec())
            .collect();
        let mut r = c.clone();
        let mut ptrs: Vec<_> = c.iter_mut().map(|r| r.as_mut_ptr()).collect();
        unsafe {
            dlsq_smoothing(
                nt as i32,
                nx as i32,
                ift as i32,
                ilt as i32,
                ifx as i32,
                ilx as i32,
                r1,
                r2,
                rw,
                ptrs.as_mut_ptr(),
            );
        }
        let c_error = error_file();
        rust::dlsq_smoothing(ift, ilt, ifx, ilx, r1, r2, rw, &mut r);
        let r_error = error_file();
        assert_eq!(r_error.len(), nt);
        assert_eq!(c_error.len(), nt);
        let scale = c_error.iter().fold(0.0_f32, |a, b| a.max(b.abs()));
        for (&a, &b) in r_error.iter().zip(&c_error) {
            // Subtraction magnifies sub-ULP differences in the smoothed
            // data when the relative RMS error is close to zero.
            assert!(
                (a - b).abs() <= 1e-5 * scale + 1e-5,
                "error file: Rust={a}, C={b}, scale={scale}"
            );
        }
        compare(&r.concat(), &c.concat());
    }
    // The C time solve accesses outside its allocation when nt=1 and r1!=0.
    // Test singleton time axes with that pass disabled.
    for nx in [1, 2, 3, 7, 16] {
        let mut c: Vec<Vec<f32>> = (0..nx).map(|i| vec![i as f32 - 2.0]).collect();
        let mut r = c.clone();
        let mut ptrs: Vec<_> = c.iter_mut().map(|r| r.as_mut_ptr()).collect();
        unsafe {
            dlsq_smoothing(
                1,
                nx as i32,
                0,
                1,
                0,
                nx as i32,
                0.0,
                2.0,
                0.0,
                ptrs.as_mut_ptr(),
            );
        }
        rust::dlsq_smoothing(0, 1, 0, nx, 0.0, 2.0, 0.0, &mut r);
        compare(&r.concat(), &c.concat());
    }
}

#[test]
fn segmented_matches_c() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(221);
    for case in 0..15000 {
        let n = size(&mut rng, case);
        let m = [1, 3, 4, 5, 6][case % 5];
        let sm = rng.range_i32(if m == 3 || m == 6 { 2 } else { 1 }, 16) as usize;
        let inc = rng.range_i32(1, 4);
        let mut x: Vec<f32> = (0..n).map(|i| -100.0 + i as f32 * inc as f32).collect();
        // The even-sized gapped regression below corrupts C's heap.
        // Odd total lengths with this split give reproducible C outputs.
        if n > 2 && n % 2 == 1 && case % 2 == 0 {
            for v in &mut x[n / 2..] {
                *v += 5.0 * inc as f32;
            }
        }
        let mut c = samples(&mut rng, n, case);
        let mut r = c.clone();
        let original_x = x.clone();
        unsafe {
            smooth_segmented_array(x.as_mut_ptr(), c.as_mut_ptr(), n as i32, sm as i32, inc, m);
        }
        rust::smooth_segmented_array(&original_x, &mut r, sm, inc, m);
        assert_eq!(x, original_x);
        compare(&r, &c);
        if m == 6 {
            assert_eq!(r, c);
        }
    }
    // Method 2 is defined for an unprocessed final singleton.
    for _ in 0..2000 {
        let mut x = vec![rng.uniform_f32(-100.0, 100.0)];
        let mut c = vec![rng.uniform_f32(-100.0, 100.0)];
        let mut r = c.clone();
        unsafe {
            smooth_segmented_array(x.as_mut_ptr(), c.as_mut_ptr(), 1, 3, 1, 2);
        }
        rust::smooth_segmented_array(&x, &mut r, 3, 1, 2);
        assert_eq!(r, c);
    }
}

/// C has no numerical result for these inputs. Run it in a separate process
/// so heap corruption cannot contaminate the other differential comparisons.
#[test]
fn segmented_c_failure_modes() {
    use std::os::unix::process::ExitStatusExt;
    const CHILD: &str = "OXWAVE_SMOOTH_C_FAILURE";
    if let Ok(mode) = std::env::var(CHILD) {
        let mut rng = XorShift::new(221);
        for _ in 0..256 {
            let (n, method) = if mode == "dlsq" { (21, 2) } else { (16, 1) };
            let mut x: Vec<f32> = (0..n)
                .map(|i| (i + if mode == "gapped" && i >= n / 2 { 5 } else { 0 }) as f32)
                .collect();
            let mut c: Vec<f32> = (0..n).map(|_| rng.uniform_f32(-10.0, 10.0)).collect();
            let mut r = c.clone();
            unsafe {
                smooth_segmented_array(x.as_mut_ptr(), c.as_mut_ptr(), n as i32, 3, 1, method);
            }
            rust::smooth_segmented_array(&x, &mut r, 3, 1, method);
            // The gapped call can return matching values before its heap
            // corruption is detected on a subsequent allocation/free.
            if mode == "gapped" {
                compare(&r, &c);
            }
        }
        if mode == "dlsq" {
            panic!("C unexpectedly survived the {mode} regression");
        }
        return;
    }
    for mode in ["dlsq", "gapped"] {
        let directory = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../.targets/lib-f_smooth/diff-c-faults");
        std::fs::create_dir_all(&directory).unwrap();
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "segmented_c_failure_modes", "--nocapture"])
            .env(CHILD, mode)
            .current_dir(directory)
            .output()
            .unwrap();
        let signal = output.status.signal();
        // Heap corruption can remain undetected with a different allocator
        // layout. A normal exit means every gapped output was compared.
        assert!(
            matches!(signal, Some(6 | 10 | 11)) || (mode == "gapped" && output.status.success()),
            "expected a C memory fault for {mode}, status={}, stderr={}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
    }
}
