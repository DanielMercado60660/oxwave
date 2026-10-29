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

//! Differential tests against cwp/lib/sge.c (LINPACK FORTRAN adaptations).
//! Credits: Dave Hale, Colorado School of Mines, 10/01/89.

use oxwave_cref::XorShift;
use oxwave_lu as rust;

extern "C" {
    fn sgefa(a: *mut *mut f32, n: i32, ipvt: *mut i32, info: *mut i32);
    fn sgeco(a: *mut *mut f32, n: i32, ipvt: *mut i32, rcond: *mut f32, z: *mut f32);
    fn sgesl(a: *mut *mut f32, n: i32, ipvt: *mut i32, b: *mut f32, job: i32);
}

// Separate column allocations exercise float ** rather than assuming a flat array.
// Every allocation has guards; even the extra column must remain untouched.
fn matrix(rng: &mut XorShift, n: usize, kind: usize, scale: f32) -> Vec<Vec<f32>> {
    let mut a = vec![vec![12345.0; n + 3]; n + 1];
    for j in 0..n {
        for i in 0..n {
            let x = rng.uniform_f32(-1.0, 1.0);
            a[j][i] = match kind {
                0 => x, // Dense matrices force partial pivoting.
                1 => {
                    if i == j {
                        x + if x < 0.0 { -2.0 } else { 2.0 }
                    } else {
                        0.0
                    }
                }
                2 => 0.0,
                3 => {
                    if i <= j {
                        x
                    } else {
                        0.0
                    }
                }
                4 => {
                    if i >= j {
                        x
                    } else {
                        0.0
                    }
                }
                5 => {
                    if j == 0 {
                        0.0
                    } else {
                        x
                    }
                } // Initial zero pivot.
                6 => {
                    if i == n - 1 {
                        0.0
                    } else {
                        x
                    }
                } // Final zero pivot.
                7 => x, // Duplicate column below.
                8 => {
                    if i == j {
                        x + n as f32
                    } else {
                        x
                    }
                }
                9 => {
                    if i == j {
                        if i % 2 == 0 {
                            1.0
                        } else {
                            1e-6
                        }
                    } else {
                        x * 1e-7
                    }
                }
                10 => {
                    if i == j {
                        1.0
                    } else {
                        0.0
                    }
                } // Swapped below.
                11 => {
                    if (i + j) % 2 == 0 {
                        1.0
                    } else {
                        -1.0
                    }
                } // Pivot ties.
                12 => {
                    if (i + j) % 2 == 0 {
                        0.0
                    } else {
                        -0.0
                    }
                }
                13 => {
                    if i == j {
                        -1.0
                    } else {
                        0.0
                    }
                }
                14 => 1.0 / (i + j + 1) as f32, // Ill-conditioned Hilbert matrix.
                _ => unreachable!(),
            } * scale;
        }
    }
    if kind == 7 && n > 1 {
        for i in 0..n {
            a[n - 1][i] = a[0][i];
        }
    }
    if kind == 10 && n > 1 {
        for j in 0..n {
            a[j].swap(0, n - 1);
        }
    }
    a
}

fn cases(mut f: impl FnMut(usize, usize, Vec<Vec<f32>>, &mut XorShift)) {
    let mut rng = XorShift::new(0x5347_454c_5546_4131);
    let sizes = [1, 2, 3, 5, 7, 8, 11, 16, 17, 23, 31, 32];
    let scales = [1.0, 1e-30, 1e30, 1e-15, 1e15, -1.0, 0.01, 1000.0];
    for case in 0..3600 {
        let n = sizes[case % sizes.len()];
        let kind = (case / sizes.len()) % 15;
        let scale = scales[rng.next_u32() as usize % scales.len()];
        let a = matrix(&mut rng, n, kind, scale);
        f(case, n, a, &mut rng);
    }
    let mut case = 3600;
    for n in [61, 64, 127, 128, 257] {
        for kind in [0, 8, 10] {
            let a = matrix(&mut rng, n, kind, 1.0);
            f(case, n, a, &mut rng);
            case += 1;
        }
    }
    // Subnormal, extreme finite, and exceptional scalars have defined IEEE
    // results in these routines. No singular or nonfinite case is skipped.
    for x in [
        0.0,
        -0.0,
        1.0,
        -1.0,
        f32::MIN_POSITIVE,
        -f32::MIN_POSITIVE,
        f32::from_bits(1),
        -f32::from_bits(1),
        f32::MAX,
        -f32::MAX,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
    ] {
        f(case, 1, vec![vec![x, 12345.0], vec![12345.0; 2]], &mut rng);
        case += 1;
    }
    for n in [2, 5, 16] {
        for scale in [f32::from_bits(1), f32::MIN_POSITIVE, 1e-38, 1e37] {
            for kind in [0, 1, 9] {
                let a = matrix(&mut rng, n, kind, scale);
                f(case, n, a, &mut rng);
                case += 1;
            }
        }
    }
}

fn pointers(a: &mut [Vec<f32>]) -> Vec<*mut f32> {
    a.iter_mut().map(|column| column.as_mut_ptr()).collect()
}

fn views(a: &mut [Vec<f32>]) -> Vec<&mut [f32]> {
    a.iter_mut().map(|column| column.as_mut_slice()).collect()
}

// Require identical bits for finite values, signed zero, and infinities.
// NaN payload/sign are not specified by C, so compare their classification.
fn same(got: &[f32], expected: &[f32], case: usize, what: &str) {
    assert_eq!(got.len(), expected.len());
    for (i, (&r, &c)) in got.iter().zip(expected).enumerate() {
        assert!(
            r.to_bits() == c.to_bits() || (r.is_nan() && c.is_nan()),
            "case {case} {what}[{i}]: Rust {r:?} ({:08x}), C {c:?} ({:08x})",
            r.to_bits(),
            c.to_bits()
        );
    }
}

fn same_matrix(got: &[Vec<f32>], expected: &[Vec<f32>], case: usize) {
    assert_eq!(got.len(), expected.len());
    for (j, (r, c)) in got.iter().zip(expected).enumerate() {
        same(r, c, case, &format!("a[{j}]"));
    }
}

#[test]
fn sgefa_matches_c() {
    cases(|case, n, mut a, _| {
        let mut c = a.clone();
        let mut cp = vec![-99; n + 2];
        let mut rp = cp.clone();
        let (mut ci, mut ri) = (-99, -99);
        let mut ptrs = pointers(&mut c);
        unsafe { sgefa(ptrs.as_mut_ptr(), n as i32, cp.as_mut_ptr(), &mut ci) };
        rust::sgefa(&mut views(&mut a), n as i32, &mut rp, &mut ri);
        assert_eq!(ri, ci, "case {case} info");
        assert_eq!(rp, cp, "case {case} pivots");
        same_matrix(&a, &c, case);
    });
}

#[test]
fn sgeco_matches_c() {
    cases(|case, n, mut a, _| {
        let mut c = a.clone();
        let mut cp = vec![-99; n + 2];
        let mut rp = cp.clone();
        let mut cz = vec![-987.0; n + 3];
        let mut rz = cz.clone();
        let (mut cc, mut rc) = (-99.0, -99.0);
        let mut ptrs = pointers(&mut c);
        unsafe {
            sgeco(
                ptrs.as_mut_ptr(),
                n as i32,
                cp.as_mut_ptr(),
                &mut cc,
                cz.as_mut_ptr(),
            )
        };
        rust::sgeco(&mut views(&mut a), n as i32, &mut rp, &mut rc, &mut rz);
        assert_eq!(rp, cp, "case {case} pivots");
        same_matrix(&a, &c, case);
        same(&[rc], &[cc], case, "rcond");
        same(&rz, &cz, case, "z");
    });
}

#[test]
fn sgesl_matches_c() {
    cases(|case, n, mut a, rng| {
        // Use C factors to isolate the solve from the Rust factorization.
        let mut ipvt = vec![-99; n + 2];
        let mut info = -99;
        let mut ptrs = pointers(&mut a);
        unsafe { sgefa(ptrs.as_mut_ptr(), n as i32, ipvt.as_mut_ptr(), &mut info) };
        let before = a.clone();
        let before_ipvt = ipvt.clone();
        let jobs = [0, 1, -1, 2, i32::MIN, i32::MAX];
        for job in [0, jobs[1 + case % (jobs.len() - 1)]] {
            // Zero RHS, signed RHS, and independent large/small RHS scales.
            // Draw independently of the size cycle so every dimension also
            // exercises nonzero right-hand sides at each magnitude.
            let scale = [0.0, 1.0, 1e-30, 1e30][rng.next_u32() as usize % 4];
            let mut cb: Vec<f32> = (0..n).map(|_| rng.uniform_f32(-1.0, 1.0) * scale).collect();
            cb.extend_from_slice(&[12345.0, -0.0]);
            let mut rb = cb.clone();
            unsafe {
                sgesl(
                    ptrs.as_mut_ptr(),
                    n as i32,
                    ipvt.as_mut_ptr(),
                    cb.as_mut_ptr(),
                    job,
                )
            };
            let cols: Vec<&[f32]> = a.iter().map(|column| column.as_slice()).collect();
            rust::sgesl(&cols, n as i32, &ipvt, &mut rb, job);
            same(&rb, &cb, case, &format!("b (job={job})"));
        }
        same_matrix(&a, &before, case);
        assert_eq!(ipvt, before_ipvt, "case {case} input pivots changed");
    });
}
