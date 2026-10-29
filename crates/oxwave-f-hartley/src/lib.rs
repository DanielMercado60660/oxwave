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

//! Port of cwp/lib/hartley.c. Credits: CENPET: Werner M. Heigl, April 2006.

/// C prototype: `int nextpow2(int n);`
/// Returns the smallest exponent m such that n <= 2^m (zero for n <= 1).
/// Requires `n <= 16777216`, the upper limit of the C lookup table.
pub fn nextpow2(n: i32) -> i32 {
    assert!(n <= 16777216);
    let pow2 = [
        1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072,
        262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216,
    ];
    let mut i: i32 = 0;
    while pow2[i as usize] < n {
        i += 1;
    }
    i
}

/// C prototype: `int nextpow4(int n);`
/// Returns the smallest exponent m such that n <= 4^m (zero for n <= 1).
/// Requires `n <= 16777216`, the upper limit of the C lookup table.
pub fn nextpow4(n: i32) -> i32 {
    assert!(n <= 16777216);
    let pow4 = [
        1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216,
    ];
    let mut i: i32 = 0;
    while pow4[i as usize] < n {
        i += 1;
    }
    i
}

macro_rules! split {
    ($name:ident,$t:ty,$doc:literal) => {
        #[doc=$doc]
        /// In-place, unnormalized split-radix Hartley transform.
        /// Requires n = 2^m, 0 <= m < 30, and at least max(n, 2) elements.
        /// For n = 1 the C last stage still transforms the first two elements;
        /// supply that padding to reproduce its behavior safely. The ordinary
        /// transform domain starts at n = 2. Remaining elements are untouched.
        pub fn $name(n: i32, m: i32, f: &mut [$t]) {
            assert!(
                n >= 1
                    && m >= 0
                    && m < 30
                    && (1i64 << m) == n as i64
                    && f.len() >= n.max(2) as usize
            );
            macro_rules! a {
                ($i:expr) => {
                    f[($i - 1) as usize]
                };
            }
            let mut n2 = n.wrapping_shl(1);
            let tpi = (2.0 * std::f64::consts::PI) as $t;
            let r2 = (2.0 as $t).sqrt();
            for _ in 1..m {
                let mut is = 1;
                let mut id = n2;
                n2 /= 2;
                let n4 = n2 / 4;
                let n8 = n2 / 8;
                let e = tpi / n2 as $t;
                loop {
                    let mut i = is;
                    while if id < 0 { i >= n } else { i <= n } {
                        let l12 = i + n4;
                        let l13 = l12 + n4;
                        let l14 = l13 + n4;
                        let mut t1 = a!(l12) - a!(l14);
                        let mut t2 = a!(i) + a!(l13);
                        let mut t3 = a!(l12) + a!(l14);
                        let mut t4 = a!(i) - a!(l13);
                        a!(l14) = t4 - t1;
                        a!(l13) = t4 + t1;
                        a!(l12) = t3;
                        a!(i) = t2;
                        if n4 != 1 {
                            let l21 = i + n8;
                            let l22 = l21 + n4;
                            let l23 = l22 + n4;
                            let l24 = l23 + n4;
                            t1 = a!(l22);
                            t2 = a!(l23);
                            t3 = a!(l24);
                            a!(l24) = (t1 - t3) * r2;
                            a!(l23) = (a!(l21) - t2) * r2;
                            a!(l22) = t1 + t3;
                            a!(l21) += t2;
                            let mut v = e;
                            for j in 2..=n8 {
                                let l11 = i + j - 1;
                                let l12 = l11 + n4;
                                let l13 = l12 + n4;
                                let l14 = l13 + n4;
                                let l21 = i + n4 - j + 1;
                                let l22 = l21 + n4;
                                let l23 = l22 + n4;
                                let l24 = l23 + n4;
                                let a3 = v * 3.0;
                                let cc1 = v.cos();
                                let ss1 = v.sin();
                                let cc3 = a3.cos();
                                let ss3 = a3.sin();
                                v = j as $t * e;
                                let mut t5 = a!(l21) - a!(l23);
                                t2 = a!(l11) - a!(l13);
                                t1 = t2 + t5;
                                t2 -= t5;
                                t5 = a!(l22) - a!(l24);
                                t4 = a!(l14) - a!(l12);
                                t3 = t4 + t5;
                                t4 -= t5;
                                a!(l11) += a!(l13);
                                a!(l12) += a!(l14);
                                a!(l21) += a!(l23);
                                a!(l22) += a!(l24);
                                a!(l13) = t1 * cc1 + t3 * ss1;
                                a!(l14) = t2 * cc3 - t4 * ss3;
                                a!(l23) = t1 * ss1 - t3 * cc1;
                                a!(l24) = t2 * ss3 + t4 * cc3;
                            }
                        }
                        i = i.wrapping_add(id);
                    }
                    is = id.wrapping_shl(1).wrapping_sub(n2).wrapping_add(1);
                    id = id.wrapping_shl(2);
                    if is >= n {
                        break;
                    }
                }
            }
            let mut is = 1;
            let mut id = 4;
            loop {
                let mut i = is;
                while if id < 0 { i >= n } else { i <= n } {
                    let j = i + 1;
                    let t = a!(i);
                    a!(i) = t + a!(j);
                    a!(j) = t - a!(j);
                    i = i.wrapping_add(id);
                }
                is = id.wrapping_shl(1).wrapping_sub(1);
                id = id.wrapping_shl(2);
                if is >= n {
                    break;
                }
            }
            let mut j = 1;
            for i in 1..n {
                if i < j {
                    f.swap((i - 1) as usize, (j - 1) as usize);
                }
                let mut k = n / 2;
                while k < j {
                    j -= k;
                    k /= 2;
                }
                j += k;
            }
        }
    };
}
split!(
    srfht,
    f32,
    "C prototype: `void srfht(int *n, int *m, float *f);`"
);
split!(
    dsrfht,
    f64,
    "C prototype: `void dsrfht(int *n, int *m, double *f);`"
);

/// C prototype: `void r4fht(int n, int m, float *f);`
/// In-place, unnormalized radix-4 Hartley transform.
/// Requires n = 4^m, 0 <= m <= 15, and at least max(n, 4) elements.
/// For n = 1 the C first stage still transforms four elements; supply that
/// padding to reproduce its behavior safely. The ordinary transform domain
/// starts at n = 4. Remaining elements are untouched.
pub fn r4fht(n: i32, m: i32, f: &mut [f32]) {
    assert!(
        n >= 1
            && m >= 0
            && m <= 15
            && (1i64 << (2 * m)) == n as i64
            && f.len() >= n.max(4) as usize
    );
    let n = n as usize;
    let n4 = n / 4;
    let r2 = 2.0f64.sqrt() as f32;
    let mut j = 1;
    let mut i = 0;
    while i < n - 1 {
        i += 1;
        if i < j {
            f.swap(j - 1, i - 1);
        }
        let mut k = n4;
        while 3 * k < j {
            j -= 3 * k;
            k /= 4;
        }
        j += k;
    }
    for i in (0..n).step_by(4) {
        let t1 = f[i] + f[i + 1];
        let t2 = f[i] - f[i + 1];
        let t3 = f[i + 2] + f[i + 3];
        let t4 = f[i + 2] - f[i + 3];
        f[i] = t1 + t3;
        f[i + 1] = t1 - t3;
        f[i + 2] = t2 + t4;
        f[i + 3] = t2 - t4;
    }
    for l in 2..=m {
        let e1 = 2.0f32.powf((l + l - 3) as f32) as usize;
        let e2 = e1 + e1;
        let e3 = e2 + e1;
        let e4 = e3 + e1;
        let e5 = e4 + e1;
        let e6 = e5 + e1;
        let e7 = e6 + e1;
        let e8 = e7 + e1;
        for j in (0..n).step_by(e8) {
            let mut t1 = f[j] + f[j + e2];
            let mut t2 = f[j] - f[j + e2];
            let mut t3 = f[j + e4] + f[j + e6];
            let mut t4 = f[j + e4] - f[j + e6];
            f[j] = t1 + t3;
            f[j + e2] = t1 - t3;
            f[j + e4] = t2 + t4;
            f[j + e6] = t2 - t4;
            t1 = f[j + e1];
            t2 = f[j + e3] * r2;
            t3 = f[j + e5];
            t4 = f[j + e7] * r2;
            f[j + e1] = t1 + t2 + t3;
            f[j + e3] = t1 - t3 + t4;
            f[j + e5] = t1 - t2 + t3;
            f[j + e7] = t1 - t3 - t4;
            for k in 1..e1 {
                let l1 = j + k;
                let l2 = l1 + e2;
                let l3 = l1 + e4;
                let l4 = l1 + e6;
                let l5 = j + e2 - k;
                let l6 = l5 + e2;
                let l7 = l5 + e4;
                let l8 = l5 + e6;
                let a1 = (std::f64::consts::PI * k as f64 / e4 as f64) as f32;
                let a2 = a1 + a1;
                let a3 = a1 + a2;
                let c1 = (a1 as f64).cos() as f32;
                let s1 = (a1 as f64).sin() as f32;
                let c2 = (a2 as f64).cos() as f32;
                let s2 = (a2 as f64).sin() as f32;
                let c3 = (a3 as f64).cos() as f32;
                let s3 = (a3 as f64).sin() as f32;
                let t5 = f[l2] * c1 + f[l6] * s1;
                let t6 = f[l3] * c2 + f[l7] * s2;
                let t7 = f[l4] * c3 + f[l8] * s3;
                let t8 = f[l6] * c1 - f[l2] * s1;
                let t9 = f[l7] * c2 - f[l3] * s2;
                let t0 = f[l8] * c3 - f[l4] * s3;
                t1 = f[l5] - t9;
                t2 = f[l5] + t9;
                t3 = -t8 - t0;
                t4 = t5 - t7;
                f[l5] = t1 + t4;
                f[l6] = t2 + t3;
                f[l7] = t1 - t4;
                f[l8] = t2 - t3;
                t1 = f[l1] + t6;
                t2 = f[l1] - t6;
                t3 = t8 - t0;
                t4 = t5 + t7;
                f[l1] = t1 + t4;
                f[l2] = t2 + t3;
                f[l3] = t1 - t4;
                f[l4] = t2 - t3;
            }
        }
    }
}
