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

//! Port of cwp/lib/bigmatrix.c.
//! Author: Dave Hale, Colorado School of Mines, 05/17/89.

/// Opaque big matrix state, holding elements as native bytes in C's page layout.
pub struct BigMatrix {
    l1: usize,
    l2: usize,
    m1: usize,
    m2: usize,
    nbpe: usize,
    n1: usize,
    n2: usize,
    p: Vec<u8>,
}

/// Allocate a paged matrix; return `None` for invalid dimensions or allocation failure.
///
/// Dimensions and element size must be positive, and the padded allocation must
/// fit in a C `int`. Storage is zeroed for Rust safety; initialize elements before
/// comparing reads with C, whose `malloc` storage is uninitialized.
///
/// C prototype: `void *bmalloc(int nbpe, int n1, int n2);`
pub fn bmalloc(nbpe: i32, n1: i32, n2: i32) -> Option<Box<BigMatrix>> {
    if nbpe <= 0 || n1 <= 0 || n2 <= 0 {
        return None;
    }
    // Keep the C int arithmetic, including the target's wrapping abs(INT_MIN).
    let mut wmin = n2.wrapping_sub(1024i32.wrapping_mul(n1)).wrapping_abs();
    let (nb, n1, n2) = (nbpe as usize, n1 as usize, n2 as usize);
    let (mut l1, mut l2) = (1usize, 1024usize);
    let (mut l1t, mut l2t) = (2usize, 512usize);
    while l1t <= 1024 {
        let w = (l1t as i32)
            .wrapping_mul(n2 as i32)
            .wrapping_sub((l2t as i32).wrapping_mul(n1 as i32))
            .wrapping_abs();
        if w < wmin {
            wmin = w;
            l1 = l1t;
            l2 = l2t;
        }
        l1t *= 2;
        l2t /= 2;
    }
    let m1 = 1 + (n1 - 1) / l1;
    let m2 = 1 + (n2 - 1) / l2;
    let size = nb
        .checked_mul(l1)?
        .checked_mul(l2)?
        .checked_mul(m1)?
        .checked_mul(m2)?;
    // The C allocation size is evaluated in signed int arithmetic.
    if size > i32::MAX as usize {
        return None;
    }
    let mut p = Vec::new();
    p.try_reserve_exact(size).ok()?;
    p.resize(size, 0);
    Some(Box::new(BigMatrix {
        l1,
        l2,
        m1,
        m2,
        nbpe: nb,
        n1,
        n2,
        p,
    }))
}

/// Print the C-style diagnostic state, including the current allocation addresses.
///
/// C prototype: `void bmdump(void *bm);`
pub fn bmdump(bm: &BigMatrix) {
    oxwave_core::printf!("big matrix state %lu:\n", bm as *const BigMatrix as u64);
    oxwave_core::printf!("l1 = %d\n", bm.l1 as i32);
    oxwave_core::printf!("l2 = %d\n", bm.l2 as i32);
    oxwave_core::printf!("m1 = %d\n", bm.m1 as i32);
    oxwave_core::printf!("m2 = %d\n", bm.m2 as i32);
    oxwave_core::printf!("nbpe = %d\n", bm.nbpe as i32);
    oxwave_core::printf!("p = %lu\n", bm.p.as_ptr() as u64);
}

/// Release the matrix and its storage.
///
/// C prototype: `void bmfree(void *bm);`
pub fn bmfree(bm: Box<BigMatrix>) {
    drop(bm);
}

impl BigMatrix {
    fn index(&self, d: i32, k1: usize, k2: usize, n: usize) -> (usize, usize, usize, usize, usize) {
        assert!(k1 < self.n1 && k2 < self.n2);
        if d == 1 {
            assert!(n <= self.n1 - k1);
        } else {
            assert!(n <= self.n2 - k2);
        }
        let (i1, i2) = (k1 % self.l1, k2 % self.l2);
        let (j01, j2) = (k1 / self.l1, k2 / self.l2);
        if d == 1 {
            let p = self.nbpe
                * (i2 * self.l1 + j01 * self.l1 * self.l2 + j2 * self.l1 * self.l2 * self.m1);
            (
                p + i1 * self.nbpe,
                p,
                self.nbpe,
                self.nbpe * self.l1 * self.l2,
                i1,
            )
        } else {
            let p = self.nbpe * (i1 + j01 * self.l1 * self.l2 + j2 * self.l1 * self.l2 * self.m1);
            (
                p + i2 * self.nbpe * self.l1,
                p,
                self.nbpe * self.l1,
                self.nbpe * self.l1 * self.l2 * self.m1,
                i2,
            )
        }
    }
}

/// Read elements into a byte slice whose length is a multiple of the element size.
/// Any direction other than 1 selects the second dimension, as in C.
///
/// C prototype: `void bmread(void *bm, int d, int k1, int k2, int n, void *v);`
pub fn bmread(bm: &BigMatrix, d: i32, k1: usize, k2: usize, v: &mut [u8]) {
    assert_eq!(v.len() % bm.nbpe, 0);
    let n = v.len() / bm.nbpe;
    let (mut q, mut p, qstep, pstep, mut i) = bm.index(d, k1, k2, n);
    let l = if d == 1 { bm.l1 } else { bm.l2 };
    for k in 0..n {
        // C's char/short/long/float/double branches all copy these exact bytes.
        // Keeping byte slices also permits arbitrary element sizes and alignment.
        v[k * bm.nbpe..(k + 1) * bm.nbpe].copy_from_slice(&bm.p[q..q + bm.nbpe]);
        i += 1;
        if i < l {
            q += qstep;
        } else {
            i = 0;
            p += pstep;
            q = p;
        }
    }
}

/// Write elements from a byte slice whose length is a multiple of the element size.
/// Any direction other than 1 selects the second dimension, as in C.
///
/// C prototype: `void bmwrite(void *bm, int d, int k1, int k2, int n, void *v);`
pub fn bmwrite(bm: &mut BigMatrix, d: i32, k1: usize, k2: usize, v: &[u8]) {
    assert_eq!(v.len() % bm.nbpe, 0);
    let n = v.len() / bm.nbpe;
    let (mut q, mut p, qstep, pstep, mut i) = bm.index(d, k1, k2, n);
    let l = if d == 1 { bm.l1 } else { bm.l2 };
    for k in 0..n {
        bm.p[q..q + bm.nbpe].copy_from_slice(&v[k * bm.nbpe..(k + 1) * bm.nbpe]);
        i += 1;
        if i < l {
            q += qstep;
        } else {
            i = 0;
            p += pstep;
            q = p;
        }
    }
}
