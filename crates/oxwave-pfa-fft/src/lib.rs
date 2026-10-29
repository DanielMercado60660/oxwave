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

//! Prime factor FFTs, ported from cwp/lib/pfafft.c.
//! Credits: Dave Hale, Colorado School of Mines, 04/27/89 (complex FFT),
//! 04/28/89 and 06/13/89 (complex size selection), 06/16/89 (real size
//! selection), 06/13/89 (real and 2-D real FFT), and 06/15/89 (multiple
//! and 2-D complex FFT); modified 08/05/91 and 11/26/89.
//! FFT cost table: Dave Hale, 08/04/91; rebuilt by Jack Cohen, 05/03/94,
//! on advice of David and Gregory Chudnovsky.
//! References: C. Temperton, Journal of Computational Physics 58 (1985),
//! 283–299, and 75 (1988), 190–198; Press et al., Numerical Recipes in C
//! (1988), p. 417.
//!
//! Transforms are unnormalized. Complex lengths must be products of mutually
//! prime factors from {2, 3, 4, 5, 7, 8, 9, 11, 13, 16}. The C documentation
//! recommends even real lengths with an admissible half-length; odd lengths
//! retain C's integer division, packing, and twiddle behavior. They are not
//! conventional odd-length DFTs. One-dimensional lengths come from
//! the real data slice, or the complex slice for pfacc. Multidimensional arrays
//! are flat, with the first dimension varying fastest. Input and output slices
//! of real transforms are separate; C's optional pointer aliasing is represented
//! by separate buffers in this safe API.
//! Invalid dimensions, insufficient storage, and unsafe radix combinations
//! terminate with an err! diagnostic instead of a Rust panic.

pub use oxwave_cwp::complex::Complex;

mod kernels;
mod table;

use table::NCTAB;
const PI: f64 = 3.141592653589793;

/// Smallest admissible complex length, saturated at 720720.
/// C prototype: int npfa(int nmin);
pub fn npfa(nmin: i32) -> i32 {
    let mut i = 0;
    while i < NCTAB.len() - 1 && NCTAB[i].0 < nmin {
        i += 1;
    }
    NCTAB[i].0
}

/// Lowest estimated cost among admissible lengths in the requested interval.
/// The upper bound is a preference; see npfa for saturation.
/// C prototype: int npfao(int nmin, int nmax);
pub fn npfao(nmin: i32, nmax: i32) -> i32 {
    let mut i = 0;
    while i < NCTAB.len() - 1 && NCTAB[i].0 < nmin {
        i += 1;
    }
    let mut j = i + 1;
    while j < NCTAB.len() - 1 && NCTAB[j].0 <= nmax {
        if NCTAB[j].1 < NCTAB[i].1 {
            i = j;
        }
        j += 1;
    }
    NCTAB[i].0
}

/// Smallest admissible even real-transform length, saturated at 1441440.
/// C prototype: int npfar(int nmin);
pub fn npfar(nmin: i32) -> i32 {
    // C's nmin+1 overflows at INT_MAX. Widen before rounding so the
    // documented saturation still applies; the half-length fits in i32.
    2 * npfa(((i64::from(nmin) + 1) / 2) as i32)
}

/// Lowest estimated cost for an even real-transform length.
/// Bounds use the same saturation as npfar, including i32::MAX.
/// C prototype: int npfaro(int nmin, int nmax);
pub fn npfaro(nmin: i32, nmax: i32) -> i32 {
    2 * npfao(
        ((i64::from(nmin) + 1) / 2) as i32,
        ((i64::from(nmax) + 1) / 2) as i32,
    )
}

// Complex is repr(C), containing exactly two adjacent f32 fields. Reborrow
// their storage exclusively as float slots, as in the C implementation.
fn floats_mut(z: &mut [Complex]) -> &mut [f32] {
    // SAFETY: Complex has the same alignment and layout as [f32; 2], all f32
    // bit patterns are valid, and the returned lifetime preserves exclusivity.
    unsafe { std::slice::from_raw_parts_mut(z.as_mut_ptr().cast(), 2 * z.len()) }
}

fn floats(z: &[Complex]) -> &[f32] {
    // SAFETY: same layout as above, with shared access for the input lifetime.
    unsafe { std::slice::from_raw_parts(z.as_ptr().cast(), 2 * z.len()) }
}

/// In-place complex FFT; the sign of isign selects the kernel exponent sign.
/// The transform length is z.len().
/// C prototype: void pfacc(int isign, int n, complex z[]);
pub fn pfacc(isign: i32, z: &mut [Complex]) {
    let n = i32::try_from(z.len())
        .unwrap_or_else(|_| oxwave_core::err!("pfafft.c: pfacc: length exceeds C int"));
    if n < 1 {
        oxwave_core::err!("pfafft.c: pfacc: invalid transform length");
    }
    kernels::pfacc(isign, n, floats_mut(z));
}

/// Multiple complex FFTs with length n, count nt, and complex-element strides
/// k (within each transform) and kt (between transforms). Slice gaps survive.
/// A zero count returns immediately without inspecting lengths or strides.
/// Invalid dimensions, negative strides, or insufficient storage exit via err!.
/// C prototype: void pfamcc(int isign, int n, int nt, int k, int kt, complex z[]);
pub fn pfamcc(isign: i32, n: i32, nt: i32, k: i32, kt: i32, z: &mut [Complex]) {
    if nt == 0 {
        return;
    }
    if n < 1 || nt < 0 || k < 0 || kt < 0 {
        oxwave_core::err!("pfafft.c: pfamcc: invalid dimensions or strides");
    }
    // Nonnegative i32 dimensions/strides fit this extent calculation in i64.
    let needed = (i64::from(n) - 1) * i64::from(k) + (i64::from(nt) - 1) * i64::from(kt) + 1;
    if (z.len() as u64) < needed as u64 {
        oxwave_core::err!("pfafft.c: pfamcc: insufficient complex storage");
    }
    kernels::pfamcc(isign, n, nt, k, kt, floats_mut(z));
}

/// Transform dimension idim (1 or 2) of a flat n2-by-n1 complex array.
/// C prototype: void pfa2cc(int isign, int idim, int n1, int n2, complex z[]);
pub fn pfa2cc(isign: i32, idim: i32, n1: i32, n2: i32, z: &mut [Complex]) {
    if n1 < 1 || n2 < 1 {
        oxwave_core::err!("pfafft.c: pfa2cc: invalid dimensions");
    }
    if z.len() < n1 as usize * n2 as usize {
        oxwave_core::err!("pfafft.c: pfa2cc: insufficient complex storage");
    }
    pfa2cc_slots(isign, idim, n1, n2, floats_mut(z));
}

fn pfa2cc_slots(isign: i32, idim: i32, n1: i32, n2: i32, z: &mut [f32]) {
    let (n, nt, k, kt);
    if idim == 1 {
        n = n1;
        nt = n2;
        k = 1;
        kt = n1;
    } else {
        n = n2;
        nt = n1;
        k = n1;
        kt = 1;
    }
    kernels::pfamcc(isign, n, nt, k, kt, z);
}

/// Complex-to-real FFT into rz; n = rz.len(), cz has at least n/2+1 elements.
/// For even n, the DC and Nyquist imaginary components are ignored, as in C.
/// Odd n uses exactly C's float-slot indexing and floor(n/2) complex kernel.
/// At n=1, C's unused write past its documented output is omitted.
/// C prototype: void pfacr(int isign, int n, complex cz[], float rz[]);
pub fn pfacr(isign: i32, cz: &[Complex], rz: &mut [f32]) {
    let n = rz.len();
    if n < 1 || n > i32::MAX as usize {
        oxwave_core::err!("pfafft.c: pfacr: invalid transform length");
    }
    if cz.len() < n / 2 + 1 {
        oxwave_core::err!("pfafft.c: pfacr: insufficient complex storage");
    }
    let z = floats(cz);
    if n == 1 {
        // C also writes rz[1], outside its documented array[n]. That slot is
        // never read again (the length-zero kernel does no work). Preserve
        // the sole in-bounds result without writing beyond the Rust slice.
        rz[0] = z[0] + z[1];
        return;
    }
    // Copy input to output and fix DC and Nyquist.
    for i in 2..n {
        rz[i] = z[i];
    }
    rz[1] = z[0] - z[n];
    rz[0] = z[0] + z[n];
    let z = rz;
    let mut theta = 2.0 * PI / n as f64;
    if isign > 0 {
        theta = -theta;
    }
    let mut wtemp = (0.5 * theta).sin();
    let wpr = -2.0 * wtemp * wtemp;
    let wpi = theta.sin();
    let mut wr = 1.0 + wpr;
    let mut wi = wpi;
    let no2 = n / 2;
    let (mut ir, mut ii, mut jr, mut ji) = (2, 3, n - 2, n - 1);
    while ir <= no2 {
        let sumr = z[ir] + z[jr];
        let sumi = z[ii] + z[ji];
        let difr = z[ir] - z[jr];
        let difi = z[ii] - z[ji];
        let tempr = (wi * difr as f64 - wr * sumi as f64) as f32;
        let tempi = (wi * sumi as f64 + wr * difr as f64) as f32;
        z[ir] = sumr + tempr;
        z[ii] = difi + tempi;
        z[jr] = sumr - tempr;
        z[ji] = tempi - difi;
        wtemp = wr;
        wr += wr * wpr - wi * wpi;
        wi += wi * wpr + wtemp * wpi;
        ir += 2;
        ii += 2;
        jr -= 2;
        ji -= 2;
    }
    kernels::pfacc(isign, (n / 2) as i32, z);
}

/// Real-to-complex FFT; n = rz.len(), cz needs ceil((n+2)/2) elements.
/// For odd n, C writes float slots n and n+1, requiring one extra complex
/// element beyond its documented n/2+1. Odd n keeps C's floor(n/2) kernel.
/// C prototype: void pfarc(int isign, int n, float rz[], complex cz[]);
pub fn pfarc(isign: i32, rz: &[f32], cz: &mut [Complex]) {
    let n = rz.len();
    if n < 1 || n > i32::MAX as usize {
        oxwave_core::err!("pfafft.c: pfarc: invalid transform length");
    }
    if cz.len() < (n + 3) / 2 {
        oxwave_core::err!("pfafft.c: pfarc: insufficient complex storage");
    }
    let z = floats_mut(cz);
    for i in 0..n {
        z[i] = (0.5 * rz[i] as f64) as f32;
    }
    kernels::pfacc(isign, (n / 2) as i32, z);
    z[n] = (2.0 * (z[0] - z[1]) as f64) as f32;
    z[0] = (2.0 * (z[0] + z[1]) as f64) as f32;
    z[n + 1] = 0.0;
    z[1] = 0.0;
    let mut theta = 2.0 * PI / n as f64;
    if isign < 0 {
        theta = -theta;
    }
    let mut wtemp = (0.5 * theta).sin();
    let wpr = -2.0 * wtemp * wtemp;
    let wpi = theta.sin();
    let mut wr = 1.0 + wpr;
    let mut wi = wpi;
    let no2 = n / 2;
    // For n=1 C's jr is -1, but the twiddle loop has no iterations.
    let (mut ir, mut ii, mut jr, mut ji) = (2, 3, n.saturating_sub(2), n - 1);
    while ir <= no2 {
        let sumr = z[ir] + z[jr];
        let sumi = z[ii] + z[ji];
        let difr = z[ir] - z[jr];
        let difi = z[ii] - z[ji];
        let tempr = (wi * difr as f64 + wr * sumi as f64) as f32;
        let tempi = (wi * sumi as f64 - wr * difr as f64) as f32;
        z[ir] = sumr + tempr;
        z[ii] = difi + tempi;
        z[jr] = sumr - tempr;
        z[ji] = tempi - difi;
        wtemp = wr;
        wr += wr * wpr - wi * wpi;
        wi += wi * wpr + wtemp * wpi;
        ir += 2;
        ii += 2;
        jr -= 2;
        ji -= 2;
    }
}

/// Complex-to-real transform of dimension idim of an n2-by-n1 real array.
/// For idim=1, cz is n2 rows of n1/2+1 complex values; for idim=2 it is
/// n2/2+1 rows of n1 complex values. rz has at least n1*n2 real values.
/// For odd transformed lengths, cz needs ceil(((n1+2)*n2-1)/2) elements for
/// idim=1, or ceil((n1*(n2+2)-1)/2) for idim=2. With odd n2 and idim=2, rz
/// needs n1*(n2+1) slots: C unmerges the final pair into the extra row.
/// For n1=1 and idim=1, rz needs n2+1 slots for C's final DC pair.
/// C prototype: void pfa2cr(int isign, int idim, int n1, int n2, complex cz[], float rz[]);
pub fn pfa2cr(isign: i32, idim: i32, n1: i32, n2: i32, cz: &[Complex], rz: &mut [f32]) {
    check_real_dims(idim, n1, n2, rz.len(), cz.len(), true);
    let (n1, n2) = (n1 as usize, n2 as usize);
    let input = floats(cz);
    let (n, nt, itmul, itinc);
    if idim == 1 {
        let (mut jt, mut kt) = (0, 0);
        for _i2 in 0..n2 {
            rz[jt + 1] = input[kt] - input[kt + n1];
            rz[jt] = input[kt] + input[kt + n1];
            let (mut j, mut k) = (jt + 2, kt + 2);
            for _i1 in 2..n1 {
                rz[j] = input[k];
                j += 1;
                k += 1;
            }
            jt += n1;
            kt += n1 + 2;
        }
        n = n1;
        nt = n2;
        itmul = 1;
        itinc = n1;
    } else {
        for i2 in 1..n2 / 2 {
            let mut j = i2 * n1 * 2;
            for _i1 in 0..n1 * 2 {
                rz[j] = input[j];
                j += 1;
            }
        }
        let mut j = n1 * n2;
        for i1 in (0..n1 * 2).step_by(2) {
            rz[i1 + 1] = input[i1] - input[j];
            rz[i1] = input[i1] + input[j];
            j += 2;
        }
        n = n2;
        nt = n1;
        itmul = n1;
        itinc = 2;
    }
    let z = rz;
    let mut theta = 2.0 * PI / n as f64;
    if isign > 0 {
        theta = -theta;
    }
    let mut wtemp = (0.5 * theta).sin();
    let wpr = -2.0 * wtemp * wtemp;
    let wpi = theta.sin();
    let mut wr = 1.0 + wpr;
    let mut wi = wpi;
    // C's unused k is -1 when the transformed length is one.
    let (mut j, mut k) = (2, n.saturating_sub(2));
    while j <= n / 2 {
        let mut jt = j * itmul;
        let mut kt = k * itmul;
        for _it in 0..nt {
            let sumr = z[jt] + z[kt];
            let sumi = z[jt + 1] + z[kt + 1];
            let difr = z[jt] - z[kt];
            let difi = z[jt + 1] - z[kt + 1];
            let tempr = (wi * difr as f64 - wr * sumi as f64) as f32;
            let tempi = (wi * sumi as f64 + wr * difr as f64) as f32;
            z[jt] = sumr + tempr;
            z[jt + 1] = difi + tempi;
            z[kt] = sumr - tempr;
            z[kt + 1] = tempi - difi;
            jt += itinc;
            kt += itinc;
        }
        wtemp = wr;
        wr += wr * wpr - wi * wpi;
        wi += wi * wpr + wtemp * wpi;
        j += 2;
        k -= 2;
    }
    if idim == 1 {
        pfa2cc_slots(isign, 1, (n1 / 2) as i32, n2 as i32, z);
    } else {
        pfa2cc_slots(isign, 2, n1 as i32, (n2 / 2) as i32, z);
        // Unmerge even and odd vectors in the same copy order as C.
        let mut temp = vec![0.0f32; n1];
        for i2 in (0..n2).step_by(2) {
            let mut j = i2 * n1 + 1;
            for i1 in 0..n1 {
                temp[i1] = z[j];
                j += 2;
            }
            let (mut j, mut k) = (i2 * n1, i2 * n1);
            for _i1 in 0..n1 {
                z[k] = z[j];
                j += 2;
                k += 1;
            }
            let mut j = (i2 + 1) * n1;
            for i1 in 0..n1 {
                z[j] = temp[i1];
                j += 1;
            }
        }
    }
}

/// Real-to-complex transform of dimension idim of an n2-by-n1 real array.
/// For idim=1, cz is n2 rows of n1/2+1 complex values; for idim=2 it is
/// n2/2+1 rows of n1 complex values. rz has at least n1*n2 real values.
/// Odd transformed lengths use the C packing: cz needs
/// ceil((n1+2)*n2/2) elements for idim=1, or ceil(n1*(n2+2)/2) for idim=2.
/// Output workspace contents can affect results for odd or unit lengths;
/// the caller's existing values are preserved until C would overwrite them.
/// C prototype: void pfa2rc(int isign, int idim, int n1, int n2, float rz[], complex cz[]);
pub fn pfa2rc(isign: i32, idim: i32, n1: i32, n2: i32, rz: &[f32], cz: &mut [Complex]) {
    check_real_dims(idim, n1, n2, rz.len(), cz.len(), false);
    let (n1, n2) = (n1 as usize, n2 as usize);
    let z = floats_mut(cz);
    let mut j = 0;
    for _i2 in 0..n2 {
        for _i1 in 0..n1 {
            z[j] = (0.5 * rz[j] as f64) as f32;
            j += 1;
        }
    }
    let (n, nt, itmul, itinc);
    if idim == 1 {
        pfa2cc_slots(isign, 1, (n1 / 2) as i32, n2 as i32, z);
        // Shift rows backwards to make room for Nyquist.
        for i2 in (1..n2).rev() {
            let jt = i2 * n1 + n1 - 1;
            let kt = jt + i2 * 2;
            // Signed indices preserve C's final decrement past zero.
            let (mut j, mut k) = (jt as i64, kt as i64);
            for _i1 in (0..n1).rev() {
                z[k as usize] = z[j as usize];
                j -= 1;
                k -= 1;
            }
        }
        n = n1;
        nt = n2;
        itmul = 1;
        itinc = n1 + 2;
    } else {
        // C stores the temporary row in the output's Nyquist space.
        let temp = n1 * n2;
        for i2 in (0..n2).step_by(2) {
            let mut j = i2 * n1;
            for i1 in 0..n1 {
                z[temp + i1] = z[j];
                j += 1;
            }
            let (mut j, mut k) = ((i2 + 1) * n1, i2 * n1 + 1);
            for _i1 in 0..n1 {
                z[k] = z[j];
                j += 1;
                k += 2;
            }
            let mut j = i2 * n1;
            for i1 in 0..n1 {
                z[j] = z[temp + i1];
                j += 2;
            }
        }
        pfa2cc_slots(isign, 2, n1 as i32, (n2 / 2) as i32, z);
        n = n2;
        nt = n1;
        itmul = n1;
        itinc = 2;
    }
    let (mut j, mut k) = (0, n * itmul);
    for _it in 0..nt {
        z[k] = (2.0 * (z[j] - z[j + 1]) as f64) as f32;
        z[j] = (2.0 * (z[j] + z[j + 1]) as f64) as f32;
        z[k + 1] = 0.0;
        z[j + 1] = 0.0;
        j += itinc;
        k += itinc;
    }
    let mut theta = 2.0 * PI / n as f64;
    if isign < 0 {
        theta = -theta;
    }
    let mut wtemp = (0.5 * theta).sin();
    let wpr = -2.0 * wtemp * wtemp;
    let wpi = theta.sin();
    let mut wr = 1.0 + wpr;
    let mut wi = wpi;
    // C's unused k is -1 when the transformed length is one.
    let (mut j, mut k) = (2, n.saturating_sub(2));
    while j <= n / 2 {
        let mut jt = j * itmul;
        let mut kt = k * itmul;
        for _it in 0..nt {
            let sumr = z[jt] + z[kt];
            let sumi = z[jt + 1] + z[kt + 1];
            let difr = z[jt] - z[kt];
            let difi = z[jt + 1] - z[kt + 1];
            let tempr = (wi * difr as f64 + wr * sumi as f64) as f32;
            let tempi = (wi * sumi as f64 - wr * difr as f64) as f32;
            z[jt] = sumr + tempr;
            z[jt + 1] = difi + tempi;
            z[kt] = sumr - tempr;
            z[kt + 1] = tempi - difi;
            jt += itinc;
            kt += itinc;
        }
        wtemp = wr;
        wr += wr * wpr - wi * wpi;
        wi += wi * wpr + wtemp * wpi;
        j += 2;
        k -= 2;
    }
}

fn check_real_dims(idim: i32, n1: i32, n2: i32, nr: usize, nc: usize, to_real: bool) {
    let name = if to_real { "pfa2cr" } else { "pfa2rc" };
    if n1 < 1 || n2 < 1 {
        oxwave_core::err!("pfafft.c: %s: invalid dimensions", name);
    }
    let (n1, n2) = (n1 as usize, n2 as usize);
    // Count the actual float slots C accesses, including its odd-length
    // workspace, rather than imposing the even-length recommendation.
    let (real_needed, complex_slots) = if idim == 1 {
        (n1 * n2 + usize::from(to_real && n1 == 1), (n1 + 2) * n2)
    } else {
        (n1 * (n2 + if to_real { n2 % 2 } else { 0 }), n1 * (n2 + 2))
    };
    let complex_slots = complex_slots - usize::from(to_real);
    if nr < real_needed || nc < (complex_slots + 1) / 2 {
        oxwave_core::err!("pfafft.c: %s: insufficient storage", name);
    }
}
