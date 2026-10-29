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

//! Port of cwp/lib/dpfafft.c.
//! Author: Dave Hale, Colorado School of Mines, 04/27/89.
//! Revised by Baoniu Han to handle double precision, 12/14/98.
//! Length selection and real, multiple, and 2-D transforms: Dave Hale,
//! 04/28/89–06/16/89; modified 08/05/91 and 11/26/89.
//! FFT cost table: Dave Hale, 08/04/91; rebuilt by Jack Cohen, 05/03/94,
//! on advice of David and Gregory Chudnovsky.
//!
//! Valid complex lengths are products of mutually prime factors from
//! {2, 3, 4, 5, 7, 8, 9, 11, 13, 16}, up to 720720. Use `npfa_d` or
//! `npfao_d` to select them. Real transforms require even lengths whose
//! half is a valid complex length; use `npfar_d` or `npfaro_d`.
//! The transforms are unnormalized, just as in C.

const NCTAB: [(i32, f64); 240] = [
    (1, 0.000052),
    (2, 0.000061),
    (3, 0.000030),
    (4, 0.000053),
    (5, 0.000066),
    (6, 0.000067),
    (7, 0.000071),
    (8, 0.000062),
    (9, 0.000079),
    (10, 0.000080),
    (11, 0.000052),
    (12, 0.000069),
    (13, 0.000103),
    (14, 0.000123),
    (15, 0.000050),
    (16, 0.000086),
    (18, 0.000108),
    (20, 0.000101),
    (21, 0.000098),
    (22, 0.000135),
    (24, 0.000090),
    (26, 0.000165),
    (28, 0.000084),
    (30, 0.000132),
    (33, 0.000158),
    (35, 0.000138),
    (36, 0.000147),
    (39, 0.000207),
    (40, 0.000156),
    (42, 0.000158),
    (44, 0.000176),
    (45, 0.000171),
    (48, 0.000185),
    (52, 0.000227),
    (55, 0.000242),
    (56, 0.000194),
    (60, 0.000215),
    (63, 0.000233),
    (65, 0.000288),
    (66, 0.000271),
    (70, 0.000248),
    (72, 0.000247),
    (77, 0.000285),
    (78, 0.000395),
    (80, 0.000285),
    (84, 0.000209),
    (88, 0.000332),
    (90, 0.000321),
    (91, 0.000372),
    (99, 0.000400),
    (104, 0.000391),
    (105, 0.000358),
    (110, 0.000440),
    (112, 0.000367),
    (117, 0.000494),
    (120, 0.000413),
    (126, 0.000424),
    (130, 0.000549),
    (132, 0.000480),
    (140, 0.000450),
    (143, 0.000637),
    (144, 0.000497),
    (154, 0.000590),
    (156, 0.000626),
    (165, 0.000654),
    (168, 0.000536),
    (176, 0.000656),
    (180, 0.000611),
    (182, 0.000730),
    (195, 0.000839),
    (198, 0.000786),
    (208, 0.000835),
    (210, 0.000751),
    (220, 0.000826),
    (231, 0.000926),
    (234, 0.000991),
    (240, 0.000852),
    (252, 0.000820),
    (260, 0.001053),
    (264, 0.000987),
    (273, 0.001152),
    (280, 0.000952),
    (286, 0.001299),
    (308, 0.001155),
    (312, 0.001270),
    (315, 0.001156),
    (330, 0.001397),
    (336, 0.001173),
    (360, 0.001259),
    (364, 0.001471),
    (385, 0.001569),
    (390, 0.001767),
    (396, 0.001552),
    (420, 0.001516),
    (429, 0.002015),
    (440, 0.001748),
    (455, 0.001988),
    (462, 0.001921),
    (468, 0.001956),
    (495, 0.002106),
    (504, 0.001769),
    (520, 0.002196),
    (528, 0.002127),
    (546, 0.002454),
    (560, 0.002099),
    (572, 0.002632),
    (585, 0.002665),
    (616, 0.002397),
    (624, 0.002711),
    (630, 0.002496),
    (660, 0.002812),
    (693, 0.002949),
    (715, 0.003571),
    (720, 0.002783),
    (728, 0.003060),
    (770, 0.003392),
    (780, 0.003553),
    (792, 0.003198),
    (819, 0.003726),
    (840, 0.003234),
    (858, 0.004354),
    (880, 0.003800),
    (910, 0.004304),
    (924, 0.003975),
    (936, 0.004123),
    (990, 0.004517),
    (1001, 0.005066),
    (1008, 0.003902),
    (1040, 0.004785),
    (1092, 0.005017),
    (1144, 0.005599),
    (1155, 0.005380),
    (1170, 0.005730),
    (1232, 0.005323),
    (1260, 0.005112),
    (1287, 0.006658),
    (1320, 0.005974),
    (1365, 0.006781),
    (1386, 0.006413),
    (1430, 0.007622),
    (1456, 0.006679),
    (1540, 0.007032),
    (1560, 0.007538),
    (1584, 0.007126),
    (1638, 0.007979),
    (1680, 0.007225),
    (1716, 0.008961),
    (1820, 0.008818),
    (1848, 0.008427),
    (1872, 0.009004),
    (1980, 0.009398),
    (2002, 0.010830),
    (2145, 0.012010),
    (2184, 0.010586),
    (2288, 0.012058),
    (2310, 0.011673),
    (2340, 0.011700),
    (2520, 0.011062),
    (2574, 0.014313),
    (2640, 0.013021),
    (2730, 0.014606),
    (2772, 0.013216),
    (2860, 0.015789),
    (3003, 0.016988),
    (3080, 0.014911),
    (3120, 0.016393),
    (3276, 0.016741),
    (3432, 0.018821),
    (3465, 0.018138),
    (3640, 0.018892),
    (3696, 0.018634),
    (3960, 0.020216),
    (4004, 0.022455),
    (4095, 0.022523),
    (4290, 0.026087),
    (4368, 0.023474),
    (4620, 0.024590),
    (4680, 0.025641),
    (5005, 0.030303),
    (5040, 0.025253),
    (5148, 0.030364),
    (5460, 0.031250),
    (5544, 0.029412),
    (5720, 0.034404),
    (6006, 0.037500),
    (6160, 0.034091),
    (6435, 0.040214),
    (6552, 0.037221),
    (6864, 0.042735),
    (6930, 0.040214),
    (7280, 0.042980),
    (7920, 0.045872),
    (8008, 0.049505),
    (8190, 0.049834),
    (8580, 0.055762),
    (9009, 0.057034),
    (9240, 0.054945),
    (9360, 0.056818),
    (10010, 0.066667),
    (10296, 0.065502),
    (10920, 0.068182),
    (11088, 0.065217),
    (11440, 0.075000),
    (12012, 0.078534),
    (12870, 0.087719),
    (13104, 0.081081),
    (13860, 0.084270),
    (15015, 0.102740),
    (16016, 0.106383),
    (16380, 0.105634),
    (17160, 0.119048),
    (18018, 0.123967),
    (18480, 0.119048),
    (20020, 0.137615),
    (20592, 0.140187),
    (21840, 0.154639),
    (24024, 0.168539),
    (25740, 0.180723),
    (27720, 0.180723),
    (30030, 0.220588),
    (32760, 0.241935),
    (34320, 0.254237),
    (36036, 0.254237),
    (40040, 0.288462),
    (45045, 0.357143),
    (48048, 0.357143),
    (51480, 0.384615),
    (55440, 0.384615),
    (60060, 0.454545),
    (65520, 0.517241),
    (72072, 0.576923),
    (80080, 0.625000),
    (90090, 0.833333),
    (102960, 0.789474),
    (120120, 1.153846),
    (144144, 1.153846),
    (180180, 1.875000),
    (240240, 2.500000),
    (360360, 3.750000),
    (720720, 7.500000),
];

/// Return the smallest valid prime-factor transform length at least `nmin`,
/// or 720720 when `nmin` exceeds the largest supported length.
/// C prototype: `int npfa_d(int nmin);`
pub fn npfa_d(nmin: i32) -> i32 {
    let mut i = 0;
    while i < 239 && NCTAB[i].0 < nmin {
        i += 1;
    }
    NCTAB[i].0
}

/// Return the least-cost valid transform length, treating `nmax` as a
/// desired upper bound. As in C, return 720720 if `nmin` exceeds it.
/// C prototype: `int npfao_d(int nmin, int nmax);`
pub fn npfao_d(nmin: i32, nmax: i32) -> i32 {
    let mut i = 0;
    while i < 239 && NCTAB[i].0 < nmin {
        i += 1;
    }
    let mut j = i + 1;
    while j < 239 && NCTAB[j].0 <= nmax {
        if NCTAB[j].1 < NCTAB[i].1 {
            i = j;
        }
        j += 1;
    }
    NCTAB[i].0
}

/// Return the smallest valid real transform length at least `nmin`.
/// C prototype: `int npfar_d(int nmin);`
pub fn npfar_d(nmin: i32) -> i32 {
    2 * npfa_d(nmin.wrapping_add(1) / 2)
}

/// Return the least-cost valid real transform length in the requested range.
/// C prototype: `int npfaro_d(int nmin, int nmax);`
pub fn npfaro_d(nmin: i32, nmax: i32) -> i32 {
    2 * npfao_d(nmin.wrapping_add(1) / 2, nmax.wrapping_add(1) / 2)
}

use oxwave_cwp::complex::DComplex;
use std::f64::consts::PI;
mod kernels;

// DComplex is repr(C), containing exactly two adjacent f64 fields. Reborrow
// their storage exclusively as double slots, as in the C implementation.
fn floats_mut(z: &mut [DComplex]) -> &mut [f64] {
    // SAFETY: DComplex has the same alignment and layout as [f64; 2], all f64
    // bit patterns are valid, and the returned lifetime preserves exclusivity.
    unsafe { std::slice::from_raw_parts_mut(z.as_mut_ptr().cast(), 2 * z.len()) }
}

fn floats(z: &[DComplex]) -> &[f64] {
    // SAFETY: same layout as above, with shared access for the input lifetime.
    unsafe { std::slice::from_raw_parts(z.as_ptr().cast(), 2 * z.len()) }
}

/// In-place complex FFT; the sign of isign selects the kernel exponent sign.
/// The transform length is z.len() and must be valid (see the module docs).
/// C prototype: void pfacc_d(int isign, int n, real_complex z[]);
pub fn pfacc_d(isign: i32, z: &mut [DComplex]) {
    let n = i32::try_from(z.len()).expect("FFT length exceeds C int");
    assert!(n >= 1);
    kernels::pfacc(isign, n, floats_mut(z));
}

/// Multiple complex FFTs with length n, count nt, and complex-element strides
/// k (within each transform) and kt (between transforms). Slice gaps survive.
/// C prototype: void pfamcc_d(int isign, int n, int nt, int k, int kt, real_complex z[]);
pub fn pfamcc_d(isign: i32, n: i32, nt: i32, k: i32, kt: i32, z: &mut [DComplex]) {
    assert!(n >= 1 && nt >= 0 && k >= 0 && kt >= 0);
    let needed = if nt == 0 {
        0
    } else {
        (n as usize - 1) * k as usize + (nt as usize - 1) * kt as usize + 1
    };
    assert!(z.len() >= needed);
    kernels::pfamcc(isign, n, nt, k, kt, floats_mut(z));
}

/// Transform dimension idim (1 or 2) of a flat n2-by-n1 complex array.
/// C prototype: void pfa2cc_d(int isign, int idim, int n1, int n2, real_complex z[]);
pub fn pfa2cc_d(isign: i32, idim: i32, n1: i32, n2: i32, z: &mut [DComplex]) {
    assert!(n1 >= 1 && n2 >= 1);
    assert!(z.len() >= n1 as usize * n2 as usize);
    pfa2cc_slots(isign, idim, n1, n2, floats_mut(z));
}

fn pfa2cc_slots(isign: i32, idim: i32, n1: i32, n2: i32, z: &mut [f64]) {
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

/// DComplex-to-real FFT into rz; n = rz.len(), cz has at least n/2+1 elements.
/// For even n, the DC and Nyquist imaginary components are ignored, as in C.
/// Odd n uses exactly C's double-slot indexing and floor(n/2) complex kernel.
/// C prototype: void pfacr_d(int isign, int n, real_complex cz[], real rz[]);
pub fn pfacr_d(isign: i32, cz: &[DComplex], rz: &mut [f64]) {
    let n = rz.len();
    assert!(n >= 1 && cz.len() >= n / 2 + 1);
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
        let tempr = (wi * difr as f64 - wr * sumi as f64) as f64;
        let tempi = (wi * sumi as f64 + wr * difr as f64) as f64;
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
/// For odd n, C writes double slots n and n+1, requiring one extra complex
/// element beyond its documented n/2+1. Odd n keeps C's floor(n/2) kernel.
/// C prototype: void pfarc_d(int isign, int n, real rz[], real_complex cz[]);
pub fn pfarc_d(isign: i32, rz: &[f64], cz: &mut [DComplex]) {
    let n = rz.len();
    assert!(n >= 1 && cz.len() >= (n + 3) / 2);
    let z = floats_mut(cz);
    for i in 0..n {
        z[i] = (0.5 * rz[i] as f64) as f64;
    }
    kernels::pfacc(isign, (n / 2) as i32, z);
    z[n] = (2.0 * (z[0] - z[1]) as f64) as f64;
    z[0] = (2.0 * (z[0] + z[1]) as f64) as f64;
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
        let tempr = (wi * difr as f64 + wr * sumi as f64) as f64;
        let tempi = (wi * sumi as f64 - wr * difr as f64) as f64;
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

/// DComplex-to-real transform of dimension idim of an n2-by-n1 real array.
/// For idim=1, cz is n2 rows of n1/2+1 complex values; for idim=2 it is
/// n2/2+1 rows of n1 complex values. rz has at least n1*n2 real values.
/// For odd transformed lengths, cz needs ceil(((n1+2)*n2-1)/2) elements for
/// idim=1, or ceil((n1*(n2+2)-1)/2) for idim=2. With odd n2 and idim=2, rz
/// needs n1*(n2+1) slots: C unmerges the final pair into the extra row.
/// For n1=1 and idim=1, rz needs n2+1 slots for C's final DC pair.
/// C prototype: void pfa2cr_d(int isign, int idim, int n1, int n2, real_complex cz[], real rz[]);
pub fn pfa2cr_d(isign: i32, idim: i32, n1: i32, n2: i32, cz: &[DComplex], rz: &mut [f64]) {
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
            let tempr = (wi * difr as f64 - wr * sumi as f64) as f64;
            let tempi = (wi * sumi as f64 + wr * difr as f64) as f64;
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
        let mut temp = vec![0.0f64; n1];
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
/// C prototype: void pfa2rc_d(int isign, int idim, int n1, int n2, real rz[], real_complex cz[]);
pub fn pfa2rc_d(isign: i32, idim: i32, n1: i32, n2: i32, rz: &[f64], cz: &mut [DComplex]) {
    check_real_dims(idim, n1, n2, rz.len(), cz.len(), false);
    let (n1, n2) = (n1 as usize, n2 as usize);
    let z = floats_mut(cz);
    let mut j = 0;
    for _i2 in 0..n2 {
        for _i1 in 0..n1 {
            z[j] = (0.5 * rz[j] as f64) as f64;
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
            let (mut j, mut k) = (jt as i32, kt as i32);
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
        z[k] = (2.0 * (z[j] - z[j + 1]) as f64) as f64;
        z[j] = (2.0 * (z[j] + z[j + 1]) as f64) as f64;
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
            let tempr = (wi * difr as f64 + wr * sumi as f64) as f64;
            let tempi = (wi * sumi as f64 - wr * difr as f64) as f64;
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
    assert!(n1 >= 1 && n2 >= 1);
    let (n1, n2) = (n1 as usize, n2 as usize);
    // Count the actual double slots C accesses, including its odd-length
    // workspace, rather than imposing the even-length recommendation.
    let (real_needed, complex_slots) = if idim == 1 {
        (n1 * n2 + usize::from(to_real && n1 == 1), (n1 + 2) * n2)
    } else {
        (n1 * (n2 + if to_real { n2 % 2 } else { 0 }), n1 * (n2 + 2))
    };
    let complex_slots = complex_slots - usize::from(to_real);
    assert!(nr >= real_needed && nc >= (complex_slots + 1) / 2);
}
