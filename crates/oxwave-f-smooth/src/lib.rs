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

//! Clean-room implementations for par/lib/smooth.c.
//! Credits: Gabriel Alvarez (1995); Zhenyue Liu (1992, 1994);
//! Balasz Nemeth, Potash Corporation (2008).

#![allow(non_snake_case)]

/// Construct a Savitzky–Golay least-squares filter.
///
/// C prototype: `void SG_smoothing_filter(int np, int nl, int nr, int ld, int m, float *filter);`
///
/// The output slice length is the C `np`. The polynomial must be identifiable
/// from the requested window, and `ld` must not exceed `m`.
pub fn SG_smoothing_filter(nl: usize, nr: usize, ld: usize, m: usize, filter: &mut [f32]) {
    assert!(filter.len() >= nl + nr + 1 && ld <= m && m <= nl + nr);
    let n = m + 1;
    let mut moments = vec![0.0_f32; 2 * m + 1];
    for (power, moment) in moments.iter_mut().enumerate() {
        let mut sum = if power == 0 { 1.0_f32 } else { 0.0_f32 };
        for k in 1..=nr {
            sum = (sum as f64 + (k as f64).powf(power as f64)) as f32;
        }
        for k in 1..=nl {
            sum = (sum as f64 + (-(k as f64)).powf(power as f64)) as f32;
        }
        *moment = sum;
    }
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = moments[i + j];
        }
    }
    let mut b = vec![0.0_f32; n];
    b[ld] = 1.0;
    let mut pivot = vec![0; n];
    let mut determinant = 1.0;
    oxwave_f_axb::LU_decomposition(&mut a, &mut pivot, &mut determinant);
    oxwave_f_axb::backward_substitution(&a, &pivot, &mut b);
    filter.fill(0.0);
    for k in -(nl as i32)..=nr as i32 {
        let mut sum = b[0];
        let mut power = 1.0_f32;
        for &coefficient in &b[1..] {
            power *= k as f32;
            sum += coefficient * power;
        }
        filter[(k + nl as i32) as usize] = sum;
    }
}

/// Construct a rectangular (`flag = 1`) or triangular (`flag = 2`) filter.
///
/// C prototype: `void rwa_smoothing_filter(int flag, int nl, int nr, float *filter);`
///
/// For a triangle with a nonempty right ramp, the C implementation includes
/// the incoming final coefficient in normalization. Initialize it to zero
/// for an ordinary triangular window.
pub fn rwa_smoothing_filter(flag: i32, nl: usize, nr: usize, filter: &mut [f32]) {
    assert!(filter.len() >= nl + nr + 1);
    assert!(flag == 1 || flag == 2);
    let n = nl + nr + 1;
    if flag == 1 {
        filter[..n].fill(1.0 / n as f32);
    } else if n == 1 {
        filter[0] = 1.0;
    } else {
        for i in 0..=nl {
            filter[i] = if nl == 0 {
                1.0
            } else {
                1.0 - (nl - i) as f32 / nl as f32
            };
        }
        // The C right ramp excludes its final endpoint. Its incoming
        // value participates in normalization; callers normally supply zero.
        for i in nl + 1..n - 1 {
            filter[i] = 1.0 - (i - nl) as f32 / nr as f32;
        }
        let mut sum = 0.0_f32;
        for &v in &filter[..n] {
            sum += v;
        }
        for v in &mut filter[..n] {
            *v /= sum;
        }
    }
}

/// Smooth a one-dimensional array with a truncated Gaussian window.
///
/// C prototype: `void gaussian1d_smoothing(int ns, int nsr, float *data);`
///
/// A zero radius leaves the input alone. If `nsr >= ns`, the C routine
/// replaces the samples with their mean instead of applying a convolution.
pub fn gaussian1d_smoothing(nsr: usize, data: &mut [f32]) {
    if nsr == 0 || data.is_empty() {
        return;
    }
    if data.len() <= nsr {
        let mut sum = 0.0_f32;
        for &v in data.iter() {
            sum += v;
        }
        let mean = sum / data.len() as f32;
        data.fill(mean);
        return;
    }
    let length = 3 * nsr + 1;
    let half = length / 2;
    let mut kernel = Vec::with_capacity(length);
    let mut sum = 0.0_f32;
    for i in 0..length {
        let x = (i as i32 - half as i32) as f64 / nsr as f64;
        let v = (-std::f64::consts::PI * x * x).exp() as f32;
        kernel.push(v);
        sum += v;
    }
    for v in &mut kernel {
        *v /= sum;
    }
    let input = data.to_vec();
    for (i, value) in data.iter_mut().enumerate() {
        let mut total = 0.0_f32;
        for (j, &weight) in kernel.iter().enumerate() {
            let k = i as isize + half as isize - j as isize;
            if k >= 0 && (k as usize) < input.len() {
                total += weight * input[k as usize];
            }
        }
        *value = total;
    }
}

/// Smooth a two-dimensional array in both dimensions.
///
/// C prototype: `void gaussian2d_smoothing(int nx, int nt, int nsx, int nst, float **data);`
pub fn gaussian2d_smoothing(nsx: usize, nst: usize, data: &mut [Vec<f32>]) {
    if data.is_empty() {
        return;
    }
    let nt = data[0].len();
    assert!(data.iter().all(|row| row.len() == nt));
    for row in data.iter_mut() {
        gaussian1d_smoothing(nst, row);
    }
    let mut column = vec![0.0; data.len()];
    for t in 0..nt {
        for (x, row) in data.iter().enumerate() {
            column[x] = row[t];
        }
        gaussian1d_smoothing(nsx, &mut column);
        for (x, row) in data.iter_mut().enumerate() {
            row[t] = column[x];
        }
    }
}

/// Smooth and normalize a histogram using a Gaussian kernel.
///
/// C prototype: `void smooth_histogram(int nintlh, float *pdf);`
pub fn smooth_histogram(pdf: &mut [f32]) {
    let mut filter = [0.0_f32; 11];
    for (i, weight) in filter.iter_mut().enumerate() {
        let r = i as f32 - 5.0;
        // C stores the exponent in float before passing it to double exp.
        let exponent = ((-r * r / 1.0 / 1.0) as f64 * 3.141592653589793_f64) as f32;
        *weight = (exponent as f64).exp() as f32;
    }
    let mut smoothed = vec![0.0; pdf.len()];
    oxwave_f_convolution::convolve_cwp(-5, &filter, 0, pdf, 0, &mut smoothed);
    let mut sum = 0.0_f32;
    for &value in &smoothed {
        sum += value;
    }
    for (out, value) in pdf.iter_mut().zip(smoothed) {
        *out = value / sum;
    }
}

/// Damped least-squares smoothing on nonuniform horizontal coordinates.
///
/// C prototype: `void smooth_1(float *x, float *z, float r, int n);`
///
/// Uses two implicit diffusion steps, each solving the symmetric
/// tridiagonal normal equations with natural boundary conditions.
/// Adjacent horizontal coordinates should be distinct.
pub fn smooth_1(x: &[f32], z: &mut [f32], r: f32) {
    assert_eq!(x.len(), z.len());
    let n = z.len();
    if n == 0 {
        return;
    }
    // Both C assignments round to float, but their unsuffixed constants
    // make the multiplication and division use double intermediates.
    let scaled = (0.5 * r as f64 * r as f64) as f32;
    let strength = (scaled as f64 / (5.19 * 5.19)) as f32;
    let mut inverse_squared = vec![0.0; n - 1];
    for i in 0..n - 1 {
        let inverse = 1.0 / (x[i + 1] - x[i]);
        inverse_squared[i] = inverse * inverse;
    }
    let mut diagonal = vec![0.0; n];
    diagonal[0] = 1.0 + strength * inverse_squared.first().copied().unwrap_or(0.0);
    for i in 1..n {
        let left = inverse_squared[i - 1];
        let right = if i + 1 < n { inverse_squared[i] } else { 0.0 };
        diagonal[i] = 1.0 + strength * (left + right);
    }
    let edge: Vec<f32> = inverse_squared.iter().map(|&v| -strength * v).collect();
    for _ in 0..2 {
        oxwave_tridiag::tripd(&mut diagonal.clone(), &mut edge.clone(), z);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn singleton_smoothing_still_solves_with_nan_strength() {
        let mut z = [2.0];
        smooth_1(&[0.0], &mut z, f32::NAN);
        assert!(z[0].is_nan());
    }

    #[test]
    fn large_strength_is_scaled_before_float_overflow() {
        let mut z = [2.0, 3.0, 4.0, 5.0];
        smooth_1(&[0.0, 0.3, 0.6, 0.9], &mut z, 2e19);
        assert!(z.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn two_point_large_strength_matches_c_infinities() {
        let mut z = [2.0, 3.0];
        smooth_1(&[0.0, 1.0], &mut z, 2e19);
        // The finite coefficient leads to a singular solve in float.
        // Premature overflow of r*r would produce NaNs instead.
        assert_eq!(z, [f32::INFINITY; 2]);
    }

    #[test]
    fn histogram_matches_c_with_fewer_bins_than_filter() {
        let mut pdf = [0.125, 0.625, 0.25];
        smooth_histogram(&mut pdf);
        // Recorded from this Mac's original smooth_histogram. Convolution
        // swaps the inputs when the histogram is shorter than the filter.
        assert_eq!(pdf.map(f32::to_bits), [0x3e1171b4, 0x3f1960b4, 0x3e8485be]);
    }

    #[test]
    fn histogram_matches_c_with_more_bins_than_filter() {
        let mut pdf = [0.0; 23];
        for (i, value) in pdf.iter_mut().enumerate() {
            *value = ((i * 17) % 31 + 1) as f32 / 256.0;
        }
        smooth_histogram(&mut pdf);
        // Original C output bits fix the accumulation order at the edges
        // and in the middle, as well as the float-rounded Gaussian weights.
        assert_eq!(
            pdf.map(f32::to_bits),
            [
                0x3b96b111, 0x3d40fff9, 0x3c70f220, 0x3d638861, 0x3cbd899a, 0x3d830853, 0x3d014d12,
                0x3d944c75, 0x3d23d557, 0x3da59072, 0x3d382c1a, 0x3c14dd40, 0x3d4c82dd, 0x3c8f7e94,
                0x3d6f0b23, 0x3cd48f1d, 0x3d88c9b3, 0x3d0ccfd3, 0x3d9a0dd6, 0x3d2f5818, 0x3dab51d3,
                0x3d43aead, 0x3c1e4741,
            ]
        );
    }

    #[test]
    fn histogram_uses_rounded_exponent_and_finite_filter() {
        let mut pdf = [0.0; 13];
        pdf[0] = 1.0;
        smooth_histogram(&mut pdf);
        let exponent = (-std::f64::consts::PI as f32) as f64;
        let weight = exponent.exp() as f32;
        let mut sum = 0.0_f32;
        for i in 0..=5 {
            let e = (-(i * i) as f64 * 3.141592653589793_f64) as f32;
            sum += (e as f64).exp() as f32;
        }
        assert_eq!(pdf[1].to_bits(), (weight / sum).to_bits());
        assert_eq!(pdf[6], 0.0);
    }
}

/// Smooth a rectangular data window by damped least squares.
///
/// C prototype: `void dlsq_smoothing(int nt, int nx, int ift, int ilt,
/// int ifx, int ilx, float r1, float r2, float rw, float **traces);`
///
/// Bounds are half-open. Rows are traces and columns are time samples.
/// The time pass keeps sample zero fixed; the spatial pass has natural
/// boundaries. A nonzero `rw` diffuses the window before smoothing the data.
/// Writes the cumulative relative RMS error through each time sample as
/// native `float`s to `error`, as does the C library. With a singleton time axis, Rust leaves
/// the time pass alone; the C time pass accesses outside its arrays.
pub fn dlsq_smoothing(
    ift: usize,
    ilt: usize,
    ifx: usize,
    ilx: usize,
    r1: f32,
    r2: f32,
    rw: f32,
    traces: &mut [Vec<f32>],
) {
    let nx = traces.len();
    if nx == 0 {
        return;
    }
    let nt = traces[0].len();
    assert!(traces.iter().all(|row| row.len() == nt));
    assert!(ift <= ilt && ilt <= nt && ifx <= ilx && ilx <= nx);
    if nt == 0 {
        return;
    }
    let original = traces.to_vec();
    let mut window = vec![vec![0.0; nt]; nx];
    for row in &mut window[ifx..ilx] {
        row[ift..ilt].fill(1.0);
    }
    let w_strength = rw * rw * 0.25;
    if rw != 0.0 {
        let mut column = vec![0.0; nx];
        for t in 0..nt {
            for x in 0..nx {
                column[x] = window[x][t];
            }
            weighted_diffusion(&mut column, &vec![1.0; nx], w_strength, false);
            for x in 0..nx {
                window[x][t] = column[x];
            }
        }
        for row in &mut window {
            weighted_diffusion(row, &vec![1.0; nt], w_strength, false);
        }
    }
    let mut column = vec![0.0; nx];
    let mut weights = vec![0.0; nx];
    for t in 0..nt {
        for x in 0..nx {
            column[x] = traces[x][t];
            weights[x] = window[x][t];
        }
        if nx == 1 {
            column[0] /= 1.0 + r2 * r2 * 0.25 * weights[0];
        } else {
            weighted_diffusion(&mut column, &weights, r2 * r2 * 0.25, false);
        }
        for x in 0..nx {
            traces[x][t] = column[x];
        }
    }
    for x in 0..nx {
        weighted_diffusion(&mut traces[x], &window[x], r1 * r1 * 0.25, true);
    }
    let mut errors = vec![0.0_f32; nt];
    let mut energy = 0.0_f32;
    let mut residual = 0.0_f32;
    // Both sums carry forward across time, matching the observable report.
    for t in 0..nt {
        for x in 0..nx {
            energy += original[x][t] * original[x][t];
            let delta = traces[x][t] - original[x][t];
            residual += delta * delta;
        }
        errors[t] = if (energy as f64) <= 0.001 {
            1.0
        } else {
            ((residual / energy) as f64).sqrt() as f32
        };
    }
    let file = oxwave_core::prelude::efopen("error", "w");
    oxwave_core::prelude::efwrite(&errors, &file);
    oxwave_core::prelude::efclose(&file);
}

fn weighted_diffusion(z: &mut [f32], weight: &[f32], strength: f32, fixed_first: bool) {
    let n = z.len();
    if n < 2 {
        return;
    }
    let start = usize::from(fixed_first);
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n - 1];
    for i in start..n {
        let left = if i > 0 { weight[i] } else { 0.0 };
        let right = if i + 1 < n { weight[i + 1] } else { 0.0 };
        d[i] = 1.0 + strength * (left + right);
    }
    for i in 0..n - 1 {
        e[i] = -strength * weight[i + 1];
    }
    // The C two-sample time solve omits the fixed sample's RHS contribution.
    if fixed_first && n > 2 {
        z[1] -= e[0] * z[0];
    }
    oxwave_tridiag::tripd(&mut d[start..], &mut e[start..], &mut z[start..]);
}

/// Smooth each contiguous run of sample indices independently.
///
/// C prototype: `void smooth_segmented_array(float *index, float *val,
/// int n, int sm, int inc, int m);`
///
/// Methods: 1 Gaussian, 2 damped least squares, 3 degree-four
/// Savitzky–Golay, 4 rectangular average, 5 triangular average, 6 median.
/// A jump larger than `inc` begins a new segment. As in the C interface,
/// a final singleton is left alone. Method 3 requires `sm >= 2`.
/// The C method-2 dispatch crashes on nontrivial segments; Rust uses the
/// documented nonuniform least-squares smoother for those segments.
/// Method 6 preserves C's selection by absolute magnitude at rank
/// `sm/2 + 1` in a window of length `2*(sm/2)+1`, rather than changing the
/// observable selection to a conventional median. For `sm < 2`, Rust leaves
/// method 6 alone instead of reading past a one-element window as C does.
pub fn smooth_segmented_array(index: &[f32], val: &mut [f32], sm: usize, inc: i32, m: i32) {
    assert_eq!(index.len(), val.len());
    assert!((1..=6).contains(&m));
    assert!(m != 3 || sm >= 2);
    let n = val.len();
    let mut start = 0;
    for i in 1..n {
        let gap = index[i] - index[i - 1] > inc as f32;
        if gap || i == n - 1 {
            let end = if gap { i } else { i + 1 };
            segmented_run(&index[start..end], &mut val[start..end], sm, m);
            start = i;
        }
    }
}

fn segmented_run(index: &[f32], val: &mut [f32], sm: usize, m: i32) {
    match m {
        1 => gaussian1d_smoothing(sm, val),
        2 => smooth_1(index, val, sm as f32),
        3..=5 => {
            let mut kernel = vec![0.0; 2 * sm + 1];
            if m == 3 {
                SG_smoothing_filter(sm, sm, 0, 4, &mut kernel);
            } else {
                rwa_smoothing_filter(m - 3, sm, sm, &mut kernel);
            }
            let input = val.to_vec();
            oxwave_f_convolution::convolve_cwp(-(sm as i32), &kernel, 0, &input, 0, val);
        }
        6 => {
            let half = sm / 2;
            if half == 0 || val.len() < 2 * half + 1 {
                return;
            }
            let input = val.to_vec();
            for i in half..val.len() - half {
                let mut window = input[i - half..=i + half].to_vec();
                window.sort_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap());
                val[i] = window[half + 1];
            }
        }
        _ => unreachable!(),
    }
}
