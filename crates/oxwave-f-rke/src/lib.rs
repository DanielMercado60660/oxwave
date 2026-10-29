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

//! Port of `par/lib/rke.c`: ordinary differential equation solver,
//! Runge-Kutta-England technique.
//!
//! Credits: Made available for SU by the author's permission, 1 Sept 2006.
//! Copyright 1988: François Pinard <pinard@iro.umontreal.ca>, 1988.
//! Modified: John Stockwell, Center for Wave Phenomena (CWP), Colorado School
//! of Mines (CSM), August 2004, for integration into the CWP/SU Seismic Unix
//! package. Some stylistic changes were made, and some SU routines were
//! incorporated into the modified code.

/// Reentrancy block and user-adjustable integration controls from `rke.h`.
pub struct RkeVariables<F> {
    pub n_equations: usize,
    pub eval_routine: F,
    pub minimum_step: f64,
    pub maximum_step: f64,
    pub current_step: f64,
    pub error_slope: f64,
    pub error_bias: f64,
    pub accepted_steps: i32,
    pub rejected_steps: i32,
}

/// Initialize a system. The callback returns zero to stop integration.
/// It receives time, mutable equation values, and a slice to fill with derivatives.
/// It may modify the values in place, including intermediate stage values.
///
/// C prototype: `rke_variables rke_init(int number, int (*routine)());`
pub fn rke_init<F>(number: usize, routine: F) -> RkeVariables<F>
where
    F: FnMut(f64, &mut [f64], &mut [f64]) -> i32,
{
    RkeVariables {
        n_equations: number,
        eval_routine: routine,
        minimum_step: 0.000001,
        maximum_step: 1000000.0,
        current_step: 1.0,
        error_slope: 0.0000001,
        error_bias: 0.00000001,
        accepted_steps: 0,
        rejected_steps: 0,
    }
}

/// Release a system's state.
///
/// C prototype: `void rke_term(rke_variables var);`
pub fn rke_term<F>(var: RkeVariables<F>) {
    drop(var);
}

/// Integrate towards `aimed_time`; returns zero on callback or convergence failure.
///
/// Controls and step counters persist between calls. The slice must contain at
/// least `var.n_equations` values; any remaining elements are untouched. On
/// failure, time retains the last completed step. Callbacks receive the actual
/// state or intermediate stage slice, so their changes persist in that array
/// until the solver overwrites it. Changes to the caller's values persist even
/// if the callback fails; midpoint changes feed the remaining stages and result.
/// Like C, the accepted counter increments before the final callback, even if
/// that callback fails.
/// Retrying a rejected trial halves `current_step`, which can exceed the
/// remaining distance to the target; the retry is not clipped to that distance.
/// Integration stops within half of `minimum_step` of the target.
///
/// C prototype: `int rke_solve(rke_variables var, double *time, double *variables, double aimed_time);`
pub fn rke_solve<F>(
    var: &mut RkeVariables<F>,
    time: &mut f64,
    variables: &mut [f64],
    aimed_time: f64,
) -> i32
where
    F: FnMut(f64, &mut [f64], &mut [f64]) -> i32,
{
    let n = var.n_equations;
    assert!(variables.len() >= n);
    let variables = &mut variables[..n];
    // C passes a byte count to ealloc1double, allocating eight times the
    // required elements, and leaks the arrays. Only the first n are accessed.
    let mut dp = vec![0.0; n];
    let mut vt = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut d = vec![0.0; n];
    let mut a1 = vec![0.0; n];
    let mut a2 = vec![0.0; n];
    let mut a3 = vec![0.0; n];
    let mut a4 = vec![0.0; n];
    let mut a5 = vec![0.0; n];
    let mut a6 = vec![0.0; n];
    let mut a7 = vec![0.0; n];
    while 2.0 * (aimed_time - *time).abs() > var.minimum_step {
        let mut whole_step = aimed_time - *time;
        if whole_step > 0.0 {
            if whole_step > var.current_step {
                whole_step = var.current_step;
            }
        } else if whole_step < -var.current_step {
            whole_step = -var.current_step;
        }
        if (var.eval_routine)(*time, variables, &mut dp) == 0 {
            return 0;
        }
        loop {
            let quarter_step = 0.25 * whole_step;
            let half_step = quarter_step + quarter_step;
            let three_quarter_step = half_step + quarter_step;
            for k in 0..n {
                a1[k] = half_step * dp[k];
                v[k] = variables[k] + 0.5 * a1[k];
            }
            if (var.eval_routine)(*time + quarter_step, &mut v, &mut d) == 0 {
                return 0;
            }
            for k in 0..n {
                a2[k] = half_step * d[k];
                v[k] = variables[k] + 0.25 * (a1[k] + a2[k]);
            }
            if (var.eval_routine)(*time + quarter_step, &mut v, &mut d) == 0 {
                return 0;
            }
            for k in 0..n {
                a3[k] = half_step * d[k];
                v[k] = variables[k] + (-a2[k] + a3[k] + a3[k]);
            }
            if (var.eval_routine)(*time + half_step, &mut v, &mut d) == 0 {
                return 0;
            }
            for k in 0..n {
                a4[k] = half_step * d[k];
                vt[k] = variables[k] + (a1[k] + 4.0 * a3[k] + a4[k]) / 6.0;
            }
            if (var.eval_routine)(*time + half_step, &mut vt, &mut d) == 0 {
                return 0;
            }
            for k in 0..n {
                a5[k] = half_step * d[k];
                v[k] = vt[k] + 0.5 * a5[k];
            }
            if (var.eval_routine)(*time + three_quarter_step, &mut v, &mut d) == 0 {
                return 0;
            }
            for k in 0..n {
                a6[k] = half_step * d[k];
                v[k] = vt[k] + 0.25 * (a5[k] + a6[k]);
            }
            if (var.eval_routine)(*time + three_quarter_step, &mut v, &mut d) == 0 {
                return 0;
            }
            for k in 0..n {
                a7[k] = half_step * d[k];
                v[k] = variables[k]
                    + (-a1[k] - 96.0 * a2[k] + 92.0 * a3[k] - 121.0 * a4[k]
                        + 144.0 * a5[k]
                        + 6.0 * a6[k]
                        - 12.0 * a7[k])
                        / 6.0;
            }
            if (var.eval_routine)(*time + whole_step, &mut v, &mut d) == 0 {
                return 0;
            }
            let mut within_tolerance = true;
            let mut all_errors_small = true;
            for k in 0..n {
                let estimated_error = ((-a1[k] + 4.0 * a3[k] + 17.0 * a4[k] - 23.0 * a5[k]
                    + 4.0 * a7[k]
                    - half_step * d[k])
                    / 90.0)
                    .abs();
                let allowable_error =
                    whole_step.abs() * (var.error_slope * vt[k].abs() + var.error_bias);
                if estimated_error > allowable_error {
                    within_tolerance = false;
                    break;
                } else if estimated_error > 0.02 * allowable_error {
                    all_errors_small = false;
                }
            }
            if within_tolerance {
                var.accepted_steps = var.accepted_steps.wrapping_add(1);
                for k in 0..n {
                    v[k] = vt[k] + (-a6[k] + a7[k] + a7[k]);
                }
                if (var.eval_routine)(*time + whole_step, &mut v, &mut d) == 0 {
                    return 0;
                }
                *time += whole_step;
                for k in 0..n {
                    variables[k] = vt[k] + (a5[k] + 4.0 * a7[k] + half_step * d[k]) / 6.0;
                }
                if all_errors_small && whole_step.abs() == var.current_step {
                    if 2.0 * var.current_step > var.maximum_step {
                        var.current_step = var.maximum_step;
                    } else {
                        var.current_step *= 2.0;
                    }
                }
                break;
            } else {
                var.rejected_steps = var.rejected_steps.wrapping_add(1);
                if whole_step.abs() > var.minimum_step {
                    if var.current_step < 2.0 * var.minimum_step {
                        var.current_step = var.minimum_step;
                    } else {
                        var.current_step *= 0.5;
                    }
                    if aimed_time > *time {
                        whole_step = var.current_step;
                    } else {
                        whole_step = -var.current_step;
                    }
                } else {
                    return 0;
                }
            }
        }
    }
    1
}
