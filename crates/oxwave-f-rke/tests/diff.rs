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

//! Differential tests for `par/lib/rke.c`.
//! Credits: Copyright 1988: François Pinard <pinard@iro.umontreal.ca>, 1988.
//! Made available for SU by the author's permission, 1 Sept 2006.
//! Modified: John Stockwell, Center for Wave Phenomena (CWP), Colorado School
//! of Mines (CSM), August 2004, for integration into the CWP/SU Seismic Unix
//! package. Some stylistic changes were made, and some SU routines were
//! incorporated into the modified code.
use oxwave_cref::XorShift;
use oxwave_f_rke::{
    rke_init as rust_init, rke_solve as rust_solve, rke_term as rust_term, RkeVariables,
};
use std::cell::Cell;

#[repr(C)]
struct CState {
    n: i32,
    callback: unsafe extern "C" fn(f64, *mut f64, *mut f64) -> i32,
    min: f64,
    max: f64,
    step: f64,
    slope: f64,
    bias: f64,
    accepted: i32,
    rejected: i32,
}
extern "C" {
    // rke_variables rke_init(int number, int (*routine)());
    fn rke_init(n: i32, f: unsafe extern "C" fn(f64, *mut f64, *mut f64) -> i32) -> *mut CState;
    // int rke_solve(rke_variables var, double *time, double *variables, double aimed_time);
    fn rke_solve(s: *mut CState, t: *mut f64, y: *mut f64, target: f64) -> i32;
    // void rke_term(rke_variables var);
    fn rke_term(s: *mut CState);
}
// Callbacks run synchronously on the calling thread. Keep each test's C
// callback parameters separate even when the test harness runs in parallel.
thread_local! {
    static N: Cell<usize> = const { Cell::new(0) };
    static CALL: Cell<usize> = const { Cell::new(0) };
    static FAIL: Cell<usize> = const { Cell::new(0) };
    static MODE: Cell<usize> = const { Cell::new(0) };
    static RATE: Cell<f64> = const { Cell::new(0.4) };
    static TRACE: std::cell::RefCell<Vec<u64>> = const { std::cell::RefCell::new(Vec::new()) };
    static STATE_TRACE: std::cell::RefCell<Vec<Vec<u64>>> = const { std::cell::RefCell::new(Vec::new()) };
    static MUTATE_ON: Cell<usize> = const { Cell::new(0) };
    static MUTATION_SCALE: Cell<f64> = const { Cell::new(1.0) };
    // a, b, c, p, K, and thresh from voltlotka's callback.
    static POPULATION_PARAMS: Cell<[f64; 6]> = const { Cell::new([0.5, 0.02, 0.75, 0.005, 1e6, 1.0]) };
    static CLAMPED: Cell<usize> = const { Cell::new(0) };
}
fn eval(t: f64, y: &mut [f64], d: &mut [f64]) -> i32 {
    TRACE.with(|trace| trace.borrow_mut().push(t.to_bits()));
    let call = CALL.with(|c| {
        let n = c.get() + 1;
        c.set(n);
        n
    });
    if FAIL.with(|f| f.get()) == call {
        return 0;
    }
    let rate = RATE.with(Cell::get);
    let mode = MODE.with(Cell::get);
    for k in 0..y.len() {
        d[k] = match mode {
            0 => -rate * y[k] + t * 0.1 + (k as f64) * 0.001,
            1 => 0.0,
            2 => rate * (k as f64 + 1.0),
            3 => -rate * y[k] / (1.0 + y[k].abs()) + 0.1 * t * t,
            // Coupled pairs of oscillator equations; odd lengths have one
            // extra decaying component.
            _ if k % 2 == 0 && k + 1 < y.len() => rate * y[k + 1],
            _ if k % 2 == 1 => -rate * y[k - 1],
            _ => -rate * y[k],
        };
    }
    1
}
unsafe extern "C" fn c_eval(t: f64, y: *mut f64, d: *mut f64) -> i32 {
    let n = N.with(Cell::get);
    eval(
        t,
        std::slice::from_raw_parts_mut(y, n),
        std::slice::from_raw_parts_mut(d, n),
    )
}
fn reset(n: usize, fail: usize) {
    N.with(|x| x.set(n));
    CALL.with(|x| x.set(0));
    FAIL.with(|x| x.set(fail));
    TRACE.with(|x| x.borrow_mut().clear());
    STATE_TRACE.with(|x| x.borrow_mut().clear());
    CLAMPED.with(|x| x.set(0));
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1e-5 * a.abs().max(b.abs()).max(1e-300),
        "{a:?} != {b:?}"
    );
    // The original libraries disable contraction. All finite cases here
    // additionally match bit for bit, including signed zero.
    assert_eq!(a.to_bits(), b.to_bits(), "{a:?} != {b:?}");
}

#[test]
fn rke_init_matches_c() {
    let mut rng = XorShift::new(1);
    for case in 0..3000 {
        let sizes = [0, 1, 2, 3, 7, 16, 127, 256, 4096];
        let n = if case % 2 == 0 {
            sizes[case % sizes.len()]
        } else {
            rng.range_i32(0, 10000) as usize
        };
        let r = rust_init(n, eval);
        unsafe {
            let c = rke_init(n as i32, c_eval);
            assert!(!c.is_null());
            let s = &*c;
            assert_eq!(s.callback as usize, c_eval as *const () as usize);
            assert_eq!(r.eval_routine as *const (), eval as *const ());
            assert_eq!(r.n_equations, s.n as usize);
            assert_eq!(
                (
                    r.minimum_step,
                    r.maximum_step,
                    r.current_step,
                    r.error_slope,
                    r.error_bias,
                    r.accepted_steps,
                    r.rejected_steps
                ),
                (s.min, s.max, s.step, s.slope, s.bias, s.accepted, s.rejected)
            );
            rke_term(c);
        }
        rust_term(r);
    }
}

#[test]
fn rke_solve_matches_c() {
    let mut rng = XorShift::new(2);
    let sizes = [0, 1, 2, 3, 7, 16, 31, 64, 127, 256, 257];
    let mut successes = 0;
    let mut failures = 0;
    let mut rejections = 0;
    for case in 0..3000 {
        let n = if case % 251 == 0 {
            4096
        } else {
            sizes[case % sizes.len()]
        };
        MODE.with(|m| m.set(case % 5));
        RATE.with(|r| r.set(rng.uniform(-2.0, 2.0)));
        let scale = [1e-9, 1.0, 1e9][case % 3];
        let mut y: Vec<f64> = (0..n).map(|_| scale * rng.uniform(-10.0, 10.0)).collect();
        if case % 11 == 0 {
            y.fill(0.0);
        }
        // Trailing elements are outside n_equations and must remain untouched.
        y.extend_from_slice(&[-0.0, 123.0]);
        let mut cy = y.clone();
        let mut t = rng.uniform(-2.0, 2.0);
        let mut ct = t;
        let target = if case % 13 == 0 {
            t
        } else {
            rng.uniform(-2.0, 2.0)
        };
        let fail = if case % 7 == 0 { 1 + case % 10 } else { 0 };
        let mut r = rust_init(n, eval);
        unsafe {
            let c = rke_init(n as i32, c_eval);
            if case % 5 == 0 {
                r.minimum_step = 0.05;
                (*c).min = r.minimum_step;
                r.error_bias = 1e-15;
                (*c).bias = r.error_bias;
                r.error_slope = 0.0;
                (*c).slope = 0.0;
            }
            if case % 3 == 0 {
                r.current_step = 0.125;
                (*c).step = r.current_step;
                r.maximum_step = 0.5;
                (*c).max = r.maximum_step;
            }
            reset(n, fail);
            let cr = rke_solve(c, &mut ct, cy.as_mut_ptr(), target);
            let calls = CALL.with(Cell::get);
            let trace = TRACE.with(|x| x.borrow().clone());
            reset(n, fail);
            let rr = rust_solve(&mut r, &mut t, &mut y, target);
            assert_eq!(rr, cr, "case {case}");
            assert_eq!(CALL.with(Cell::get), calls, "case {case}");
            TRACE.with(|x| assert_eq!(*x.borrow(), trace, "case {case}"));
            close(t, ct);
            for (a, b) in y.iter().zip(&cy) {
                close(*a, *b);
            }
            close(r.current_step, (*c).step);
            assert_eq!(
                (r.accepted_steps, r.rejected_steps),
                ((*c).accepted, (*c).rejected)
            );
            successes += usize::from(cr != 0);
            failures += usize::from(cr == 0);
            rejections += (*c).rejected as usize;

            // Continue using the same state and accumulated counters, changing
            // direction and controls between solves as permitted by rke.h.
            let target = t + rng.uniform(-0.5, 0.5);
            r.error_slope = 1e-7;
            (*c).slope = r.error_slope;
            r.error_bias = 1e-8;
            (*c).bias = r.error_bias;
            reset(n, 0);
            let cr = rke_solve(c, &mut ct, cy.as_mut_ptr(), target);
            let calls = CALL.with(Cell::get);
            let trace = TRACE.with(|x| x.borrow().clone());
            reset(n, 0);
            assert_eq!(rust_solve(&mut r, &mut t, &mut y, target), cr);
            assert_eq!(CALL.with(Cell::get), calls);
            TRACE.with(|x| assert_eq!(*x.borrow(), trace));
            close(t, ct);
            for (a, b) in y.iter().zip(&cy) {
                close(*a, *b);
            }
            close(r.current_step, (*c).step);
            assert_eq!(
                (r.accepted_steps, r.rejected_steps),
                ((*c).accepted, (*c).rejected)
            );
            rke_term(c);
        }
        rust_term(r);
    }
    assert!(successes > 1000 && failures > 100 && rejections > 100);
}

#[test]
fn rke_solve_step_boundaries_match_c() {
    let mut rng = XorShift::new(4);
    MODE.with(|m| m.set(1)); // Zero derivatives make exact step decisions visible.
    for case in 0..2000 {
        let n = 1 + case % 7;
        let mut r = rust_init(n, eval);
        // Binary fractions make the half-minimum stopping boundary exact.
        r.minimum_step = 2.0_f64.powf(rng.range_i32(-20, -4) as f64);
        r.current_step = 2.0 * r.minimum_step;
        r.maximum_step = 3.0 * r.minimum_step;
        let sign = if case % 2 == 0 { 1.0 } else { -1.0 };
        let target = sign * [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 8.0][case % 7] * r.minimum_step;
        let mut y = vec![rng.uniform(-100.0, 100.0); n];
        let mut cy = y.clone();
        let (mut t, mut ct) = (0.0, 0.0);
        unsafe {
            let c = rke_init(n as i32, c_eval);
            (*c).min = r.minimum_step;
            (*c).step = r.current_step;
            (*c).max = r.maximum_step;
            reset(n, 0);
            assert_eq!(rke_solve(c, &mut ct, cy.as_mut_ptr(), target), 1);
            let calls = CALL.with(Cell::get);
            reset(n, 0);
            assert_eq!(rust_solve(&mut r, &mut t, &mut y, target), 1);
            assert_eq!(CALL.with(Cell::get), calls);
            if target.abs() <= 0.5 * r.minimum_step {
                assert_eq!(calls, 0);
            }
            close(t, ct);
            close(r.current_step, (*c).step);
            assert_eq!(y, cy);
            assert_eq!(
                (r.accepted_steps, r.rejected_steps),
                ((*c).accepted, (*c).rejected)
            );
            rke_term(c);
        }
        rust_term(r);
    }
}

#[test]
fn rke_solve_callback_failure_at_every_stage_matches_c() {
    MODE.with(|m| m.set(2));
    RATE.with(|r| r.set(-0.5));
    // Nine calls per accepted step, including the call after accepted_steps
    // increments. Also fail during the second step after time has advanced.
    for fail in 1..=18 {
        let mut r = rust_init(2, eval);
        r.current_step = 0.125;
        r.maximum_step = 0.125;
        let (mut t, mut ct) = (0.0, 0.0);
        let mut y = [1.0, -1.0];
        let mut cy = y;
        unsafe {
            let c = rke_init(2, c_eval);
            (*c).step = r.current_step;
            (*c).max = r.maximum_step;
            reset(2, fail);
            assert_eq!(rke_solve(c, &mut ct, cy.as_mut_ptr(), 1.0), 0);
            assert_eq!(CALL.with(Cell::get), fail);
            let trace = TRACE.with(|x| x.borrow().clone());
            reset(2, fail);
            assert_eq!(rust_solve(&mut r, &mut t, &mut y, 1.0), 0);
            assert_eq!(CALL.with(Cell::get), fail);
            TRACE.with(|x| assert_eq!(*x.borrow(), trace));
            close(t, ct);
            for k in 0..2 {
                close(y[k], cy[k]);
            }
            assert_eq!(
                (r.accepted_steps, r.rejected_steps),
                ((*c).accepted, (*c).rejected)
            );
            rke_term(c);
        }
        rust_term(r);
    }
}

#[test]
fn rke_solve_rejected_trials_match_c() {
    let mut rng = XorShift::new(5);
    MODE.with(|m| m.set(1)); // Zero derivatives give an exactly zero error estimate.
    for case in 0..2048 {
        let sizes = [1, 2, 3, 7, 16, 31];
        let n = sizes[case % sizes.len()];
        let mut r = rust_init(n, eval);
        r.minimum_step = 2.0_f64.powf(rng.range_i32(-20, -4) as f64);
        let exponent = rng.range_i32(5, 10);
        r.current_step = r.minimum_step * 2.0_f64.powf(exponent as f64);
        // C does not validate tolerances. A negative bias forces every trial
        // to reject, without relying on roundoff or a particular ODE stiffness.
        r.error_slope = 0.0;
        r.error_bias = -rng.uniform(1e-12, 1e-6);
        r.accepted_steps = rng.range_i32(0, 1000);
        r.rejected_steps = rng.range_i32(0, 1000);
        let initial_accepted = r.accepted_steps;
        let initial_rejected = r.rejected_steps;
        let sign = if case % 2 == 0 { 1.0 } else { -1.0 };
        let mut t = rng.range_i32(-1024, 1024) as f64 * r.minimum_step;
        let mut ct = t;
        let initial_time = t;
        let mut distance = 1.5 * r.minimum_step;
        let fail = if case % 3 == 2 {
            // Seven callbacks per rejected trial after the initial derivative.
            // Exercise every callback on each of the first two retries.
            9 + (case / 3) % 14
        } else {
            0
        };
        if case % 3 == 1 {
            // Exercise current_step < 2 * minimum_step: clamp directly to
            // the minimum instead of halving below it.
            r.current_step = 1.5 * r.minimum_step;
            distance = 1.25 * r.minimum_step;
        }
        // In the other cases, rejection expands the initially clipped trial
        // past the target because C halves current_step, not whole_step.
        let target = t + sign * distance;
        let scale = [1e-12, 1.0, 1e12][(case / 3) % 3];
        let mut y: Vec<f64> = (0..n).map(|_| scale * rng.uniform(-10.0, 10.0)).collect();
        if case % 11 == 0 {
            y.fill(0.0);
        }
        y.extend_from_slice(&[-0.0, -123.0]);
        let initial_y = y.clone();
        let mut cy = y.clone();
        unsafe {
            let c = rke_init(n as i32, c_eval);
            (*c).min = r.minimum_step;
            (*c).step = r.current_step;
            (*c).slope = r.error_slope;
            (*c).bias = r.error_bias;
            (*c).accepted = r.accepted_steps;
            (*c).rejected = r.rejected_steps;
            reset(n, fail);
            assert_eq!(rke_solve(c, &mut ct, cy.as_mut_ptr(), target), 0);
            let calls = CALL.with(Cell::get);
            let trace = TRACE.with(|x| x.borrow().clone());
            reset(n, fail);
            assert_eq!(rust_solve(&mut r, &mut t, &mut y, target), 0);
            assert_eq!(CALL.with(Cell::get), calls, "case {case}");
            TRACE.with(|x| assert_eq!(*x.borrow(), trace, "case {case}"));
            close(t, ct);
            close(t, initial_time);
            close(r.current_step, (*c).step);
            for k in 0..y.len() {
                close(y[k], cy[k]);
                close(y[k], initial_y[k]);
            }
            assert_eq!(
                (r.accepted_steps, r.rejected_steps),
                ((*c).accepted, (*c).rejected)
            );
            assert_eq!(r.accepted_steps, initial_accepted);
            let rejected = r.rejected_steps - initial_rejected;
            if fail != 0 {
                assert_eq!(calls, fail);
                assert_eq!(rejected as usize, (fail - 2) / 7);
            } else {
                close(r.current_step, r.minimum_step);
                assert_eq!(rejected, if case % 3 == 1 { 2 } else { exponent + 1 });
                assert_eq!(calls, 1 + 7 * rejected as usize);
            }
            rke_term(c);
        }
        rust_term(r);
    }
}

// Record the actual inputs at every callback so that copying a stage before
// calling user code cannot hide lost mutations behind matching derivatives.
fn record_state(t: f64, y: &[f64]) -> usize {
    TRACE.with(|x| x.borrow_mut().push(t.to_bits()));
    STATE_TRACE.with(|x| x.borrow_mut().push(y.iter().map(|v| v.to_bits()).collect()));
    CALL.with(|x| {
        let call = x.get() + 1;
        x.set(call);
        call
    })
}

fn mutating_eval(t: f64, y: &mut [f64], d: &mut [f64]) -> i32 {
    let call = record_state(t, y);
    let mutate_on = MUTATE_ON.with(Cell::get);
    let scale = MUTATION_SCALE.with(Cell::get);
    for k in 0..y.len() {
        if mutate_on == 0 || mutate_on == call {
            y[k] = 0.5 * y[k] + scale * (k as f64 + 1.0);
        }
        d[k] = 0.0;
    }
    // Mutate before reporting failure, as a C callback is allowed to do.
    i32::from(FAIL.with(Cell::get) != call)
}

unsafe extern "C" fn c_mutating_eval(t: f64, y: *mut f64, d: *mut f64) -> i32 {
    let n = N.with(Cell::get);
    mutating_eval(
        t,
        std::slice::from_raw_parts_mut(y, n),
        std::slice::from_raw_parts_mut(d, n),
    )
}

fn population_eval(t: f64, y: &mut [f64], d: &mut [f64]) -> i32 {
    record_state(t, y);
    let [a, b, c, p, capacity, threshold] = POPULATION_PARAMS.with(Cell::get);
    // Same state changes and equation order as the callback in
    // par/main/ordinary_differential_equations/voltlotka.c.
    for k in 0..2 {
        if y[k] <= threshold {
            if y[k] != 0.0 {
                CLAMPED.with(|x| x.set(x.get() + 1));
            }
            y[k] = 0.0;
        }
    }
    d[0] = a * y[0] * (1.0 - y[0] / capacity) - b * y[0] * y[1];
    d[1] = -c * y[1] + p * y[0] * y[1];
    1
}

unsafe extern "C" fn c_population_eval(t: f64, y: *mut f64, d: *mut f64) -> i32 {
    population_eval(
        t,
        std::slice::from_raw_parts_mut(y, 2),
        std::slice::from_raw_parts_mut(d, 2),
    )
}

unsafe fn compare_mutating_solve<F>(
    r: &mut RkeVariables<F>,
    c: *mut CState,
    t: &mut f64,
    y: &mut [f64],
    target: f64,
    fail: usize,
) -> i32
where
    F: FnMut(f64, &mut [f64], &mut [f64]) -> i32,
{
    let mut ct = *t;
    let mut cy = y.to_vec();
    reset(r.n_equations, fail);
    let cr = rke_solve(c, &mut ct, cy.as_mut_ptr(), target);
    let calls = CALL.with(Cell::get);
    let trace = TRACE.with(|x| x.borrow().clone());
    let states = STATE_TRACE.with(|x| x.borrow().clone());
    let clamped = CLAMPED.with(Cell::get);
    reset(r.n_equations, fail);
    let rr = rust_solve(r, t, y, target);
    assert_eq!(rr, cr);
    assert_eq!(CALL.with(Cell::get), calls);
    assert_eq!(CLAMPED.with(Cell::get), clamped);
    TRACE.with(|x| assert_eq!(*x.borrow(), trace));
    STATE_TRACE.with(|x| assert_eq!(*x.borrow(), states));
    close(*t, ct);
    for k in 0..y.len() {
        close(y[k], cy[k]);
    }
    close(r.current_step, (*c).step);
    assert_eq!(
        (r.accepted_steps, r.rejected_steps),
        ((*c).accepted, (*c).rejected)
    );
    rr
}

#[test]
fn rke_solve_state_mutations_at_every_stage_match_c() {
    let mut rng = XorShift::new(6);
    let sizes = [1, 2, 3, 7, 16, 31, 64, 127, 256, 257];
    for case in 0..3072 {
        let n = if case % 251 == 0 {
            4096
        } else {
            sizes[case % sizes.len()]
        };
        let scale = [1e-12, 1.0, 1e12][(case / 4) % 3];
        MUTATION_SCALE.with(|x| x.set(scale * rng.uniform(-2.0, 2.0)));
        // Isolate each of the nine evaluation sites over two accepted steps,
        // then mutate every site on failures and rejected trials.
        MUTATE_ON.with(|x| {
            x.set(if case % 4 == 0 {
                1 + (case / 4) % 18
            } else {
                0
            })
        });
        let fail = match case % 4 {
            1 => 1 + (case / 4) % 18,
            3 => 9 + (case / 4) % 14,
            _ => 0,
        };
        let mut r = rust_init(n, mutating_eval);
        r.minimum_step = 0.0078125;
        r.current_step = 0.125;
        r.maximum_step = 0.125;
        r.error_slope = 0.0;
        // Zero derivatives give zero estimated error. A negative tolerance
        // forces retries and convergence failure independently of mutations.
        r.error_bias = if case % 4 >= 2 { -1.0 } else { 1.0 };
        let mut y: Vec<f64> = (0..n).map(|_| scale * rng.uniform(-10.0, 10.0)).collect();
        if case % 11 == 0 {
            y.fill(0.0);
        }
        y.extend_from_slice(&[-0.0, -123.0]);
        let mut t = rng.range_i32(-16, 16) as f64 * 0.125;
        let sign = if (case / 4) % 2 == 0 { 1.0 } else { -1.0 };
        let target = t + sign * 0.25;
        unsafe {
            let c = rke_init(n as i32, c_mutating_eval);
            (*c).min = r.minimum_step;
            (*c).max = r.maximum_step;
            (*c).step = r.current_step;
            (*c).slope = r.error_slope;
            (*c).bias = r.error_bias;
            let result = compare_mutating_solve(&mut r, c, &mut t, &mut y, target, fail);
            assert_eq!(result, i32::from(case % 4 == 0));
            if fail != 0 {
                assert_eq!(CALL.with(Cell::get), fail);
            } else if case % 4 == 0 {
                assert_eq!(CALL.with(Cell::get), 18);
                assert_eq!(r.accepted_steps, 2);
            } else {
                assert_eq!(r.rejected_steps, 5);
                close(r.current_step, r.minimum_step);
            }
            close(y[n], -0.0);
            close(y[n + 1], -123.0);
            rke_term(c);
        }
        rust_term(r);
    }
}

#[test]
fn rke_solve_midpoint_clamping_matches_c() {
    for sign in [-1.0, 1.0] {
        POPULATION_PARAMS.with(|x| x.set([0.0, 0.0, sign, 0.0, 1e6, 1.0]));
        let mut r = rust_init(2, population_eval);
        r.current_step = 0.125;
        r.error_slope = 0.0;
        r.error_bias = 1.0;
        let mut t = 0.0;
        // The predator starts above the threshold, then crosses it during
        // the first half-step. Its midpoint must be changed in place.
        let mut y = [0.0, 1.04];
        unsafe {
            let c = rke_init(2, c_population_eval);
            (*c).step = r.current_step;
            (*c).slope = r.error_slope;
            (*c).bias = r.error_bias;
            assert_eq!(
                compare_mutating_solve(&mut r, c, &mut t, &mut y, sign * 0.125, 0),
                1
            );
            assert_eq!((r.accepted_steps, r.rejected_steps), (1, 0));
            STATE_TRACE.with(|x| {
                let states = x.borrow();
                assert_eq!(states.len(), 9);
                let midpoint = f64::from_bits(states[4][1]);
                assert!(midpoint > 0.0 && midpoint <= 1.0);
                assert_eq!(states[5][1], 0.0_f64.to_bits());
            });
            close(y[1], 0.0);
            rke_term(c);
        }
        rust_term(r);
    }
}

#[test]
fn rke_solve_population_clamping_matches_c() {
    let mut rng = XorShift::new(7);
    let mut clamped_cases = 0;
    let mut rejected_cases = 0;
    for case in 0..3000 {
        let threshold = [0.0, 1e-9, 1.0, 10.0][case % 4];
        POPULATION_PARAMS.with(|x| {
            x.set([
                rng.uniform(-0.5, 1.0),
                rng.uniform(0.0, 0.05),
                rng.uniform(0.0, 1.5),
                rng.uniform(0.0, 0.01),
                rng.uniform(100.0, 1e6),
                threshold,
            ])
        });
        let mut y = [rng.uniform(-5.0, 100.0), rng.uniform(-5.0, 50.0), -0.0];
        // Exact threshold, immediately above it, negative values and zeros.
        if case % 3 == 0 {
            y[0] = threshold;
        } else if case % 3 == 1 {
            y[1] = threshold * (1.0 + rng.uniform(1e-10, 0.05));
        }
        let mut t = rng.uniform(-2.0, 2.0);
        let sign = if case % 2 == 0 { 1.0 } else { -1.0 };
        let target = t + sign * rng.uniform(0.125, 0.5);
        let mut r = rust_init(2, population_eval);
        r.minimum_step = 1e-4;
        r.current_step = 0.125;
        r.maximum_step = 0.25;
        unsafe {
            let c = rke_init(2, c_population_eval);
            (*c).min = r.minimum_step;
            (*c).step = r.current_step;
            (*c).max = r.maximum_step;
            compare_mutating_solve(&mut r, c, &mut t, &mut y, target, 0);
            clamped_cases += usize::from(CLAMPED.with(Cell::get) > 0);
            rejected_cases += usize::from(r.rejected_steps > 0);
            // Reuse the mutated state and adaptive controls, reversing time.
            let target = t - sign * 0.125;
            compare_mutating_solve(&mut r, c, &mut t, &mut y, target, 0);
            close(y[2], -0.0);
            rke_term(c);
        }
        rust_term(r);
    }
    assert!(clamped_cases > 1000);
    assert!(rejected_cases > 100);
}

#[test]
fn rke_term_matches_c() {
    let mut rng = XorShift::new(3);
    for _ in 0..3000 {
        let n = rng.range_i32(0, 1000) as usize;
        unsafe {
            rke_term(rke_init(n as i32, c_eval));
        }
        let marker = std::rc::Rc::new(());
        let captured = marker.clone();
        let r = rust_init(n, move |_, _: &mut [f64], _: &mut [f64]| {
            let _ = &captured;
            1
        });
        rust_term(r);
        assert_eq!(std::rc::Rc::strong_count(&marker), 1);
    }
}
