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

//! Differential tests for Refl/lib/auxil.c, compute.c, convolve.c,
//! psvreflect.c, shreflect.c and synthetic.c.
//! Credits: Copyright (c) Colorado School of Mines, 2011; the original
//! reflectivity sources' Credits sections name no individual authors.
use oxwave_cref::XorShift;
use oxwave_cwp::complex::Complex;
use oxwave_lib_refl as refl;

#[link(name = "reflect", kind = "static")]
extern "C" {
    fn compute_block_size(np: i32, b: *mut i32, n: *mut i32, l: *mut i32);
    fn source_receiver_type(
        n: i32,
        nl: i32,
        ls: i32,
        obs: *const i32,
        cl: *const f32,
        ct: *const f32,
        a: *mut i32,
        f: *mut i32,
    );
    fn apply_earth_flattening(
        n: i32,
        z: f32,
        cl: *mut f32,
        ct: *mut f32,
        rho: *mut f32,
        t: *mut f32,
    );
    fn compute_Hanning_window(i: i32, w: i32, f1: i32, f2: i32, nf: i32, win: *mut Complex);
    fn compute_prs(
        n: i32,
        l: *const i32,
        al: *const Complex,
        at: *const Complex,
        rho: *const f32,
        prs: *mut Complex,
    );
    fn compute_alpha_betha(
        ls: i32,
        n: i32,
        ats: Complex,
        rho: f32,
        gl: *const Complex,
        gt: *const Complex,
        a: *mut Complex,
        b: *mut Complex,
    );
    fn red_vel_factor(x: f32, vel: f32, lag: f32, wp: Complex, ws: Complex, out: *mut Complex);
    fn compute_pwin(
        wtype: i32,
        n: i32,
        bp: f32,
        dp: f32,
        w: f32,
        decay: f32,
        pw1: f32,
        pw2: f32,
        pw3: f32,
        pw4: f32,
        pwin: *mut Complex,
        p: *mut Complex,
        pp: *mut Complex,
    );
    fn compute_gl_gt_gam(
        wtype: i32,
        nl: i32,
        n: i32,
        al: *const Complex,
        at: *const Complex,
        rho: *const f32,
        pp: *const Complex,
        gl: *mut Complex,
        gt: *mut Complex,
        gam: *mut Complex,
    );
    fn compute_moment_tensor(
        wtype: i32,
        phi: f32,
        lambda: f32,
        delta: f32,
        phis: f32,
        m0: f32,
        m1: *mut f32,
        m2: *mut f32,
        m3: *mut f32,
    );

    fn parameter_interpolation(
        n: i32,
        layers: *const i32,
        ns: *const i32,
        th: *const f32,
        cl: *mut f32,
        ql: *mut f32,
        ct: *mut f32,
        qt: *mut f32,
        rho: *mut f32,
        t: *mut f32,
    );

    fn compute_slowness(
        wtype: i32,
        lsource: i32,
        w: f32,
        p2w: f32,
        np: *mut i32,
        fp: *mut f32,
        fs: f32,
        e4: f32,
        dk: f32,
        decay: f32,
        xmax: f32,
        dp: *mut f32,
        pw1: *mut f32,
        pw2: *mut f32,
        pw3: *mut f32,
        pw4: *mut f32,
        cdp: *mut Complex,
        cl: *const f32,
        ct: *const f32,
        t: *const f32,
    );

}

fn close(a: f32, b: f32) {
    assert!(
        a == b
            || (a.is_nan() && b.is_nan())
            || (a - b).abs() <= 1e-5 * a.abs().max(b.abs()).max(1e-30),
        "{a} != {b}"
    );
}

// Include singleton, even, prime, power-of-two and large allocations without
// discarding any generated input. Other iterations vary sizes independently.
fn sample_size(rng: &mut XorShift, i: usize, maximum: u32) -> usize {
    if i % 8 == 0 {
        [1, 2, 3, 7, 16, 31, 64, 257, 1024][(i / 8) % 9]
    } else {
        (rng.next_u32() % maximum + 1) as usize
    }
}
#[test]
fn block_size() {
    let mut rng = XorShift::new(9);
    for i in 0..10000 {
        let np = if i < 5 {
            [0, 1, 100, 101, -1][i]
        } else {
            (rng.next_u32() % 100000) as i32
        };
        let (mut b, mut n, mut l) = (0, 0, 0);
        let (mut cb, mut cn, mut cl) = (0, 0, 0);
        refl::compute_block_size(np, &mut b, &mut n, &mut l);
        unsafe { compute_block_size(np, &mut cb, &mut cn, &mut cl) };
        assert_eq!((b, n, l), (cb, cn, cl));
    }
}
#[test]
fn receiver_type() {
    let mut rng = XorShift::new(19);
    for case in 0..3000 {
        let n = sample_size(&mut rng, case, 30);
        let nr = sample_size(&mut rng, case + 8, 30);
        let ls = (rng.next_u32() % n as u32 + 1) as i32;
        let cl: Vec<_> = (0..n).map(|_| (rng.next_u32() % 4) as f32).collect();
        let ct: Vec<_> = (0..n).map(|_| (rng.next_u32() % 4) as f32).collect();
        let obs: Vec<_> = (0..nr)
            .map(|_| (rng.next_u32() % n as u32 + 1) as i32)
            .collect();
        let (mut a, mut f, mut ca, mut cf) = (vec![0; nr], vec![0; nr], vec![0; nr], vec![0; nr]);
        refl::source_receiver_type(ls, &obs, &cl, &ct, &mut a, &mut f);
        unsafe {
            source_receiver_type(
                nr as i32,
                n as i32,
                ls,
                obs.as_ptr(),
                cl.as_ptr(),
                ct.as_ptr(),
                ca.as_mut_ptr(),
                cf.as_mut_ptr(),
            )
        };
        assert_eq!((a, f), (ca, cf));
    }
}
#[test]
fn earth_flattening() {
    let mut rng = XorShift::new(22);
    for case in 0..3000 {
        let n = sample_size(&mut rng, case, 50);
        let z = rng.uniform_f32(-100.0, 100.0);
        let mut cl = rng.vec_f32(n, -3.0, 8.0);
        let mut ct = rng.vec_f32(n, 0.0, 4.0);
        let mut rho = rng.vec_f32(n, -1.0, 4.0);
        let mut t: Vec<_> = (0..n)
            .map(|_| (rng.next_u32() % 100) as f32 / 100.0)
            .collect();
        let (mut cc, mut cs, mut cr, mut tt) = (cl.clone(), ct.clone(), rho.clone(), t.clone());
        refl::apply_earth_flattening(z, &mut cl, &mut ct, &mut rho, &mut t);
        unsafe {
            apply_earth_flattening(
                n as i32,
                z,
                cc.as_mut_ptr(),
                cs.as_mut_ptr(),
                cr.as_mut_ptr(),
                tt.as_mut_ptr(),
            )
        };
        for (a, b) in cl
            .iter()
            .chain(&ct)
            .chain(&rho)
            .chain(&t)
            .zip(cc.iter().chain(&cs).chain(&cr).chain(&tt))
        {
            close(*a, *b);
        }
    }
}
#[test]
fn hanning_window() {
    let mut rng = XorShift::new(29);
    for case in 0..10000 {
        let nf = sample_size(&mut rng, case, 200) as i32;
        let f1 = (rng.next_u32() % nf as u32 + 1) as i32;
        let f2 = f1 + (rng.next_u32() % (nf - f1 + 1) as u32) as i32;
        let iw = (rng.next_u32() % (nf + 5) as u32) as i32;
        let i = (rng.next_u32() % 3) as i32;
        let (mut a, mut b) = (Complex { r: 0.0, i: 0.0 }, Complex { r: 0.0, i: 0.0 });
        refl::compute_Hanning_window(i, iw, f1, f2, nf, &mut a);
        unsafe { compute_Hanning_window(i, iw, f1, f2, nf, &mut b) };
        close(a.r, b.r);
        close(a.i, b.i);
    }
}

fn z(r: &mut XorShift) -> Complex {
    Complex {
        r: (r.next_u32() % 2000) as f32 / 1000.0 + 0.1,
        i: (r.next_u32() % 2000) as f32 / 1000.0 - 1.0,
    }
}
fn closez(a: Complex, b: Complex) {
    let scale =
        a.r.abs()
            .max(a.i.abs())
            .max(b.r.abs())
            .max(b.i.abs())
            .max(1e-30);
    assert!((a.r - b.r).abs() <= 1e-5 * scale && (a.i - b.i).abs() <= 1e-5 * scale);
}
#[test]
fn prs() {
    let mut rng = XorShift::new(37);
    for i in 0..4000 {
        let n = sample_size(&mut rng, i, 31);
        let al: Vec<_> = (0..n).map(|_| z(&mut rng)).collect();
        let at: Vec<_> = (0..n)
            .map(|j| {
                let mut a = z(&mut rng);
                if (i + j) % 2 == 0 {
                    a.r = 0.0;
                }
                a
            })
            .collect();
        let obs: Vec<_> = (0..n).map(|_| rng.range_i32(1, n as i32)).collect();
        let rho = rng.vec_f32(n, -3.0, 4.0);
        let (mut a, mut b) = (vec![Complex::default(); n], vec![Complex::default(); n]);
        refl::compute_prs(&obs, &al, &at, &rho, &mut a);
        unsafe {
            compute_prs(
                n as i32,
                obs.as_ptr(),
                al.as_ptr(),
                at.as_ptr(),
                rho.as_ptr(),
                b.as_mut_ptr(),
            );
        }
        for j in 0..n {
            closez(a[j], b[j]);
        }
    }
}
#[test]
fn alpha_betha() {
    let mut rng = XorShift::new(42);
    for i in 0..4000 {
        let n = sample_size(&mut rng, i, 31);
        let ats = if i % 2 == 0 {
            Complex { r: 0.0, i: 0.0 }
        } else {
            z(&mut rng)
        };
        let gl: Vec<_> = (0..n * 2).map(|_| z(&mut rng)).collect();
        let gt: Vec<_> = (0..n * 2).map(|_| z(&mut rng)).collect();
        let (mut a, mut b, mut ca, mut cb) =
            (vec![ats; n], vec![ats; n], vec![ats; n], vec![ats; n]);
        refl::compute_alpha_betha(2, ats, 2.0, &gl, &gt, &mut a, &mut b);
        unsafe {
            compute_alpha_betha(
                2,
                n as i32,
                ats,
                2.0,
                gl.as_ptr(),
                gt.as_ptr(),
                ca.as_mut_ptr(),
                cb.as_mut_ptr(),
            );
        }
        for j in 0..n {
            closez(a[j], ca[j]);
            closez(b[j], cb[j]);
        }
    }
}
#[test]
fn red_velocity() {
    let mut rng = XorShift::new(43);
    for i in 0..4000 {
        let x = rng.uniform_f32(-10.0, 10.0);
        let vel = [0.0, 2.0, -2.0][i % 3];
        let lag = rng.uniform_f32(-1.0, 1.0);
        let wp = z(&mut rng);
        let ws = z(&mut rng);
        let (mut a, mut b) = (ws, ws);
        refl::red_vel_factor(x, vel, lag, wp, ws, &mut a);
        unsafe {
            red_vel_factor(x, vel, lag, wp, ws, &mut b);
        }
        closez(a, b);
    }
}

#[test]
fn moment_tensor() {
    let mut rng = XorShift::new(701);
    for i in 0..5000 {
        let mut vals: Vec<f32> = (0..5)
            .map(|_| (rng.next_u32() % 72000) as f32 / 100.0 - 360.0)
            .collect();
        if i % 11 == 0 {
            vals.fill(0.0);
        }
        let wtype = (i % 3) as i32;
        let (mut a, mut b, mut c) = (1.0, 2.0, 3.0);
        let (mut x, mut y, mut z) = (a, b, c);
        refl::compute_moment_tensor(
            wtype, vals[0], vals[1], vals[2], vals[3], vals[4], &mut a, &mut b, &mut c,
        );
        unsafe {
            compute_moment_tensor(
                wtype, vals[0], vals[1], vals[2], vals[3], vals[4], &mut x, &mut y, &mut z,
            );
        }
        close(a, x);
        close(b, y);
        close(c, z);
    }
}

#[test]
fn pwin() {
    let mut rng = XorShift::new(702);
    for i in 0..3000 {
        let n = sample_size(&mut rng, i as usize, 50);
        let bp = (rng.next_u32() % 100) as f32 / 100.0 - 0.2;
        let dp = (rng.next_u32() % 100 + 1) as f32 / 10000.0;
        let zero = Complex { r: 0.0, i: 0.0 };
        let (mut a, mut b, mut c) = (vec![zero; n], vec![zero; n], vec![zero; n]);
        let (mut x, mut y, mut z) = (a.clone(), b.clone(), c.clone());
        refl::compute_pwin(
            1 + i % 2,
            bp,
            dp,
            10.0,
            0.2,
            0.0,
            0.1,
            0.5,
            0.8,
            &mut a,
            &mut b,
            &mut c,
        );
        unsafe {
            compute_pwin(
                1 + i % 2,
                n as i32,
                bp,
                dp,
                10.0,
                0.2,
                0.0,
                0.1,
                0.5,
                0.8,
                x.as_mut_ptr(),
                y.as_mut_ptr(),
                z.as_mut_ptr(),
            );
        }
        for j in 0..n {
            closez(a[j], x[j]);
            closez(b[j], y[j]);
            closez(c[j], z[j]);
        }
    }
}

#[test]
fn gl_gt_gam() {
    let mut rng = XorShift::new(703);
    for i in 0..3000 {
        let n = sample_size(&mut rng, i, 20);
        let layers = (rng.next_u32() % 5 + 1) as usize;
        let al: Vec<_> = (0..layers).map(|_| z(&mut rng)).collect();
        let at: Vec<_> = (0..layers)
            .map(|j| {
                if (i + j) % 3 == 0 {
                    Complex { r: 0.0, i: 0.0 }
                } else {
                    z(&mut rng)
                }
            })
            .collect();
        let pp: Vec<_> = (0..n).map(|_| z(&mut rng)).collect();
        let rho = vec![2.0; layers];
        let zero = Complex { r: 0.0, i: 0.0 };
        let (mut a, mut b, mut c) = (
            vec![zero; n * layers],
            vec![zero; n * layers],
            vec![zero; n * layers],
        );
        let (mut x, mut y, mut z) = (a.clone(), b.clone(), c.clone());
        refl::compute_gl_gt_gam(
            1 + (i % 2) as i32,
            &al,
            &at,
            &rho,
            &pp,
            &mut a,
            &mut b,
            &mut c,
        );
        unsafe {
            compute_gl_gt_gam(
                1 + (i % 2) as i32,
                layers as i32,
                n as i32,
                al.as_ptr(),
                at.as_ptr(),
                rho.as_ptr(),
                pp.as_ptr(),
                x.as_mut_ptr(),
                y.as_mut_ptr(),
                z.as_mut_ptr(),
            );
        }
        for j in 0..n * layers {
            closez(a[j], x[j]);
            closez(b[j], y[j]);
            closez(c[j], z[j]);
        }
    }
}

#[test]
fn interpolation() {
    let mut rng = XorShift::new(704);
    for k in 0..3000 {
        let n = sample_size(&mut rng, k, 30);
        let nint = if n > 2 {
            (rng.next_u32() % (n - 2) as u32 + 1) as usize
        } else {
            0
        };
        let layers = [if n > 2 { 1 } else { -1 }, -1];
        let ns = [nint as i32, 0];
        let th = [if k % 2 == 0 { 2.0 } else { -4.0 }, 0.0];
        let mut values: Vec<Vec<f32>> = (0..6)
            .map(|_| {
                (0..n)
                    .map(|_| (rng.next_u32() % 1000) as f32 / 100.0 + 0.1)
                    .collect()
            })
            .collect();
        if k % 2 == 0 && n > 2 {
            values[2][nint + 1] = -1.0;
            values[3][nint + 1] = -1.0;
            values[4][nint] = -1.0;
        }
        let mut expected = values.clone();
        let [cl, ql, ct, qt, rho, t] = &mut values[..] else {
            unreachable!()
        };
        refl::parameter_interpolation(&layers, &ns, &th, cl, ql, ct, qt, rho, t);
        let [cl, ql, ct, qt, rho, t] = &mut expected[..] else {
            unreachable!()
        };
        unsafe {
            parameter_interpolation(
                n as i32,
                layers.as_ptr(),
                ns.as_ptr(),
                th.as_ptr(),
                cl.as_mut_ptr(),
                ql.as_mut_ptr(),
                ct.as_mut_ptr(),
                qt.as_mut_ptr(),
                rho.as_mut_ptr(),
                t.as_mut_ptr(),
            );
        }
        for (a, b) in values.iter().flatten().zip(expected.iter().flatten()) {
            close(*a, *b);
        }
    }
}

#[test]
fn slowness() {
    let mut rng = XorShift::new(705);
    for i in 0..5000 {
        let n = sample_size(&mut rng, i, 31).max(2);
        let cl = rng.vec_f32(n, 2.0, 8.0);
        let ct: Vec<_> = (0..n)
            .map(|j| {
                if (i + j) % 3 == 0 {
                    0.0
                } else {
                    rng.uniform_f32(1.0, 4.0)
                }
            })
            .collect();
        let t = rng.vec_f32(n, 0.01, 10.0);
        let ls = rng.range_i32(2, n as i32);
        let wtype = 1 + (i % 2) as i32;
        let w = 1.0 + (rng.next_u32() % 10000) as f32 / 100.0;
        let (mut np, mut fp, mut dp, mut p1, mut p2, mut p3, mut p4, mut c) =
            (0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Complex { r: 0.0, i: 0.0 });
        let (mut n2, mut f2, mut d2, mut q1, mut q2, mut q3, mut q4, mut c2) =
            (np, fp, dp, p1, p2, p3, p4, c);
        refl::compute_slowness(
            wtype, ls, w, 1.0, &mut np, &mut fp, 0.1, 0.4, 1.0, 0.2, 100.0, &mut dp, &mut p1,
            &mut p2, &mut p3, &mut p4, &mut c, &cl, &ct, &t,
        );
        unsafe {
            compute_slowness(
                wtype,
                ls,
                w,
                1.0,
                &mut n2,
                &mut f2,
                0.1,
                0.4,
                1.0,
                0.2,
                100.0,
                &mut d2,
                &mut q1,
                &mut q2,
                &mut q3,
                &mut q4,
                &mut c2,
                cl.as_ptr(),
                ct.as_ptr(),
                t.as_ptr(),
            );
        }
        assert_eq!(np, n2);
        for (a, b) in [(fp, f2), (dp, d2), (p1, q1), (p2, q2), (p3, q3), (p4, q4)] {
            close(a, b);
        }
        closez(c, c2);
    }
}
// Additional cases.

extern "C" {
    fn apply_filters(
        min_phase: *mut i32,
        nfilters: i32,
        nw: i32,
        tsec: f32,
        f1: *mut f32,
        f2: *mut f32,
        dbpo: *mut f32,
        filter_type: *mut i32,
        refw: *mut Complex,
    );
    fn compute_al_at(
        wtype: i32,
        nlayers: i32,
        layern: i32,
        eps: f32,
        epsp: f32,
        epss: f32,
        sigma: f32,
        sigp: f32,
        sigs: f32,
        wrefp: *mut f32,
        wrefs: *mut f32,
        fw: f32,
        cl: *mut f32,
        ct: *mut f32,
        ql: *mut f32,
        qt: *mut f32,
        rho: *mut f32,
        wpie: Complex,
        al: *mut Complex,
        at: *mut Complex,
    );
    fn compute_p_aux_arrays(
        wtype: i32,
        nlayers: i32,
        lsource: i32,
        block_size: i32,
        bp: f32,
        dp: f32,
        w: f32,
        decay: f32,
        pw1: f32,
        pw2: f32,
        pw3: f32,
        pw4: f32,
        pwin: *mut Complex,
        m1: f32,
        m2: f32,
        m3: f32,
        h1: f32,
        h2: f32,
        wpie: Complex,
        divfac: *mut Complex,
        p: *mut Complex,
        pp: *mut Complex,
        al: *mut Complex,
        at: *mut Complex,
        rho: *mut f32,
        gl: *mut Complex,
        gt: *mut Complex,
        gam: *mut Complex,
        alpha: *mut Complex,
        betha: *mut Complex,
        sigmad1: *mut Complex,
        sigmad2: *mut Complex,
        sigmau1: *mut Complex,
        sigmau2: *mut Complex,
    );
    fn compute_reflectivities(
        int_type: i32,
        verbose: i32,
        wtype: i32,
        wfield: i32,
        vsp: i32,
        flt: i32,
        win: i32,
        nx: i32,
        nt: i32,
        ntc: i32,
        nor: i32,
        nf: i32,
        nlayers: i32,
        lsource: i32,
        layern: i32,
        nfilters: i32,
        filters_phase: *mut i32,
        nw: i32,
        np: i32,
        bp: f32,
        tlag: f32,
        red_vel: f32,
        w1: f32,
        w2: f32,
        fx: f32,
        dx: f32,
        bx: f32,
        fs: f32,
        decay: f32,
        p2w: f32,
        tsec: f32,
        fw: f32,
        wrefp: f32,
        wrefs: f32,
        epsp: f32,
        epss: f32,
        sigp: f32,
        sigs: f32,
        pw1: f32,
        pw2: f32,
        pw3: f32,
        pw4: f32,
        h1: f32,
        h2: f32,
        m1: f32,
        m2: f32,
        m3: f32,
        fref: f32,
        lobs: *mut i32,
        filters_type: *mut i32,
        dbpo: *mut f32,
        f1: *mut f32,
        f2: *mut f32,
        cl: *mut f32,
        ct: *mut f32,
        ql: *mut f32,
        qt: *mut f32,
        rho: *mut f32,
        t: *mut f32,
        wavefield1: *mut *mut f32,
        wavefield2: *mut *mut f32,
        wavefield3: *mut *mut f32,
        outfp: *mut std::ffi::c_void,
    );
    fn compute_sigmas(
        wtype: i32,
        block_size: i32,
        lsource: i32,
        m1: f32,
        m2: f32,
        meu: Complex,
        alpha: *mut Complex,
        betha: *mut Complex,
        p: *mut Complex,
        gl: *mut Complex,
        gt: *mut Complex,
        gam: *mut Complex,
        a2: Complex,
        a3: Complex,
        pwin: *mut Complex,
        s1: Complex,
        s2: Complex,
        s4: Complex,
        sigmau1: *mut Complex,
        sigmau2: *mut Complex,
        sigmad1: *mut Complex,
        sigmad2: *mut Complex,
    );
    fn compute_synthetics(
        verbose: i32,
        nt: i32,
        ntc: i32,
        nx: i32,
        nor: i32,
        nw: i32,
        nlayers: i32,
        lsource: i32,
        nf: i32,
        flt: i32,
        vsp: i32,
        win: i32,
        wtype: i32,
        tlag: f32,
        red_vel: f32,
        decay: f32,
        tsec: f32,
        bx: f32,
        dx: f32,
        w1: f32,
        w2: f32,
        fil_phase: *mut i32,
        nfilters: i32,
        fil_type: *mut i32,
        dbpo: *mut f32,
        f1: *mut f32,
        f2: *mut f32,
        lobs: *mut i32,
        cl: *mut f32,
        t: *mut f32,
        response: *mut *mut *mut Complex,
        reflectivity: *mut *mut f32,
        outfp: *mut std::ffi::c_void,
    );
    fn compute_w_aux_arrays(
        wtype: i32,
        layern: i32,
        nlayers: i32,
        nor: i32,
        lsource: i32,
        np: *mut i32,
        block_size: *mut i32,
        nblock: *mut i32,
        left: *mut i32,
        fref: f32,
        wrefp: f32,
        wrefs: f32,
        lobs: *mut i32,
        tsec: f32,
        p2w: f32,
        fs: f32,
        xmax: f32,
        w: f32,
        decay: f32,
        epsp: f32,
        epss: f32,
        sigp: f32,
        sigs: f32,
        pw1: *mut f32,
        pw2: *mut f32,
        pw3: *mut f32,
        pw4: *mut f32,
        dp: *mut f32,
        fp: *mut f32,
        wpie: Complex,
        cdp: *mut Complex,
        cl: *mut f32,
        ct: *mut f32,
        ql: *mut f32,
        qt: *mut f32,
        rho: *mut f32,
        t: *mut f32,
        al: *mut Complex,
        at: *mut Complex,
        prs: *mut Complex,
    );
    fn construct_tx_trace(
        nw: i32,
        nt: i32,
        ntc: i32,
        ntfft: i32,
        nfilters: i32,
        filters_phase: *mut i32,
        tsec: f32,
        unexp: f32,
        sphrd: f32,
        refw: *mut Complex,
        fil_type: *mut i32,
        f1: *mut f32,
        f2: *mut f32,
        dbpo: *mut f32,
        reft: *mut f32,
    );
    fn convolve_wavelet(
        wavelet_type: i32,
        nx: i32,
        nt: i32,
        dt: f32,
        fpeak: f32,
        wfield: *mut *mut f32,
    );
    fn psv_reflectivities(
        int_type: i32,
        verbose: i32,
        wtype: i32,
        nw: i32,
        nlayers: i32,
        nx: i32,
        layern: i32,
        nor: i32,
        np: i32,
        bp1: f32,
        m1: f32,
        m2: f32,
        m3: f32,
        h1: f32,
        h2: f32,
        lsource: i32,
        bx: f32,
        dx: f32,
        xmax: f32,
        decay: f32,
        fref: f32,
        wrefp: f32,
        wrefs: f32,
        tsec: f32,
        p2w: f32,
        fs: f32,
        epsp: f32,
        epss: f32,
        sigp: f32,
        sigs: f32,
        pw1: f32,
        pw2: f32,
        pw3: f32,
        pw4: f32,
        acoustic: *mut i32,
        flag: *mut i32,
        lobs: *mut i32,
        rho: *mut f32,
        t: *mut f32,
        cl: *mut f32,
        ct: *mut f32,
        ql: *mut f32,
        qt: *mut f32,
        response1: *mut *mut *mut Complex,
        response2: *mut *mut *mut Complex,
        response3: *mut *mut *mut Complex,
        outfp: *mut std::ffi::c_void,
    );
    fn random_velocity_layers(
        nlayers: *mut i32,
        lsource: *mut i32,
        nrand_layers: i32,
        sdcl: f32,
        sdct: f32,
        layer: f32,
        zlayer: f32,
        cl: *mut f32,
        ql: *mut f32,
        ct: *mut f32,
        qt: *mut f32,
        rho: *mut f32,
        t: *mut f32,
    );
    fn sh_reflectivities(
        int_type: i32,
        verbose: i32,
        wtype: i32,
        nw: i32,
        nlayers: i32,
        nx: i32,
        layern: i32,
        nor: i32,
        np: i32,
        bp1: f32,
        m1: f32,
        m2: f32,
        m3: f32,
        h1: f32,
        h2: f32,
        lsource: i32,
        bx: f32,
        dx: f32,
        xmax: f32,
        decay: f32,
        fref: f32,
        wrefp: f32,
        wrefs: f32,
        tsec: f32,
        p2w: f32,
        fs: f32,
        epsp: f32,
        epss: f32,
        sigp: f32,
        sigs: f32,
        pw1: f32,
        pw2: f32,
        pw3: f32,
        pw4: f32,
        flag: *mut i32,
        lobs: *mut i32,
        rho: *mut f32,
        t: *mut f32,
        cl: *mut f32,
        ct: *mut f32,
        ql: *mut f32,
        qt: *mut f32,
        response1: *mut *mut *mut Complex,
        outfp: *mut std::ffi::c_void,
    );
}

// Serialize the tests which manipulate the C RNG or silence its unconditional
// progress messages. C diagnostics use fd 2, outside the Rust test capture.
static C_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
extern "C" {
    fn sranuni(seed: i32);
    fn alloc1complex(n: usize) -> *mut Complex;
    fn alloc1float(n: usize) -> *mut f32;
    fn free1complex(p: *mut Complex);
    fn free1float(p: *mut f32);
    fn dup(fd: i32) -> i32;
    fn dup2(src: i32, dst: i32) -> i32;
    fn fdopen(fd: i32, mode: *const std::ffi::c_char) -> *mut std::ffi::c_void;
    fn fclose(file: *mut std::ffi::c_void) -> i32;
    #[link_name = "close"]
    fn c_close(fd: i32) -> i32;
}
struct Quiet {
    fd: i32,
    _lock: std::sync::MutexGuard<'static, ()>,
}
impl Quiet {
    fn new() -> Self {
        use std::os::fd::AsRawFd;
        let lock = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let file = std::fs::OpenOptions::new()
            .write(true)
            .open("/dev/null")
            .unwrap();
        let fd = unsafe { dup(2) };
        assert!(fd >= 0);
        assert_eq!(unsafe { dup2(file.as_raw_fd(), 2) }, 2);
        Self { fd, _lock: lock }
    }
}
impl Drop for Quiet {
    fn drop(&mut self) {
        unsafe {
            dup2(self.fd, 2);
            c_close(self.fd);
        }
    }
}
struct COutput(*mut std::ffi::c_void);
impl COutput {
    fn new() -> Self {
        let file = unsafe { fdopen(dup(2), b"w\0".as_ptr().cast()) };
        assert!(!file.is_null());
        Self(file)
    }
}
impl Drop for COutput {
    fn drop(&mut self) {
        unsafe {
            fclose(self.0);
        }
    }
}

// C leaves the initial synthesis spectrum and the tail of the Larner wavelet
// unwritten. Prime allocator blocks with zeros so their defined computations
// can be compared to Rust's initialized workspace. Several blocks cover both
// scratch allocations even when they share a malloc size class.
unsafe fn prime_complex(n: usize) {
    let mut blocks = Vec::new();
    for _ in 0..8 {
        let p = alloc1complex(n);
        std::ptr::write_bytes(p, 0, n);
        blocks.push(p);
    }
    for p in blocks {
        free1complex(p);
    }
}
unsafe fn prime_float(n: usize) {
    let mut blocks = Vec::new();
    for _ in 0..8 {
        let p = alloc1float(n);
        std::ptr::write_bytes(p, 0, n);
        blocks.push(p);
    }
    for p in blocks {
        free1float(p);
    }
}
fn compare_floats(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    let scale = a
        .iter()
        .chain(b)
        .filter(|x| x.is_finite())
        .fold(1e-30f32, |s, x| s.max(x.abs()));
    for (j, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert!(
            x == y || (x.is_nan() && y.is_nan()) || (x - y).abs() <= 1e-5 * scale,
            "element {j}: Rust {x} C {y}, scale {scale}"
        );
    }
}
fn compare_complex(a: &[Complex], b: &[Complex]) {
    compare_floats(
        &a.iter().flat_map(|z| [z.r, z.i]).collect::<Vec<_>>(),
        &b.iter().flat_map(|z| [z.r, z.i]).collect::<Vec<_>>(),
    );
}
fn check_matrix(a: &[Vec<f32>], b: &[Vec<f32>]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        compare_floats(x, y);
    }
}
fn check_cube(a: &[Vec<Vec<Complex>>], b: &[Vec<Vec<Complex>>]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        for (u, v) in x.iter().zip(y) {
            compare_complex(u, v);
        }
    }
}

#[derive(Clone)]
struct Case {
    min_phase: Vec<i32>,
    nfilters: i32,
    nw: i32,
    tsec: f32,
    f1: Vec<f32>,
    f2: Vec<f32>,
    dbpo: Vec<f32>,
    filter_type: Vec<i32>,
    refw: Vec<Complex>,
    wtype: i32,
    nlayers: i32,
    layern: i32,
    eps: f32,
    epsp: f32,
    epss: f32,
    sigma: f32,
    sigp: f32,
    sigs: f32,
    wrefp: f32,
    wrefs: f32,
    fw: f32,
    cl: Vec<f32>,
    ct: Vec<f32>,
    ql: Vec<f32>,
    qt: Vec<f32>,
    rho: Vec<f32>,
    wpie: Complex,
    al: Vec<Complex>,
    at: Vec<Complex>,
    lsource: i32,
    block_size: i32,
    bp: f32,
    dp: f32,
    w: f32,
    decay: f32,
    pw1: f32,
    pw2: f32,
    pw3: f32,
    pw4: f32,
    pwin: Vec<Complex>,
    m1: f32,
    m2: f32,
    m3: f32,
    h1: f32,
    h2: f32,
    divfac: Complex,
    p: Vec<Complex>,
    pp: Vec<Complex>,
    gl: Vec<Complex>,
    gt: Vec<Complex>,
    gam: Vec<Complex>,
    alpha: Vec<Complex>,
    betha: Vec<Complex>,
    sigmad1: Vec<Complex>,
    sigmad2: Vec<Complex>,
    sigmau1: Vec<Complex>,
    sigmau2: Vec<Complex>,
    int_type: i32,
    verbose: i32,
    wfield: Vec<Vec<f32>>,
    vsp: i32,
    flt: i32,
    win: i32,
    nx: i32,
    nt: i32,
    ntc: i32,
    nor: i32,
    nf: i32,
    filters_phase: Vec<i32>,
    np: i32,
    tlag: f32,
    red_vel: f32,
    w1: f32,
    w2: f32,
    fx: f32,
    dx: f32,
    bx: f32,
    fs: f32,
    p2w: f32,
    fref: f32,
    lobs: Vec<i32>,
    filters_type: Vec<i32>,
    t: Vec<f32>,
    wavefield1: Vec<Vec<f32>>,
    wavefield2: Vec<Vec<f32>>,
    wavefield3: Vec<Vec<f32>>,
    meu: Complex,
    a2: Complex,
    a3: Complex,
    s1: Complex,
    s2: Complex,
    s4: Complex,
    fil_phase: Vec<i32>,
    fil_type: Vec<i32>,
    response: Vec<Vec<Vec<Complex>>>,
    reflectivity: Vec<Vec<f32>>,
    nblock: i32,
    left: i32,
    xmax: f32,
    fp: f32,
    cdp: Complex,
    prs: Vec<Complex>,
    ntfft: i32,
    unexp: f32,
    sphrd: f32,
    reft: Vec<f32>,
    wavelet_type: i32,
    dt: f32,
    fpeak: f32,
    bp1: f32,
    acoustic: Vec<i32>,
    flag: Vec<i32>,
    response1: Vec<Vec<Vec<Complex>>>,
    response2: Vec<Vec<Vec<Complex>>>,
    response3: Vec<Vec<Vec<Complex>>>,
    nrand_layers: i32,
    sdcl: f32,
    sdct: f32,
    layer: f32,
    zlayer: f32,
}

impl Case {
    fn new(fn_name: &str, iteration: usize, rng: &mut XorShift) -> Self {
        let sizes = [1usize, 2, 3, 7, 8, 11, 16, 31, 64, 127, 256, 1024];
        let points = if ["compute_p_aux_arrays", "compute_sigmas"].contains(&fn_name) {
            sizes[iteration % sizes.len()]
        } else {
            1
        };
        let nlayers = if fn_name == "compute_al_at" {
            sizes[iteration % sizes.len()]
        } else if ["compute_w_aux_arrays", "random_velocity_layers"].contains(&fn_name) {
            sizes[iteration % sizes.len()].max(2)
        } else {
            5 + rng.next_u32() as usize % 5
        };
        let nt = if [
            "convolve_wavelet",
            "construct_tx_trace",
            "compute_synthetics",
        ]
        .contains(&fn_name)
        {
            sizes[iteration % sizes.len()]
        } else if iteration % 32 == 31 {
            if iteration % 64 == 31 {
                127
            } else {
                256
            }
        } else {
            [2usize, 3, 4, 7, 8, 16, 32][iteration % 7]
        };
        let ntfft = oxwave_pfa_fft::npfa(nt as i32) as usize;
        let nw = if fn_name == "apply_filters" {
            [1usize, 2, 3, 4, 5, 7, 8, 11, 13, 16, 30, 64, 128][iteration % 13]
        } else if fn_name == "construct_tx_trace" {
            ntfft / 2
        } else {
            1 + iteration % 3
        };
        let nx = if fn_name == "compute_reflectivities" {
            ntfft
        } else {
            1 + iteration % 3
        };
        let nor = if fn_name == "psv_reflectivities" || fn_name == "compute_reflectivities" {
            1
        } else {
            1 + iteration % 3
        };
        let complex_size = if fn_name == "apply_filters" {
            2 * nw
        } else {
            ntfft
        };
        let zero = Complex { r: 0.0, i: 0.0 };
        let filter_count = if fn_name == "apply_filters" {
            1 + iteration % 3
        } else if ntfft >= 2
            && [
                "construct_tx_trace",
                "compute_synthetics",
                "compute_reflectivities",
            ]
            .contains(&fn_name)
        {
            iteration % 3
        } else {
            0
        };
        let al = (0..nlayers)
            .map(|_| Complex {
                r: rng.uniform_f32(0.1, 0.5),
                i: rng.uniform_f32(0.001, 0.04),
            })
            .collect();
        let at = (0..nlayers)
            .map(|j| {
                if iteration % 4 == 0 && j == 0 {
                    zero
                } else {
                    Complex {
                        r: rng.uniform_f32(0.3, 0.8),
                        i: rng.uniform_f32(0.001, 0.03),
                    }
                }
            })
            .collect();
        Self {
            min_phase: vec![0; filter_count],
            nfilters: filter_count as i32,
            nw: nw as i32,
            tsec: 2.0,
            f1: (0..filter_count)
                .map(|_| rng.uniform_f32(3.0, 20.0))
                .collect(),
            f2: vec![25.0; filter_count],
            dbpo: (0..filter_count)
                .map(|_| rng.uniform_f32(6.0, 24.0))
                .collect(),
            filter_type: (0..filter_count)
                .map(|j| 1 + ((iteration + j) % 2) as i32)
                .collect(),
            refw: (0..complex_size).map(|_| z(rng)).collect(),
            wtype: if fn_name == "sh_reflectivities" {
                2
            } else if fn_name == "psv_reflectivities" {
                1
            } else {
                (1 + iteration % 2) as i32
            },
            nlayers: nlayers as i32,
            layern: if iteration % 3 == 0 { 0 } else { 2 },
            eps: 0.001,
            epsp: rng.uniform_f32(0.0001, 0.1),
            epss: rng.uniform_f32(0.0001, 0.1),
            sigma: 0.1,
            sigp: rng.uniform_f32(0.05, 0.5),
            sigs: rng.uniform_f32(0.05, 0.5),
            wrefp: 1.0,
            wrefs: 1.0,
            fw: 1.0,
            cl: (0..nlayers).map(|_| rng.uniform_f32(3.0, 7.0)).collect(),
            ct: (0..nlayers)
                .map(|j| {
                    if iteration % 4 == 0 && j == 0 && fn_name != "sh_reflectivities" {
                        0.0
                    } else {
                        rng.uniform_f32(1.4, 2.8)
                    }
                })
                .collect(),
            ql: rng.vec_f32(nlayers, 20.0, 200.0),
            qt: rng.vec_f32(nlayers, 20.0, 200.0),
            rho: rng.vec_f32(nlayers, 1.5, 3.5),
            wpie: Complex {
                r: rng.uniform_f32(1.0, 100.0),
                i: 0.1,
            },
            al: al,
            at: at,
            lsource: 2,
            block_size: points as i32,
            bp: rng.uniform_f32(-0.02, 0.03),
            dp: rng.uniform_f32(0.0001, 0.002),
            w: rng.uniform_f32(1.0, 100.0),
            decay: 0.1,
            pw1: 0.0,
            pw2: 0.005,
            pw3: 0.03,
            pw4: 0.04,
            pwin: (0..points).map(|_| z(rng)).collect(),
            m1: rng.uniform_f32(-2.0, 2.0),
            m2: rng.uniform_f32(-2.0, 2.0),
            m3: rng.uniform_f32(-2.0, 2.0),
            h1: rng.uniform_f32(-1.0, 1.0),
            h2: rng.uniform_f32(-1.0, 1.0),
            divfac: zero,
            p: (0..points)
                .map(|j| if j == 0 { zero } else { z(rng) })
                .collect(),
            pp: vec![zero; points],
            gl: (0..points * nlayers).map(|_| z(rng)).collect(),
            gt: (0..points * nlayers)
                .map(|j| {
                    if iteration % 3 == 0 && j % 2 == 0 {
                        zero
                    } else {
                        z(rng)
                    }
                })
                .collect(),
            gam: (0..points * nlayers).map(|_| z(rng)).collect(),
            alpha: (0..points).map(|_| z(rng)).collect(),
            betha: (0..points).map(|_| z(rng)).collect(),
            sigmad1: vec![zero; points],
            sigmad2: vec![zero; points],
            sigmau1: vec![zero; points],
            sigmau2: vec![zero; points],
            int_type: (1 + iteration % 2) as i32,
            verbose: 0,
            wfield: (0..nx).map(|_| rng.vec_f32(nt, -1.0, 1.0)).collect(),
            vsp: (iteration % 2) as i32,
            flt: (iteration % 2) as i32,
            win: (iteration % 2) as i32,
            nx: nx as i32,
            nt: nt as i32,
            ntc: ntfft as i32,
            nor: nor as i32,
            nf: nw as i32,
            filters_phase: vec![0; filter_count],
            np: 0,
            tlag: rng.uniform_f32(-0.1, 0.1),
            red_vel: if iteration % 3 == 0 {
                0.0
            } else if iteration % 3 == 1 {
                -1.0
            } else {
                3.0
            },
            w1: 0.0,
            w2: 0.0,
            fx: 1.0,
            dx: 0.1,
            bx: rng.uniform_f32(0.05, 0.2),
            fs: 0.1,
            p2w: rng.uniform_f32(0.008, 0.04),
            fref: 1.0,
            lobs: (0..nor)
                .map(|j| (2 + j as i32).min(nlayers as i32))
                .collect(),
            filters_type: (0..filter_count)
                .map(|j| 1 + ((iteration + j) % 2) as i32)
                .collect(),
            t: (0..nlayers)
                .map(|j| {
                    if j == nlayers - 1 {
                        10000.0
                    } else {
                        rng.uniform_f32(0.1, 0.9)
                    }
                })
                .collect(),
            wavefield1: vec![vec![0.0; nt]; nx],
            wavefield2: vec![vec![0.0; nt]; nx],
            wavefield3: vec![vec![0.0; nt]; nx],
            meu: z(rng),
            a2: z(rng),
            a3: z(rng),
            s1: z(rng),
            s2: z(rng),
            s4: z(rng),
            fil_phase: vec![0; filter_count],
            fil_type: (0..filter_count)
                .map(|j| 1 + ((iteration + j) % 2) as i32)
                .collect(),
            response: (0..nw)
                .map(|_| {
                    (0..nx)
                        .map(|_| (0..nor).map(|_| z(rng)).collect())
                        .collect()
                })
                .collect(),
            reflectivity: vec![vec![0.0; nt]; nx],
            nblock: 0,
            left: 0,
            xmax: 10.0,
            fp: 0.0,
            cdp: zero,
            prs: vec![zero; nor],
            ntfft: ntfft as i32,
            unexp: rng.uniform_f32(-0.5, 0.5),
            sphrd: rng.uniform_f32(0.5, 2.0),
            reft: vec![0.0; nt],
            wavelet_type: if nt == 1 {
                2 + 2 * ((iteration / 12) % 2) as i32
            } else {
                (1 + iteration % 4) as i32
            },
            dt: 0.004,
            fpeak: rng.uniform_f32(5.0, 80.0),
            bp1: rng.uniform_f32(-0.01, 0.02),
            acoustic: vec![0; nor],
            flag: vec![0; nor],
            response1: (0..nw)
                .map(|_| {
                    (0..nx)
                        .map(|_| (0..nor).map(|_| z(rng)).collect())
                        .collect()
                })
                .collect(),
            response2: (0..nw)
                .map(|_| {
                    (0..nx)
                        .map(|_| (0..nor).map(|_| z(rng)).collect())
                        .collect()
                })
                .collect(),
            response3: (0..nw)
                .map(|_| {
                    (0..nx)
                        .map(|_| (0..nor).map(|_| z(rng)).collect())
                        .collect()
                })
                .collect(),
            nrand_layers: 2,
            sdcl: rng.uniform_f32(-10.0, 10.0),
            sdct: rng.uniform_f32(-10.0, 10.0),
            layer: 1.0,
            zlayer: 1.0,
        }
    }
}

#[test]
fn diff_apply_filters() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1000);
    for iteration in 0..2048 {
        let mut a = Case::new("apply_filters", iteration, &mut rng);

        // Exercise both minimum phase paths and all three filter designs.
        if iteration % 4 == 0 && oxwave_pfa_fft::npfa(2 * a.nw) == 2 * a.nw {
            a.min_phase[0] = 1;
            a.f1.fill(10.0);
        }
        if a.nw >= 30 && iteration % 3 == 0 {
            let delta = if a.min_phase[0] == 1 {
                (2.0 * std::f64::consts::PI / a.tsec as f64) as f32
            } else {
                0.0
            };
            a.filter_type[0] = 3;
            a.f1[0] = 1.0 + delta;
            a.f2[0] = (a.nw - 3) as f32 / a.tsec - delta;
        }
        let mut b = a.clone();
        refl::apply_filters(
            &a.min_phase,
            a.nw,
            a.tsec,
            &a.f1,
            &a.f2,
            &a.dbpo,
            &a.filter_type,
            &mut a.refw,
        );
        unsafe {
            apply_filters(
                b.min_phase.as_mut_ptr(),
                b.nfilters,
                b.nw,
                b.tsec,
                b.f1.as_mut_ptr(),
                b.f2.as_mut_ptr(),
                b.dbpo.as_mut_ptr(),
                b.filter_type.as_mut_ptr(),
                b.refw.as_mut_ptr(),
            );
        }
        compare_complex(&a.refw, &b.refw);
    }
}

#[test]
fn diff_compute_al_at() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1001);
    for iteration in 0..2048 {
        let mut a = Case::new("compute_al_at", iteration, &mut rng);
        if iteration % 7 == 0 {
            a.ql[0] = -a.ql[0];
            a.qt[0] = -a.qt[0];
        }
        let mut b = a.clone();
        refl::compute_al_at(
            a.wtype,
            a.layern,
            a.eps,
            a.epsp,
            a.epss,
            a.sigma,
            a.sigp,
            a.sigs,
            &mut a.wrefp,
            &mut a.wrefs,
            a.fw,
            &a.cl,
            &a.ct,
            &a.ql,
            &a.qt,
            &a.rho,
            a.wpie,
            &mut a.al,
            &mut a.at,
        );
        unsafe {
            compute_al_at(
                b.wtype,
                b.nlayers,
                b.layern,
                b.eps,
                b.epsp,
                b.epss,
                b.sigma,
                b.sigp,
                b.sigs,
                &mut b.wrefp,
                &mut b.wrefs,
                b.fw,
                b.cl.as_mut_ptr(),
                b.ct.as_mut_ptr(),
                b.ql.as_mut_ptr(),
                b.qt.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.wpie,
                b.al.as_mut_ptr(),
                b.at.as_mut_ptr(),
            );
        }
        compare_complex(&a.at, &b.at);
        compare_complex(&a.al, &b.al);
        close(a.wrefp, b.wrefp);
        close(a.wrefs, b.wrefs);
    }
}

#[test]
fn diff_compute_p_aux_arrays() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1002);
    for iteration in 0..2048 {
        let mut a = Case::new("compute_p_aux_arrays", iteration, &mut rng);
        let mut b = a.clone();
        refl::compute_p_aux_arrays(
            a.wtype,
            a.lsource,
            a.bp,
            a.dp,
            a.w,
            a.decay,
            a.pw1,
            a.pw2,
            a.pw3,
            a.pw4,
            &mut a.pwin,
            a.m1,
            a.m2,
            a.m3,
            a.h1,
            a.h2,
            a.wpie,
            &mut a.divfac,
            &mut a.p,
            &mut a.pp,
            &a.al,
            &a.at,
            &a.rho,
            &mut a.gl,
            &mut a.gt,
            &mut a.gam,
            &mut a.alpha,
            &mut a.betha,
            &mut a.sigmad1,
            &mut a.sigmad2,
            &mut a.sigmau1,
            &mut a.sigmau2,
        );
        unsafe {
            compute_p_aux_arrays(
                b.wtype,
                b.nlayers,
                b.lsource,
                b.block_size,
                b.bp,
                b.dp,
                b.w,
                b.decay,
                b.pw1,
                b.pw2,
                b.pw3,
                b.pw4,
                b.pwin.as_mut_ptr(),
                b.m1,
                b.m2,
                b.m3,
                b.h1,
                b.h2,
                b.wpie,
                &mut b.divfac,
                b.p.as_mut_ptr(),
                b.pp.as_mut_ptr(),
                b.al.as_mut_ptr(),
                b.at.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.gl.as_mut_ptr(),
                b.gt.as_mut_ptr(),
                b.gam.as_mut_ptr(),
                b.alpha.as_mut_ptr(),
                b.betha.as_mut_ptr(),
                b.sigmad1.as_mut_ptr(),
                b.sigmad2.as_mut_ptr(),
                b.sigmau1.as_mut_ptr(),
                b.sigmau2.as_mut_ptr(),
            );
        }
        compare_complex(&a.sigmau1, &b.sigmau1);
        compare_complex(&a.gam, &b.gam);
        compare_complex(&a.alpha, &b.alpha);
        compare_complex(&a.pwin, &b.pwin);
        compare_complex(&a.pp, &b.pp);
        compare_complex(&a.gl, &b.gl);
        compare_complex(&a.gt, &b.gt);
        compare_complex(&a.sigmau2, &b.sigmau2);
        compare_complex(&a.sigmad2, &b.sigmad2);
        compare_complex(&a.p, &b.p);
        compare_complex(&a.betha, &b.betha);
        compare_complex(&a.sigmad1, &b.sigmad1);
        closez(a.divfac, b.divfac);
    }
}

#[test]
fn diff_compute_reflectivities() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1003);
    for iteration in 0..2048 {
        let mut a = Case::new("compute_reflectivities", iteration, &mut rng);
        let output = COutput::new();
        a.verbose = (iteration % 4) as i32;
        a.lsource = 2 + (iteration % (a.nlayers as usize - 1)) as i32;
        refl::source_receiver_type(
            a.lsource,
            &a.lobs,
            &a.cl,
            &a.ct,
            &mut a.acoustic,
            &mut a.flag,
        );

        // Include many 100-slowness blocks, their exact boundary, and a remainder.
        if iteration % 8 == 0 {
            a.p2w = 0.202;
        } else if iteration % 8 == 1 {
            a.p2w = 0.4;
        }
        if iteration % 13 == 0 {
            a.m1 = 0.0;
            a.m2 = 0.0;
            a.m3 = 0.0;
            a.h1 = 0.0;
            a.h2 = 0.0;
        }

        if a.wtype == 1 {
            match iteration % 8 {
                0 => {
                    a.ct.fill(0.0);
                    a.lsource = 2;
                } // liquid/liquid
                2 => {
                    a.ct[1] = 0.0;
                    a.ct[0] = 0.0;
                    a.lsource = 2;
                    a.lobs[0] = 3;
                } // liquid/solid
                4 => {
                    a.ct[2] = 0.0;
                    a.lsource = 2;
                    a.lobs[0] = 3;
                } // solid/liquid
                6 => {
                    a.cl[2] = a.cl[1];
                    a.ct[2] = a.ct[1];
                    a.ql[2] = a.ql[1];
                    a.qt[2] = a.qt[1];
                } // equal velocities
                _ => {}
            }
        }
        refl::source_receiver_type(
            a.lsource,
            &a.lobs,
            &a.cl,
            &a.ct,
            &mut a.acoustic,
            &mut a.flag,
        );
        let mut b = a.clone();
        b.lobs.push(i32::MAX); // readable end sentinel for the C receiver scan
        refl::compute_reflectivities(
            a.int_type,
            a.verbose,
            a.wtype,
            a.wtype,
            a.vsp,
            a.flt,
            a.win,
            a.ntc,
            a.nf,
            a.lsource,
            a.layern,
            &a.filters_phase,
            a.nw,
            a.np,
            a.bp,
            a.tlag,
            a.red_vel,
            a.w1,
            a.w2,
            a.fx,
            a.dx,
            a.bx,
            a.fs,
            a.decay,
            a.p2w,
            a.tsec,
            a.fw,
            a.wrefp,
            a.wrefs,
            a.epsp,
            a.epss,
            a.sigp,
            a.sigs,
            a.pw1,
            a.pw2,
            a.pw3,
            a.pw4,
            a.h1,
            a.h2,
            a.m1,
            a.m2,
            a.m3,
            a.fref,
            &a.lobs,
            &a.filters_type,
            &a.dbpo,
            &a.f1,
            &a.f2,
            &a.cl,
            &a.ct,
            &a.ql,
            &a.qt,
            &a.rho,
            &a.t,
            &mut a.wavefield1,
            &mut a.wavefield2,
            &mut a.wavefield3,
            &oxwave_core::prelude::stderr(),
        );
        let mut wavefield1_rows: Vec<_> = b.wavefield1.iter_mut().map(|r| r.as_mut_ptr()).collect();
        let mut wavefield2_rows: Vec<_> = b.wavefield2.iter_mut().map(|r| r.as_mut_ptr()).collect();
        let mut wavefield3_rows: Vec<_> = b.wavefield3.iter_mut().map(|r| r.as_mut_ptr()).collect();
        unsafe {
            prime_complex(b.ntc as usize);
        }
        unsafe {
            compute_reflectivities(
                b.int_type,
                b.verbose,
                b.wtype,
                b.wtype,
                b.vsp,
                b.flt,
                b.win,
                b.nx,
                b.nt,
                b.ntc,
                b.nor,
                b.nf,
                b.nlayers,
                b.lsource,
                b.layern,
                b.nfilters,
                b.filters_phase.as_mut_ptr(),
                b.nw,
                b.np,
                b.bp,
                b.tlag,
                b.red_vel,
                b.w1,
                b.w2,
                b.fx,
                b.dx,
                b.bx,
                b.fs,
                b.decay,
                b.p2w,
                b.tsec,
                b.fw,
                b.wrefp,
                b.wrefs,
                b.epsp,
                b.epss,
                b.sigp,
                b.sigs,
                b.pw1,
                b.pw2,
                b.pw3,
                b.pw4,
                b.h1,
                b.h2,
                b.m1,
                b.m2,
                b.m3,
                b.fref,
                b.lobs.as_mut_ptr(),
                b.filters_type.as_mut_ptr(),
                b.dbpo.as_mut_ptr(),
                b.f1.as_mut_ptr(),
                b.f2.as_mut_ptr(),
                b.cl.as_mut_ptr(),
                b.ct.as_mut_ptr(),
                b.ql.as_mut_ptr(),
                b.qt.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.t.as_mut_ptr(),
                wavefield1_rows.as_mut_ptr(),
                wavefield2_rows.as_mut_ptr(),
                wavefield3_rows.as_mut_ptr(),
                output.0,
            );
        }
        check_matrix(&a.wavefield3, &b.wavefield3);
        check_matrix(&a.wavefield2, &b.wavefield2);
        check_matrix(&a.wavefield1, &b.wavefield1);
    }
}

#[test]
fn diff_compute_sigmas() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1004);
    for iteration in 0..2048 {
        let mut a = Case::new("compute_sigmas", iteration, &mut rng);
        let mut b = a.clone();
        refl::compute_sigmas(
            a.wtype,
            a.lsource,
            a.m1,
            a.m2,
            a.meu,
            &a.alpha,
            &a.betha,
            &mut a.p,
            &a.gl,
            &mut a.gt,
            &a.gam,
            a.a2,
            a.a3,
            &a.pwin,
            a.s1,
            a.s2,
            a.s4,
            &mut a.sigmau1,
            &mut a.sigmau2,
            &mut a.sigmad1,
            &mut a.sigmad2,
        );
        unsafe {
            compute_sigmas(
                b.wtype,
                b.block_size,
                b.lsource,
                b.m1,
                b.m2,
                b.meu,
                b.alpha.as_mut_ptr(),
                b.betha.as_mut_ptr(),
                b.p.as_mut_ptr(),
                b.gl.as_mut_ptr(),
                b.gt.as_mut_ptr(),
                b.gam.as_mut_ptr(),
                b.a2,
                b.a3,
                b.pwin.as_mut_ptr(),
                b.s1,
                b.s2,
                b.s4,
                b.sigmau1.as_mut_ptr(),
                b.sigmau2.as_mut_ptr(),
                b.sigmad1.as_mut_ptr(),
                b.sigmad2.as_mut_ptr(),
            );
        }
        compare_complex(&a.sigmau1, &b.sigmau1);
        compare_complex(&a.sigmau2, &b.sigmau2);
        compare_complex(&a.sigmad2, &b.sigmad2);
        compare_complex(&a.gt, &b.gt);
        compare_complex(&a.p, &b.p);
        compare_complex(&a.sigmad1, &b.sigmad1);
    }
}

#[test]
fn diff_compute_synthetics() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1005);
    for iteration in 0..2048 {
        let mut a = Case::new("compute_synthetics", iteration, &mut rng);
        let output = COutput::new();
        a.verbose = (iteration % 4) as i32;
        a.wtype = (iteration % 4) as i32;
        let mut b = a.clone();
        refl::compute_synthetics(
            a.verbose,
            a.ntc,
            a.lsource,
            a.nf,
            a.flt,
            a.vsp,
            a.win,
            a.wtype,
            a.tlag,
            a.red_vel,
            a.decay,
            a.tsec,
            a.bx,
            a.dx,
            a.w1,
            a.w2,
            &a.fil_phase,
            &a.fil_type,
            &a.dbpo,
            &a.f1,
            &a.f2,
            &a.lobs,
            &a.cl,
            &a.t,
            &mut a.response,
            &mut a.reflectivity,
            &oxwave_core::prelude::stderr(),
        );
        let mut response_rows: Vec<Vec<_>> = b
            .response
            .iter_mut()
            .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
            .collect();
        let mut response_planes: Vec<_> =
            response_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
        let mut reflectivity_rows: Vec<_> =
            b.reflectivity.iter_mut().map(|r| r.as_mut_ptr()).collect();
        unsafe {
            prime_complex(b.ntc as usize);
        }
        unsafe {
            compute_synthetics(
                b.verbose,
                b.nt,
                b.ntc,
                b.nx,
                b.nor,
                b.nw,
                b.nlayers,
                b.lsource,
                b.nf,
                b.flt,
                b.vsp,
                b.win,
                b.wtype,
                b.tlag,
                b.red_vel,
                b.decay,
                b.tsec,
                b.bx,
                b.dx,
                b.w1,
                b.w2,
                b.fil_phase.as_mut_ptr(),
                b.nfilters,
                b.fil_type.as_mut_ptr(),
                b.dbpo.as_mut_ptr(),
                b.f1.as_mut_ptr(),
                b.f2.as_mut_ptr(),
                b.lobs.as_mut_ptr(),
                b.cl.as_mut_ptr(),
                b.t.as_mut_ptr(),
                response_planes.as_mut_ptr(),
                reflectivity_rows.as_mut_ptr(),
                output.0,
            );
        }
        check_matrix(&a.reflectivity, &b.reflectivity);
        check_cube(&a.response, &b.response);
    }
}

#[test]
fn diff_compute_w_aux_arrays() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1006);
    for iteration in 0..2048 {
        let mut a = Case::new("compute_w_aux_arrays", iteration, &mut rng);
        let mut b = a.clone();
        refl::compute_w_aux_arrays(
            a.wtype,
            a.layern,
            a.lsource,
            &mut a.np,
            &mut a.block_size,
            &mut a.nblock,
            &mut a.left,
            a.fref,
            a.wrefp,
            a.wrefs,
            &a.lobs,
            a.tsec,
            a.p2w,
            a.fs,
            a.xmax,
            a.w,
            a.decay,
            a.epsp,
            a.epss,
            a.sigp,
            a.sigs,
            &mut a.pw1,
            &mut a.pw2,
            &mut a.pw3,
            &mut a.pw4,
            &mut a.dp,
            &mut a.fp,
            a.wpie,
            &mut a.cdp,
            &a.cl,
            &a.ct,
            &a.ql,
            &a.qt,
            &a.rho,
            &a.t,
            &mut a.al,
            &mut a.at,
            &mut a.prs,
        );
        unsafe {
            compute_w_aux_arrays(
                b.wtype,
                b.layern,
                b.nlayers,
                b.nor,
                b.lsource,
                &mut b.np,
                &mut b.block_size,
                &mut b.nblock,
                &mut b.left,
                b.fref,
                b.wrefp,
                b.wrefs,
                b.lobs.as_mut_ptr(),
                b.tsec,
                b.p2w,
                b.fs,
                b.xmax,
                b.w,
                b.decay,
                b.epsp,
                b.epss,
                b.sigp,
                b.sigs,
                &mut b.pw1,
                &mut b.pw2,
                &mut b.pw3,
                &mut b.pw4,
                &mut b.dp,
                &mut b.fp,
                b.wpie,
                &mut b.cdp,
                b.cl.as_mut_ptr(),
                b.ct.as_mut_ptr(),
                b.ql.as_mut_ptr(),
                b.qt.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.t.as_mut_ptr(),
                b.al.as_mut_ptr(),
                b.at.as_mut_ptr(),
                b.prs.as_mut_ptr(),
            );
        }
        compare_complex(&a.prs, &b.prs);
        compare_complex(&a.al, &b.al);
        compare_complex(&a.at, &b.at);
        assert_eq!(a.np, b.np);
        assert_eq!(a.block_size, b.block_size);
        assert_eq!(a.nblock, b.nblock);
        assert_eq!(a.left, b.left);
        close(a.pw1, b.pw1);
        close(a.pw2, b.pw2);
        close(a.pw3, b.pw3);
        close(a.pw4, b.pw4);
        close(a.dp, b.dp);
        close(a.fp, b.fp);
        closez(a.cdp, b.cdp);
    }
}

#[test]
fn diff_construct_tx_trace() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1007);
    for iteration in 0..2048 {
        let mut a = Case::new("construct_tx_trace", iteration, &mut rng);
        let mut b = a.clone();
        refl::construct_tx_trace(
            a.nw,
            a.ntc,
            &a.filters_phase,
            a.tsec,
            a.unexp,
            a.sphrd,
            &mut a.refw,
            &a.fil_type,
            &a.f1,
            &a.f2,
            &a.dbpo,
            &mut a.reft,
        );
        unsafe {
            construct_tx_trace(
                b.nw,
                b.nt,
                b.ntc,
                b.ntfft,
                b.nfilters,
                b.filters_phase.as_mut_ptr(),
                b.tsec,
                b.unexp,
                b.sphrd,
                b.refw.as_mut_ptr(),
                b.fil_type.as_mut_ptr(),
                b.f1.as_mut_ptr(),
                b.f2.as_mut_ptr(),
                b.dbpo.as_mut_ptr(),
                b.reft.as_mut_ptr(),
            );
        }
        compare_complex(&a.refw, &b.refw);
        compare_floats(&a.reft, &b.reft);
    }
}

#[test]
fn diff_convolve_wavelet() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1008);
    for iteration in 0..2048 {
        let mut a = Case::new("convolve_wavelet", iteration, &mut rng);
        let mut b = a.clone();
        refl::convolve_wavelet(a.wavelet_type, a.dt, a.fpeak, &mut a.wfield);
        let mut wfield_rows: Vec<_> = b.wfield.iter_mut().map(|r| r.as_mut_ptr()).collect();
        unsafe {
            prime_float(b.nt as usize);
        }
        unsafe {
            convolve_wavelet(
                b.wavelet_type,
                b.nx,
                b.nt,
                b.dt,
                b.fpeak,
                wfield_rows.as_mut_ptr(),
            );
        }
        check_matrix(&a.wfield, &b.wfield);
    }
}

#[test]
fn diff_psv_reflectivities() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1009);
    for iteration in 0..2048 {
        let mut a = Case::new("psv_reflectivities", iteration, &mut rng);
        let output = COutput::new();
        a.verbose = (iteration % 4) as i32;
        a.lsource = 2 + (iteration % (a.nlayers as usize - 1)) as i32;
        refl::source_receiver_type(
            a.lsource,
            &a.lobs,
            &a.cl,
            &a.ct,
            &mut a.acoustic,
            &mut a.flag,
        );

        // Include many 100-slowness blocks, their exact boundary, and a remainder.
        if iteration % 8 == 0 {
            a.p2w = 0.202;
        } else if iteration % 8 == 1 {
            a.p2w = 0.4;
        }
        if iteration % 13 == 0 {
            a.m1 = 0.0;
            a.m2 = 0.0;
            a.m3 = 0.0;
            a.h1 = 0.0;
            a.h2 = 0.0;
        }

        if a.wtype == 1 {
            match iteration % 8 {
                0 => {
                    a.ct.fill(0.0);
                    a.lsource = 2;
                } // liquid/liquid
                2 => {
                    a.ct[1] = 0.0;
                    a.ct[0] = 0.0;
                    a.lsource = 2;
                    a.lobs[0] = 3;
                } // liquid/solid
                4 => {
                    a.ct[2] = 0.0;
                    a.lsource = 2;
                    a.lobs[0] = 3;
                } // solid/liquid
                6 => {
                    a.cl[2] = a.cl[1];
                    a.ct[2] = a.ct[1];
                    a.ql[2] = a.ql[1];
                    a.qt[2] = a.qt[1];
                } // equal velocities
                _ => {}
            }
        }
        refl::source_receiver_type(
            a.lsource,
            &a.lobs,
            &a.cl,
            &a.ct,
            &mut a.acoustic,
            &mut a.flag,
        );
        let mut b = a.clone();
        b.lobs.push(i32::MAX); // readable end sentinel for the C receiver scan
        refl::psv_reflectivities(
            a.int_type,
            a.verbose,
            a.wtype,
            a.layern,
            a.np,
            a.bp1,
            a.m1,
            a.m2,
            a.m3,
            a.h1,
            a.h2,
            a.lsource,
            a.bx,
            a.dx,
            a.xmax,
            a.decay,
            a.fref,
            a.wrefp,
            a.wrefs,
            a.tsec,
            a.p2w,
            a.fs,
            a.epsp,
            a.epss,
            a.sigp,
            a.sigs,
            a.pw1,
            a.pw2,
            a.pw3,
            a.pw4,
            &a.acoustic,
            &mut a.flag,
            &a.lobs,
            &a.rho,
            &a.t,
            &a.cl,
            &a.ct,
            &a.ql,
            &a.qt,
            &mut a.response1,
            &mut a.response2,
            &mut a.response3,
            &oxwave_core::prelude::stderr(),
        );
        let mut response1_rows: Vec<Vec<_>> = b
            .response1
            .iter_mut()
            .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
            .collect();
        let mut response1_planes: Vec<_> =
            response1_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
        let mut response2_rows: Vec<Vec<_>> = b
            .response2
            .iter_mut()
            .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
            .collect();
        let mut response2_planes: Vec<_> =
            response2_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
        let mut response3_rows: Vec<Vec<_>> = b
            .response3
            .iter_mut()
            .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
            .collect();
        let mut response3_planes: Vec<_> =
            response3_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
        unsafe {
            psv_reflectivities(
                b.int_type,
                b.verbose,
                b.wtype,
                b.nw,
                b.nlayers,
                b.nx,
                b.layern,
                b.nor,
                b.np,
                b.bp1,
                b.m1,
                b.m2,
                b.m3,
                b.h1,
                b.h2,
                b.lsource,
                b.bx,
                b.dx,
                b.xmax,
                b.decay,
                b.fref,
                b.wrefp,
                b.wrefs,
                b.tsec,
                b.p2w,
                b.fs,
                b.epsp,
                b.epss,
                b.sigp,
                b.sigs,
                b.pw1,
                b.pw2,
                b.pw3,
                b.pw4,
                b.acoustic.as_mut_ptr(),
                b.flag.as_mut_ptr(),
                b.lobs.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.t.as_mut_ptr(),
                b.cl.as_mut_ptr(),
                b.ct.as_mut_ptr(),
                b.ql.as_mut_ptr(),
                b.qt.as_mut_ptr(),
                response1_planes.as_mut_ptr(),
                response2_planes.as_mut_ptr(),
                response3_planes.as_mut_ptr(),
                output.0,
            );
        }
        assert_eq!(a.flag, b.flag);
        check_cube(&a.response2, &b.response2);
        check_cube(&a.response3, &b.response3);
        check_cube(&a.response1, &b.response1);
    }
}

#[test]
fn diff_random_velocity_layers() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1010);
    for iteration in 0..2048 {
        let mut a = Case::new("random_velocity_layers", iteration, &mut rng);
        if iteration % 2 == 0 && a.nlayers > 2 {
            a.t[0] = rng.uniform_f32(2.0, 2.9);
        }
        let mut b = a.clone();
        oxwave_f_franuni::sranuni(iteration as i32 + 1);
        refl::random_velocity_layers(
            &mut a.nlayers,
            &mut a.lsource,
            a.nrand_layers,
            a.sdcl,
            a.sdct,
            a.layer,
            a.zlayer,
            &mut a.cl,
            &mut a.ql,
            &mut a.ct,
            &mut a.qt,
            &mut a.rho,
            &mut a.t,
        );
        unsafe {
            sranuni(iteration as i32 + 1);
        }
        unsafe {
            random_velocity_layers(
                &mut b.nlayers,
                &mut b.lsource,
                b.nrand_layers,
                b.sdcl,
                b.sdct,
                b.layer,
                b.zlayer,
                b.cl.as_mut_ptr(),
                b.ql.as_mut_ptr(),
                b.ct.as_mut_ptr(),
                b.qt.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.t.as_mut_ptr(),
            );
        }
        compare_floats(&a.ct, &b.ct);
        compare_floats(&a.cl, &b.cl);
        compare_floats(&a.rho, &b.rho);
        compare_floats(&a.qt, &b.qt);
        compare_floats(&a.t, &b.t);
        compare_floats(&a.ql, &b.ql);
        assert_eq!(a.nlayers, b.nlayers);
        assert_eq!(a.lsource, b.lsource);
    }
}

#[test]
fn diff_sh_reflectivities() {
    let _quiet = Quiet::new();
    let mut rng = XorShift::new(1011);
    for iteration in 0..2048 {
        let mut a = Case::new("sh_reflectivities", iteration, &mut rng);
        let output = COutput::new();
        a.verbose = (iteration % 4) as i32;
        a.lsource = 2 + (iteration % (a.nlayers as usize - 1)) as i32;
        refl::source_receiver_type(
            a.lsource,
            &a.lobs,
            &a.cl,
            &a.ct,
            &mut a.acoustic,
            &mut a.flag,
        );

        // Include many 100-slowness blocks, their exact boundary, and a remainder.
        if iteration % 8 == 0 {
            a.p2w = 0.202;
        } else if iteration % 8 == 1 {
            a.p2w = 0.4;
        }
        if iteration % 13 == 0 {
            a.m1 = 0.0;
            a.m2 = 0.0;
            a.m3 = 0.0;
            a.h1 = 0.0;
            a.h2 = 0.0;
        }
        refl::source_receiver_type(
            a.lsource,
            &a.lobs,
            &a.cl,
            &a.ct,
            &mut a.acoustic,
            &mut a.flag,
        );
        let mut b = a.clone();
        b.lobs.push(i32::MAX); // readable end sentinel for the C receiver scan
        refl::sh_reflectivities(
            a.int_type,
            a.verbose,
            a.wtype,
            a.layern,
            a.np,
            a.bp1,
            a.m1,
            a.m2,
            a.m3,
            a.h1,
            a.h2,
            a.lsource,
            a.bx,
            a.dx,
            a.xmax,
            a.decay,
            a.fref,
            a.wrefp,
            a.wrefs,
            a.tsec,
            a.p2w,
            a.fs,
            a.epsp,
            a.epss,
            a.sigp,
            a.sigs,
            a.pw1,
            a.pw2,
            a.pw3,
            a.pw4,
            &mut a.flag,
            &a.lobs,
            &a.rho,
            &a.t,
            &a.cl,
            &a.ct,
            &a.ql,
            &a.qt,
            &mut a.response1,
            &oxwave_core::prelude::stderr(),
        );
        let mut response1_rows: Vec<Vec<_>> = b
            .response1
            .iter_mut()
            .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
            .collect();
        let mut response1_planes: Vec<_> =
            response1_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
        unsafe {
            sh_reflectivities(
                b.int_type,
                b.verbose,
                b.wtype,
                b.nw,
                b.nlayers,
                b.nx,
                b.layern,
                b.nor,
                b.np,
                b.bp1,
                b.m1,
                b.m2,
                b.m3,
                b.h1,
                b.h2,
                b.lsource,
                b.bx,
                b.dx,
                b.xmax,
                b.decay,
                b.fref,
                b.wrefp,
                b.wrefs,
                b.tsec,
                b.p2w,
                b.fs,
                b.epsp,
                b.epss,
                b.sigp,
                b.sigs,
                b.pw1,
                b.pw2,
                b.pw3,
                b.pw4,
                b.flag.as_mut_ptr(),
                b.lobs.as_mut_ptr(),
                b.rho.as_mut_ptr(),
                b.t.as_mut_ptr(),
                b.cl.as_mut_ptr(),
                b.ct.as_mut_ptr(),
                b.ql.as_mut_ptr(),
                b.qt.as_mut_ptr(),
                response1_planes.as_mut_ptr(),
                output.0,
            );
        }
        assert_eq!(a.flag, b.flag);
        check_cube(&a.response1, &b.response1);
    }
}
