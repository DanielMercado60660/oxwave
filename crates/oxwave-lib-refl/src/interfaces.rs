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

//! Slice interfaces for Refl/lib/auxil.c, compute.c, convolve.c,
//! psvreflect.c, shreflect.c and synthetic.c.
//! Credits: Copyright (c) Colorado School of Mines, 2011.
//! The C sources' Credits sections name no individual authors.
use crate::bounds::BoundScope;
use crate::kernels::*;
use oxwave_core::prelude::CFile;
use oxwave_cwp::complex::Complex;

fn same_layers(n: i32, lengths: &[usize]) {
    assert!(n >= 0 && lengths.iter().all(|&len| len >= n as usize));
}
fn model(
    ls: i32,
    obs: &[i32],
    cl: &[f32],
    ct: &[f32],
    ql: &[f32],
    qt: &[f32],
    rho: &[f32],
    t: &[f32],
) {
    assert!(cl.len() >= 2 && ls >= 2 && ls as usize <= cl.len());
    same_layers(
        cl.len() as i32,
        &[ct.len(), ql.len(), qt.len(), rho.len(), t.len()],
    );
    assert!(obs.iter().all(|&l| l >= 1 && l as usize <= cl.len()));
}
fn matrix<T>(a: &[Vec<T>], n: i32, m: i32) {
    assert!(n >= 0 && m >= 0 && a.len() >= n as usize);
    assert!(a[..n as usize].iter().all(|r| r.len() >= m as usize));
}
fn cube<T>(a: &[Vec<Vec<T>>], n: i32, m: i32, k: i32) {
    assert!(n >= 0 && a.len() >= n as usize);
    for plane in &a[..n as usize] {
        matrix(plane, m, k);
    }
}

/// C prototype: `void apply_filters(int * min_phase, int nfilters, int nw, float tsec, float * f1, float * f2, float * dbpo, int * filter_type, complex * refw);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn apply_filters(
    min_phase: &[i32],
    nw: i32,
    tsec: f32,
    f1: &[f32],
    f2: &[f32],
    dbpo: &[f32],
    filter_type: &[i32],
    refw: &mut [Complex],
) {
    let nfilters = i32::try_from(f1.len()).expect("array too large for C dimensions");
    assert!(nw >= 1 && tsec > 0.0);
    assert!(
        min_phase.len() >= f1.len()
            && f2.len() >= f1.len()
            && dbpo.len() >= f1.len()
            && filter_type.len() >= f1.len()
    );
    assert!(refw.len() >= 2 * nw as usize);
    for i in 0..f1.len() {
        if filter_type[i] == 3 {
            let delta = if min_phase[..=i].contains(&1) {
                (2.0 * std::f64::consts::PI / tsec as f64) as f32
            } else {
                0.0
            };
            let first = ((f1[i] - delta) * tsec + 1.0) as i32;
            let last = ((f2[i] + delta) * tsec + 2.0) as i32;
            assert!(last - first < 16 || (first >= 1 && last <= nw + 1));
        }
    }
    if min_phase.contains(&1) {
        assert_eq!(oxwave_pfa_fft::npfa(2 * nw), 2 * nw);
    }

    let mut bounds = BoundScope::default();
    bounds.track_slice(min_phase);
    bounds.track_slice(f1);
    bounds.track_slice(f2);
    bounds.track_slice(dbpo);
    bounds.track_slice(filter_type);
    bounds.track_slice(refw);
    unsafe {
        raw_apply_filters(
            min_phase.as_ptr() as *mut i32,
            nfilters,
            nw,
            tsec,
            f1.as_ptr() as *mut f32,
            f2.as_ptr() as *mut f32,
            dbpo.as_ptr() as *mut f32,
            filter_type.as_ptr() as *mut i32,
            refw.as_mut_ptr(),
        );
    }
}

/// C prototype: `void compute_al_at(int wtype, int nlayers, int layern, float eps, float epsp, float epss, float sigma, float sigp, float sigs, float * wrefp, float * wrefs, float fw, float * cl, float * ct, float * ql, float * qt, float * rho, complex wpie, complex * al, complex * at);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn compute_al_at(
    wtype: i32,
    layern: i32,
    eps: f32,
    epsp: f32,
    epss: f32,
    sigma: f32,
    sigp: f32,
    sigs: f32,
    wrefp: &mut f32,
    wrefs: &mut f32,
    fw: f32,
    cl: &[f32],
    ct: &[f32],
    ql: &[f32],
    qt: &[f32],
    rho: &[f32],
    wpie: Complex,
    al: &mut [Complex],
    at: &mut [Complex],
) {
    let nlayers = i32::try_from(cl.len()).expect("array too large for C dimensions");
    assert!(layern >= 0 && layern <= nlayers);
    same_layers(
        nlayers,
        &[ct.len(), ql.len(), qt.len(), rho.len(), al.len(), at.len()],
    );

    let mut bounds = BoundScope::default();
    bounds.track_slice(cl);
    bounds.track_slice(ct);
    bounds.track_slice(ql);
    bounds.track_slice(qt);
    bounds.track_slice(rho);
    bounds.track_slice(al);
    bounds.track_slice(at);
    unsafe {
        raw_compute_al_at(
            wtype,
            nlayers,
            layern,
            eps,
            epsp,
            epss,
            sigma,
            sigp,
            sigs,
            wrefp,
            wrefs,
            fw,
            cl.as_ptr() as *mut f32,
            ct.as_ptr() as *mut f32,
            ql.as_ptr() as *mut f32,
            qt.as_ptr() as *mut f32,
            rho.as_ptr() as *mut f32,
            wpie,
            al.as_mut_ptr(),
            at.as_mut_ptr(),
        );
    }
}

/// C prototype: `void compute_p_aux_arrays(int wtype, int nlayers, int lsource, int block_size, float bp, float dp, float w, float decay, float pw1, float pw2, float pw3, float pw4, complex * pwin, float m1, float m2, float m3, float h1, float h2, complex wpie, complex * divfac, complex * p, complex * pp, complex * al, complex * at, float * rho, complex * gl, complex * gt, complex * gam, complex * alpha, complex * betha, complex * sigmad1, complex * sigmad2, complex * sigmau1, complex * sigmau2);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn compute_p_aux_arrays(
    wtype: i32,
    lsource: i32,
    bp: f32,
    dp: f32,
    w: f32,
    decay: f32,
    pw1: f32,
    pw2: f32,
    pw3: f32,
    pw4: f32,
    pwin: &mut [Complex],
    m1: f32,
    m2: f32,
    m3: f32,
    h1: f32,
    h2: f32,
    wpie: Complex,
    divfac: &mut Complex,
    p: &mut [Complex],
    pp: &mut [Complex],
    al: &[Complex],
    at: &[Complex],
    rho: &[f32],
    gl: &mut [Complex],
    gt: &mut [Complex],
    gam: &mut [Complex],
    alpha: &mut [Complex],
    betha: &mut [Complex],
    sigmad1: &mut [Complex],
    sigmad2: &mut [Complex],
    sigmau1: &mut [Complex],
    sigmau2: &mut [Complex],
) {
    let nlayers = i32::try_from(al.len()).expect("array too large for C dimensions");
    let block_size = i32::try_from(pwin.len()).expect("array too large for C dimensions");
    assert!(lsource >= 1 && lsource <= nlayers);
    same_layers(nlayers, &[at.len(), rho.len()]);
    same_layers(
        block_size,
        &[
            p.len(),
            pp.len(),
            alpha.len(),
            betha.len(),
            sigmad1.len(),
            sigmad2.len(),
            sigmau1.len(),
            sigmau2.len(),
        ],
    );
    same_layers(nlayers * block_size, &[gl.len(), gt.len(), gam.len()]);

    let mut bounds = BoundScope::default();
    bounds.track_slice(pwin);
    bounds.track_slice(p);
    bounds.track_slice(pp);
    bounds.track_slice(al);
    bounds.track_slice(at);
    bounds.track_slice(rho);
    bounds.track_slice(gl);
    bounds.track_slice(gt);
    bounds.track_slice(gam);
    bounds.track_slice(alpha);
    bounds.track_slice(betha);
    bounds.track_slice(sigmad1);
    bounds.track_slice(sigmad2);
    bounds.track_slice(sigmau1);
    bounds.track_slice(sigmau2);
    unsafe {
        raw_compute_p_aux_arrays(
            wtype,
            nlayers,
            lsource,
            block_size,
            bp,
            dp,
            w,
            decay,
            pw1,
            pw2,
            pw3,
            pw4,
            pwin.as_mut_ptr(),
            m1,
            m2,
            m3,
            h1,
            h2,
            wpie,
            divfac,
            p.as_mut_ptr(),
            pp.as_mut_ptr(),
            al.as_ptr() as *mut Complex,
            at.as_ptr() as *mut Complex,
            rho.as_ptr() as *mut f32,
            gl.as_mut_ptr(),
            gt.as_mut_ptr(),
            gam.as_mut_ptr(),
            alpha.as_mut_ptr(),
            betha.as_mut_ptr(),
            sigmad1.as_mut_ptr(),
            sigmad2.as_mut_ptr(),
            sigmau1.as_mut_ptr(),
            sigmau2.as_mut_ptr(),
        );
    }
}

/// C prototype: `void compute_reflectivities(int int_type, int verbose, int wtype, int wfield, int vsp, int flt, int win, int nx, int nt, int ntc, int nor, int nf, int nlayers, int lsource, int layern, int nfilters, int * filters_phase, int nw, int np, float bp, float tlag, float red_vel, float w1, float w2, float fx, float dx, float bx, float fs, float decay, float p2w, float tsec, float fw, float wrefp, float wrefs, float epsp, float epss, float sigp, float sigs, float pw1, float pw2, float pw3, float pw4, float h1, float h2, float m1, float m2, float m3, float fref, int * lobs, int * filters_type, float * dbpo, float * f1, float * f2, float * cl, float * ct, float * ql, float * qt, float * rho, float * t, float ** wavefield1, float ** wavefield2, float ** wavefield3, FILE * outfp);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn compute_reflectivities(
    int_type: i32,
    verbose: i32,
    wtype: i32,
    wfield: i32,
    vsp: i32,
    flt: i32,
    win: i32,
    ntc: i32,
    nf: i32,
    lsource: i32,
    layern: i32,
    filters_phase: &[i32],
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
    lobs: &[i32],
    filters_type: &[i32],
    dbpo: &[f32],
    f1: &[f32],
    f2: &[f32],
    cl: &[f32],
    ct: &[f32],
    ql: &[f32],
    qt: &[f32],
    rho: &[f32],
    t: &[f32],
    wavefield1: &mut [Vec<f32>],
    wavefield2: &mut [Vec<f32>],
    wavefield3: &mut [Vec<f32>],
    outfp: &CFile,
) {
    let nx = i32::try_from(wavefield1.len()).expect("array too large for C dimensions");
    let nt = i32::try_from(wavefield1.first().map_or(0, Vec::len))
        .expect("array too large for C dimensions");
    let nor = i32::try_from(lobs.len()).expect("array too large for C dimensions");
    let nlayers = i32::try_from(cl.len()).expect("array too large for C dimensions");
    let nfilters = i32::try_from(f1.len()).expect("array too large for C dimensions");
    model(lsource, lobs, cl, ct, ql, qt, rho, t);
    assert!(nor <= 3); // C uses fixed 300-element receiver/block scratch arrays.

    matrix(wavefield1, nx, nt);
    matrix(wavefield2, nx, nt);
    matrix(wavefield3, nx, nt);
    same_layers(
        nfilters,
        &[
            filters_phase.len(),
            filters_type.len(),
            f2.len(),
            dbpo.len(),
        ],
    );
    assert!(ntc > 0 && nt <= oxwave_pfa_fft::npfa(ntc));
    // Original SH branch passes nx and ntc in the opposite order.
    assert!(wtype != 2 || nx == ntc);

    let mut wavefield1_rows: Vec<_> = wavefield1.iter_mut().map(|r| r.as_mut_ptr()).collect();
    let mut wavefield2_rows: Vec<_> = wavefield2.iter_mut().map(|r| r.as_mut_ptr()).collect();
    let mut wavefield3_rows: Vec<_> = wavefield3.iter_mut().map(|r| r.as_mut_ptr()).collect();
    // C scans one receiver past the last match while finishing the layers.
    let mut lobs_with_end = lobs.to_vec();
    lobs_with_end.push(i32::MAX);
    let mut bounds = BoundScope::default();
    bounds.track_slice(&lobs_with_end);
    bounds.track_slice(filters_phase);
    bounds.track_slice(lobs);
    bounds.track_slice(filters_type);
    bounds.track_slice(dbpo);
    bounds.track_slice(f1);
    bounds.track_slice(f2);
    bounds.track_slice(cl);
    bounds.track_slice(ct);
    bounds.track_slice(ql);
    bounds.track_slice(qt);
    bounds.track_slice(rho);
    bounds.track_slice(t);
    bounds.track_slice(&wavefield1_rows);
    for row in wavefield1.iter() {
        bounds.track_slice(row);
    }
    bounds.track_slice(&wavefield2_rows);
    for row in wavefield2.iter() {
        bounds.track_slice(row);
    }
    bounds.track_slice(&wavefield3_rows);
    for row in wavefield3.iter() {
        bounds.track_slice(row);
    }
    unsafe {
        raw_compute_reflectivities(
            int_type,
            verbose,
            wtype,
            wfield,
            vsp,
            flt,
            win,
            nx,
            nt,
            ntc,
            nor,
            nf,
            nlayers,
            lsource,
            layern,
            nfilters,
            filters_phase.as_ptr() as *mut i32,
            nw,
            np,
            bp,
            tlag,
            red_vel,
            w1,
            w2,
            fx,
            dx,
            bx,
            fs,
            decay,
            p2w,
            tsec,
            fw,
            wrefp,
            wrefs,
            epsp,
            epss,
            sigp,
            sigs,
            pw1,
            pw2,
            pw3,
            pw4,
            h1,
            h2,
            m1,
            m2,
            m3,
            fref,
            lobs_with_end.as_mut_ptr(),
            filters_type.as_ptr() as *mut i32,
            dbpo.as_ptr() as *mut f32,
            f1.as_ptr() as *mut f32,
            f2.as_ptr() as *mut f32,
            cl.as_ptr() as *mut f32,
            ct.as_ptr() as *mut f32,
            ql.as_ptr() as *mut f32,
            qt.as_ptr() as *mut f32,
            rho.as_ptr() as *mut f32,
            t.as_ptr() as *mut f32,
            wavefield1_rows.as_mut_ptr(),
            wavefield2_rows.as_mut_ptr(),
            wavefield3_rows.as_mut_ptr(),
            outfp as *const CFile as *mut CFile,
        );
    }
}

/// C prototype: `void compute_sigmas(int wtype, int block_size, int lsource, float m1, float m2, complex meu, complex * alpha, complex * betha, complex * p, complex * gl, complex * gt, complex * gam, complex a2, complex a3, complex * pwin, complex s1, complex s2, complex s4, complex * sigmau1, complex * sigmau2, complex * sigmad1, complex * sigmad2);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn compute_sigmas(
    wtype: i32,
    lsource: i32,
    m1: f32,
    m2: f32,
    meu: Complex,
    alpha: &[Complex],
    betha: &[Complex],
    p: &mut [Complex],
    gl: &[Complex],
    gt: &mut [Complex],
    gam: &[Complex],
    a2: Complex,
    a3: Complex,
    pwin: &[Complex],
    s1: Complex,
    s2: Complex,
    s4: Complex,
    sigmau1: &mut [Complex],
    sigmau2: &mut [Complex],
    sigmad1: &mut [Complex],
    sigmad2: &mut [Complex],
) {
    let block_size = i32::try_from(alpha.len()).expect("array too large for C dimensions");
    assert!(lsource >= 1);
    same_layers(
        block_size,
        &[
            betha.len(),
            p.len(),
            pwin.len(),
            sigmau1.len(),
            sigmau2.len(),
            sigmad1.len(),
            sigmad2.len(),
        ],
    );
    same_layers(lsource * block_size, &[gl.len(), gt.len(), gam.len()]);

    let mut bounds = BoundScope::default();
    bounds.track_slice(alpha);
    bounds.track_slice(betha);
    bounds.track_slice(p);
    bounds.track_slice(gl);
    bounds.track_slice(gt);
    bounds.track_slice(gam);
    bounds.track_slice(pwin);
    bounds.track_slice(sigmau1);
    bounds.track_slice(sigmau2);
    bounds.track_slice(sigmad1);
    bounds.track_slice(sigmad2);
    unsafe {
        raw_compute_sigmas(
            wtype,
            block_size,
            lsource,
            m1,
            m2,
            meu,
            alpha.as_ptr() as *mut Complex,
            betha.as_ptr() as *mut Complex,
            p.as_mut_ptr(),
            gl.as_ptr() as *mut Complex,
            gt.as_mut_ptr(),
            gam.as_ptr() as *mut Complex,
            a2,
            a3,
            pwin.as_ptr() as *mut Complex,
            s1,
            s2,
            s4,
            sigmau1.as_mut_ptr(),
            sigmau2.as_mut_ptr(),
            sigmad1.as_mut_ptr(),
            sigmad2.as_mut_ptr(),
        );
    }
}

/// C prototype: `void compute_synthetics(int verbose, int nt, int ntc, int nx, int nor, int nw, int nlayers, int lsource, int nf, int flt, int vsp, int win, int wtype, float tlag, float red_vel, float decay, float tsec, float bx, float dx, float w1, float w2, int * fil_phase, int nfilters, int * fil_type, float * dbpo, float * f1, float * f2, int * lobs, float * cl, float * t, complex *** response, float ** reflectivity, FILE * outfp);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn compute_synthetics(
    verbose: i32,
    ntc: i32,
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
    fil_phase: &[i32],
    fil_type: &[i32],
    dbpo: &[f32],
    f1: &[f32],
    f2: &[f32],
    lobs: &[i32],
    cl: &[f32],
    t: &[f32],
    response: &mut [Vec<Vec<Complex>>],
    reflectivity: &mut [Vec<f32>],
    outfp: &CFile,
) {
    let nt = i32::try_from(reflectivity.first().map_or(0, Vec::len))
        .expect("array too large for C dimensions");
    let nx = i32::try_from(reflectivity.len()).expect("array too large for C dimensions");
    let nor = i32::try_from(lobs.len()).expect("array too large for C dimensions");
    let nw = i32::try_from(response.len()).expect("array too large for C dimensions");
    let nlayers = i32::try_from(cl.len()).expect("array too large for C dimensions");
    let nfilters = i32::try_from(f1.len()).expect("array too large for C dimensions");
    assert!(ntc > 0 && nt <= oxwave_pfa_fft::npfa(ntc));
    assert!(lsource >= 1 && lsource <= nlayers && t.len() >= cl.len());
    assert!(lobs.iter().all(|&l| l >= 1 && l <= nlayers));
    same_layers(
        nfilters,
        &[fil_phase.len(), fil_type.len(), f2.len(), dbpo.len()],
    );
    cube(response, nw, nx, nor);
    matrix(reflectivity, nx, nt);
    let mut response_rows: Vec<Vec<_>> = response
        .iter_mut()
        .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
        .collect();
    let mut response_planes: Vec<_> = response_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
    let mut reflectivity_rows: Vec<_> = reflectivity.iter_mut().map(|r| r.as_mut_ptr()).collect();
    let mut bounds = BoundScope::default();
    bounds.track_slice(fil_phase);
    bounds.track_slice(fil_type);
    bounds.track_slice(dbpo);
    bounds.track_slice(f1);
    bounds.track_slice(f2);
    bounds.track_slice(lobs);
    bounds.track_slice(cl);
    bounds.track_slice(t);
    bounds.track_slice(&response_planes);
    for rows in &response_rows {
        bounds.track_slice(rows);
    }
    for plane in response.iter() {
        for row in plane {
            bounds.track_slice(row);
        }
    }
    bounds.track_slice(&reflectivity_rows);
    for row in reflectivity.iter() {
        bounds.track_slice(row);
    }
    unsafe {
        raw_compute_synthetics(
            verbose,
            nt,
            ntc,
            nx,
            nor,
            nw,
            nlayers,
            lsource,
            nf,
            flt,
            vsp,
            win,
            wtype,
            tlag,
            red_vel,
            decay,
            tsec,
            bx,
            dx,
            w1,
            w2,
            fil_phase.as_ptr() as *mut i32,
            nfilters,
            fil_type.as_ptr() as *mut i32,
            dbpo.as_ptr() as *mut f32,
            f1.as_ptr() as *mut f32,
            f2.as_ptr() as *mut f32,
            lobs.as_ptr() as *mut i32,
            cl.as_ptr() as *mut f32,
            t.as_ptr() as *mut f32,
            response_planes.as_mut_ptr(),
            reflectivity_rows.as_mut_ptr(),
            outfp as *const CFile as *mut CFile,
        );
    }
}

/// C prototype: `void compute_w_aux_arrays(int wtype, int layern, int nlayers, int nor, int lsource, int * np, int * block_size, int * nblock, int * left, float fref, float wrefp, float wrefs, int * lobs, float tsec, float p2w, float fs, float xmax, float w, float decay, float epsp, float epss, float sigp, float sigs, float * pw1, float * pw2, float * pw3, float * pw4, float * dp, float * fp, complex wpie, complex * cdp, float * cl, float * ct, float * ql, float * qt, float * rho, float * t, complex * al, complex * at, complex * prs);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn compute_w_aux_arrays(
    wtype: i32,
    layern: i32,
    lsource: i32,
    np: &mut i32,
    block_size: &mut i32,
    nblock: &mut i32,
    left: &mut i32,
    fref: f32,
    wrefp: f32,
    wrefs: f32,
    lobs: &[i32],
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
    pw1: &mut f32,
    pw2: &mut f32,
    pw3: &mut f32,
    pw4: &mut f32,
    dp: &mut f32,
    fp: &mut f32,
    wpie: Complex,
    cdp: &mut Complex,
    cl: &[f32],
    ct: &[f32],
    ql: &[f32],
    qt: &[f32],
    rho: &[f32],
    t: &[f32],
    al: &mut [Complex],
    at: &mut [Complex],
    prs: &mut [Complex],
) {
    let nlayers = i32::try_from(cl.len()).expect("array too large for C dimensions");
    let nor = i32::try_from(lobs.len()).expect("array too large for C dimensions");
    model(lsource, lobs, cl, ct, ql, qt, rho, t);
    same_layers(nlayers, &[al.len(), at.len()]);
    same_layers(nor, &[prs.len()]);

    let mut bounds = BoundScope::default();
    bounds.track_slice(lobs);
    bounds.track_slice(cl);
    bounds.track_slice(ct);
    bounds.track_slice(ql);
    bounds.track_slice(qt);
    bounds.track_slice(rho);
    bounds.track_slice(t);
    bounds.track_slice(al);
    bounds.track_slice(at);
    bounds.track_slice(prs);
    unsafe {
        raw_compute_w_aux_arrays(
            wtype,
            layern,
            nlayers,
            nor,
            lsource,
            np,
            block_size,
            nblock,
            left,
            fref,
            wrefp,
            wrefs,
            lobs.as_ptr() as *mut i32,
            tsec,
            p2w,
            fs,
            xmax,
            w,
            decay,
            epsp,
            epss,
            sigp,
            sigs,
            pw1,
            pw2,
            pw3,
            pw4,
            dp,
            fp,
            wpie,
            cdp,
            cl.as_ptr() as *mut f32,
            ct.as_ptr() as *mut f32,
            ql.as_ptr() as *mut f32,
            qt.as_ptr() as *mut f32,
            rho.as_ptr() as *mut f32,
            t.as_ptr() as *mut f32,
            al.as_mut_ptr(),
            at.as_mut_ptr(),
            prs.as_mut_ptr(),
        );
    }
}

/// C prototype: `void construct_tx_trace(int nw, int nt, int ntc, int ntfft, int nfilters, int * filters_phase, float tsec, float unexp, float sphrd, complex * refw, int * fil_type, float * f1, float * f2, float * dbpo, float * reft);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn construct_tx_trace(
    nw: i32,
    ntc: i32,
    filters_phase: &[i32],
    tsec: f32,
    unexp: f32,
    sphrd: f32,
    refw: &mut [Complex],
    fil_type: &[i32],
    f1: &[f32],
    f2: &[f32],
    dbpo: &[f32],
    reft: &mut [f32],
) {
    let nt = i32::try_from(reft.len()).expect("array too large for C dimensions");
    let ntfft = i32::try_from(refw.len()).expect("array too large for C dimensions");
    let nfilters = i32::try_from(f1.len()).expect("array too large for C dimensions");
    assert!(nw >= 0 && ntc > 0 && nt <= ntfft);
    assert_eq!(oxwave_pfa_fft::npfa(ntfft), ntfft);
    same_layers(
        nfilters,
        &[filters_phase.len(), fil_type.len(), f2.len(), dbpo.len()],
    );
    assert!(nfilters == 0 || 2 * nw <= ntfft);

    let mut bounds = BoundScope::default();
    bounds.track_slice(filters_phase);
    bounds.track_slice(refw);
    bounds.track_slice(fil_type);
    bounds.track_slice(f1);
    bounds.track_slice(f2);
    bounds.track_slice(dbpo);
    bounds.track_slice(reft);
    unsafe {
        raw_construct_tx_trace(
            nw,
            nt,
            ntc,
            ntfft,
            nfilters,
            filters_phase.as_ptr() as *mut i32,
            tsec,
            unexp,
            sphrd,
            refw.as_mut_ptr(),
            fil_type.as_ptr() as *mut i32,
            f1.as_ptr() as *mut f32,
            f2.as_ptr() as *mut f32,
            dbpo.as_ptr() as *mut f32,
            reft.as_mut_ptr(),
        );
    }
}

/// C prototype: `void convolve_wavelet(int wavelet_type, int nx, int nt, float dt, float fpeak, float ** wfield);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn convolve_wavelet(wavelet_type: i32, dt: f32, fpeak: f32, wfield: &mut [Vec<f32>]) {
    let nx = i32::try_from(wfield.len()).expect("array too large for C dimensions");
    let nt = i32::try_from(wfield.first().map_or(0, Vec::len))
        .expect("array too large for C dimensions");
    assert!(nt >= 1 && (nt >= 2 || wavelet_type == 2 || wavelet_type == 4));
    matrix(wfield, nx, nt);
    let mut wfield_rows: Vec<_> = wfield.iter_mut().map(|r| r.as_mut_ptr()).collect();
    let mut bounds = BoundScope::default();
    bounds.track_slice(&wfield_rows);
    for row in wfield.iter() {
        bounds.track_slice(row);
    }
    unsafe {
        raw_convolve_wavelet(wavelet_type, nx, nt, dt, fpeak, wfield_rows.as_mut_ptr());
    }
}

/// C prototype: `void psv_reflectivities(int int_type, int verbose, int wtype, int nw, int nlayers, int nx, int layern, int nor, int np, float bp1, float m1, float m2, float m3, float h1, float h2, int lsource, float bx, float dx, float xmax, float decay, float fref, float wrefp, float wrefs, float tsec, float p2w, float fs, float epsp, float epss, float sigp, float sigs, float pw1, float pw2, float pw3, float pw4, int * acoustic, int * flag, int * lobs, float * rho, float * t, float * cl, float * ct, float * ql, float * qt, complex *** response1, complex *** response2, complex *** response3, FILE * outfp);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn psv_reflectivities(
    int_type: i32,
    verbose: i32,
    wtype: i32,
    layern: i32,
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
    acoustic: &[i32],
    flag: &mut [i32],
    lobs: &[i32],
    rho: &[f32],
    t: &[f32],
    cl: &[f32],
    ct: &[f32],
    ql: &[f32],
    qt: &[f32],
    response1: &mut [Vec<Vec<Complex>>],
    response2: &mut [Vec<Vec<Complex>>],
    response3: &mut [Vec<Vec<Complex>>],
    outfp: &CFile,
) {
    let nw = i32::try_from(response1.len()).expect("array too large for C dimensions");
    let nlayers = i32::try_from(cl.len()).expect("array too large for C dimensions");
    let nx = i32::try_from(response1.first().map_or(0, Vec::len))
        .expect("array too large for C dimensions");
    let nor = i32::try_from(lobs.len()).expect("array too large for C dimensions");
    model(lsource, lobs, cl, ct, ql, qt, rho, t);
    assert!(nor <= 3); // C uses fixed 300-element receiver/block scratch arrays.
    same_layers(nor, &[flag.len()]);
    same_layers(nor, &[acoustic.len()]);
    cube(response1, nw, nx, nor);
    cube(response2, nw, nx, nor);
    cube(response3, nw, nx, nor);
    let mut response1_rows: Vec<Vec<_>> = response1
        .iter_mut()
        .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
        .collect();
    let mut response1_planes: Vec<_> = response1_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
    let mut response2_rows: Vec<Vec<_>> = response2
        .iter_mut()
        .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
        .collect();
    let mut response2_planes: Vec<_> = response2_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
    let mut response3_rows: Vec<Vec<_>> = response3
        .iter_mut()
        .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
        .collect();
    let mut response3_planes: Vec<_> = response3_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
    // C scans one receiver past the last match while finishing the layers.
    let mut lobs_with_end = lobs.to_vec();
    lobs_with_end.push(i32::MAX);
    let mut bounds = BoundScope::default();
    bounds.track_slice(&lobs_with_end);
    bounds.track_slice(acoustic);
    bounds.track_slice(flag);
    bounds.track_slice(lobs);
    bounds.track_slice(rho);
    bounds.track_slice(t);
    bounds.track_slice(cl);
    bounds.track_slice(ct);
    bounds.track_slice(ql);
    bounds.track_slice(qt);
    bounds.track_slice(&response1_planes);
    for rows in &response1_rows {
        bounds.track_slice(rows);
    }
    for plane in response1.iter() {
        for row in plane {
            bounds.track_slice(row);
        }
    }
    bounds.track_slice(&response2_planes);
    for rows in &response2_rows {
        bounds.track_slice(rows);
    }
    for plane in response2.iter() {
        for row in plane {
            bounds.track_slice(row);
        }
    }
    bounds.track_slice(&response3_planes);
    for rows in &response3_rows {
        bounds.track_slice(rows);
    }
    for plane in response3.iter() {
        for row in plane {
            bounds.track_slice(row);
        }
    }
    unsafe {
        raw_psv_reflectivities(
            int_type,
            verbose,
            wtype,
            nw,
            nlayers,
            nx,
            layern,
            nor,
            np,
            bp1,
            m1,
            m2,
            m3,
            h1,
            h2,
            lsource,
            bx,
            dx,
            xmax,
            decay,
            fref,
            wrefp,
            wrefs,
            tsec,
            p2w,
            fs,
            epsp,
            epss,
            sigp,
            sigs,
            pw1,
            pw2,
            pw3,
            pw4,
            acoustic.as_ptr() as *mut i32,
            flag.as_mut_ptr(),
            lobs_with_end.as_mut_ptr(),
            rho.as_ptr() as *mut f32,
            t.as_ptr() as *mut f32,
            cl.as_ptr() as *mut f32,
            ct.as_ptr() as *mut f32,
            ql.as_ptr() as *mut f32,
            qt.as_ptr() as *mut f32,
            response1_planes.as_mut_ptr(),
            response2_planes.as_mut_ptr(),
            response3_planes.as_mut_ptr(),
            outfp as *const CFile as *mut CFile,
        );
    }
}

/// C prototype: `void random_velocity_layers(int * nlayers, int * lsource, int nrand_layers, float sdcl, float sdct, float layer, float zlayer, float * cl, float * ql, float * ct, float * qt, float * rho, float * t);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn random_velocity_layers(
    nlayers: &mut i32,
    lsource: &mut i32,
    nrand_layers: i32,
    sdcl: f32,
    sdct: f32,
    layer: f32,
    zlayer: f32,
    cl: &mut [f32],
    ql: &mut [f32],
    ct: &mut [f32],
    qt: &mut [f32],
    rho: &mut [f32],
    t: &mut [f32],
) {
    assert!(*nlayers >= 1 && *lsource >= 2 && *lsource <= *nlayers);
    assert!(layer >= 1.0 && layer <= *nlayers as f32 && layer.fract() == 0.0 && zlayer > 0.0);
    let original = *nlayers as usize;
    same_layers(
        *nlayers,
        &[cl.len(), ct.len(), ql.len(), qt.len(), rho.len(), t.len()],
    );
    let limit = (nrand_layers as f32 + *nlayers as f32 - layer + 1.0) as i32;
    assert!(limit >= 0);
    let mut generated = 0i32;
    for il in layer as usize - 1..original - 1 {
        let count = if t[il] <= zlayer {
            1
        } else {
            (t[il] / zlayer) as i32
        };
        assert!(count >= 1);
        generated += count;
        if generated > limit {
            generated = limit;
            break;
        }
    }
    // The C scratch allocations hold only the original number of layers.
    assert!(generated <= *nlayers && generated > 0);
    let final_layers = generated + layer as i32 - 1;
    same_layers(
        final_layers,
        &[cl.len(), ct.len(), ql.len(), qt.len(), rho.len(), t.len()],
    );

    let mut bounds = BoundScope::default();
    bounds.track_slice(cl);
    bounds.track_slice(ql);
    bounds.track_slice(ct);
    bounds.track_slice(qt);
    bounds.track_slice(rho);
    bounds.track_slice(t);
    unsafe {
        raw_random_velocity_layers(
            nlayers,
            lsource,
            nrand_layers,
            sdcl,
            sdct,
            layer,
            zlayer,
            cl.as_mut_ptr(),
            ql.as_mut_ptr(),
            ct.as_mut_ptr(),
            qt.as_mut_ptr(),
            rho.as_mut_ptr(),
            t.as_mut_ptr(),
        );
    }
}

/// C prototype: `void sh_reflectivities(int int_type, int verbose, int wtype, int nw, int nlayers, int nx, int layern, int nor, int np, float bp1, float m1, float m2, float m3, float h1, float h2, int lsource, float bx, float dx, float xmax, float decay, float fref, float wrefp, float wrefs, float tsec, float p2w, float fs, float epsp, float epss, float sigp, float sigs, float pw1, float pw2, float pw3, float pw4, int * flag, int * lobs, float * rho, float * t, float * cl, float * ct, float * ql, float * qt, complex *** response1, FILE * outfp);`
///
/// Array dimensions are inferred from slices; explicit frequency and transform
/// lengths retain their C meaning. Model and output slices must cover the C dimensions.
pub fn sh_reflectivities(
    int_type: i32,
    verbose: i32,
    wtype: i32,
    layern: i32,
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
    flag: &mut [i32],
    lobs: &[i32],
    rho: &[f32],
    t: &[f32],
    cl: &[f32],
    ct: &[f32],
    ql: &[f32],
    qt: &[f32],
    response1: &mut [Vec<Vec<Complex>>],
    outfp: &CFile,
) {
    let nw = i32::try_from(response1.len()).expect("array too large for C dimensions");
    let nlayers = i32::try_from(cl.len()).expect("array too large for C dimensions");
    let nx = i32::try_from(response1.first().map_or(0, Vec::len))
        .expect("array too large for C dimensions");
    let nor = i32::try_from(lobs.len()).expect("array too large for C dimensions");
    model(lsource, lobs, cl, ct, ql, qt, rho, t);
    assert!(nor <= 3); // C uses fixed 300-element receiver/block scratch arrays.
    same_layers(nor, &[flag.len()]);
    cube(response1, nw, nx, nor);
    let mut response1_rows: Vec<Vec<_>> = response1
        .iter_mut()
        .map(|p| p.iter_mut().map(|r| r.as_mut_ptr()).collect())
        .collect();
    let mut response1_planes: Vec<_> = response1_rows.iter_mut().map(|p| p.as_mut_ptr()).collect();
    // C scans one receiver past the last match while finishing the layers.
    let mut lobs_with_end = lobs.to_vec();
    lobs_with_end.push(i32::MAX);
    let mut bounds = BoundScope::default();
    bounds.track_slice(&lobs_with_end);
    bounds.track_slice(flag);
    bounds.track_slice(lobs);
    bounds.track_slice(rho);
    bounds.track_slice(t);
    bounds.track_slice(cl);
    bounds.track_slice(ct);
    bounds.track_slice(ql);
    bounds.track_slice(qt);
    bounds.track_slice(&response1_planes);
    for rows in &response1_rows {
        bounds.track_slice(rows);
    }
    for plane in response1.iter() {
        for row in plane {
            bounds.track_slice(row);
        }
    }
    unsafe {
        raw_sh_reflectivities(
            int_type,
            verbose,
            wtype,
            nw,
            nlayers,
            nx,
            layern,
            nor,
            np,
            bp1,
            m1,
            m2,
            m3,
            h1,
            h2,
            lsource,
            bx,
            dx,
            xmax,
            decay,
            fref,
            wrefp,
            wrefs,
            tsec,
            p2w,
            fs,
            epsp,
            epss,
            sigp,
            sigs,
            pw1,
            pw2,
            pw3,
            pw4,
            flag.as_mut_ptr(),
            lobs_with_end.as_mut_ptr(),
            rho.as_ptr() as *mut f32,
            t.as_ptr() as *mut f32,
            cl.as_ptr() as *mut f32,
            ct.as_ptr() as *mut f32,
            ql.as_ptr() as *mut f32,
            qt.as_ptr() as *mut f32,
            response1_planes.as_mut_ptr(),
            outfp as *const CFile as *mut CFile,
        );
    }
}
