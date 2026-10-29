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

//! Anisotropic Gaussian beam migration. Port of su/main/migration_inversion/sumiggbzoan.c (Revision 1.5).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: Tariq Alkhalifah, based on MIGGBZO by Dave Hale.
//! CWP: repackaged as an SU program by John Stockwell, April 2006.

#![allow(
    unused_mut,
    unused_assignments,
    unused_variables,
    non_snake_case,
    unused_parens
)]
use super::*;
pub(super) fn cellTimeAmp(cells: &mut Cells, jx: i32, jz: i32) -> i32 {
    let mut lx: i32 = cells.lx;
    let mut lz: i32 = cells.lz;
    let mut dx: f32 = cells.dx;
    let mut fx: f32 = cells.fx;
    let mut dz: f32 = cells.dz;
    let mut fz: f32 = cells.fz;
    let mut wmin: f32 = cells.wmin;
    let mut lmin: f32 = cells.lmin;
    let ray = &mut cells.ray;
    let cell = &mut cells.cell;
    let mut irs: i32 = 0;
    let mut kmah: i32 = 0;
    let mut tmax: f32 = 0.0;
    let mut xc: f32 = 0.0;
    let mut zc: f32 = 0.0;
    let mut t0: f32 = 0.0;
    let mut x0: f32 = 0.0;
    let mut z0: f32 = 0.0;
    let mut x: f32 = 0.0;
    let mut z: f32 = 0.0;
    let mut tr: f32 = 0.0;
    let mut ti: f32 = 0.0;
    let mut ar: f32 = 0.0;
    let mut ai: f32 = 0.0;
    let mut v: f32 = 0.0;
    let mut q1: f32 = 0.0;
    let mut p1: f32 = 0.0;
    let mut q2: f32 = 0.0;
    let mut p2: f32 = 0.0;
    let mut e: f32 = 0.0;
    let mut es: f32 = 0.0;
    let mut scale: f32 = 0.0;
    let mut phase: f32 = 0.0;
    let mut vvmr: f32 = 0.0;
    let mut vvmi: f32 = 0.0;
    let mut c: f32 = 0.0;
    let mut s: f32 = 0.0;
    let mut ov: f32 = 0.0;
    let mut px: f32 = 0.0;
    let mut pz: f32 = 0.0;
    let mut pzpz: f32 = 0.0;
    let mut pxpx: f32 = 0.0;
    let mut pxpz: f32 = 0.0;
    let mut dvdx: f32 = 0.0;
    let mut dvdz: f32 = 0.0;
    let mut dvds: f32 = 0.0;
    let mut dvdn: f32 = 0.0;
    let mut wxxr: f32 = 0.0;
    let mut wxzr: f32 = 0.0;
    let mut wzzr: f32 = 0.0;
    let mut wxxi: f32 = 0.0;
    let mut wxzi: f32 = 0.0;
    let mut wzzi: f32 = 0.0;
    /* maximum time */
    tmax = ray.rs[(ray.rs.len() as i32 - 1i32) as usize].t;
    /* cell coordinates */
    xc = (fx + (((jx * lx) as f32) * dx));
    zc = (fz + (((jz * lz) as f32) * dz));
    /* ray step nearest to cell */
    irs = nearestRayStep(ray, xc, zc);
    let rs = ray.rs[(irs) as usize];
    /* ray time and coordinates */
    t0 = rs.t;
    x0 = rs.x;
    z0 = rs.z;
    /* real and imaginary parts of v*v*m = v*v*p/q */
    v = rs.v;
    q1 = rs.q1;
    p1 = rs.p1;
    q2 = rs.q2;
    p2 = rs.p2;
    e = ((wmin * lmin) * lmin);
    es = (e * e);
    scale = ((v * v) / ((q2 * q2) + ((es * q1) * q1)));
    vvmr = (((p2 * q2) + ((es * p1) * q1)) * scale);
    vvmi = ((-e) * scale);
    /* components of slowness vector, px and pz */
    c = rs.c;
    s = rs.s;
    ov = ((1i32 as f32) / v);
    px = (s * ov);
    pz = (c * ov);
    pzpz = (pz * pz);
    pxpx = (px * px);
    pxpz = (px * pz);
    /* velocity derivatives along tangent and normal */
    dvdx = rs.dvdx;
    dvdz = rs.dvdz;
    dvds = ((s * dvdx) + (c * dvdz));
    dvdn = ((c * dvdx) - (s * dvdz));
    /* real part of W matrix */
    wxxr = (((((pzpz * vvmr) as f64) - ((2_f64 * (pxpz as f64)) * (dvdn as f64)))
        - ((pxpx * dvds) as f64)) as f32);
    wxzr = (((pxpx - pzpz) * dvdn) - (pxpz * (vvmr + dvds)));
    wzzr = (((((pxpx * vvmr) as f64) + ((2_f64 * (pxpz as f64)) * (dvdn as f64)))
        - ((pzpz * dvds) as f64)) as f32);
    /* imaginary part of W matrix */
    wxxi = (pzpz * vvmi);
    wxzi = ((-pxpz) * vvmi);
    wzzi = (pxpx * vvmi);
    /* vector from ray to cell */
    x = (xc - x0);
    z = (zc - z0);
    /* real and imaginary parts of complex time */
    tr = ((((t0 as f64) + (((px as f64) + ((0.5_f64 * (wxxr as f64)) * (x as f64))) * (x as f64)))
        + ((((pz + (wxzr * x)) as f64) + ((0.5_f64 * (wzzr as f64)) * (z as f64))) * (z as f64)))
        as f32);
    ti = (((((0.5_f64 * (wxxi as f64)) * (x as f64)) * (x as f64))
        + ((((wxzi * x) as f64) + ((0.5_f64 * (wzzi as f64)) * (z as f64))) * (z as f64)))
        as f32);
    /* real and imaginary parts of complex amplitude */
    kmah = rs.kmah;
    scale = (((es / ((q2 * q2) + ((es * q1) * q1))) as f64).powf(0.25_f64) as f32);
    phase = ((0.5_f64
        * ((q2 as f64).atan2(((q1 * e) as f64))
            + ((2_f64 * (3.1415926535897931_f64))
                * (((kmah + 1i32).checked_div(2i32).unwrap_or(0)) as f64)))) as f32);
    ar = (((scale as f64) * (phase as f64).cos()) as f32);
    ai = (((scale as f64) * (phase as f64).sin()) as f32);
    /*fprintf(stderr,"tr=%f ti=%f ar=%f ai=%f\n",tr,ti,ar,ai);*/
    /* set cell parameters */
    cell[(jx) as usize][(jz) as usize].tr = tr;
    cell[(jx) as usize][(jz) as usize].ti = ti;
    cell[(jx) as usize][(jz) as usize].ar = ar;
    cell[(jx) as usize][(jz) as usize].ai = ai;
    /* return 1 if Gaussian amplitude is significant, 0 otherwise */
    return (if ((((wmin * ti) as f64) > (-5_f64)) && (tr <= tmax)) {
        1i32
    } else {
        0i32
    });
}
pub(super) fn cellBeam(cells: &mut Cells, jx: i32, jz: i32) {
    let mut lx: i32 = cells.lx;
    let mut lz: i32 = cells.lz;
    let mut nx: i32 = cells.nx;
    let mut nz: i32 = cells.nz;
    let g = &mut cells.g;
    let cell = &cells.cell;
    let bd = &cells.bd;
    let mut ntr: i32 = bd.ntr;
    let mut nti: i32 = bd.nti;
    let mut dtr: f32 = bd.dtr;
    let mut ftr: f32 = bd.ftr;
    let mut dti: f32 = bd.dti;
    let mut fti: f32 = bd.fti;
    let cf = &bd.cf;
    let mut kxlo: i32 = 0;
    let mut kxhi: i32 = 0;
    let mut kzlo: i32 = 0;
    let mut kzhi: i32 = 0;
    let mut kx: i32 = 0;
    let mut kz: i32 = 0;
    let mut itr: i32 = 0;
    let mut iti: i32 = 0;
    let mut odtr: f32 = 0.0;
    let mut odti: f32 = 0.0;
    let mut t00r: f32 = 0.0;
    let mut t01r: f32 = 0.0;
    let mut t10r: f32 = 0.0;
    let mut t11r: f32 = 0.0;
    let mut t00i: f32 = 0.0;
    let mut t01i: f32 = 0.0;
    let mut t10i: f32 = 0.0;
    let mut t11i: f32 = 0.0;
    let mut a00r: f32 = 0.0;
    let mut a01r: f32 = 0.0;
    let mut a10r: f32 = 0.0;
    let mut a11r: f32 = 0.0;
    let mut a00i: f32 = 0.0;
    let mut a01i: f32 = 0.0;
    let mut a10i: f32 = 0.0;
    let mut a11i: f32 = 0.0;
    let mut tx0r: f32 = 0.0;
    let mut tx1r: f32 = 0.0;
    let mut tx0i: f32 = 0.0;
    let mut tx1i: f32 = 0.0;
    let mut ax0r: f32 = 0.0;
    let mut ax1r: f32 = 0.0;
    let mut ax0i: f32 = 0.0;
    let mut ax1i: f32 = 0.0;
    let mut txzr: f32 = 0.0;
    let mut txzi: f32 = 0.0;
    let mut axzr: f32 = 0.0;
    let mut axzi: f32 = 0.0;
    let mut dtx0r: f32 = 0.0;
    let mut dtx0i: f32 = 0.0;
    let mut dtx1r: f32 = 0.0;
    let mut dtx1i: f32 = 0.0;
    let mut dax0r: f32 = 0.0;
    let mut dax0i: f32 = 0.0;
    let mut dax1r: f32 = 0.0;
    let mut dax1i: f32 = 0.0;
    let mut dtxzr: f32 = 0.0;
    let mut dtxzi: f32 = 0.0;
    let mut daxzr: f32 = 0.0;
    let mut daxzi: f32 = 0.0;
    let mut xdelta: f32 = 0.0;
    let mut zdelta: f32 = 0.0;
    let mut trn: f32 = 0.0;
    let mut tin: f32 = 0.0;
    let mut trfrac: f32 = 0.0;
    let mut mtrfrac: f32 = 0.0;
    let mut tifrac: f32 = 0.0;
    let mut mtifrac: f32 = 0.0;
    let mut cf0r: f32 = 0.0;
    let mut cf0i: f32 = 0.0;
    let mut cf1r: f32 = 0.0;
    let mut cf1i: f32 = 0.0;
    let mut cfr: f32 = 0.0;
    let mut cfi: f32 = 0.0;
    /* inverse of time sampling intervals */
    odtr = ((1_f64 / (dtr as f64)) as f32);
    odti = ((1_f64 / (dti as f64)) as f32);
    /* complex time and amplitude for each corner */
    t00r = cell[(jx) as usize][(jz) as usize].tr;
    t01r = cell[(jx) as usize][(jz + 1i32) as usize].tr;
    t10r = cell[(jx + 1i32) as usize][(jz) as usize].tr;
    t11r = cell[(jx + 1i32) as usize][(jz + 1i32) as usize].tr;
    t00i = cell[(jx) as usize][(jz) as usize].ti;
    t01i = cell[(jx) as usize][(jz + 1i32) as usize].ti;
    t10i = cell[(jx + 1i32) as usize][(jz) as usize].ti;
    t11i = cell[(jx + 1i32) as usize][(jz + 1i32) as usize].ti;
    a00r = cell[(jx) as usize][(jz) as usize].ar;
    a01r = cell[(jx) as usize][(jz + 1i32) as usize].ar;
    a10r = cell[(jx + 1i32) as usize][(jz) as usize].ar;
    a11r = cell[(jx + 1i32) as usize][(jz + 1i32) as usize].ar;
    a00i = cell[(jx) as usize][(jz) as usize].ai;
    a01i = cell[(jx) as usize][(jz + 1i32) as usize].ai;
    a10i = cell[(jx + 1i32) as usize][(jz) as usize].ai;
    a11i = cell[(jx + 1i32) as usize][(jz + 1i32) as usize].ai;
    /* x and z samples for cell */
    kxlo = (jx * lx);
    kxhi = (kxlo + lx);
    if (kxhi > nx) {
        kxhi = nx;
    }
    kzlo = (jz * lz);
    kzhi = (kzlo + lz);
    if (kzhi > nz) {
        kzhi = nz;
    }
    /* fractional increments for linear interpolation */
    xdelta = ((1_f64 / (lx as f64)) as f32);
    zdelta = ((1_f64 / (lz as f64)) as f32);
    /* increments for times and amplitudes at top and bottom of cell */
    dtx0r = ((t10r - t00r) * xdelta);
    dtx1r = ((t11r - t01r) * xdelta);
    dtx0i = ((t10i - t00i) * xdelta);
    dtx1i = ((t11i - t01i) * xdelta);
    dax0r = ((a10r - a00r) * xdelta);
    dax1r = ((a11r - a01r) * xdelta);
    dax0i = ((a10i - a00i) * xdelta);
    dax1i = ((a11i - a01i) * xdelta);
    /* times and amplitudes at top-left and bottom-left of cell */
    tx0r = t00r;
    tx1r = t01r;
    tx0i = t00i;
    tx1i = t01i;
    ax0r = a00r;
    ax1r = a01r;
    ax0i = a00i;
    ax1i = a01i;
    /* loop over x samples */
    {
        kx = kxlo;
        while (kx < kxhi) {
            /* increments for time and amplitude */
            dtxzr = ((tx1r - tx0r) * zdelta);
            dtxzi = ((tx1i - tx0i) * zdelta);
            daxzr = ((ax1r - ax0r) * zdelta);
            daxzi = ((ax1i - ax0i) * zdelta);
            /* time and amplitude at top of cell */
            txzr = tx0r;
            txzi = tx0i;
            axzr = ax0r;
            axzi = ax0i;
            /* loop over z samples */
            {
                kz = kzlo;
                while (kz < kzhi) {
                    /* index of imaginary time */
                    iti = ({
                        tin = ((txzi - fti) * odti);
                        tin
                    } as i32);
                    if ((iti < 0i32) || (iti >= (nti - 1i32))) {
                        kz += 1;
                        continue;
                    }
                    /* pointers to left and right imaginary time samples */
                    let cf0 = &cf[(iti) as usize];
                    let cf1 = &cf[(iti + 1i32) as usize];
                    /* imaginary time linear interpolation coefficients */
                    tifrac = (tin - (iti as f32));
                    mtifrac = ((1_f64 - (tifrac as f64)) as f32);
                    /* index of real time */
                    itr = ({
                        trn = ((txzr - ftr) * odtr);
                        trn
                    } as i32);
                    if ((itr < 0i32) || (itr >= (ntr - 1i32))) {
                        kz += 1;
                        continue;
                    }
                    /* real time linear interpolation coefficients */
                    trfrac = (trn - (itr as f32));
                    mtrfrac = ((1_f64 - (trfrac as f64)) as f32);
                    /* real and imaginary parts of complex beam data */
                    cf0r = ((mtrfrac * cf0[(itr) as usize].r)
                        + (trfrac * cf0[(itr + 1i32) as usize].r));
                    cf1r = ((mtrfrac * cf1[(itr) as usize].r)
                        + (trfrac * cf1[(itr + 1i32) as usize].r));
                    cfr = ((mtifrac * cf0r) + (tifrac * cf1r));
                    cf0i = ((mtrfrac * cf0[(itr) as usize].i)
                        + (trfrac * cf0[(itr + 1i32) as usize].i));
                    cf1i = ((mtrfrac * cf1[(itr) as usize].i)
                        + (trfrac * cf1[(itr + 1i32) as usize].i));
                    cfi = ((mtifrac * cf0i) + (tifrac * cf1i));
                    /* accumulate beam */
                    g[(kx) as usize][(kz) as usize] =
                        (g[(kx) as usize][(kz) as usize] + ((axzr * cfr) - (axzi * cfi)));
                    /* increment time and amplitude */
                    txzr = (txzr + dtxzr);
                    txzi = (txzi + dtxzi);
                    axzr = (axzr + daxzr);
                    axzi = (axzi + daxzi);
                    kz += 1;
                }
            }
            /* increment times and amplitudes at top and bottom of cell */
            tx0r = (tx0r + dtx0r);
            tx1r = (tx1r + dtx1r);
            tx0i = (tx0i + dtx0i);
            tx1i = (tx1i + dtx1i);
            ax0r = (ax0r + dax0r);
            ax1r = (ax1r + dax1r);
            ax0i = (ax0i + dax0i);
            ax1i = (ax1i + dax1i);
            kx += 1;
        }
    }
}
