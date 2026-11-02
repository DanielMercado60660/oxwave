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
// Copyright (c) Colorado School of Mines, 2011. All rights reserved.

//! SURELAN: port of su/main/velocity_analysis/surelan.c (revision 1.8).
//! Reference: Zhenyue Liu, Migration Velocity Analysis, Ph.D. thesis, CWP report #168, 1995.
//! Credits: CWP: Zhenyue Liu, 1995.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SURELAN - compute residual-moveout semblance for cdp gathers based\t",
    "\ton z(h)*z(h) = z(0)*z(0) + r*h*h where z depth and h offset.\t",
    "\t\t\t\t\t\t\t\t\t",
    " surelan <stdin >stdout   [optional parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nr=51\t\t\tnumber of r-parameter samples   \t\t",
    " dr=0.01               r-parameter sampling interval\t\t\t",
    " fr=-0.25               first value of b-parameter\t\t\t",
    " smute=1.5             samples with RMO stretch exceeding smute are zeroed",
    " dzratio=5             ratio of output to input depth sampling intervals",
    " nsmooth=dzratio*2+1   length of semblance num and den smoothing window",
    " verbose=0             =1 for diagnostic print on stderr\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: \t\t\t\t\t\t\t\t",
    " 1. This program is part of Zhenyue Liu's velocity analysis technique.\t",
    " 2. Input migrated traces should be sorted by cdp - surelan outputs a \t",
    "    group of semblance traces every time cdp changes.  Therefore, the  ",
    "    output will be useful only if cdp gathers are input.  \t\t",
    " 3. The parameter r may take negative values. The range of r can be \t",
    "     controlled by maximum of (z(h)*z(h)-z(0)*z(0))/(h*h)   \t\t",
];

fn alloc_accumulator(nz: usize, nr: i32) -> Vec<Vec<f32>> {
    if nr < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let mut rows = Vec::new();
    if rows.try_reserve_exact(nr as usize).is_err() {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    for _ in 0..nr {
        let mut row = Vec::new();
        if row.try_reserve_exact(nz).is_err() {
            err!("ealloc.c: malloc failed (Cannot allocate memory)");
        }
        row.resize(nz, 0.0);
        rows.push(row);
    }
    rows
}

fn main() {
    let argv = args(); let mut tr = Trace::new();
    initargs(&argv); requestdoc(0, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nz = tr.ns as i32; let dz = tr.d1; let fz = tr.f1;
    let mut cdp = tr.cdp; let mut offset = tr.offset as f32;
    let mut nr = 51; let mut dr = 0.01f32; let mut fr = -0.25f32;
    let mut smute = 1.5f32; let mut dzratio = 5; let mut nsmooth = 0; let mut verbose = 0;
    getparint("nr", &mut nr); getparfloat("dr", &mut dr); getparfloat("fr", &mut fr);
    getparfloat("smute", &mut smute);
    if smute <= 1.0 { err!("smute must be greater than 1.0"); }
    getparint("dzratio", &mut dzratio);
    if getparint("nsmooth", &mut nsmooth) == 0 { nsmooth = dzratio.wrapping_mul(2).wrapping_add(1); }
    getparint("verbose", &mut verbose); checkpars();
    // ARM64 C integer division by zero yields zero.
    let nzout = 1 + (nz - 1).checked_div(dzratio).unwrap_or(0);
    CHECK_NT("nzout", nzout); let dzout = dz * dzratio as f32;
    if verbose != 0 {
        warn!("number of depth samples is %d\n", nzout);
        warn!("output depth sampling interval is %g\n", dzout);
        warn!("output depth of first sample is %g\n", fz);
    }
    let mut data = vec![0.0f32; nz as usize];
    let mut num = alloc_accumulator(nz as usize, nr);
    let mut den = alloc_accumulator(nz as usize, nr); let mut nnz = alloc_accumulator(nz as usize, nr);
    if nzout < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut sem = vec![0.0f32; nzout as usize]; let mut gottrace = true; let mut cdpprev = tr.cdp;
    loop {
        if gottrace { offset = tr.offset as f32; cdp = tr.cdp; data.copy_from_slice(&tr.data[..nz as usize]); }
        if cdp != cdpprev || !gottrace {
            tr.offset = 0; tr.cdp = cdpprev; tr.ns = nzout as u16; tr.d1 = dzout;
            for ir in 0..nr as usize {
                for izout in 0..nzout {
                    let iz = izout * dzratio; let mut ismin = iz - nsmooth / 2; let mut ismax = iz + nsmooth / 2;
                    if ismin < 0 { ismin = 0; } if ismax > nz - 1 { ismax = nz - 1; }
                    let mut nsum = 0.0f32; let mut dsum = 0.0f32;
                    for is in ismin..=ismax { let i = is as usize; nsum += num[ir][i] * num[ir][i]; dsum += nnz[ir][i] * den[ir][i]; }
                    sem[izout as usize] = if dsum != 0.0 { nsum / dsum } else { 0.0 };
                }
                tr.data[..nzout as usize].copy_from_slice(&sem); puttr(&mut tr);
                num[ir].fill(0.0); den[ir].fill(0.0); nnz[ir].fill(0.0);
            }
            if verbose != 0 { warn!("tsemblance output for cdp=%d\n", cdpprev); }
            if !gottrace { break; } cdpprev = cdp;
        }
        let mut r = fr;
        for ir in 0..nr as usize {
            let roffs2 = r * offset * offset;
            let znmute = ((roffs2 as f64) / ((smute * smute) as f64 - 1.0)) as f32;
            let znmute = if znmute > 0.0 { (znmute as f64).sqrt() } else { (-(znmute as f64)).sqrt() } as f32;
            let mut iz = ((znmute - fz) / dz) as i32; if iz < 0 { iz = 0; }
            let mut zn = fz + iz as f32 * dz;
            while iz < nz {
                let mut temp = zn * zn + roffs2;
                let zi = if temp > fz * fz { (((temp as f64).sqrt() - fz as f64) / dz as f64) as f32 } else { 0.0 };
                let izi = zi as i32;
                if izi < nz - 1 {
                    if izi < 0 { err!("invalid interpolation index"); }
                    let frac = zi - izi as f32;
                    temp = ((1.0 - frac as f64) * data[izi as usize] as f64 + (frac * data[izi as usize + 1]) as f64) as f32;
                    if temp != 0.0 { let i = iz as usize; num[ir][i] += temp; den[ir][i] += temp * temp; nnz[ir][i] += 1.0; }
                }
                iz += 1; zn += dz;
            }
            r += dr;
        }
        if gettr(&mut tr) == 0 { gottrace = false; }
    }
    exit(CWP_Exit());
}
