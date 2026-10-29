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

//! SUVELAN - port of su/main/velocity_analysis/suvelan.c (Revision 1.20).
//! Credits: CWP: Dave Hale and Bjoern Rommel; SINTEF, IKU Petroleumsforskning:
//! Bjoorn Rommel (power-of-semblance function).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t     ",
    " SUVELAN - compute stacking velocity semblance for cdp gathers\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " suvelan <stdin >stdout [optional parameters]\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Optional Parameters:\t\t\t\t\t\t\t     ",
    " nv=50                   number of velocities\t\t\t\t     ",
    " dv=50.0                 velocity sampling interval\t\t\t     ",
    " fv=1500.0               first velocity\t\t\t\t     ",
    " anis1=0.0               quartic term, numerator of an extended quartic term",
    " anis2=0.0               in denominator of an extended quartic term         ",
    " smute=1.5               samples with NMO stretch exceeding smute are zeroed",
    " dtratio=5               ratio of output to input time sampling intervals   ",
    " nsmooth=dtratio*2+1     length of semblance num and den smoothing window   ",
    " verbose=0               =1 for diagnostic print on stderr\t\t     ",
    " pwr=1.0                 semblance value to the power      \t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Notes:\t\t\t\t\t\t\t\t     ",
    " Velocity analysis is usually a two-dimensional screen for optimal values of",
    " the vertical two-way traveltime and stacking velocity. But if the travel-  ",
    " time curve is no longer close to a hyperbola, the quartic term of the      ",
    " traveltime series should be considered. In its easiest form (with anis2=0) ",
    " the optimizion of all parameters requires a three-dimensional screen. This ",
    " is done by a repetition of the conventional two-dimensional screen with a  ",
    " variation of the quartic term. The extended quartic term is more accurate, ",
    " though the function is no more a polynomial. When screening for optimal    ",
    " values the theoretical dependencies between these paramters can be taken   ",
    " into account. The traveltime function is defined by                        ",
    "                                                                            ",
    "                1            anis1                                          ",
    " t^2 = t_0^2 + --- x^2 + ------------- x^4                                  ",
    "               v^2       1 + anis2 x^2                                      ",
    "                                                                            ",
    " The coefficients anis1, anis2 are assumed to be small, that means the non- ",
    " hyperbolicity is assumed to be small. Triplications cannot be handled.     ",
    "                                                                            ",
    " Semblance is defined by the following quotient:\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    "                 n-1                 \t\t\t\t\t     ",
    "               [ sum q(t,j) ]^2      \t\t\t\t\t     ",
    "                 j=0                 \t\t\t\t\t     ",
    "       s(t) = ------------------     \t\t\t\t\t     ",
    "                 n-1                 \t\t\t\t\t     ",
    "               n sum [q(t,j)]^2      \t\t\t\t\t     ",
    "                 j=0                 \t\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " where n is the number of non-zero samples after muting.\t\t     ",
    " Smoothing (nsmooth) is applied separately to the numerator and denominator ",
    " before computing this semblance quotient.\t\t\t\t     ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Then, the semblance is set to the power of the parameter pwr. With pwr > 1 ",
    " the difference between semblance values is stretched in the upper half of  ",
    " the range of semblance values [0,1], but compressed in the lower half of   ",
    " it; thus, the few large semblance values are enhanced. With pwr < 1 the    ",
    " many small values are enhanced, thus more discernible against background   ",
    " noise. Of course, always at the expanse of the respective other feature.   ",
    "\t\t\t\t\t\t\t\t\t     ",
    " Input traces should be sorted by cdp - suvelan outputs a group of\t     ",
    " semblance traces every time cdp changes.  Therefore, the output will\t     ",
    " be useful only if cdp gathers are input.\t\t\t\t     ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(0, SDOC);
    if gettr(&mut tr) == 0 { err!("can't get first trace"); }
    let nt = tr.ns as i32;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let ft = (tr.delrt as f64 / 1000.0) as f32;
    let mut nv = 50i32;
    let mut dv = 50.0f32;
    let mut fv = 1500.0f32;
    let mut anis1 = 0.0f32;
    let mut anis2 = 0.0f32;
    let mut smute = 1.5f32;
    let mut dtratio = 5i32;
    let mut verbose = 0i32;
    let mut pwr = 1.0f32;
    getparint("nv", &mut nv);
    getparfloat("dv", &mut dv);
    getparfloat("fv", &mut fv);
    getparfloat("anis1", &mut anis1);
    getparfloat("anis2", &mut anis2);
    getparfloat("smute", &mut smute);
    if smute <= 1.0 { err!("smute must be greater than 1.0"); }
    getparint("dtratio", &mut dtratio);
    let mut nsmooth = dtratio.wrapping_mul(2).wrapping_add(1);
    getparint("nsmooth", &mut nsmooth);
    getparint("verbose", &mut verbose);
    getparfloat("pwr", &mut pwr);
    if pwr < 0.0 { err!("we are not looking for noise: pwr < 0"); }
    if pwr == 0.0 { err!("we are creating an all-white semblance: pwr = 0"); }
    checkpars();
    // arm64 C integer division by zero produces zero.
    let ntout = 1 + (nt - 1).checked_div(dtratio).unwrap_or(0);
    CHECK_NT("ntout", ntout);
    let dtout = dt * dtratio as f32;
    if verbose != 0 {
        fprintf!(stderr(), "\tnumber of output time samples is %d\n", ntout);
        fprintf!(stderr(), "\toutput time sampling interval is %g\n", dtout);
        fprintf!(stderr(), "\toutput time of first sample is %g\n", ft);
    }
    let mut data = vec![0.0f32; nt as usize];
    // The C allocator rejects negative dimensions after conversion to size_t.
    if nv < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut num = vec![vec![0.0f32; nt as usize]; nv as usize];
    let mut den = vec![vec![0.0f32; nt as usize]; nv as usize];
    let mut nnz = vec![vec![0.0f32; nt as usize]; nv as usize];
    if ntout < 0 { err!("ealloc.c: malloc failed (Cannot allocate memory)"); }
    let mut sem = vec![0.0f32; ntout as usize];
    let mut gottrace = true;
    let mut cdpprev = tr.cdp;
    let mut cdp = tr.cdp;
    let mut offset = tr.offset as f32;
    let mut offan = 0.0f32;
    loop {
        if gottrace {
            offset = tr.offset as f32;
            cdp = tr.cdp;
            if 1.0 + (offset * offset * anis2) as f64 <= 0.0 { err!("anis2 too small"); }
            offan = ((offset * offset * offset * offset * anis1) as f64 /
                (1.0 + (offset * offset * anis2) as f64)) as f32;
            data.copy_from_slice(&tr.data[..nt as usize]);
        }
        if cdp != cdpprev || !gottrace {
            tr.offset = 0;
            tr.cdp = cdpprev;
            tr.ns = ntout as u16;
            tr.dt = (dtout as f64 * 1000000.0) as u16;
            for iv in 0..nv as usize {
                for itout in 0..ntout {
                    let it = itout * dtratio;
                    let ismin = (it - nsmooth / 2).max(0);
                    let ismax = (it + nsmooth / 2).min(nt - 1);
                    let mut nsum = 0.0f32;
                    let mut dsum = 0.0f32;
                    for is in ismin..ismax {
                        let i = is as usize;
                        nsum += num[iv][i] * num[iv][i];
                        dsum += nnz[iv][i] * den[iv][i];
                    }
                    sem[itout as usize] = if dsum != 0.0 { nsum / dsum } else { 0.0 };
                }
                if pwr != 1.0 {
                    for s in &mut sem { *s = (*s as f64).powf(pwr as f64) as f32; }
                }
                tr.data[..ntout as usize].copy_from_slice(&sem);
                puttr(&mut tr);
                num[iv].fill(0.0);
                den[iv].fill(0.0);
                nnz[iv].fill(0.0);
            }
            if verbose != 0 { warn!("semblance output for cdp=%d", cdpprev); }
            if !gottrace { break; }
            cdpprev = cdp;
        }
        let mut v = fv;
        for iv in 0..nv as usize {
            let offovs = (offset * offset) / (v * v) + offan;
            if offovs < 0.0 { warn!("no moveout; check anis1 and anis2"); }
            let tnmute = ((offovs as f64 / ((smute * smute) as f64 - 1.0)).sqrt()) as f32;
            let itmute = if tnmute > ft { ((tnmute - ft) / dt) as i32 } else { 0 };
            let mut tn = ft + itmute as f32 * dt;
            for it in itmute..nt {
                let ti = ((((tn * tn + offovs) as f64).sqrt() - ft as f64) / dt as f64) as f32;
                let iti = ti as i32;
                if iti < nt - 1 {
                    if iti < 0 { err!("interpolation index out of bounds"); }
                    let frac = ti - iti as f32;
                    let temp = ((1.0 - frac as f64) * data[iti as usize] as f64 +
                        (frac * data[(iti + 1) as usize]) as f64) as f32;
                    if temp != 0.0 {
                        num[iv][it as usize] += temp;
                        den[iv][it as usize] += temp * temp;
                        nnz[iv][it as usize] += 1.0;
                    }
                }
                tn += dt;
            }
            v += dv;
        }
        if gettr(&mut tr) == 0 { gottrace = false; }
    }
    exit(CWP_Exit());
}
