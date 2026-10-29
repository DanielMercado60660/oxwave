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

//! UNI2TRI - convert uniformly sampled model to a triangulated model.
//! Port of tri/main/uni2tri.c (Revision 1.7).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Author: Craig Artley, Colorado School of Mines, 03/31/94.
//! NOTE: After a program outlined by Dave Hale, 12/27/90.

#![allow(non_snake_case)]

use oxwave_core::prelude::*;
use oxwave_lib_tri::*;
use oxwave_sort::*;
use std::ffi::c_void;
use std::ptr::null_mut;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " UNI2TRI - convert UNIformly sampled model to a TRIangulated model\t",
    "\t\t\t\t\t\t\t\t\t",
    " uni2tri <slothfile >modelfile n2= n1= [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=                    number of samples in first (fast) dimension\t",
    " n2=                    number of samples in second dimension\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " d1=1.0                 sampling interval in dimension 1\t\t",
    " d2=1.0                 sampling interval in dimension 2\t\t",
    " f1=0.0                 first value in dimension 1\t\t\t",
    " f2=0.0                 first value in dimesion 2\t\t\t",
    " ifile=                 triangulated model file used as initial model\t",
    " errmax=                maximum sloth error (see notes below)\t\t",
    " verbose=1              =0 for silence\t\t\t\t\t",
    "                        =1 to report maximum error at each stage to stderr",
    "                        =2 to also write the normalized error to efile\t",
    " efile=emax.dat         filename for error file (for verbose=2)\t",
    " mm=0\t\t\toutput every mm-th intermediate model (0 for none)",
    " mfile=intmodel        intermediate models written to intmodel%d\t",
    " method=3              =1 add 1 vertex at maximum error\t\t",
    "                       =2 add vertex to every triangle that exceeds errmax",
    "                       =3 method 2, but avoid closely spaced vertices\t",
    " tol=10                closeness criterion for (in samples)\t\t",
    " sfill=                 x, z, x0, z0, s00, dsdx, dsdz to fill a region\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Triangles are constructed until the maximum error is\t\t\t",
    " not greater than the user-specified errmax.  The default errmax\t",
    " is 1% of the maximum value in the sampled input file.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " After the uniform values have been triangulated, the optional sfill\t",
    " parameters are used to fill closed regions bounded by fixed edges.\t",
    " Let (x,z) denote any point inside a closed region.  Values inside\t",
    " this region is determined by s(x,z) = s00+(x-x0)*dsdx+(z-z0)*dsdz.\t",
    " The (x,z) component of the sfill parameter is used to identify a\t",
    " closed region.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The uniformly sampled quantity is assumed to be sloth=1/v^2.\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// Triangle I/O uses native FILE pointers. Duplicate descriptors so the native
// streams can be closed without closing the corresponding core streams.
extern "C" {
    fn fdopen(fd: i32, mode: *const i8) -> *mut c_void;
    fn dup(fd: i32) -> i32;
    #[link_name = "fclose"]
    fn native_fclose(fp: *mut c_void) -> i32;
    fn calloc(n: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn signal(sig: i32, handler: usize) -> usize;
    fn raise(sig: i32) -> i32;
}

#[repr(C)]
struct VertexAttributes {
    s: f32,
}
#[repr(C)]
struct FaceAttributes {
    s00: f32,
    dsdx: f32,
    dsdz: f32,
    dens: f32,
    qfac: f32,
}

// The sampled-grid arguments of the C file's private helpers.
struct Grid {
    n2: i32,
    d2: f32,
    f2: f32,
    n1: i32,
    d1: f32,
    f1: f32,
    s: Vec<f32>,
}
impl Grid {
    fn sample(&self, ix: i32, iz: i32) -> f32 {
        let i = ix * self.n1 + iz;
        // For n1=0 the reference reads zero beyond a zero-size allocation
        // (including [-1]). Use zero explicitly for this undefined read.
        // n2=0 also has an uninitialized row pointer in C and is unstable.
        if self.n1 == 0 || self.n2 == 0 {
            return 0.0;
        }
        self.s[i as usize]
    }
}

unsafe fn attributes<T>() -> *mut T {
    unsafe { calloc(1, std::mem::size_of::<T>()).cast() }
}

unsafe fn write_model(m: *mut Model, fp: &CFile) {
    unsafe {
        efflush(fp);
        let native = fdopen(dup(fileno(fp)), c"w".as_ptr());
        writeModel(m, native);
        native_fclose(native);
    }
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC);
    let (mut n2, mut n1) = (0i32, 0i32);
    if getparint("n2", &mut n2) == 0 {
        err!("must specify n2!");
    }
    if getparint("n1", &mut n1) == 0 {
        err!("must specify n1!");
    }

    // Read native-endian sloth before parsing optional parameters.
    if n1 < 0 || n2 < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    let count = n1.wrapping_mul(n2) as usize;
    let mut s = vec![0.0f32; count];
    if fread(&mut s, &stdin()) != count {
        err!("error reading input file!");
    }
    // The observed n1=0 reference runs read zero from their unused block.
    let mut smax = if count == 0 { 0.0 } else { s[0] };
    for ix in 0..n2 {
        for iz in 0..n1 {
            if s[(ix * n1 + iz) as usize] > smax {
                smax = s[(ix * n1 + iz) as usize];
            }
        }
    }
    let (mut d2, mut d1, mut f2, mut f1) = (1.0f32, 1.0f32, 0.0f32, 0.0f32);
    getparfloat("d2", &mut d2);
    getparfloat("d1", &mut d1);
    getparfloat("f2", &mut f2);
    getparfloat("f1", &mut f1);
    let mut errmax = (0.01 * smax as f64) as f32;
    getparfloat("errmax", &mut errmax);
    let mut efile = "emax.dat".to_string();
    getparstring("efile", &mut efile);
    let mut mm = 0i32;
    getparint("mm", &mut mm);
    mm = if mm > 0 { mm } else { 0 };
    let mut mfile = "intmodel".to_string();
    getparstring("mfile", &mut mfile);
    let mut method = 3i32;
    getparint("method", &mut method);
    let mut tol = 10i32;
    getparint("tol", &mut tol);
    let mut verbose = 1i32;
    getparint("verbose", &mut verbose);

    let nr = countparname("sfill");
    let mut sr = vec![vec![0.0f32; 7]; nr as usize];
    for ir in 0..nr {
        if getnparfloat(ir + 1, "sfill", &mut sr[ir as usize]) != 7 {
            err!("7 values must be specified in sfill parameter");
        }
    }
    if verbose != 0 {
        fprintf!(
            stderr(),
            "n2=%d d2=%g f2=%g\nn1=%d d1=%g f1=%g\nerrmax=%g\n",
            n2,
            d2,
            f2,
            n1,
            d1,
            f1,
            errmax
        );
    }
    let g = Grid {
        n2,
        d2,
        f2,
        n1,
        d1,
        f1,
        s,
    };
    unsafe {
        let m = match method {
            1 => makeMod1(errmax, &g, smax, &efile, mm, &mfile, verbose),
            2 => makeMod2(errmax, &g, smax, &efile, mm, &mfile, verbose),
            3 => makeMod3(errmax, &g, smax, &efile, mm, &mfile, tol, verbose),
            _ => err!("Unknown method=%d!\n", method),
        };
        fillsloth(m, &sr);
        write_model(m, &stdout());
    }
    exit(CWP_Exit());
}

// Identical initialization block in makeMod1, makeMod2 and makeMod3.
unsafe fn initialModel(g: &Grid) -> *mut Model {
    unsafe {
        let Grid {
            n2,
            d2,
            f2,
            n1,
            d1,
            f1,
            ..
        } = *g;
        let xmin = f2;
        let zmin = f1;
        let xmax = f2 + (n2 - 1) as f32 * d2;
        let zmax = f1 + (n1 - 1) as f32 * d1;
        let mut ifile = String::new();
        let m;
        if getparstring("ifile", &mut ifile) != 0 {
            let ifp = efopen(&ifile, "r");
            let native = fdopen(dup(fileno(&ifp)), c"r".as_ptr());
            m = readModel(native);
            native_fclose(native);
            efclose(&ifp);
            // Empty initial files return NULL in the C library, which the
            // original program immediately dereferences.
            if m.is_null() {
                signal(11, 0);
                raise(11);
                exit(139);
            }
            if ((xmin - (*m).ymin).abs() as f64) > 0.01 * d2 as f64 {
                err!("xmin!=m->ymin %g!=%g", xmin, (*m).ymin);
            }
            if ((xmax - (*m).ymax).abs() as f64) > 0.01 * d2 as f64 {
                err!("xmax!=m->ymax %g!=%g", xmax, (*m).ymax);
            }
            if ((zmin - (*m).xmin).abs() as f64) > 0.01 * d1 as f64 {
                err!("zmin!=m->xmin %g!=%g", zmin, (*m).xmin);
            }
            if ((zmax - (*m).xmax).abs() as f64) > 0.01 * d1 as f64 {
                err!("zmax!=m->xmax %g!=%g", zmax, (*m).xmax);
            }
            let mut f = (*m).f;
            loop {
                free((*f).fa);
                (*f).fa = null_mut();
                setSlothOfVertex((*(*(*f).eu).vu).v, g);
                setSlothOfVertex((*(*(*(*f).eu).euCW).vu).v, g);
                setSlothOfVertex((*(*(*(*f).eu).euCCW).vu).v, g);
                f = (*f).fNext;
                if f == (*m).f {
                    break;
                }
            }
        } else {
            m = makeModel(zmin, xmin, zmax, xmax);
            (*m).sfa = std::mem::size_of::<FaceAttributes>();
            (*m).eps = (0.01 * (if d2 < d1 { d2 } else { d1 }) as f64) as f32;
            let v = nearestVertexInModel(m, null_mut(), zmin, xmin);
            let va = attributes::<VertexAttributes>();
            (*va).s = g.sample(0, 0);
            (*v).va = va.cast();
            let v = nearestVertexInModel(m, null_mut(), zmax, xmin);
            let va = attributes::<VertexAttributes>();
            (*va).s = g.sample(0, n1 - 1);
            (*v).va = va.cast();
            let v = nearestVertexInModel(m, null_mut(), zmin, xmax);
            let va = attributes::<VertexAttributes>();
            (*va).s = g.sample(n2 - 1, 0);
            (*v).va = va.cast();
            let v = nearestVertexInModel(m, null_mut(), zmax, xmax);
            let va = attributes::<VertexAttributes>();
            (*va).s = g.sample(n2 - 1, n1 - 1);
            (*v).va = va.cast();
        }
        m
    }
}

// Final face-attribute block shared verbatim by the three C methods.
unsafe fn finalSloth(m: *mut Model) {
    unsafe {
        let mut f = (*m).f;
        loop {
            let fa = attributes::<FaceAttributes>();
            (*f).fa = fa.cast();
            makeSlothForTri(f, &mut (*fa).s00, &mut (*fa).dsdx, &mut (*fa).dsdz);
            (*fa).dens = f32::MAX;
            (*fa).qfac = f32::MAX;
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
    }
}

unsafe fn makeMod1(
    errmax: f32,
    g: &Grid,
    smax: f32,
    efile: &str,
    mm: i32,
    mfile: &str,
    verbose: i32,
) -> *mut Model {
    unsafe {
        let m = initialModel(g);
        let efp = if verbose == 2 {
            Some(fopen(efile, "w").unwrap_or_else(|| err!("Could not open efile=%s!\n", efile)))
        } else {
            None
        };
        if verbose != 0 {
            fprintf!(stderr(), "iterations faces  max errror    norm error\n");
            fprintf!(stderr(), "--------   -----  ----------    ----------\n");
        }
        let mut im = 0i32;
        let (mut ixmax, mut izmax) = (0i32, 0i32);
        let (mut ix, mut iz) = (0i32, 0i32);
        loop {
            let mut emax = 0.0f32;
            let mut f = (*m).f;
            let mut nf = 0i32;
            loop {
                let mut e = 0.0f32;
                maxErrorInTri(f, g, &mut e, &mut ix, &mut iz);
                if e >= emax {
                    emax = e;
                    ixmax = ix;
                    izmax = iz;
                }
                f = (*f).fNext;
                nf += 1;
                if f == (*m).f {
                    break;
                }
            }
            let errfrac = emax / smax;
            if verbose != 0 {
                fprintf!(
                    stderr(),
                    " %5d     %5d   % .6g     % .6g\n",
                    im,
                    nf,
                    emax,
                    errfrac
                );
            }
            if let Some(ref fp) = efp {
                fwrite(&[errfrac], fp);
                if im % 50 == 0 {
                    efflush(fp);
                }
            }
            if mm > 0 && im % mm == 0 {
                writeIntermediateModel(m, im, mfile);
            }
            if emax > errmax {
                let v = addVertexToModel(m, g.f1 + izmax as f32 * g.d1, g.f2 + ixmax as f32 * g.d2);
                if v.is_null() {
                    err!("This can't happen!");
                }
                let va = attributes::<VertexAttributes>();
                (*va).s = g.sample(ixmax, izmax);
                (*v).va = va.cast();
            }
            im += 1;
            if !(emax > errmax) {
                break;
            }
        }
        if let Some(fp) = efp {
            fclose(&fp);
        }
        finalSloth(m);
        m
    }
}

unsafe fn makeMod2(
    errmax: f32,
    g: &Grid,
    smax: f32,
    efile: &str,
    mm: i32,
    mfile: &str,
    verbose: i32,
) -> *mut Model {
    unsafe {
        let m = initialModel(g);
        let efp = if verbose == 2 {
            Some(fopen(efile, "w").unwrap_or_else(|| err!("Could not open efile=%s!\n", efile)))
        } else {
            None
        };
        if verbose != 0 {
            fprintf!(stderr(), "iterations faces  max errror    norm error\n");
            fprintf!(stderr(), "--------   -----  ----------    ----------\n");
        }
        let mut nfalloc = 1024;
        let mut e = vec![0.0f32; nfalloc];
        let mut ix = vec![0i32; nfalloc];
        let mut iz = vec![0i32; nfalloc];
        let mut im = 0i32;
        loop {
            let mut emax = 0.0f32;
            let mut f = (*m).f;
            let mut jf = 0;
            loop {
                maxErrorInTri(f, g, &mut e[jf], &mut ix[jf], &mut iz[jf]);
                if e[jf] >= emax {
                    emax = e[jf];
                }
                f = (*f).fNext;
                jf += 1;
                if f == (*m).f {
                    break;
                }
            }
            let nf = jf;
            let errfrac = emax / smax;
            if verbose != 0 {
                fprintf!(
                    stderr(),
                    " %5d     %5d   % .6g     % .6g\n",
                    im,
                    nf,
                    emax,
                    errfrac
                );
            }
            if let Some(ref fp) = efp {
                fwrite(&[errfrac], fp);
                if im % 50 == 0 {
                    fflush(fp);
                }
            }
            if mm > 0 && im % mm == 0 {
                writeIntermediateModel(m, im, mfile);
            }
            let mut iv = 0;
            for jf in 0..nf {
                if e[jf] > errmax {
                    let v = addVertexToModel(
                        m,
                        g.f1 + iz[jf] as f32 * g.d1,
                        g.f2 + ix[jf] as f32 * g.d2,
                    );
                    if v.is_null() {
                        continue;
                    }
                    let va = attributes::<VertexAttributes>();
                    (*va).s = g.sample(ix[jf], iz[jf]);
                    (*v).va = va.cast();
                    iv += 1;
                }
            }
            if iv == 0 && emax > errmax {
                err!("added 0 vertices!");
            }
            if nf + 2 * iv > nfalloc {
                nfalloc *= 2;
                e = vec![0.0; nfalloc];
                ix = vec![0; nfalloc];
                iz = vec![0; nfalloc];
            }
            im += 1;
            if !(emax > errmax) {
                break;
            }
        }
        if let Some(fp) = efp {
            fclose(&fp);
        }
        finalSloth(m);
        m
    }
}

unsafe fn makeMod3(
    errmax: f32,
    g: &Grid,
    smax: f32,
    efile: &str,
    mm: i32,
    mfile: &str,
    tol: i32,
    verbose: i32,
) -> *mut Model {
    unsafe {
        let m = initialModel(g);
        checkpars();
        let efp = if verbose == 2 {
            Some(fopen(efile, "w").unwrap_or_else(|| err!("Could not open efile=%s!\n", efile)))
        } else {
            None
        };
        if verbose != 0 {
            fprintf!(stderr(), "iterations faces  max errror    norm error\n");
            fprintf!(stderr(), "--------   -----  ----------    ----------\n");
        }
        let mut nfalloc = 1024;
        let mut e = vec![0.0f32; nfalloc];
        let mut ix = vec![0i32; nfalloc];
        let mut iz = vec![0i32; nfalloc];
        let mut ind2 = vec![0i32; nfalloc];
        let mut used = vec![0i32; nfalloc];
        let (mut im, mut nv) = (0i32, 0i32);
        loop {
            let mut emax = 0.0f32;
            let mut f = (*m).f;
            let mut jf = 0;
            loop {
                maxErrorInTri(f, g, &mut e[jf], &mut ix[jf], &mut iz[jf]);
                if e[jf] >= emax {
                    emax = e[jf];
                }
                f = (*f).fNext;
                jf += 1;
                if f == (*m).f {
                    break;
                }
            }
            let nf = jf;
            let errfrac = emax / smax;
            if verbose != 0 {
                fprintf!(
                    stderr(),
                    " %5d     %5d   % .6g     % .6g\n",
                    im,
                    nf,
                    emax,
                    errfrac
                );
            }
            if let Some(ref fp) = efp {
                fwrite(&[errfrac], fp);
                if im % 50 == 0 {
                    fflush(fp);
                }
            }
            if mm > 0 && im % mm == 0 {
                writeIntermediateModel(m, im, mfile);
            }
            for jf in 0..nf {
                ind2[jf] = jf as i32;
                used[jf] = 0;
            }
            qkisort(&e[..nf], &mut ind2[..nf]);
            for jf in 0..nf / 2 {
                let kf = nf - 1 - jf;
                ind2.swap(jf, kf);
            }
            let mut iv = 0;
            for jf in 0..nf {
                let j = ind2[jf] as usize;
                if e[j] > errmax {
                    let mut skip = false;
                    for kf in 0..jf {
                        // The C deliberately also considers candidates skipped
                        // earlier in this iteration (its used[] test is disabled).
                        let k = ind2[kf] as usize;
                        let dix = ix[k] - ix[j];
                        let diz = iz[k] - iz[j];
                        if dix.abs() <= tol && diz.abs() <= tol {
                            skip = true;
                            break;
                        }
                    }
                    if skip {
                        continue;
                    }
                    let v =
                        addVertexToModel(m, g.f1 + iz[j] as f32 * g.d1, g.f2 + ix[j] as f32 * g.d2);
                    if v.is_null() {
                        continue;
                    }
                    let va = attributes::<VertexAttributes>();
                    (*va).s = g.sample(ix[j], iz[j]);
                    (*v).va = va.cast();
                    used[jf] = 1;
                    iv += 1;
                    nv += 1;
                }
            }
            if iv == 0 && emax > errmax {
                err!("added 0 vertices!");
            }
            if nf + 2 * iv > nfalloc {
                nfalloc *= 2;
                e = vec![0.0; nfalloc];
                ix = vec![0; nfalloc];
                iz = vec![0; nfalloc];
                ind2 = vec![0; nfalloc];
                used = vec![0; nfalloc];
            }
            im += 1;
            if !(emax > errmax) {
                break;
            }
        }
        if let Some(fp) = efp {
            fclose(&fp);
        }
        if verbose != 0 {
            fprintf!(stderr(), "added %d vertices\x08", nv);
        }
        finalSloth(m);
        m
    }
}

// Dave Hale; modified by Craig Artley to handle nonzero model origins.
unsafe fn makeSlothForTri(t: *mut Tri, s00: &mut f32, dsdx: &mut f32, dsdz: &mut f32) {
    unsafe {
        let eu = (*t).eu;
        let v1 = (*(*eu).vu).v;
        let v2 = (*(*(*eu).euCW).vu).v;
        let v3 = (*(*(*eu).euCCW).vu).v;
        let (x1, z1) = ((*v1).y, (*v1).x);
        let s1 = (*(*v1).va.cast::<VertexAttributes>()).s;
        let (x2, z2) = ((*v2).y, (*v2).x);
        let s2 = (*(*v2).va.cast::<VertexAttributes>()).s;
        let (x3, z3) = ((*v3).y, (*v3).x);
        let s3 = (*(*v3).va.cast::<VertexAttributes>()).s;
        let (x2mx1, z2mz1, s2ms1) = (x2 - x1, z2 - z1, s2 - s1);
        let (x3mx1, z3mz1, s3ms1) = (x3 - x1, z3 - z1, s3 - s1);
        let det = z3mz1 * x2mx1 - z2mz1 * x3mx1;
        let a = z3mz1 * s2ms1 - z2mz1 * s3ms1;
        let b = s3ms1 * x2mx1 - s2ms1 * x3mx1;
        *dsdx = a / det;
        *dsdz = b / det;
        *s00 = s2 - (*dsdx) * x2 - (*dsdz) * z2;
    }
}

// Dave Hale; modified by Craig Artley for triangles with vertical edges.
unsafe fn maxErrorInTri(t: *mut Tri, g: &Grid, emax: &mut f32, ixmax: &mut i32, izmax: &mut i32) {
    unsafe {
        let Grid { d2, f2, d1, f1, .. } = *g;
        let (mut s00, mut dsdx, mut dsdz) = (0.0f32, 0.0f32, 0.0f32);
        makeSlothForTri(t, &mut s00, &mut dsdx, &mut dsdz);
        let eu = (*t).eu;
        let v1 = (*(*eu).vu).v;
        let v2 = (*(*(*eu).euCW).vu).v;
        let v3 = (*(*(*eu).euCCW).vu).v;
        let (mut x1, mut z1) = ((*v1).y, (*v1).x);
        let (mut x2, mut z2) = ((*v2).y, (*v2).x);
        let (mut x3, mut z3) = ((*v3).y, (*v3).x);
        if x1 > x2 {
            std::mem::swap(&mut x1, &mut x2);
            std::mem::swap(&mut z1, &mut z2);
        }
        if x1 > x3 {
            std::mem::swap(&mut x1, &mut x3);
            std::mem::swap(&mut z1, &mut z3);
        }
        if x2 > x3 {
            std::mem::swap(&mut x2, &mut x3);
            std::mem::swap(&mut z2, &mut z3);
        }
        // The double literal in each conditional promotes both alternatives.
        let slope21 = ((z2 - z1) as f64
            / if x2 != x1 {
                (x2 - x1) as f64
            } else {
                1e-10 * d2 as f64
            }) as f32;
        let slope31 = ((z3 - z1) as f64
            / if x3 != x1 {
                (x3 - x1) as f64
            } else {
                1e-10 * d2 as f64
            }) as f32;
        let slope32 = ((z3 - z2) as f64
            / if x3 != x2 {
                (x3 - x2) as f64
            } else {
                1e-10 * d2 as f64
            }) as f32;
        let mut ixlo = ((x1 - f2) / d2) as i32;
        let ixhi = ((x3 - f2) / d2) as i32;
        if f2 + ixlo as f32 * d2 < x1 {
            ixlo += 1;
        }
        let low = if slope21 < slope31 { slope21 } else { slope31 };
        let high = if slope21 > slope31 { slope21 } else { slope31 };
        let mut zlo = z1 + (f2 + ixlo as f32 * d2 - x1) * low;
        let mut zhi = z1 + (f2 + ixlo as f32 * d2 - x1) * high;
        let mut d1lo = d2 * low;
        let mut d1hi = d2 * high;
        if x1 == x2 && z1 < z2 {
            zhi = z2 + (f2 + ixlo as f32 * d2 - x2) * slope32;
            d1hi = d2 * slope32;
        }
        if x1 == x2 && z1 > z2 {
            zlo = z2 + (f2 + ixlo as f32 * d2 - x2) * slope32;
            d1lo = d2 * slope32;
        }
        *emax = 0.0;
        for ix in ixlo..=ixhi {
            let x = f2 + ix as f32 * d2;
            if x > x2 {
                if slope21 < slope31 {
                    zlo = z2 + (x - x2) * slope32;
                    d1lo = d2 * slope32;
                } else {
                    zhi = z2 + (x - x2) * slope32;
                    d1hi = d2 * slope32;
                }
            }
            let mut izlo = ((zlo - f1) / d1) as i32;
            if f1 + izlo as f32 * d1 < zlo {
                izlo += 1;
            }
            let izhi = ((zhi - f1) / d1) as i32;
            for iz in izlo..=izhi {
                let z = f1 + iz as f32 * d1;
                let mut e = s00 + x * dsdx + z * dsdz - g.sample(ix, iz);
                e = if e < 0.0 { -e } else { e };
                if e >= *emax {
                    *emax = e;
                    *ixmax = ix;
                    *izmax = iz;
                }
            }
            zlo += d1lo;
            zhi += d1hi;
        }
    }
}

// Craig Artley.
unsafe fn writeIntermediateModel(m: *mut Model, im: i32, mfile: &str) {
    unsafe {
        let mut f = (*m).f;
        loop {
            let fa = attributes::<FaceAttributes>();
            (*f).fa = fa.cast();
            makeSlothForTri(f, &mut (*fa).s00, &mut (*fa).dsdx, &mut (*fa).dsdz);
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
        let mexfile = sprintf!("%s%.4d", mfile, im);
        let mfp = efopen(&mexfile, "w");
        write_model(m, &mfp);
        efclose(&mfp);
        let mut f = (*m).f;
        loop {
            free((*f).fa);
            (*f).fa = null_mut();
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
    }
}

// Craig Artley.
unsafe fn setSlothOfVertex(v: *mut Vertex, g: &Grid) {
    unsafe {
        let xi = ((*v).y - g.f2) / g.d2;
        let zi = ((*v).x - g.f1) / g.d1;
        let (mut ix, mut iz) = (xi as i32, zi as i32);
        let (mut fracx, mut fracz) = (xi - ix as f32, zi - iz as f32);
        if ix < 0 {
            ix = 0;
            fracx = 0.0;
        }
        if ix > g.n2 - 2 {
            ix = g.n2 - 2;
            fracx = 1.0;
        }
        if iz < 0 {
            iz = 0;
            fracz = 0.0;
        }
        if iz > g.n1 - 2 {
            iz = g.n1 - 2;
            fracz = 1.0;
        }
        if (*v).va.is_null() {
            (*v).va = attributes::<VertexAttributes>().cast();
        }
        let va = (*v).va.cast::<VertexAttributes>();
        (*va).s = ((1.0 - fracx as f64)
            * ((1.0 - fracz as f64) * g.sample(ix, iz) as f64
                + (fracz * g.sample(ix, iz + 1)) as f64)
            + fracx as f64
                * ((1.0 - fracz as f64) * g.sample(ix + 1, iz) as f64
                    + (fracz * g.sample(ix + 1, iz + 1)) as f64)) as f32;
    }
}

unsafe fn fillsloth(m: *mut Model, sr: &[Vec<f32>]) {
    unsafe {
        for r in sr {
            let (x, z, x0, z0) = (r[0], r[1], r[2], r[3]);
            let (mut s00, dsdx, dsdz) = (r[4], r[5], r[6]);
            s00 -= x0 * dsdx + z0 * dsdz;
            let t = insideTriInModel(m, null_mut(), z, x);
            setsloth(t, s00, dsdx, dsdz);
        }
    }
}

unsafe fn setsloth(t: *mut Tri, s00: f32, dsdx: f32, dsdz: f32) {
    unsafe {
        let mut fa = (*t).fa.cast::<FaceAttributes>();
        if !fa.is_null() && (*fa).s00 == s00 && (*fa).dsdx == dsdx && (*fa).dsdz == dsdz {
            return;
        }
        if fa.is_null() {
            fa = attributes();
            (*t).fa = fa.cast();
        }
        (*fa).s00 = s00;
        (*fa).dsdx = dsdx;
        (*fa).dsdz = dsdz;
        let mut eu = (*t).eu;
        loop {
            if (*(*eu).e).fixed == 0 {
                setsloth((*(*eu).euMate).f, s00, dsdx, dsdz);
            }
            eu = (*eu).euCW;
            if eu == (*t).eu {
                break;
            }
        }
    }
}
