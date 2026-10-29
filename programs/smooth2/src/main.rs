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

//! SMOOTH2 - smooth a uniformly sampled 2d array via damped least squares.
//! Port of par/main/smoothing/smooth2.c (Revision 1.17).
//! Credits: CWP: Zhen-yue Liu; adapted for par/main by John Stockwell,
//! 1 Oct 92. Windowing feature added by Zliu on 16 Nov 1992.

use oxwave_core::prelude::*;
use oxwave_tridiag::tripd;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SMOOTH2 --- SMOOTH a uniformly sampled 2d array of data, within a user-",
    "\t\tdefined window, via a damped least squares technique\t",
    "\t\t\t\t\t\t\t\t\t",
    " smooth2 < stdin n1= n2= [optional parameters ] > stdout\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=\t\t\tnumber of samples in the 1st (fast) dimension\t",
    " n2=\t\t\tnumber of samples in the 2nd (slow) dimension\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " r1=0\t\t\tsmoothing parameter in the 1 direction\t\t",
    " r2=0\t\t\tsmoothing parameter in the 2 direction\t\t",
    " win=0,n1,0,n2\t\tarray for window range\t\t\t\t",
    " rw=0\t\t\tsmoothing parameter for window function\t\t",
    " efile=                 =efilename if set write relative error(x1) to\t",
    "\t\t\t\tefilename\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Larger r1 and r2 result in a smoother data. Recommended ranges of r1 \t",
    " and r2 are from 1 to 20.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The file verror gives the relative error between the original velocity ",
    " and the smoothed one, as a function of depth. If the error is\t\t",
    " between 0.01 and 0.1, the smoothing parameters are suitable. Otherwise,",
    " consider increasing or decreasing the smoothing parameter values.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Smoothing can be implemented in a selected window. The range of 1st   ",
    " dimension for window is from win[0] to win[1]; the range of 2nd   \t",
    " dimension is from win[2] to win[3]. \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Smoothing the window function (i.e. blurring the edges of the window)\t",
    " may be done by setting a nonzero value for rw, otherwise the edges\t",
    " of the window will be sharp.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut n1 = 0i32;
    let mut n2 = 0i32;
    let mut r1 = 0.0f32;
    let mut r2 = 0.0f32;
    let mut rw = 0.0f32;
    let infp = stdin();
    let outfp = stdout();
    let mut efile = String::new();

    /* hook up getpar */
    initargs(&argv);
    requestdoc(0, SDOC);

    /* get parameters */
    if getparint("n1", &mut n1) == 0 {
        err!("must specify n1!\n");
    }
    if getparint("n2", &mut n2) == 0 {
        err!("must specify n2!\n");
    }
    if getparfloat("r1", &mut r1) == 0 {
        r1 = 0.0;
    }
    if getparfloat("r2", &mut r2) == 0 {
        r2 = 0.0;
    }

    /* scale the smoothing parameter */
    r1 = ((r1 * r1) as f64 * 0.25) as f32;
    r2 = ((r2 * r2) as f64 * 0.25) as f32;

    // These dimensions make the C program dereference a failed/empty
    // alloc2float result. Preserve its SIGSEGV instead of a Rust panic.
    if n1 < 0 || n2 < 0 || (n2 == 0 && n1 > 0) {
        unsafe extern "C" {
            fn signal(sig: i32, handler: usize) -> usize;
            fn raise(sig: i32) -> i32;
        }
        unsafe {
            signal(11, 0);
            raise(11);
        }
        exit(128 + 11);
    }

    /* allocate space (the C two-dimensional arrays are contiguous) */
    let nmax = if n1 < n2 { n2 } else { n1 } as usize;
    let n1 = n1 as usize;
    let n2 = n2 as usize;
    let mut win = vec![0i32; 4];
    let mut v = vec![0.0f32; n1 * n2];
    let mut v0 = vec![0.0f32; n1 * n2];
    let mut w = vec![0.0f32; n1 * n2];
    let mut errz = vec![0.0f32; nmax];
    let mut d = vec![0.0f32; nmax];
    let mut e = vec![0.0f32; nmax];
    let mut f = vec![0.0f32; nmax];

    /* read velocities */
    fread(&mut v, &infp);

    /* save the original velocity */
    for ix in 0..n2 {
        for iz in 0..n1 {
            v0[ix * n1 + iz] = v[ix * n1 + iz];
        }
    }

    /* get parameters for window function */
    if getparint("win", &mut win) == 0 {
        win[0] = 0;
        win[1] = n1 as i32;
        win[2] = 0;
        win[3] = n2 as i32;
    }
    if getparfloat("rw", &mut rw) == 0 {
        rw = 0.0;
    }
    rw = ((rw * rw) as f64 * 0.25) as f32;

    /* define the window function */
    for ix in 0..n2 {
        for iz in 0..n1 {
            w[ix * n1 + iz] = 0.0;
        }
    }
    for ix in win[2]..win[3] {
        for iz in win[0]..win[1] {
            w[ix as usize * n1 + iz as usize] = 1.0;
        }
    }

    if win[0] > 0 || win[1] < n1 as i32 || win[2] > 0 || win[3] < n2 as i32 {
        /* smooth the window function */
        for iz in 0..n1 {
            for ix in 0..n2 {
                d[ix] = (1.0 + 2.0 * rw as f64) as f32;
                e[ix] = -rw;
                f[ix] = w[ix * n1 + iz];
            }
            d[0] -= rw;
            d[n2 - 1] -= rw;
            tripd(&mut d[..n2], &mut e, &mut f);
            for ix in 0..n2 {
                w[ix * n1 + iz] = f[ix];
            }
        }
        // With n1 == 0 the C only touches memory outside its arrays;
        // it copies no window or velocity samples back.
        if n1 > 0 {
            for ix in 0..n2 {
                for iz in 0..n1 {
                    d[iz] = (1.0 + 2.0 * rw as f64) as f32;
                    e[iz] = -rw;
                    f[iz] = w[ix * n1 + iz];
                }
                d[0] -= rw;
                d[n1 - 1] -= rw;
                tripd(&mut d[..n1], &mut e, &mut f);
                for iz in 0..n1 {
                    w[ix * n1 + iz] = f[iz];
                }
            }
        }
    }

    /* solving for the smoothing velocity */
    for iz in 0..n1 {
        for ix in 0..n2 - 1 {
            d[ix] = (1.0 + (r2 * (w[ix * n1 + iz] + w[(ix + 1) * n1 + iz])) as f64) as f32;
            e[ix] = -r2 * w[(ix + 1) * n1 + iz];
            f[ix] = v[ix * n1 + iz];
        }
        d[0] -= r2 * w[iz];
        d[n2 - 1] = (1.0 + (r2 * w[(n2 - 1) * n1 + iz]) as f64) as f32;
        f[n2 - 1] = v[(n2 - 1) * n1 + iz];
        tripd(&mut d[..n2], &mut e, &mut f);
        for ix in 0..n2 {
            v[ix * n1 + iz] = f[ix];
        }
    }
    // For n1 <= 1 the original fast-axis block writes outside its scratch
    // arrays (tripd gets n <= 0), but copies no samples into v. Omit those
    // invalid scratch writes while retaining the observed output.
    if n1 > 1 {
        for ix in 0..n2 {
            for iz in 0..n1 - 2 {
                d[iz] = (1.0 + (r1 * (w[ix * n1 + iz + 1] + w[ix * n1 + iz + 2])) as f64) as f32;
                e[iz] = -r1 * w[ix * n1 + iz + 2];
                f[iz] = v[ix * n1 + iz + 1];
            }
            f[0] += r1 * w[ix * n1 + 1] * v[ix * n1];
            d[n1 - 2] = (1.0 + (r1 * w[ix * n1 + n1 - 1]) as f64) as f32;
            f[n1 - 2] = v[ix * n1 + n1 - 1];
            tripd(&mut d[..n1 - 1], &mut e, &mut f);
            for iz in 0..n1 - 1 {
                v[ix * n1 + iz + 1] = f[iz];
            }
        }
    }

    /* write smoothed data */
    fwrite(&v, &outfp);

    /* if the user specifies the name of an error file */
    getparstring("efile", &mut efile);
    checkpars();
    if !efile.is_empty() {
        let errorfp = efopen(&efile, "w");

        /* calculate the RMS error of velocity */
        let mut err0 = 0.0f32;
        let mut vrms = 0.0f32;
        for iz in 0..n1 {
            for ix in 0..n2 {
                err0 += (v0[ix * n1 + iz] - v[ix * n1 + iz]) * (v0[ix * n1 + iz] - v[ix * n1 + iz]);
                vrms += v0[ix * n1 + iz] * v0[ix * n1 + iz];
            }
            errz[iz] = if vrms != 0.0 {
                ((err0 / vrms) as f64).sqrt() as f32
            } else {
                ((err0 / f32::MIN_POSITIVE) as f64).sqrt() as f32
            };
        }
        fwrite(&errz[..n1], &errorfp);
        fclose(&errorfp);
    }
    exit(CWP_Exit());
}
