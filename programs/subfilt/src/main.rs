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

//! SUBFILT - apply Butterworth bandpass filter, port of su/main/filters/subfilt.c (Revision 1.22).
//! Original source copyright (c) Colorado School of Mines, 2011. All rights reserved.
//!
//! Credits (from the C source):
//! CWP: Dave Hale c. 1993 for bf.c subs and test drivers
//! CWP: Jack K. Cohen for su wrapper c. 1993
//! SEAM Project: Bruce Verwest 2009 added explicit pole option
//!               in a program called "subfiltpole"
//! CWP: John Stockwell (2012) combined Bruce Verwests changes
//!      into the original subfilt.
//!
//! Caveat: zerophase will not do good if trace has a spike near the end.
//! Trace header fields accessed: ns, dt, trid.

use oxwave_core::prelude::*;
use oxwave_f_butterworth::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUBFILT - apply Butterworth bandpass filter \t\t\t",
    " \t\t\t\t\t\t\t\t",
    " subfilt <stdin >stdout [optional parameters]\t\t\t",
    " \t\t\t\t\t\t\t        ",
    " Required parameters:\t\t\t\t\t\t",
    " \tif dt is not set in header, then dt is mandatory\t",
    " \t\t\t\t\t\t\t        ",
    " Optional parameters: (nyquist calculated internally)\t\t",
    " \tzerophase=1\t\t=0 for minimum phase filter \t",
    " \tlocut=1\t\t\t=0 for no low cut filter \t",
    " \thicut=1\t\t\t=0 for no high cut filter \t",
    " \tfstoplo=0.10*(nyq)\tfreq(Hz) in low cut stop band\t",
    " \tastoplo=0.05\t\tupper bound on amp at fstoplo \t",
    " \tfpasslo=0.15*(nyq)\tfreq(Hz) in low cut pass band\t",
    " \tapasslo=0.95\t\tlower bound on amp at fpasslo \t",
    " \tfpasshi=0.40*(nyq)\tfreq(Hz) in high cut pass band\t",
    " \tapasshi=0.95\t\tlower bound on amp at fpasshi \t",
    " \tfstophi=0.55*(nyq)\tfreq(Hz) in high cut stop band\t",
    " \tastophi=0.05\t\tupper bound on amp at fstophi \t",
    " \tverbose=0\t\t=1 for filter design info \t",
    " \tdt = (from header)\ttime sampling interval (sec)\t",
    " \t\t\t\t\t\t\t        ",
    " ... or  set filter by defining  poles and 3db cutoff frequencies",
    "\tnpoleselo=calculated     number of poles of the lo pass band",
    "\tnpolesehi=calculated     number of poles of the lo pass band",
    "\tf3dblo=calculated\tfrequency of 3db cutoff frequency",
    "\tf3dbhi=calculated\tfrequency of 3db cutoff frequency",
    " \t\t\t\t\t\t\t        ",
    " Notes:\t\t\t\t\t\t        ",
    " Butterworth filters were originally of interest because they  ",
    " can be implemented in hardware form through the combination of",
    " inductors, capacitors, and an amplifier. Such a filter can be ",
    " constructed in such a way as to have very small oscillations\t",
    " in the flat portion of the bandpass---a desireable attribute.\t",
    " Because the filters are composed of LC circuits, the impulse  ",
    " response is an ordinary differential equation, which translates",
    " into a polynomial in the transform domain. The filter is expressed",
    " as the division by this polynomial. Hence the poles of the filter",
    " are of interest.\t\t\t\t\t        ",
    " \t\t\t\t\t\t\t        ",
    " The user may define low pass, high pass, and band pass filters",
    " that are either minimum phase or are zero phase.  The default\t",
    " is to let the program calculate the optimal number of poles in",
    " low and high cut bands. \t\t\t\t\t",
    " \t\t\t\t\t\t\t        ",
    " Alternately the user may manually define the filter by the 3db",
    " frequency and by the number of poles in the low and or high\t",
    " cut region. \t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t        ",
    " The advantage of using the alternate method is that the user  ",
    " can control the smoothness of the filter. Greater smoothness  ",
    " through a larger pole number results in a more bell shaped    ",
    " amplitude spectrum.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t        ",
    " For simple zero phase filtering with sin squared tapering use ",
    " \"sufilter\".\t\t\t\t\t\t        ",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Get info from first trace */
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    if !seismic {
        warn!("input is not seismic data, trid=%d", tr.trid);
    }
    let nt = tr.ns as i32;
    let mut dt = 0.0f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = (tr.dt as f64 / 1000000.0) as f32;
    }
    if dt == 0.0 {
        err!("dt field is zero and not getparred");
    }
    let nyq = (0.5f64 / dt as f64) as f32;

    /* Get design frequencies and amplitudes */
    let mut verbose = 0;
    let mut zerophase = 1;
    let mut locut = 1;
    let mut hicut = 1;
    getparint("verbose", &mut verbose);
    getparint("zerophase", &mut zerophase);
    getparint("locut", &mut locut);
    getparint("hicut", &mut hicut);
    let mut fstoplo = (0.10 * nyq as f64) as f32;
    let mut fpasslo = (0.15 * nyq as f64) as f32;
    let mut fpasshi = (0.40 * nyq as f64) as f32;
    let mut fstophi = (0.55 * nyq as f64) as f32;
    getparfloat("fstoplo", &mut fstoplo);
    getparfloat("fpasslo", &mut fpasslo);
    getparfloat("fpasshi", &mut fpasshi);
    getparfloat("fstophi", &mut fstophi);
    if locut != 0 {
        if fstoplo <= 0.0 {
            err!("fstoplo must be positive");
        }
        if fstoplo > fpasslo {
            err!("fstoplo must be < fpasslo");
        }
    }
    if hicut != 0 {
        if fpasshi > fstophi {
            err!("fpasshi must be < fstophi");
        }
        if fstophi > nyq {
            err!("fstophi must be < nyquist (%f)", nyq);
        }
    }
    let mut astoplo = 0.05f32;
    let mut apasslo = 0.95f32;
    let mut apasshi = 0.95f32;
    let mut astophi = 0.05f32;
    getparfloat("astoplo", &mut astoplo);
    getparfloat("apasslo", &mut apasslo);
    getparfloat("apasshi", &mut apasshi);
    getparfloat("astophi", &mut astophi);
    if astoplo > apasslo || apasshi < astophi {
        err!("Bad amplitude parameters");
    }

    /* Normalize frequencies to [0, 0.5] for bfdesign */
    fstoplo *= dt;
    fpasslo *= dt;
    fstophi *= dt;
    fpasshi *= dt;

    /* Adapt user frequencies if zerophase selected */
    if zerophase != 0 {
        astoplo = (astoplo as f64).sqrt() as f32;
        apasslo = (apasslo as f64).sqrt() as f32;
        astophi = (astophi as f64).sqrt() as f32;
        apasshi = (apasshi as f64).sqrt() as f32;
    }

    let mut npoleslo = 0;
    let mut npoleshi = 0;
    let mut f3dblo = 0.0f32;
    let mut f3dbhi = 0.0f32;
    let is_npoleslo = getparint("npoleslo", &mut npoleslo) != 0;
    let is_npoleshi = getparint("npoleshi", &mut npoleshi) != 0;
    if is_npoleslo {
        if getparfloat("f3dblo", &mut f3dblo) == 0 {
            f3dblo = (0.15 * nyq as f64) as f32;
        }
        f3dblo *= dt;
    } else if locut != 0 {
        /* Use bfdesign to make lo cut filters */
        bfdesign(
            fpasslo,
            apasslo,
            fstoplo,
            astoplo,
            &mut npoleslo,
            &mut f3dblo,
        );
    }
    if is_npoleshi {
        if getparfloat("f3dbhi", &mut f3dbhi) == 0 {
            f3dbhi = (0.40 * nyq as f64) as f32;
        }
        f3dbhi *= dt;
    } else if hicut != 0 {
        /* Use bfdesign to make hi cut filters */
        bfdesign(
            fpasshi,
            apasshi,
            fstophi,
            astophi,
            &mut npoleshi,
            &mut f3dbhi,
        );
    }

    /* Give verbose info if requested */
    if verbose != 0 && locut != 0 {
        if zerophase != 0 {
            warn!(
                "low-cut filter: npoles = %d, 3db point = %f(Hz)",
                npoleslo.wrapping_mul(2),
                f3dblo / dt
            );
        } else {
            warn!(
                "low-cut filter: npoles = %d, 3db point = %f(Hz)",
                npoleslo,
                f3dblo / dt
            );
        }
    }
    if verbose != 0 && hicut != 0 {
        if zerophase != 0 {
            warn!(
                "high-cut filter: npoles = %d, 3db point = %f(Hz)",
                npoleshi.wrapping_mul(2),
                f3dbhi / dt
            );
        } else {
            warn!(
                "high-cut filter: npoles = %d, 3db point = %f(Hz)",
                npoleshi,
                f3dbhi / dt
            );
        }
    }

    // C aliases the input and output of each filter call. The Rust library
    // takes separate slices, so preserve the pass input in a reusable buffer.
    let mut input = vec![0.0f32; nt as usize];

    // subfilt.c does not call checkpars().
    /* Main loop over traces */
    loop {
        /* low-cut (high pass) filter */
        if locut != 0 {
            input.copy_from_slice(&tr.data[..nt as usize]);
            bfhighpass(npoleslo, f3dblo, &input, &mut tr.data[..nt as usize]);
            if zerophase != 0 {
                for i in 0..nt / 2 {
                    /* reverse trace in place */
                    let tmp = tr.data[i as usize];
                    tr.data[i as usize] = tr.data[(nt - 1 - i) as usize];
                    tr.data[(nt - 1 - i) as usize] = tmp;
                }
                input.copy_from_slice(&tr.data[..nt as usize]);
                bfhighpass(npoleslo, f3dblo, &input, &mut tr.data[..nt as usize]);
                for i in 0..nt / 2 {
                    /* flip trace back */
                    let tmp = tr.data[i as usize];
                    tr.data[i as usize] = tr.data[(nt - 1 - i) as usize];
                    tr.data[(nt - 1 - i) as usize] = tmp;
                }
            }
        }

        /* high-cut (low pass) filter */
        if hicut != 0 {
            input.copy_from_slice(&tr.data[..nt as usize]);
            bflowpass(npoleshi, f3dbhi, &input, &mut tr.data[..nt as usize]);
            if zerophase != 0 {
                for i in 0..nt / 2 {
                    /* reverse trace */
                    let tmp = tr.data[i as usize];
                    tr.data[i as usize] = tr.data[(nt - 1 - i) as usize];
                    tr.data[(nt - 1 - i) as usize] = tmp;
                }
                input.copy_from_slice(&tr.data[..nt as usize]);
                bflowpass(npoleshi, f3dbhi, &input, &mut tr.data[..nt as usize]);
                for i in 0..nt / 2 {
                    /* flip trace back */
                    let tmp = tr.data[i as usize];
                    tr.data[i as usize] = tr.data[(nt - 1 - i) as usize];
                    tr.data[(nt - 1 - i) as usize] = tmp;
                }
            }
        }

        puttr(&mut tr);
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}
