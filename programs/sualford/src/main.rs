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

//! SUALFORD - trace by trace Alford Rotation of shear wave data volumes.
//! Port of su/main/multicomponent/sualford.c (Revision 1.6).
//! Original copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: CWP: Rodrigo Felicio Fuck. Code translated and adapted from the
//! original version in Fortran by Ted Schuck (1993).
//! Schuck, E. L., 1993, Multicomponent, three dimensional seismic
//! characterization of a fractured coalbed methane reservoir,
//! Cedar Hill Field, San Juan County, New Mexico, Ph.D. Thesis,
//! Colorado School of Mines.

use oxwave_core::prelude::*;
use oxwave_correlation::*;
use oxwave_cwp::cwp::{nint, PI};
use oxwave_sinc_interp::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUALFORD - trace by trace Alford Rotation of shear wave data volumes  ",
    " \t\t\t\t\t\t\t\t\t",
    " sualford inS11=file1 inS22=file2 inS12=file3 inS21=file4\t\t",
    " outS11=file5 outS22=file6 outS12=file7 outS21=file8 [optional         ",
    " parameters]                                                           ",
    "                                                                       ",
    " Required Parameters:                                                  ",
    " inS11=\tinput data volume for the 11 component\t\t\t",
    " inS12=\tinput data volume for the 12 component\t\t\t",
    " inS21=\tinput data volume for the 21 component\t\t\t",
    " inS22=\tinput data volume for the 22 component\t\t\t",
    " outS11=\toutput data volume for the 11 component\t\t\t",
    " outS12=\toutput data volume for the 11 component\t\t\t",
    " outS21=\toutput data volume for the 11 component\t\t\t",
    " outS22=\toutput data volume for the 11 component\t\t\t",
    "                                                                       ",
    " Optional parameters:                                                  ",
    " angle_inc=               sets the increment to the angle by which\t",
    "                         the data sets are rotated. The minimum is     ",
    "                         set to be 1 degree and default is 5.          ",
    " Az_key=                  to set the header field storing the azimuths\t",
    "                         for the fast shear wave on the output volumes ",
    " Q_key=                   to set the header field storing the quality\t",
    "                         factors of performed optimum rotations        ",
    " lag_key=                 to set the header field storing the lag in\t",
    "                         miliseconds the fast and slow shear components",
    " xcorr_key=               to set the header field containing the maxi-\t",
    "                         mum normalized cross-correlation between the\t",
    "                         and slow shear waves.                         ",
    " taper=\t\t  2*taper+1 is the length of the sample overlap ",
    "\t\t\t  between the unrotated data with the rotated   ",
    "\t\t\t  data on the traces. The boundary between them ",
    " \t\t\t  is defined by time windowning.                ",
    "\t\t\t\ttaper = -1, for no-overlap\t\t",
    "\t\t\t\ttaper = 0, for overlap of one sample\t",
    "\t\t\t\ttaper =>1, for use of cosine scale to   ",
    "\t\t\t\t\t   to interpolate between the \t",
    "\t\t\t\t\t   unrotated and rotated parts\t",
    "\t\t\t\t\t   of the traces\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " taperwin=               another taper used to taper the data within   ",
    "\t\t\t  the window of analysis to diminish the effect ",
    "                         of data near the window edges.In this way one ",
    "                         can focus on a given reflector. Also given in ",
    "                         number of samples                             ",
    "\t\t\t\t\t\t\t\t\t",
    " maxlag=\t\t  maximum limit in ms for the lag between fast  ",
    " \t\t\t  and slow shear waves. If this threshold is \t",
    "\t\t\t  attained or surpassed, the quality factor for\t",
    "\t\t\t  the rotation is zeroed as well as all the     ",
    "\t\t\t  parameters found for that certain rotation \t",
    "\t\t\t \t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " ntraces=\t\t  number of traces to be used per computation   ",
    "\t\t\t  ntraces=3 will use three adjacent traces to   ",
    "\t\t          compute the angle of rotation                                                                                        ",
    " Notes:                                                                ",
    "                                                                       ",
    " The Alford Rotation is a method to rotate the four components         ",
    " of a shear wave survey into its natural coordinate system, where      ",
    " the fast and slow shear correspond to the inline to inline shear (S11)",
    " and xline to xline (S22) volumes, respectively.                       ",
    "                                                                       ",
    " This Alford Rotation code tries to maximize the energy in the         ",
    " diagonal volumes, i.e., S11 and S22, while minimizing the energy      ",
    " in the off-diagonals, i.e., in volumes S12 and S21, in a trace by     ",
    " trace manner. It then returns the new rotated volumes, saving the     ",
    " the quality factor for the rotation and azimuth angle of the fast     ",
    " shear wave direction for each trace headers of the new rotated S11    ",
    " volume.                                                               ",
    "                                                                       ",
    " The fields in the header containing the Azimuth and Quality factor    ",
    " and the sample lag between fast and slow shear are otrav, grnolf and  ",
    " grnofr, respectively, by default. The values are multiplied by ten in ",
    " the case of the angles and by a thousand for quality factors. To      ",
    " change this defaults use the optional parameters Az_key, Q_key and    ",
    " lag_key                                                             \t",
    "                                                                       ",
    "\t\t\t   \t\t\t\t\t\t",
    " modified header fields:                                               ",
    " the ones specified by Az_key, Q_key, lag_key and xcorr_key. By default",
    " these are otrav, grnlof, tstat and grnors, respectively.            \t",
];

// Component order throughout is S11, S12, S21, S22.
type Components = [Vec<f32>; 4];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut ntraces = 5i32;
    getparint("ntraces", &mut ntraces);
    let mut names: [String; 4] = std::array::from_fn(|_| String::new());
    if getparstring("inS11", &mut names[0]) == 0 {
        err!("no file specified for S11 volume. Aborting...");
    }
    if getparstring("inS12", &mut names[1]) == 0 {
        err!("no file specified for S12 volume. Aborting...");
    }
    if getparstring("inS21", &mut names[2]) == 0 {
        err!("no file specified for S21 volume. Aborting...");
    }
    if getparstring("inS22", &mut names[3]) == 0 {
        err!("no file specified for S22 volume. Aborting...");
    }
    let mut outnames: [String; 4] = std::array::from_fn(|_| String::new());
    for (k, component) in ["S11", "S12", "S21", "S22"].iter().enumerate() {
        if getparstring(&format!("out{}", component), &mut outnames[k]) == 0 {
            outnames[k] = format!("{}_out", component);
        }
    }
    let mut input: [Trace; 4] = std::array::from_fn(|_| Trace::new());
    let mut ns = [0usize; 4];
    let mut dt = [0.0f32; 4];
    let files: [CFile; 4] =
        std::array::from_fn(|k| open_vol(&names[k], &mut ns[k], &mut dt[k], &mut input[k]));
    if ntraces < 0 {
        err!("ealloc.c: malloc failed (Cannot allocate memory)");
    }
    // ealloc2float stores each panel contiguously.
    // C copies the first trace before validating lengths, even for a zero
    // mixing count. Reserve that row and copy space so these writes do not
    // corrupt memory before the existing length checks run.
    let mut s: [Vec<Vec<f32>>; 4] =
        std::array::from_fn(|k| vec![vec![0.0; ns[k].max(ns[0])]; (ntraces as usize).max(1)]);
    for k in 0..4 {
        s[k][0][..ns[0]].copy_from_slice(&input[k].data[..ns[0]]);
    }

    let mut keys: [String; 5] = std::array::from_fn(|_| String::new());
    for (k, (name, default)) in [
        ("Az_key", "otrav"),
        ("Q_key", "grnlof"),
        ("lag_key", "tstat"),
        ("xcorr_key", "grnors"),
        ("INPUTQ_key", "gaps"),
    ]
    .iter()
    .enumerate()
    {
        if getparstring(name, &mut keys[k]) == 0 {
            keys[k] = default.to_string();
        }
    }
    // Preserve the C order of header lookups.
    let az_type = hdtype(&keys[0]);
    let q_type = hdtype(&keys[1]);
    let inputq_type = hdtype(&keys[4]);
    let lag_type = hdtype(&keys[2]);
    let xcorr_type = hdtype(&keys[3]);
    let types = [az_type, q_type, lag_type, xcorr_type, inputq_type];
    let index_az = getindex(&keys[0]);
    let index_q = getindex(&keys[1]);
    let index_inputq = getindex(&keys[4]);
    let index_lag = getindex(&keys[2]);
    let index_xcorr = getindex(&keys[3]);
    let indices = [index_az, index_q, index_lag, index_xcorr, index_inputq];

    let mut taper = 3i32;
    getparint("taper", &mut taper);
    let mut taper2 = 10i32;
    getparint("taperwin", &mut taper2);
    let mut angle_inc = 5.0f32;
    getparfloat("angle_inc", &mut angle_inc);
    let mut total = 0i32;
    getparint("Total", &mut total);
    let mut maxlag = 40.0f32;
    getparfloat("maxlag", &mut maxlag);
    let mut tmin = 0.0f32;
    let itmin;
    if getparfloat("tmin", &mut tmin) != 0 {
        itmin = nint((tmin as f64 - input[0].delrt as f64 / 1000.0) / dt[0] as f64);
    } else {
        itmin = 0;
    }
    let mut tmax = 0.0f32;
    let itmax;
    if getparfloat("tmax", &mut tmax) != 0 {
        itmax = nint((tmax as f64 - input[0].delrt as f64 / 1000.0) / dt[0] as f64);
    } else {
        itmax = ns[0] as i32 - 1;
    }
    let nt = itmax - itmin + 1;
    checkpars();
    if itmin < 0 {
        err!("itmin = %d should be positive", itmin);
    }
    if nt > SU_NFLTS as i32 {
        err!(" nt = %d exceeds SU_NFLTS=%d", nt, SU_NFLTS);
    }
    if itmin > itmax {
        err!("itmin=%d, itmax = %d conflict", itmin, itmax);
    }
    let outputs: [CFile; 4] = std::array::from_fn(|k| check_output(&outnames[k]));
    if ns[0] != ns[3] || ns[0] != ns[1] || ns[0] != ns[2] {
        err!("trace length for S11 is not equal to other traces \n");
    }
    if ns[3] != ns[1] || ns[3] != ns[2] {
        err!("trace length for S22 is not equal to other traces\n");
    }
    if ns[1] != ns[2] {
        err!("trace length for off diagonals is not the same\n");
    }
    // C only checks the window length against SU_NFLTS, then indexes panels
    // allocated for ns samples. Reject its undefined out-of-bounds reads.
    if itmax >= ns[0] as i32 {
        err!(
            "itmax = %d exceeds last trace sample = %d",
            itmax,
            ns[0] as i32 - 1
        );
    }

    let mix = vec![1.0f32 / ntraces as f32; ntraces as usize];
    let mut mixed: Components = std::array::from_fn(|k| vec![0.0; ns[k]]);
    let mut output: [Trace; 4] = std::array::from_fn(|_| Trace::new());
    let mut counter = 0i32;
    loop {
        if counter != 0 {
            s[3][0][..ns[3]].copy_from_slice(&input[3].data[..ns[3]]);
        }
        for k in 0..4 {
            mixed[k].fill(0.0);
        }
        for i in itmin..=itmax {
            for j in 0..ntraces as usize {
                for k in 0..4 {
                    mixed[k][i as usize] += s[k][j][i as usize] * mix[j];
                }
            }
        }
        let (mut angle, mut quality, inputquality, mut lag, xcorr) = find_angle(
            &mut mixed,
            angle_inc,
            ns[0],
            dt[0],
            nt as usize,
            itmin,
            itmax,
            taper2,
        );
        if lag * 1000.0 > maxlag || quality <= 0.0 {
            quality = 0.0;
            lag = 0.0;
            angle = 0.0;
        }
        if nt as usize != ns[0] {
            for k in 0..4 {
                mixed[k].copy_from_slice(&s[k][0]);
            }
            final_rot(&mut mixed, angle, ns[0], itmin, taper);
        }
        for k in 0..4 {
            output[k].data[..ns[k]].copy_from_slice(&mixed[k]);
            output[k].hdr = input[k].hdr.clone();
        }
        for (field, number) in [
            (0, angle * 10.0),
            (1, quality * 1000.0),
            (4, inputquality * 1000.0),
        ] {
            let value = type_conv(nint(number as f64), types[field]);
            for tr in &mut output {
                puthval(tr, indices[field], &value);
            }
        }
        let value = type_conv(nint((lag * 1000.0) as f64), lag_type);
        puthval(&mut output[0], index_lag, &value);
        puthval(&mut output[3], index_lag, &value);
        let value = type_conv(nint((lag * 500.0) as f64), lag_type);
        puthval(&mut output[1], index_lag, &value);
        puthval(&mut output[2], index_lag, &value);
        let value = type_conv(nint((xcorr * 100.0) as f64), xcorr_type);
        for tr in &mut output {
            puthval(tr, index_xcorr, &value);
        }
        for k in 0..4 {
            fputtr(&outputs[k], &mut output[k]);
        }
        for i in (1..ntraces as usize).rev() {
            for j in 0..ns[0] {
                for k in 0..4 {
                    s[k][i][j] = s[k][i - 1][j];
                }
            }
        }
        if fgettr(&files[0], &mut input[0]) == 0 {
            break;
        }
        if fgettr(&files[1], &mut input[1]) == 0 {
            break;
        }
        if fgettr(&files[2], &mut input[2]) == 0 {
            break;
        }
        for k in 0..3 {
            s[k][0].copy_from_slice(&input[k].data[..ns[k]]);
        }
        counter += 1;
        if total != 0 {
            let f_counter = counter as f32;
            printf!("%2.f %% done\r", 100.0f32 * (f_counter / total as f32));
        }
        if fgettr(&files[3], &mut input[3]) == 0 {
            break;
        }
    }
    exit(CWP_Exit());
}

fn rms_amp(tr: &[f32], nt: usize) -> f32 {
    let mut amp = 0.0f32;
    let mut counter = 0i32;
    for i in 0..nt {
        amp += tr[i] * tr[i];
        counter += 1;
    }
    ((amp / counter as f32) as f64).sqrt() as f32
}

fn quality_factor(s: &Components, nt: usize, mode: i32) -> f32 {
    let noise = 1.0e-30f32;
    let amp_s11 = rms_amp(&s[0], nt) + noise;
    let amp_s12 = rms_amp(&s[1], nt) + noise;
    let amp_s21 = rms_amp(&s[2], nt) + noise;
    let amp_s22 = rms_amp(&s[3], nt) + noise;
    let (max, min) = if amp_s11 >= amp_s22 {
        (amp_s11, amp_s22)
    } else {
        (amp_s22, amp_s11)
    };
    let num = amp_s12 + amp_s21;
    let den = amp_s11 + amp_s22;
    if mode == 0 {
        (1.0 - (num / den) as f64) as f32
    } else {
        ((min / max) as f64 * (1.0 - (num / den) as f64)) as f32
    }
}

fn rot4tr(s: &Components, out: &mut Components, theta: f32, itmin: i32, itmax: i32, taper: i32) {
    let rad = (theta as f64 * (PI / 180.0)) as f32;
    let cos2 = (rad as f64).cos().powf(2.0) as f32;
    let sin2 = (rad as f64).sin().powf(2.0) as f32;
    let sincos = ((rad as f64).cos() * (rad as f64).sin()) as f32;
    let nt = itmax - itmin + 1;
    for i in itmin..=itmax {
        let j = (i - itmin) as usize;
        let i = i as usize;
        let temp = sincos * (s[1][i] + s[2][i]);
        out[0][j] = s[0][i] * cos2 + s[3][i] * sin2 + temp;
        out[3][j] = s[0][i] * sin2 + s[3][i] * cos2 - temp;
        let temp = sincos * (s[3][i] - s[0][i]);
        out[1][j] = s[1][i] * cos2 - s[2][i] * sin2 + temp;
        out[2][j] = -s[1][i] * sin2 + s[2][i] * cos2 + temp;
    }
    if taper > 1 {
        // For windows shorter than taper, the C also touches negative sample
        // indices and samples beyond nt. Put those accesses in guard storage;
        // retain the exact sequence of multiplications on actual samples.
        if taper as usize > nt as usize {
            let offset = taper as usize;
            let mut padded: Components = std::array::from_fn(|k| {
                let mut v = vec![0.0; nt as usize + 2 * offset + 1];
                v[offset..offset + nt as usize].copy_from_slice(&out[k][..nt as usize]);
                v
            });
            for i in 0..taper {
                let f = i as f32 / taper as f32;
                let env = (0.5 * (1.0 - (PI * f as f64).cos())) as f32;
                for k in [0, 3, 1, 2] {
                    padded[k][offset + i as usize] *= env;
                }
                for k in [0, 3, 1, 2] {
                    padded[k][(taper + nt - i) as usize] *= env;
                }
            }
            for k in 0..4 {
                out[k][..nt as usize].copy_from_slice(&padded[k][offset..offset + nt as usize]);
            }
            return;
        }
        for i in 0..taper {
            let f = i as f32 / taper as f32;
            let env = (0.5 * (1.0 - (PI * f as f64).cos())) as f32;
            for k in [0, 3, 1, 2] {
                out[k][i as usize] *= env;
            }
            // C addresses sample nt on its first taper iteration. This guard
            // sample preserves the write without corrupting the allocation.
            for k in [0, 3, 1, 2] {
                out[k][(nt - i) as usize] *= env;
            }
        }
    }
}

fn final_rot(s: &mut Components, theta: f32, nt: usize, itmin: i32, taper: i32) {
    // The C overlap can extend before sample zero or beyond the trace.
    // These guard cells keep those writes separate from the actual samples.
    let n = if taper == -1 { 1 } else { 0 };
    let prefix = (taper + n - itmin).max(0);
    let end = (nt as i32).max(itmin + taper + 1).max(itmin - taper - n);
    let source: Components = std::array::from_fn(|k| {
        let mut v = vec![0.0; (prefix + end) as usize];
        v[prefix as usize..prefix as usize + nt].copy_from_slice(&s[k][..nt]);
        v
    });
    let mut out: Components = std::array::from_fn(|_| vec![0.0; (prefix + end) as usize]);
    let rad = (theta as f64 * (PI / 180.0)) as f32;
    let cos2 = (rad as f64).cos().powf(2.0) as f32;
    let sin2 = (rad as f64).sin().powf(2.0) as f32;
    let sincos = ((rad as f64).cos() * (rad as f64).sin()) as f32;
    for i in 0..=itmin - taper - 1 - n {
        for k in 0..4 {
            out[k][(prefix + i) as usize] = source[k][(prefix + i) as usize];
        }
    }
    for i in itmin - taper - n..nt as i32 {
        let i = (prefix + i) as usize;
        let temp = sincos * (source[1][i] + source[2][i]);
        out[0][i] = source[0][i] * cos2 + source[3][i] * sin2 + temp;
        out[3][i] = source[0][i] * sin2 + source[3][i] * cos2 - temp;
        let temp = sincos * (source[3][i] - source[0][i]);
        out[1][i] = source[1][i] * cos2 - source[2][i] * sin2 + temp;
        out[2][i] = -source[1][i] * sin2 + source[2][i] * cos2 + temp;
    }
    let overlap = 2 * taper + 1;
    if overlap == 1 {
        for k in [0, 1, 2, 3] {
            out[k][(prefix + itmin) as usize] = (0.5
                * (out[k][(prefix + itmin) as usize] + source[k][(prefix + itmin) as usize]) as f64)
                as f32;
        }
    }
    if overlap > 1 {
        for i in 0..overlap {
            // The integer division in the original cosine scale is deliberate.
            let c1 = ((PI / 2.0) * (i / (overlap - 1)) as f64).cos() as f32;
            let c2 = 1.0f32 - c1;
            let j = (prefix + itmin - taper + i) as usize;
            for k in 0..4 {
                out[k][j] = c1 * source[k][j] + c2 * out[k][j];
            }
        }
    }
    for k in 0..4 {
        s[k][..nt].copy_from_slice(&out[k][prefix as usize..prefix as usize + nt]);
    }
}

fn find_angle(
    s: &mut Components,
    mut ang_inc: f32,
    ns: usize,
    dt: f32,
    nt: usize,
    itmin: i32,
    itmax: i32,
    taper: i32,
) -> (f32, f32, f32, f32, f32) {
    let mut out: Components = std::array::from_fn(|_| vec![0.0; nt + 1]);
    if ang_inc < 1.0 {
        ang_inc = 5.0;
    }
    let ns_angle = (1.0 + ((90.0f32 / ang_inc) as f64).floor()) as usize;
    // The periodic peak interpolation reads two neighboring Q entries.
    // With fewer than two search angles those reads are undefined in C.
    if ns_angle < 2 {
        err!("angle_inc = %g gives fewer than two search angles", ang_inc);
    }
    let mut angles = vec![0.0f32; ns_angle];
    let mut q = vec![0.0f32; ns_angle];
    for i in 1..ns_angle {
        angles[i] = angles[i - 1] + ang_inc;
    }
    for i in 0..ns_angle {
        rot4tr(s, &mut out, angles[i], itmin, itmax, taper);
        q[i] = quality_factor(&out, nt, 0);
    }
    let inputquality = q[0];
    let mut qmax = -1.0e30f32;
    let mut imax = 0usize;
    for i in 0..ns_angle {
        if q[i] > qmax {
            qmax = q[i];
            imax = i;
        }
    }
    let y = if imax != 0 && imax != ns_angle - 1 {
        [q[imax - 1], q[imax], q[imax + 1]]
    } else if imax == 0 {
        [q[ns_angle - 2], q[imax], q[imax + 1]]
    } else {
        [q[imax - 1], q[imax], q[1]]
    };
    let angle_frac = quad_fit(y, &mut qmax);
    angles[imax] += ang_inc * angle_frac;
    rot4tr(s, &mut out, angles[imax], itmin, itmax, taper);
    let ntres = 4 * nt;
    let dtres = dt / 4.0f32;
    let mut s1res = vec![0.0f32; ntres];
    let mut s2res = vec![0.0f32; ntres];
    interp(&out[0][..nt], &out[3][..nt], &mut s1res, &mut s2res);
    let ntemp = 2 * ntres - 1;
    let mut temp = vec![0.0f32; ntemp];
    nxcor(&s1res, &s2res, -(ntres as i32) + 1, &mut temp);
    let (peak_index, xcorr) = find_peak(&temp);
    let zero_lag = ntres as i32 - 1;
    let mut lag = dtres * (peak_index - zero_lag) as f32;
    if peak_index < zero_lag {
        angles[imax] += 90.0;
        if angles[imax] >= 180.0 {
            angles[imax] -= 180.0;
        }
        lag *= -1.0;
    }
    let amps1 = rms_amp(&out[0], nt);
    let amps2 = rms_amp(&out[3], nt);
    if peak_index == zero_lag && amps2 > amps1 {
        angles[imax] += 90.0;
        if angles[imax] >= 180.0 {
            angles[imax] -= 180.0;
        }
    }
    let quality = quality_factor(&out, nt, 0);
    if nt == ns {
        az_swap(&mut out, nt);
        for k in [0, 3, 1, 2] {
            s[k][..nt].copy_from_slice(&out[k][..nt]);
        }
    }
    (angles[imax], quality, inputquality, lag, xcorr)
}

#[allow(dead_code)]
fn total_quality(s11: &[f32], s22: &[f32], qmin: f32, qmax: f32, nt: usize) -> f32 {
    let noise = 1.0e-30f32;
    let amp_s11 = rms_amp(s11, nt) + noise;
    let amp_s22 = rms_amp(s22, nt) + noise;
    let (max, min) = if amp_s11 >= amp_s22 {
        (amp_s11, amp_s22)
    } else {
        (amp_s22, amp_s11)
    };
    ((min / max) as f64 * (1.0 - (qmin / qmax) as f64)) as f32
}

fn az_swap(s: &mut Components, nt: usize) {
    s.swap(0, 3);
    for i in 0..nt {
        s[1][i] *= -1.0;
        s[2][i] *= -1.0;
    }
    s.swap(1, 2);
}

fn find_peak(y: &[f32]) -> (i32, f32) {
    let mut max = -1.0e30f32;
    let mut index = 0i32;
    for i in 0..y.len() {
        if y[i] > max {
            max = y[i];
            index = i as i32;
        }
    }
    (index, max)
}

fn quad_fit(y: [f32; 3], max: &mut f32) -> f32 {
    let c = y[1];
    let b = ((y[2] - y[0]) as f64 * 0.5) as f32;
    let a = y[0] + b - c;
    let x = (-b as f64 / (2.0 * a as f64)) as f32;
    *max = c + x * (b + x * a);
    x
}

// TypeConv, adapted in the C source from John Stockwell's suazimuth.
fn type_conv(dval_out: i32, type_in: &str) -> Value {
    let mut val = Value::default();
    match type_in.as_bytes().first().copied().unwrap_or(0) {
        b's' => err!("can't change char header word"),
        b'h' => val.set_h(dval_out as i16),
        b'u' => val.set_u(dval_out as u16),
        b'l' => val.set_l(dval_out as i64),
        b'v' => val.set_v(dval_out as u64),
        b'i' => val.set_i(dval_out),
        b'p' => val.set_p(dval_out as u32),
        b'f' => val.set_f(dval_out as f32),
        b'd' => val.set_d(dval_out as f64),
        _ => err!("unknown type %s", type_in),
    }
    val
}

fn nxcor(x: &[f32], y: &[f32], z1: i32, z: &mut [f32]) {
    let mut sum_x = 0.0f32;
    let mut sum_y = 0.0f32;
    xcor(0, x, 0, y, z1, z);
    if y.len() == x.len() {
        for i in 0..x.len() {
            sum_x += x[i] * x[i];
            sum_y += y[i] * y[i];
        }
    } else {
        for i in 0..x.len() {
            sum_x += x[i] * x[i];
        }
        for i in 0..y.len() {
            sum_y += y[i] * y[i];
        }
    }
    sum_x = (sum_x as f64).sqrt() as f32;
    sum_y = (sum_y as f64).sqrt() as f32;
    let n = sum_x * sum_y;
    for i in 0..z.len() {
        z[i] /= n;
    }
}

fn interp(s1: &[f32], s2: &[f32], s1res: &mut [f32], s2res: &mut [f32]) {
    let mut t = vec![0.0f32; s1res.len()];
    let mut tvalue = 0.0f32;
    for itime in 0..t.len() {
        t[itime] = tvalue;
        tvalue = (tvalue as f64 + 0.25) as f32;
    }
    ints8r(1.0, 0.0, s1, 0.0, 0.0, &t, s1res);
    ints8r(1.0, 0.0, s2, 0.0, 0.0, &t, s2res);
}

fn open_vol(name: &str, ns: &mut usize, dt: &mut f32, tr: &mut Trace) -> CFile {
    let fp = match fopen(name, "r") {
        Some(fp) => fp,
        None => err!("cannot open volume %s\n", name),
    };
    let mut tr2 = Trace::new();
    if fgettr(&fp, &mut tr2) == 0 {
        err!("cannot read first trace in %s volume\n", name);
    }
    *dt = (tr2.dt as f64 / 1000000.0) as f32;
    *ns = tr2.ns as usize;
    // The C tests the non-null dt pointer, rather than the header value.
    *tr = tr2;
    fp
}

fn check_output(filename: &str) -> CFile {
    match fopen(filename, "w") {
        Some(fp) => fp,
        None => err!("cannot open output volume %s file \n", filename),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::CommandExt;
    use std::process::{Command, Stdio};

    // Run after cargo build --release -p sualford with cargo test --release.
    // These inputs invoke undefined reads in C, so golden comparison cannot
    // establish the required clean failure. Exercise the Rust CLI directly.
    #[test]
    fn invalid_analysis_bounds_and_single_angle_exit_cleanly() {
        let exe = std::env::current_exe().unwrap();
        let target = exe.parent().unwrap().parent().unwrap();
        let binary = target.join("sualford");
        assert!(binary.is_file(), "build sualford before running this test");
        let work = target.join(format!("sualford-boundaries-{}", std::process::id()));
        std::fs::create_dir(&work).unwrap();
        let input = work.join("input.su");
        let fp = efopen(input.to_str().unwrap(), "w");
        let mut tr = Trace::new();
        tr.ns = 65;
        tr.dt = 4000;
        tr.trid = 1;
        fputtr(&fp, &mut tr);
        efclose(&fp);

        for (options, message) in [
            (
                vec!["tmax=0.26"],
                "itmax = 65 exceeds last trace sample = 64",
            ),
            (
                vec!["tmin=0.004", "tmax=0.26"],
                "itmax = 65 exceeds last trace sample = 64",
            ),
            (
                vec!["tmin=0.26", "tmax=0.26"],
                "itmax = 65 exceeds last trace sample = 64",
            ),
            (
                vec!["angle_inc=90.00001"],
                "angle_inc = 90 gives fewer than two search angles",
            ),
            (
                vec!["angle_inc=91"],
                "angle_inc = 91 gives fewer than two search angles",
            ),
            (
                vec!["angle_inc=180"],
                "angle_inc = 180 gives fewer than two search angles",
            ),
            (vec!["angle_inc=inf"], "atopkge.c: eatof: overflow"),
            (
                vec!["angle_inc=nan"],
                "angle_inc = nan gives fewer than two search angles",
            ),
        ] {
            let result = Command::new(&binary)
                .arg0("sualford")
                .args([
                    "inS11=input.su",
                    "inS12=input.su",
                    "inS21=input.su",
                    "inS22=input.su",
                ])
                .args(&options)
                .current_dir(&work)
                .stdin(Stdio::null())
                .output()
                .unwrap();
            assert_eq!(result.status.code(), Some(1), "{options:?}: {result:?}");
            assert!(result.stdout.is_empty(), "{options:?}");
            assert_eq!(
                String::from_utf8(result.stderr).unwrap(),
                format!("\nsualford: {message}\n"),
                "{options:?}",
            );
        }
        std::fs::remove_dir_all(work).unwrap();
    }
}
