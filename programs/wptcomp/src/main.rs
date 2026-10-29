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

//! WPTCOMP - Compression by Wavelet Packet Compression.
//! Port of comp/dct/main/wptcomp.c (Revision 1.6).
//! Credits: CWP: Tong Chen, Dec 1995.

use oxwave_core::prelude::*;
use oxwave_lib_comp_dct::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " WPTCOMP - Compression by Wavelet Packet Compression \t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   wptcomp < stdin n1= n2=   [optional parameter] > sdtout\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " n1=\t\t\tnumber of samples in the fast (first) dimension\t",
    " n2=\t\t\tnumber of samples in the slow (second) dimension",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " nfilter1=11\t\tnumber of filters in direction 1\t\t",
    " nfilter2=11\t\tnumber of filters in direction 2\t\t",
    " nstage1=\t\tfilter stages (automatically set based on n1)\t",
    " nstage2=\t\tfilter stages (automatically set based on n2)\t",
    " error=0.01\t\tacceptable error\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    run(&argv);
}

fn run(argv: &[String]) {
    initargs(argv);
    requestdoc(1, SDOC);

    /* get the parameters */
    let (mut n1, mut n2) = (0i32, 0i32);
    let (mut nfilter1, mut nfilter2) = (11i32, 11i32);
    let mut error = 0.01f32;
    if getparint("n1", &mut n1) == 0 {
        err!("Must specify n1\n");
    }
    if getparint("n2", &mut n2) == 0 {
        err!("Must specify n2\n");
    }
    if getparint("nfilter1", &mut nfilter1) == 0 {
        nfilter1 = 11;
    }
    if getparint("nfilter2", &mut nfilter2) == 0 {
        nfilter2 = 11;
    }
    if getparfloat("error", &mut error) == 0 {
        error = 0.01;
    }

    /* regular sizes */
    let mut npow1 = 0i32;
    let npad1;
    if n1 == 1 {
        npad1 = 1;
    } else {
        // C's arm64 shifts mask the count; nonpositive sizes never terminate.
        while n1.wrapping_sub(1).wrapping_shr(npow1 as u32) != 0 {
            npow1 = npow1.wrapping_add(1);
        }
        npad1 = 1i32.wrapping_shl(npow1 as u32);
    }
    let mut npow2 = 0i32;
    let npad2;
    if n2 == 1 {
        npad2 = 1;
    } else {
        while n2.wrapping_sub(1).wrapping_shr(npow2 as u32) != 0 {
            npow2 = npow2.wrapping_add(1);
        }
        npad2 = 1i32.wrapping_shl(npow2 as u32);
    }

    let (mut nstage1, mut nstage2) = (0i32, 0i32);
    if getparint("nstage1", &mut nstage1) == 0 {
        nstage1 = if npow1 > 4 {
            3
        } else if npow1 > 1 {
            2
        } else {
            npow1
        };
    }
    if getparint("nstage2", &mut nstage2) == 0 {
        nstage2 = if npow2 > 4 {
            3
        } else if npow2 > 1 {
            2
        } else {
            npow2
        };
    }
    // Negative stages cause undefined shifts in C. Reject them, as directed,
    // before allocation or input instead of casting them to huge usize values.
    if nstage1 < 0 {
        err!("nstage1 must be nonnegative");
    }
    if nstage2 < 0 {
        err!("nstage2 must be nonnegative");
    }
    checkpars();

    /* allocate space */
    let mut f = vec![vec![0.0f32; npad1 as usize]; npad2 as usize];
    let mut g = vec![vec![0.0f32; npad1 as usize]; npad2 as usize];
    let mut qx = vec![0i32; npad1.wrapping_mul(npad2) as usize];

    /* filter to use */
    let filter1 = waveGetfilter(nfilter1);
    let filter2 = waveGetfilter(nfilter2);

    /* input data */
    for i2 in 0..n2 as usize {
        fread(&mut f[i2][..n1 as usize], &stdin());
        /* pad with zeroes */
        for i1 in n1 as usize..npad1 as usize {
            f[i2][i1] = 0.0;
        }
    }
    /* pad with dead traces */
    for i2 in n2 as usize..npad2 as usize {
        for i1 in 0..npad1 as usize {
            f[i2][i1] = 0.0;
        }
    }

    /* perform the transform */
    // C dereferences an unsupported filter only when that dimension is active.
    if nstage1 <= npow1 && nstage2 <= npow2 {
        if (filter1.is_none() && nstage1 != 0) || (filter2.is_none() && nstage2 != 0) {
            unsafe extern "C" {
                fn signal(sig: i32, handler: usize) -> usize;
                fn raise(sig: i32) -> i32;
            }
            unsafe {
                // Restore SIG_DFL: Rust installs a SIGSEGV stack-overflow handler.
                signal(11, 0);
                raise(11);
            }
            exit(139);
        }
    }
    // A filter is unused at stage zero or when the stage check rejects the transform.
    let filter1 = filter1.unwrap_or_else(|| waveGetfilter(0).unwrap());
    let filter2 = filter2.unwrap_or_else(|| waveGetfilter(0).unwrap());
    wavePack_2(
        &mut f,
        &mut g,
        &filter1,
        &filter2,
        npow1 as usize,
        npow2 as usize,
        nstage1 as usize,
        nstage2 as usize,
        0,
    );

    /* quantization */
    let mut step = -1.0f32;
    let mut ave = 0.0f32;
    // The C allocation is contiguous: flatten rows in their original order.
    let g: Vec<f32> = g.into_iter().flatten().collect();
    uniQuant(&g, error, &mut ave, &mut step, &mut qx);

    /* prefix encoding */
    let mut ibuff = pEncode(&qx);
    /* allocate out buffer */
    let mut obuff = buffAlloc1(ibuff.code.len());
    /* rewind the in buffer */
    ibuff.pos = 0;
    /* run-length coding */
    codeSilence(&mut ibuff, &mut obuff);
    /* Huffman coding */
    let pos = obuff.pos;
    buffRealloc1(&mut obuff, pos);
    obuff.pos = 0;
    ibuff.pos = 0;
    let nsize = huffCompress(&mut obuff, &mut ibuff);

    fprintf!(stderr(), "size after compression = %d bytes\n", nsize);
    fprintf!(
        stderr(),
        "compression ratio = %f\n",
        (n1 as f32 * n2 as f32 * std::mem::size_of::<i32>() as f32 / nsize as f32) as f64
    );

    /* Write out compressed data */
    fwrite(&[nsize], &stdout());
    fwrite(&[n1], &stdout());
    fwrite(&[n2], &stdout());
    fwrite(&[nfilter1], &stdout());
    fwrite(&[nfilter2], &stdout());
    fwrite(&[nstage1], &stdout());
    fwrite(&[nstage2], &stdout());
    fwrite(&[ave], &stdout());
    fwrite(&[step], &stdout());
    fwrite(&ibuff.code[..ibuff.pos], &stdout());
    exit(CWP_Exit());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::process::{Command, Stdio};

    #[test]
    fn reject_negative_stages() {
        const CHILD_ARGS: &str = "WPTCOMP_NEGATIVE_STAGE_TEST_ARGS";
        const MARKER: &str = "BEGIN_WPTCOMP\n";
        if let Ok(parameters) = std::env::var(CHILD_ARGS) {
            let mut argv: Vec<String> = ["wptcomp", "n1=16", "n2=8", "nfilter1=15", "nfilter2=15"]
                .into_iter()
                .map(String::from)
                .collect();
            argv.extend(parameters.split_whitespace().map(String::from));
            // Separate the test harness preamble from the program's output.
            print!("{MARKER}");
            std::io::stdout().flush().unwrap();
            run(&argv);
            unreachable!("run must exit");
        }

        // Exercise each dimension, signed limits, and validation order. Invalid
        // filters would fault if the program reached the transform first.
        for (parameters, rejected) in [
            ("nstage1=-1", "nstage1"),
            ("nstage2=-1", "nstage2"),
            ("nstage1=-2147483648", "nstage1"),
            ("nstage2=-2147483648", "nstage2"),
            ("nstage1=-1 nstage2=-1", "nstage1"),
            ("nstage1=5 nstage2=-1", "nstage2"),
            ("nstage1=-1 nstage2=4", "nstage1"),
        ] {
            let output = Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "tests::reject_negative_stages", "--nocapture"])
                .env(CHILD_ARGS, parameters)
                .stdin(Stdio::null())
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(1), "{parameters}");
            let marker_end = output
                .stdout
                .windows(MARKER.len())
                .position(|window| window == MARKER.as_bytes())
                .expect("child must enter run")
                + MARKER.len();
            assert!(
                output.stdout[marker_end..].is_empty(),
                "{parameters}: compressed output emitted"
            );
            assert_eq!(
                output.stderr,
                format!("\nwptcomp: {rejected} must be nonnegative\n").as_bytes(),
                "{parameters}"
            );
        }
    }
}
