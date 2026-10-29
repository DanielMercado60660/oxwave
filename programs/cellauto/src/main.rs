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

//! CELLAUTO - Two-dimensional CELLular AUTOmata.
//! Port of par/main/cellular_automata/cellauto.c (Revision 1.4).
//! Credits (from the C source): UHouston: Chris Liner.

use oxwave_core::err::strerror;
use oxwave_core::prelude::*;
use oxwave_cwp::cwp::isodd;
use oxwave_f_franuni::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " CELLAUTO - Two-dimensional CELLular AUTOmata\t\t\t  \t",
    "\t\t\t\t\t\t\t\t\t",
    "   cellauto > stdout [optional params]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " n1=500\toutput dimensions of image (n1 x n1 pixels)\t \t",
    " rule=30\tCA rule (Wolfram classification)\t\t\t",
    " \t\tOthers: 54,60,62,90,94,102,110,122,126\t\t\t",
    "                       150,158,182,188,190,220,222,225,226,250\t\t",
    " fill=0\tDon't fill image (=1 fill image)\t\t\t",
    " f0=330\tfill zero values with f0\t\t\t\t",
    " f1=3000\tfill non-zero values with f1\t\t\t\t",
    " ic=1\t\tinitial condition for centered unit value at t=0\t",
    "               = 2 for multiple random units\t\t\t\t",
    " nc=20\t\tnumber of random units (if ic=2)\t\t\t",
    " tc=1\t\trandom initial units at t=0 (if ic=2)\t\t\t",
    "               = 2 for initial units at random (t,x)\t\t\t",
    " verbose=0\tsilent operation\t\t\t\t\t",
    "               = 1 echos 'porosity' of the CA in bottom half of image\t",
    " seed=from_clock    \trandom number seed (integer)            \t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " This program generates a select set of Wolframs fundamental cellular\t",
    " automata. This may be useful for constructing rough, vuggy wavespeed\t",
    " profiles. The numbering scheme follows Stephen Wolfram's.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Example: \t\t\t\t\t\t\t\t",
    "  cellauto rule=110 ic=2 nc=100 fill=1 f1=3000 | ximage n1=500 nx=500 &",
    "\t\t\t\t\t\t\t\t\t",
    " Here we simulate a complex near surface with air-filled \t\t",
    " vugs in hard country rock, with smoothing applied via smooth2 \t",
    "\t\t\t\t\t\t\t\t\t",
    "  cellauto rule=110 ic=2 nc=100 fill=1 f1=3000 n1=500 |\t\t",
    "  smooth2 n1=500 n2=500 r1=5 r2=5 > vfile.bin\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    let mut n1 = 0i32;
    let mut code = 0i32;
    let mut rule = 0i32;
    let mut ic = 0i32;
    let mut nc = 0i32;
    let mut tc = 0i32;
    let mut fill = 0i32;
    let mut f0 = 0.0f32;
    let mut f1 = 0.0f32;
    let mut seed = 0u32;

    /* hook up getpar */
    initargs(&argv);
    requestdoc(0, SDOC);

    /* get parameters */
    if getparint("n1", &mut n1) == 0 {
        n1 = 500;
    }
    let nz = if isodd(n1) { n1 } else { n1.wrapping_add(1) };
    /* set output slow dimension */
    let n2 = n1;
    let nx = 2i32.wrapping_mul(nz).wrapping_add(1);

    if getparint("rule", &mut rule) == 0 {
        rule = 30;
    }
    if getparint("ic", &mut ic) == 0 {
        ic = 1;
    }
    if getparint("nc", &mut nc) == 0 {
        nc = 20;
    }
    if getparint("tc", &mut tc) == 0 {
        tc = 1;
    }
    if getparint("fill", &mut fill) == 0 {
        fill = 0;
    }
    if getparfloat("f0", &mut f0) == 0 {
        f0 = 330.0;
    }
    if getparfloat("f1", &mut f1) == 0 {
        f1 = 3000.0;
    }

    /* allocate workspace */
    // Negative dimensions become enormous size_t allocations in C and fail
    // with ENOMEM. Avoid a Rust capacity-overflow panic for the same input.
    if nz < 0 || nx < 0 {
        err!("%s: malloc failed (%s)", "ealloc.c", strerror(12));
    }
    let mut data = vec![0.0f32; nz as usize];
    // C leaves this allocation uninitialized. Fresh allocations in the SU
    // executable are zero-filled in the reproducible cases; make that state
    // explicit rather than reading uninitialized Rust memory.
    let mut ca = vec![vec![0.0f32; nz as usize]; nx as usize];

    /* Set seed */
    if getparuint("seed", &mut seed) == 0 {
        seed = match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
            Ok(t) => t.as_secs() as u32,
            Err(e) => (0u64.wrapping_sub(e.duration().as_secs())) as u32,
        };
        if seed == u32::MAX {
            err!("time() failed to set seed");
        }
    }
    sranuni(seed as i32);

    /* set CA seed */
    if ic == 1 {
        ca[((nx - 1) / 2) as usize][0] = 1.0;
    }
    if ic == 2 {
        for _i in 0..nc {
            let ix = ((nx / 4) as f32 + franuni() * nx as f32 / 2.0f32) as i32;
            let it = if tc == 2 {
                (franuni() * nz as f32 / 2.0f32) as i32
            } else {
                0
            };
            ca[ix as usize][it as usize] = 1.0;
        }
    }

    /* build ca matrix */
    for it in 1..nz as usize {
        for ix in 1..(nx - 1) as usize {
            /* get 3 vals centered above this location */
            let a = ca[ix - 1][it - 1] as i32;
            let b = ca[ix][it - 1] as i32;
            let c = ca[ix + 1][it - 1] as i32;

            /* map the pattern into an integer code
             * 111,110,101,100,011,010,001,000  <-- input pattern
             * 111,110,101,100, 11, 10,  1,  0  <-- code
             */
            if a != 0 && b != 0 && c != 0 {
                code = 111;
            }
            if a != 0 && b != 0 && c == 0 {
                code = 110;
            }
            if a != 0 && b == 0 && c != 0 {
                code = 101;
            }
            if a != 0 && b == 0 && c == 0 {
                code = 100;
            }
            if a == 0 && b != 0 && c != 0 {
                code = 11;
            }
            if a == 0 && b != 0 && c == 0 {
                code = 10;
            }
            if a == 0 && b == 0 && c != 0 {
                code = 1;
            }
            if a == 0 && b == 0 && c == 0 {
                code = 0;
            }

            /* apply ca rule to get value at this location */
            let val = rule_ca(rule, code);

            /* update ca matrix */
            if val == 1.0 {
                ca[ix][it] += val;
            }
        }
    }

    /* pull off traces and output */
    let mut wht = 0.0f32;
    let mut blk = 0.0f32;
    let mut count = 0i32;
    for ix in nx / 4..3 * nx / 4 {
        for it in 0..nz {
            if it > nz / 2 && ca[ix as usize][it as usize] == 0.0 {
                wht = (wht as f64 + 1.0) as f32;
            }
            if it > nz / 2 && ca[ix as usize][it as usize] != 0.0 {
                blk = (blk as f64 + 1.0) as f32;
            }
            if fill == 0 {
                data[it as usize] = ca[ix as usize][it as usize];
            } else {
                if ca[ix as usize][it as usize] == 0.0 {
                    data[it as usize] = f0;
                }
                if ca[ix as usize][it as usize] != 0.0 {
                    data[it as usize] = f1;
                }
            }
        }
        if count < n2 {
            // C uses fwrite(data,n1,sizeof(float),stdout): the byte count is
            // still n1*sizeof(float), and the return value is ignored.
            fwrite(&data[..n1 as usize], &stdout());
        }
        count += 1;
    }

    // Deliberately unconditional: verbose is documented but never read in C.
    // The original does not call checkpars(), either.
    warn!("bottom half porosity = %f", 100.0f32 * wht / (wht + blk));
    exit(CWP_Exit());
}

/// Cellular automaton value for a pattern code (Chris Liner, 10/21/2006).
fn rule_ca(rule: i32, code: i32) -> f32 {
    if rule == 30 {
        match code {
            111 => return 0.0,
            110 => return 0.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 54 {
        match code {
            111 => return 0.0,
            110 => return 0.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 0.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 60 {
        match code {
            111 => return 0.0,
            110 => return 0.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 0.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 62 {
        match code {
            111 => return 0.0,
            110 => return 0.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 90 {
        match code {
            111 => return 0.0,
            110 => return 1.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 0.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 94 {
        match code {
            111 => return 0.0,
            110 => return 1.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 102 {
        match code {
            111 => return 0.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 0.0,
            11 => return 0.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 110 {
        match code {
            111 => return 0.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 0.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 122 {
        match code {
            111 => return 0.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 0.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 126 {
        match code {
            111 => return 0.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 150 {
        match code {
            111 => return 1.0,
            110 => return 0.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 0.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 158 {
        match code {
            111 => return 1.0,
            110 => return 0.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 182 {
        match code {
            111 => return 1.0,
            110 => return 0.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 0.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 188 {
        match code {
            111 => return 1.0,
            110 => return 0.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 0.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 190 {
        match code {
            111 => return 1.0,
            110 => return 0.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 220 {
        match code {
            111 => return 1.0,
            110 => return 1.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 0.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 222 {
        match code {
            111 => return 1.0,
            110 => return 1.0,
            101 => return 0.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 225 {
        match code {
            111 => return 1.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 0.0,
            11 => return 0.0,
            10 => return 0.0,
            1 => return 0.0,
            0 => return 1.0,
            _ => {}
        }
    }
    if rule == 226 {
        match code {
            111 => return 1.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 0.0,
            11 => return 0.0,
            10 => return 0.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 250 {
        match code {
            111 => return 1.0,
            110 => return 1.0,
            101 => return 1.0,
            100 => return 1.0,
            11 => return 1.0,
            10 => return 0.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    if rule == 260 {
        // x(it+1,ix) = x(it,ix-1) - x(it,ix) + x(it,ix+1)
        match code {
            111 => return 1.0,
            110 => return 0.0,
            101 => return 2.0,
            100 => return 1.0,
            11 => return 0.0,
            10 => return -1.0,
            1 => return 1.0,
            0 => return 0.0,
            _ => {}
        }
    }
    err!("Bad rule number");
}
