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

//! SURANGE - report nonzero header ranges. Port of su/main/headers/surange.c (Revision 1.19).
//! Credits: Arkansas: Chris Liner; Stanford: Stewart A. Levin; Geocon: Garry Perratt;
//! SEP: Stew Levin; CWP: Jack K. Cohen.

use oxwave_core::prelude::*;
use std::sync::atomic::{AtomicBool, Ordering};

static INPUT_CLOSED: AtomicBool = AtomicBool::new(false);

unsafe extern "C" {
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn open(path: *const std::ffi::c_char, flags: i32, ...) -> i32;
    fn dup2(oldfd: i32, newfd: i32) -> i32;
    fn close(fd: i32) -> i32;
}

// C closeinput reopens stdin on /dev/null. Use async-signal-safe descriptor
// operations: the core input stream may already be borrowed by gettr here.
extern "C" fn closeinput(_sig: i32) {
    INPUT_CLOSED.store(true, Ordering::Relaxed);
    unsafe {
        let fd = open(c"/dev/null".as_ptr(), 0);
        if fd >= 0 {
            dup2(fd, 0);
            if fd != 0 {
                close(fd);
            }
        }
    }
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SURANGE - get max and min values for non-zero header entries\t",
    " \t\t\t\t\t\t\t\t",
    " surange <stdin\t \t\t\t\t\t",
    "\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t",
    "\tkey=\t\tHeader key(s) to range (default=all)\t",
    "\tdim=0\t\tdim seismic flag\t",
    "\t    \t\t0 = not dim, 1 = coord in ft, 2 = coord in m\t",
    " \t\t\t\t\t\t\t\t",
    " Note: Gives partial results if interrupted\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Output is: \t\t\t\t\t\t\t",
    " number of traces \t\t\t\t\t\t",
    " keyword min max (first - last) \t\t\t\t",
    " north-south-east-west limits of shot, receiver and midpoint   ",
    " if dim then also midpoint interval and line length   ",
    " \t\t\t\t\t\t\t\t",
];

fn printrange(min: &Trace, max: &Trace, first: &Trace, last: &Trace) {
    for i in 0..SU_NKEYS as i32 {
        let key = getkey(i).unwrap();
        let typ = hdtype(key);
        let (mut a, mut b, mut c, mut d) = (
            Value::default(),
            Value::default(),
            Value::default(),
            Value::default(),
        );
        gethval(min, i, &mut a);
        gethval(max, i, &mut b);
        gethval(first, i, &mut c);
        gethval(last, i, &mut d);
        let da = vtod(typ, a);
        let db = vtod(typ, b);
        if da != 0.0 || db != 0.0 {
            printf!("%-8s ", key);
            printfval(typ, a);
            if da < db {
                printf!(" ");
                printfval(typ, b);
                printf!(" (");
                printfval(typ, c);
                printf!(" - ");
                printfval(typ, d);
                printf!(")");
            }
            putchar(b'\n' as i32);
        }
    }
}

// North, south, east and west each store the complete coordinate pair.
fn update_limits(limits: &mut [[f64; 2]; 4], x: f64, y: f64) {
    if limits[0][1] < y {
        limits[0] = [x, y];
    }
    if limits[1][1] > y {
        limits[1] = [x, y];
    }
    if limits[2][0] < x {
        limits[2] = [x, y];
    }
    if limits[3][0] > x {
        limits[3] = [x, y];
    }
}
fn print_limits(name: &str, a: &[[f64; 2]; 4]) {
    if a[0][1] != 0.0 || a[1][1] != 0.0 || a[2][0] != 0.0 || a[3][0] != 0.0 {
        printf!(
            "\n%s coordinate limits:\n\tNorth(%g,%g) South(%g,%g) East(%g,%g) West(%g,%g)\n",
            name,
            a[0][0],
            a[0][1],
            a[1][0],
            a[1][1],
            a[2][0],
            a[2][1],
            a[3][0],
            a[3][1]
        );
    }
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let nkeys = countparval("key");
    let mut keys = vec![String::new(); SU_NKEYS];
    if nkeys != 0 {
        getparstringarray("key", &mut keys);
    }
    let mut dim = 0;
    if getparint("dim", &mut dim) == 0 {
        dim = 0;
    }
    checkpars();
    let (mut tr, mut trmin, mut trmax, mut trfirst, mut trlast) = (
        Trace::new(),
        Trace::new(),
        Trace::new(),
        Trace::new(),
        Trace::new(),
    );
    // Set up closing commands, as in the C program (SIGINT and SIGTERM).
    unsafe {
        signal(2, closeinput);
        signal(15, closeinput);
    }
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let mut shot = [[0.0; 2]; 4];
    let mut rec = shot;
    let mut cmp = shot;
    let mut scale = 1.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut gx = 0.0;
    let mut gy = 0.0;
    let indices: Vec<i32> = if nkeys == 0 {
        (0..SU_NKEYS as i32).collect()
    } else {
        keys[..nkeys as usize].iter().map(|k| getindex(k)).collect()
    };
    let mut val = Value::default();
    for &i in &indices {
        gethval(&tr, i, &mut val);
        puthval(&mut trmin, i, &val);
        puthval(&mut trmax, i, &val);
        puthval(&mut trfirst, i, &val);
        if nkeys == 0 {
            if i == 20 {
                let coscal = if val.h() == 0 { 1 } else { val.h() as i32 };
                scale = if coscal > 0 {
                    coscal as f64
                } else {
                    -1.0 / coscal as f64
                };
            }
            if i == 21 {
                sx = val.i() as f64 * scale;
            }
            if i == 22 {
                sy = val.i() as f64 * scale;
            }
            if i == 23 {
                gx = val.i() as f64 * scale;
            }
            if i == 24 {
                gy = val.i() as f64 * scale;
            }
        }
    }
    if nkeys == 0 {
        shot = [[sx, sy]; 4];
        rec = [[gx, gy]; 4];
        cmp = [[0.5 * (sx + gx), 0.5 * (sy + gy)]; 4];
    }
    let mut ntr: i32 = 1;
    let (mut mx1, mut my1, mut dmin, mut dmax, mut davg) = (0.0, 0.0, 0.0, 0.0, 0.0);
    // freopen also discards read-ahead. Stop before consuming any headers
    // still buffered by the core stream when closeinput ran.
    while !INPUT_CLOSED.load(Ordering::Relaxed) && gettr(&mut tr) != 0 {
        sx = 0.0;
        sy = 0.0;
        gx = 0.0;
        gy = 0.0;
        for &i in &indices {
            let typ = hdtype(getkey(i).unwrap());
            let (mut vmin, mut vmax) = (Value::default(), Value::default());
            gethval(&tr, i, &mut val);
            gethval(&trmin, i, &mut vmin);
            gethval(&trmax, i, &mut vmax);
            if valcmp(typ, val, vmin) < 0 {
                puthval(&mut trmin, i, &val);
            }
            if valcmp(typ, val, vmax) > 0 {
                puthval(&mut trmax, i, &val);
            }
            puthval(&mut trlast, i, &val);
            if nkeys == 0 {
                if i == 20 {
                    let coscal = if val.h() == 0 { 1 } else { val.h() as i32 };
                    scale = if coscal > 0 {
                        coscal as f64
                    } else {
                        -1.0 / coscal as f64
                    };
                }
                if i == 21 {
                    sx = val.i() as f64 * scale;
                }
                if i == 22 {
                    sy = val.i() as f64 * scale;
                }
                if i == 23 {
                    gx = val.i() as f64 * scale;
                }
                if i == 24 {
                    gy = val.i() as f64 * scale;
                }
            }
        }
        if nkeys == 0 {
            update_limits(&mut shot, sx, sy);
            update_limits(&mut rec, gx, gy);
            update_limits(&mut cmp, 0.5 * (sx + gx), 0.5 * (sy + gy));
        }
        let mx2 = 0.5 * tr.sx.wrapping_add(tr.gx) as f64;
        let my2 = 0.5 * tr.sy.wrapping_add(tr.gy) as f64;
        if ntr == 1 {
            mx1 = mx2;
            my1 = my2;
        } else {
            let dm = ((mx1 - mx2) * (mx1 - mx2) + (my1 - my2) * (my1 - my2)).sqrt();
            if ntr == 2 {
                dmin = dm;
                dmax = dm;
                davg = (dmin + dmax) / 2.0;
            } else {
                if dm < dmin {
                    dmin = dm;
                }
                if dm > dmax {
                    dmax = dm;
                }
                davg = (davg + (dmin + dmax) / 2.0) / 2.0;
            }
            mx1 = mx2;
            my1 = my2;
        }
        ntr = ntr.wrapping_add(1);
    }
    printf!("%d traces:\n", ntr);
    printrange(&trmin, &trmax, &trfirst, &trlast);
    if nkeys == 0 {
        print_limits("Shot", &shot);
        print_limits("Receiver", &rec);
        print_limits("Midpoint", &cmp);
    }
    if dim == 1 {
        printf!("\n2D line: \n");
        printf!("Min CMP interval = %g ft\n", dmin);
        printf!("Max CMP interval = %g ft\n", dmax);
        printf!(
            "Line length = %g miles (using avg CMP interval of %g ft)\n",
            davg * ntr as f64 / 5280.0,
            davg
        );
    } else if dim == 2 {
        printf!("ddim line: \n");
        printf!("Min CMP interval = %g m\n", dmin);
        printf!("Max CMP interval = %g m\n", dmax);
        printf!(
            "Line length = %g km (using avg CMP interval of %g m)\n",
            davg * ntr as f64 / 1000.0,
            davg
        );
    }
    exit(CWP_Exit());
}
