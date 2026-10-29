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

//! SWAPBYTES - swap the bytes of various data types.
//! Port of par/main/data_conversion/swapbytes.c (Revision 1.12).
//! Credits (from the C source): CWP: John Stockwell (Jan 1994).
//! Institut fur Geophysik, Hamburg: Jens Hartmann supplied byte swapping subroutines.

use oxwave_core::prelude::*;
use oxwave_core::segy::{swap_float_4, swap_double_8, swap_short_2, swap_u_short_2,
    swap_long_4, swap_u_long_4, swap_int_4};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SWAPBYTES - SWAP the BYTES of various  data types\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " swapbytes <stdin [optional parameters]  >stdout \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " \tnone\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " in=float\tinput type\t(float)\t\t\t\t\t",
    " \t\t=double\t\t(double)\t\t\t\t",
    " \t\t=short\t\t(short)\t\t\t\t\t",
    " \t\t=ushort\t\t(unsigned short)\t\t\t",
    " \t\t=long\t\t(long)\t\t\t\t\t",
    " \t\t=ulong\t\t(unsigned long)\t\t\t\t",
    " \t\t=int\t\t(int)\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " outpar=/dev/tty\t\toutput parameter file, contains the\t",
    "\t\t\t\tnumber of values (n1=)\t\t\t",
    " \t\t\tother choices for outpar are: /dev/tty,\t\t",
    " \t\t\t/dev/stderr, or a name of a disk file\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The byte order of the mantissa of binary data values on PC's and DEC's",
    " is the reverse of so called \"big endian\" machines (IBM RS6000, SUN,etc.)",
    " hence the need for byte-swapping capability. The subroutines in this code",
    " have been tested for swapping between PCs and\tbig endian machines, but",
    " have not been tested for DEC products.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Caveat:\t\t\t\t\t\t\t\t",
    " 2 byte short, 4 byte long, 4 byte float, 4 byte int,\t\t\t",
    " and 8 bit double assumed.\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);

    let mut outpar = String::new();
    if getparstring("outpar", &mut outpar) == 0 { outpar = "/dev/tty".into(); }
    let outparfp = efopen(&outpar, "w");
    let mut input = String::new();
    if getparstring("in", &mut input) == 0 { input = "float".into(); }
    checkpars();

    if !matches!(input.as_str(), "float" | "double" | "short" | "ushort" | "long" | "ulong" | "int") {
        err!("%s is an unsupported type", input);
    }

    let mut n: i32 = 0;
    macro_rules! convert {
        ($ty:ty, $swap:ident) => {{
            let mut x: [$ty; 1] = [0 as $ty];
            while efread(&mut x, &stdin()) != 0 {
                n = n.wrapping_add(1);
                $swap(&mut x[0]);
                efwrite(&x, &stdout());
            }
        }};
    }
    match input.as_str() {
        "float" => convert!(f32, swap_float_4),
        "double" => convert!(f64, swap_double_8),
        "short" => convert!(i16, swap_short_2),
        "ushort" => convert!(u16, swap_u_short_2),
        "long" => convert!(i64, swap_long_4),
        "ulong" => convert!(u64, swap_u_long_4),
        "int" => convert!(i32, swap_int_4),
        _ => unreachable!(),
    }
    fprintf!(outparfp, "n=%d\n", n);
    exit(CWP_Exit());
}
