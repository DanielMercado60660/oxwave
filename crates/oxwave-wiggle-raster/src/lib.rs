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

//! Ports of cwp/lib/rfwtva.c and cwp/lib/rfwtvaint.c.
//! Credits (from the C sources):
//! Author:  Dave Hale, Colorado School of Mines, 07/01/89
//! Modified:  Craig Artley, Colorado School of Mines, 04/14/92
//!            Fixed bug in computing yoffset.  Previously, when zmin==zmax
//!            the rasterized trace was shifted to the left by one trace.
//! rfwtvaint.c additionally credits:
//! Memorial University of Newfoundland: Tony Kocurko, Sept 1995.
//! Added sinc interpolation.
//! MODIFIED:  Paul Michaels, Boise State University, 29 December 2000
//!            Added solid/grey shading scheme, wiggle>=2 option for peaks/troughs
//! MODIFIED:  Kris Vanneste, Royal Observatory of Belgium, July 2005
//!            Added missing wiggle trace bits in y (amplitude) direction
//!
//! The input slice length replaces C's `n` and must fit in `int`.
//! Projected amplitudes must fit in `int`. `nbpr` must be positive,
//! the clipped horizontal interval must be nonempty, and the raster slice
//! must cover every addressed byte. Grey shading preserves C's byte/bit
//! stepping, which can touch one byte beyond a row (including the last row).
//! An attempted write outside the supplied slice exits through `err!`;
//! provide the extra byte when this shading is needed. Invalid sample counts,
//! row widths, scan-line counts, and scratch allocation failures also use
//! `err!` instead of panicking. An empty [`rfwtva`] trace is a no-op;
//! [`rfwtvaint`] requires at least one sample.

use oxwave_sinc_interp::ints8r;

// Keep C's byte and bit cursors separate: grey shading resets the bit to 7
// on crossing a byte, rather than recomputing it from y.
fn set_bit(bits: &mut [u8], byte: i32, bit: i32, endian: i32) {
    if endian == 0 || endian == 1 {
        // C can step past the last row during grey fill. Preserve the write
        // when the caller supplied that byte, and reject memory corruption.
        let target = bits.get_mut(byte as usize).unwrap_or_else(|| {
            oxwave_core::err!("wiggle raster: raster byte outside supplied buffer")
        });
        *target |= 1 << if endian == 0 { -bit + 7 } else { bit };
    } else {
        oxwave_core::fprintf!(
            oxwave_core::cio::stderr(),
            "endian must equal either 0 or 1\n"
        );
    }
}

fn interpolation_scratch(len: usize) -> Vec<f32> {
    let mut scratch = Vec::new();
    scratch.try_reserve_exact(len).unwrap_or_else(|_| {
        oxwave_core::err!("rfwtvaint: cannot allocate interpolation scratch buffer")
    });
    // Capacity was reserved fallibly, so resize cannot allocate.
    scratch.resize(len, 0.0);
    scratch
}

/// Rasterize samples as wiggle trace / variable area, using linear segments.
///
/// C prototype:
/// `void rfwtva(int n, float z[], float zmin, float zmax, float zbase,
/// int yzmin, int yzmax, int xfirst, int xlast, int wiggle, int nbpr,
/// unsigned char *bits, int endian);`
///
/// Coordinates increase down (`x`) and right (`y`). Amplitudes map from
/// `zmin..zmax` to `yzmin..yzmax`; reversing either interval reverses polarity.
/// `wiggle == 0` fills only, nonzero draws the trace, and values >= 2 also
/// shade troughs (2 light grey through 5 black). `endian` selects low-bit-first
/// (0) or high-bit-first (1); other values print C's diagnostic for each bit.
/// Existing raster bits are retained. See the module's buffer requirements.
/// An empty sample slice leaves the raster unchanged.
#[allow(clippy::too_many_arguments)]
pub fn rfwtva(
    z: &[f32],
    zmin: f32,
    zmax: f32,
    zbase: f32,
    yzmin: i32,
    yzmax: i32,
    xfirst: i32,
    xlast: i32,
    mut wiggle: i32,
    nbpr: i32,
    bits: &mut [u8],
    endian: i32,
) {
    if z.is_empty() {
        return;
    }
    let n = i32::try_from(z.len())
        .unwrap_or_else(|_| oxwave_core::err!("rfwtva: sample count exceeds int range"));
    if nbpr <= 0 {
        oxwave_core::err!("rfwtva: nbpr must be positive");
    }
    let igrey = if wiggle >= 2 {
        let igrey = wiggle.wrapping_abs();
        wiggle = 1;
        igrey
    } else {
        0
    };
    let mut ymin = if yzmin < yzmax { yzmin } else { yzmax };
    let mut ymax = if yzmax > yzmin { yzmax } else { yzmin };
    ymin = if ymin > 0 { ymin } else { 0 };
    let right = nbpr.wrapping_mul(8).wrapping_sub(1);
    ymax = if ymax < right { ymax } else { right };
    let iscale = n - 1;
    let yscale = if zmax != zmin {
        yzmax.wrapping_sub(yzmin) as f32 / (zmax - zmin)
    } else {
        1.0
    };
    let yoffset = if zmax != zmin {
        yzmin as f32 - zmin * yscale
    } else {
        (0.5 * yzmin.wrapping_add(yzmax) as f64) as f32
    };
    let xscale = if n > 1 { xlast.wrapping_sub(xfirst) } else { 0 };
    let xstep = if xlast > xfirst { 1 } else { -1 };
    let mut ybase = (yoffset + zbase * yscale) as i32;
    ybase = if ybase > ymin { ybase } else { ymin };
    ybase = if ybase < ymax { ybase } else { ymax };
    let mut znext = z[0];
    let mut ynext = (yoffset + znext * yscale) as i32;
    let mut xnext = xfirst;

    for i in 0..n {
        let xthis = xnext;
        xnext = if i < iscale {
            xfirst.wrapping_add((i + 1).wrapping_mul(xscale) / iscale)
        } else {
            xthis.wrapping_add(xstep)
        };
        let dx = xnext.wrapping_sub(xthis);
        let zthis = znext;
        znext = if i < n - 1 {
            z[(i + 1) as usize]
        } else {
            zthis
        };
        let ythis = ynext;
        ynext = (yoffset + znext * yscale) as i32;

        // Coincident scan lines still draw the intervening amplitudes.
        if dx == 0 {
            if wiggle != 0 {
                let ystep = if ythis < ynext { 1 } else { -1 };
                let mut y2 = ythis;
                while y2 != ynext {
                    let mut clip = false;
                    if y2 < ymin {
                        y2 = ymin;
                        clip = true;
                    }
                    if y2 > ymax {
                        y2 = ymax;
                        clip = true;
                    }
                    let bit = 7i32.wrapping_sub(y2) & 7;
                    let byte = xthis.wrapping_mul(nbpr).wrapping_add(y2 >> 3);
                    set_bit(bits, byte, bit, endian);
                    if clip {
                        break;
                    }
                    y2 = y2.wrapping_add(ystep);
                }
            } else {
                continue;
            }
        }
        let dy = ynext.wrapping_sub(ythis);
        let mut x = xthis;
        let mut y = ythis;
        while x != xnext {
            // dx cannot be zero in this loop. wrapping_div also matches
            // arm64's signed MIN / -1 result if coordinates overflow.
            if wiggle != 0 {
                let yxthis = y;
                let yxnext = ythis.wrapping_add(
                    x.wrapping_sub(xthis)
                        .wrapping_add(xstep)
                        .wrapping_mul(dy)
                        .wrapping_div(dx),
                );
                let ystep = if yxthis < yxnext { 1 } else { -1 };
                let mut y2 = yxthis;
                while y2 != yxnext.wrapping_add(ystep) {
                    let mut clip = false;
                    if y2 < ymin {
                        y2 = ymin;
                        clip = true;
                    }
                    if y2 > ymax {
                        y2 = ymax;
                        clip = true;
                    }
                    let bit = 7i32.wrapping_sub(y2) & 7;
                    let byte = x.wrapping_mul(nbpr).wrapping_add(y2 >> 3);
                    set_bit(bits, byte, bit, endian);
                    if clip {
                        break;
                    }
                    y2 = y2.wrapping_add(ystep);
                }
            }
            if y < ymin {
                y = ymin;
            }
            if y > ymax {
                y = ymax;
            }
            let mut bit = 7i32.wrapping_sub(y) & 7;
            let mut byte = x.wrapping_mul(nbpr).wrapping_add(y >> 3);
            if y > ybase {
                set_bit(bits, byte, bit, endian);
            }
            while y > ybase {
                y = y.wrapping_sub(1);
                bit += 1;
                if bit >= 8 {
                    byte = byte.wrapping_sub(1);
                    bit = 0;
                }
                set_bit(bits, byte, bit, endian);
            }
            if igrey > 0 {
                let mut ideci = 6i32.wrapping_sub(igrey);
                if ideci < 1 {
                    ideci = 1;
                }
                while y < ybase {
                    y = y.wrapping_add(ideci);
                    bit -= ideci;
                    if bit < 0 {
                        byte = byte.wrapping_add(1);
                        bit = 7;
                    }
                    set_bit(bits, byte, bit, endian);
                }
            }
            x = x.wrapping_add(xstep);
            y = ythis.wrapping_add(x.wrapping_sub(xthis).wrapping_mul(dy).wrapping_div(dx));
        }
    }
}

/// Rasterize samples after eight-point sinc interpolation.
///
/// C prototype:
/// `void rfwtvaint(int n, float z[], float zmin, float zmax, float zbase,
/// int yzmin, int yzmax, int xfirst, int xlast, int wiggle, int nbpr,
/// unsigned char *bits, int endian);`
///
/// Amplitude mapping, shading and bit order match [`rfwtva`]. Like C, this
/// writes `abs(xlast - xfirst + 1)` rows starting at row zero, regardless of
/// `xfirst`; reversed coordinates also reverse the interpolated samples.
/// A zero row count leaves the raster unchanged. Scratch arrays are local
/// instead of C's static reusable buffers. See the module's requirements.
#[allow(clippy::too_many_arguments)]
pub fn rfwtvaint(
    z: &[f32],
    zmin: f32,
    zmax: f32,
    zbase: f32,
    yzmin: i32,
    yzmax: i32,
    xfirst: i32,
    xlast: i32,
    mut wiggle: i32,
    nbpr: i32,
    bits: &mut [u8],
    endian: i32,
) {
    let n = i32::try_from(z.len())
        .unwrap_or_else(|_| oxwave_core::err!("rfwtvaint: sample count exceeds int range"));
    if n == 0 {
        oxwave_core::err!("rfwtvaint: sample array must not be empty");
    }
    if nbpr <= 0 {
        oxwave_core::err!("rfwtvaint: nbpr must be positive");
    }
    let igrey = if wiggle >= 2 {
        let igrey = wiggle.wrapping_abs();
        wiggle = 1;
        igrey
    } else {
        0
    };
    // Check each C signed operation: overflowing subtraction, addition, or
    // ABS(INT_MIN) is undefined, and must not turn into a huge allocation.
    let nxout = xlast
        .checked_sub(xfirst)
        .and_then(|span| span.checked_add(1))
        .and_then(i32::checked_abs)
        .unwrap_or_else(|| oxwave_core::err!("rfwtvaint: scan-line count exceeds int range"));
    let mut xout = interpolation_scratch(nxout as usize);
    let mut yout = interpolation_scratch(nxout as usize);
    let dxin = 1.0f32;
    let fxin = 0.0f32;
    let yin = z;
    let yinl = z[0];
    let yinr = z[(n - 1) as usize];
    let mut deltax = (n - 1) as f32 / nxout.wrapping_sub(1) as f32;
    let mut x0 = 0.0f32;
    if xfirst > xlast {
        x0 = (n - 1) as f32;
        deltax = -deltax;
    }
    for i in 0..nxout {
        xout[i as usize] = x0 + i as f32 * deltax;
    }
    ints8r(dxin, fxin, yin, yinl, yinr, &xout, &mut yout);

    let mut ymin = if yzmin < yzmax { yzmin } else { yzmax };
    let mut ymax = if yzmax > yzmin { yzmax } else { yzmin };
    ymin = if ymin > 0 { ymin } else { 0 };
    let right = nbpr.wrapping_mul(8).wrapping_sub(1);
    ymax = if ymax < right { ymax } else { right };
    let yscale = if zmax != zmin {
        yzmax.wrapping_sub(yzmin) as f32 / (zmax - zmin)
    } else {
        1.0
    };
    let yoffset = if zmax != zmin {
        yzmin as f32 - zmin * yscale
    } else {
        (0.5 * yzmin.wrapping_add(yzmax) as f64) as f32
    };
    let mut ybase = (yoffset + zbase * yscale) as i32;
    ybase = if ybase > ymin { ybase } else { ymin };
    ybase = if ybase < ymax { ybase } else { ymax };

    for i in 0..nxout {
        let zthis = yout[i as usize];
        let mut y = (yoffset + zthis * yscale) as i32;
        if y < ymin {
            y = ymin;
        }
        if y > ymax {
            y = ymax;
        }
        if wiggle != 0 && i > 0 {
            let mut y_prev = (yoffset + yout[(i - 1) as usize] * yscale) as i32;
            if y_prev < ymin {
                y_prev = ymin;
            } else if y_prev > ymax {
                y_prev = ymax;
            }
            let ystep = if y < y_prev { 1 } else { -1 };
            if y != y_prev {
                let mut y2 = y.wrapping_add(ystep);
                while y2 != y_prev {
                    let bit = 7i32.wrapping_sub(y2) & 7;
                    let byte = i.wrapping_mul(nbpr).wrapping_add(y2 >> 3);
                    set_bit(bits, byte, bit, endian);
                    y2 = y2.wrapping_add(ystep);
                }
            }
        }
        let mut bit = 7i32.wrapping_sub(y) & 7;
        let mut byte = i.wrapping_mul(nbpr).wrapping_add(y >> 3);
        if wiggle != 0 || y > ybase {
            set_bit(bits, byte, bit, endian);
        }
        while y > ybase {
            y = y.wrapping_sub(1);
            bit += 1;
            if bit >= 8 {
                byte = byte.wrapping_sub(1);
                bit = 0;
            }
            set_bit(bits, byte, bit, endian);
        }
        if igrey > 0 {
            let mut ideci = 6i32.wrapping_sub(igrey);
            if ideci < 1 {
                ideci = 1;
            }
            while y < ybase {
                y = y.wrapping_add(ideci);
                bit -= ideci;
                if bit < 0 {
                    byte = byte.wrapping_add(1);
                    bit = 7;
                }
                set_bit(bits, byte, bit, endian);
            }
        }
    }
}
