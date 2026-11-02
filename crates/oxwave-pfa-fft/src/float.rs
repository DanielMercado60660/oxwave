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

//! Ordered floating operations for the C FFT expressions.
//!
//! Rust permits commutative operands to be exchanged and, for example,
//! multiplication by -1 to become negation. These preserve numerical results
//! but can select a different NaN payload or flip its sign. The reference
//! Darwin/arm64 build executes the source expressions in operand order.

use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy)]
#[repr(transparent)]
pub(super) struct Single(pub f32);

#[derive(Clone, Copy)]
#[repr(transparent)]
pub(super) struct Double(pub f64);

macro_rules! operations {
    ($ty:ident, $reg:tt) => {
        operations!(@op $ty, $reg, Add, add, "fadd", +);
        operations!(@op $ty, $reg, Sub, sub, "fsub", -);
        operations!(@op $ty, $reg, Mul, mul, "fmul", *);
        impl Neg for $ty {
            type Output = Self;
            #[inline(always)]
            fn neg(self) -> Self {
                Self(-self.0)
            }
        }
    };
    (@op $ty:ident, $reg:tt, $trait:ident, $method:ident, $instruction:literal, $op:tt) => {
        impl $trait for $ty {
            type Output = Self;
            #[inline(always)]
            fn $method(self, rhs: Self) -> Self {
                #[cfg(target_arch = "aarch64")]
                {
                    let result;
                    // SAFETY: register-only scalar arithmetic, with no memory
                    // access or stack use. Keep both operands in C order.
                    unsafe {
                        std::arch::asm!(
                            concat!($instruction, " {out:", $reg, "}, {a:", $reg,
                                    "}, {b:", $reg, "}"),
                            out = lateout(vreg) result,
                            a = in(vreg) self.0,
                            b = in(vreg) rhs.0,
                            options(pure, nomem, nostack, preserves_flags),
                        );
                    }
                    Self(result)
                }
                #[cfg(not(target_arch = "aarch64"))]
                {
                    Self(self.0 $op rhs.0)
                }
            }
        }
    };
}

operations!(Single, "s");
operations!(Double, "d");

impl Single {
    #[inline(always)]
    pub(super) fn wide(self) -> Double {
        #[cfg(target_arch = "aarch64")]
        {
            let result;
            // SAFETY: register-only conversion. Rust constant folding of
            // `as f64` is allowed to discard a NaN sign/payload as well.
            unsafe {
                std::arch::asm!(
                    "fcvt {out:d}, {a:s}",
                    out = lateout(vreg) result,
                    a = in(vreg) self.0,
                    options(pure, nomem, nostack, preserves_flags),
                );
            }
            Double(result)
        }
        #[cfg(not(target_arch = "aarch64"))]
        {
            Double(self.0 as f64)
        }
    }
}

impl Double {
    #[inline(always)]
    pub(super) fn narrow(self) -> Single {
        #[cfg(target_arch = "aarch64")]
        {
            let result;
            // SAFETY: register-only conversion, preserving native quieting
            // and payload truncation rather than Rust's unspecified NaN bits.
            unsafe {
                std::arch::asm!(
                    "fcvt {out:s}, {a:d}",
                    out = lateout(vreg) result,
                    a = in(vreg) self.0,
                    options(pure, nomem, nostack, preserves_flags),
                );
            }
            Single(result)
        }
        #[cfg(not(target_arch = "aarch64"))]
        {
            Single(self.0 as f32)
        }
    }
}

pub(super) fn slots(z: &mut [f32]) -> &mut [Single] {
    // SAFETY: Single is transparent over f32, with identical validity and
    // alignment. The exclusive reborrow retains the slice's extent/lifetime.
    unsafe { std::slice::from_raw_parts_mut(z.as_mut_ptr().cast(), z.len()) }
}
