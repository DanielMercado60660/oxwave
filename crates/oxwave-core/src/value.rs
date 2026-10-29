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

//! `value`: SU's `Value` union (su.h) and the routines of `valpkge.c`.
//!
//! `Value` is an 8-byte union in C (`s,h,u,l,v,i,p,f,d,U,P` all at offset 0).
//! Here it is an 8-byte little-endian cell with typed accessors, so reading
//! a different member than the one written behaves exactly like the C union
//! on this machine (e.g. `zeroval.l = 0` then comparing as any type).
//!
//! | C                    | Rust                          |
//! |----------------------|-------------------------------|
//! | `val.h = x`          | `val.set_h(x)`                |
//! | `val.h`              | `val.h()`                     |
//! | `vtoi(type, val)`    | `vtoi(type, val)` (type: `&str`) |
//! | `vtol/vtof/vtod`     | same                          |
//! | `valcmp(type,a,b)`   | `valcmp(type, a, b) -> i32`   |
//! | `printfval(type,v)`  | `printfval(type, v)`          |
//! | `fprintfval(fp,t,v)` | `fprintfval(&fp, t, v)`       |
//! | `atoval(t, s, &v)`   | `atoval(t, s, &mut v)`        |
//! | `getparval(n,t,k,v)` | `getparval(n, t, k, &mut v[..])` |
//! | `valtoabs(t, v)`     | `valtoabs(t, v)`              |
//! | `printftype(key)`    | `printftype(key)`             |

use crate::cio::{stdout, CFile};
use crate::par;

/// The C `Value` union: 8 bytes, members aliased at offset 0 (little-endian).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Value {
    pub bytes: [u8; 8],
}

macro_rules! member {
    ($get:ident, $set:ident, $t:ty, $n:expr) => {
        #[inline]
        pub fn $get(&self) -> $t {
            let mut b = [0u8; $n];
            b.copy_from_slice(&self.bytes[..$n]);
            <$t>::from_le_bytes(b)
        }
        #[inline]
        pub fn $set(&mut self, x: $t) {
            self.bytes[..$n].copy_from_slice(&x.to_le_bytes());
        }
    };
}

impl Value {
    member!(h, set_h, i16, 2);
    member!(u, set_u, u16, 2);
    member!(l, set_l, i64, 8);
    member!(v, set_v, u64, 8);
    member!(i, set_i, i32, 4);
    member!(p, set_p, u32, 4);
    member!(f, set_f, f32, 4);
    member!(d, set_d, f64, 8);
    /// `unsigned int U:16` bit-field (low 16 bits of the first word).
    #[allow(non_snake_case)]
    pub fn U(&self) -> u32 {
        self.p() & 0xFFFF
    }
    /// `unsigned int P:32` bit-field.
    #[allow(non_snake_case)]
    pub fn P(&self) -> u32 {
        self.p()
    }
    /// `char s[8]` member, up to (not including) the first NUL.
    pub fn s(&self) -> &[u8] {
        let n = self.bytes.iter().position(|&c| c == 0).unwrap_or(8);
        &self.bytes[..n]
    }
    /// `strcpy(val.s, str)` (truncated to the 8 bytes of the union).
    pub fn set_s(&mut self, s: &str) {
        let b = s.as_bytes();
        let n = b.len().min(8);
        self.bytes[..n].copy_from_slice(&b[..n]);
        if n < 8 {
            self.bytes[n] = 0;
        }
    }
}

fn tc(type_: &str) -> u8 {
    type_.as_bytes().first().copied().unwrap_or(0)
}

pub fn vtoi(type_: &str, val: Value) -> i32 {
    match tc(type_) {
        b's' => val.bytes[0] as i8 as i32,
        b'h' => val.h() as i32,
        b'u' => val.u() as i32,
        b'l' => val.l() as i32,
        b'v' => val.v() as i32,
        b'i' => val.i(),
        b'p' => val.p() as i32,
        b'f' => val.f() as i32,
        b'U' => val.U() as i32,
        b'P' => val.P() as i32,
        _ => crate::err!("vtoi: unknown type %s", type_),
    }
}

pub fn vtol(type_: &str, val: Value) -> i64 {
    match tc(type_) {
        b's' => val.bytes[0] as i8 as i64,
        b'h' => val.h() as i64,
        b'u' => val.u() as i64,
        b'l' => val.l(),
        b'v' => val.v() as i64,
        b'i' => val.i() as i64,
        b'p' => val.p() as i64,
        b'f' => val.f() as i64,
        b'd' => val.d() as i64,
        b'U' => val.U() as i64,
        b'P' => val.P() as i64,
        _ => crate::err!("vtol: unknown type %s", type_),
    }
}

pub fn vtof(type_: &str, val: Value) -> f32 {
    match tc(type_) {
        b's' => val.bytes[0] as i8 as f32,
        b'h' => val.h() as f32,
        b'u' => val.u() as f32,
        b'l' => val.l() as f32,
        b'v' => val.v() as f32,
        b'i' => val.i() as f32,
        b'p' => val.p() as f32,
        b'f' => val.f(),
        b'd' => val.d() as f32,
        b'U' => val.U() as f32,
        b'P' => val.P() as f32,
        _ => crate::err!("vtof: unknown type %s", type_),
    }
}

pub fn vtod(type_: &str, val: Value) -> f64 {
    match tc(type_) {
        b's' => val.bytes[0] as i8 as f64,
        b'h' => val.h() as f64,
        b'u' => val.u() as f64,
        b'l' => val.l() as f64,
        b'v' => val.v() as f64,
        b'i' => val.i() as f64,
        b'p' => val.p() as f64,
        b'f' => val.f() as f64,
        b'd' => val.d(),
        b'U' => val.U() as f64,
        b'P' => val.P() as f64,
        _ => crate::err!("vtod: unknown type %s", type_),
    }
}

fn cmp3<T: PartialOrd>(a: T, b: T) -> i32 {
    if a < b {
        -1
    } else if a > b {
        1
    } else {
        0
    }
}

/// C `strcmp` on the 8-byte `s` member.
fn strcmp8(a: &Value, b: &Value) -> i32 {
    for k in 0..8 {
        let (x, y) = (a.bytes[k], b.bytes[k]);
        if x != y {
            return x as i32 - y as i32;
        }
        if x == 0 {
            return 0;
        }
    }
    0
}

pub fn valcmp(type_: &str, val1: Value, val2: Value) -> i32 {
    match tc(type_) {
        b's' => strcmp8(&val1, &val2),
        b'h' => cmp3(val1.h(), val2.h()),
        b'u' => cmp3(val1.u(), val2.u()),
        b'l' => cmp3(val1.l(), val2.l()),
        b'v' => cmp3(val1.v(), val2.v()),
        b'i' => cmp3(val1.i(), val2.i()),
        b'p' => cmp3(val1.p(), val2.p()),
        b'f' => cmp3(val1.f(), val2.f()),
        b'd' => cmp3(val1.d(), val2.d()),
        b'U' => cmp3(val1.U(), val2.U()),
        b'P' => cmp3(val1.P(), val2.P()),
        _ => crate::err!("valcmp: unknown type %s", type_),
    }
}

/// The text `printfval` would print for this value.
pub fn sprintfval(type_: &str, val: Value) -> String {
    match tc(type_) {
        b's' => String::from_utf8_lossy(val.s()).into_owned(),
        b'h' => crate::sprintf!("%d", val.h()),
        b'u' => crate::sprintf!("%d", val.u()),
        b'i' => crate::sprintf!("%d", val.i()),
        b'p' => crate::sprintf!("%d", val.p()),
        b'l' => crate::sprintf!("%ld", val.l()),
        b'v' => crate::sprintf!("%ld", val.v()),
        b'f' => crate::sprintf!("%f", val.f()),
        b'd' => crate::sprintf!("%f", val.d()),
        b'U' => crate::sprintf!("%d", val.U()),
        b'P' => crate::sprintf!("%d", val.P()),
        _ => crate::err!("printfval: unknown type %s", type_),
    }
}

/// printfval: print a value on stdout in SU's per-type format.
pub fn printfval(type_: &str, val: Value) {
    let s = sprintfval(type_, val);
    crate::cio::fputs(&s, &stdout());
}

/// fprintfval: print a value on a stream in SU's per-type format.
pub fn fprintfval(stream: &CFile, type_: &str, val: Value) {
    let s = match tc(type_) {
        b's' | b'h' | b'u' | b'i' | b'p' | b'l' | b'v' | b'f' | b'd' | b'U' | b'P' => {
            sprintfval(type_, val)
        }
        _ => crate::err!("fprintfval: unknown type %s", type_),
    };
    crate::cio::fputs(&s, stream);
}

/// printftype: print the C type name of a header key on stdout.
pub fn printftype(key: &str) {
    let t = crate::header::hdtype(key);
    let s = match tc(t) {
        b's' => "char\n",
        b'h' => "short\n",
        b'u' => "ushort\n",
        b'i' => "int\n",
        b'p' => "unsigned int\n",
        b'l' => "long\n",
        b'v' => "ulong\n",
        b'f' => "float\n",
        b'd' => "double\n",
        b'U' => "unsigned int:16\n",
        b'P' => "unsigned int:32\n",
        _ => crate::err!("printftype: unknown type %s", t),
    };
    crate::cio::fputs(s, &stdout());
}

/// atoval: convert ascii to Value according to the key type (eatoX family).
pub fn atoval(type_: &str, keyval: &str, valp: &mut Value) {
    match tc(type_) {
        b's' => valp.set_s(keyval),
        b'h' => valp.set_h(par::eatoh(keyval)),
        b'u' => valp.set_u(par::eatou(keyval)),
        b'i' => valp.set_i(par::eatoi(keyval)),
        b'p' => valp.set_p(par::eatop(keyval)),
        b'l' => valp.set_l(par::eatol(keyval)),
        b'v' => valp.set_v(par::eatov(keyval)),
        b'f' => valp.set_f(par::eatof(keyval)),
        b'd' => valp.set_d(par::eatod(keyval)),
        _ => crate::err!("%s: %d: mysterious data type: %s", "valpkge.c", 403, keyval),
    }
}

/// getparval: getpar into an array of Values by type ("h","u","i","p","l","v","f","d").
/// Elements beyond the number of values given are left unchanged... except that
/// C copies `n` elements from an uninitialised buffer; here they are zero.
pub fn getparval(name: &str, type_: &str, n: usize, valp: &mut [Value]) {
    match tc(type_) {
        b'h' => {
            let mut h = vec![0i16; n];
            par::getparshort(name, &mut h[..]);
            for k in 0..n {
                valp[k].set_h(h[k]);
            }
        }
        b'u' => {
            let mut u = vec![0u16; n];
            par::getparushort(name, &mut u[..]);
            for k in 0..n {
                valp[k].set_u(u[k]);
            }
        }
        b'i' => {
            let mut x = vec![0i32; n];
            par::getparint(name, &mut x[..]);
            for k in 0..n {
                valp[k].set_i(x[k]);
            }
        }
        b'p' => {
            let mut x = vec![0u32; n];
            par::getparuint(name, &mut x[..]);
            for k in 0..n {
                valp[k].set_p(x[k]);
            }
        }
        b'l' => {
            let mut x = vec![0i64; n];
            par::getparlong(name, &mut x[..]);
            for k in 0..n {
                valp[k].set_l(x[k]);
            }
        }
        b'v' => {
            let mut x = vec![0u64; n];
            par::getparulong(name, &mut x[..]);
            for k in 0..n {
                valp[k].set_v(x[k]);
            }
        }
        b'f' => {
            let mut x = vec![0f32; n];
            par::getparfloat(name, &mut x[..]);
            for k in 0..n {
                valp[k].set_f(x[k]);
            }
        }
        b'd' => {
            let mut x = vec![0f64; n];
            par::getpardouble(name, &mut x[..]);
            for k in 0..n {
                valp[k].set_d(x[k]);
            }
        }
        _ => crate::err!("getparval: %d: mysterious type %s", 457, type_),
    }
}

/// valtoabs: absolute value for a typed value (C `ABS` macro semantics:
/// `x < 0 ? -x : x`, so -0.0 stays -0.0 and SHRT_MIN wraps).
pub fn valtoabs(type_: &str, mut val: Value) -> Value {
    match tc(type_) {
        b'u' | b'v' | b'p' | b'U' | b'P' => {}
        b'h' => {
            let x = val.h();
            val.set_h(if x < 0 { x.wrapping_neg() } else { x })
        }
        b'i' => {
            let x = val.i();
            val.set_i(if x < 0 { x.wrapping_neg() } else { x })
        }
        b'l' => {
            let x = val.l();
            val.set_l(if x < 0 { x.wrapping_neg() } else { x })
        }
        b'f' => {
            let x = val.f();
            val.set_f(if x < 0.0 { -x } else { x })
        }
        b'd' => {
            let x = val.d();
            val.set_d(if x < 0.0 { -x } else { x })
        }
        _ => crate::err!("valtoabs: %d: mysterious type %s", 478, type_),
    }
    val
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn union_aliasing() {
        let mut v = Value::default();
        v.set_l(0);
        assert_eq!(v.h(), 0);
        v.set_i(-1);
        assert_eq!(v.u(), 65535);
        assert_eq!(v.U(), 65535);
        v.set_f(1.0);
        assert_eq!(v.p(), 0x3f80_0000);
    }

    #[test]
    fn formatting() {
        let mut v = Value::default();
        v.set_h(-3);
        assert_eq!(sprintfval("h", v), "-3");
        v.set_u(65535);
        assert_eq!(sprintfval("u", v), "65535");
        v.set_f(0.004);
        assert_eq!(sprintfval("f", v), "0.004000");
        v.set_p(3_000_000_000);
        assert_eq!(sprintfval("p", v), "-1294967296");
        assert_eq!(vtof("p", v), 3e9);
    }

    #[test]
    fn compare_and_abs() {
        let mut a = Value::default();
        let mut b = Value::default();
        a.set_f(-1.0);
        b.set_f(2.0);
        assert_eq!(valcmp("f", a, b), -1);
        assert_eq!(valtoabs("f", a).f(), 1.0);
        a.set_h(i16::MIN);
        assert_eq!(valtoabs("h", a).h(), i16::MIN);
        a.set_f(f32::NAN);
        assert_eq!(valcmp("f", a, b), 0);
        assert_eq!(vtoi("f", a), 0);
    }
}
