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

//! `cfmt`: C `printf` formatting, so ported programs can keep their C format
//! strings verbatim (`%d`, `%6s`, `%g`, `%a`, `%-10.3e`, `%ld`, `%hd`, `*`, `%1$.6g`, ...).
//! Output is identical to macOS libc `printf` (tested against `snprintf`).
//! The uninitialized-byte tie in Darwin's `%.12a` uses nearest-even rounding here.
//!
//! Macros (exported at the crate root):
//! * `sprintf!(fmt, args..) -> String`   (C `sprintf(buf, fmt, ..)`)
//! * `printf!(fmt, args..)`              (C `printf`, buffered stdout)
//! * `fprintf!(fp, fmt, args..)`         (C `fprintf`, `fp: &CFile`, e.g. `&stderr()`)
//!
//! Arguments are any Rust integer, `f32`/`f64` (f32 is promoted to double as in
//! C varargs), `&str`/`String`, `char`, `bool`.  Integer conversions honour the
//! length modifier exactly like C on LP64: `%d` truncates to `int`, `%hd` to
//! `short`, `%ld` uses 64 bits, `%u` reinterprets as unsigned.  A float passed
//! to `%d` (or an integer to `%f`) is a porting bug and panics.

use std::fmt::Write as _;

/// One printf argument.
#[derive(Clone, Copy, Debug)]
pub enum Arg<'a> {
    /// Any integer, sign- or zero-extended to 64 bits (two's complement bits).
    Int(u64),
    /// Any float, promoted to double.
    Float(f64),
    /// A string (`%s`).
    Str(&'a str),
}

/// Conversion of Rust values to printf arguments.
pub trait CArg {
    fn to_arg(&self) -> Arg<'_>;
}

macro_rules! carg_int {
    ($($t:ty),*) => {$(
        impl CArg for $t {
            fn to_arg(&self) -> Arg<'_> { Arg::Int(*self as i64 as u64) }
        }
    )*};
}
carg_int!(i8, i16, i32, i64, isize, u8, u16, u32);
impl CArg for u64 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Int(*self)
    }
}
impl CArg for usize {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Int(*self as u64)
    }
}
impl CArg for bool {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Int(*self as u64)
    }
}
impl CArg for char {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Int(*self as u64)
    }
}
impl CArg for f32 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Float(*self as f64)
    }
}
impl CArg for f64 {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Float(*self)
    }
}
impl CArg for str {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Str(self)
    }
}
impl CArg for String {
    fn to_arg(&self) -> Arg<'_> {
        Arg::Str(self.as_str())
    }
}
impl<T: CArg + ?Sized> CArg for &T {
    fn to_arg(&self) -> Arg<'_> {
        (**self).to_arg()
    }
}
impl<T: CArg + ?Sized> CArg for &mut T {
    fn to_arg(&self) -> Arg<'_> {
        (**self).to_arg()
    }
}

/// Format like C `sprintf(fmt, args...)`.
pub fn format(fmt: &str, args: &[Arg]) -> String {
    let mut out = Vec::with_capacity(fmt.len() + 16);
    format_bytes(&mut out, fmt.as_bytes(), args);
    // Arguments and format are UTF-8; %c may emit a lone byte >= 0x80.
    match String::from_utf8(out) {
        Ok(s) => s,
        Err(e) => String::from_utf8_lossy(e.as_bytes()).into_owned(),
    }
}

/// Format into raw bytes (use when `%c` may produce non-UTF-8 bytes).
pub fn format_bytes(out: &mut Vec<u8>, fmt: &[u8], args: &[Arg]) {
    let mut ai = 0usize;
    let mut next = |position: Option<usize>, what: &str| -> Arg {
        let index = position.unwrap_or_else(|| {
            let index = ai;
            ai += 1;
            index
        });
        if index >= args.len() {
            panic!("printf: too few arguments for {:?} ({})", String::from_utf8_lossy(fmt), what);
        }
        args[index]
    };
    let mut i = 0;
    while i < fmt.len() {
        let c = fmt[i];
        if c != b'%' {
            out.push(c);
            i += 1;
            continue;
        }
        i += 1;
        if i >= fmt.len() {
            break;
        }
        let position = arg_position(fmt, &mut i);
        // flags
        let (mut minus, mut plus, mut space, mut alt, mut zero) = (false, false, false, false, false);
        while i < fmt.len() {
            match fmt[i] {
                b'-' => minus = true,
                b'+' => plus = true,
                b' ' => space = true,
                b'#' => alt = true,
                b'0' => zero = true,
                _ => break,
            }
            i += 1;
        }
        // width
        let mut width: usize = 0;
        if i < fmt.len() && fmt[i] == b'*' {
            i += 1;
            let w = int_of(next(arg_position(fmt, &mut i), "*"), true, 32);
            if w < 0 {
                minus = true;
                width = w.unsigned_abs() as usize;
            } else {
                width = w as usize;
            }
        } else {
            while i < fmt.len() && fmt[i].is_ascii_digit() {
                width = width * 10 + (fmt[i] - b'0') as usize;
                i += 1;
            }
        }
        // precision
        let mut prec: Option<usize> = None;
        if i < fmt.len() && fmt[i] == b'.' {
            i += 1;
            if i < fmt.len() && fmt[i] == b'*' {
                i += 1;
                let p = int_of(next(arg_position(fmt, &mut i), ".*"), true, 32);
                prec = if p < 0 { None } else { Some(p as usize) };
            } else {
                let mut p = 0usize;
                while i < fmt.len() && fmt[i].is_ascii_digit() {
                    p = p * 10 + (fmt[i] - b'0') as usize;
                    i += 1;
                }
                prec = Some(p);
            }
        }
        // length modifier -> integer width in bits
        let mut bits = 32u32;
        while i < fmt.len() {
            match fmt[i] {
                b'h' => {
                    bits = if bits == 16 { 8 } else { 16 };
                }
                b'l' | b'L' | b'q' | b'j' | b'z' | b't' => bits = 64,
                _ => break,
            }
            i += 1;
        }
        if i >= fmt.len() {
            break;
        }
        let conv = fmt[i];
        i += 1;
        let mut body: Vec<u8> = Vec::new(); // digits/text without sign
        let mut sign: &[u8] = b"";
        let mut prefix: &[u8] = b"";
        let mut numeric = true;
        let mut zero_ok = zero && !minus;
        match conv {
            b'%' => {
                out.push(b'%');
                continue;
            }
            b'd' | b'i' => {
                let v = int_of(next(position, "%d"), true, bits);
                if v < 0 {
                    sign = b"-";
                } else if plus {
                    sign = b"+";
                } else if space {
                    sign = b" ";
                }
                body = int_digits(v.unsigned_abs() as u64, 10, false, prec);
                if prec.is_some() {
                    zero_ok = false;
                }
            }
            b'u' | b'o' | b'x' | b'X' => {
                let v = int_of(next(position, "%u"), false, bits) as u64;
                let (base, upper) = match conv {
                    b'o' => (8, false),
                    b'x' => (16, false),
                    b'X' => (16, true),
                    _ => (10, false),
                };
                body = int_digits(v, base, upper, prec);
                if alt {
                    if conv == b'o' && body.first() != Some(&b'0') {
                        body.insert(0, b'0');
                    } else if conv == b'x' && v != 0 {
                        prefix = b"0x";
                    } else if conv == b'X' && v != 0 {
                        prefix = b"0X";
                    }
                }
                if prec.is_some() {
                    zero_ok = false;
                }
            }
            b'c' => {
                let v = int_of(next(position, "%c"), true, 32);
                body.push(v as u8);
                numeric = false;
            }
            b's' => {
                let s = match next(position, "%s") {
                    Arg::Str(s) => s.as_bytes(),
                    other => panic!("printf: %s given {:?}", other),
                };
                let n = prec.map_or(s.len(), |p| p.min(s.len()));
                body.extend_from_slice(&s[..n]);
                numeric = false;
            }
            b'f' | b'F' | b'e' | b'E' | b'g' | b'G' | b'a' | b'A' => {
                let x = match next(position, "%f") {
                    Arg::Float(x) => x,
                    other => panic!("printf: %{} given {:?}", conv as char, other),
                };
                let upper = conv.is_ascii_uppercase();
                if x.is_nan() {
                    // macOS libc prints "nan" with no sign, whatever the flags
                } else if x.is_sign_negative() {
                    sign = b"-";
                } else if plus {
                    sign = b"+";
                } else if space {
                    sign = b" ";
                }
                if !x.is_finite() {
                    let t: &[u8] = match (x.is_nan(), upper) {
                        (true, false) => b"nan",
                        (true, true) => b"NAN",
                        (false, false) => b"inf",
                        (false, true) => b"INF",
                    };
                    body.extend_from_slice(t);
                    zero_ok = false;
                } else {
                    if conv == b'a' || conv == b'A' {
                        prefix = if upper { b"0X" } else { b"0x" };
                        body = hex_float_body(x.abs(), prec, alt);
                    } else {
                        let p = prec.unwrap_or(6);
                        body = float_body(x.abs(), conv.to_ascii_lowercase(), p, alt);
                    }
                    if upper {
                        body.make_ascii_uppercase();
                    }
                }
            }
            b'p' => {
                let v = int_of(next(position, "%p"), false, 64) as u64;
                prefix = b"0x";
                body = int_digits(v, 16, false, None);
            }
            other => panic!("printf: unsupported conversion %{}", other as char),
        }
        let len = sign.len() + prefix.len() + body.len();
        let pad = width.saturating_sub(len);
        if minus {
            out.extend_from_slice(sign);
            out.extend_from_slice(prefix);
            out.extend_from_slice(&body);
            out.extend(std::iter::repeat(b' ').take(pad));
        } else if zero_ok && numeric {
            out.extend_from_slice(sign);
            out.extend_from_slice(prefix);
            out.extend(std::iter::repeat(b'0').take(pad));
            out.extend_from_slice(&body);
        } else {
            out.extend(std::iter::repeat(b' ').take(pad));
            out.extend_from_slice(sign);
            out.extend_from_slice(prefix);
            out.extend_from_slice(&body);
        }
    }
}

/// POSIX argument number, present only when the digits are followed by `$`.
/// Widths such as `%12g` must leave the parser at the first digit.
fn arg_position(fmt: &[u8], i: &mut usize) -> Option<usize> {
    let mut end = *i;
    let mut n = 0usize;
    while end < fmt.len() && fmt[end].is_ascii_digit() {
        n = n.saturating_mul(10).saturating_add((fmt[end] - b'0') as usize);
        end += 1;
    }
    if end > *i && fmt.get(end) == Some(&b'$') && n > 0 {
        *i = end + 1;
        Some(n - 1)
    } else {
        None
    }
}

/// Interpret an integer argument as C would after the default promotions and
/// the length modifier (`bits` = 8, 16, 32 or 64).
fn int_of(a: Arg, signed: bool, bits: u32) -> i128 {
    let raw = match a {
        Arg::Int(v) => v,
        other => panic!("printf: integer conversion given {:?}", other),
    };
    match (signed, bits) {
        (true, 8) => raw as i8 as i128,
        (true, 16) => raw as i16 as i128,
        (true, 32) => raw as i32 as i128,
        (true, _) => raw as i64 as i128,
        (false, 8) => raw as u8 as i128,
        (false, 16) => raw as u16 as i128,
        (false, 32) => raw as u32 as i128,
        (false, _) => raw as i128,
    }
}

fn int_digits(v: u64, base: u64, upper: bool, prec: Option<usize>) -> Vec<u8> {
    let digits: &[u8] = if upper { b"0123456789ABCDEF" } else { b"0123456789abcdef" };
    let mut d = Vec::new();
    let mut x = v;
    while x > 0 {
        d.push(digits[(x % base) as usize]);
        x /= base;
    }
    let p = prec.unwrap_or(1);
    while d.len() < p {
        d.push(b'0');
    }
    d.reverse();
    d
}

/// Digits for %f / %e / %g of a finite non-negative `x` (lowercase).
fn float_body(x: f64, conv: u8, p: usize, alt: bool) -> Vec<u8> {
    match conv {
        b'f' => {
            let mut s = format!("{:.*}", p, x);
            if alt && p == 0 {
                s.push('.');
            }
            s.into_bytes()
        }
        b'e' => exp_form(x, p, alt, false),
        _ => {
            // %g
            let p = if p == 0 { 1 } else { p };
            let e = format!("{:.*e}", p - 1, x);
            let xexp: i64 = if x == 0.0 {
                0
            } else {
                e[e.find('e').unwrap() + 1..].parse().unwrap()
            };
            let mut s = if (p as i64) > xexp && xexp >= -4 {
                let fp = (p as i64 - 1 - xexp) as usize;
                let mut s = format!("{:.*}", fp, x);
                if alt && fp == 0 {
                    s.push('.');
                }
                s.into_bytes()
            } else {
                exp_form(x, p - 1, alt, false)
            };
            if !alt && !g_retains_zeros(x, p, xexp) {
                strip_zeros(&mut s);
            }
            s
        }
    }
}

/// Darwin dtoa's small-integer path (decimal exponent at most 14) returns
/// all requested digits on an exact halfway tie rounded down. Its other
/// paths trim zeros, as does rounding up through a run of nines.
fn g_retains_zeros(x: f64, p: usize, exp: i64) -> bool {
    if !(1..=14).contains(&exp) || p > exp as usize {
        return false;
    }
    let scale = 10u64.pow((exp as usize + 1 - p) as u32);
    let n = x as u64;
    x == n as f64 && n % scale == scale / 2 && (n / scale) % 10 == 0
}

/// Darwin %a normalizes even subnormal doubles to a leading 1. Rounding
/// can produce a leading 2 without changing the binary exponent.
fn hex_float_body(x: f64, prec: Option<usize>, alt: bool) -> Vec<u8> {
    let bits = x.to_bits();
    let fraction = bits & ((1u64 << 52) - 1);
    let exponent = (bits >> 52) as i32;
    let (mut mantissa, exp) = if exponent != 0 {
        (fraction | (1u64 << 52), exponent - 1023)
    } else if fraction != 0 {
        let shift = fraction.leading_zeros() - 11;
        (fraction << shift, -1022 - shift as i32)
    } else {
        (0, 0)
    };
    let p = prec.unwrap_or_else(|| {
        let mut p = 13;
        let mut m = mantissa;
        while p > 0 && m & 15 == 0 {
            p -= 1;
            m >>= 4;
        }
        p
    });
    if p < 13 {
        let shift = 52 - 4 * p;
        let mut rounded = mantissa >> shift;
        let guard = (mantissa >> (shift - 4)) & 15;
        // Darwin __hdtoa tests the low bit of the *second discarded*
        // nibble when the first is 8, rather than the retained digit or
        // all remaining bits. Preserve this libc behaviour. At p=12 that
        // second nibble is beyond its initialized digit buffer; use the
        // defined nearest-even rule instead of reproducing that read.
        let odd = if shift >= 8 {
            (mantissa >> (shift - 8)) & 1 != 0
        } else {
            rounded & 1 != 0
        };
        if guard > 8 || (guard == 8 && odd) {
            rounded += 1;
        }
        mantissa = rounded << shift;
    }
    let digits = b"0123456789abcdef";
    let mut s = String::new();
    s.push(digits[(mantissa >> 52) as usize] as char);
    if p > 0 || alt {
        s.push('.');
    }
    for j in 0..p.min(13) {
        s.push(digits[((mantissa >> (48 - 4 * j)) & 15) as usize] as char);
    }
    for _ in 13..p {
        s.push('0');
    }
    let _ = write!(s, "p{:+}", exp);
    s.into_bytes()
}

fn exp_form(x: f64, p: usize, alt: bool, _upper: bool) -> Vec<u8> {
    let e = format!("{:.*e}", p, x);
    let epos = e.find('e').unwrap();
    let mant = &e[..epos];
    let ex: i64 = e[epos + 1..].parse().unwrap();
    let mut s = String::with_capacity(mant.len() + 6);
    s.push_str(mant);
    if alt && p == 0 {
        s.push('.');
    }
    s.push('e');
    s.push(if ex < 0 { '-' } else { '+' });
    let _ = write!(s, "{:02}", ex.abs());
    s.into_bytes()
}

/// %g: remove trailing zeros of the fraction, and the point if no fraction.
fn strip_zeros(s: &mut Vec<u8>) {
    let epos = s.iter().position(|&c| c == b'e').unwrap_or(s.len());
    let (mant, rest) = s.split_at(epos);
    if !mant.contains(&b'.') {
        return;
    }
    let mut m = mant.to_vec();
    while m.last() == Some(&b'0') {
        m.pop();
    }
    if m.last() == Some(&b'.') {
        m.pop();
    }
    m.extend_from_slice(rest);
    *s = m;
}

/// `sprintf!(fmt, args..) -> String` with C printf semantics.
#[macro_export]
macro_rules! sprintf {
    ($fmt:expr $(, $arg:expr)* $(,)?) => {
        $crate::cfmt::format($fmt, &[$($crate::cfmt::CArg::to_arg(&$arg)),*])
    };
}

/// `printf!(fmt, args..)`: C printf to the (buffered) global stdout.
#[macro_export]
macro_rules! printf {
    ($fmt:expr $(, $arg:expr)* $(,)?) => {{
        let __s = $crate::cfmt::format($fmt, &[$($crate::cfmt::CArg::to_arg(&$arg)),*]);
        $crate::cio::fputs(&__s, &$crate::cio::stdout());
    }};
}

/// `fprintf!(fp, fmt, args..)`: C fprintf to a `CFile` (e.g. `&stderr()`).
#[macro_export]
macro_rules! fprintf {
    ($fp:expr, $fmt:expr $(, $arg:expr)* $(,)?) => {{
        let __s = $crate::cfmt::format($fmt, &[$($crate::cfmt::CArg::to_arg(&$arg)),*]);
        $crate::cio::fputs(&__s, &$fp);
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::{c_char, c_int, CStr, CString};

    extern "C" {
        fn snprintf(buf: *mut c_char, n: usize, fmt: *const c_char, ...) -> c_int;
    }

    fn c_f(fmt: &str, x: f64) -> String {
        let mut buf = vec![0 as c_char; 2048];
        let f = CString::new(fmt).unwrap();
        unsafe { snprintf(buf.as_mut_ptr(), buf.len(), f.as_ptr(), x) };
        unsafe { CStr::from_ptr(buf.as_ptr()) }.to_str().unwrap().to_string()
    }
    fn c_i(fmt: &str, x: i64) -> String {
        let mut buf = vec![0 as c_char; 512];
        let f = CString::new(fmt).unwrap();
        unsafe { snprintf(buf.as_mut_ptr(), buf.len(), f.as_ptr(), x) };
        unsafe { CStr::from_ptr(buf.as_ptr()) }.to_str().unwrap().to_string()
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    const FFMTS: &[&str] = &[
        "%f", "%e", "%g", "%G", "%E", "%.0f", "%.1f", "%.3f", "%.10f", "%12.4f", "%-12.4f|",
        "%+f", "% e", "%#g", "%#.0f", "%#.0e", "%.0e", "%.0g", "%.1g", "%.3g", "%.10g", "%.17g",
        "%015.3e", "%-+15.3g|", "%010f", "%.20e", "%#x-%g", "%11.4e", "%lf", "%10.2g",
    ];

    #[test]
    fn floats_match_libc() {
        let mut r = Rng(0x9E3779B97F4A7C15);
        let mut vals: Vec<f64> = vec![
            0.0, -0.0, 1.0, -1.0, 0.5, 1.5, 2.5, 0.125, 0.25, 1e-5, 1e-4, 9.9999995e-5, 123456.0,
            1234567.0, 999999.5, 9.5, 0.05, 1e300, -1e-300, 5e-324, f64::MAX, f64::MIN_POSITIVE,
            f64::INFINITY, -f64::INFINITY, f64::NAN, -f64::NAN, 3.0f32 as f64, 0.1f32 as f64,
            1e22, 1e23, 4.35, 0.004, 2.0e-3,
        ];
        for _ in 0..3000 {
            vals.push(f64::from_bits(r.next()));
            vals.push(f32::from_bits(r.next() as u32) as f64);
            let m = (r.next() % 100000) as f64 / 1000.0;
            vals.push(m);
            vals.push(-m / 7.0);
        }
        for &x in &vals {
            for f in FFMTS {
                if f.contains("%#x") {
                    continue;
                }
                let want = c_f(f, x);
                let got = format(f, &[Arg::Float(x)]);
                assert_eq!(got, want, "fmt {} x {:e} bits {:x}", f, x, x.to_bits());
            }
        }
    }

    #[test]
    fn hex_floats_match_libc() {
        assert_eq!(sprintf!("%a|%A|%#a", 1.5, 0.125, 0.0), "0x1.8p+0|0X1P-3|0x0.p+0");
        assert_eq!(sprintf!("%a|%.0a|%.1a", f64::from_bits(1), 1.5, 1.96875),
            "0x1p-1074|0x1p+0|0x1.fp+0");
        assert_eq!(sprintf!("%.1a|%.1a", 1.0 + 129.0 / 4096.0, 1.0 + 136.0 / 4096.0),
            "0x1.1p+0|0x1.0p+0");
        let fmts = [
            "%a", "%A", "%#a", "%+a", "% a", "%020a", "%+025.3A", "%-#25.0a|",
            "%.0a", "%.1a", "%.2a", "%.3a", "%.6a", "%.12a", "%.13a", "%.14a", "%.20a",
        ];
        let mut vals = vec![
            0.0, -0.0, 1.0, -1.0, 1.5, 2.5, 1.96875, 1.90625, 1.03125, 1.09375,
            f64::from_bits(1), f64::from_bits(15), f64::from_bits((1u64 << 52) - 1),
            f64::MIN_POSITIVE, f64::MAX, f64::INFINITY, -f64::INFINITY, f64::NAN, -f64::NAN,
        ];
        // Exercise subnormal normalization and exact halfway hexadecimal ties
        // at every precision that discards any mantissa bits.
        for bit in 0..52 {
            vals.push(f64::from_bits(1u64 << bit));
        }
        for p in 0..13 {
            let half = 1u64 << (51 - 4 * p);
            for n in [0, 1, 2, 14, 15] {
                let bits = 1.0f64.to_bits() + n * (half * 2) + half;
                for bits in [bits - 1, bits, bits + 1] {
                    vals.push(f64::from_bits(bits));
                }
            }
        }
        let mut r = Rng(0xd1b54a32d192ed03);
        for _ in 0..3000 {
            vals.push(f64::from_bits(r.next()));
        }
        for x in vals {
            let bits = x.abs().to_bits();
            let mantissa = if bits != 0 && bits >> 52 == 0 {
                bits << (bits.leading_zeros() - 11)
            } else {
                bits
            };
            for f in fmts {
                // __hdtoa reads an uninitialized byte for this final-nibble
                // tie. Its answer depends on the allocator's previous use.
                if f == "%.12a" && x.is_finite() && mantissa & 15 == 8 {
                    continue;
                }
                assert_eq!(format(f, &[Arg::Float(x)]), c_f(f, x),
                    "fmt {} bits {:x}", f, x.to_bits());
            }
        }
    }

    #[test]
    fn hex_final_nibble_ties_are_deterministic() {
        let even = f64::from_bits(1.0f64.to_bits() + 8);
        let odd = f64::from_bits(1.0f64.to_bits() + 24);
        assert_eq!(sprintf!("%.12a|%.12a", even, odd),
            "0x1.000000000000p+0|0x1.000000000002p+0");
        assert_eq!(sprintf!("%.12A|%.12A", -even, -odd),
            "-0X1.000000000000P+0|-0X1.000000000002P+0");
    }

    #[test]
    fn positional_arguments_match_libc() {
        assert_eq!(sprintf!("%1$.6g|%1$A", 12.375), "12.375|0X1.8CP+3");
        assert_eq!(sprintf!("%3$s|%1$04d|%2$hhd|%3$.2s", 7, 255, "abc"), "abc|0007|-1|ab");
        for x in [0.0, 1000005.0, -0.125, f64::from_bits(1), f64::INFINITY] {
            assert_eq!(sprintf!("%1$.6g", x), c_f("%1$.6g", x));
            for w in [-24, 0, 24] {
                for p in [-1, 0, 3, 14] {
                    for f in ["%3$*1$.*2$f|%3$g", "%3$#*1$.*2$A", "%2$d|%3$*1$g"] {
                        let cf = CString::new(f).unwrap();
                        let mut buf = vec![0 as c_char; 2048];
                        unsafe { snprintf(buf.as_mut_ptr(), buf.len(), cf.as_ptr(), w, p, x) };
                        let want = unsafe { CStr::from_ptr(buf.as_ptr()) }.to_str().unwrap();
                        assert_eq!(sprintf!(f, w, p, x), want, "fmt {} w {} p {} x {}", f, w, p, x);
                    }
                }
            }
        }
        // Sequential stars still consume arguments before the conversion.
        assert_eq!(sprintf!("%*.*a|%d", 12, 2, 1.5, 7), "   0x1.80p+0|7");
    }

    #[test]
    fn diagnostic_strings_match_libc() {
        // These are the only conversions reached before elamodel's model
        // construction. Keep their sequential argument handling unchanged
        // when extending the positional and floating conversion parser.
        let program = CString::new("elamodel").unwrap();
        for argument in ["fill=0.1,0.01,2,1,2", "fill=0.1,0.99,3,1,3"] {
            let arg = CString::new(argument).unwrap();
            let fmt = CString::new("Unknown %s argument %s\n").unwrap();
            let mut buf = vec![0 as c_char; 256];
            unsafe {
                snprintf(buf.as_mut_ptr(), buf.len(), fmt.as_ptr(), program.as_ptr(), arg.as_ptr());
            }
            let want = unsafe { CStr::from_ptr(buf.as_ptr()) }.to_str().unwrap();
            assert_eq!(sprintf!("Unknown %s argument %s\n", "elamodel", argument), want);
        }
        let fmt = CString::new("getpar() call after checkpars(): %s\n").unwrap();
        let name = CString::new("fill").unwrap();
        let mut buf = vec![0 as c_char; 256];
        unsafe { snprintf(buf.as_mut_ptr(), buf.len(), fmt.as_ptr(), name.as_ptr()) };
        let want = unsafe { CStr::from_ptr(buf.as_ptr()) }.to_str().unwrap();
        assert_eq!(sprintf!("getpar() call after checkpars(): %s\n", "fill"), want);
    }

    #[test]
    fn g_halfway_zeros_match_darwin() {
        assert_eq!(sprintf!("%g|%g|%g", 1000005.0, 6217305000000.0, -2930105.0),
            "1.00000e+06|6.21730e+12|-2.93010e+06");
        // The small-integer path stops at exponent 14. Exact numbers, nearby
        // non-ties, and ties rounded up must still have their zeros stripped.
        for exp in 2..=16u32 {
            for p in 2..=exp.min(15) {
                let scale = 10u64.pow(exp + 1 - p);
                for leading in [10u64.pow(p - 1), 12 * 10u64.pow(p - 2)] {
                    for last in [0, 5, 9] {
                        let x = ((leading + last) * scale) as f64 + (scale / 2) as f64;
                        for x in [x, f64::from_bits(x.to_bits() - 1), f64::from_bits(x.to_bits() + 1), -x] {
                            for f in [format!("%.{}g", p), format!("%#.{}g", p), format!("%+025.{}G", p)] {
                                assert_eq!(format(&f, &[Arg::Float(x)]), c_f(&f, x),
                                    "fmt {} x {} bits {:x}", f, x, x.to_bits());
                            }
                        }
                    }
                }
            }
        }
        for x in [1000000.0, 999999.5, 1000095.0, 100000.5, 1.25, 0.0001000005] {
            assert_eq!(sprintf!("%g", x), c_f("%g", x));
        }
    }

    #[test]
    fn ints_match_libc() {
        let fmts = [
            "%d", "%5d", "%-5d|", "%05d", "%+d", "% d", "%.3d", "%8.3d", "%u", "%x", "%X", "%#x",
            "%#o", "%o", "%ld", "%lu", "%lx", "%hd", "%hu", "%hhd", "%c", "%.0d", "%#.0o", "%i",
            "%-08d|",
        ];
        let mut r = Rng(12345);
        let mut vals: Vec<i64> = vec![0, 1, -1, 65, 255, 32767, 32768, -32768, 70000, i32::MAX as i64,
            i32::MIN as i64, i64::MAX, i64::MIN, 4294967295];
        for _ in 0..2000 {
            vals.push(r.next() as i64);
            vals.push((r.next() as i64) % 100000);
        }
        for &v in &vals {
            for f in fmts {
                if f == "%c" && !(32..127).contains(&(v as i32)) {
                    continue;
                }
                assert_eq!(format(f, &[Arg::Int(v as u64)]), c_i(f, v), "fmt {} v {}", f, v);
            }
        }
    }

    #[test]
    fn strings_and_misc() {
        assert_eq!(sprintf!("%6s=", "tracl"), " tracl=");
        assert_eq!(sprintf!("%-4s|%.2s|%%", "ab", "xyz"), "ab  |xy|%");
        assert_eq!(sprintf!("%*d|%-*d|%.*f", 5, 42, 4, 7, 2, 3.14159f32), "   42|7   |3.14");
        assert_eq!(sprintf!("%d %s %g", 3i16, String::from("x"), 0.1f32), "3 x 0.1");
        assert_eq!(sprintf!("%u", -1i32), "4294967295");
        assert_eq!(sprintf!("%ld", -1i32), "-1");
    }
}
