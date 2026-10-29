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
//! strings verbatim (`%d`, `%6s`, `%g`, `%-10.3e`, `%ld`, `%hd`, `*`, ...).
//! Output is identical to macOS libc `printf` (tested against `snprintf`).
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
    let mut next = |what: &str| -> Arg {
        if ai >= args.len() {
            panic!("printf: too few arguments for {:?} ({})", String::from_utf8_lossy(fmt), what);
        }
        ai += 1;
        args[ai - 1]
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
            let w = int_of(next("*"), true, 32);
            if w < 0 {
                minus = true;
                width = w.unsigned_abs() as usize;
            } else {
                width = w as usize;
            }
            i += 1;
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
                let p = int_of(next(".*"), true, 32);
                prec = if p < 0 { None } else { Some(p as usize) };
                i += 1;
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
                let v = int_of(next("%d"), true, bits);
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
                let v = int_of(next("%u"), false, bits) as u64;
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
                let v = int_of(next("%c"), true, 32);
                body.push(v as u8);
                numeric = false;
            }
            b's' => {
                let s = match next("%s") {
                    Arg::Str(s) => s.as_bytes(),
                    other => panic!("printf: %s given {:?}", other),
                };
                let n = prec.map_or(s.len(), |p| p.min(s.len()));
                body.extend_from_slice(&s[..n]);
                numeric = false;
            }
            b'f' | b'F' | b'e' | b'E' | b'g' | b'G' => {
                let x = match next("%f") {
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
                    let p = prec.unwrap_or(6);
                    body = float_body(x.abs(), conv.to_ascii_lowercase(), p, alt);
                    if upper {
                        body.make_ascii_uppercase();
                    }
                }
            }
            b'p' => {
                let v = int_of(next("%p"), false, 64) as u64;
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
            if !alt {
                strip_zeros(&mut s);
            }
            s
        }
    }
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
