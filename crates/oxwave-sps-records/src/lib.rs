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

//! Port of su/lib/getSPSfield.c.
//! Author: Fernando Roxo <petro@roxo.org>.

use std::io::Write;

/// Destination for the C field's `%c`, `%s`, `%d`, `%f`, or `%lf` conversion.
/// String destinations receive the terminating NUL as in C.
pub enum SPSField<'a> {
    Char(&'a mut u8),
    String(&'a mut [u8]),
    Int(&'a mut i32),
    Float(&'a mut f32),
    Double(&'a mut f64),
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Char,
    String,
    Int,
    Float,
    Double,
}
use Kind::*;
const P0: [(usize, usize, Kind); 15] = [
    (0, 1, Char),
    (1, 16, Float),
    (17, 8, Float),
    (25, 1, Int),
    (26, 2, Int),
    (28, 4, Int),
    (32, 4, Float),
    (36, 4, Int),
    (40, 2, Int),
    (42, 4, Float),
    (46, 9, Double),
    (55, 10, Double),
    (65, 6, Float),
    (71, 3, Int),
    (74, 6, Int),
];
const P2: [(usize, usize, Kind); 15] = [
    (0, 1, Char),
    (1, 10, Float),
    (11, 10, Float),
    (23, 1, Int),
    (24, 2, Int),
    (26, 4, Int),
    (30, 4, Float),
    (34, 4, Int),
    (38, 2, Int),
    (40, 6, Float),
    (46, 9, Double),
    (55, 10, Double),
    (65, 6, Float),
    (71, 3, Int),
    (74, 6, Int),
];
const X0: [(usize, usize, Kind); 15] = [
    (0, 1, Char),
    (1, 6, String),
    (7, 4, Int),
    (11, 1, Int),
    (12, 1, Char),
    (13, 16, Float),
    (29, 8, Float),
    (37, 1, Int),
    (38, 4, Int),
    (42, 4, Int),
    (46, 1, Int),
    (47, 16, Float),
    (63, 8, Float),
    (71, 8, Float),
    (79, 1, Int),
];
const X2: [(usize, usize, Kind); 15] = [
    (0, 1, Char),
    (1, 6, String),
    (7, 8, Int),
    (15, 1, Int),
    (16, 1, Char),
    (17, 10, Float),
    (27, 10, Float),
    (37, 1, Int),
    (38, 5, Int),
    (43, 5, Int),
    (48, 1, Int),
    (49, 10, Float),
    (59, 10, Float),
    (69, 10, Float),
    (79, 1, Int),
];

/// Extract one SPS record field, returning the C `sscanf` result (1, 0, or -1).
///
/// C prototype: `int getSPSfield(char *Line, const int rev, int const idx, void *addr)`.
/// Input ends at the first NUL or at the end of the slice. `idx` must be
/// 0..15, and `addr` must match the field's C conversion type. A string
/// destination must fit the converted token plus its terminating NUL.
/// Numeric parsing and whitespace follow the C locale. Unsupported revisions
/// or record labels print the original diagnostic, leave `addr` untouched,
/// and return 1, as the C implementation does.
#[allow(non_snake_case)]
pub fn getSPSfield(line: &[u8], rev: i32, idx: i32, addr: SPSField<'_>) -> i32 {
    if rev != 0 && rev != 2 {
        eprintln!("****  Unsupported SPS format {}!!\n", rev);
        return 1;
    }
    let label = line.first().copied().unwrap_or(0);
    let table = match (label, rev) {
        (b'R' | b'S', 0) => &P0,
        (b'R' | b'S', 2) => &P2,
        (b'X', 0) => &X0,
        (b'X', 2) => &X2,
        _ => {
            let _ = std::io::stderr().lock().write_all(
                &[
                    b"****  Unsupported SPS record ".as_slice(),
                    &[label],
                    b"!!\n\n",
                ]
                .concat(),
            );
            return 1;
        }
    };
    let (start, len, kind) = table[usize::try_from(idx).expect("negative SPS field index")];
    let line_len = line.iter().position(|&c| c == 0).unwrap_or(line.len());
    if line_len < start {
        return -1;
    }
    let end = (start + len).min(line_len);
    let field = &line[start..end];
    let first = field.iter().position(|&c| !scanf_space(c));
    match addr {
        SPSField::Char(v) if kind == Char => {
            if let Some(&c) = field.first() {
                *v = c;
                1
            } else {
                -1
            }
        }
        SPSField::String(v) if kind == String => {
            let Some(start) = first else { return -1 };
            let token = &field[start..];
            let n = token
                .iter()
                .position(|&c| scanf_space(c))
                .unwrap_or(token.len());
            assert!(v.len() > n, "SPS string destination is too small");
            v[..n].copy_from_slice(&token[..n]);
            v[n] = 0;
            1
        }
        SPSField::Int(v) if kind == Int => {
            let Some(start) = first else { return -1 };
            let token = &field[start..];
            let mut i = usize::from(matches!(token.first(), Some(b'+' | b'-')));
            let digits = i;
            let mut value = 0i64;
            while i < token.len() && token[i].is_ascii_digit() {
                value = value
                    .saturating_mul(10)
                    .saturating_add((token[i] - b'0') as i64);
                i += 1;
            }
            if i == digits {
                return 0;
            }
            if token[0] == b'-' {
                value = -value;
            }
            *v = value as i32;
            1
        }
        SPSField::Float(v) if kind == Float => {
            let Some(start) = first else { return -1 };
            if let Some(x) = special_float(&field[start..], true) {
                *v = x as f32;
                return 1;
            }
            match parse_float(&field[start..]).and_then(|s| s.parse::<f32>().ok()) {
                Some(x) => {
                    *v = x;
                    1
                }
                None => 0,
            }
        }
        SPSField::Double(v) if kind == Double => {
            let Some(start) = first else { return -1 };
            if let Some(x) = special_float(&field[start..], false) {
                *v = x;
                return 1;
            }
            match parse_float(&field[start..]).and_then(|s| s.parse::<f64>().ok()) {
                Some(x) => {
                    *v = x;
                    1
                }
                None => 0,
            }
        }
        _ => panic!("SPS field destination has the wrong type"),
    }
}

// C's isspace includes vertical tab, unlike u8::is_ascii_whitespace.
fn scanf_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

fn special_float(token: &[u8], single: bool) -> Option<f64> {
    let negative = token.first() == Some(&b'-');
    let sign_len = usize::from(matches!(token.first(), Some(b'+' | b'-')));
    let s = &token[sign_len..];
    let sign = if negative { -1.0 } else { 1.0 };
    if s.len() >= 3 && s[..3].eq_ignore_ascii_case(b"inf") {
        return Some(sign * f64::INFINITY);
    }
    if s.len() >= 3 && s[..3].eq_ignore_ascii_case(b"nan") {
        return Some(sign * f64::NAN);
    }
    if s.len() < 2 || s[0] != b'0' || !matches!(s[1], b'x' | b'X') {
        return None;
    }
    // The longest field is 16 bytes, so at most 14 hexadecimal digits
    // follow 0x. Preserve them exactly until rounding to the target type.
    let mut significand = 0u64;
    let mut i = 2;
    let mut fraction = false;
    let mut digits = 0;
    let mut fractional_digits = 0;
    while i < s.len() {
        if s[i] == b'.' && !fraction {
            fraction = true;
            i += 1;
            continue;
        }
        let digit = (s[i] as char).to_digit(16);
        let Some(digit) = digit else { break };
        significand = significand * 16 + digit as u64;
        digits += 1;
        if fraction {
            fractional_digits += 1;
        }
        i += 1;
    }
    if digits == 0 {
        return None;
    }
    let mut exponent: i32 = -4 * fractional_digits;
    if matches!(s.get(i), Some(b'p' | b'P')) {
        i += 1;
        let neg_exp = s.get(i) == Some(&b'-');
        if matches!(s.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let mut e: i32 = 0;
        while i < s.len() && s[i].is_ascii_digit() {
            e = e.saturating_mul(10).saturating_add((s[i] - b'0') as i32);
            i += 1;
        }
        exponent = exponent.saturating_add(if neg_exp { -e } else { e });
    }
    Some(hex_float(significand, exponent, negative, single))
}

// Convert significand * 2^exponent directly to the target IEEE format.
// Rounding through f64 first can double-round a hexadecimal %f field;
// powi also underflows prematurely at f64's smallest subnormal exponents.
fn hex_float(significand: u64, exponent: i32, negative: bool, single: bool) -> f64 {
    let (precision, min_exp, max_exp) = if single {
        (24, -126, 127)
    } else {
        (53, -1022, 1023)
    };
    let infinity = ((2 * max_exp + 1) as u64) << (precision - 1);
    let mut bits = 0;
    if significand != 0 {
        let nbits = 64 - significand.leading_zeros() as i32;
        let mut e = exponent.saturating_add(nbits - 1);
        if e > max_exp {
            bits = infinity;
        } else if e < min_exp {
            // Subnormals have fixed spacing; rounding can reach the smallest
            // normal, whose encoding is contiguous with the subnormals.
            bits = round_hex(
                significand,
                (min_exp - precision + 1).saturating_sub(exponent),
            );
        } else {
            let mut rounded = round_hex(significand, nbits - precision);
            if rounded == 1u64 << precision {
                rounded >>= 1;
                e += 1;
            }
            bits = if e > max_exp {
                infinity
            } else {
                (((e + max_exp) as u64) << (precision - 1))
                    | (rounded & ((1u64 << (precision - 1)) - 1))
            };
        }
    }
    if single {
        f32::from_bits(bits as u32 | ((negative as u32) << 31)) as f64
    } else {
        f64::from_bits(bits | ((negative as u64) << 63))
    }
}

fn round_hex(significand: u64, shift: i32) -> u64 {
    if shift <= 0 {
        return significand << -shift;
    }
    // The input has at most 56 significant bits, so this rounds to zero.
    if shift >= 64 {
        return 0;
    }
    let rounded = significand >> shift;
    let remainder = significand & ((1u64 << shift) - 1);
    let half = 1u64 << (shift - 1);
    rounded + u64::from(remainder > half || (remainder == half && rounded & 1 != 0))
}

// Identify the initial decimal floating-point subject sequence accepted by scanf.
fn parse_float(token: &[u8]) -> Option<&str> {
    let mut i = usize::from(matches!(token.first(), Some(b'+' | b'-')));
    let start = i;
    while i < token.len() && token[i].is_ascii_digit() {
        i += 1;
    }
    if token.get(i) == Some(&b'.') {
        i += 1;
        while i < token.len() && token[i].is_ascii_digit() {
            i += 1;
        }
    }
    if !token[start..i].iter().any(u8::is_ascii_digit) {
        return None;
    }
    let mantissa_end = i;
    if matches!(token.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(token.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let exp_start = i;
        while i < token.len() && token[i].is_ascii_digit() {
            i += 1;
        }
        if i == exp_start {
            i = mantissa_end;
        }
    }
    std::str::from_utf8(&token[..i]).ok()
}
