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

//! `par`: command-line parameters (par/lib/getpars.c) and checked ascii
//! conversions (par/lib/atopkge.c), with SU's exact semantics.
//!
//! | C                                   | Rust                                              |
//! |-------------------------------------|---------------------------------------------------|
//! | `main(int argc, char **argv)`       | `let argv = args();` (`Vec<String>`, argv[0] as invoked) |
//! | `initargs(argc, argv)`              | `initargs(&argv)`                                 |
//! | `xargc`, `xargv`                    | `xargc()`, `xargv()`                              |
//! | `getparint("n", &n)`                | `getparint("n", &mut n)` -> i32 (number of values) |
//! | `getparfloat("v", v)` (array)       | `getparfloat("v", &mut v[..])` or `&mut v` (Vec)  |
//! | `getparuint/short/ushort/long/ulong/double` | same names (u32/i16/u16/i64/u64/f64)      |
//! | `getparstring("s", &s)`             | `getparstring("s", &mut s)` (`s: String`)         |
//! | `getparstringarray("k", k)`         | `getparstringarray("k", &mut k)` (`k: Vec<String>`) |
//! | `getnparint(n, "x", &x)` etc.       | `getnparint(n, "x", &mut x)` etc.                 |
//! | `countparname/countparval/countnparval` | same names, return i32                        |
//! | `checkpars()`                       | `checkpars()`                                     |
//! | `MUSTGETPARINT("n", &n)`            | `mustgetparint("n", &mut n)` (also float/double/string) |
//! | `eatoi(s)` etc. (`eatoh/u/i/p/l/v/f/d`) | same names, `&str` in                         |
//! | `strtol(s,NULL,10)`, `strtoul`, `strtod`, `atoi`, `atof` | `strtol(s).0`, `strtoul(s).0`, `strtod(s).0`, `atoi(s)`, `atof(s)` |
//!
//! Semantics (all verified against the C code):
//! * Numeric getpars return the number of comma-separated values in the
//!   chosen occurrence and store them in order. A scalar destination keeps
//!   only the first; a slice keeps what fits; a `Vec` is overwritten from
//!   index 0 and grows if needed (never truncated).
//! * `n == 0` (getparX) means the LAST occurrence; `getnparX(n>=1)` the n'th,
//!   counting `par=` file entries first, then the command line.
//! * `x=` (empty) returns 0 for numeric types but 1 (empty string) for strings.
//!   `"1,,2"` -> [1,0,2]; `"1,"` -> [1]; leading garbage converts to 0 with
//!   no error, trailing garbage is ignored (C `strtol`/`strtod`); integers are
//!   base 10 only; floats accept `inf`, `nan`, hex floats.
//! * Range errors exit like SU: `atopkge.c: eatoi: overflow`, ... including
//!   macOS `strtod` underflow (e.g. `1e-310` -> `eatod: overflow`).
//! * `checkpars()` warns `Unknown <argv0> argument <arg>` for every `name=`
//!   argument whose name was never asked for, prints
//!   `Using Big Endian SU data format w/ XDR.` when `verbose=1`, and any
//!   later getpar prints `getpar() call after checkpars(): <name>`.

use std::cell::{Cell, RefCell};

const PAR_NAMES_MAX: usize = 512;

#[derive(Default)]
struct ParState {
    xargv: Vec<String>,
    tabled: bool,
    argtbl: Vec<(Vec<u8>, Vec<u8>)>,
    par_names: Vec<String>,
    parcheck: bool,
}

thread_local! {
    static PAR: RefCell<ParState> = RefCell::new(ParState::default());
    static ERRNO: Cell<i32> = const { Cell::new(0) };
}

/// ERANGE on this platform (macOS).
pub const ERANGE: i32 = 34;
/// EINVAL on this platform.
pub const EINVAL: i32 = 22;

/// The emulated C `errno`, as set by `strtol`/`strtoul`/`strtod` here.
/// Like C, it is never reset by those functions.
pub fn errno() -> i32 {
    ERRNO.with(|e| e.get())
}
pub fn set_errno(e: i32) {
    ERRNO.with(|c| c.set(e))
}

/// The process arguments (argv[0] as invoked, lossily converted to UTF-8).
pub fn args() -> Vec<String> {
    std::env::args_os().map(|a| a.to_string_lossy().into_owned()).collect()
}

/// initargs: make the command line available to getpar (re-entrant).
/// Also restores default SIGPIPE handling like a C program.
pub fn initargs(argv: &[String]) {
    crate::cio::restore_sigpipe();
    PAR.with(|p| {
        let mut p = p.borrow_mut();
        p.par_names = vec!["par".to_string(), "lheader".to_string()];
        p.xargv = argv.to_vec();
        p.tabled = false;
        p.argtbl.clear();
        p.parcheck = false;
    });
}

/// `xargc` global.
pub fn xargc() -> usize {
    PAR.with(|p| p.borrow().xargv.len())
}
/// `xargv` global (a copy).
pub fn xargv() -> Vec<String> {
    PAR.with(|p| p.borrow().xargv.clone())
}

fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}
fn c_isgraph(c: u8) -> bool {
    (0x21..=0x7e).contains(&c)
}

fn getparinit(p: &mut ParState) {
    p.tabled = true;
    if p.xargv.is_empty() {
        drop_and_err_xargc();
    }
    // par file name: last "par=..." with something after "par="
    let mut pfname: Option<String> = None;
    for a in p.xargv.iter().skip(1).rev() {
        if a.starts_with("par=") && a.len() != 4 {
            pfname = Some(a[4..].to_string());
            break;
        }
    }
    let mut tokens: Vec<Vec<u8>> = Vec::new();
    if let Some(name) = pfname {
        let data = match std::fs::read(&name) {
            Ok(d) => d,
            Err(_) => {
                // efopen failure inside getparinit
                PAR_ERR.with(|e| e.set(true));
                return;
            }
        };
        let mut s: Vec<u8> = data.iter().map(|b| b & 0x7F).collect();
        s.push(0);
        let nread = data.len();
        let mut starts = Vec::new();
        let mut start = true;
        let mut quote = false;
        for i in 0..nread {
            if start {
                if c_isgraph(s[i]) {
                    starts.push(i);
                    start = false;
                } else {
                    s[i] = 0;
                }
            } else if !quote && c_isspace(s[i]) {
                s[i] = 0;
                start = true;
            }
            if s[i] == b'\'' || s[i] == b'"' {
                quote = !quote;
            }
        }
        for st in starts {
            let end = s[st..].iter().position(|&c| c == 0).map_or(s.len(), |e| st + e);
            let tok = s[st..end].to_vec();
            if tok.contains(&b'=') {
                tokens.push(tok);
            }
        }
    }
    for a in p.xargv.iter().skip(1) {
        tokens.push(a.as_bytes().to_vec());
    }
    p.argtbl.clear();
    for t in tokens {
        if let Some(eq) = t.iter().position(|&c| c == b'=') {
            p.argtbl.push((t[..eq].to_vec(), t[eq + 1..].to_vec()));
        }
    }
}

thread_local! {
    static PAR_ERR: Cell<bool> = const { Cell::new(false) };
}

fn drop_and_err_xargc() -> ! {
    crate::err!("%s: xargc=%d -- not initiated in main", "getpars.c", 0)
}

/// Ensure the table is built (errors are raised outside the RefCell borrow).
fn ensure_tabled() {
    PAR.with(|p| {
        let mut p = p.borrow_mut();
        if !p.tabled {
            getparinit(&mut p);
        }
    });
    if PAR_ERR.with(|e| e.replace(false)) {
        crate::err!("%s: efopen: fopen failed", "subcalls.c");
    }
}

/// Index of the n'th occurrence (n==0: last) of name, or None.
fn getparindex(p: &ParState, n: i32, name: &str) -> Option<usize> {
    let nb = name.as_bytes();
    if n == 0 {
        (0..p.argtbl.len()).rev().find(|&i| p.argtbl[i].0 == nb)
    } else {
        let mut n = n;
        for i in 0..p.argtbl.len() {
            if p.argtbl[i].0 == nb {
                n -= 1;
                if n == 0 {
                    return Some(i);
                }
            }
        }
        None
    }
}

/// The bookkeeping + lookup half of C getnpar: returns the ascii value of
/// the selected occurrence, or None (C returns 0).
fn getnpar_ascii(n: i32, name: &str) -> Option<Vec<u8>> {
    let (parcheck, full, argv0) = PAR.with(|p| {
        let mut p = p.borrow_mut();
        let parcheck = p.parcheck;
        if !p.par_names.iter().any(|s| s == name) && p.par_names.len() < PAR_NAMES_MAX {
            p.par_names.push(name.to_string());
        }
        let full = p.par_names.len() == PAR_NAMES_MAX;
        (parcheck, full, p.xargv.first().cloned().unwrap_or_default())
    });
    if parcheck && name != "lheader" {
        crate::fprintf!(crate::cio::stderr(), "getpar() call after checkpars(): %s\n", name);
    }
    if full {
        crate::fprintf!(crate::cio::stderr(), " %s exceeded PAR_NAMES_MAX %d \n", argv0, PAR_NAMES_MAX);
    }
    if xargc() == 1 {
        return None;
    }
    ensure_tabled();
    PAR.with(|p| {
        let p = p.borrow();
        getparindex(&p, n, name).map(|i| p.argtbl[i].1.clone())
    })
}

/// Split an ascii value into the C getnpar fields: each field starts at the
/// current position (the converter reads up to the comma by itself).
fn fields(aval: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut pos = 0;
    while pos < aval.len() {
        out.push(&aval[pos..]);
        // while (*aval++ != ',') { if (!*aval) break; }
        loop {
            let c = aval[pos];
            pos += 1;
            if c == b',' {
                break;
            }
            if pos >= aval.len() {
                break;
            }
        }
    }
    out
}

fn b2s(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// Destination of a numeric getpar: a scalar, a slice or a Vec.
pub trait ParDest<T> {
    /// Store value number `k` (0-based) of the occurrence.
    fn put(&mut self, k: usize, v: T);
}

macro_rules! pardest {
    ($($t:ty),*) => {$(
        impl ParDest<$t> for $t {
            fn put(&mut self, k: usize, v: $t) { if k == 0 { *self = v; } }
        }
        impl ParDest<$t> for [$t] {
            fn put(&mut self, k: usize, v: $t) { if k < self.len() { self[k] = v; } }
        }
        impl<const N: usize> ParDest<$t> for [$t; N] {
            fn put(&mut self, k: usize, v: $t) { if k < N { self[k] = v; } }
        }
        impl ParDest<$t> for Vec<$t> {
            fn put(&mut self, k: usize, v: $t) {
                if k < self.len() { self[k] = v; } else { self.push(v); }
            }
        }
    )*};
}
pardest!(i32, u32, i16, u16, i64, u64, f32, f64);

macro_rules! getpar_fns {
    ($($get:ident, $getn:ident, $t:ty, $conv:ident);* $(;)?) => {$(
        #[doc = concat!("getpar for `", stringify!($t), "` (last occurrence). Returns the number of values.")]
        pub fn $get<D: ParDest<$t> + ?Sized>(name: &str, ptr: &mut D) -> i32 {
            $getn(0, name, ptr)
        }
        #[doc = concat!("getnpar for `", stringify!($t), "` (n'th occurrence, 0 = last).")]
        pub fn $getn<D: ParDest<$t> + ?Sized>(n: i32, name: &str, ptr: &mut D) -> i32 {
            let aval = match getnpar_ascii(n, name) { Some(a) => a, None => return 0 };
            let mut nval = 0;
            for f in fields(&aval) {
                let v = $conv(&b2s(f));
                ptr.put(nval, v);
                nval += 1;
            }
            nval as i32
        }
    )*};
}

getpar_fns!(
    getparint, getnparint, i32, eatoi;
    getparuint, getnparuint, u32, eatop;
    getparshort, getnparshort, i16, eatoh;
    getparushort, getnparushort, u16, eatou;
    getparlong, getnparlong, i64, eatol;
    getparulong, getnparulong, u64, eatov;
    getparfloat, getnparfloat, f32, eatof;
    getpardouble, getnpardouble, f64, eatod;
);

/// getparstring: last occurrence as a string (may contain commas). Returns 1/0.
pub fn getparstring(name: &str, ptr: &mut String) -> i32 {
    getnparstring(0, name, ptr)
}
/// getnparstring: n'th occurrence as a string.
pub fn getnparstring(n: i32, name: &str, ptr: &mut String) -> i32 {
    match getnpar_ascii(n, name) {
        Some(a) => {
            *ptr = b2s(&a);
            1
        }
        None => 0,
    }
}

/// getparstringarray: comma-separated strings (C `strchop` per field).
/// Overwrites `ptr[0..nval]`, pushing if the Vec is shorter. Returns nval.
pub fn getparstringarray(name: &str, ptr: &mut Vec<String>) -> i32 {
    getnparstringarray(0, name, ptr)
}
/// getnparstringarray: n'th occurrence of a string array.
pub fn getnparstringarray(n: i32, name: &str, ptr: &mut Vec<String>) -> i32 {
    let aval = match getnpar_ascii(n, name) {
        Some(a) => a,
        None => return 0,
    };
    let mut nval = 0;
    for f in fields(&aval) {
        let end = f.iter().position(|&c| c == b',').unwrap_or(f.len());
        let s = b2s(&f[..end]);
        if nval < ptr.len() {
            ptr[nval] = s;
        } else {
            ptr.push(s);
        }
        nval += 1;
    }
    nval as i32
}

/// countparname: number of occurrences of a parameter name.
pub fn countparname(name: &str) -> i32 {
    if xargc() == 1 {
        return 0;
    }
    ensure_tabled();
    PAR.with(|p| p.borrow().argtbl.iter().filter(|(k, _)| k == name.as_bytes()).count() as i32)
}

/// countnparval: number of values (commas + 1) in the n'th occurrence.
pub fn countnparval(n: i32, name: &str) -> i32 {
    if xargc() == 1 {
        return 0;
    }
    ensure_tabled();
    PAR.with(|p| {
        let p = p.borrow();
        match getparindex(&p, n, name) {
            Some(i) => p.argtbl[i].1.iter().filter(|&&c| c == b',').count() as i32 + 1,
            None => 0,
        }
    })
}

/// countparval: number of values in the last occurrence.
pub fn countparval(name: &str) -> i32 {
    countnparval(0, name)
}

/// checkpars: warn about unknown `name=` arguments (see module docs).
pub fn checkpars() {
    let mut i = 0i32;
    if getparint("verbose", &mut i) != 0 && i == 1 {
        crate::fprintf!(crate::cio::stderr(), "Using Big Endian SU data format w/ XDR.\n");
    }
    let (argv, names) = PAR.with(|p| {
        let p = p.borrow();
        (p.xargv.clone(), p.par_names.clone())
    });
    for a in argv.iter().skip(1) {
        let known = names.iter().any(|n| a.as_bytes().starts_with(format!("{}=", n).as_bytes()));
        if !known && a.contains('=') {
            crate::fprintf!(crate::cio::stderr(), "Unknown %s argument %s\n", argv[0], a);
        }
    }
    PAR.with(|p| p.borrow_mut().parcheck = true);
}

macro_rules! mustget {
    ($($f:ident, $g:ident, $t:ty);*) => {$(
        #[doc = concat!("MUSTGETPAR: `", stringify!($g), "` or err(\"must specify %s=\").")]
        pub fn $f(name: &str, ptr: &mut $t) {
            if $g(name, ptr) == 0 {
                crate::err!("must specify %s=", name);
            }
        }
    )*};
}
mustget!(mustgetparint, getparint, i32; mustgetparfloat, getparfloat, f32;
         mustgetpardouble, getpardouble, f64; mustgetparstring, getparstring, String);

// ------------------------------------------------------------------------
// C conversions: strtol / strtoul / strtod (macOS semantics) and atopkge.c
// ------------------------------------------------------------------------

fn skip_ws(b: &[u8]) -> usize {
    let mut i = 0;
    while i < b.len() && c_isspace(b[i]) {
        i += 1;
    }
    i
}

/// Parse sign and decimal digits. Returns (negative, magnitude-overflowed,
/// magnitude (saturated u64), bytes consumed (0 if no digits)).
fn parse_dec(s: &str) -> (bool, bool, u64, usize) {
    let b = s.as_bytes();
    let mut i = skip_ws(b);
    let mut neg = false;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        neg = b[i] == b'-';
        i += 1;
    }
    let d0 = i;
    let mut mag: u64 = 0;
    let mut ovf = false;
    while i < b.len() && b[i].is_ascii_digit() {
        let d = (b[i] - b'0') as u64;
        match mag.checked_mul(10).and_then(|m| m.checked_add(d)) {
            Some(m) => mag = m,
            None => ovf = true,
        }
        i += 1;
    }
    if i == d0 {
        return (false, false, 0, 0);
    }
    (neg, ovf, mag, i)
}

/// C `strtol(s, &end, 10)` -> (value, bytes consumed). Sets errno like macOS.
pub fn strtol(s: &str) -> (i64, usize) {
    let (neg, ovf, mag, used) = parse_dec(s);
    if used == 0 {
        set_errno(EINVAL);
        return (0, 0);
    }
    if neg {
        if ovf || mag > (i64::MAX as u64) + 1 {
            set_errno(ERANGE);
            (i64::MIN, used)
        } else {
            ((mag as i64).wrapping_neg(), used)
        }
    } else if ovf || mag > i64::MAX as u64 {
        set_errno(ERANGE);
        (i64::MAX, used)
    } else {
        (mag as i64, used)
    }
}

/// C `strtoul(s, &end, 10)` -> (value, bytes consumed). A leading '-'
/// negates in unsigned arithmetic (so "-1" is ULONG_MAX, no error).
pub fn strtoul(s: &str) -> (u64, usize) {
    let (neg, ovf, mag, used) = parse_dec(s);
    if used == 0 {
        set_errno(EINVAL);
        return (0, 0);
    }
    if ovf {
        set_errno(ERANGE);
        return (u64::MAX, used);
    }
    (if neg { mag.wrapping_neg() } else { mag }, used)
}

fn eq_ci(b: &[u8], i: usize, word: &[u8]) -> bool {
    b.len() >= i + word.len() && b[i..i + word.len()].eq_ignore_ascii_case(word)
}

/// C `strtod(s, &end)` with macOS semantics -> (value, bytes consumed).
/// Sets errno=ERANGE on overflow and underflow. Decimal subnormal results
/// report underflow; hexadecimal subnormals report it only when inexact.
pub fn strtod(s: &str) -> (f64, usize) {
    let b = s.as_bytes();
    let mut i = skip_ws(b);
    let mut neg = false;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        neg = b[i] == b'-';
        i += 1;
    }
    let sgn = |x: f64| if neg { -x } else { x };
    // inf / infinity
    if eq_ci(b, i, b"inf") {
        let n = if eq_ci(b, i, b"infinity") { 8 } else { 3 };
        return (sgn(f64::INFINITY), i + n);
    }
    // nan, nan(n-char-seq)
    if eq_ci(b, i, b"nan") {
        let mut j = i + 3;
        if j < b.len() && b[j] == b'(' {
            let mut k = j + 1;
            while k < b.len() && (b[k].is_ascii_alphanumeric() || b[k] == b'_') {
                k += 1;
            }
            if k < b.len() && b[k] == b')' {
                j = k + 1;
            }
        }
        return (sgn(f64::NAN), j);
    }
    // hex float
    if i + 1 < b.len() && b[i] == b'0' && (b[i + 1] == b'x' || b[i + 1] == b'X') {
        if let Some((x, used)) = parse_hexfloat(b, i + 2) {
            return (sgn(x), used);
        }
        // "0x" without hex digits: just the "0"
        return (sgn(0.0), i + 1);
    }
    // decimal
    let start = i;
    let mut j = i;
    let mut ndig = 0;
    let mut nonzero = false;
    while j < b.len() && b[j].is_ascii_digit() {
        nonzero |= b[j] != b'0';
        j += 1;
        ndig += 1;
    }
    if j < b.len() && b[j] == b'.' {
        j += 1;
        while j < b.len() && b[j].is_ascii_digit() {
            nonzero |= b[j] != b'0';
            j += 1;
            ndig += 1;
        }
    }
    if ndig == 0 {
        // Darwin strtod leaves errno unchanged when no conversion occurs.
        return (0.0, 0);
    }
    if j < b.len() && (b[j] == b'e' || b[j] == b'E') {
        let mut k = j + 1;
        if k < b.len() && (b[k] == b'+' || b[k] == b'-') {
            k += 1;
        }
        if k < b.len() && b[k].is_ascii_digit() {
            while k < b.len() && b[k].is_ascii_digit() {
                k += 1;
            }
            j = k;
        }
    }
    // Rust accepts the same decimal significands (including a leading dot).
    // Parse the consumed prefix directly without allocating another string.
    let text = &s[start..j];
    // Rust's parser is correctly rounded like macOS strtod; it does not
    // accept huge exponents gracefully in all cases, so clamp them.
    let x: f64 = match text.parse::<f64>() {
        Ok(v) => v,
        Err(_) => clamp_parse(text),
    };
    if x.is_infinite() || (nonzero && (x == 0.0 || x.is_subnormal())) {
        set_errno(ERANGE);
    }
    (sgn(x), j)
}

fn clamp_parse(text: &str) -> f64 {
    // Exponent too large to parse: decide by its sign.
    let epos = text.find(['e', 'E']).unwrap_or(text.len());
    let neg_exp = text[epos..].contains('-');
    if neg_exp {
        0.0
    } else {
        f64::INFINITY
    }
}

fn parse_hexfloat(b: &[u8], mut i: usize) -> Option<(f64, usize)> {
    let mut mant: u64 = 0;
    let mut sticky = false;
    let mut exp: i64 = 0;
    let mut ndig = 0;
    let mut seen_dot = false;
    loop {
        if i < b.len() && b[i].is_ascii_hexdigit() {
            let d = (b[i] as char).to_digit(16).unwrap() as u64;
            if mant >> 56 == 0 {
                mant = (mant << 4) | d;
                if seen_dot {
                    exp -= 4;
                }
            } else {
                sticky |= d != 0;
                if !seen_dot {
                    exp += 4;
                }
            }
            ndig += 1;
            i += 1;
        } else if i < b.len() && b[i] == b'.' && !seen_dot {
            seen_dot = true;
            i += 1;
        } else {
            break;
        }
    }
    if ndig == 0 {
        return None;
    }
    if i < b.len() && (b[i] == b'p' || b[i] == b'P') {
        let mut k = i + 1;
        let mut eneg = false;
        if k < b.len() && (b[k] == b'+' || b[k] == b'-') {
            eneg = b[k] == b'-';
            k += 1;
        }
        if k < b.len() && b[k].is_ascii_digit() {
            let mut e: i64 = 0;
            while k < b.len() && b[k].is_ascii_digit() {
                e = e.saturating_mul(10).saturating_add((b[k] - b'0') as i64);
                k += 1;
            }
            exp = exp.saturating_add(if eneg { -e } else { e });
            i = k;
        }
    }
    if mant == 0 {
        return Some((0.0, i));
    }

    // mant holds the leading (at most 60) bits; sticky records any nonzero
    // digits beyond them. Round directly to the destination's spacing so a
    // subnormal is not rounded once to 53 bits and then again during scaling.
    let top = exp.saturating_add(63 - mant.leading_zeros() as i64);
    if top > 1023 {
        set_errno(ERANGE);
        return Some((f64::INFINITY, i));
    }
    let mut unit_exp = top.saturating_sub(52).max(-1074);
    let shift = unit_exp.saturating_sub(exp);
    let (mut sig, inexact) = if shift >= 64 {
        (0, true)
    } else if shift > 0 {
        let rem = mant & ((1u64 << shift) - 1);
        let half = 1u64 << (shift - 1);
        let sig = mant >> shift;
        let up = rem > half || (rem == half && (sticky || sig & 1 != 0));
        (sig + u64::from(up), rem != 0 || sticky)
    } else {
        (mant << -shift, false)
    };
    if sig == 1u64 << 53 {
        sig >>= 1;
        unit_exp += 1;
    }
    let bits = if unit_exp == -1074 {
        sig // Includes rounding up to the smallest normal.
    } else {
        ((unit_exp + 1075) as u64) << 52 | (sig & ((1u64 << 52) - 1))
    };
    let x = f64::from_bits(bits);
    // Darwin detects hex underflow after rounding, and accepts exact
    // subnormals. Successful conversions leave any previous errno intact.
    if x.is_infinite() || (x < f64::MIN_POSITIVE && inexact) {
        set_errno(ERANGE);
    }
    Some((x, i))
}

/// C `atoi` (strtol truncated to int, no errno semantics needed).
pub fn atoi(s: &str) -> i32 {
    strtol(s).0 as i32
}
/// C `atof` (= strtod).
pub fn atof(s: &str) -> f64 {
    strtod(s).0
}

/// eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX}
pub fn eatoh(s: &str) -> i16 {
    let n = strtol(s).0;
    if n > i16::MAX as i64 || n < i16::MIN as i64 || errno() == ERANGE {
        crate::err!("%s: eatoh: overflow", "atopkge.c");
    }
    n as i16
}
/// eatou - convert string s to unsigned short integer {0:USHRT_MAX}
pub fn eatou(s: &str) -> u16 {
    let n = strtoul(s).0;
    if n > u16::MAX as u64 || errno() == ERANGE {
        crate::err!("%s: eatou: overflow", "atopkge.c");
    }
    n as u16
}
/// eatoi - convert string s to integer {INT_MIN:INT_MAX}
pub fn eatoi(s: &str) -> i32 {
    let n = strtol(s).0;
    if n > i32::MAX as i64 || n < i32::MIN as i64 || errno() == ERANGE {
        crate::err!("%s: eatoi: overflow", "atopkge.c");
    }
    n as i32
}
/// eatop - convert string s to unsigned integer {0:UINT_MAX}
pub fn eatop(s: &str) -> u32 {
    let n = strtoul(s).0;
    if n > u32::MAX as u64 || errno() == ERANGE {
        crate::err!("%s: eatop: overflow", "atopkge.c");
    }
    n as u32
}
/// eatol - convert string s to long integer {LONG_MIN:LONG_MAX}
pub fn eatol(s: &str) -> i64 {
    let n = strtol(s).0;
    if errno() == ERANGE {
        crate::err!("%s: eatol: overflow", "atopkge.c");
    }
    n
}
/// eatov - convert string s to unsigned long {0:ULONG_MAX}
pub fn eatov(s: &str) -> u64 {
    let n = strtoul(s).0;
    if errno() == ERANGE {
        crate::err!("%s: eatov: overflow", "atopkge.c");
    }
    n
}
/// eatof - convert string s to float {-FLT_MAX:FLT_MAX}
/// (the double from strtod is narrowed to float *before* the range test).
pub fn eatof(s: &str) -> f32 {
    let x = strtod(s).0 as f32;
    if x > f32::MAX || x < -f32::MAX || errno() == ERANGE {
        crate::err!("%s: eatof: overflow", "atopkge.c");
    }
    x
}
/// eatod - convert string s to double {-DBL_MAX:DBL_MAX}
pub fn eatod(s: &str) -> f64 {
    let x = strtod(s).0;
    if errno() == ERANGE || x > f64::MAX || x < -f64::MAX {
        crate::err!("%s: eatod: overflow", "atopkge.c");
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup(args: &[&str]) {
        let v: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        initargs(&v);
        set_errno(0);
    }

    #[test]
    fn strtod_matches_macos_probe() {
        let cases: &[(&str, f64, usize, bool)] = &[
            ("4.9406564584124654e-324", 5e-324, 23, true),
            ("0x1p-1074", 5e-324, 9, false),
            ("2.2250738585072014e-308", 2.2250738585072014e-308, 23, false),
            ("1e-400", 0.0, 6, true),
            ("1.7976931348623159e308", f64::INFINITY, 22, true),
            ("0x1.8p1", 3.0, 7, false),
            ("0x", 0.0, 1, false),
            ("1e", 1.0, 1, false),
            ("1e+", 1.0, 1, false),
            ("-.5e-1x", -0.05, 6, false),
            ("+.e5", 0.0, 0, false),
            ("  \t\n12", 12.0, 6, false),
            ("infinit", f64::INFINITY, 3, false),
            ("INFINITY", f64::INFINITY, 8, false),
            ("0x.8", 0.5, 4, false),
            ("0x1P+3", 8.0, 6, false),
            ("1e5000", f64::INFINITY, 6, true),
            ("0e-999", 0.0, 6, false),
        ];
        for &(s, v, used, erange) in cases {
            set_errno(0);
            let (x, u) = strtod(s);
            assert_eq!((x, u), (v, used), "{}", s);
            assert_eq!(errno() == ERANGE, erange, "{}", s);
        }
        set_errno(0);
        let (x, u) = strtod("nan(abc)");
        assert!(x.is_nan() && u == 8);
        assert_eq!(strtod("nAn(abc").1, 3);
    }

    #[test]
    fn strtol_strtoul() {
        assert_eq!(strtol("  +7x"), (7, 4));
        assert_eq!(strtol("0x10"), (0, 1));
        assert_eq!(strtoul("-1"), (u64::MAX, 2));
        set_errno(0);
        assert_eq!(strtoul("-9223372036854775809").0, 9223372036854775807);
        assert_eq!(errno(), 0);
        assert_eq!(strtol("99999999999999999999").0, i64::MAX);
        assert_eq!(errno(), ERANGE);
    }

    #[test]
    fn hexfloat_rounding_and_range() {
        let cases: &[(&str, u64, bool)] = &[
            ("0x1.1p-1074", 1, true),
            ("0x1.8p-1074", 2, true),
            ("0x1.0000000000001p-1023", 1 << 51, true),
            ("0x1p-1074", 1, false),
            ("0x1.8p-1073", 3, false),
            ("0x1p-1023", 1 << 51, false),
            ("0x0.fffffffffffffp-1022", (1 << 52) - 1, false),
            ("0x1p-1075", 0, true),
            // Bits beyond the retained mantissa must break a midpoint tie.
            ("0x1.0000000000000001p-1075", 1, true),
            ("0x1.000000000000000000000000000001p-1075", 1, true),
            ("0x2.fffffffffffffffep-1075", 1, true),
            ("0x3p-1075", 2, true),
            ("0x0.fffffffffffff7p-1022", (1 << 52) - 1, true),
            // Tininess is detected after rounding, including at the tie.
            ("0x0.fffffffffffff8p-1022", 1 << 52, false),
            ("0x0.fffffffffffff9p-1022", 1 << 52, false),
            ("0x1.00000000000008p0", 1.0f64.to_bits(), false),
            ("0x1.00000000000008000001p0", 1.0f64.to_bits() + 1, false),
            ("0x1.00000000000018p0", 1.0f64.to_bits() + 2, false),
            ("0x1.fffffffffffff4p1023", f64::MAX.to_bits(), false),
            ("0x1.fffffffffffff8p1023", f64::INFINITY.to_bits(), true),
            ("0x1p9999999999999999999999", f64::INFINITY.to_bits(), true),
            ("0x1p-9999999999999999999999", 0, true),
            ("0x.1p-9999999999999999999999", 0, true),
            ("0x0p9999999999999999999999", 0, false),
            ("0x0p-9999999999999999999999", 0, false),
        ];
        for &(s, bits, range) in cases {
            for sign in ["", "-"] {
                let text = format!(" \t{}{}!", sign, s);
                // A successful conversion must preserve even a nonzero errno.
                set_errno(EINVAL);
                let (x, used) = strtod(&text);
                let sign_bit = if sign == "-" { 1u64 << 63 } else { 0 };
                assert_eq!(x.to_bits(), bits | sign_bit, "{}", text);
                assert_eq!(used, text.len() - 1, "{}", text);
                assert_eq!(errno(), if range { ERANGE } else { EINVAL }, "{}", text);
            }
        }
    }

    #[test]
    fn decimal_conversion_compatibility() {
        // Decimal coordinates from the elamodel regression: preserve both
        // strtod's double result and eatof's subsequent float rounding.
        let cases = [
            ("0.5280995857712361", 0x3fe0e6311a3d6538, 0x3f073189),
            ("0.5816681028436204", 0x3fe29d066cdadef0, 0x3f14e833),
            ("0.385852903967372", 0x3fd8b1d060e6d0c1, 0x3ec58e83),
            ("0.6403706677278277", 0x3fe47deaa066adb1, 0x3f23ef55),
            ("0.3497277229821637", 0x3fd661f0632da0d2, 0x3eb30f83),
            ("0.5721368953465087", 0x3fe24ef208cb24b2, 0x3f127790),
            ("0.45123948844005607", 0x3fdce11b9760e247, 0x3ee708dd),
            ("0.7170173854732829", 0x3fe6f1ce71a8ad78, 0x3f378e74),
            ("0.528327681141805", 0x3fe0e80f74026b0e, 0x3f07407c),
            ("0.652199668607643", 0x3fe4ded1d6e4383c, 0x3f26f68f),
            ("0.2328092997054671", 0x3fcdccb1f438463e, 0x3e6e6590),
        ];
        for (s, double_bits, float_bits) in cases {
            set_errno(EINVAL);
            let text = format!("{},ignored", s);
            let (x, used) = strtod(&text);
            assert_eq!(used, s.len());
            assert_eq!(x.to_bits(), double_bits, "{}", s);
            assert_eq!(eatod(&text).to_bits(), double_bits, "{}", s);
            assert_eq!(eatof(&text).to_bits(), float_bits, "{}", s);
            assert_eq!(errno(), EINVAL);
        }
    }

    #[test]
    fn strtod_no_conversion_preserves_errno() {
        for s in ["", " ", "+", "-", ".", "+.e1", "garbage"] {
            for previous in [0, EINVAL, ERANGE] {
                set_errno(previous);
                assert_eq!(strtod(s), (0.0, 0), "{}", s);
                assert_eq!(errno(), previous, "{}", s);
            }
        }
    }

    #[test]
    fn checked_conversion_range_errors() {
        const CHILD: &str = "OXWAVE_TEST_CHECKED_RANGE";
        if let Ok(conv) = std::env::var(CHILD) {
            setup(&["range-probe"]);
            if conv == "eatof" {
                eatof("0x1.8p-1074");
            } else if conv == "eatod" {
                eatod("0x1.0000000000001p-1023");
            } else {
                // Like atopkge.c, every eato* also rejects a stale ERANGE.
                atof("0x1.1p-1074");
                match conv.as_str() {
                    "eatoh" => { eatoh("1"); }
                    "eatou" => { eatou("1"); }
                    "eatoi" => { eatoi("1"); }
                    "eatop" => { eatop("1"); }
                    "eatol" => { eatol("1"); }
                    "eatov" => { eatov("1"); }
                    "eatof-invalid" => { eatof("garbage"); }
                    "eatod-invalid" => { eatod("."); }
                    _ => unreachable!(),
                }
            }
            panic!("{} accepted ERANGE", conv);
        }
        for conv in ["eatoh", "eatou", "eatoi", "eatop", "eatol", "eatov", "eatof", "eatod",
                     "eatof-invalid", "eatod-invalid"] {
            let out = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "par::tests::checked_conversion_range_errors", "--nocapture"])
                .env(CHILD, conv)
                .output()
                .unwrap();
            assert_eq!(out.status.code(), Some(1), "{}: {:?}", conv, out);
            assert_eq!(String::from_utf8(out.stderr).unwrap(),
                format!("\nrange-probe: atopkge.c: {}: overflow\n", conv.split('-').next().unwrap()));
        }
        set_errno(0);
        assert_eq!(eatof("0x1p-1074"), 0.0);
        assert_eq!(eatod("0x1p-1074").to_bits(), 1);
        assert_eq!(eatof("0x1p-149").to_bits(), 1);
        assert_eq!(errno(), 0);
    }

    #[test]
    fn getpar_semantics() {
        setup(&["prog", "n=1", "v=1,,2", "w=3,", "n=5", "s=a,b", "e=", "f=2.5", "x", "key=sx,,gx,"]);
        let mut n = 0i32;
        assert_eq!(getparint("n", &mut n), 1);
        assert_eq!(n, 5);
        assert_eq!(getnparint(1, "n", &mut n), 1);
        assert_eq!(n, 1);
        let mut v = vec![9i32; 1];
        assert_eq!(getparint("v", &mut v), 3);
        assert_eq!(v, vec![1, 0, 2]);
        let mut w = [0i32; 4];
        assert_eq!(getparint("w", &mut w[..]), 1);
        assert_eq!(countparval("w"), 2);
        assert_eq!(countparval("e"), 1);
        assert_eq!(getparint("e", &mut n), 0);
        let mut s = String::new();
        assert_eq!(getparstring("e", &mut s), 1);
        assert_eq!(s, "");
        assert_eq!(getparstring("s", &mut s), 1);
        assert_eq!(s, "a,b");
        let mut f = 0f32;
        assert_eq!(getparfloat("f", &mut f), 1);
        assert_eq!(f, 2.5);
        let mut keys: Vec<String> = Vec::new();
        assert_eq!(getparstringarray("key", &mut keys), 3);
        assert_eq!(keys, vec!["sx", "", "gx"]);
        assert_eq!(countparname("n"), 2);
        assert_eq!(countparname("x"), 0);
        assert_eq!(getparint("nothere", &mut n), 0);
    }

    #[test]
    fn parfile() {
        let dir = std::env::temp_dir().join(format!("oxw_par_{}", std::process::id()));
        std::fs::write(&dir, b"a=1 b=2\n c='x y' junk\n\ta=7").unwrap();
        let pf = format!("par={}", dir.display());
        setup(&["prog", "a=3", &pf]);
        let mut a = 0i32;
        assert_eq!(getparint("a", &mut a), 1);
        assert_eq!(a, 3);
        assert_eq!(getnparint(1, "a", &mut a), 1);
        assert_eq!(a, 1);
        assert_eq!(getnparint(2, "a", &mut a), 1);
        assert_eq!(a, 7);
        let mut c = String::new();
        getparstring("c", &mut c);
        assert_eq!(c, "'x y'");
        assert_eq!(countparname("a"), 3);
        let mut b = 0i32;
        getparint("b", &mut b);
        assert_eq!(b, 2);
        std::fs::remove_file(&dir).unwrap();
    }

    #[test]
    fn conversions() {
        set_errno(0);
        assert_eq!(eatoi(" -12abc"), -12);
        assert_eq!(eatoi("abc"), 0);
        assert_eq!(eatov("-1"), u64::MAX);
        assert_eq!(eatof("3.4028235e38"), f32::MAX);
        assert_eq!(eatof("1e-50"), 0.0);
        assert_eq!(eatoh("-32768"), i16::MIN);
    }
}
