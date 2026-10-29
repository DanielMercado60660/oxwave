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

//! Differential tests for su/lib/getSPSfield.c.
//! Author: Fernando Roxo <petro@roxo.org>.

use oxwave_cref::XorShift;
use oxwave_sps_records::{getSPSfield, SPSField};
use std::ffi::{c_char, c_int, c_void};
use std::sync::Mutex;

// The C routine selects its parse table through a process-global pointer.
static C_LOCK: Mutex<()> = Mutex::new(());

mod original {
    use super::{c_char, c_int, c_void};

    extern "C" {
        pub fn getSPSfield(line: *mut c_char, rev: c_int, idx: c_int, addr: *mut c_void) -> c_int;
    }
}

#[derive(Clone, Copy, Debug)]
enum Kind {
    Char,
    String,
    Int,
    Float,
    Double,
}

// Independent transcription of the original C record layouts.
fn fields(rev: i32, label: u8) -> [(usize, usize, Kind); 15] {
    use Kind::*;
    match (rev, label) {
        (0, b'R' | b'S') => [
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
        ],
        (2, b'R' | b'S') => [
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
        ],
        (0, b'X') => [
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
        ],
        (2, b'X') => [
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
        ],
        _ => unreachable!(),
    }
}

fn compare(line: &[u8], rev: i32, idx: i32, kind: Kind) {
    // C needs writable, NUL-terminated storage; Rust accepts either a slice
    // with a NUL or the slice's end. Include guard bytes to check restoration.
    let mut c_line = line.to_vec();
    c_line.extend_from_slice(&[0, 0xa5, 0x5a]);
    let before = c_line.clone();
    let context = || format!("rev={rev} idx={idx} kind={kind:?} line={line:?}");

    macro_rules! numeric {
        ($initial:expr, $variant:ident, $equal:expr) => {{
            let mut c = $initial;
            let mut rust = $initial;
            // SAFETY: layout selects this exact C type, input is writable and
            // terminated, and idx is in the original table's 0..15 range.
            let rc = unsafe {
                original::getSPSfield(
                    c_line.as_mut_ptr().cast(),
                    rev,
                    idx,
                    std::ptr::from_mut(&mut c).cast(),
                )
            };
            let rr = getSPSfield(line, rev, idx, SPSField::$variant(&mut rust));
            assert_eq!(rc, rr, "{}", context());
            assert!($equal(c, rust), "{}: C={c:?} Rust={rust:?}", context());
        }};
    }

    match kind {
        Kind::Char => numeric!(0xa5u8, Char, |a: u8, b: u8| a == b),
        Kind::Int => numeric!(-123456789i32, Int, |a: i32, b: i32| a == b),
        Kind::Float => numeric!(-123.5f32, Float, |a: f32, b: f32| a.to_bits()
            == b.to_bits()
            || (a.is_nan() && b.is_nan())),
        Kind::Double => numeric!(-123.5f64, Double, |a: f64, b: f64| a.to_bits()
            == b.to_bits()
            || (a.is_nan() && b.is_nan())),
        Kind::String => {
            // Maximum %s field width is six bytes plus its terminating NUL.
            // Compare the whole destination, including untouched guard bytes.
            let mut c = [0xa5u8; 32];
            let mut rust = c;
            let rc = unsafe {
                original::getSPSfield(c_line.as_mut_ptr().cast(), rev, idx, c.as_mut_ptr().cast())
            };
            let rr = getSPSfield(line, rev, idx, SPSField::String(&mut rust));
            assert_eq!((rc, c), (rr, rust), "{}", context());
        }
    }
    assert_eq!(c_line, before, "C must restore the input: {}", context());
}

fn record(label: u8, start: usize, width: usize, token: &[u8]) -> Vec<u8> {
    let mut line = vec![b' '; 80];
    line[0] = label;
    let n = width.min(token.len());
    line[start..start + n].copy_from_slice(&token[..n]);
    line
}

#[test]
fn get_sps_field_matches_c_random_records() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x12345678);
    let alphabet = b"  0123456789-+.eExpPinNaAfF()xYZ\t\r\n\x0b\x0c";
    for rev in [0, 2] {
        for label in [b'R', b'S', b'X'] {
            for (idx, &(_, _, kind)) in fields(rev, label).iter().enumerate() {
                for trial in 0..1000 {
                    let n = match trial % 12 {
                        0 => 1,
                        1 => 2,
                        2 => 17,
                        3 => 32,
                        4 => 79,
                        5 => 80,
                        6 => 127,
                        7 => 256,
                        8 => 4096,
                        _ => (rng.next_u32() % 100 + 1) as usize,
                    };
                    let mut line = vec![0u8; n];
                    line[0] = label;
                    for ch in &mut line[1..] {
                        *ch = alphabet[(rng.next_u32() as usize) % alphabet.len()];
                    }
                    // Some inputs contain an early NUL and nonzero bytes after it.
                    if trial % 7 == 0 && n > 1 {
                        line[1 + rng.next_u32() as usize % (n - 1)] = 0;
                    }
                    compare(&line, rev, idx as i32, kind);
                }
            }
        }
    }
}

#[test]
fn get_sps_field_matches_c_generated_numbers() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut rng = XorShift::new(0x5a5073);
    for rev in [0, 2] {
        for label in [b'R', b'S', b'X'] {
            for (idx, &(start, width, kind)) in fields(rev, label).iter().enumerate().skip(1) {
                for trial in 0..1000 {
                    let signed = rng.next_u32() as i32;
                    let token = match kind {
                        Kind::Char => vec![rng.next_u32() as u8],
                        Kind::String => format!("T{:05}", rng.next_u32() % 100000).into_bytes(),
                        Kind::Int => match trial % 5 {
                            0 => b"0".to_vec(),
                            1 => b"-0".to_vec(),
                            2 => format!("{signed:+}").into_bytes(),
                            3 => format!("\t{}", signed % 1000).into_bytes(),
                            _ => {
                                format!("{}", signed % 10i32.pow(width.min(8) as u32)).into_bytes()
                            }
                        },
                        Kind::Float | Kind::Double => {
                            let exponent = rng.next_u32() as i32 % 1100;
                            match trial % 8 {
                                0 => format!("{}", (signed % 100000) as f64 / 10.0).into_bytes(),
                                1 => format!("{signed}e{exponent}").into_bytes(),
                                2 => format!(
                                    "{}.{:06}e{exponent}",
                                    signed % 10,
                                    rng.next_u32() % 1000000
                                )
                                .into_bytes(),
                                3 => format!("0x{:x}p{exponent}", rng.next_u32()).into_bytes(),
                                4 => format!(
                                    "-0x{:x}.{:x}p{exponent}",
                                    rng.next_u32() % 16,
                                    rng.next_u32() % 256
                                )
                                .into_bytes(),
                                5 => format!(
                                    "0x{:x}{:08x}",
                                    rng.next_u32() % 0x100000,
                                    rng.next_u32()
                                )
                                .into_bytes(),
                                6 => format!(
                                    "{}e{exponent}",
                                    if trial % 3 == 0 { "-0" } else { "0" }
                                )
                                .into_bytes(),
                                _ => format!("\x0b{:+}", signed % 1000).into_bytes(),
                            }
                        }
                    };
                    let line = record(label, start, width, &token);
                    compare(&line, rev, idx as i32, kind);
                }
            }
        }
    }
}

#[test]
fn get_sps_field_matches_c_edges() {
    let _guard = C_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let tokens: &[&[u8]] = &[
        b"",
        b"0",
        b"-0",
        b"+0",
        b"1.25e-3",
        b"-9.9E+12",
        b"1e+",
        b"0x1.8p2",
        b"-0x1.2p-3",
        b"inf",
        b"nan",
        b".",
        b"+",
        b"-",
        b" ",
        b"\x0b1",
        b"\x0c-2",
        b"\t 3",
        b"1e9999",
        b"0e9999",
        b"-0e9999",
        b"1e-9999",
        b"1e-45",
        b"3.402824e38",
        b"0x1p128",
        b"0x1p-149",
        b"0x1p-150",
        b"0x1p-151",
        b"0x3p-150",
        b"0x5p-150",
        b"0x7p-150",
        b"-0x1p-150",
        b"0x0p9999",
        b"-0x0p9999",
        b"0x1p9999",
        b"0x1p-9999",
        b"0x",
        b"0x.",
        b"0xg",
        b"0xp",
        b"0x1p",
        b"0x1p+",
        b"0x1p-",
        b"in",
        b"infi",
        b"infin",
        b"infinit",
        b"infinity",
        b"nan(1)",
        b"-nan",
        b"nan(foo)",
        b"1e",
        b"1e-",
        b"1.2.3",
        b"-.5",
        b"+.5",
        b"1_2",
        b"12\0garbage",
        b"\xff12",
        b"\x8512",
        b"\xa012",
        b"000000",
        b"1234567890123456",
        b"-123456789012345",
        b"99999999",
        b"-9999999",
        b"1e-308",
        b"1e309",
        b"-1e309",
        b"0x1p-1022",
        b"0x1p-1023",
        b"0x1p-1074",
        b"0x1p-1075",
        b"0x3p-1075",
        b"0x5p-1075",
        b"-0x1p-1074",
        b"0x1p1023",
        b"0x1p1024",
        b"-0x1p-9999",
        b"0x10000010000001",
        b"0x1000000fffffff",
        b"0x80000080000001",
        b"0x8000007fffffff",
        b"0x1.fffffep127",
        b"0x1.ffffffp127",
        b"\t\r\n\x0b\x0c ",
        b"abc def",
        b"abc\x0bdef",
        b"0XABC.P+0",
        b"+INF",
    ];
    for rev in [0, 2] {
        for label in [b'R', b'S', b'X'] {
            for (idx, &(start, width, kind)) in fields(rev, label).iter().enumerate() {
                if idx > 0 {
                    for &token in tokens {
                        compare(&record(label, start, width, token), rev, idx as i32, kind);
                    }
                }
                // Empty field, one byte before/after the start and end, with
                // whitespace, numeric text, and NUL-terminated shorter records.
                for fill in [b' ', b'0', b'9', b'\x0b'] {
                    for n in [
                        start.saturating_sub(1).max(1),
                        start.max(1),
                        start + 1,
                        start + width - 1,
                        start + width,
                        start + width + 1,
                    ] {
                        let n = n.max(1);
                        let mut line = vec![fill; n];
                        line[0] = label;
                        compare(&line, rev, idx as i32, kind);
                        line.extend_from_slice(b"\0ignored");
                        compare(&line, rev, idx as i32, kind);
                    }
                }
            }
        }
    }
}

#[test]
fn get_sps_field_matches_c_errors() {
    // Capture both implementations' diagnostics through child processes so
    // no other tests' stderr or libc FILE state needs to be redirected.
    const CHILD: &str = "OXWAVE_SPS_RECORDS_ERROR_CHILD";
    if let Ok(side) = std::env::var(CHILD) {
        let rev: i32 = std::env::var("OXWAVE_SPS_RECORDS_REV")
            .unwrap()
            .parse()
            .unwrap();
        let label: u8 = std::env::var("OXWAVE_SPS_RECORDS_LABEL")
            .unwrap()
            .parse()
            .unwrap();
        let mut output = 0xa5u8;
        let rc = if side == "C" {
            let mut line = [label, 0];
            unsafe {
                original::getSPSfield(
                    line.as_mut_ptr().cast(),
                    rev,
                    0,
                    std::ptr::from_mut(&mut output).cast(),
                )
            }
        } else {
            getSPSfield(&[label], rev, 0, SPSField::Char(&mut output))
        };
        assert_eq!(rc, 1);
        assert_eq!(output, 0xa5);
        return;
    }

    for (rev, label) in [
        (i32::MIN, b'R'),
        (-1, b'X'),
        (1, b'H'),
        (3, b'S'),
        (i32::MAX, 0),
        (0, b'H'),
        (2, b'r'),
        (0, b'x'),
        (2, 0),
        (0, 0xff),
    ] {
        let expected = if rev != 0 && rev != 2 {
            format!("****  Unsupported SPS format {rev}!!\n\n").into_bytes()
        } else {
            let mut message = b"****  Unsupported SPS record ".to_vec();
            message.push(label);
            message.extend_from_slice(b"!!\n\n");
            message
        };
        for side in ["C", "Rust"] {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "get_sps_field_matches_c_errors", "--nocapture"])
                .env(CHILD, side)
                .env("OXWAVE_SPS_RECORDS_REV", rev.to_string())
                .env("OXWAVE_SPS_RECORDS_LABEL", label.to_string())
                .output()
                .unwrap();
            assert!(output.status.success(), "{side}: {output:?}");
            assert_eq!(output.stderr, expected, "{side}: rev={rev} label={label}");
        }
    }
}
