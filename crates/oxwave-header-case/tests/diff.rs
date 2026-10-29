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

//! Differential tests against su/lib/headcase.c in the original SU library.
//! Credits: the original SU source contains no named author credits.

use oxwave_core::header::{SegyHeader, HDR, HDRBYTES, MAXSEGY, SU_NFLTS};
use oxwave_cref::XorShift;
use oxwave_header_case as rust;
use std::ffi::{c_char, CString};

// C segy includes inline sample storage; Trace uses a Vec and is not ABI-compatible.
#[repr(C)]
#[derive(Clone, Copy)]
struct CSegy {
    hdr: SegyHeader,
    data: [f32; SU_NFLTS],
}

extern "C" {
    // int GetCase(char* cbuf);
    #[allow(non_snake_case)]
    fn GetCase(cbuf: *mut c_char) -> i32;
    // double fromhead(segy tr, int k);
    fn fromhead(tr: CSegy, k: i32) -> f64;
    // void tohead(segy *tr, int k, double dvalue);
    fn tohead(tr: *mut CSegy, k: i32, dvalue: f64);
}

const SIZES: &[u16] = &[0, 1, 2, 3, 7, 16, 31, 64, 257, 1024, 4093, 65535];

fn check_name(name: &str) {
    // Include the whole Rust string, including any embedded NUL and its suffix.
    let mut bytes = name.as_bytes().to_vec();
    bytes.push(0);
    let c = unsafe { GetCase(bytes.as_mut_ptr().cast()) };
    assert_eq!(rust::GetCase(name), c, "name {name:?}");
}

#[test]
fn get_case_matches_c() {
    for (i, key) in HDR.iter().enumerate() {
        let ckey = CString::new(key.key).unwrap();
        assert_eq!(
            unsafe { GetCase(ckey.as_ptr().cast_mut()) },
            i as i32 + 1
        );
        assert_eq!(rust::GetCase(key.key), i as i32 + 1);
        for name in [
            key.key.to_uppercase(),
            format!("{}x", key.key),
            format!(" {}", key.key),
            format!("{} ", key.key),
            format!("{}\0garbage", key.key),
        ] {
            check_name(&name);
        }
    }
    for name in [
        "",
        "n",
        "nu",
        "nul",
        "null",
        "nullsuffix",
        "null-100",
        "numb",
        "numb+",
        "numb-",
        "numbxyz",
        "numb  -12trailing",
        "numb0x10",
        "numb+00123",
        "numb2147483647",
        "numb-2147483648",
        "numb2147483648",
        "numb-2147483649",
        "numb9223372036854775807",
        "numb-9223372036854775808",
        "numb99999999999999999999999999",
        "numb-99999999999999999999999999",
        "numb\0-42",
        "null\0x",
        "é",
        "numbé",
        "tracl\0junk",
    ] {
        check_name(name);
    }
    let mut rng = XorShift::new(0x43617365);
    for _ in 0..20_000 {
        let n = rng.next_u32() as i32;
        let whitespace = ["", " ", "\t", "\n", "\r", "\x0b", "\x0c"][rng.next_u32() as usize % 7];
        check_name(&format!("numb{whitespace}{n}suffix"));
        check_name(&format!("null{n}"));
        let key = HDR[rng.next_u32() as usize % HDR.len()].key;
        check_name(key);
        check_name(&format!("{key}{n}"));
    }
    // Empty/short strings, primes, powers of two, and very long unrecognized names.
    for &size in SIZES {
        for _ in 0..16 {
            let mut name = String::with_capacity(size as usize);
            for _ in 0..size {
                name.push((b'a' + (rng.next_u32() % 26) as u8) as char);
            }
            check_name(&name);
            check_name(&format!("null{name}"));
            check_name(&format!("numb{name}"));
        }
    }
}

fn random_header(rng: &mut XorShift) -> SegyHeader {
    let mut hdr = SegyHeader::zeroed();
    for byte in hdr.as_bytes_mut() {
        *byte = rng.next_u32() as u8;
    }
    hdr
}

fn new_c_trace() -> Box<CSegy> {
    assert_eq!(std::mem::size_of::<SegyHeader>(), HDRBYTES);
    assert_eq!(std::mem::size_of::<CSegy>(), MAXSEGY);
    assert_eq!(std::mem::offset_of!(CSegy, data), HDRBYTES);
    Box::new(CSegy {
        hdr: SegyHeader::zeroed(),
        data: [0.0; SU_NFLTS],
    })
}

#[test]
fn fromhead_matches_c() {
    // Passing C segy by value requires several full trace copies in debug mode.
    std::thread::Builder::new()
        .stack_size(8 * 1024 * 1024)
        .spawn(|| {
            let mut rng = XorShift::new(0x46726f6d);
            let mut tr = new_c_trace();
            let float_edges = [
                0.0f32,
                -0.0,
                1.0,
                -1.0,
                f32::MIN_POSITIVE,
                f32::from_bits(1),
                f32::MAX,
                -f32::MAX,
                f32::INFINITY,
                f32::NEG_INFINITY,
                f32::NAN,
                f32::from_bits(0x7f80_0001),
                f32::from_bits(0xffc1_2345),
            ];
            for i in 0..2048 {
                tr.hdr = random_header(&mut rng);
                if i < 16 {
                    for key in &HDR {
                        let bytes = tr.hdr.as_bytes_mut();
                        match key.type_ {
                            "i" => {
                                let v = [i32::MIN, i32::MAX, -1, 0, 1, 16_777_217][i % 6];
                                bytes[key.offs..key.offs + 4].copy_from_slice(&v.to_ne_bytes());
                            }
                            "h" | "u" => {
                                let v = [0u16, 1, 0x7fff, 0x8000, 0xffff][i % 5];
                                bytes[key.offs..key.offs + 2].copy_from_slice(&v.to_ne_bytes());
                            }
                            "f" => {
                                bytes[key.offs..key.offs + 4].copy_from_slice(
                                    &float_edges[i % float_edges.len()].to_ne_bytes(),
                                );
                            }
                            _ => unreachable!(),
                        }
                    }
                }
                tr.hdr.ns = SIZES[i % SIZES.len()];
                let before = *tr.hdr.as_bytes();
                for k in 1..=80 {
                    let c = unsafe { fromhead(*tr, k) };
                    let r = rust::fromhead(&tr.hdr, k);
                    assert_eq!(
                        r.to_bits(),
                        c.to_bits(),
                        "iteration {i}, case {k}: {r:?} vs {c:?}"
                    );
                }
                assert_eq!(tr.hdr.as_bytes(), &before);
            }
        })
        .unwrap()
        .join()
        .unwrap();
}

fn check_write(tr: &mut CSegy, k: i32, value: f64) {
    let mut hdr = tr.hdr;
    unsafe { tohead(tr, k, value) };
    rust::tohead(&mut hdr, k, value);
    assert_eq!(
        hdr.as_bytes(),
        tr.hdr.as_bytes(),
        "case {k}, value {value:?}"
    );
}

#[test]
fn tohead_matches_c() {
    let mut rng = XorShift::new(0x546f68656164);
    let mut tr = new_c_trace();
    for (i, x) in tr.data.iter_mut().enumerate() {
        *x = (i as i32 - 32768) as f32;
    }
    let edges = [
        0.0,
        -0.0,
        0.5,
        -0.5,
        1.5,
        -1.5,
        2.5,
        -2.5,
        32766.5,
        32767.5,
        -32768.5,
        -32769.5,
        65534.5,
        65535.5,
        -65535.5,
        2147483646.5,
        2147483647.5,
        -2147483648.5,
        -2147483649.5,
        4294967295.5,
        -4294967296.5,
        9_223_372_036_854_774_784.0,
        -9_223_372_036_854_775_808.0,
        9_223_372_036_854_775_808.0,
        -9_223_372_036_854_777_856.0,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        f32::MIN_POSITIVE as f64,
        f32::from_bits(1) as f64,
        f32::MAX as f64,
        -(f32::MAX as f64),
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ];
    for k in 1..=80 {
        for &value in &edges {
            tr.hdr = random_header(&mut rng);
            check_write(&mut tr, k, value);
        }
        for i in 0..4096 {
            tr.hdr = random_header(&mut rng);
            tr.hdr.ns = SIZES[i % SIZES.len()];
            let value = match i % 8 {
                0 => rng.uniform(-32768.0, 32768.0),
                1 => rng.uniform(-65536.0, 65536.0),
                2 => rng.uniform(-3e9, 3e9),
                3 => rng.uniform(-1e18, 1e18),
                4 => rng.range_i32(-1_000_000, 1_000_000) as f64 + 0.5,
                5 => rng.uniform(-1e-30, 1e-30),
                6 => rng.uniform(-1e38, 1e38),
                _ => f64::from_bits(rng.next_u64()),
            };
            check_write(&mut tr, k, value);
        }
    }
    // Invalid cases, including pseudo-keys, must leave all header bytes intact.
    for k in [-1, 0, 81, 1000, 1234, i32::MIN, i32::MAX] {
        for _ in 0..2048 {
            tr.hdr = random_header(&mut rng);
            check_write(&mut tr, k, f64::from_bits(rng.next_u64()));
        }
    }
    for (i, &x) in tr.data.iter().enumerate() {
        assert_eq!(x.to_bits(), ((i as i32 - 32768) as f32).to_bits());
    }

    // The public API accepts a core Trace through its header dereference.
    let mut trace = oxwave_core::trace::Trace::new();
    rust::tohead(&mut trace, 1, 2.5);
    assert_eq!(rust::fromhead(&trace, 1), 2.0);
}
