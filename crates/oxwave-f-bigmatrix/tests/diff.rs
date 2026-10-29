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

//! Differential tests for cwp/lib/bigmatrix.c.
//! Author: Dave Hale, Colorado School of Mines, 05/17/89.

use oxwave_cref::XorShift;
use oxwave_f_bigmatrix as rust;
use std::ffi::{c_int, c_void};

extern "C" {
    fn bmalloc(nbpe: c_int, n1: c_int, n2: c_int) -> *mut c_void;
    fn bmdump(bm: *mut c_void);
    fn bmfree(bm: *mut c_void);
    fn bmread(bm: *mut c_void, d: c_int, k1: c_int, k2: c_int, n: c_int, v: *mut c_void);
    fn bmwrite(bm: *mut c_void, d: c_int, k1: c_int, k2: c_int, n: c_int, v: *mut c_void);
    fn fflush(stream: *mut c_void) -> c_int;
}

#[repr(C)]
#[derive(Clone, Copy)]
struct CState {
    l1: i32,
    l2: i32,
    m1: i32,
    m2: i32,
    nbpe: i32,
    p: *mut u8,
}

fn dimensions(r: &mut XorShift, i: usize) -> (usize, usize, usize) {
    let dims = [1, 2, 3, 7, 16, 31, 32, 33, 63, 97, 128, 257, 1024];
    let widths = [1, 2, 3, 4, 5, 8, 12, 16];
    let edges = [
        (1, 1),
        (2, 2),
        (257, 263),
        (1024, 1024),
        (1, 4093),
        (4093, 1),
        (1025, 33),
        (33, 1025),
    ];
    let (n1, n2) = if i < edges.len() {
        edges[i]
    } else if i % 11 == 0 {
        (r.range_i32(1, 513) as usize, r.range_i32(1, 513) as usize)
    } else {
        (
            dims[r.next_u32() as usize % dims.len()],
            dims[r.next_u32() as usize % dims.len()],
        )
    };
    (n1, n2, widths[i % widths.len()])
}

fn bytes(r: &mut XorShift, len: usize, case: usize) -> Vec<u8> {
    let mut data: Vec<_> = (0..len).map(|_| r.next_u32() as u8).collect();
    match case % 11 {
        0 => data.fill(0),
        1 => data.fill(0xff), // Signed -1 and NaN bit patterns.
        2 => {
            // Native f32 bits: negative numbers, signed zero, subnormal,
            // infinities and quiet/signaling NaNs. Every bit must be preserved.
            let words = [
                0xbf800000u32,
                0x80000000,
                0,
                1,
                0xff800000,
                0x7f800000,
                0x7fc12345,
                0x7f812345,
            ];
            let pattern: Vec<_> = words.iter().flat_map(|x| x.to_ne_bytes()).collect();
            for (i, byte) in data.iter_mut().enumerate() {
                *byte = pattern[i % pattern.len()];
            }
        }
        _ => {}
    }
    data
}

// Independent address calculation by page number and position within the page.
fn offset(s: CState, k1: usize, k2: usize) -> usize {
    let (l1, l2, m1, nb) = (s.l1 as usize, s.l2 as usize, s.m1 as usize, s.nbpe as usize);
    let page = (k2 / l2) * m1 + k1 / l1;
    let within = (k2 % l2) * l1 + k1 % l1;
    (page * l1 * l2 + within) * nb
}

fn check_read(
    c: *mut c_void,
    matrix: &rust::BigMatrix,
    expected: &[u8],
    n1: usize,
    nb: usize,
    d: i32,
    k1: usize,
    k2: usize,
    n: usize,
) {
    // Eight-byte guards preserve alignment for C's typed branches.
    let mut cv = vec![0xa5; n * nb + 16];
    let mut rv = cv.clone();
    unsafe {
        bmread(
            c,
            d,
            k1 as i32,
            k2 as i32,
            n as i32,
            cv.as_mut_ptr().add(8).cast(),
        )
    };
    rust::bmread(matrix, d, k1, k2, &mut rv[8..8 + n * nb]);
    let mut want = vec![0xa5; n * nb + 16];
    for i in 0..n {
        let (x, y) = if d == 1 { (k1 + i, k2) } else { (k1, k2 + i) };
        let start = (y * n1 + x) * nb;
        want[8 + i * nb..8 + (i + 1) * nb].copy_from_slice(&expected[start..start + nb]);
    }
    assert_eq!(cv, want, "C read: d={d}, k1={k1}, k2={k2}, n={n}, nb={nb}");
    assert_eq!(
        rv, want,
        "Rust read: d={d}, k1={k1}, k2={k2}, n={n}, nb={nb}"
    );
}

fn exercise(seed: u64, writes: bool) {
    let mut r = XorShift::new(seed);
    for case in 0..2500 {
        let (n1, n2, nb) = dimensions(&mut r, case);
        let c = unsafe { bmalloc(nb as i32, n1 as i32, n2 as i32) };
        assert!(!c.is_null());
        let state = unsafe { *(c as *const CState) };
        let mut matrix = rust::bmalloc(nb as i32, n1 as i32, n2 as i32).unwrap();
        let mut expected = bytes(&mut r, n1 * n2 * nb, case);
        // Initialize C storage directly, independently of both bmwrite ports.
        // Never read the indeterminate padding left by C malloc.
        for y in 0..n2 {
            let row = &expected[y * n1 * nb..(y + 1) * n1 * nb];
            for x in 0..n1 {
                unsafe {
                    std::ptr::copy_nonoverlapping(
                        row.as_ptr().add(x * nb),
                        state.p.add(offset(state, x, y)),
                        nb,
                    );
                }
            }
            rust::bmwrite(&mut matrix, 1, 0, y, row);
        }
        for op in 0..8 {
            let d = [1, 2, 0, -1, i32::MIN, i32::MAX][(case + op) % 6];
            let (k1, k2) = match op {
                0 | 2 => (0, 0),
                1 => (n1 - 1, n2 - 1),
                _ => (r.next_u32() as usize % n1, r.next_u32() as usize % n2),
            };
            let remaining = if d == 1 { n1 - k1 } else { n2 - k2 };
            let n = match op {
                0 => 0,
                1 => 1,
                2 => remaining,
                _ => r.next_u32() as usize % (remaining + 1),
            };
            if writes {
                let data = bytes(&mut r, n * nb, case + op);
                let mut cdata = data.clone();
                unsafe {
                    bmwrite(
                        c,
                        d,
                        k1 as i32,
                        k2 as i32,
                        n as i32,
                        cdata.as_mut_ptr().cast(),
                    )
                };
                rust::bmwrite(&mut matrix, d, k1, k2, &data);
                assert_eq!(cdata, data, "bmwrite must preserve its input");
                for i in 0..n {
                    let (x, y) = if d == 1 { (k1 + i, k2) } else { (k1, k2 + i) };
                    let start = (y * n1 + x) * nb;
                    expected[start..start + nb].copy_from_slice(&data[i * nb..(i + 1) * nb]);
                }
            }
            check_read(c, &matrix, &expected, n1, nb, d, k1, k2, n);
            let other = if d == 1 { 2 } else { 1 };
            let len = if other == 1 { n1 - k1 } else { n2 - k2 };
            check_read(c, &matrix, &expected, n1, nb, other, k1, k2, len);
        }
        if writes {
            // Check every coordinate, including elements outside the written
            // vectors, against C's raw storage and an independent flat oracle.
            let mut row = vec![0u8; n1 * nb];
            for y in 0..n2 {
                rust::bmread(&matrix, 1, 0, y, &mut row);
                assert_eq!(row, expected[y * n1 * nb..(y + 1) * n1 * nb]);
                for x in 0..n1 {
                    let actual =
                        unsafe { std::slice::from_raw_parts(state.p.add(offset(state, x, y)), nb) };
                    assert_eq!(
                        actual,
                        &row[x * nb..(x + 1) * nb],
                        "case={case}, x={x}, y={y}"
                    );
                }
            }
        }
        unsafe { bmfree(c) };
        rust::bmfree(matrix);
    }
}

#[test]
fn bmread_matches_c() {
    exercise(123, false);
}
#[test]
fn bmwrite_matches_c() {
    exercise(456, true);
}
#[test]
fn bmfree_matches_c() {
    let mut r = XorShift::new(789);
    for i in 0..2500 {
        let (n1, n2, nb) = dimensions(&mut r, i);
        let c = unsafe { bmalloc(nb as i32, n1 as i32, n2 as i32) };
        assert!(!c.is_null());
        let matrix = rust::bmalloc(nb as i32, n1 as i32, n2 as i32).unwrap();
        // Free one matrix while another remains live, then verify and free it.
        let survivor = unsafe { bmalloc(nb as i32, 2, 3) };
        assert!(!survivor.is_null());
        let mut survivor_r = rust::bmalloc(nb as i32, 2, 3).unwrap();
        let data = bytes(&mut r, 6 * nb, i);
        for y in 0..3 {
            let mut row = data[y * 2 * nb..(y + 1) * 2 * nb].to_vec();
            unsafe { bmwrite(survivor, 1, 0, y as i32, 2, row.as_mut_ptr().cast()) };
            rust::bmwrite(&mut survivor_r, 1, 0, y, &row);
        }
        unsafe { bmfree(c) };
        rust::bmfree(matrix);
        check_read(survivor, &survivor_r, &data, 2, nb, 2, 1, 0, 3);
        unsafe { bmfree(survivor) };
        rust::bmfree(survivor_r);
    }
}

// The diagnostic exposes the allocation layout, so both allocation and dump
// tests compare all its fields against C. Isolate stdout in a child process.
fn check_dumps(test_name: &str, seed: u64) {
    if std::env::var_os("OXWAVE_DUMP_CHILD").is_some() {
        let mut r = XorShift::new(seed);
        for i in 0..2500 {
            let (n1, n2, nb) = dimensions(&mut r, i);
            let c = unsafe { bmalloc(nb as i32, n1 as i32, n2 as i32) };
            assert!(!c.is_null());
            let matrix = rust::bmalloc(nb as i32, n1 as i32, n2 as i32).unwrap();
            let s = unsafe { &*(c as *const CState) };
            assert_eq!(s.nbpe, nb as i32);
            assert_eq!(s.l1 * s.l2, 1024);
            assert_eq!(s.m1, 1 + (n1 as i32 - 1) / s.l1);
            assert_eq!(s.m2, 1 + (n2 as i32 - 1) / s.l2);
            unsafe {
                bmdump(c);
                assert_eq!(fflush(std::ptr::null_mut()), 0);
            }
            rust::bmdump(&matrix);
            assert_eq!(oxwave_core::cio::fflush(&oxwave_core::cio::stdout()), 0);
            unsafe { bmfree(c) };
            rust::bmfree(matrix);
        }
        return;
    }
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", test_name, "--nocapture"])
        .env("OXWAVE_DUMP_CHILD", "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<_> = text.lines().collect();
    let mut count = 0;
    for (i, line) in lines.iter().enumerate() {
        if let Some(addr) = line.strip_prefix("big matrix state ") {
            assert!(addr.strip_suffix(':').unwrap().parse::<usize>().unwrap() != 0);
            assert!(
                lines[i + 6]
                    .strip_prefix("p = ")
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
                    != 0
            );
            if count % 2 == 0 {
                assert!(lines[i + 7].starts_with("big matrix state "));
                assert_eq!(&lines[i + 1..i + 6], &lines[i + 8..i + 13]);
            }
            count += 1;
        }
    }
    assert_eq!(count, 5000);
}

#[test]
fn bmalloc_matches_c() {
    check_dumps("bmalloc_matches_c", 42);
    // Outside the positive-size Rust contract; do not invoke undefined C accesses.
    for (nb, n1, n2) in [
        (0, 1, 1),
        (-1, 1, 1),
        (1, 0, 1),
        (1, 1, -1),
        (i32::MAX, 1, 1),
    ] {
        assert!(rust::bmalloc(nb, n1, n2).is_none());
    }
}

#[test]
fn bmdump_matches_c() {
    check_dumps("bmdump_matches_c", 987);
}
