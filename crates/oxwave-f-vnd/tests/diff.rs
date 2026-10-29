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

//! Differential tests for par/lib/VND.c.
//! Credits: J.E. Anderson (CSM Mobil Visiting Scientist, 1993).
use oxwave_cref::XorShift;
use oxwave_f_vnd::*;
use std::ffi::c_void;

#[repr(C)]
struct CVnd {
    num_dim: i64,
    bytes_node: i64,
    blocks_panel: i64,
    files: i64,
    panels: i64,
    open_panels: i64,
    bytes_mem: i64,
    bytes_block: i64,
    bytes_first: i64,
    nodes_block: i64,
    blocks_file: i64,
    current_file: i32,
    mode: *mut c_void,
    fd: *mut c_void,
    names_len: *mut c_void,
    names: *mut c_void,
    first_data: *mut c_void,
    file_len: *mut c_void,
    n: *mut i64,
    nodes: *mut i64,
    blocks: *mut i64,
    byte_incr: *mut i64,
    block_incr: *mut i64,
    current_panel: *mut i64,
    current_dimension: *mut i32,
    current_block: *mut i64,
    modified: *mut i32,
    w: *mut u8,
}
mod c {
    use super::*;
    extern "C" {
        pub fn VNDGetDimensions(v: *mut CVnd, max: i64, ich: i32) -> i32;
        pub fn VNDr2c(v: *mut CVnd);
        pub fn VNDc2r(v: *mut CVnd);
        pub fn vndr2c_(v: *mut CVnd);
        pub fn vndc2r_(v: *mut CVnd);
        pub fn VNDemalloc(n: usize, msg: *const i8) -> *mut u8;
        pub fn VNDmemchk(p: *mut u8, msg: *const i8);
        pub fn VNDfree(p: *mut u8, msg: *const i8);
        pub fn VNDtotalmem() -> i64;
        pub fn VNDtempname(root: *const i8) -> *mut i8;
        pub fn VNDerr(msg: *const i8);
        pub fn VNDlook(w: *const i32, lw: i32);
        pub fn VNDop(
            mode: i32,
            lwmax: i64,
            ndim: i32,
            n: *const i64,
            npanels: i64,
            nbytes: i64,
            file: *const i8,
            ndir: i32,
            dir: *const *const i8,
            noppanels: i32,
        ) -> *mut CVnd;
        pub fn V2Dop(
            mode: i32,
            lwmax: i64,
            nbytes: i64,
            file: *const i8,
            n0: i64,
            n1: i64,
        ) -> *mut CVnd;
        pub fn VNDrw(
            iop: i8,
            iopanel: i32,
            v: *mut CVnd,
            idim: i32,
            key: *const i64,
            ipanel: i64,
            values: *mut u8,
            start: i64,
            incr: i64,
            nvalues: i64,
            msgnum: i32,
            msg: *const i8,
        );
        pub fn VNDrwslice(
            v: *mut CVnd,
            iopanel: i32,
            mode: i8,
            idim: i32,
            iblock: i64,
            ipanel: i64,
        ) -> i32;
        pub fn VNDseek(v: *mut CVnd, ib: i64) -> i32;
        pub fn VNDflush(v: *mut CVnd);
        pub fn vndflush_(v: *mut *mut CVnd);
        pub fn VNDcl(v: *mut CVnd, mode: i32);
        pub fn vndcl_(v: *mut *mut CVnd, mode: *const i64);
        pub fn VNDdump(v: *mut CVnd, fp: *mut c_void);
        pub fn V2Dr0(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i32);
        pub fn V2Dr1(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i32);
        pub fn V2Dw0(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i32);
        pub fn V2Dw1(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i32);
        pub fn V2Dr1_(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i64);
        pub fn v2dr2_(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i64);
        pub fn v2dw1_(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i64);
        pub fn v2dw2_(v: *mut CVnd, key: i64, buf: *mut u8, msgnum: i64);
        pub fn fflush(fp: *mut c_void) -> i32;
        pub fn fclose(fp: *mut c_void) -> i32;
        pub fn tmpfile() -> *mut c_void;
        pub fn fseeko(fp: *mut c_void, offset: i64, whence: i32) -> i32;
        pub fn ftello(fp: *mut c_void) -> i64;
        pub fn fread(buf: *mut u8, size: usize, n: usize, fp: *mut c_void) -> usize;
    }
}

fn compare(which: usize) {
    let _guard = MEMORY_TEST.lock().unwrap();
    let mut rng = XorShift::new(91 + which as u64);
    for iteration in 0..4000 {
        let ndim = 1 + (rng.next_u32() % 4) as usize;
        let mut n = (0..ndim)
            .map(|_| 1 + (rng.next_u32() % 37) as i64)
            .collect::<Vec<_>>();
        if iteration % 7 == 0 {
            n[0] = [1, 2, 17, 32, 101][iteration % 5];
        }
        let mut bytes = 1 + (rng.next_u32() % 16) as i64;
        let mut max = (rng.next_u32() % 200_000) as i64;
        // Geometry arithmetic is defined for these degenerate signed inputs;
        // force the nonpartitioning branch so no allocation sizes are involved.
        if iteration % 8 == 0 {
            n[0] = [0, -1, -2, -17][iteration / 8 % 4];
            bytes = [0, -1, -8, 4][iteration / 8 % 4];
            max = 1_000_000_000;
        }
        let ich = (rng.next_u32() % 2) as i32;
        let mut geom = VNDGeometry::new(&n, bytes);
        let mut nodes = vec![0; ndim];
        let mut blocks = vec![0; ndim];
        let mut byte_incr = vec![0; ndim];
        let mut block_incr = vec![0; ndim];
        // The trailing fields are not accessed by these routines.
        let mut c: CVnd = unsafe { std::mem::zeroed() };
        c.num_dim = ndim as i64;
        c.bytes_node = bytes;
        c.n = n.as_mut_ptr();
        c.nodes = nodes.as_mut_ptr();
        c.blocks = blocks.as_mut_ptr();
        c.byte_incr = byte_incr.as_mut_ptr();
        c.block_incr = block_incr.as_mut_ptr();
        let expected = unsafe { c::VNDGetDimensions(&mut c, max, ich) };
        assert_eq!(VNDGetDimensions(&mut geom, max, ich), expected);
        match which {
            1 => {
                VNDr2c(&mut geom);
                unsafe { c::VNDr2c(&mut c) };
            }
            2 => {
                VNDc2r(&mut geom);
                unsafe { c::VNDc2r(&mut c) };
            }
            3 => {
                vndr2c_(&mut geom);
                unsafe { c::vndr2c_(&mut c) };
            }
            4 => {
                vndc2r_(&mut geom);
                unsafe { c::vndc2r_(&mut c) };
            }
            _ => {}
        }
        assert_eq!(geom.n, n);
        assert_eq!(geom.num_bytes_per_node, c.bytes_node);
        assert_eq!(geom.num_blocks_per_panel, c.blocks_panel);
        assert_eq!(geom.num_bytes_mem_buf, c.bytes_mem);
        assert_eq!(geom.num_bytes_per_block, c.bytes_block);
        assert_eq!(geom.num_bytes_1st_dim, c.bytes_first);
        assert_eq!(geom.num_nodes_per_block, c.nodes_block);
        assert_eq!(geom.nodes_per_block, nodes);
        assert_eq!(geom.blocks, blocks);
        assert_eq!(geom.byte_incr, byte_incr);
        assert_eq!(geom.block_incr, block_incr);
    }
}
#[test]
fn dimensions() {
    compare(0);
}
#[test]
fn r2c() {
    compare(1);
}
#[test]
fn c2r() {
    compare(2);
}
#[test]
fn fortran_r2c() {
    compare(3);
}
#[test]
fn fortran_c2r() {
    compare(4);
}

static MEMORY_TEST: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn compare_memory(seed: u64) {
    let _guard = MEMORY_TEST.lock().unwrap();
    let mut rng = XorShift::new(seed);
    let msg = b"differential test\0";
    for i in 0..4000 {
        let n = if i < 8 {
            [0, 1, 2, 7, 8, 9, 17, 4096][i]
        } else {
            1 + (rng.next_u32() as usize % 4096)
        };
        unsafe {
            let before_c = c::VNDtotalmem();
            let before_r = VNDtotalmem();
            let c = c::VNDemalloc(n, msg.as_ptr().cast());
            let r = VNDemalloc(n, "differential test");
            assert_eq!(c::VNDtotalmem() - before_c, VNDtotalmem() - before_r);
            for j in 0..n {
                let value = rng.next_u32() as u8;
                c.add(j).write(value);
                r.add(j).write(value);
            }
            c::VNDmemchk(c, msg.as_ptr().cast());
            VNDmemchk(r, "differential test");
            assert_eq!(
                std::slice::from_raw_parts(c, n),
                std::slice::from_raw_parts(r, n)
            );
            c::VNDfree(c, msg.as_ptr().cast());
            VNDfree(r, "differential test");
            assert_eq!(c::VNDtotalmem(), before_c);
            assert_eq!(VNDtotalmem(), before_r);
        }
    }
}

#[test]
fn emalloc() {
    compare_memory(101);
}
#[test]
fn memchk() {
    compare_memory(102);
}
#[test]
fn free() {
    compare_memory(103);
}
#[test]
fn totalmem() {
    compare_memory(104);
}

#[test]
fn tempname() {
    let _guard = MEMORY_TEST.lock().unwrap();
    let mut rng = XorShift::new(105);
    for i in 0..4000 {
        let maxroot = 79 - std::process::id().to_string().len() - (i + 1).to_string().len();
        let length = [0, 1, 2, 7, 31, 63, maxroot][i % 7].min(maxroot);
        let root: String = (0..length)
            .map(|_| (b'a' + (rng.next_u32() % 26) as u8) as char)
            .collect();
        let croot = std::ffi::CString::new(root.clone()).unwrap();
        unsafe {
            let c = c::VNDtempname(croot.as_ptr());
            let r = VNDtempname(&root);
            assert_eq!(
                std::ffi::CStr::from_ptr(c).to_bytes(),
                std::ffi::CStr::from_ptr(r.cast()).to_bytes()
            );
            c::VNDfree(c.cast(), b"test\0".as_ptr().cast());
            VNDfree(r, "test");
        }
    }
}

fn workspace() -> std::path::PathBuf {
    let p = std::path::PathBuf::from("../../.targets/lib-f_vnd/vnd-diff");
    std::fs::create_dir_all(&p).unwrap();
    p
}
unsafe fn assert_state(r: &VND, c: *mut CVnd) {
    let c = &*c;
    assert_eq!(r.n.len() as i64, c.num_dim);
    assert_eq!(r.num_bytes_per_node, c.bytes_node);
    assert_eq!(r.num_blocks_per_panel, c.blocks_panel);
    assert_eq!(r.num_bytes_mem_buf, c.bytes_mem);
    assert_eq!(r.num_bytes_per_block, c.bytes_block);
    assert_eq!(r.num_bytes_1st_dim, c.bytes_first);
    assert_eq!(r.num_nodes_per_block, c.nodes_block);
    assert_eq!(r.num_blocks_per_file, c.blocks_file);
    assert_eq!(r.num_files, c.files);
    assert_eq!(r.num_panels, c.panels);
    assert_eq!(r.num_open_panels, c.open_panels);
    assert_eq!(r.current_file, c.current_file);
    let nd = r.n.len();
    let np = r.num_open_panels as usize;
    for (a, p) in [
        (&r.n, c.n),
        (&r.nodes_per_block, c.nodes),
        (&r.blocks, c.blocks),
        (&r.byte_incr, c.byte_incr),
        (&r.block_incr, c.block_incr),
    ] {
        assert_eq!(a.as_slice(), std::slice::from_raw_parts(p, nd));
    }
    assert_eq!(
        r.first_data_byte,
        std::slice::from_raw_parts(c.first_data.cast::<i64>(), r.num_files as usize)
    );
    assert_eq!(
        r.len_file,
        std::slice::from_raw_parts(c.file_len.cast::<i64>(), r.num_files as usize)
    );
    assert_eq!(
        r.current_panel,
        std::slice::from_raw_parts(c.current_panel, np)
    );
    assert_eq!(
        r.current_dimension,
        std::slice::from_raw_parts(c.current_dimension, np)
    );
    assert_eq!(
        r.current_block,
        std::slice::from_raw_parts(c.current_block, np)
    );
    assert_eq!(r.modified, std::slice::from_raw_parts(c.modified, np));
    for j in 0..np {
        // C initializes only the first panel buffer at allocation time.
        if r.current_panel[j] >= 0 || j == 0 {
            let off = j * r.num_bytes_mem_buf as usize;
            let len = if r.current_dimension[j] < 0 {
                r.num_bytes_mem_buf as usize
            } else {
                (r.blocks[r.current_dimension[j] as usize] * r.num_bytes_per_block) as usize
            };
            let expected = std::slice::from_raw_parts(c.w.add(off), len);
            let difference = r.w[off..off + len]
                .iter()
                .zip(expected)
                .position(|(a, b)| a != b);
            assert!(difference.is_none(),"cache buffer {j} differs at {difference:?}, geometry {:?}, panel {:?}, dims {:?}, blocks {:?}",r.geometry,r.current_panel,r.current_dimension,r.current_block);
        }
    }
}
unsafe fn read_stream(fp: *mut c_void) -> Vec<u8> {
    assert_eq!(c::fflush(fp), 0);
    let size = c::ftello(fp);
    assert!(size >= 0);
    assert_eq!(c::fseeko(fp, 0, 0), 0);
    let mut bytes = vec![0; size as usize];
    assert_eq!(
        c::fread(bytes.as_mut_ptr(), 1, bytes.len(), fp),
        bytes.len()
    );
    bytes
}
unsafe fn dump_pair(r: &VND, c: *mut CVnd, cname: &str, rname: &str) {
    let fp = c::tmpfile();
    assert!(!fp.is_null());
    c::VNDdump(c, fp);
    let text = String::from_utf8(read_stream(fp)).unwrap();
    assert_eq!(c::fclose(fp), 0);
    let mut mine = Vec::new();
    VNDdump(r, &mut mine);
    assert_eq!(String::from_utf8(mine).unwrap(), text.replace(cname, rname));
}
fn random_bytes(rng: &mut XorShift, len: usize, zero: bool) -> Vec<u8> {
    (0..len)
        .map(|_| if zero { 0 } else { rng.next_u32() as u8 })
        .collect()
}

// Every entry invokes its named API on 2,048 independently generated shapes.
// Byte equality covers f32/f64, signed integers, NaNs, infinities and zeroes
// without introducing a numeric tolerance into an opaque byte storage API.
fn engine_case(which: &str) {
    let _guard = MEMORY_TEST.lock().unwrap();
    let root = workspace();
    let mut rng = XorShift::new(0x93ac712 ^ which.bytes().map(u64::from).sum::<u64>());
    for iteration in 0..2048 {
        let wrappers = which.starts_with("V2D") || which.starts_with("v2d");
        let ndim = if wrappers {
            2
        } else {
            1 + (rng.next_u32() % 4) as usize
        };
        let mut n: Vec<i64> = (0..ndim)
            .map(|_| [1, 2, 7, 16, 31, 37][rng.next_u32() as usize % 6])
            .collect();
        if iteration % 128 == 0 {
            n[0] = [127, 256, 257, 1024][iteration / 128 % 4];
        }
        let bytes = [1, 2, 3, 4, 8, 16][rng.next_u32() as usize % 6];
        let max = [-1, 0, 1024, 4096, 8192, 65536, 1_000_000][iteration % 7];
        let panels = if wrappers {
            1
        } else {
            1 + (rng.next_u32() % 3) as i64
        };
        let opens = if wrappers {
            1
        } else {
            1 + (rng.next_u32() % 2) as i32
        };
        let cname = root.join(format!("{which}-c"));
        let rname = root.join(format!("{which}-r"));
        let cname = cname.to_str().unwrap();
        let rname = rname.to_str().unwrap();
        assert!(cname.len() + 6 < 120);
        let name = std::ffi::CString::new(cname).unwrap();
        let before_c = unsafe { c::VNDtotalmem() };
        let before_r = VNDtotalmem();
        let mut c = unsafe {
            if which == "V2Dop" {
                c::V2Dop(2, max, bytes, name.as_ptr(), n[0], n[1])
            } else {
                c::VNDop(
                    2,
                    max,
                    ndim as i32,
                    n.as_ptr(),
                    panels,
                    bytes,
                    name.as_ptr(),
                    0,
                    std::ptr::null(),
                    opens,
                )
            }
        };
        assert!(!c.is_null());
        let mut r = if which == "V2Dop" {
            V2Dop(2, max, bytes, rname, n[0], n[1]).unwrap()
        } else {
            VNDop(2, max, &n, panels, bytes, rname, 0, &[], opens).unwrap()
        };
        unsafe {
            assert_state(&r, c);
        }
        assert_eq!(
            VNDtotalmem() - before_r,
            unsafe { c::VNDtotalmem() } - before_c
        );
        for op in 0..5 {
            let d =
                if which == "V2Dr0" || which == "V2Dr1_" || which == "V2Dw0" || which == "v2dw1_" {
                    0
                } else if wrappers && which != "V2Dop" {
                    1
                } else {
                    rng.next_u32() as usize % ndim
                };
            let key: Vec<i64> = n.iter().map(|&x| rng.next_u32() as i64 % x).collect();
            let p = if wrappers {
                0
            } else {
                rng.next_u32() as i64 % panels
            };
            // Different buffers own different panels; never modify the same
            // panel through two buffers (explicitly forbidden by VND.h).
            let opened = if opens > 1 {
                (p % opens as i64) as i32
            } else {
                0
            };
            let mut data = random_bytes(&mut rng, (n[d] * bytes) as usize, iteration % 11 == 0);
            let mut cd = data.clone();
            unsafe {
                match which {
                    "V2Dw0" => {
                        V2Dw0(&mut r, key[1], &mut data, op);
                        c::V2Dw0(c, key[1], cd.as_mut_ptr(), op);
                    }
                    "V2Dw1" => {
                        V2Dw1(&mut r, key[0], &mut data, op);
                        c::V2Dw1(c, key[0], cd.as_mut_ptr(), op);
                    }
                    "v2dw1_" => {
                        v2dw1_(&mut r, key[1] + 1, &mut data, op as i64);
                        c::v2dw1_(c, key[1] + 1, cd.as_mut_ptr(), op as i64);
                    }
                    "v2dw2_" => {
                        v2dw2_(&mut r, key[0] + 1, &mut data, op as i64);
                        c::v2dw2_(c, key[0] + 1, cd.as_mut_ptr(), op as i64);
                    }
                    _ => {
                        VNDrw(
                            b'w', opened, &mut r, d as i32, &key, p, &mut data, 0, 1, n[d], op,
                            "diff",
                        );
                        c::VNDrw(
                            b'w' as i8,
                            opened,
                            c,
                            d as i32,
                            key.as_ptr(),
                            p,
                            cd.as_mut_ptr(),
                            0,
                            1,
                            n[d],
                            op,
                            b"diff\0".as_ptr().cast(),
                        );
                    }
                }
                assert_state(&r, c);
                data.fill(0xa5);
                cd.fill(0x5a);
                match which {
                    "V2Dr0" => {
                        V2Dr0(&mut r, key[1], &mut data, op);
                        c::V2Dr0(c, key[1], cd.as_mut_ptr(), op);
                    }
                    "V2Dr1" => {
                        V2Dr1(&mut r, key[0], &mut data, op);
                        c::V2Dr1(c, key[0], cd.as_mut_ptr(), op);
                    }
                    "V2Dr1_" => {
                        V2Dr1_(&mut r, key[1] + 1, &mut data, op as i64);
                        c::V2Dr1_(c, key[1] + 1, cd.as_mut_ptr(), op as i64);
                    }
                    "v2dr2_" => {
                        v2dr2_(&mut r, key[0] + 1, &mut data, op as i64);
                        c::v2dr2_(c, key[0] + 1, cd.as_mut_ptr(), op as i64);
                    }
                    _ => {
                        VNDrw(
                            b'r', opened, &mut r, d as i32, &key, p, &mut data, 0, 1, n[d], op,
                            "diff",
                        );
                        c::VNDrw(
                            b'r' as i8,
                            opened,
                            c,
                            d as i32,
                            key.as_ptr(),
                            p,
                            cd.as_mut_ptr(),
                            0,
                            1,
                            n[d],
                            op,
                            b"diff\0".as_ptr().cast(),
                        );
                    }
                }
                assert_eq!(data, cd, "{which} iteration {iteration}");
                // Exercise both signed strides and the zero-element case.
                if which == "VNDrw" {
                    let stride: i64 = if op % 2 == 0 { -2 } else { 3 };
                    let start = if stride < 0 { n[d] - 1 } else { 0 };
                    let count = if op == 0 {
                        0
                    } else {
                        1 + (n[d] - 1) / stride.abs()
                    };
                    let len = (count * bytes) as usize;
                    let mut data = random_bytes(&mut rng, len, false);
                    let mut cd = data.clone();
                    VNDrw(
                        b'w', opened, &mut r, d as i32, &key, p, &mut data, start, stride, count,
                        op, "stride",
                    );
                    c::VNDrw(
                        b'w' as i8,
                        opened,
                        c,
                        d as i32,
                        key.as_ptr(),
                        p,
                        cd.as_mut_ptr(),
                        start,
                        stride,
                        count,
                        op,
                        b"stride\0".as_ptr().cast(),
                    );
                    data.fill(0xa5);
                    cd.fill(0x5a);
                    VNDrw(
                        b'r', opened, &mut r, d as i32, &key, p, &mut data, start, stride, count,
                        op, "stride",
                    );
                    c::VNDrw(
                        b'r' as i8,
                        opened,
                        c,
                        d as i32,
                        key.as_ptr(),
                        p,
                        cd.as_mut_ptr(),
                        start,
                        stride,
                        count,
                        op,
                        b"stride\0".as_ptr().cast(),
                    );
                    assert_eq!(data, cd);
                }
                if which == "vndflush_" {
                    vndflush_(&mut r);
                    c::vndflush_(&mut c);
                } else {
                    VNDflush(&mut r);
                    c::VNDflush(c);
                }
                assert_state(&r, c);
                if which == "VNDseek" {
                    let ib = rng.next_u32() as i64 % (r.num_blocks_per_panel * panels);
                    assert_eq!(VNDseek(&mut r, ib), c::VNDseek(c, ib));
                    assert_eq!(
                        c::ftello((*c).fd),
                        r.first_data_byte[0] + ib * r.num_bytes_per_block
                    );
                }
                if which == "VNDrwslice" {
                    let iblock: i64 = (0..ndim)
                        .filter(|&j| j != d)
                        .map(|j| key[j] / r.nodes_per_block[j] * r.block_incr[j])
                        .sum();
                    assert_eq!(
                        VNDrwslice(&mut r, opened, b'r', d as i32, iblock, p),
                        c::VNDrwslice(c, opened, b'r' as i8, d as i32, iblock, p)
                    );
                    assert_state(&r, c);
                    let base = opened as usize * r.num_bytes_mem_buf as usize;
                    let count = (r.blocks[d] * r.num_bytes_per_block) as usize;
                    let pattern = random_bytes(&mut rng, count, op == 0);
                    r.w[base..base + count].copy_from_slice(&pattern);
                    std::ptr::copy_nonoverlapping(pattern.as_ptr(), (*c).w.add(base), count);
                    assert_eq!(
                        VNDrwslice(&mut r, opened, b'w', d as i32, iblock, p),
                        c::VNDrwslice(c, opened, b'w' as i8, d as i32, iblock, p)
                    );
                }
                if which == "VNDdump" {
                    dump_pair(&r, c, cname, rname);
                }
            }
        }
        unsafe {
            c::fflush((*c).fd);
        }
        let rd = std::fs::read(format!("{rname}.VND0")).unwrap();
        let cd = std::fs::read(format!("{cname}.VND0")).unwrap();
        let offset = r.first_data_byte[0] as usize;
        assert_eq!(
            &rd[offset..],
            &cd[offset..],
            "disk payload: {which} iteration {iteration}"
        );
        // The scalar header and geometry arrays are persistent and exact;
        // process-local pointer slots and the two filenames differ by design.
        assert_eq!(&rd[..88], &cd[..88]);
        assert_eq!(&rd[224..224 + ndim * 40], &cd[224..224 + ndim * 40]);
        let keep = (which == "VNDcl" || which == "vndcl_") && iteration % 2 == 0;
        let mode = if keep { 0 } else { 1 };
        if which == "vndcl_" {
            vndcl_(r, mode as i64);
            unsafe {
                c::vndcl_(&mut c, &(mode as i64));
            }
        } else {
            VNDcl(r, mode);
            unsafe {
                c::VNDcl(c, mode);
            }
        }
        assert_eq!(
            std::path::Path::new(&format!("{rname}.VND0")).exists(),
            keep
        );
        assert_eq!(
            std::path::Path::new(&format!("{cname}.VND0")).exists(),
            keep
        );
        let leaked = if keep {
            (1 + (cname.len() + 6 - 1) / 8 + 3) * 8
        } else {
            0
        } as i64;
        // VNDcl(mode=0) leaks its filename allocation in the original C.
        assert_eq!(VNDtotalmem() - before_r, leaked);
        assert_eq!(unsafe { c::VNDtotalmem() } - before_c, leaked);
        if keep {
            std::fs::remove_file(format!("{rname}.VND0")).unwrap();
            std::fs::remove_file(format!("{cname}.VND0")).unwrap();
        }
    }
}

#[test]
fn vndop() {
    engine_case("VNDop");
}
#[test]
fn v2dop() {
    engine_case("V2Dop");
}
#[test]
fn v2dr0() {
    engine_case("V2Dr0");
}
#[test]
fn v2dr1() {
    engine_case("V2Dr1");
}
#[test]
fn v2dw0() {
    engine_case("V2Dw0");
}
#[test]
fn v2dw1() {
    engine_case("V2Dw1");
}
#[test]
fn v2dr1_fortran() {
    engine_case("V2Dr1_");
}
#[test]
fn v2dr2_fortran() {
    engine_case("v2dr2_");
}
#[test]
fn v2dw1_fortran() {
    engine_case("v2dw1_");
}
#[test]
fn v2dw2_fortran() {
    engine_case("v2dw2_");
}
#[test]
fn vndrw() {
    engine_case("VNDrw");
}
#[test]
fn vndrwslice() {
    engine_case("VNDrwslice");
}
#[test]
fn vndseek() {
    engine_case("VNDseek");
}
#[test]
fn vnddump() {
    engine_case("VNDdump");
}
#[test]
fn vndflush() {
    engine_case("VNDflush");
}
#[test]
fn vndflush_fortran() {
    engine_case("vndflush_");
}
#[test]
fn vndcl() {
    engine_case("VNDcl");
}
#[test]
fn vndcl_fortran() {
    engine_case("vndcl_");
}

extern "C" {
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn _exit(status: i32) -> !;
}
extern "C" fn caught_abort(_: i32) {
    unsafe {
        _exit(134);
    }
}

#[test]
fn vnderr() {
    if let Ok(message) = std::env::var("VND_ERROR_MESSAGE") {
        // Catch SIGABRT only to avoid OS crash-report/core-file side effects;
        // exit 134 proves the C and Rust routines both raised SIGABRT.
        unsafe {
            signal(6, caught_abort);
        }
        if std::env::var("VND_ORACLE").unwrap() == "c" {
            let msg = std::ffi::CString::new(message).unwrap();
            unsafe {
                c::VNDerr(msg.as_ptr());
            }
            panic!("C VNDerr returned");
        } else {
            VNDerr(&message);
        }
    }
    let mut rng = XorShift::new(0x512a);
    for i in 0..2048 {
        let len = [0, 1, 2, 7, 31, 64, 257, 1024][i % 8];
        let message: String = (0..len)
            .map(|_| match rng.next_u32() % 70 {
                0 => '\n',
                1 => '\t',
                k => (b' ' + k as u8) as char,
            })
            .collect();
        let run = |oracle: &str| {
            std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "vnderr", "--nocapture"])
                .env("VND_ERROR_MESSAGE", &message)
                .env("VND_ORACLE", oracle)
                .output()
                .unwrap()
        };
        let c = run("c");
        let r = run("r");
        assert_eq!(c.status.code(), Some(134));
        assert_eq!(r.status, c.status);
        assert_eq!(r.stderr, c.stderr, "error case {i}");
        assert_eq!(
            r.stderr,
            format!("VNDerr: fatal error\n{message}\n").as_bytes()
        );
    }
}

#[test]
fn vndlook() {
    if let Ok(oracle) = std::env::var("VND_LOOK_ORACLE") {
        let mut rng = XorShift::new(0x741ad);
        for i in 0..4000 {
            let len = [0, 1, 2, 7, 16, 31, 257, 1024][i % 8];
            let mut words: Vec<i32> = (0..len).map(|_| rng.next_u32() as i32).collect();
            if len > 0 {
                words[0] = [0, -1, 1, i32::MIN, i32::MAX][i % 5];
            }
            if oracle == "c" {
                unsafe {
                    c::VNDlook(words.as_ptr(), words.len() as i32);
                }
            } else {
                VNDlook(&words);
            }
        }
        unsafe {
            c::fflush(std::ptr::null_mut());
        }
        std::io::Write::flush(&mut std::io::stdout()).unwrap();
        std::process::exit(0);
    }
    let run = |oracle: &str| {
        std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "vndlook", "--nocapture"])
            .env("VND_LOOK_ORACLE", oracle)
            .output()
            .unwrap()
    };
    let c = run("c");
    let r = run("r");
    assert!(c.status.success());
    assert_eq!(r.status, c.status);
    assert_eq!(r.stdout, c.stdout);
    assert_eq!(r.stderr, c.stderr);
}

#[test]
fn memory_faults() {
    if let Ok(kind) = std::env::var("VND_MEMORY_FAULT") {
        unsafe {
            signal(6, caught_abort);
        }
        let c = std::env::var("VND_ORACLE").unwrap() == "c";
        let n: usize = std::env::var("VND_FAULT_SIZE").unwrap().parse().unwrap();
        let msg = std::env::var("VND_FAULT_MESSAGE").unwrap();
        let cm = std::ffi::CString::new(msg.clone()).unwrap();
        let p = if c {
            unsafe { c::VNDemalloc(n, cm.as_ptr()) }
        } else {
            VNDemalloc(n, &msg)
        };
        unsafe {
            let base = p.cast::<i64>().sub(2);
            let words = base.add(1).read() as usize;
            if kind == "begin" || kind == "free" {
                base.write(-123);
            } else {
                base.add(words + 2).write(0);
            }
            if kind == "free" {
                if c {
                    c::VNDfree(p, cm.as_ptr());
                } else {
                    VNDfree(p, &msg);
                }
            } else if c {
                c::VNDmemchk(p, cm.as_ptr());
            } else {
                VNDmemchk(p, &msg);
            }
        }
        panic!("damaged guard was accepted");
    }
    let mut rng = XorShift::new(0x144a);
    for i in 0..64 {
        let n = if i == 0 {
            0
        } else {
            1 + rng.next_u32() as usize % 8192
        };
        let msg = format!("guard test {:08x}\nline two", rng.next_u32());
        let kind = ["begin", "end", "free"][i % 3];
        let run = |oracle: &str| {
            std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "memory_faults", "--nocapture"])
                .env("VND_MEMORY_FAULT", kind)
                .env("VND_ORACLE", oracle)
                .env("VND_FAULT_SIZE", n.to_string())
                .env("VND_FAULT_MESSAGE", &msg)
                .output()
                .unwrap()
        };
        let c = run("c");
        let r = run("r");
        assert_eq!(c.status.code(), Some(134));
        assert_eq!(r.status, c.status);
        assert_eq!(r.stderr, c.stderr);
    }
}

#[test]
fn io_errors() {
    if let Ok(kind) = std::env::var("VND_IO_ERROR") {
        unsafe {
            signal(6, caught_abort);
        }
        let c = std::env::var("VND_ORACLE").unwrap() == "c";
        let root = workspace();
        let name = root.join("io-error");
        let name = name.to_str().unwrap();
        let cn = std::ffi::CString::new(name).unwrap();
        let mut r = None;
        let mut cv = std::ptr::null_mut();
        if c {
            cv = unsafe {
                c::VNDop(
                    2,
                    8192,
                    2,
                    [8, 7].as_ptr(),
                    2,
                    4,
                    cn.as_ptr(),
                    0,
                    std::ptr::null(),
                    2,
                )
            };
        } else {
            r = VNDop(2, 8192, &[8, 7], 2, 4, name, 0, &[], 2);
        }
        let mut iop = b'r';
        let mut iopanel = 0;
        let idim = 0;
        let mut key = [0, -1];
        let mut ipanel = 0;
        let mut start = 0;
        let mut count = 8;
        match kind.as_str() {
            "panel-low" => ipanel = -1,
            "panel-high" => ipanel = 2,
            "buffer-low" => iopanel = -1,
            "buffer-high" => iopanel = 2,
            "key-low" => key[1] = -1,
            "key-high" => key[1] = 8,
            "read-index" => {
                key[1] = 0;
                start = -1;
                count = 1;
            }
            "write-index" => {
                key[1] = 0;
                start = 8;
                count = 1;
                iop = b'w';
            }
            "invalid-mode" => {
                key[1] = 0;
                iop = b'x';
            }
            "slice-eof" => {
                let offset = if c {
                    unsafe { *(*cv).first_data.cast::<i64>() }
                } else {
                    r.as_ref().unwrap().first_data_byte[0]
                };
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(format!("{name}.VND0"))
                    .unwrap()
                    .set_len(offset as u64)
                    .unwrap();
                if c {
                    unsafe {
                        c::VNDrwslice(cv, 0, b'r' as i8, 0, 0, 0);
                    }
                } else {
                    VNDrwslice(r.as_mut().unwrap(), 0, b'r', 0, 0, 0);
                }
                panic!("EOF accepted");
            }
            _ => panic!("bad kind"),
        }
        let msg = std::env::var("VND_FAULT_MESSAGE").unwrap();
        let cm = std::ffi::CString::new(msg.clone()).unwrap();
        let mut values = vec![0; 32];
        if c {
            unsafe {
                c::VNDrw(
                    iop as i8,
                    iopanel,
                    cv,
                    idim,
                    key.as_ptr(),
                    ipanel,
                    values.as_mut_ptr(),
                    start,
                    1,
                    count,
                    -7,
                    cm.as_ptr(),
                );
            }
        } else {
            VNDrw(
                iop,
                iopanel,
                r.as_mut().unwrap(),
                idim,
                &key,
                ipanel,
                &mut values,
                start,
                1,
                count,
                -7,
                &msg,
            );
        }
        panic!("invalid access accepted");
    }
    let kinds = [
        "panel-low",
        "panel-high",
        "buffer-low",
        "buffer-high",
        "key-low",
        "key-high",
        "read-index",
        "write-index",
        "invalid-mode",
        "slice-eof",
    ];
    let mut rng = XorShift::new(0x8aa2);
    for i in 0..40 {
        let message = format!("bad input {:08x}", rng.next_u32());
        let run = |oracle: &str| {
            std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "io_errors", "--nocapture"])
                .env("VND_IO_ERROR", kinds[i % kinds.len()])
                .env("VND_ORACLE", oracle)
                .env("VND_FAULT_MESSAGE", &message)
                .output()
                .unwrap()
        };
        let c = run("c");
        let r = run("r");
        assert_eq!(
            c.status.code(),
            Some(134),
            "{}: {:?}",
            kinds[i % kinds.len()],
            c
        );
        assert_eq!(
            r.status,
            c.status,
            "{}: {}",
            kinds[i % kinds.len()],
            String::from_utf8_lossy(&r.stderr)
        );
        assert_eq!(r.stderr, c.stderr, "{}", kinds[i % kinds.len()]);
    }
    std::fs::remove_file(workspace().join("io-error.VND0")).unwrap();
}

#[test]
fn multiple_files() {
    if std::env::var_os("VND_MULTIFILE_CHILD").is_none() {
        let child = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "multiple_files", "--nocapture"])
            .env("VND_MULTIFILE_CHILD", "1")
            .output()
            .unwrap();
        assert!(
            child.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&child.stdout),
            String::from_utf8_lossy(&child.stderr)
        );
        return;
    }
    let _guard = MEMORY_TEST.lock().unwrap();
    let root = workspace();
    let mut rng = XorShift::new(0x74109);
    let dirs: Vec<_> = (0..3).map(|i| root.join(format!("dir{i}"))).collect();
    for d in &dirs {
        std::fs::create_dir_all(d).unwrap();
    }
    let ds: Vec<&str> = dirs.iter().map(|p| p.to_str().unwrap()).collect();
    let cds: Vec<_> = ds
        .iter()
        .map(|p| std::ffi::CString::new(*p).unwrap())
        .collect();
    let pointers: Vec<_> = cds.iter().map(|p| p.as_ptr()).collect();
    for i in 0..2048 {
        let n = [
            1 + (rng.next_u32() % 65) as i64,
            1 + (rng.next_u32() % 47) as i64,
        ];
        let bytes = [1, 4, 8][i % 3];
        let nf = 1 + i % 3;
        let max = [0, 4096, 65536][i % 3];
        let panels = 3;
        let name = std::ffi::CString::new("multi-c").unwrap();
        let c = unsafe {
            c::VNDop(
                2,
                max,
                2,
                n.as_ptr(),
                panels,
                bytes,
                name.as_ptr(),
                nf as i32,
                pointers.as_ptr(),
                1,
            )
        };
        assert!(!c.is_null());
        let mut r = VNDop(
            2,
            max,
            &n,
            panels,
            bytes,
            "multi-r",
            nf as i32,
            &ds[..nf],
            1,
        )
        .unwrap();
        unsafe {
            assert_state(&r, c);
        }
        for j in 0..9 {
            let d = j % 2;
            let key = [rng.next_u32() as i64 % n[0], rng.next_u32() as i64 % n[1]];
            let panel = (j % 3) as i64;
            let mut data = random_bytes(&mut rng, (n[d] * bytes) as usize, j == 0);
            let mut cd = data.clone();
            VNDrw(
                b'w', 0, &mut r, d as i32, &key, panel, &mut data, 0, 1, n[d], j as i32, "multi",
            );
            unsafe {
                c::VNDrw(
                    b'w' as i8,
                    0,
                    c,
                    d as i32,
                    key.as_ptr(),
                    panel,
                    cd.as_mut_ptr(),
                    0,
                    1,
                    n[d],
                    j as i32,
                    b"multi\0".as_ptr().cast(),
                );
            }
            data.fill(0);
            cd.fill(0);
            VNDrw(
                b'r', 0, &mut r, d as i32, &key, panel, &mut data, 0, 1, n[d], j as i32, "multi",
            );
            unsafe {
                c::VNDrw(
                    b'r' as i8,
                    0,
                    c,
                    d as i32,
                    key.as_ptr(),
                    panel,
                    cd.as_mut_ptr(),
                    0,
                    1,
                    n[d],
                    j as i32,
                    b"multi\0".as_ptr().cast(),
                );
            }
            assert_eq!(data, cd);
        }
        VNDflush(&mut r);
        unsafe {
            c::VNDflush(c);
            c::fflush((*c).fd);
            assert_state(&r, c);
        }
        for j in 0..nf {
            let rn = dirs[j].join(format!("multi-r.VND{j}"));
            let cn = dirs[j].join(format!("multi-c.VND{j}"));
            if rn.exists() {
                let rd = std::fs::read(rn).unwrap();
                let cd = std::fs::read(cn).unwrap();
                let offset = r.first_data_byte[j] as usize;
                assert_eq!(&rd[offset..], &cd[offset..]);
            }
        }
        VNDcl(r, 1);
        unsafe {
            c::VNDcl(c, 1);
        }
    }
}

#[test]
fn geometry_only_and_restart() {
    let _guard = MEMORY_TEST.lock().unwrap();
    let mut rng = XorShift::new(0x74b9);
    let root = workspace();
    for i in 0..2048 {
        let n = [
            2 * (1 + (rng.next_u32() % 32) as i64),
            1 + (rng.next_u32() % 31) as i64,
        ];
        let bytes = [1, 4, 8][i % 3];
        let max = [0, 1024, 8192, 65536][i % 4];
        // VNDop(3) reads an uninitialized NumFiles in C, so use the original
        // geometry function as the reference for this explicitly documented UB.
        let r = VNDop(3, max, &n, 1, bytes, "unused", 0, &[], 1).unwrap();
        let mut cn = n;
        let mut nodes = [0; 2];
        let mut blocks = [0; 2];
        let mut bi = [0; 2];
        let mut bli = [0; 2];
        let mut cv: CVnd = unsafe { std::mem::zeroed() };
        cv.num_dim = 2;
        cv.bytes_node = bytes;
        cv.n = cn.as_mut_ptr();
        cv.nodes = nodes.as_mut_ptr();
        cv.blocks = blocks.as_mut_ptr();
        cv.byte_incr = bi.as_mut_ptr();
        cv.block_incr = bli.as_mut_ptr();
        assert_eq!(unsafe { c::VNDGetDimensions(&mut cv, max, 1) }, 0);
        assert_eq!(r.nodes_per_block, nodes);
        assert_eq!(r.blocks, blocks);
        assert_eq!(r.num_bytes_mem_buf, cv.bytes_mem);
        VNDcl(r, 1);
        // Read files produced by the original C. C's LP64 reopen path cannot
        // be an oracle: it fread(sizeof(long)) into LenFileNames (int*).
        let name = root.join("restart");
        let name = name.to_str().unwrap();
        let cname = std::ffi::CString::new(name).unwrap();
        let c = unsafe { c::V2Dop(2, max, bytes, cname.as_ptr(), n[0], n[1]) };
        assert!(!c.is_null());
        let row = rng.next_u32() as i64 % n[1];
        let mut data = random_bytes(&mut rng, (n[0] * bytes) as usize, i % 11 == 0);
        unsafe {
            c::V2Dw0(c, row, data.as_mut_ptr(), -1);
            c::VNDcl(c, 0);
        }
        let mut r = VNDop(0, max, &n, 1, bytes, name, 0, &[], 1).unwrap();
        let mut got = vec![0; data.len()];
        V2Dr0(&mut r, row, &mut got, -1);
        assert_eq!(got, data);
        VNDcl(r, 0);
        let mut r = VNDop(1, max, &n, 1, bytes, name, 0, &[], 1).unwrap();
        let mut new = random_bytes(&mut rng, (n[1] * bytes) as usize, false);
        V2Dw1(&mut r, 0, &mut new, 0);
        VNDflush(&mut r);
        let mut check = vec![0; new.len()];
        V2Dr1(&mut r, 0, &mut check, 0);
        assert_eq!(check, new);
        VNDcl(r, 1);
    }
}

#[test]
fn directory_defaults() {
    if std::env::var_os("VND_DEFAULT_CHILD").is_none() {
        let root = workspace().join("defaults");
        std::fs::create_dir_all(&root).unwrap();
        let root = root.canonicalize().unwrap();
        let child = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "directory_defaults", "--nocapture"])
            .env("VND_DEFAULT_CHILD", "1")
            .current_dir(root)
            .output()
            .unwrap();
        assert!(
            child.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&child.stdout),
            String::from_utf8_lossy(&child.stderr)
        );
        return;
    }
    let mut rng = XorShift::new(0xa519);
    std::fs::create_dir_all("disk-a").unwrap();
    std::fs::create_dir_all("disk-b").unwrap();
    for i in 0..2048 {
        if i % 2 == 0 {
            std::fs::write(".VND", "2\ndisk-a\ndisk-b\n").unwrap();
        } else {
            std::fs::remove_file(".VND").unwrap();
        }
        let n = [
            1 + rng.next_u32() as i64 % 31,
            1 + rng.next_u32() as i64 % 37,
        ];
        let c = unsafe {
            c::VNDop(
                2,
                4096,
                2,
                n.as_ptr(),
                2,
                4,
                b"default-c\0".as_ptr().cast(),
                -1,
                std::ptr::null(),
                1,
            )
        };
        let mut r = VNDop(2, 4096, &n, 2, 4, "default-r", -1, &[], 1).unwrap();
        assert!(!c.is_null());
        unsafe {
            assert_state(&r, c);
        }
        let mut data = random_bytes(&mut rng, (n[0] * 4) as usize, i % 3 == 0);
        let mut cd = data.clone();
        V2Dw0(&mut r, 0, &mut data, 0);
        unsafe {
            c::V2Dw0(c, 0, cd.as_mut_ptr(), 0);
        }
        VNDflush(&mut r);
        unsafe {
            c::VNDflush(c);
        }
        data.fill(0);
        cd.fill(0);
        V2Dr0(&mut r, 0, &mut data, 0);
        unsafe {
            c::V2Dr0(c, 0, cd.as_mut_ptr(), 0);
        }
        assert_eq!(data, cd);
        VNDcl(r, 1);
        unsafe {
            c::VNDcl(c, 1);
        }
    }
}
