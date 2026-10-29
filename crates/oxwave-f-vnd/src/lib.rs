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

//! Port of par/lib/VND.c. Credits: J.E. Anderson (CSM Mobil Visiting Scientist, 1993).
#![allow(non_snake_case)]

/// Geometry used by the VND block partitioner.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VNDGeometry {
    pub n: Vec<i64>,
    pub num_bytes_per_node: i64,
    pub num_blocks_per_panel: i64,
    pub num_bytes_mem_buf: i64,
    pub num_bytes_per_block: i64,
    pub num_bytes_1st_dim: i64,
    pub num_nodes_per_block: i64,
    pub nodes_per_block: Vec<i64>,
    pub blocks: Vec<i64>,
    pub byte_incr: Vec<i64>,
    pub block_incr: Vec<i64>,
}
impl VNDGeometry {
    /// Construct the input geometry for VNDGetDimensions.
    pub fn new(n: &[i64], num_bytes_per_node: i64) -> Self {
        let z = vec![0; n.len()];
        Self {
            n: n.to_vec(),
            num_bytes_per_node,
            num_blocks_per_panel: 0,
            num_bytes_mem_buf: 0,
            num_bytes_per_block: 0,
            num_bytes_1st_dim: 0,
            num_nodes_per_block: 0,
            nodes_per_block: z.clone(),
            blocks: z.clone(),
            byte_incr: z.clone(),
            block_incr: z,
        }
    }
}

/// C: `int VNDGetDimensions(VND *vnd, long lwmax, int ich)`.
pub fn VNDGetDimensions(vnd: &mut VNDGeometry, lwmax: i64, ich: i32) -> i32 {
    #[cfg(target_os = "macos")]
    const SECTOR: i64 = 1024;
    #[cfg(not(target_os = "macos"))]
    const SECTOR: i64 = 8192;
    let mut total_nodes = 1i64;
    let mut max_n = 0;
    for &n in &vnd.n {
        total_nodes = total_nodes.wrapping_mul(n);
        if n > max_n {
            max_n = n;
        }
    }
    let mut total_bytes = total_nodes.wrapping_mul(vnd.num_bytes_per_node);
    if vnd.n.len() > 1 {
        total_bytes = SECTOR.wrapping_mul(1 + (total_bytes - 1) / SECTOR);
    }
    if total_bytes <= lwmax || vnd.n.len() == 1 {
        vnd.num_blocks_per_panel = 1;
        vnd.num_bytes_mem_buf = total_bytes;
        vnd.num_bytes_per_block = total_bytes;
        vnd.num_bytes_1st_dim = total_bytes;
        vnd.num_nodes_per_block = total_nodes;
        let mut byte = vnd.num_bytes_per_node;
        let mut block = 1i64;
        for j in 0..vnd.n.len() {
            vnd.nodes_per_block[j] = vnd.n[j];
            vnd.blocks[j] = 1;
            vnd.byte_incr[j] = byte;
            vnd.block_incr[j] = block;
            byte = byte.wrapping_mul(vnd.nodes_per_block[j]);
            block = block.wrapping_mul(vnd.blocks[j]);
        }
        return 0;
    }
    for k in 1..=max_n {
        let mut byte = vnd.num_bytes_per_node;
        let mut block = 1i64;
        let mut nodes = 1i64;
        for j in 0..vnd.n.len() {
            vnd.blocks[j] = if k > vnd.n[j] { vnd.n[j] } else { k };
            // arm64 C integer division by zero returns zero.
            vnd.nodes_per_block[j] = 1 + vnd.n[j]
                .wrapping_sub(1)
                .checked_div(vnd.blocks[j])
                .unwrap_or(0);
            if j == 0 && ich == 1 {
                let mut i = vnd.nodes_per_block[j];
                if i != 2 * (i / 2) {
                    i += 1;
                }
                vnd.nodes_per_block[j] = i;
                vnd.blocks[j] = 1 + (vnd.n[j] - 1) / i;
            }
            vnd.byte_incr[j] = byte;
            vnd.block_incr[j] = block;
            byte = byte.wrapping_mul(vnd.nodes_per_block[j]);
            block = block.wrapping_mul(vnd.blocks[j]);
            nodes = nodes.wrapping_mul(vnd.nodes_per_block[j]);
        }
        vnd.num_blocks_per_panel = block;
        byte = SECTOR.wrapping_mul(1 + (byte - 1) / SECTOR);
        vnd.num_bytes_per_block = byte;
        vnd.num_bytes_mem_buf = k.wrapping_mul(byte);
        vnd.num_bytes_1st_dim = vnd.blocks[0].wrapping_mul(byte);
        vnd.num_nodes_per_block = nodes;
        if byte == SECTOR {
            return 0;
        }
        if k.wrapping_mul(byte) <= lwmax && byte <= SECTOR * 8 {
            return 0;
        }
    }
    -1
}

/// C: `void VNDr2c(VND *vnd)`.
pub fn VNDr2c(vnd: &mut VNDGeometry) {
    vnd.num_bytes_per_node = vnd.num_bytes_per_node.wrapping_mul(2);
    vnd.n[0] /= 2;
    vnd.nodes_per_block[0] /= 2;
    vnd.byte_incr[0] = vnd.byte_incr[0].wrapping_mul(2);
}
/// C: `void VNDc2r(VND *vnd)`.
pub fn VNDc2r(vnd: &mut VNDGeometry) {
    vnd.num_bytes_per_node /= 2;
    vnd.n[0] = vnd.n[0].wrapping_mul(2);
    vnd.nodes_per_block[0] = vnd.nodes_per_block[0].wrapping_mul(2);
    vnd.byte_incr[0] /= 2;
}
/// C: `void vndr2c_(VND *vnd)`.
pub fn vndr2c_(vnd: &mut VNDGeometry) {
    VNDr2c(vnd);
}
/// C: `void vndc2r_(VND *vnd)`.
pub fn vndc2r_(vnd: &mut VNDGeometry) {
    VNDc2r(vnd);
}

use std::alloc::{alloc, dealloc, Layout};
use std::sync::atomic::{AtomicI32, AtomicI64, Ordering};

static TOTAL_MEM: AtomicI64 = AtomicI64::new(0);
static TEMP_COUNT: AtomicI32 = AtomicI32::new(0);

/// C: `void *VNDemalloc(size_t n, char *s)`.
/// The allocation must be released with `VNDfree`.
pub fn VNDemalloc(n: usize, message: &str) -> *mut u8 {
    let words = 1usize.wrapping_add(n.wrapping_sub(1) / std::mem::size_of::<i64>());
    // size_t arithmetic wraps in C. In particular n=0 allocates 24 bytes:
    // the huge logical word count wraps back to a guard at byte offset 16.
    let size = words.wrapping_add(3).wrapping_mul(8);
    let layout = Layout::from_size_align(size, 8).ok();
    let Some(layout) = layout else {
        eprintln!("Error allocating memory in VNDemalloc\n{message}");
        std::process::abort();
    };
    let p = unsafe { alloc(layout) as *mut i64 };
    if p.is_null() {
        eprintln!("Error allocating memory in VNDemalloc\n{message}");
        std::process::abort();
    }
    unsafe {
        p.write(6789);
        p.add(1).write((words as i32) as i64);
        (p as *mut u8)
            .wrapping_add(words.wrapping_add(2).wrapping_mul(8))
            .cast::<i64>()
            .write(1234);
    }
    TOTAL_MEM.fetch_add(layout.size() as i64, Ordering::SeqCst);
    unsafe { p.add(2) as *mut u8 }
}

/// C: `void VNDmemchk(void *p, char *msg)`.
///
/// # Safety
/// `p` must be a live pointer returned by `VNDemalloc`.
pub unsafe fn VNDmemchk(p: *mut u8, msg: &str) {
    let j = (p as *mut i64).sub(2);
    if j.read() != 6789 {
        eprintln!("Beginning memory overrun detected by VNDmemchk\n{msg}");
        std::process::abort();
    }
    if j.add(j.add(1).read() as usize + 2).read() != 1234 {
        eprintln!("Ending memory overrun detected by VNDmemchk\n{msg}");
        std::process::abort();
    }
}

/// C: `long VNDtotalmem(void)`.
pub fn VNDtotalmem() -> i64 {
    TOTAL_MEM.load(Ordering::SeqCst)
}

/// C: `void VNDfree(void *p, char *msg)`.
///
/// # Safety
/// `p` must be a live pointer returned by `VNDemalloc` and freed only once.
pub unsafe fn VNDfree(p: *mut u8, msg: &str) {
    VNDmemchk(p, msg);
    let j = (p as *mut i64).sub(2);
    let size = (j.add(1).read() as usize + 3) * 8;
    TOTAL_MEM.fetch_sub((size as i32) as i64, Ordering::SeqCst);
    dealloc(j as *mut u8, Layout::from_size_align_unchecked(size, 8));
}

/// C: `void VNDerr(char *s)`.
pub fn VNDerr(s: &str) -> ! {
    eprintln!("VNDerr: fatal error\n{s}");
    std::process::abort();
}

/// C: `void VNDlook(int *w, int lw)`.
pub fn VNDlook(w: &[i32]) {
    for (j, value) in w.iter().enumerate() {
        println!(" j={j} w[j]={value} ");
    }
}

/// C: `char *VNDtempname(char *root)`.
/// The returned pointer must be released with `VNDfree`.
pub fn VNDtempname(root: &str) -> *mut u8 {
    let count = TEMP_COUNT.fetch_add(1, Ordering::SeqCst) + 1;
    let name = format!("{}{}{}", root, std::process::id(), count);
    assert!(name.len() < 80, "VNDtempname buffer overflow");
    let p = VNDemalloc(80, "VNDtempname");
    unsafe {
        std::ptr::copy_nonoverlapping(name.as_ptr(), p, name.len());
        p.add(name.len()).write(0);
    }
    p
}

mod engine;
pub use engine::*;
