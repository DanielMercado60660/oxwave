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

//! File and slice operations from par/lib/VND.c.
//! Credits: J.E. Anderson (CSM Mobil Visiting Scientist, 1993).
use super::*;
use std::fs::{File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::ops::{Deref, DerefMut};
#[cfg(target_os = "macos")]
const SECTOR: i64 = 1024;
#[cfg(not(target_os = "macos"))]
const SECTOR: i64 = 8192;
// LP64 sizeof(VND): 11 longs, int and padding, 16 pointers.
const CONTROL_SIZE: usize = 224;

/// VND-managed byte storage, accessible as a slice. Its allocation has the
/// same leading and trailing guards as the original C buffer.
#[derive(Debug)]
pub struct VNDBuffer {
    p: *mut u8,
    len: usize,
}
impl Deref for VNDBuffer {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.p, self.len) }
    }
}
impl DerefMut for VNDBuffer {
    fn deref_mut(&mut self) -> &mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(self.p, self.len) }
    }
}
impl Drop for VNDBuffer {
    fn drop(&mut self) {
        if self.len != 0 {
            unsafe {
                VNDfree(self.p, "VNDcl: freeing vnd->w");
            }
        }
    }
}

/// Owned VND control structure. Dimension zero varies fastest.
/// Dereferences to geometry; w contains the open panel buffers.
#[derive(Debug)]
pub struct VND {
    pub geometry: VNDGeometry,
    pub num_files: i64,
    pub num_panels: i64,
    pub num_open_panels: i64,
    pub num_blocks_per_file: i64,
    pub current_file: i32,
    pub mode: String,
    pub file_names: Vec<String>,
    pub first_data_byte: Vec<i64>,
    pub len_file: Vec<i64>,
    pub current_panel: Vec<i64>,
    pub current_dimension: Vec<i32>,
    pub current_block: Vec<i64>,
    pub modified: Vec<i32>,
    pub w: VNDBuffer,
    file: Option<File>,
    // Guard allocations preserve VNDtotalmem accounting. Typed vectors own
    // Rust data independently, so process-local C pointers never escape.
    allocations: Vec<*mut u8>,
}
impl Deref for VND {
    type Target = VNDGeometry;
    fn deref(&self) -> &VNDGeometry {
        &self.geometry
    }
}
impl DerefMut for VND {
    fn deref_mut(&mut self) -> &mut VNDGeometry {
        &mut self.geometry
    }
}
impl Drop for VND {
    fn drop(&mut self) {
        for p in self.allocations.drain(..) {
            unsafe {
                VNDfree(p, "VNDcl: freeing vnd");
            }
        }
    }
}
fn open_file(name: &str, mode: &str) -> std::io::Result<File> {
    let mut op = OpenOptions::new();
    match mode {
        "w" | "w+" => {
            op.write(true).create(true).truncate(true);
            if mode == "w+" {
                op.read(true);
            }
        }
        "r+" => {
            op.read(true).write(true);
        }
        _ => {
            op.read(true);
        }
    }
    op.open(name)
}
fn put_long(w: &mut impl Write, n: i64) {
    w.write_all(&n.to_ne_bytes())
        .unwrap_or_else(|_| VNDerr("VNDopen: write error"));
}
fn get_long(r: &mut impl Read) -> i64 {
    let mut b = [0; 8];
    r.read_exact(&mut b)
        .unwrap_or_else(|_| VNDerr("VNDopen: read error"));
    i64::from_ne_bytes(b)
}
fn charge(v: &mut VND, n: usize) {
    v.allocations
        .push(VNDemalloc(n, "VNDopen: allocating memory"));
}

/// C: `VND *VNDop(int mode, long lwmax, int ndim, long *N, long npanels,
/// long nbytes, char *file, int ndir, char **dir, int noppanels)`.
/// N.len() replaces ndim. None means that the first file cannot be opened.
/// Mode 3 computes geometry without files (C leaves its file count undefined).
#[allow(clippy::too_many_arguments)]
pub fn VNDop(
    mode: i32,
    lwmax: i64,
    n: &[i64],
    npanels: i64,
    nbytes: i64,
    file: &str,
    ndir: i32,
    dir: &[&str],
    noppanels: i32,
) -> Option<VND> {
    assert!(!n.is_empty() && noppanels > 0);
    let control = VNDemalloc(CONTROL_SIZE, "VNDop: allocating VND");
    let mut dirs: Vec<String> = dir.iter().map(|s| s.to_string()).collect();
    if ndir < 0 {
        dirs.clear();
        if let Ok(text) = std::fs::read_to_string(".VND") {
            let mut words = text.split_whitespace();
            let count: usize = words.next().unwrap_or("0").parse().unwrap_or(0);
            dirs.extend(words.take(count).map(str::to_owned));
        }
    } else if ndir == 0 {
        dirs.clear();
    } else {
        dirs.truncate(ndir as usize);
    }
    let num_files = if dirs.is_empty() { 1 } else { dirs.len() };
    let names: Vec<String> = (0..num_files)
        .map(|j| {
            if dirs.is_empty() {
                format!("{file}.VND{j}")
            } else {
                format!("{}/{file}.VND{j}", dirs[j])
            }
        })
        .collect();
    let mut fp = if mode == 3 {
        None
    } else {
        match open_file(&names[0], if mode <= 1 { "r" } else { "w+" }) {
            Ok(f) => Some(f),
            Err(_) => {
                if mode <= 1 {
                    eprintln!("Could not open file {} ", names[0]);
                } else {
                    eprintln!("Could not create file {} ", names[0]);
                }
                return None;
            }
        }
    };
    let mut v = VND {
        geometry: VNDGeometry::new(n, nbytes),
        num_files: num_files as i64,
        num_panels: npanels,
        num_open_panels: noppanels as i64,
        num_blocks_per_file: 0,
        current_file: 0,
        mode: String::new(),
        file_names: names,
        first_data_byte: vec![0; num_files],
        len_file: vec![0; num_files],
        current_panel: vec![-1; noppanels as usize],
        current_dimension: vec![-1; noppanels as usize],
        current_block: vec![-1; noppanels as usize],
        modified: vec![0; noppanels as usize],
        w: VNDBuffer {
            p: std::ptr::NonNull::dangling().as_ptr(),
            len: 0,
        },
        file: None,
        allocations: vec![control],
    };
    if mode <= 1 {
        let f = fp.as_mut().unwrap();
        let mut header = [0; CONTROL_SIZE];
        f.read_exact(&mut header)
            .unwrap_or_else(|_| VNDerr(&v.file_names[0]));
        let vals: Vec<i64> = header[..88]
            .chunks_exact(8)
            .map(|b| i64::from_ne_bytes(b.try_into().unwrap()))
            .collect();
        v.geometry = VNDGeometry::new(&vec![0; vals[0] as usize], vals[1]);
        v.geometry.num_blocks_per_panel = vals[2];
        v.num_files = vals[3];
        v.num_panels = vals[4];
        v.geometry.num_bytes_mem_buf = vals[6];
        v.geometry.num_bytes_per_block = vals[7];
        v.geometry.num_bytes_1st_dim = vals[8];
        v.geometry.num_nodes_per_block = vals[9];
        v.num_blocks_per_file = vals[10];
        for a in [
            &mut v.geometry.n,
            &mut v.geometry.nodes_per_block,
            &mut v.geometry.blocks,
            &mut v.geometry.byte_incr,
            &mut v.geometry.block_incr,
        ] {
            for x in a {
                *x = get_long(f);
            }
        }
        let nf = v.num_files as usize;
        // C fread(sizeof(long)) into an int array here corrupts the following
        // reads. Decode the int elements actually written by VNDop instead.
        let mut lens = vec![0; nf];
        for x in &mut lens {
            let mut b = [0; 4];
            f.read_exact(&mut b)
                .unwrap_or_else(|_| VNDerr("VNDopen: read error"));
            *x = i32::from_ne_bytes(b);
        }
        v.first_data_byte = (0..nf).map(|_| get_long(f)).collect();
        v.len_file = (0..nf).map(|_| get_long(f)).collect();
        v.file_names.clear();
        for len in lens {
            let mut name = vec![0; len as usize];
            f.read_exact(&mut name)
                .unwrap_or_else(|_| VNDerr("VNDopen: read error"));
            name.pop();
            v.file_names.push(String::from_utf8(name).unwrap());
        }
    } else {
        VNDGetDimensions(&mut v.geometry, lwmax / noppanels as i64, 1);
    }
    let nf = v.num_files as usize;
    let nd = v.n.len();
    for len in [(nf + 1) * 4, (nf + 1) * 8, (nf + 1) * 8, (nf + 1) * 8] {
        charge(&mut v, len);
    }
    if mode != 3 {
        for j in 0..nf {
            let len = v.file_names[j].len() + 1;
            charge(&mut v, len);
        }
    }
    for width in [8, 4, 4, 8] {
        charge(&mut v, noppanels as usize * width);
    }
    for _ in 0..5 {
        charge(&mut v, nd * 8);
    }
    if mode == 3 {
        return Some(v);
    }
    let len = (v.num_bytes_mem_buf * noppanels as i64) as usize;
    v.w = VNDBuffer {
        p: VNDemalloc(len, "VNDopen: allocating w"),
        len,
    };
    // C only clears the first panel's buffer. Clear all allocated bytes here
    // so slices never expose uninitialized memory in Rust.
    v.w.fill(0);
    if mode == 2 {
        for name in &v.file_names {
            v.first_data_byte[0] += (name.len() + 1) as i64;
        }
        v.first_data_byte[0] += (CONTROL_SIZE + nd * 8 * 5 + nf * 8 * 2 + nf * 4) as i64;
        if nd > 1 {
            v.first_data_byte[0] = SECTOR * (1 + (v.first_data_byte[0] - 1) / SECTOR);
        }
        v.geometry.num_bytes_1st_dim = v.num_bytes_per_block * v.blocks[0];
        v.num_blocks_per_file = (v.num_blocks_per_panel * v.num_panels) / v.blocks[0];
        v.num_blocks_per_file = 1 + (v.num_blocks_per_file - 1) / v.num_files;
        v.num_blocks_per_file *= v.blocks[0];
        for j in 0..nf {
            v.len_file[j] = v.first_data_byte[j] + v.num_blocks_per_file * v.num_bytes_per_block;
        }
        let f = fp.as_mut().unwrap();
        for x in [
            nd as i64,
            v.num_bytes_per_node,
            v.num_blocks_per_panel,
            v.num_files,
            v.num_panels,
            v.num_open_panels,
            v.num_bytes_mem_buf,
            v.num_bytes_per_block,
            v.num_bytes_1st_dim,
            v.num_nodes_per_block,
            v.num_blocks_per_file,
        ] {
            put_long(f, x);
        }
        // C stores process-local pointers here. Zero these nonpersistent
        // slots; all persistent information follows in typed arrays.
        f.write_all(&vec![0; CONTROL_SIZE - 88])
            .unwrap_or_else(|_| VNDerr("VNDopen: write error"));
        for a in [
            &v.n,
            &v.nodes_per_block,
            &v.blocks,
            &v.byte_incr,
            &v.block_incr,
        ] {
            for &x in a {
                put_long(f, x);
            }
        }
        for name in &v.file_names {
            f.write_all(&((name.len() + 1) as i32).to_ne_bytes())
                .unwrap_or_else(|_| VNDerr("VNDopen: write error"));
        }
        for &x in &v.first_data_byte {
            put_long(f, x);
        }
        for &x in &v.len_file {
            put_long(f, x);
        }
        for name in &v.file_names {
            f.write_all(name.as_bytes())
                .unwrap_or_else(|_| VNDerr("VNDopen: write error"));
            f.write_all(&[0])
                .unwrap_or_else(|_| VNDerr("VNDopen: write error"));
        }
        let padding = v.first_data_byte[0] - f.stream_position().unwrap() as i64;
        if padding > 0 {
            f.write_all(&vec![0; padding as usize]).unwrap();
        }
        v.file = fp.take();
        v.mode = "w".into();
        for ipanel in 0..v.num_panels {
            let mut iblock = 0;
            while iblock < v.num_blocks_per_panel {
                VNDrwslice(&mut v, 0, b'w', 0, iblock, ipanel);
                iblock += v.blocks[0];
            }
        }
        v.file.take();
    }
    v.mode = if mode == 0 { "r" } else { "r+" }.into();
    v.file = open_file(&v.file_names[0], &v.mode).ok();
    Some(v)
}

/// C: `int VNDseek(VND *vnd, long ib)`.
pub fn VNDseek(v: &mut VND, mut ib: i64) -> i32 {
    let ifile = (ib / v.num_blocks_per_file) as i32;
    if v.current_file == -1 || ifile != v.current_file {
        let reopen = v.current_file != -1;
        v.file.take();
        match open_file(&v.file_names[ifile as usize], &v.mode) {
            Ok(f) => {
                v.file = Some(f);
                v.current_file = ifile;
            }
            Err(_) => {
                eprint!(
                    "{}",
                    if reopen {
                        "freopen failed in VNDseek"
                    } else {
                        "fopen failed in VNDseek"
                    }
                );
                return -1;
            }
        }
    }
    ib -= ifile as i64 * v.num_blocks_per_file;
    ib = ib * v.num_bytes_per_block + v.first_data_byte[ifile as usize];
    if ib < 0
        || v.file
            .as_mut()
            .unwrap()
            .seek(SeekFrom::Start(ib as u64))
            .is_err()
    {
        eprint!("fseeko failed in VNDseek");
        return -2;
    }
    0
}

// fread/fwrite report complete elements. VND slices use byte elements.
fn io_bytes(f: &mut File, bytes: &mut [u8], mode: u8) -> usize {
    let mut done = 0;
    while done < bytes.len() {
        let result = if mode == b'r' {
            f.read(&mut bytes[done..])
        } else {
            f.write(&bytes[done..])
        };
        match result {
            Ok(0) => break,
            Ok(n) => done += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(_) => break,
        }
    }
    done
}
fn slice_context(v: &VND, mode: u8, iopanel: i32, idim: i32, iblock: i64, ipanel: i64, dump: bool) {
    eprintln!(
        " mode = {} iopanel = {iopanel} idim = {idim} iblock = {iblock} ipanel = {ipanel} ",
        mode as char
    );
    if dump {
        VNDdump(v, &mut std::io::stderr());
    }
}

/// C: `int VNDrwslice(VND *vnd, int iopanel, char mode, int idim, long iblock, long ipanel)`.
pub fn VNDrwslice(v: &mut VND, iopanel: i32, mode: u8, idim: i32, iblock: i64, ipanel: i64) -> i32 {
    let base = iopanel as usize * v.num_bytes_mem_buf as usize;
    let mut ib = ipanel * v.num_blocks_per_panel + iblock;
    let contiguous = idim == 0 && v.num_bytes_mem_buf < SECTOR * 8;
    let n2read = if contiguous {
        v.num_bytes_1st_dim
    } else {
        v.num_bytes_per_block
    } as i32;
    let count = if contiguous {
        1
    } else {
        v.blocks[idim as usize]
    };
    for j in 0..count {
        if VNDseek(v, ib) != 0 {
            eprint!(
                "{}",
                if contiguous {
                    "VNDrwslice: seek failed on dim 0 \n"
                } else {
                    "VNDrwslice: seek failed"
                }
            );
            slice_context(v, mode, iopanel, idim, iblock, ipanel, contiguous);
            std::process::abort();
        }
        let offset = base + j as usize * n2read as usize;
        let bytes = &mut v.w[offset..offset + n2read as usize];
        let n = io_bytes(v.file.as_mut().unwrap(), bytes, mode) as i32;
        if n != n2read {
            if contiguous {
                eprintln!("VNDrwslice: wrong number of values");
                if mode == b'r' {
                    eprint!(" Read {n} bytes out of {n2read} on 1st dim");
                } else {
                    eprint!(" Wrote {n} bytes out of {n2read} on 1st dim");
                }
            } else {
                eprint!("VNDrwslice wrong number of values");
                if mode == b'r' {
                    eprintln!(" Read {n} bytes out of {n2read} on dim {idim} ");
                } else {
                    eprintln!(" Wrote {n} bytes out of {n2read} on dim {idim} ");
                }
            }
            slice_context(v, mode, iopanel, idim, iblock, ipanel, true);
            std::process::abort();
        }
        ib += v.block_incr[idim as usize];
    }
    0
}

/// C: `void VNDflush(VND *vnd)`.
pub fn VNDflush(v: &mut VND) {
    for j in 0..v.num_open_panels as usize {
        if v.modified[j] == 1 {
            VNDrwslice(
                v,
                j as i32,
                b'w',
                v.current_dimension[j],
                v.current_block[j],
                v.current_panel[j],
            );
        }
    }
}

/// C: `void VNDrw(char iop, int iopanel, VND *vnd, int idim, long *key,
/// long ipanel, char *values, long start, long incr, long nvalues, int msgnum, char *msg)`.
#[allow(clippy::too_many_arguments)]
pub fn VNDrw(
    iop: u8,
    iopanel: i32,
    v: &mut VND,
    idim: i32,
    key: &[i64],
    ipanel: i64,
    values: &mut [u8],
    start: i64,
    incr: i64,
    nvalues: i64,
    msgnum: i32,
    msg: &str,
) {
    let mut iblock = 0;
    let mut ibyte = 0;
    if ipanel < 0 || ipanel > v.num_panels - 1 {
        rw_error(
            "VNDrw: ipanel is out of range \n",
            iopanel,
            idim,
            key,
            ipanel,
            msgnum,
            msg,
        );
    }
    let d = idim as usize;
    if v.n.len() == 1 && start == 0 && incr == 1 && nvalues == v.n[0] {
        if VNDseek(v, ipanel) != 0 {
            eprintln!("VNDrw: Seek error on trace {ipanel} ");
            eprintln!("VNDrw: Message # {msgnum} {msg} ");
            std::process::abort();
        }
        let count = (v.num_bytes_per_node * v.n[0]) as usize;
        let nbytes = v.num_bytes_per_node as usize;
        let n = if iop == b'r' || iop == b'w' {
            io_bytes(v.file.as_mut().unwrap(), &mut values[..count], iop)
                .checked_div(nbytes)
                .unwrap_or(0) as i32
        } else {
            0
        };
        if n as i64 != v.n[0] {
            eprintln!("VNDrw: I/O error on trace {ipanel} ");
            eprintln!("       Only read/wrote {n} values out of {} ", v.n[0]);
            eprintln!("VNDrw: Message # {msgnum} {msg} ");
            std::process::abort();
        }
        return;
    }
    if iopanel < 0 || iopanel as i64 > v.num_open_panels - 1 {
        rw_error(
            "VNDWR: iopanel is out of range \n",
            iopanel,
            idim,
            key,
            ipanel,
            msgnum,
            msg,
        );
    }
    for j in 0..v.n.len() {
        if j != d {
            if key[j] < 0 || key[j] > v.n[j] {
                rw_error(
                    &format!(
                        "VNDrw: key of {} for dimension {j} is out of range \n",
                        key[j]
                    ),
                    iopanel,
                    idim,
                    key,
                    ipanel,
                    msgnum,
                    msg,
                );
            }
            let ib = key[j] / v.nodes_per_block[j];
            let ik = key[j] - ib * v.nodes_per_block[j];
            iblock += ib * v.block_incr[j];
            ibyte += ik * v.byte_incr[j];
        }
    }
    ibyte += v.num_bytes_mem_buf * iopanel as i64;
    let p = iopanel as usize;
    let mut incore = ipanel == v.current_panel[p];
    if v.num_blocks_per_panel > 1 {
        if idim != v.current_dimension[p] {
            incore = false;
        }
        if iblock != v.current_block[p] {
            incore = false;
        }
    }
    if !incore {
        if v.modified[p] == 1 {
            VNDrwslice(
                v,
                iopanel,
                b'w',
                v.current_dimension[p],
                v.current_block[p],
                v.current_panel[p],
            );
        }
        VNDrwslice(v, iopanel, b'r', idim, iblock, ipanel);
        v.current_dimension[p] = idim;
        v.current_block[p] = iblock;
        v.current_panel[p] = ipanel;
        v.modified[p] = 0;
    }
    if iop != b'r' && iop != b'w' {
        eprint!("VNDrw: Invalid option.  Not  r  or  w   \n");
        std::process::abort();
    }
    if iop == b'w' {
        v.modified[p] = 1;
    }
    if start == 0 && incr == 1 && nvalues == v.n[d] {
        let mut k = 0usize;
        let mut jj = 0;
        if d == 0 {
            let mut nbytes = v.num_bytes_per_node * v.nodes_per_block[0];
            let nbytesm = v.num_bytes_per_node * v.n[0];
            for _ in 0..v.blocks[0] {
                nbytes = nbytes.min(nbytesm - jj);
                for j in 0..nbytes {
                    copy_byte(v, iop, values, k, (ibyte + j) as usize);
                    k += 1;
                }
                ibyte += v.num_bytes_per_block;
                jj += v.nodes_per_block[0] * v.num_bytes_per_node;
            }
        } else {
            for _ in 0..v.blocks[d] {
                let nodes = (v.n[d] - jj).min(v.nodes_per_block[d]);
                let mut jjj = ibyte as i32;
                for _ in 0..nodes {
                    for kk in 0..v.num_bytes_per_node {
                        copy_byte(v, iop, values, k, (jjj as i64 + kk) as usize);
                        k += 1;
                    }
                    jjj = (jjj as i64 + v.byte_incr[d]) as i32;
                }
                ibyte += v.num_bytes_per_block;
                jj += v.nodes_per_block[d];
            }
        }
    } else {
        let mut iv = 0;
        for j in 0..nvalues {
            let k = (start + j * incr) as i32;
            if k < 0 || k as i64 >= v.n[d] {
                eprintln!(
                    "VNDrw: index out of range on {}",
                    if iop == b'r' { "read" } else { "write" }
                );
                eprintln!("VNDrw: Message # {msgnum} {msg} ");
                eprintln!("VNDrw: dimension (0 based for C, add 1 for Fortran) = {idim}");
                eprintln!("VNDrw: number of nodes = {}", v.n[d]);
                eprintln!("VNDrw: index (0 based for C, add 1 for Fortran) = {k}");
                eprintln!(
                    "VNDrw: erroneous values for start, incr, or nvalues might cause this error"
                );
                eprintln!("VNDrw: start (0 based for C, add 1 for Fortran) = {start}");
                eprintln!("VNDrw: incr (0 based for C, add 1 for Fortran) = {incr}");
                eprintln!("VNDrw: nvalues = {nvalues}");
                rw_tail(iopanel, idim, key, ipanel);
                std::process::abort();
            }
            let ib = k as i64 / v.nodes_per_block[d];
            let mut iw = k as i64 - ib * v.nodes_per_block[d];
            iw = ibyte + iw * v.byte_incr[d] + ib * v.num_bytes_per_block;
            for _ in 0..v.num_bytes_per_node {
                copy_byte(v, iop, values, iv, iw as usize);
                iw += 1;
                iv += 1;
            }
        }
    }
}
fn copy_byte(v: &mut VND, iop: u8, values: &mut [u8], iv: usize, iw: usize) {
    if iop == b'r' {
        values[iv] = v.w[iw];
    } else {
        v.w[iw] = values[iv];
    }
}
fn rw_tail(iopanel: i32, idim: i32, key: &[i64], ipanel: i64) {
    eprintln!("VNDrw: iopanel = {iopanel} idim = {idim} ipanel = {ipanel}");
    for (jj, k) in key.iter().enumerate() {
        eprintln!("      Dimension = {jj}  key = {k} ");
    }
}
fn rw_error(
    text: &str,
    iopanel: i32,
    idim: i32,
    key: &[i64],
    ipanel: i64,
    msgnum: i32,
    msg: &str,
) -> ! {
    eprint!("{text}");
    eprintln!("VNDrw: Message # {msgnum} {msg} ");
    rw_tail(iopanel, idim, key, ipanel);
    std::process::abort();
}

/// C: `void VNDcl(VND *vnd, int mode)`.
/// Consumes the handle; mode 1 also removes all physical files.
pub fn VNDcl(mut v: VND, mode: i32) {
    if v.file.is_some() {
        VNDflush(&mut v);
    }
    v.file.take();
    if mode != 1 && !v.w.is_empty() {
        // The original retains the filename allocations when files are kept.
        // Preserve that externally observable VNDtotalmem behavior.
        let nf = v.num_files as usize;
        v.allocations.drain(5..5 + nf);
    }
    if mode == 1 {
        for name in &v.file_names {
            let _ = std::fs::remove_file(name);
        }
    }
}

/// C: `void VNDdump(VND *vnd, FILE *fp)`.
pub fn VNDdump(v: &VND, fp: &mut impl Write) {
    writeln!(fp, "\nDump of VND structure").unwrap();
    for (name, x) in [
        ("NumDim", v.n.len() as i64),
        ("NumBytesPerNode", v.num_bytes_per_node),
        ("NumBlocksPerPanel", v.num_blocks_per_panel),
        ("NumFiles", v.num_files),
        ("NumPanels", v.num_panels),
        ("NumOpenPanels", v.num_open_panels),
        ("NumBytesMemBuf", v.num_bytes_mem_buf),
        ("NumBytesPerBlock", v.num_bytes_per_block),
        ("NumBytes1stDim", v.num_bytes_1st_dim),
        ("NumNodesPerBlock", v.num_nodes_per_block),
        ("NumBlocksPerFile", v.num_blocks_per_file),
    ] {
        writeln!(fp, "{name} = {x}").unwrap();
    }
    writeln!(
        fp,
        "CurrentFile = {}\nMode = {}\nFileDescriptor will not be dumped",
        v.current_file, v.mode
    )
    .unwrap();
    writeln!(fp, "\nFor each block matrix transpose dimension:").unwrap();
    for j in 0..v.n.len() {
        writeln!(
            fp,
            "Dimension = {j} N = {} NNodesPerBlock = {} NBlocks = {} ByteIncr = {} BlockIncr = {}",
            v.n[j], v.nodes_per_block[j], v.blocks[j], v.byte_incr[j], v.block_incr[j]
        )
        .unwrap();
    }
    writeln!(fp, "\nFor each simultaneously open panel buffer:").unwrap();
    for j in 0..v.num_open_panels as usize {
        writeln!(fp,"Open panel buffer = {j} Current Panel = {} CurrentDimension = {} CurrentBlock = {} Modified = {} ",
            v.current_panel[j],v.current_dimension[j],v.current_block[j],v.modified[j]).unwrap();
    }
    writeln!(fp, "\nFor each VND file:").unwrap();
    for j in 0..v.num_files as usize {
        writeln!(
            fp,
            "File name = {}\n   LenFileNames = {} FirstDataByte = {} LenFile = {} ",
            v.file_names[j],
            v.file_names[j].len() + 1,
            v.first_data_byte[j],
            v.len_file[j]
        )
        .unwrap();
    }
    write!(fp, "\n \n").unwrap();
}

/// C: `VND *V2Dop(int mode, long lwmax, long nbytes, char *file, long n0, long n1)`.
pub fn V2Dop(mode: i32, lwmax: i64, nbytes: i64, file: &str, n0: i64, n1: i64) -> Option<VND> {
    VNDop(mode, lwmax, &[n0, n1], 1, nbytes, file, 0, &[], 1)
}

/// C: `void V2Dr0(VND *vnd, long key, char *buf, int msgnum)`.
pub fn V2Dr0(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i32) {
    VNDrw(
        b'r',
        0,
        v,
        0,
        &[0, key],
        0,
        buf,
        0,
        1,
        v.n[0],
        msgnum as i32,
        "v2dr0",
    );
}

/// C: `void V2Dr1(VND *vnd, long key, char *buf, int msgnum)`.
pub fn V2Dr1(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i32) {
    VNDrw(
        b'r',
        0,
        v,
        1,
        &[key, 0],
        0,
        buf,
        0,
        1,
        v.n[1],
        msgnum as i32,
        "v2dr1",
    );
}

/// C: `void V2Dw0(VND *vnd, long key, char *buf, int msgnum)`.
pub fn V2Dw0(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i32) {
    VNDrw(
        b'w',
        0,
        v,
        0,
        &[0, key],
        0,
        buf,
        0,
        1,
        v.n[0],
        msgnum as i32,
        "v2dw0",
    );
}

/// C: `void V2Dw1(VND *vnd, long key, char *buf, int msgnum)`.
pub fn V2Dw1(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i32) {
    VNDrw(
        b'w',
        0,
        v,
        1,
        &[key, 0],
        0,
        buf,
        0,
        1,
        v.n[1],
        msgnum as i32,
        "v2dw1",
    );
}

/// C: `void V2Dr1_(VND *vnd, FORTINT key, void *buf, FORTINT msgnum)`.
pub fn V2Dr1_(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i64) {
    VNDrw(
        b'r',
        0,
        v,
        0,
        &[0, key - 1],
        0,
        buf,
        0,
        1,
        v.n[0],
        msgnum as i32,
        "V2DR1",
    );
}

/// C: `void v2dr2_(VND *vnd, FORTINT key, void *buf, FORTINT msgnum)`.
pub fn v2dr2_(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i64) {
    VNDrw(
        b'r',
        0,
        v,
        1,
        &[key - 1, 0],
        0,
        buf,
        0,
        1,
        v.n[1],
        msgnum as i32,
        "V2DR2",
    );
}

/// C: `void v2dw1_(VND *vnd, FORTINT key, void *buf, FORTINT msgnum)`.
pub fn v2dw1_(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i64) {
    VNDrw(
        b'w',
        0,
        v,
        0,
        &[0, key - 1],
        0,
        buf,
        0,
        1,
        v.n[0],
        msgnum as i32,
        "V2DW1",
    );
}

/// C: `void v2dw2_(VND *vnd, FORTINT key, void *buf, FORTINT msgnum)`.
pub fn v2dw2_(v: &mut VND, key: i64, buf: &mut [u8], msgnum: i64) {
    VNDrw(
        b'w',
        0,
        v,
        1,
        &[key - 1, 0],
        0,
        buf,
        0,
        1,
        v.n[1],
        msgnum as i32,
        "V2DW2",
    );
}

/// C: `void vndcl_(VND **vnd, FORTINT *fort_mode)`.
pub fn vndcl_(v: VND, mode: i64) {
    VNDcl(v, mode as i32);
}
/// C: `void vndflush_(VND **vnd)`.
pub fn vndflush_(v: &mut VND) {
    VNDflush(v);
}
