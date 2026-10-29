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

//! `segy`: what segyread/segywrite and friends need from SEG-Y: the binary
//! reel header `bhed` (segy.h) with the `bhdr` key table (bhdr.h), the
//! 3200-byte textual header conversions, IBM <-> IEEE float conversion
//! bit-exact with segyread.c/segywrite.c, and cwp/lib/swapbyte.c.
//!
//! | C                                   | Rust                                    |
//! |-------------------------------------|-----------------------------------------|
//! | `bhed bh;`                          | `let mut bh = Bhed::zeroed();`          |
//! | `tapehdr[]`, `tapebhdr[]`, `gettapehval/puttapehval/gettapebhval/puttapebhval` | `TAPEHDR`, `TAPEBHDR`, same names on raw byte buffers |
//! | `bhdr[i].key/type/offs`, BHED_NKEYS | `BHDR[i]`, `BHED_NKEYS`                 |
//! | `getbhval/putbhval/swapbhval`       | same names                              |
//! | `xdrbhdrsub` (encode/decode)        | `bh.to_xdr()` / `Bhed::from_xdr(&buf)`  |
//! | `EBCBYTES, BNYBYTES`                | same                                    |
//! | `dd conv=ascii` (segyread)          | `dd_conv_ascii(&mut buf)` / `DD_EBCDIC_TO_ASCII` |
//! | `dd conv=ebcdic`                    | `dd_conv_ebcdic(&mut buf)` / `DD_ASCII_TO_EBCDIC` |
//! | `zebc(in, out, n)` (segywrite)      | `zebc(&input, &mut output, n)`          |
//! | `ibm_to_float(from,to,n,endian,v)`  | `ibm_to_float(&mut data[..], n, endian)` (in place, f32 slots hold raw words) |
//! | `float_to_ibm(from,to,n,endian)`    | `float_to_ibm(&mut data[..], n, endian)` |
//! | `swap_short_2(&x)` etc.             | `swap_short_2(&mut x)` etc.             |
//!
//! Note: segyread's text header goes through `dd conv=ascii`, segywrite's
//! through `zebc`; the tables differ (zebc maps everything outside
//! 32..=127 to an ASCII blank 0x20).

use crate::value::Value;

/// Bytes in the card image EBCDIC block.
pub const EBCBYTES: usize = 3200;
/// Bytes in the binary coded block.
pub const BNYBYTES: usize = 400;
/// Number of mandated binary fields.
pub const BHED_NKEYS: usize = 27;

/// SEG-Y binary reel header, identical to `bhed` in segy.h (400 bytes).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bhed {
    pub jobid: i32,
    pub lino: i32,
    pub reno: i32,
    pub ntrpr: i16,
    pub nart: i16,
    pub hdt: u16,
    pub dto: u16,
    pub hns: u16,
    pub nso: u16,
    pub format: i16,
    pub fold: i16,
    pub tsort: i16,
    pub vscode: i16,
    pub hsfs: i16,
    pub hsfe: i16,
    pub hslen: i16,
    pub hstyp: i16,
    pub schn: i16,
    pub hstas: i16,
    pub hstae: i16,
    pub htatyp: i16,
    pub hcorr: i16,
    pub bgrcv: i16,
    pub rcvm: i16,
    pub mfeet: i16,
    pub polyt: i16,
    pub vpol: i16,
    pub hunass: [i16; 170],
}

const _: () = assert!(core::mem::size_of::<Bhed>() == BNYBYTES);
const _: () = assert!(core::mem::offset_of!(Bhed, ntrpr) == 12);
const _: () = assert!(core::mem::offset_of!(Bhed, format) == 24);
const _: () = assert!(core::mem::offset_of!(Bhed, vpol) == 58);
const _: () = assert!(core::mem::offset_of!(Bhed, hunass) == 60);

/// One entry of the `bhdr[]` table of bhdr.h.
#[derive(Clone, Copy, Debug)]
pub struct BhdrKey {
    pub key: &'static str,
    pub type_: &'static str,
    pub offs: usize,
}

const fn bk(key: &'static str, type_: &'static str, offs: usize) -> BhdrKey {
    BhdrKey { key, type_, offs }
}

/// bhdr.h table (note: hdt/dto/hns/nso are typed "h" here, as in SU).
pub static BHDR: [BhdrKey; BHED_NKEYS] = [
    bk("jobid", "i", 0),
    bk("lino", "i", 4),
    bk("reno", "i", 8),
    bk("ntrpr", "h", 12),
    bk("nart", "h", 14),
    bk("hdt", "h", 16),
    bk("dto", "h", 18),
    bk("hns", "h", 20),
    bk("nso", "h", 22),
    bk("format", "h", 24),
    bk("fold", "h", 26),
    bk("tsort", "h", 28),
    bk("vscode", "h", 30),
    bk("hsfs", "h", 32),
    bk("hsfe", "h", 34),
    bk("hslen", "h", 36),
    bk("hstyp", "h", 38),
    bk("schn", "h", 40),
    bk("hstas", "h", 42),
    bk("hstae", "h", 44),
    bk("htatyp", "h", 46),
    bk("hcorr", "h", 48),
    bk("bgrcv", "h", 50),
    bk("rcvm", "h", 52),
    bk("mfeet", "h", 54),
    bk("polyt", "h", 56),
    bk("vpol", "h", 58),
];

impl Default for Bhed {
    fn default() -> Self {
        Bhed::zeroed()
    }
}

impl Bhed {
    pub const fn zeroed() -> Self {
        // SAFETY: plain integers; all-zero is valid.
        unsafe { core::mem::zeroed() }
    }
    pub fn as_bytes(&self) -> &[u8; BNYBYTES] {
        // SAFETY: repr(C), 400 bytes, no padding.
        unsafe { &*(self as *const Bhed as *const [u8; BNYBYTES]) }
    }
    pub fn as_bytes_mut(&mut self) -> &mut [u8; BNYBYTES] {
        // SAFETY: as above; all bit patterns valid.
        unsafe { &mut *(self as *mut Bhed as *mut [u8; BNYBYTES]) }
    }
    /// Big-endian (XDR / SEG-Y tape) encoding, as `xdrbhdrsub` writes it.
    pub fn to_xdr(&self) -> [u8; BNYBYTES] {
        let s = self.as_bytes();
        let mut o = [0u8; BNYBYTES];
        for k in 0..3 {
            for b in 0..4 {
                o[4 * k + b] = s[4 * k + 3 - b];
            }
        }
        for k in (12..BNYBYTES).step_by(2) {
            o[k] = s[k + 1];
            o[k + 1] = s[k];
        }
        o
    }
    /// Decode big-endian bytes (inverse of `to_xdr`).
    pub fn from_xdr(buf: &[u8]) -> Bhed {
        let mut h = Bhed::zeroed();
        let x = {
            let mut t = [0u8; BNYBYTES];
            t.copy_from_slice(&buf[..BNYBYTES]);
            t
        };
        let tmp = Bhed::from_native(&x).to_xdr(); // byte permutation is an involution
        h.as_bytes_mut().copy_from_slice(&tmp);
        h
    }
    fn from_native(b: &[u8; BNYBYTES]) -> Bhed {
        let mut h = Bhed::zeroed();
        h.as_bytes_mut().copy_from_slice(b);
        h
    }
}

fn bh_index(index: i32) -> usize {
    if index < 0 || index as usize >= BHED_NKEYS {
        panic!("hdrpkge: binary header index {} out of range", index);
    }
    index as usize
}

/// getbhval: get a binary header word by index.
pub fn getbhval(bh: &Bhed, index: i32, valp: &mut Value) {
    let k = &BHDR[bh_index(index)];
    let b = bh.as_bytes();
    let o = k.offs;
    match k.type_ {
        "h" => valp.set_h(i16::from_ne_bytes([b[o], b[o + 1]])),
        _ => valp.set_i(i32::from_ne_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])),
    }
}

/// putbhval: put a binary header word by index.
pub fn putbhval(bh: &mut Bhed, index: i32, valp: &Value) {
    let k = &BHDR[bh_index(index)];
    let o = k.offs;
    let b = bh.as_bytes_mut();
    match k.type_ {
        "h" => b[o..o + 2].copy_from_slice(&valp.h().to_ne_bytes()),
        _ => b[o..o + 4].copy_from_slice(&valp.i().to_ne_bytes()),
    }
}

/// swapbhval: byte-swap a binary header word in place.
pub fn swapbhval(bh: &mut Bhed, index: i32) {
    let k = &BHDR[bh_index(index)];
    let w = if k.type_ == "h" { 2 } else { 4 };
    bh.as_bytes_mut()[k.offs..k.offs + w].reverse();
}

// ------------------------------------------------------- tape headers --
// tapehdr.h / tapebhdr.h tables and hdrpkge.c gettapehval & co: raw SEG-Y
// trace (240-byte) and binary (400-byte) headers as bytes, type "U" = 2-byte
// and "P" = 4-byte words read/written in NATIVE byte order (segyread swaps
// separately, as in C).

/// tapehdr.h `tapehdr[]` (SEG-Y trace header as bytes).
pub static TAPEHDR: [BhdrKey; 71] = [
    bk("tracl", "P", 0),
    bk("tracr", "P", 4),
    bk("fldr", "P", 8),
    bk("tracf", "P", 12),
    bk("ep", "P", 16),
    bk("cdp", "P", 20),
    bk("cdpt", "P", 24),
    bk("trid", "U", 28),
    bk("nvs", "U", 30),
    bk("nhs", "U", 32),
    bk("duse", "U", 34),
    bk("offset", "P", 36),
    bk("gelev", "P", 40),
    bk("selev", "P", 44),
    bk("sdepth", "P", 48),
    bk("gdel", "P", 52),
    bk("sdel", "P", 56),
    bk("swdep", "P", 60),
    bk("gwdep", "P", 64),
    bk("scalel", "U", 68),
    bk("scalco", "U", 70),
    bk("sx", "P", 72),
    bk("sy", "P", 76),
    bk("gx", "P", 80),
    bk("gy", "P", 84),
    bk("counit", "U", 88),
    bk("wevel", "U", 90),
    bk("swevel", "U", 92),
    bk("sut", "U", 94),
    bk("gut", "U", 96),
    bk("sstat", "U", 98),
    bk("gstat", "U", 100),
    bk("tstat", "U", 102),
    bk("laga", "U", 104),
    bk("lagb", "U", 106),
    bk("delrt", "U", 108),
    bk("muts", "U", 110),
    bk("mute", "U", 112),
    bk("ns", "U", 114),
    bk("dt", "U", 116),
    bk("gain", "U", 118),
    bk("igc", "U", 120),
    bk("igi", "U", 122),
    bk("corr", "U", 124),
    bk("sfs", "U", 126),
    bk("sfe", "U", 128),
    bk("slen", "U", 130),
    bk("styp", "U", 132),
    bk("stas", "U", 134),
    bk("stae", "U", 136),
    bk("tatyp", "U", 138),
    bk("afilf", "U", 140),
    bk("afils", "U", 142),
    bk("nofilf", "U", 144),
    bk("nofils", "U", 146),
    bk("lcf", "U", 148),
    bk("hcf", "U", 150),
    bk("lcs", "U", 152),
    bk("hcs", "U", 154),
    bk("year", "U", 156),
    bk("day", "U", 158),
    bk("hour", "U", 160),
    bk("minute", "U", 162),
    bk("sec", "U", 164),
    bk("timbas", "U", 166),
    bk("trwf", "U", 168),
    bk("grnors", "U", 170),
    bk("grnofr", "U", 172),
    bk("grnlof", "U", 174),
    bk("gaps", "U", 176),
    bk("otrav", "U", 178),
];

/// tapehdr.h `tapebhdr[]` (SEG-Y binary header as bytes).
pub static TAPEBHDR: [BhdrKey; 27] = [
    bk("jobid", "P", 0),
    bk("lino", "P", 4),
    bk("reno", "P", 8),
    bk("ntrpr", "U", 12),
    bk("nart", "U", 14),
    bk("hdt", "U", 16),
    bk("dto", "U", 18),
    bk("hns", "U", 20),
    bk("nso", "U", 22),
    bk("format", "U", 24),
    bk("fold", "U", 26),
    bk("tsort", "U", 28),
    bk("vscode", "U", 30),
    bk("hsfs", "U", 32),
    bk("hsfe", "U", 34),
    bk("hslen", "U", 36),
    bk("hstyp", "U", 38),
    bk("schn", "U", 40),
    bk("hstas", "U", 42),
    bk("hstae", "U", 44),
    bk("htatyp", "U", 46),
    bk("hcorr", "U", 48),
    bk("bgrcv", "U", 50),
    bk("rcvm", "U", 52),
    bk("mfeet", "U", 54),
    bk("polyt", "U", 56),
    bk("vpol", "U", 58),
];

fn tape_get(buf: &[u8], k: &BhdrKey, valp: &mut Value) {
    let o = k.offs;
    match k.type_ {
        "U" => valp.set_h(i16::from_ne_bytes([buf[o], buf[o + 1]])),
        "P" => valp.set_i(i32::from_ne_bytes([buf[o], buf[o + 1], buf[o + 2], buf[o + 3]])),
        _ => crate::err!("%s: %s: mysterious data type", "hdrpkge.c", "?"),
    }
}

fn tape_put(buf: &mut [u8], k: &BhdrKey, valp: &Value) {
    let o = k.offs;
    match k.type_ {
        "U" => buf[o..o + 2].copy_from_slice(&valp.h().to_ne_bytes()),
        "P" => buf[o..o + 4].copy_from_slice(&valp.i().to_ne_bytes()),
        _ => crate::err!("%s: %s: mysterious data type", "hdrpkge.c", "?"),
    }
}

/// gettapehval: get a tape trace header word (`tapetr`: the 240 raw bytes).
pub fn gettapehval(tapetr: &[u8], index: i32, valp: &mut Value) {
    tape_get(tapetr, &TAPEHDR[index as usize], valp)
}
/// puttapehval: put a tape trace header word.
pub fn puttapehval(tapetr: &mut [u8], index: i32, valp: &Value) {
    tape_put(tapetr, &TAPEHDR[index as usize], valp)
}
/// gettapebhval: get a tape binary header word (`tapebh`: the 400 raw bytes).
pub fn gettapebhval(tapebh: &[u8], index: i32, valp: &mut Value) {
    tape_get(tapebh, &TAPEBHDR[index as usize], valp)
}
/// puttapebhval: put a tape binary header word.
pub fn puttapebhval(tapebh: &mut [u8], index: i32, valp: &Value) {
    tape_put(tapebh, &TAPEBHDR[index as usize], valp)
}

// ---------------------------------------------------------------- text --

include!("ebcdic_tables.rs");

/// In-place `dd conv=ascii` (EBCDIC -> ASCII), what segyread runs.
pub fn dd_conv_ascii(buf: &mut [u8]) {
    for c in buf.iter_mut() {
        *c = DD_EBCDIC_TO_ASCII[*c as usize];
    }
}

/// In-place `dd conv=ebcdic` (ASCII -> EBCDIC).
pub fn dd_conv_ebcdic(buf: &mut [u8]) {
    for c in buf.iter_mut() {
        *c = DD_ASCII_TO_EBCDIC[*c as usize];
    }
}

/// cwp/lib/zebc.c ASCII -> EBCDIC table: bytes 32..=127 are mapped, every
/// other byte becomes an ASCII blank (0x20), exactly like the f2c code.
pub static ZEBC_TABLE: [u8; 256] = {
    const T2: [u8; 96] = [
        64, 79, 127, 123, 91, 108, 80, 125, 77, 93, 92, 78, 107, 96, 75, 97, 240, 241, 242, 243,
        244, 245, 246, 247, 248, 249, 122, 94, 76, 126, 110, 111, 124, 193, 194, 195, 196, 197,
        198, 199, 200, 201, 209, 210, 211, 212, 213, 214, 215, 216, 217, 226, 227, 228, 229, 230,
        231, 232, 233, 74, 224, 90, 95, 109, 64, 129, 130, 131, 132, 133, 134, 135, 136, 137, 145,
        146, 147, 148, 149, 150, 151, 152, 153, 162, 163, 164, 165, 166, 167, 168, 169, 192, 106,
        208, 161, 64,
    ];
    let mut t = [0x20u8; 256];
    let mut i = 32;
    while i < 128 {
        t[i] = T2[i - 32];
        i += 1;
    }
    t
};

/// zebc: translate `nchar` characters from ASCII to EBCDIC (segywrite).
pub fn zebc(ainput: &[u8], aoutput: &mut [u8], nchar: usize) {
    for i in 0..nchar {
        aoutput[i] = ZEBC_TABLE[ainput[i] as usize];
    }
}

// ------------------------------------------------------------ IBM float --

fn ibm_tables() -> &'static ([u32; 512], [f32; 512]) {
    use std::sync::OnceLock;
    static T: OnceLock<([u32; 512], [f32; 512])> = OnceLock::new();
    T.get_or_init(|| {
        let mut m1 = [0u32; 512];
        let mut r1 = [0f32; 512];
        let r_infinity = f32::from_bits(2139095039); // 0x7F7FFFFF = FLT_MAX
        for i in 0..512u32 {
            let i1 = i & 255;
            let eibm = i1 >> 1;
            let mhibm = i & 1;
            let et3e = ((eibm << 2) as i32) - 130;
            let iu = i as usize;
            if et3e > 0 && et3e <= 255 {
                let x = (et3e as u32) ^ (i & 255);
                m1[iu] = x.wrapping_shl(23);
                if mhibm == 1 {
                    r1[iu] = 0.0;
                } else {
                    let y = (et3e as u32) | (i & 256);
                    r1[iu] = -f32::from_bits(y.wrapping_shl(23));
                }
            } else if et3e <= 0 {
                m1[iu] = i.wrapping_shl(23);
                r1[iu] = 0.0;
            } else {
                m1[iu] = i.wrapping_shl(23);
                r1[iu] = if i < 256 { r_infinity } else { -r_infinity };
            }
        }
        (m1, r1)
    })
}

/// Convert one 32-bit IBM word (as loaded natively from memory) to IEEE,
/// exactly like segyread.c `ibm_to_float` (endian 0 = byte-reverse first).
pub fn ibm_to_float_word(word: u32, endian: i32) -> f32 {
    let (m1, r1) = ibm_tables();
    let w = if endian == 0 { word.swap_bytes() } else { word };
    let jj = (w >> 23) as usize;
    let bits = m1[jj] ^ w;
    f32::from_bits(bits) + r1[jj]
}

/// segyread.c `ibm_to_float`, in place: each f32 slot of `data[..n]` holds
/// a raw IBM word on input and the IEEE float on output.
pub fn ibm_to_float(data: &mut [f32], n: usize, endian: i32) {
    for x in data[..n].iter_mut() {
        *x = ibm_to_float_word(x.to_bits(), endian);
    }
}

/// Convert one IEEE float's bits to an IBM word, exactly like segywrite.c
/// `float_to_ibm` (endian 0 = byte-reverse the result).
pub fn float_to_ibm_word(from: u32, endian: i32) -> u32 {
    let mut fconv = from;
    if fconv != 0 {
        let mut fmant = (0x007f_ffff & fconv) | 0x0080_0000;
        let mut t: i32 = ((0x7f80_0000 & fconv) >> 23) as i32 - 126;
        while t & 0x3 != 0 {
            t += 1;
            fmant >>= 1;
        }
        fconv = (0x8000_0000 & fconv) | (((t >> 2) + 64) as u32).wrapping_shl(24) | fmant;
    }
    if endian == 0 {
        fconv = fconv.swap_bytes();
    }
    fconv
}

/// segywrite.c `float_to_ibm`, in place on `data[..n]` (result words are
/// stored as f32 bit patterns).
pub fn float_to_ibm(data: &mut [f32], n: usize, endian: i32) {
    for x in data[..n].iter_mut() {
        *x = f32::from_bits(float_to_ibm_word(x.to_bits(), endian));
    }
}

// ------------------------------------------------------------- swapbyte --

pub fn swap_short_2(x: &mut i16) {
    *x = x.swap_bytes();
}
pub fn swap_u_short_2(x: &mut u16) {
    *x = x.swap_bytes();
}
pub fn swap_int_4(x: &mut i32) {
    *x = x.swap_bytes();
}
pub fn swap_u_int_4(x: &mut u32) {
    *x = x.swap_bytes();
}
/// swap_long_4 on LP64: swaps the low 4 bytes and clears the high 4 (C masks).
pub fn swap_long_4(x: &mut i64) {
    *x = ((*x as u32).swap_bytes()) as i64;
}
/// swap_u_long_4 on LP64: as swap_long_4.
pub fn swap_u_long_4(x: &mut u64) {
    *x = ((*x as u32).swap_bytes()) as u64;
}
pub fn swap_float_4(x: &mut f32) {
    *x = f32::from_bits(x.to_bits().swap_bytes());
}
pub fn swap_double_8(x: &mut f64) {
    *x = f64::from_bits(x.to_bits().swap_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ibm_be(x: u32) -> f32 {
        // word as it sits in a big-endian file, loaded natively on LE
        ibm_to_float_word(x.swap_bytes(), 0)
    }

    #[test]
    fn ibm_known_values() {
        assert_eq!(ibm_be(0x0000_0000), 0.0);
        assert_eq!(ibm_be(0x4110_0000), 1.0);
        assert_eq!(ibm_be(0xC110_0000), -1.0);
        assert_eq!(ibm_be(0x4220_0000), 32.0);
        assert_eq!(ibm_be(0x4080_0000), 0.5);
        assert_eq!(ibm_be(0xC276_A000), -118.625);
        // too large -> FLT_MAX (not inf)
        assert_eq!(ibm_be(0x7FFF_FFFF), f32::MAX);
        assert_eq!(ibm_be(0xFFFF_FFFF), -f32::MAX);
        // tiny IBM values pass through the table unchanged (as in C)
        assert_eq!(ibm_be(0x0010_0000), f32::from_bits(0x0010_0000));
    }

    #[test]
    fn float_to_ibm_known() {
        assert_eq!(float_to_ibm_word(1.0f32.to_bits(), 1), 0x4110_0000);
        assert_eq!(float_to_ibm_word((-118.625f32).to_bits(), 1), 0xC276_A000);
        assert_eq!(float_to_ibm_word(0, 1), 0);
        assert_eq!(float_to_ibm_word(1.0f32.to_bits(), 0), 0x0000_1041);
        // round trip of exactly representable values
        for &x in &[1.0f32, -2.5, 0.15625, 1024.0, -3.0e-5, 7.0e20] {
            let w = float_to_ibm_word(x.to_bits(), 1);
            let back = ibm_to_float_word(w, 1);
            assert!((back - x).abs() <= x.abs() * 1e-6, "{} {}", x, back);
        }
    }

    #[test]
    fn bhed_xdr() {
        let mut bh = Bhed::zeroed();
        bh.jobid = 1;
        bh.hns = 1000;
        bh.format = 1;
        bh.hunass[169] = 0x0102;
        let x = bh.to_xdr();
        assert_eq!(&x[0..4], &[0, 0, 0, 1]);
        assert_eq!(&x[20..22], &1000u16.to_be_bytes());
        assert_eq!(&x[398..400], &[1, 2]);
        assert_eq!(Bhed::from_xdr(&x), bh);
        let mut v = Value::default();
        getbhval(&bh, 9, &mut v);
        assert_eq!(v.h(), 1);
    }

    #[test]
    fn tape_tables() {
        for (t, h) in TAPEHDR.iter().zip(crate::header::HDR.iter()) {
            assert_eq!((t.key, t.offs), (h.key, h.offs));
        }
        for (t, b) in TAPEBHDR.iter().zip(BHDR.iter()) {
            assert_eq!((t.key, t.offs), (b.key, b.offs));
        }
        let mut raw = [0u8; 240];
        raw[114..116].copy_from_slice(&1000i16.to_ne_bytes());
        let mut v = Value::default();
        gettapehval(&raw, 38, &mut v);
        assert_eq!(v.h(), 1000);
        v.set_i(-5);
        puttapehval(&mut raw, 0, &v);
        assert_eq!(&raw[0..4], &(-5i32).to_ne_bytes());
    }

    #[test]
    fn text_tables() {
        let mut s = *b"C 1 CLIENT";
        dd_conv_ebcdic(&mut s);
        let e = s;
        dd_conv_ascii(&mut s);
        assert_eq!(&s, b"C 1 CLIENT");
        let mut z = [0u8; 10];
        zebc(b"C 1 CLIENT", &mut z, 10);
        assert_eq!(z, e);
        assert_eq!(ZEBC_TABLE[b'\n' as usize], 0x20);
    }
}
