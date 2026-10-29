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

//! `header`: the SU trace header (`segy.h`), the key table (`hdr.h`) and the
//! name/index accessors of `hdrpkge.c`.
//!
//! C -> Rust:
//! | C                         | Rust                                   |
//! |---------------------------|----------------------------------------|
//! | `segy tr; tr.offset`      | `Trace` derefs to `SegyHeader`: `tr.offset` |
//! | `hdr[i].key/.type/.offs`  | `HDR[i].key / .type_ / .offs`          |
//! | `gethval(&tr, i, &val)`   | `gethval(&tr, i, &mut val)`            |
//! | `puthval(&tr, i, &val)`   | `puthval(&mut tr, i, &val)`            |
//! | `gethdval(&tr, key, &val)`| `gethdval(&tr, key, &mut val)`         |
//! | `puthdval(&tr, key, &val)`| `puthdval(&mut tr, key, &val)`         |
//! | `hdtype(key)`             | `hdtype(key) -> &'static str`          |
//! | `getkey(i)`               | `getkey(i) -> Option<&'static str>`    |
//! | `getindex(key)`           | `getindex(key) -> i32` (-1 if absent)  |
//! | `swaphval(&tr, i)`        | `swaphval(&mut tr, i)`                 |
//! | `printheader(&tr)`        | `printheader(&tr)`                     |
//! | `SU_NKEYS, HDRBYTES, SU_NFLTS` | same names                        |

use crate::value::{printfval, valcmp, Value};

/// Arbitrary limit on data array size (segy.h).
pub const SU_NFLTS: usize = 65535;
/// Number of key header words (header.h).
pub const SU_NKEYS: usize = 80;
/// Bytes in the trace header (header.h).
pub const HDRBYTES: usize = 240;
/// header.h MAXSEGY.
pub const MAXSEGY: usize = 262380;

/// SU trace identification header, field-for-field identical to `segy.h`
/// (the non-SLTSU branch). 240 bytes, `#[repr(C)]`.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct SegyHeader {
    pub tracl: i32,
    pub tracr: i32,
    pub fldr: i32,
    pub tracf: i32,
    pub ep: i32,
    pub cdp: i32,
    pub cdpt: i32,
    pub trid: i16,
    pub nvs: i16,
    pub nhs: i16,
    pub duse: i16,
    pub offset: i32,
    pub gelev: i32,
    pub selev: i32,
    pub sdepth: i32,
    pub gdel: i32,
    pub sdel: i32,
    pub swdep: i32,
    pub gwdep: i32,
    pub scalel: i16,
    pub scalco: i16,
    pub sx: i32,
    pub sy: i32,
    pub gx: i32,
    pub gy: i32,
    pub counit: i16,
    pub wevel: i16,
    pub swevel: i16,
    pub sut: i16,
    pub gut: i16,
    pub sstat: i16,
    pub gstat: i16,
    pub tstat: i16,
    pub laga: i16,
    pub lagb: i16,
    pub delrt: i16,
    pub muts: i16,
    pub mute: i16,
    pub ns: u16,
    pub dt: u16,
    pub gain: i16,
    pub igc: i16,
    pub igi: i16,
    pub corr: i16,
    pub sfs: i16,
    pub sfe: i16,
    pub slen: i16,
    pub styp: i16,
    pub stas: i16,
    pub stae: i16,
    pub tatyp: i16,
    pub afilf: i16,
    pub afils: i16,
    pub nofilf: i16,
    pub nofils: i16,
    pub lcf: i16,
    pub hcf: i16,
    pub lcs: i16,
    pub hcs: i16,
    pub year: i16,
    pub day: i16,
    pub hour: i16,
    pub minute: i16,
    pub sec: i16,
    pub timbas: i16,
    pub trwf: i16,
    pub grnors: i16,
    pub grnofr: i16,
    pub grnlof: i16,
    pub gaps: i16,
    pub otrav: i16,
    pub d1: f32,
    pub f1: f32,
    pub d2: f32,
    pub f2: f32,
    pub ungpow: f32,
    pub unscale: f32,
    pub ntr: i32,
    pub mark: i16,
    pub shortpad: i16,
    pub unass: [i16; 14],
}

/// One entry of the `hdr[]` table of hdr.h.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HdrKey {
    pub key: &'static str,
    /// SU type code: "i" int, "h" short, "u" unsigned short, "f" float
    /// (the table also allows "p","l","v","d","s" though segy.h uses none).
    pub type_: &'static str,
    pub offs: usize,
}

const fn k(key: &'static str, type_: &'static str, offs: usize) -> HdrKey {
    HdrKey { key, type_, offs }
}

/// The `hdr[]` table of hdr.h, in the same order (index = SU key index).
pub static HDR: [HdrKey; SU_NKEYS] = [
    k("tracl", "i", 0),
    k("tracr", "i", 4),
    k("fldr", "i", 8),
    k("tracf", "i", 12),
    k("ep", "i", 16),
    k("cdp", "i", 20),
    k("cdpt", "i", 24),
    k("trid", "h", 28),
    k("nvs", "h", 30),
    k("nhs", "h", 32),
    k("duse", "h", 34),
    k("offset", "i", 36),
    k("gelev", "i", 40),
    k("selev", "i", 44),
    k("sdepth", "i", 48),
    k("gdel", "i", 52),
    k("sdel", "i", 56),
    k("swdep", "i", 60),
    k("gwdep", "i", 64),
    k("scalel", "h", 68),
    k("scalco", "h", 70),
    k("sx", "i", 72),
    k("sy", "i", 76),
    k("gx", "i", 80),
    k("gy", "i", 84),
    k("counit", "h", 88),
    k("wevel", "h", 90),
    k("swevel", "h", 92),
    k("sut", "h", 94),
    k("gut", "h", 96),
    k("sstat", "h", 98),
    k("gstat", "h", 100),
    k("tstat", "h", 102),
    k("laga", "h", 104),
    k("lagb", "h", 106),
    k("delrt", "h", 108),
    k("muts", "h", 110),
    k("mute", "h", 112),
    k("ns", "u", 114),
    k("dt", "u", 116),
    k("gain", "h", 118),
    k("igc", "h", 120),
    k("igi", "h", 122),
    k("corr", "h", 124),
    k("sfs", "h", 126),
    k("sfe", "h", 128),
    k("slen", "h", 130),
    k("styp", "h", 132),
    k("stas", "h", 134),
    k("stae", "h", 136),
    k("tatyp", "h", 138),
    k("afilf", "h", 140),
    k("afils", "h", 142),
    k("nofilf", "h", 144),
    k("nofils", "h", 146),
    k("lcf", "h", 148),
    k("hcf", "h", 150),
    k("lcs", "h", 152),
    k("hcs", "h", 154),
    k("year", "h", 156),
    k("day", "h", 158),
    k("hour", "h", 160),
    k("minute", "h", 162),
    k("sec", "h", 164),
    k("timbas", "h", 166),
    k("trwf", "h", 168),
    k("grnors", "h", 170),
    k("grnofr", "h", 172),
    k("grnlof", "h", 174),
    k("gaps", "h", 176),
    k("otrav", "h", 178),
    k("d1", "f", 180),
    k("f1", "f", 184),
    k("d2", "f", 188),
    k("f2", "f", 192),
    k("ungpow", "f", 196),
    k("unscale", "f", 200),
    k("ntr", "i", 204),
    k("mark", "h", 208),
    k("shortpad", "h", 210),
];

// ---- static layout checks against hdr.h ----------------------------------
macro_rules! assert_offs {
    ($($f:ident = $o:expr),* $(,)?) => {
        $( const _: () = assert!(core::mem::offset_of!(SegyHeader, $f) == $o); )*
    };
}
const _: () = assert!(core::mem::size_of::<SegyHeader>() == HDRBYTES);
assert_offs!(
    tracl = 0, tracr = 4, fldr = 8, tracf = 12, ep = 16, cdp = 20, cdpt = 24,
    trid = 28, nvs = 30, nhs = 32, duse = 34, offset = 36, gelev = 40,
    selev = 44, sdepth = 48, gdel = 52, sdel = 56, swdep = 60, gwdep = 64,
    scalel = 68, scalco = 70, sx = 72, sy = 76, gx = 80, gy = 84, counit = 88,
    wevel = 90, swevel = 92, sut = 94, gut = 96, sstat = 98, gstat = 100,
    tstat = 102, laga = 104, lagb = 106, delrt = 108, muts = 110, mute = 112,
    ns = 114, dt = 116, gain = 118, igc = 120, igi = 122, corr = 124,
    sfs = 126, sfe = 128, slen = 130, styp = 132, stas = 134, stae = 136,
    tatyp = 138, afilf = 140, afils = 142, nofilf = 144, nofils = 146,
    lcf = 148, hcf = 150, lcs = 152, hcs = 154, year = 156, day = 158,
    hour = 160, minute = 162, sec = 164, timbas = 166, trwf = 168,
    grnors = 170, grnofr = 172, grnlof = 174, gaps = 176, otrav = 178,
    d1 = 180, f1 = 184, d2 = 188, f2 = 192, ungpow = 196, unscale = 200,
    ntr = 204, mark = 208, shortpad = 210, unass = 212,
);

impl SegyHeader {
    /// An all-zero header (like a static C `segy`).
    pub const fn zeroed() -> Self {
        // SAFETY: all fields are plain integers/floats; all-zero is valid.
        unsafe { core::mem::zeroed() }
    }

    /// Raw native-endian bytes of the header (the C struct memory).
    pub fn as_bytes(&self) -> &[u8; HDRBYTES] {
        // SAFETY: repr(C), size 240, no padding (checked above), any byte valid.
        unsafe { &*(self as *const SegyHeader as *const [u8; HDRBYTES]) }
    }

    /// Mutable raw native-endian bytes of the header.
    pub fn as_bytes_mut(&mut self) -> &mut [u8; HDRBYTES] {
        // SAFETY: as above; every bit pattern is a valid SegyHeader.
        unsafe { &mut *(self as *mut SegyHeader as *mut [u8; HDRBYTES]) }
    }

    /// Encode as XDR (big-endian per field), exactly what `xdrhdrsub`
    /// produces in this SUXDR build.
    pub fn to_xdr(&self) -> [u8; HDRBYTES] {
        let src = self.as_bytes();
        let mut out = [0u8; HDRBYTES];
        for kk in HDR.iter() {
            let w = if kk.type_ == "h" || kk.type_ == "u" { 2 } else { 4 };
            for b in 0..w {
                out[kk.offs + b] = src[kk.offs + w - 1 - b];
            }
        }
        for j in 0..14 {
            let o = 212 + 2 * j;
            out[o] = src[o + 1];
            out[o + 1] = src[o];
        }
        out
    }

    /// Decode from XDR (big-endian) bytes; inverse of `to_xdr`.
    pub fn from_xdr(buf: &[u8]) -> Self {
        let mut h = SegyHeader::zeroed();
        let dst = h.as_bytes_mut();
        for kk in HDR.iter() {
            let w = if kk.type_ == "h" || kk.type_ == "u" { 2 } else { 4 };
            for b in 0..w {
                dst[kk.offs + b] = buf[kk.offs + w - 1 - b];
            }
        }
        for j in 0..14 {
            let o = 212 + 2 * j;
            dst[o] = buf[o + 1];
            dst[o + 1] = buf[o];
        }
        h
    }
}

fn index_ok(index: i32) -> usize {
    if index < 0 || index as usize >= SU_NKEYS {
        panic!("hdrpkge: header index {} out of range", index);
    }
    index as usize
}

fn get_at(tr: &SegyHeader, idx: usize, valp: &mut Value) {
    let b = tr.as_bytes();
    let o = HDR[idx].offs;
    match HDR[idx].type_.as_bytes()[0] {
        b'h' => valp.set_h(i16::from_ne_bytes([b[o], b[o + 1]])),
        b'u' => valp.set_u(u16::from_ne_bytes([b[o], b[o + 1]])),
        b'i' => valp.set_i(i32::from_ne_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])),
        b'f' => valp.set_f(f32::from_ne_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])),
        _ => crate::err!("%s: %s: mysterious data type", "hdrpkge.c", "?"),
    }
}

fn put_at(tr: &mut SegyHeader, idx: usize, valp: &Value) {
    let o = HDR[idx].offs;
    let b = tr.as_bytes_mut();
    match HDR[idx].type_.as_bytes()[0] {
        b'h' => b[o..o + 2].copy_from_slice(&valp.h().to_ne_bytes()),
        b'u' => b[o..o + 2].copy_from_slice(&valp.u().to_ne_bytes()),
        b'i' => b[o..o + 4].copy_from_slice(&valp.i().to_ne_bytes()),
        b'f' => b[o..o + 4].copy_from_slice(&valp.f().to_ne_bytes()),
        _ => crate::err!("%s: %s: mysterious data type", "hdrpkge.c", "?"),
    }
}

/// gethval: get a trace header word by index.
pub fn gethval(tr: &SegyHeader, index: i32, valp: &mut Value) {
    get_at(tr, index_ok(index), valp)
}

/// puthval: put a trace header word by index.
pub fn puthval(tr: &mut SegyHeader, index: i32, valp: &Value) {
    put_at(tr, index_ok(index), valp)
}

/// gethdval: get a trace header word by name.
/// Errors exactly like SU: `hdrpkge.c: key word not in segy.h: 'xx'`.
pub fn gethdval(tr: &SegyHeader, key: &str, valp: &mut Value) {
    let index = getindex(key);
    if index == -1 {
        crate::err!("%s: key word not in segy.h: '%s'", "hdrpkge.c", key);
    }
    get_at(tr, index as usize, valp)
}

/// puthdval: put a trace header word by name.
pub fn puthdval(tr: &mut SegyHeader, key: &str, valp: &Value) {
    let index = getindex(key);
    if index == -1 {
        crate::err!("%s: key word not in segy.h: '%s'", "hdrpkge.c", key);
    }
    put_at(tr, index as usize, valp)
}

/// hdtype: data type string of a header word by name ("i", "h", "u", "f").
pub fn hdtype(key: &str) -> &'static str {
    let index = getindex(key);
    if index == -1 {
        crate::err!("%s: key word not in segy.h: '%s'", "hdrpkge.c", key);
    }
    HDR[index as usize].type_
}

/// getkey: name of a header word from its index (None for out of range, like NULL).
pub fn getkey(index: i32) -> Option<&'static str> {
    if index >= 0 && (index as usize) < SU_NKEYS {
        Some(HDR[index as usize].key)
    } else {
        None
    }
}

/// getindex: index of a header word from its name, -1 if not found.
pub fn getindex(key: &str) -> i32 {
    for (i, kk) in HDR.iter().enumerate() {
        if kk.key == key {
            return i as i32;
        }
    }
    -1
}

/// swaphval: byte-swap the header word with this index in place.
pub fn swaphval(tr: &mut SegyHeader, index: i32) {
    let idx = index_ok(index);
    let o = HDR[idx].offs;
    let w = match HDR[idx].type_.as_bytes()[0] {
        b'h' | b'u' => 2,
        b'i' | b'p' | b'f' => 4,
        _ => crate::err!("%s: %s: unsupported data type", "hdrpkge.c", "?"),
    };
    tr.as_bytes_mut()[o..o + w].reverse();
}

/// printheader: display non-null header field values on stdout
/// (" key=val", six per line, then a newline).
pub fn printheader(tp: &SegyHeader) {
    let mut j = 0;
    let mut zeroval = Value::default();
    zeroval.set_l(0);
    for i in 0..SU_NKEYS as i32 {
        let mut val = Value::default();
        gethval(tp, i, &mut val);
        let key = getkey(i).unwrap();
        let type_ = hdtype(key);
        if valcmp(type_, val, zeroval) != 0 {
            crate::printf!(" %s=", key);
            printfval(type_, val);
            j += 1;
            if j % 6 == 0 {
                crate::printf!("\n");
            }
        }
    }
    crate::printf!("\n");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_matches_struct_offsets() {
        let mut h = SegyHeader::zeroed();
        h.ns = 0x1234;
        h.offset = -400;
        h.d1 = 1.5;
        let mut v = Value::default();
        gethdval(&h, "ns", &mut v);
        assert_eq!(v.u(), 0x1234);
        gethdval(&h, "offset", &mut v);
        assert_eq!(v.i(), -400);
        gethdval(&h, "d1", &mut v);
        assert_eq!(v.f(), 1.5);
        assert_eq!(getindex("shortpad"), 79);
        assert_eq!(getindex("nope"), -1);
        assert_eq!(getkey(80), None);
        assert_eq!(hdtype("dt"), "u");
        // contiguous coverage of 0..212
        let mut end = 0;
        for kk in HDR.iter() {
            assert_eq!(kk.offs, end, "{}", kk.key);
            end += if kk.type_ == "i" || kk.type_ == "f" { 4 } else { 2 };
        }
        assert_eq!(end, 212);
    }

    #[test]
    fn xdr_roundtrip_every_field() {
        let mut h = SegyHeader::zeroed();
        for (i, b) in h.as_bytes_mut().iter_mut().enumerate() {
            *b = (i * 7 + 3) as u8;
        }
        let x = h.to_xdr();
        let back = SegyHeader::from_xdr(&x);
        assert_eq!(back.as_bytes(), h.as_bytes());
        let mut h2 = SegyHeader::zeroed();
        h2.offset = 400;
        h2.ns = 64;
        h2.d1 = 2.0;
        h2.unass[13] = -2;
        let x2 = h2.to_xdr();
        assert_eq!(&x2[36..40], &[0, 0, 1, 0x90]);
        assert_eq!(&x2[114..116], &[0, 64]);
        assert_eq!(&x2[180..184], &2.0f32.to_be_bytes());
        assert_eq!(&x2[238..240], &[0xff, 0xfe]);
    }

    #[test]
    fn put_and_swap() {
        let mut h = SegyHeader::zeroed();
        let mut v = Value::default();
        v.set_h(-5);
        puthdval(&mut h, "trid", &v);
        assert_eq!(h.trid, -5);
        v.set_i(0x01020304);
        puthval(&mut h, 0, &v);
        swaphval(&mut h, 0);
        assert_eq!(h.tracl, 0x04030201);
    }
}
