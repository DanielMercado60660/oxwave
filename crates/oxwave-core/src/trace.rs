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

//! `trace`: SU trace I/O (su/lib/fgettr.c, su/lib/fputtr.c), XDR branch.
//!
//! This SU build defines `SUXDR`: every SU trace stream is **big-endian**
//! (header fields at their segy.h offsets, samples as big-endian IEEE
//! floats), even on this little-endian machine. Native-endian input is
//! rejected with "fgettr.c: no data on first trace" (or misread).
//!
//! | C                          | Rust                                  |
//! |----------------------------|---------------------------------------|
//! | `segy tr;` (global/static) | `let mut tr = Trace::new();` (data = SU_NFLTS zeros) |
//! | `tr.ns`, `tr.offset`, ...  | `tr.ns`, `tr.offset` (Deref to `SegyHeader`), or `tr.hdr.ns` |
//! | `tr.data[i]`               | `tr.data[i]` (`Vec<f32>`, length SU_NFLTS) |
//! | `(char*)tr.data`, `(short*)tr.data` | `tr.data_bytes()`, `tr.data_bytes_mut()` |
//! | `gettr(&tr)`               | `gettr(&mut tr)` -> i32 bytes read, 0 at end |
//! | `fgettr(fp, &tr)`          | `fgettr(&fp, &mut tr)`                |
//! | `vgettr` / `fvgettr`       | `vgettr` / `fvgettr` (variable ns)    |
//! | `puttr(&tr)`               | `puttr(&mut tr)` (&mut: SHORTPACK data is byte-swapped in place, as in C) |
//! | `fputtr(fp, &tr)`          | `fputtr(&fp, &mut tr)`                |
//! | `vputtr` / `fvputtr`       | `vputtr` / `fvputtr`                  |
//! | `gettra(&tr, itr)`         | `gettra(&mut tr, itr)` -> ntr         |
//! | `fgettra(fp, &tr, itr)`    | `fgettra(&fp, &mut tr, itr)`          |
//!
//! Semantics (identical to the C):
//! * The sample count of the **first** trace fixes the record length of
//!   the stream. With gettr, every later trace reads nsfirst samples and
//!   then errs if its header ns differs:
//!   `fgettr.c: on trace #<itr>, number of samples in header (<ns>) differs from number for first trace (<nsfirst>)`
//!   (itr = number of traces already read).
//! * Empty input or a short first header: `fgettr.c: bad first header`.
//!   First trace with missing/short data: `fgettr.c: no data on first trace`.
//! * Later traces: a short header ends input silently (returns 0); a short
//!   float data block returns 240 with the previous samples left in place
//!   (and the next call returns 0).
//! * stdin a tty: `fgettr.c: segy input can't be tty`; stdout a tty:
//!   `fputtr.c: segy output can't be tty`.
//! * fputtr uses nsfirst (first trace's ns) for data output, whatever later
//!   headers say; fvputtr uses each trace's own ns.
//!   The first trace's trid also fixes bytes per sample: later float traces
//!   are XDR encoded in full, but only that byte count's prefix is written.
//! * The state is per stream (a `CFile`), separately for reading and
//!   writing, so `fputtr(tmp) .. erewind(tmp) .. fgettr(tmp)` works.
//!
//! Reopening a file starts fresh state here. SU never removes its FILE-pointer
//! entries after fclose, so Darwin libc reusing a FILE slot can make a newly
//! opened file inherit an unrelated file's trace state. We do not emulate that
//! slot allocation: CFile has its own identity, and neither a reused file
//! descriptor nor a reused Rust allocation identifies a Darwin FILE slot.

use crate::cio::{efread, efseeko, efwrite, eftello, fileno, filestat, stdin, stdout, CFile, FileType, TrIn, TrOut, SEEK_END, SEEK_SET};
use crate::header::{SegyHeader, HDRBYTES, SU_NFLTS};

/// trid of byte-packed traces (supack1).
pub const CHARPACK: i16 = 201;
/// trid of short-packed traces (supack2).
pub const SHORTPACK: i16 = 202;
/// BYTES_PER_XDR_UNIT.
const XDR_UNIT: u16 = 4;

/// An SU trace: header plus SU_NFLTS samples (like C `segy`).
#[derive(Clone, Debug, PartialEq)]
pub struct Trace {
    pub hdr: SegyHeader,
    pub data: Vec<f32>,
}

impl Default for Trace {
    fn default() -> Self {
        Trace::new()
    }
}

impl Trace {
    /// A zeroed trace with `SU_NFLTS` samples (the C `segy tr;`).
    pub fn new() -> Trace {
        Trace { hdr: SegyHeader::zeroed(), data: vec![0.0; SU_NFLTS] }
    }

    /// The data array viewed as bytes (C `(char *) tr.data`).
    pub fn data_bytes(&self) -> &[u8] {
        crate::cio::as_bytes(&self.data[..])
    }

    /// The data array viewed as mutable bytes.
    pub fn data_bytes_mut(&mut self) -> &mut [u8] {
        crate::cio::as_bytes_mut(&mut self.data[..])
    }

    fn ensure_bytes(&mut self, nbytes: usize) {
        let nf = nbytes.div_ceil(4);
        if self.data.len() < nf {
            self.data.resize(nf, 0.0);
        }
    }
}

impl std::ops::Deref for Trace {
    type Target = SegyHeader;
    fn deref(&self) -> &SegyHeader {
        &self.hdr
    }
}
impl std::ops::DerefMut for Trace {
    fn deref_mut(&mut self) -> &mut SegyHeader {
        &mut self.hdr
    }
}

/// C `swab` on the first n bytes (pairs swapped).
fn swab(b: &mut [u8], n: usize) {
    let n = n.min(b.len()) & !1;
    for k in (0..n).step_by(2) {
        b.swap(k, k + 1);
    }
}

fn bytes_per(trid: i16) -> u16 {
    match trid {
        CHARPACK => 1,
        SHORTPACK => 2,
        _ => XDR_UNIT,
    }
}

/// dataread() of fgettr.c (XDR branch). Returns bytes read.
fn dataread(fp: &CFile, st: &TrIn, tp: &mut Trace, fixed_length: bool) -> i32 {
    let nsread: u32 = if fixed_length { st.nsfirst as u32 } else { tp.hdr.ns as u32 };
    let databytes = st.bytesper as usize * nsread as usize;
    tp.ensure_bytes(databytes);
    let mut nread: i32;
    match tp.hdr.trid {
        CHARPACK | SHORTPACK => {
            if tp.hdr.trid == CHARPACK {
                // C falls through from CHARPACK into SHORTPACK: data read twice.
                efread(&mut tp.data_bytes_mut()[..databytes], fp);
            }
            nread = efread(&mut tp.data_bytes_mut()[..databytes], fp) as i32;
            swab(tp.data_bytes_mut(), databytes);
        }
        _ => {
            let mut buf = vec![0u8; databytes];
            nread = efread(&mut buf[..], fp) as i32;
            if nread as usize != databytes {
                nread = 0;
            } else {
                for (k, c) in buf.chunks_exact(4).enumerate() {
                    tp.data[k] = f32::from_bits(u32::from_be_bytes([c[0], c[1], c[2], c[3]]));
                }
                nread = databytes as i32;
            }
        }
    }
    if nread > 0 && nread as usize != databytes {
        crate::err!(
            "%s: on trace #%ld, tried to read %d bytes, read %d bytes",
            "fgettr.c",
            st.itr + 1,
            databytes as i32,
            nread
        );
    }
    nread
}

fn read_header(fp: &CFile) -> Option<SegyHeader> {
    let mut buf = [0u8; HDRBYTES];
    let nread = efread(&mut buf[..], fp);
    if nread != HDRBYTES {
        None
    } else {
        Some(SegyHeader::from_xdr(&buf))
    }
}

fn fgettr_internal(fp: &CFile, tp: &mut Trace, fixed_length: bool) -> i32 {
    let st = fp.with(|i| i.trin.clone());
    let mut nread: i32;
    match st {
        None => {
            let ftype = filestat(fileno(fp));
            match ftype {
                FileType::Directory => crate::err!("%s: segy input can't be a directory", "fgettr.c"),
                FileType::Tty => crate::err!("%s: segy input can't be tty", "fgettr.c"),
                _ => {}
            }
            match read_header(fp) {
                Some(h) => tp.hdr = h,
                None => crate::err!("%s: bad first header", "fgettr.c"),
            }
            let nsfirst = tp.hdr.ns as i32;
            if nsfirst > SU_NFLTS as i32 {
                crate::err!("%s: unable to handle %d > %d samples per trace", "fgettr.c", nsfirst, SU_NFLTS as i32);
            }
            let bytesper = bytes_per(tp.hdr.trid);
            let databytes = bytesper as u32 * tp.hdr.ns as u32;
            let mut st = TrIn {
                itr: 0,
                nsfirst,
                bytesper,
                nsegy: databytes as usize + HDRBYTES,
                ntr: -1,
                ftype,
            };
            fp.with(|i| i.trin = Some(st.clone()));
            nread = dataread(fp, &st, tp, fixed_length);
            match nread {
                0 => crate::err!("%s: no data on first trace", "fgettr.c"),
                _ => {
                    if nread as u32 != databytes {
                        crate::err!(
                            "%s: first trace: tried to read %d bytes read %d bytes",
                            "fgettr.c",
                            databytes,
                            nread
                        );
                    } else {
                        nread += HDRBYTES as i32;
                    }
                }
            }
            if ftype == FileType::Disk {
                let curoff = eftello(fp);
                efseeko(fp, 0, SEEK_END);
                st.ntr = eftello(fp) / st.nsegy as i64;
                efseeko(fp, curoff, SEEK_SET);
            }
            st.itr += 1;
            fp.with(|i| i.trin = Some(st));
            nread
        }
        Some(mut st) => {
            match read_header(fp) {
                Some(h) => {
                    tp.hdr = h;
                    nread = HDRBYTES as i32 + dataread(fp, &st, tp, fixed_length);
                }
                None => nread = 0,
            }
            if nread <= 0 {
                return 0;
            }
            if fixed_length && tp.hdr.ns as i32 != st.nsfirst {
                crate::err!(
                    "%s: on trace #%ld, number of samples in header (%d) differs from number for first trace (%d)",
                    "fgettr.c",
                    st.itr,
                    tp.hdr.ns,
                    st.nsfirst
                );
            }
            st.itr += 1;
            fp.with(|i| i.trin = Some(st));
            nread
        }
    }
}

/// fgettr: get a fixed-length trace from a stream. Returns bytes read, 0 at end.
pub fn fgettr(fp: &CFile, tp: &mut Trace) -> i32 {
    fgettr_internal(fp, tp, true)
}
/// fvgettr: get a variable-length trace from a stream.
pub fn fvgettr(fp: &CFile, tp: &mut Trace) -> i32 {
    fgettr_internal(fp, tp, false)
}
/// gettr: fgettr(stdin).
pub fn gettr(tp: &mut Trace) -> i32 {
    fgettr(&stdin(), tp)
}
/// vgettr: fvgettr(stdin).
pub fn vgettr(tp: &mut Trace) -> i32 {
    fvgettr(&stdin(), tp)
}

/// fgettra: get trace number `itr` (0-based) from a disk file; returns ntr.
pub fn fgettra(fp: &CFile, tp: &mut Trace, itr: i32) -> i32 {
    if fp.with(|i| i.trin.is_none()) {
        if 0 >= fgettr(fp, tp) {
            return 0;
        }
        let ftype = fp.with(|i| i.trin.as_ref().unwrap().ftype);
        match ftype {
            FileType::Tty => crate::warn!("stdin not redirected"),
            FileType::Disk => {}
            _ => crate::err!("%s: input must be disk file", "fgettr.c"),
        }
        efseeko(fp, 0, SEEK_END);
        let end = eftello(fp);
        fp.with(|i| {
            let s = i.trin.as_mut().unwrap();
            s.ntr = end / s.nsegy as i64;
        });
    }
    let (ntr, nsegy, nsfirst) = fp.with(|i| {
        let s = i.trin.as_ref().unwrap();
        (s.ntr, s.nsegy, s.nsfirst)
    });
    if itr as i64 >= ntr {
        crate::err!("%s: trying to read off end of file", "fgettr.c");
    }
    if 0 > efseeko(fp, itr as i64 * nsegy as i64, SEEK_SET) {
        crate::err!("%s: unable to seek xdr disk file to trace %d", "fgettr.c", itr);
    }
    let nread = fgettr(fp, tp);
    if nread as usize != nsegy {
        crate::err!("%s: read %d bytes with %d bytes in trace", "fgettr.c", nread, nsegy as i32);
    }
    if tp.hdr.ns as i32 != nsfirst {
        crate::warn!("%s: header ns field = %d differs from first trace = %d", "fgettr.c", tp.hdr.ns, nsfirst);
    }
    ntr as i32
}

/// gettra: fgettra(stdin).
pub fn gettra(tp: &mut Trace, itr: i32) -> i32 {
    fgettra(&stdin(), tp, itr)
}

fn datawrite(fp: &CFile, st: &TrOut, tp: &mut Trace, fixed_length: bool) -> i32 {
    let nstobewritten: u32 = if fixed_length { st.nsfirst } else { tp.hdr.ns as u32 };
    let databytes = st.bytesper as usize * nstobewritten as usize;
    tp.ensure_bytes(databytes);
    let mut nwritten: i32;
    match tp.hdr.trid {
        CHARPACK | SHORTPACK => {
            if tp.hdr.trid == CHARPACK {
                // C falls through from CHARPACK into SHORTPACK: written twice.
                efwrite(&tp.data_bytes()[..databytes], fp);
            }
            swab(tp.data_bytes_mut(), databytes);
            nwritten = efwrite(&tp.data_bytes()[..databytes], fp) as i32;
        }
        _ => {
            let mut buf = Vec::with_capacity(XDR_UNIT as usize * nstobewritten as usize);
            for k in 0..nstobewritten as usize {
                buf.extend_from_slice(&tp.data[k].to_bits().to_be_bytes());
            }
            // C encodes all samples, then writes the byte count fixed by the
            // first trace, even when its packed trid selected fewer bytes.
            nwritten = efwrite(&buf[..databytes], fp) as i32;
            if nwritten as usize != databytes {
                nwritten = 0;
            }
        }
    }
    nwritten
}

fn fputtr_internal(fp: &CFile, tp: &mut Trace, fixed_length: bool) {
    let st = fp.with(|i| i.trout.clone());
    let mut st = match st {
        Some(s) => s,
        None => {
            match filestat(fileno(fp)) {
                FileType::Directory => crate::err!("%s: segy output can't be a directory", "fputtr.c"),
                FileType::Tty => crate::err!("%s: segy output can't be tty", "fputtr.c"),
                _ => {}
            }
            let nsfirst = tp.hdr.ns as u32;
            if nsfirst > SU_NFLTS as u32 {
                crate::err!("%s: unable to handle %d > %d samples per trace", "fputtr.c", nsfirst, SU_NFLTS as i32);
            }
            let s = TrOut { itr: 0, nsfirst, bytesper: bytes_per(tp.hdr.trid) };
            fp.with(|i| i.trout = Some(s.clone()));
            s
        }
    };
    let hb = tp.hdr.to_xdr();
    let nwritten = efwrite(&hb[..], fp);
    if nwritten != HDRBYTES {
        crate::err!("%s: unable to write header on trace #%ld", "fputtr.c", st.itr + 1);
    }
    let databytes = st.bytesper as u32 * if fixed_length { st.nsfirst } else { tp.hdr.ns as u32 };
    let nwritten = datawrite(fp, &st, tp, fixed_length);
    if nwritten as u32 != databytes {
        crate::err!(
            "%s: on trace #%ld, tried to write %d bytes, wrote %d bytes",
            "fputtr.c",
            st.itr + 1,
            databytes,
            nwritten
        );
    }
    st.itr += 1;
    fp.with(|i| i.trout = Some(st));
}

/// fputtr: write a trace (header + nsfirst samples) to a stream.
pub fn fputtr(fp: &CFile, tp: &mut Trace) {
    fputtr_internal(fp, tp, true)
}
/// fvputtr: write a trace with its own ns.
pub fn fvputtr(fp: &CFile, tp: &mut Trace) {
    fputtr_internal(fp, tp, false)
}
/// puttr: fputtr(stdout).
pub fn puttr(tp: &mut Trace) {
    fputtr(&stdout(), tp)
}
/// vputtr: fvputtr(stdout).
pub fn vputtr(tp: &mut Trace) {
    fvputtr(&stdout(), tp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cio::{efclose, efopen, erewind, etmpfile};

    #[test]
    fn reopened_file_starts_with_fresh_input_state() {
        // Keep the closed handles alive: a new CFile, even for the same path,
        // must not inherit nsfirst, itr or the disk trace count from them.
        let path = std::env::temp_dir().join(format!("oxwave_trace_reopen_{}", std::process::id()));
        let name = path.to_str().unwrap();
        let mut tr = Trace::new();
        tr.ns = 3;
        tr.data[..4].copy_from_slice(&[1.0, 2.0, 3.0, 4.0]);
        let out = efopen(name, "w");
        fputtr(&out, &mut tr);
        efclose(&out);
        let first = efopen(name, "r");
        let mut actual = Trace::new();
        assert_eq!(fgettr(&first, &mut actual), 252);
        assert_eq!(fgettr(&first, &mut actual), 0);
        efclose(&first);

        tr.ns = 4;
        let out = efopen(name, "w");
        fputtr(&out, &mut tr);
        fputtr(&out, &mut tr);
        efclose(&out);
        let second = efopen(name, "r");
        assert_eq!(fgettra(&second, &mut actual, 1), 2);
        assert_eq!(actual.ns, 4);
        assert_eq!(&actual.data[..4], &[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(second.with(|i| i.trin.as_ref().unwrap().itr), 2);
        efclose(&second);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn roundtrip_tmpfile() {
        let fp = etmpfile();
        let mut tr = Trace::new();
        tr.ns = 4;
        tr.tracl = 1;
        tr.offset = -400;
        tr.data[..4].copy_from_slice(&[1.0, -2.5, f32::MIN_POSITIVE, 3e38]);
        fputtr(&fp, &mut tr);
        tr.tracl = 2;
        tr.data[0] = 7.0;
        fputtr(&fp, &mut tr);
        erewind(&fp);
        let mut raw = vec![0u8; 2 * (240 + 16)];
        assert_eq!(efread(&mut raw[..], &fp), raw.len());
        assert_eq!(&raw[0..4], &[0, 0, 0, 1]);
        assert_eq!(&raw[36..40], &(-400i32).to_be_bytes());
        assert_eq!(&raw[114..116], &[0, 4]);
        assert_eq!(&raw[244..248], &(-2.5f32).to_be_bytes());
        erewind(&fp);
        let mut t2 = Trace::new();
        assert_eq!(fgettr(&fp, &mut t2), 256);
        assert_eq!(t2.tracl, 1);
        assert_eq!(&t2.data[..4], &[1.0, -2.5, f32::MIN_POSITIVE, 3e38]);
        assert_eq!(fgettr(&fp, &mut t2), 256);
        assert_eq!(t2.data[0], 7.0);
        assert_eq!(fgettr(&fp, &mut t2), 0);
        assert_eq!(fgettra(&fp, &mut t2, 0), 2);
        assert_eq!(t2.tracl, 1);
    }

    #[test]
    fn short_data_keeps_previous_samples() {
        let fp = etmpfile();
        let mut tr = Trace::new();
        tr.ns = 2;
        tr.data[0] = 5.0;
        tr.data[1] = 6.0;
        fputtr(&fp, &mut tr);
        let mut h = tr.hdr;
        h.tracl = 9;
        let hb = h.to_xdr();
        efwrite(&hb[..], &fp);
        efwrite(&[1u8, 2, 3][..], &fp);
        erewind(&fp);
        let mut t = Trace::new();
        assert_eq!(fgettr(&fp, &mut t), 248);
        assert_eq!(fgettr(&fp, &mut t), 240);
        assert_eq!(t.tracl, 9);
        assert_eq!(&t.data[..2], &[5.0, 6.0]);
        assert_eq!(fgettr(&fp, &mut t), 0);
    }

    fn mixed_output_prefix(first_trid: i16, fixed_length: bool) {
        let fp = etmpfile();
        let mut tr = Trace::new();
        tr.ns = 7;
        tr.trid = first_trid;
        tr.data_bytes_mut()[..14].copy_from_slice(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]);
        let mut expected = tr.hdr.to_xdr().to_vec();
        match first_trid {
            CHARPACK => expected.extend_from_slice(&[1, 2, 3, 4, 5, 6, 7, 2, 1, 4, 3, 6, 5, 7]),
            SHORTPACK => expected.extend_from_slice(&[2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13]),
            _ => unreachable!(),
        }
        if fixed_length {
            fputtr(&fp, &mut tr);
        } else {
            fvputtr(&fp, &mut tr);
        }

        tr.ns = 5;
        tr.trid = 1;
        tr.data[..7].copy_from_slice(&[1.0, -2.5, 3.0, 4.0, 5.0, 6.0, 7.0]);
        expected.extend_from_slice(&tr.hdr.to_xdr());
        let xdr = [
            0x3f, 0x80, 0, 0, 0xc0, 0x20, 0, 0, 0x40, 0x40, 0, 0, 0x40, 0x80, 0, 0,
            0x40, 0xa0, 0, 0, 0x40, 0xc0, 0, 0, 0x40, 0xe0, 0, 0,
        ];
        let prefix_len = match (first_trid, fixed_length) {
            (CHARPACK, true) => 7,
            (CHARPACK, false) => 5,
            (SHORTPACK, true) => 14,
            (SHORTPACK, false) => 10,
            _ => unreachable!(),
        };
        expected.extend_from_slice(&xdr[..prefix_len]);
        if fixed_length {
            fputtr(&fp, &mut tr);
        } else {
            fvputtr(&fp, &mut tr);
        }
        assert_eq!(&tr.data[..7], &[1.0, -2.5, 3.0, 4.0, 5.0, 6.0, 7.0]);
        assert_eq!(eftello(&fp), expected.len() as i64);
        erewind(&fp);
        let mut actual = vec![0; expected.len() + 1];
        assert_eq!(efread(&mut actual, &fp), expected.len());
        assert_eq!(&actual[..expected.len()], &expected);
    }

    #[test]
    fn fixed_output_float_after_shortpack_writes_xdr_prefix() {
        mixed_output_prefix(SHORTPACK, true);
    }

    #[test]
    fn fixed_output_float_after_charpack_writes_xdr_prefix() {
        mixed_output_prefix(CHARPACK, true);
    }

    #[test]
    fn variable_output_float_after_packed_uses_current_ns() {
        mixed_output_prefix(CHARPACK, false);
        mixed_output_prefix(SHORTPACK, false);
    }
}
