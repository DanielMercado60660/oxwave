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

//! Gather I/O, ported from `su/lib/fgetgthr.c` and `su/lib/fputgthr.c`
//! (Revision 1.3, 2011/11/11).
//! Credits: Potash Corporation Sascatchewan, Balasz Nemeth; given to CWP in 2008.
//!
//! Gather keys must be contiguous in the input. `nt` and `dt` are fixed by
//! the first trace, even when later headers differ. As in C, `n_val` on
//! return describes the lookahead trace, not necessarily the returned gather.
//! The C reader's EOF quirk is preserved: a final singleton after an earlier
//! gather is discarded. Read until `None`; do not call again without resetting
//! `first` to zero. Temporary files close automatically, including on EOF.

use oxwave_core::prelude::*;

/// The two static lookahead files in each C reader, owned by its caller here.
/// Use a separate state for each input stream and initialize `first` to zero.
#[derive(Default, Debug)]
pub struct GatherState {
    ntracefp: Option<CFile>,
    nheaderfp: Option<CFile>,
}

/// Read one gather from `fp`, returning `None` at the C reader's EOF condition.
///
/// C prototype:
/// `segy **fget_gather(FILE *fp, cwp_String *key, cwp_String *type, Value *n_val, int *nt, int *ntr, float *dt, int *first);`
///
/// `state` replaces C's static files. `key` is borrowed (C never changes it);
/// `type_` receives the header type. Returned traces own exactly `nt` samples.
/// Later shorter traces retain the preceding trace's samples up to `nt`,
/// and longer traces are truncated to `nt`, exactly like C's temporary files.
#[allow(clippy::too_many_arguments)]
pub fn fget_gather(
    fp: &CFile,
    state: &mut GatherState,
    key: &str,
    type_: &mut &'static str,
    n_val: &mut Value,
    nt: &mut i32,
    ntr: &mut i32,
    dt: &mut f32,
    first: &mut i32,
) -> Option<Vec<Trace>> {
    let mut tr = Trace::new();
    let mut ntrr = 0i32;
    *type_ = hdtype(key);
    let indx = getindex(key);
    *ntr = 0;

    if *first == 0 {
        let nsegy = fvgettr(fp, &mut tr);
        if nsegy == 0 {
            err!("can't get first trace");
        }
        *nt = tr.ns as i32;
        *dt = ((tr.dt as f32 as f64) / 1000000.0) as f32;
        ntrr = ntrr.wrapping_add(1);
        gethval(&tr, indx, n_val);
        state.ntracefp = Some(etmpfile());
        state.nheaderfp = Some(etmpfile());
        *first = 1;
    } else {
        let nheaderfp = state
            .nheaderfp
            .as_ref()
            .expect("gather state is not initialized");
        let ntracefp = state
            .ntracefp
            .as_ref()
            .expect("gather state is not initialized");
        erewind(nheaderfp);
        erewind(ntracefp);
        fread(std::slice::from_mut(&mut tr.hdr), nheaderfp);
        fread(&mut tr.data[..*nt as usize], ntracefp);
        gethval(&tr, indx, n_val);
    }

    let tracefp = etmpfile();
    let headerfp = etmpfile();
    let nsegy = loop {
        efwrite(std::slice::from_ref(&tr.hdr), &headerfp);
        efwrite(&tr.data[..*nt as usize], &tracefp);
        *ntr = ntr.wrapping_add(1);
        let val = *n_val;
        let nsegy = fvgettr(fp, &mut tr);
        if nsegy != 0 {
            ntrr = ntrr.wrapping_add(1);
        }
        gethval(&tr, indx, n_val);
        if nsegy == 0 || valcmp(type_, val, *n_val) != 0 {
            break nsegy;
        }
    };

    if nsegy == 0 && ntrr == 0 {
        *ntr = 0;
        state.nheaderfp = None;
        state.ntracefp = None;
        return None;
    } else {
        let nheaderfp = state.nheaderfp.as_ref().unwrap();
        let ntracefp = state.ntracefp.as_ref().unwrap();
        erewind(nheaderfp);
        erewind(ntracefp);
        efwrite(std::slice::from_ref(&tr.hdr), nheaderfp);
        efwrite(&tr.data[..*nt as usize], ntracefp);
    }

    let mut rec = Vec::with_capacity(*ntr as usize);
    for _ in 0..*ntr {
        rec.push(Trace {
            hdr: SegyHeader::zeroed(),
            data: vec![0.0; *nt as usize],
        });
    }
    erewind(&headerfp);
    erewind(&tracefp);
    for trace in &mut rec {
        fread(std::slice::from_mut(&mut trace.hdr), &headerfp);
    }
    efclose(&headerfp);
    for trace in &mut rec {
        fread(&mut trace.data, &tracefp);
    }
    efclose(&tracefp);
    Some(rec)
}

/// Read one gather from standard input, with the same semantics as `fget_gather`.
///
/// C prototype:
/// `segy **get_gather(cwp_String *key, cwp_String *type, Value *n_val, int *nt, int *ntr, float *dt, int *first);`
///
/// `state` replaces this function's independent C static lookahead files.
#[allow(clippy::too_many_arguments)]
pub fn get_gather(
    state: &mut GatherState,
    key: &str,
    type_: &mut &'static str,
    n_val: &mut Value,
    nt: &mut i32,
    ntr: &mut i32,
    dt: &mut f32,
    first: &mut i32,
) -> Option<Vec<Trace>> {
    fget_gather(&stdin(), state, key, type_, n_val, nt, ntr, dt, first)
}

/// Write all traces in a gather using variable-length trace output.
///
/// C prototype:
/// `segy **fput_gather(FILE *fp, segy **rec, int *nt, int *ntr);`
///
/// The slice supplies `ntr`; `nt` is the number of samples copied per trace,
/// while each header's `ns` determines the samples written. Each trace must
/// contain at least `nt` samples and have `ns <= nt`, as required to avoid
/// uninitialized samples in C. Rust borrows the gather instead of freeing its
/// rows; the caller controls its lifetime. C always returns NULL, mapped to `()`.
pub fn fput_gather(fp: &CFile, rec: &[Trace], nt: i32) {
    let mut tr = Trace::new();
    for trace in rec {
        tr.hdr = trace.hdr;
        tr.data[..nt as usize].copy_from_slice(&trace.data[..nt as usize]);
        fvputtr(fp, &mut tr);
    }
}

/// Write a gather to standard output, borrowing its traces like `fput_gather`.
///
/// C prototype: `segy **put_gather(segy **rec, int *nt, int *ntr);`
/// The slice supplies `ntr`; C's unconditional NULL return is mapped to `()`.
pub fn put_gather(rec: &[Trace], nt: i32) {
    fput_gather(&stdout(), rec, nt)
}
