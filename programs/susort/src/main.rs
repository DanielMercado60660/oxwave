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

//! SUSORT - sort on any segy header keywords.
//! Port of su/main/windowing_sorting_muting/susort.c (Revision 1.28).
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits (from the C source):
//! SEP: Einar Kjartansson, Stew Levin.
//! CWP: Shuki Ronen, Jack K. Cohen.
//! Earlier SEP tape support was removed in version 1.16; version 1.15 is
//! in the Portability directory for users of SU data stored on tape.
//! Trace header fields modified: tracl, tracr.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;
use std::cell::Cell;
use std::ffi::{c_char, c_void, CString};
use std::sync::atomic::{AtomicI32, AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t",
    " SUSORT - sort on any segy header keywords\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " susort <stdin >stdout [[+-]key1 [+-]key2 ...]\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Susort supports any number of (secondary) keys with either\t",
    " ascending (+, the default) or descending (-) directions for \t",
    " each.  The default sort key is cdp.\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Note:\tOnly the following types of input/output are supported\t",
    "\tDisk input --> any output\t\t\t\t",
    "\tPipe input --> Disk output\t\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Caveat:  On some Linux systems Pipe input and or output often ",
    "\t\tfails\t\t\t\t\t\t",
    "\tDisk input ---> Disk output is recommended\t\t",
    " \t\t\t\t\t\t\t\t",
    " Note: If the the CWP_TMPDIR environment variable is set use\t",
    "\tits value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t",
    " Example:\t\t\t\t\t\t\t",
    " To sort traces by cdp gather and within each gather\t\t",
    " by offset with both sorts in ascending order:\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " \tsusort <INDATA >OUTDATA cdp offset\t\t\t",
    " \t\t\t\t\t\t\t\t",
    " Caveat: In the case of Pipe input a temporary file is made\t",
    "\tto hold the ENTIRE data set.  This temporary is\t\t",
    "\teither an actual disk file (usually in /tmp) or in some\t",
    "\timplementations, a memory buffer.  It is left to the\t",
    "\tuser to be SENSIBLE about how big a file to pipe into\t",
    "\tsusort relative to the user's computer.\t\t\t",
    " \t\t\t\t\t\t\t\t",
];

const NTRSTEP: usize = 1024;

// Call the same platform library routines as C. In particular, qsort's
// permutation of equal keys is observable in the resulting trace order.
unsafe extern "C" {
    fn qsort(
        base: *mut c_void,
        n: usize,
        size: usize,
        compare: unsafe extern "C" fn(*const c_void, *const c_void) -> i32,
    );
    fn access(path: *const c_char, mode: i32) -> i32;
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
    fn close(fd: i32) -> i32;
    fn unlink(path: *const c_char) -> i32;
    fn _exit(status: i32) -> !;
}

thread_local! {
    // The C comparator's globals refer to main's immutable type list only
    // during the synchronous qsort call.
    static TYPES: Cell<*const Vec<&'static str>> = const { Cell::new(std::ptr::null()) };
}

// Signal cleanup needs only the descriptor and immutable filename. Avoid
// borrowing a CFile that might already be borrowed by interrupted trace I/O.
static TRACE_FD: AtomicI32 = AtomicI32::new(-1);
static TRACEFILE: AtomicPtr<c_char> = AtomicPtr::new(std::ptr::null_mut());

fn main() {
    let argv = args();
    let argc = argv.len();
    let mut tr = Trace::new();

    /* Initialize */
    initargs(&argv);
    requestdoc(1, SDOC);

    /* Look for user-supplied tmpdir from environment */
    let tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    if !tmpdir.is_empty() {
        let directory = CString::new(tmpdir.as_bytes()).unwrap();
        // SAFETY: directory is NUL terminated; WRITE_OK == W_OK == 2.
        if unsafe { access(directory.as_ptr(), 2) } != 0 {
            err!("you can't write in %s (or it doesn't exist)", tmpdir);
        }
    }

    /* The algorithm requires stdin to be rewound and random access to
     * either stdin or stdout. */
    let ftypein = filestat(STDIN);
    let ftypeout = filestat(STDOUT);
    let isdisk = ftypein == FileType::Disk;
    let ispipe = ftypein == FileType::Pipe || ftypein == FileType::Fifo;
    if ispipe && ftypeout != FileType::Disk {
        err!("OK IO types are: DISK->ANYTHING, PIPE->DISK");
    }
    /* If pipe, prepare temporary file to hold data. */
    let mut tracefile = String::new();
    let tracefile_c;
    let tracefp = if ispipe {
        if tmpdir.is_empty() {
            Some(etmpfile())
        } else {
            tracefile = temporary_filename(&tmpdir);
            tracefile_c = CString::new(tracefile.as_bytes()).unwrap();
            TRACEFILE.store(tracefile_c.as_ptr().cast_mut(), Ordering::Relaxed);
            // SAFETY: tracefile_c stays alive until exit; closefiles has the
            // signal-handler ABI. These are SIGINT, SIGQUIT, SIGHUP, SIGTERM.
            unsafe {
                for sig in [2, 3, 1, 15] {
                    signal(sig, closefiles);
                }
            }
            let fp = efopen(&tracefile, "w+");
            TRACE_FD.store(fileno(&fp), Ordering::Relaxed);
            Some(fp)
        }
    } else {
        None
    };
    /* Set number of sort keys and allocate lists. */
    let nkey = if argc == 1 { 1 } else { argc - 1 };
    let mut index = vec![0; nkey];
    let mut up = vec![false; nkey];
    let mut types = vec![""; nkey];

    /* Initialize index, type and up. */
    if argc == 1 {
        index[0] = getindex("cdp");
        up[0] = true;
        types[0] = "l";
    } else {
        for i in 0..nkey {
            let arg = argv[i + 1].as_str();
            let key;
            if let Some(k) = arg.strip_prefix('+') {
                up[i] = true;
                key = k;
            } else if let Some(k) = arg.strip_prefix('-') {
                up[i] = false;
                key = k;
            } else {
                up[i] = true;
                key = arg;
            }
            index[i] = getindex(key);
            types[i] = hdtype(key);
        }
    }

    /* Allocate list of trace numbers + key values. Keep the C record width
     * and layout because the system qsort may use them when choosing swaps. */
    let ngroup = nkey + 1;
    let nvsize = ngroup * std::mem::size_of::<Value>();
    let mut nv = NTRSTEP;
    let mut val_list = vec![Value::default(); nv * ngroup];

    /* Run through traces once to collect header values. */
    let mut ntr = 0i64;
    let nsegy = gettr(&mut tr);
    if nsegy == 0 {
        err!("can't get first trace");
    }
    // The C Value is reused across keys and traces; narrow header writes
    // preserve the other union bytes. In particular the default uses type l.
    let mut val = Value::default();
    loop {
        let itr = ntr;
        ntr += 1;
        if ntr % NTRSTEP as i64 == 0 {
            nv += NTRSTEP;
            val_list.resize(nv * ngroup, Value::default());
        }

        /* Enter trace index in list and then key values. */
        let mut vptr = itr as usize * ngroup;
        val_list[vptr].set_l(itr);
        vptr += 1;
        for i in 0..nkey {
            gethval(&tr, index[i], &mut val);
            val_list[vptr] = if up[i] { val } else { negval(types[i], val) };
            vptr += 1;
        }
        if let Some(fp) = &tracefp {
            // C writes the in-memory segy, including native-endian samples.
            efwrite(tr.hdr.as_bytes(), fp);
            efwrite(&tr.data_bytes()[..nsegy as usize - HDRBYTES], fp);
        }
        if gettr(&mut tr) == 0 {
            break;
        }
    }
    if let Some(fp) = &tracefp {
        erewind(fp);
    } else {
        erewind(&stdin());
    }

    /* Sort the header values. */
    TYPES.with(|ptr| ptr.set(&types));
    // SAFETY: there are ntr initialized records of nvsize bytes. Value is
    // Copy, all its bit patterns are valid, and cmp_list only reads records.
    // The immutable types list remains alive and in place during qsort.
    unsafe { qsort(val_list.as_mut_ptr().cast(), ntr as usize, nvsize, cmp_list) };
    TYPES.with(|ptr| ptr.set(std::ptr::null()));

    if isdisk {
        /* Run through sorted list and write output sequentially. */
        for i in 0..ntr as usize {
            let itr = val_list[i * ngroup].l();
            gettra(&mut tr, itr as i32);
            tr.tracl = (i + 1) as i32;
            tr.tracr = tr.tracl;
            puttr(&mut tr);
        }
    } else if let Some(fp) = &tracefp {
        /* Invert permutation and read input sequentially. */
        for i in 0..ntr as usize {
            let itr = val_list[i * ngroup].l();
            val_list[itr as usize * ngroup + 1].set_l(i as i64);
        }
        for i in 0..ntr as usize {
            let itr = val_list[i * ngroup + 1].l();
            efread(tr.hdr.as_bytes_mut(), fp);
            efread(&mut tr.data_bytes_mut()[..nsegy as usize - HDRBYTES], fp);
            tr.tracl = (itr + 1) as i32;
            tr.tracr = tr.tracl;
            efseeko(&stdout(), itr * nsegy as i64, SEEK_SET);
            efwrite(tr.hdr.as_bytes(), &stdout());
            efwrite(&tr.data_bytes()[..nsegy as usize - HDRBYTES], &stdout());
        }
    }

    /* Clean up. */
    if let Some(fp) = &tracefp {
        efclose(fp);
        TRACE_FD.store(-1, Ordering::Relaxed);
        if !tracefile.is_empty() {
            eremove(&tracefile);
        }
    }
    exit(CWP_Exit());
}

/* Comparison routine for qsort. */
unsafe extern "C" fn cmp_list(a: *const c_void, b: *const c_void) -> i32 {
    TYPES.with(|ptr| {
        // SAFETY: main publishes types for the duration of qsort, which
        // supplies pointers to initialized groups of nkey + 1 Values.
        let types = unsafe { &*ptr.get() };
        let pa = a.cast::<Value>();
        let pb = b.cast::<Value>();
        for i in 0..types.len() {
            let va = unsafe { *pa.add(i + 1) };
            let vb = unsafe { *pb.add(i + 1) };
            let compare = valcmp(types[i], va, vb);
            if compare != 0 {
                return compare;
            }
        }
        0
    })
}

/* Reverse sign of value. */
fn negval(t: &str, mut val: Value) -> Value {
    match t.as_bytes()[0] {
        b'h' => val.set_h(val.h().wrapping_neg()),
        b'u' => val.set_u(u16::MAX - val.u()),
        b'l' => val.set_l(val.l().wrapping_neg()),
        b'v' => val.set_v(u64::MAX - val.v()),
        b'i' => val.set_i(val.i().wrapping_neg()),
        b'p' => val.set_p(val.p().wrapping_neg()),
        b'f' => val.set_f(-val.f()),
        b'd' => val.set_d(-val.d()),
        _ => err!("%d: mysterious type %s", 304, t),
    }
    val
}

/* For graceful interrupt termination. */
extern "C" fn closefiles(_sig: i32) {
    // SAFETY: main publishes the open descriptor and a NUL-terminated name
    // that stays alive until exit. Atomic loads are lock-free; close, unlink,
    // and _exit are async-signal-safe. Do not use the core exit/flush routine:
    // a signal can interrupt its thread-local state, allocations, or streams.
    unsafe {
        close(TRACE_FD.load(Ordering::Relaxed));
        unlink(TRACEFILE.load(Ordering::Relaxed));
        _exit(EXIT_FAILURE);
    }
}
