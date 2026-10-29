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

//! SUQUANTILE - display quantiles or ranks of a data set.
//! Port of su/main/attributes_parameter_estimation/suquantile.c (Revision 1.9).
//! Credits: CWP: Jack K. Cohen.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;
use std::cell::Cell;
use std::ffi::{c_char, c_void, CString};
use std::sync::atomic::{AtomicI32, AtomicPtr, Ordering};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUQUANTILE - display some quantiles or ranks of a data set            ",
    "                                                                       ",
    " suquantile <stdin >stdout [optional parameters]\t\t\t",
    "                                                                       ",
    " Required parameters:                                                  ",
    "       none (no-op)                                                    ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "\tpanel=1\t\tflag; 0 = do trace by trace (vs. whole data set)",
    "\tquantiles=1\tflag; 0 = give ranks instead of quantiles\t",
    " \tverbose=0\tverbose = 1 echoes information\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " \ttmpdir= \t if non-empty, use the value as a directory path",
    "\t\t\t prefix for storing temporary files; else if the",
    "\t\t         the CWP_TMPDIR environment variable is set use\t",
    "\t\t         its value for the path; else use tmpfile()\t",
    " \t\t\t\t\t\t\t\t\t",
];

// These are platform libc/OS functions, not SU library routines. In particular,
// qsort's permutation of equal values determines the reported sample indices.
unsafe extern "C" {
    fn qsort(
        base: *mut c_void,
        count: usize,
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
    // The original global data pointer is only needed during synchronous qsort.
    static DATA: Cell<*const f32> = const { Cell::new(std::ptr::null()) };
}

// Comparison function for qsort. Keep subtraction in float, including NaNs.
unsafe extern "C" fn cmp_indirect(r: *const c_void, s: *const c_void) -> i32 {
    DATA.with(|data| {
        // SAFETY: qsort passes pointers into rank; its indices reference data,
        // which remains alive and immutable throughout the sort.
        let diff = unsafe {
            *data.get().add(*r.cast::<i32>() as usize) - *data.get().add(*s.cast::<i32>() as usize)
        };
        if diff > 0.0 {
            1
        } else if diff < 0.0 {
            -1
        } else {
            0
        }
    })
}

fn sort_ranks(rank: &mut [i32], data: &[f32]) {
    for i in 0..data.len() {
        rank[i] = i as i32;
    }
    DATA.with(|ptr| ptr.set(data.as_ptr()));
    // SAFETY: rank contains count initialized i32s. FSIZE is also sizeof(i32),
    // as assumed by the original qsort call; the callback cannot unwind.
    unsafe {
        qsort(
            rank.as_mut_ptr().cast(),
            data.len(),
            std::mem::size_of::<f32>(),
            cmp_indirect,
        );
    }
    DATA.with(|ptr| ptr.set(std::ptr::null()));
}

// Globals used by the named temporary file's interrupt cleanup. Publish the
// immutable, NUL-terminated filename before installing the handlers. Use OS
// close/unlink/_exit and lock-free atomics only: the handler must not allocate,
// borrow thread-local state, flush streams, or run Rust/libc exit callbacks.
static TRACE_FD: AtomicI32 = AtomicI32::new(-1);
static TRACEFILE: AtomicPtr<c_char> = AtomicPtr::new(std::ptr::null_mut());

extern "C" fn closefiles(_sig: i32) {
    let fd = TRACE_FD.load(Ordering::Relaxed);
    let path = TRACEFILE.load(Ordering::Relaxed);
    // SAFETY: the CString lives until process exit. close tolerates an unopened
    // or already closed descriptor. These POSIX calls are async-signal-safe.
    // Reports go to unbuffered stderr and stdout is unused; only the discarded
    // temporary file can have buffered writes. Normal termination uses exit().
    unsafe {
        close(fd);
        unlink(path);
        _exit(EXIT_FAILURE);
    }
}

// GET_RANKS - fprint some sample ranks.
fn get_ranks(rank: &mut [i32], data: &[f32]) {
    let nt = data.len();
    sort_ranks(rank, data);

    // Could make table of desired i's and loop.
    let mut i = 0;
    fprintf!(
        stderr(),
        " rank[%d] = %8.2e",
        (i + 1) as i32,
        data[rank[i] as usize]
    );
    i = nt / 20;
    fprintf!(
        stderr(),
        " rank[%d] = %8.2e",
        (i + 1) as i32,
        data[rank[i] as usize]
    );
    i = nt / 2 - i;
    fprintf!(
        stderr(),
        " rank[%d] = %8.2e",
        (i + 1) as i32,
        data[rank[i] as usize]
    );
    fprintf!(stderr(), "\n");
    i = nt - 1 - i;
    fprintf!(
        stderr(),
        " rank[%d] = %8.2e",
        (i + 1) as i32,
        data[rank[i] as usize]
    );
    i = nt - 1 - nt / 20;
    fprintf!(
        stderr(),
        " rank[%d] = %8.2e",
        (i + 1) as i32,
        data[rank[i] as usize]
    );
    i = nt - 1;
    fprintf!(
        stderr(),
        " rank[%d] = %8.2e",
        (i + 1) as i32,
        data[rank[i] as usize]
    );
    fprintf!(
        stderr(),
        "\nmin is at sample %d,  max at %d\n",
        rank[0] + 1,
        rank[nt - 1] + 1
    );
}

// GET_QUANTILES - fprint some sample quantiles.
fn get_quantiles(rank: &mut [i32], data: &[f32]) {
    let nt = data.len();
    sort_ranks(rank, data);

    // Round to qth quantile (100 qth percentile):
    // thus (q*nt - 1) + .5 (-1 for zero basing).
    // Cast through i32 to preserve C's truncation toward zero for small nt.
    let mut i = 1;
    let mut iq = (0.01 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        " %dst percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    i = 5;
    iq = (0.05 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        " %dth percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    i = 25;
    iq = (0.25 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        "%dth percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    i = 50;
    iq = (0.50 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        "%dth percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    i = 75;
    iq = (0.75 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        "%dth percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    i = 95;
    iq = (0.95 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        "%dth percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    i = 99;
    iq = (0.99 * nt as f64 - 0.5) as i32;
    fprintf!(
        stderr(),
        "%dth percentile is %8.2e\n",
        i + 1,
        data[rank[iq as usize] as usize]
    );
    fprintf!(
        stderr(),
        "min at sample %d equals %8.2e\n",
        rank[0] + 1,
        data[rank[0] as usize]
    );
    fprintf!(
        stderr(),
        "max at sample %d equals %8.2e\n",
        rank[nt - 1] + 1,
        data[rank[nt - 1] as usize]
    );
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();

    // Initialize.
    initargs(&argv);
    requestdoc(1, SDOC);

    // Get nt from first trace.
    if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let nt = tr.ns as usize;

    // Get parameters.
    let mut panel = 1;
    let mut quantiles = 1;
    let mut verbose = 0;
    getparint("panel", &mut panel);
    getparint("quantiles", &mut quantiles);
    getparint("verbose", &mut verbose);

    // Look for user-supplied tmpdir.
    let mut tmpdir = String::new();
    if getparstring("tmpdir", &mut tmpdir) == 0 {
        tmpdir = std::env::var("CWP_TMPDIR").unwrap_or_default();
    }
    let tmpdir_c = CString::new(tmpdir.as_str()).expect("getpar strings contain no NUL");
    // WRITE_OK is 2 on the target POSIX platform. Match access() directly,
    // including writable regular files and directories with ACLs.
    if !tmpdir.is_empty() && unsafe { access(tmpdir_c.as_ptr(), 2) } != 0 {
        err!("you can't write in %s (or it doesn't exist)", tmpdir);
    }
    checkpars();

    // Main loop over traces.
    if panel == 0 {
        let mut rank = vec![0i32; nt];
        let mut data = vec![0.0f32; nt];
        fprintf!(stderr(), "trace length: %d\n", nt as i32);
        loop {
            fprintf!(stderr(), "trace: %d\n", tr.tracl);
            data.copy_from_slice(&tr.data[..nt]);
            if tr.mark == 0 {
                if quantiles == 0 {
                    get_ranks(&mut rank, &data);
                } else {
                    get_quantiles(&mut rank, &data);
                }
                tr.data[..nt].copy_from_slice(&data);
            }
            if gettr(&mut tr) == 0 {
                break;
            }
        }
    } else {
        let (tracefp, tracefile) = if tmpdir.is_empty() {
            let fp = etmpfile();
            if verbose != 0 {
                warn!("using tmpfile() call");
            }
            (fp, String::new())
        } else {
            let name = temporary_filename(&tmpdir);
            let filename = CString::new(name.as_str()).expect("temporary filename contains no NUL");
            // Handlers and their immutable filename live until exit(), like the
            // original globals. The temporary filename library owns creation.
            TRACEFILE.store(filename.into_raw(), Ordering::Relaxed);
            for sig in [2, 3, 1, 15] {
                unsafe {
                    signal(sig, closefiles);
                }
            }
            let fp = efopen(&name, "w+");
            TRACE_FD.store(fileno(&fp), Ordering::Relaxed);
            if verbose != 0 {
                warn!("putting temporary files in %s", tmpdir);
            }
            (fp, name)
        };
        let mut ntr = 0i32;
        loop {
            ntr = ntr.wrapping_add(1);
            efwrite(&tr.data[..nt], &tracefp);
            if gettr(&mut tr) == 0 {
                break;
            }
        }
        erewind(&tracefp);
        let ndata = (nt as i32).wrapping_mul(ntr);
        fprintf!(stderr(), "trace length: %d\n", nt as i32);
        fprintf!(stderr(), "number of traces: %d\n", ntr);
        fprintf!(stderr(), "number of samples: %d\n", ndata);
        let mut rank = vec![0i32; ndata as usize];
        let mut data = vec![0.0f32; ndata as usize];
        efread(&mut data, &tracefp);
        efclose(&tracefp);
        if quantiles == 0 {
            get_ranks(&mut rank, &data);
        } else {
            get_quantiles(&mut rank, &data);
        }
        if !tracefile.is_empty() {
            eremove(&tracefile);
        }
    }
    exit(CWP_Exit());
}
