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

//! SUEDIT - examine segy diskfiles and edit headers.
//! Port of su/main/headers/suedit.c (revision 1.38).
//! Credits: SEP: Einar Kjartansson, Shuki Ronen, Stew Levin;
//! CWP: Jack K. Cohen; Unocal: Reg Beardsley.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;
use oxwave_trace_display::tabplot;
use std::ffi::CString;
use std::os::raw::{c_char,c_int};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUEDIT - examine segy diskfiles and edit headers\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suedit diskfile  (open for possible header modification if writable)\t",
    " suedit <diskfile  (open read only)\t\t\t\t\t",
    " \t\t\t\t\t\t\t        \t",
    " The following commands are recognized:\t\t\t\t",
    " number\tread in that trace and print nonzero header words\t",
    " <CR>\t\tgo to trace one step away (step is initially -1)\t",
    " +\t\tread in next trace (step is set to +1)\t\t\t",
    " -\t\tread in previous trace (step is set to -1)\t\t",
    " dN\t\tadvance N traces (step is set to N)\t\t\t",
    " %\t\tprint some percentiles of the trace data\t\t",
    " r\t\tprint some ranks (rank[j] = jth smallest datum) \t",
    " p [n1 [n2]]  \ttab plot sample n1 to n2 on current trace\t\t",
    " g [tr1 tr2]  \tximage plot the trace [traces tr1 to tr2]\t",
    " w [tr1 tr2]  \txwigb plot the trace [traces tr1 to tr2]\t",
    " f [tr1 tr2]   ximage plot the amplitude spectra of the trace\t\t",
    " u [tr1 tr2]   apply user pipeline to specified traces ",
    " ! key=val  \tchange a value in a field (e.g. ! tracr=101)\t\t",
    " ?\t\tprint help file\t\t\t\t\t\t",
    " q\t\tquit\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " NOTE: sample numbers are 1-based (first sample is 1).\t\t\t",
    "                                                                       ",
    " 'u 1000000  1000100 suwind >subset.su' will quickly extract a few     ",
    " traces from the middle of a large dataset                             ",
    " \t\t\t\t\t\t\t\t\t",
];

const HELP: &[&str] = &[
    "\t\t\t\t\t", " n\t\tread in trace #n\t", " <CR>\t\tstep\t\t\t",
    " +\t\tnext trace;   step -> +1", " -\t\tprev trace;   step -> -1",
    " dN\t\tadv N traces; step -> N\t", " %\t\tpercentiles\t\t", " r\t\tranks\t\t\t",
    " p [n1 [n2]]  \ttabplot\t\t\t", " w [tr1 tr2]  \twiggle plot\t\t",
    " g [tr1 tr2]  \tgreyscale plot\t\t", " f [tr1 tr2]   amplitude spectra plot  ",
    " u [tr1 tr2]   user specified pipeline ", " ! key=val  \tmodify field\t\t",
    " ?\t\tprint this file\t\t", " q\t\tquit\t\t\t", "\t\t\t\t\t",
];
unsafe extern "C" {
    fn access(p: *const c_char, m: c_int) -> c_int;
    fn sscanf(s: *const c_char, f: *const c_char, ...) -> c_int;
    fn system(s: *const c_char) -> c_int;
    fn qsort(base: *mut std::ffi::c_void, n: usize, size: usize, cmp: unsafe extern "C" fn(*const std::ffi::c_void, *const std::ffi::c_void) -> c_int);
}
static mut SORT_DATA: *const f32 = std::ptr::null();
unsafe extern "C" fn cmp_indirect(r: *const std::ffi::c_void, s: *const std::ffi::c_void) -> c_int {
    unsafe {
        let diff = *SORT_DATA.add(*(r as *const i32) as usize) - *SORT_DATA.add(*(s as *const i32) as usize);
        if diff > 0.0 { 1 } else if diff < 0.0 { -1 } else { 0 }
    }
}
fn scan(s: &str, a: &mut i32, b: &mut i32, user: &mut [u8; 8192]) -> i32 {
    let s = CString::new(s).unwrap_or_default();
    // fgets converts invalid UTF-8 lossily, so its result can exceed 8191 bytes.
    unsafe { sscanf(s.as_ptr(), c"%d %d %8191[^\n]".as_ptr(), a, b, user.as_mut_ptr()) }
}

#[cfg(test)]
mod tests {
    use super::scan;

    #[test]
    fn scan_bounds_lossy_input() {
        let input = format!("1 2 {}", "\u{fffd}".repeat(8191));
        let (mut a, mut b) = (0, 0);
        let mut user = [0u8; 8192];
        assert_eq!(scan(&input, &mut a, &mut b, &mut user), 3);
        assert_eq!((a, b), (1, 2));
        assert_eq!(user[8191], 0);
    }
}
fn help() { for line in HELP { fprintf!(stderr(), "%s\n", *line); } }
fn edxplot(mode: i32, input: &str, infp: &CFile, tr: &mut Trace) {
    let directory = temporary_filename("");
    let tmpfp = efopen(&directory, "w+");
    let (mut n1, mut n2) = (0, 0);
    let mut user = [0u8; 8192];
    match scan(input, &mut n1, &mut n2, &mut user) {
        1 => { fgettra(infp, tr, n1 - 1); fputtr(&tmpfp, tr); }
        2 | 3 => {
            if n2 < n1 { warn!("must specify n1=%d < n2=%d", n1, n2); return; }
            for i in n1 - 1..=n2 - 1 { fgettra(infp, tr, i); fputtr(&tmpfp, tr); }
        }
        _ => fputtr(&tmpfp, tr),
    }
    efclose(&tmpfp);
    let user = String::from_utf8_lossy(&user[..user.iter().position(|x| *x == 0).unwrap_or(user.len())]);
    let cmd = match mode {
        0 => format!("(suximage perc=99 <{directory} ;sleep 10;rm {directory})&"),
        1 => format!("(suxwigb <{directory} perc=99;sleep 10;rm {directory})  &"),
        2 => format!("(suspecfx <{directory} | suximage legend=1 cmap=hue perc=99;sleep 10;rm {directory} )&"),
        _ => format!("(<{directory} {user} ;sleep 10;rm {directory})&"),
    };
    if let Ok(cmd) = CString::new(cmd) { unsafe { system(cmd.as_ptr()); } }
    printheader(tr);
}
fn main() {
    let argv = args();
    initargs(&argv); requestdoc(1, SDOC);
    if argv.len() > 2 { err!("only one filename argument is allowed"); }
    let tty = efopen("/dev/tty", "r");
    let write_ok;
    let infp = if !isatty(STDIN) {
        write_ok = false; warn!("! examine only (no header editing from STDIN)\n"); stdin()
    } else {
        let name = &argv[1]; let path = CString::new(name.as_str()).unwrap_or_default();
        if unsafe { access(path.as_ptr(), 6) } == 0 {
            write_ok = true; fopen(name, "r+").unwrap_or_else(|| efopen(name, "r+"))
        } else if unsafe { access(path.as_ptr(), 4) } == 0 {
            write_ok = false; warn!("! %s is readonly (no header editing)\n", name);
            fopen(name, "r").unwrap_or_else(|| efopen(name, "r"))
        } else { err!("can't open %s for reading", name); }
    };
    fseeko(&infp, 0, SEEK_CUR);
    let mut tr = Trace::new();
    let ntr = fgettra(&infp, &mut tr, 0); let nt = tr.ns as i32;
    let mut rank: Vec<i32> = (0..nt).collect();
    printf!("%d traces in input file\n", ntr);
    let mut itr = ntr - 1; fgettra(&infp, &mut tr, itr); printheader(&tr);
    printf!("> "); efflush(&stdout());
    let mut step = -1; let (mut p1, mut p2) = (0, 0);
    while let Some(userin) = fgets(8192, &tty) {
        let uptr = userin.trim_start_matches([' ', '\t']);
        let cmd = uptr.as_bytes().first().copied().unwrap_or(0);
        match cmd {
            b'0'..=b'9' => { itr = eatoi(uptr) - 1; fgettra(&infp, &mut tr, itr); printheader(&tr); }
            b'q' => exit(CWP_Exit()),
            b'p' => {
                let mut user = [0u8; 8192];
                if scan(&uptr[1..], &mut p1, &mut p2, &mut user) >= 2 {
                    if p2 < p1 { warn!("need p1=%d < p2=%d", p1, p2); exit(CWP_Exit()); }
                    p2 = p2.min(nt);
                } else { p2 = (p1 + 19).min(nt); }
                if p1 >= nt || p1 <= 0 { p1 = 1; }
                tabplot(&tr.data, p1 - 1, p2 - 1); p1 = p2;
            }
            b'g' | b'w' | b'f' | b'u' => {
                let mode = match cmd { b'g' => 0, b'w' => 1, b'f' => 2, _ => 3 };
                edxplot(mode, &uptr[1..], &infp, &mut tr);
            }
            b'+' => { itr += 1; if itr > ntr - 1 { itr = ntr - 2; printf!("\nBounced off end of data:\n\n"); } fgettra(&infp, &mut tr, itr); printheader(&tr); step = 1; }
            b'-' => { itr -= 1; if itr < 0 { itr = 2; printf!("\nBounced off end of data:\n\n"); } fgettra(&infp, &mut tr, itr); printheader(&tr); step = -1; }
            b'\n' | b'd' => {
                if cmd == b'd' { step = eatoi(&uptr[1..]); }
                itr += step;
                if itr < 0 || itr > ntr - 1 { step *= -1; itr += 2 * step; printf!("\nBounced off end of data:\n\n"); }
                fgettra(&infp, &mut tr, itr); printheader(&tr);
            }
            b'r' | b'%' => {
                unsafe {
                    SORT_DATA = tr.data.as_ptr();
                    qsort(rank.as_mut_ptr().cast(), nt as usize, 4, cmp_indirect);
                }
                if cmd == b'r' {
                    let mut i = 0usize; printf!(" rank[%d] = %8.2e", i as i32 + 1, tr.data[rank[i] as usize]);
                    i = nt as usize / 20; printf!(" rank[%d] = %8.2e", i as i32 + 1, tr.data[rank[i] as usize]);
                    i = nt as usize / 2 - i; printf!(" rank[%d] = %8.2e", i as i32 + 1, tr.data[rank[i] as usize]); printf!("\n");
                    i = nt as usize - 1 - i; printf!(" rank[%d] = %8.2e", i as i32 + 1, tr.data[rank[i] as usize]);
                    i = nt as usize - 1 - nt as usize / 20; printf!(" rank[%d] = %8.2e", i as i32 + 1, tr.data[rank[i] as usize]);
                    i = nt as usize - 1; printf!(" rank[%d] = %8.2e", i as i32 + 1, tr.data[rank[i] as usize]);
                    printf!("\nmin is at sample %d,  max at %d\n", rank[0] as i32 + 1, rank[nt as usize - 1] as i32 + 1);
                } else {
                    for (i, q) in [(1, 0.01), (5, 0.05), (25, 0.25), (50, 0.50), (75, 0.75), (95, 0.95), (99, 0.99)] {
                        let iq = (q * nt as f64 - 0.5) as usize;
                        if i == 1 { printf!(" %dst percentile is %8.2e\n", i + 1, tr.data[rank[iq] as usize]); }
                        else if i == 5 { printf!(" %dth percentile is %8.2e\n", i + 1, tr.data[rank[iq] as usize]); }
                        else { printf!("%dth percentile is %8.2e\n", i + 1, tr.data[rank[iq] as usize]); }
                    }
                    printf!("min at sample %d equals %8.2e\n", rank[0] as i32 + 1, tr.data[rank[0] as usize]);
                    printf!("max at sample %d equals %8.2e\n", rank[nt as usize - 1] as i32 + 1, tr.data[rank[nt as usize - 1] as usize]);
                }
            }
            b'!' => {
                if write_ok {
                    if let Some((key, value)) = userin[1..].trim_start().split_once('=') {
                        let key = key.trim_end(); let mut val = Value::default();
                        atoval(hdtype(key), value, &mut val); puthdval(&mut tr, key, &val);
                        let nsegy = nt as i64 * 4 + HDRBYTES as i64;
                        fseeko(&infp, -nsegy, SEEK_CUR); fwrite(tr.hdr.as_bytes(), &infp); fseeko(&infp, nsegy, SEEK_CUR); printheader(&tr);
                    }
                } else { warn!("file not writable"); }
            }
            b'?' => help(),
            _ => { warn!("bad key %s\n%s", uptr, ""); help(); }
        }
        printf!("> "); efflush(&stdout());
    }
    exit(CWP_Exit());
}
