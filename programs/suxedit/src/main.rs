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

//! SUXEDIT - examine segy diskfiles and edit headers. Port of su/main/headers/suxedit.c (revision 1.29).
//! Credits: SEP: Einar Kjartansson, Shuki Ronen, Stew Levin; CWP: Jack K. Cohen.

use oxwave_core::prelude::*;
use oxwave_tempfile::temporary_filename;
use oxwave_trace_display::tabplot;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUXEDIT - examine segy diskfiles and edit headers\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " suxedit diskfile  (open for possible header modification if writable)\t",
    " suxedit <diskfile  (open read only)\t\t\t\t\t",
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
    " g [tr1 tr2] [\"opts\"] \twiggle plot (graph) the trace\t\t",
    "\t\t\t\t[traces tr1 to tr2]\t\t\t",
    " f\t\twiggle plot the Fourier transform of the trace\t\t",
    " ! key=val  \tchange a value in a field (e.g. ! tracr=101)\t\t",
    " ?\t\tprint help file\t\t\t\t\t\t",
    " q\t\tquit\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " NOTE: sample numbers are 1-based (first sample is 1).\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

const HELP: &[&str] = &[
    "\t\t\t\t\t", " n\t\tread in trace #n\t", " <CR>\t\tstep\t\t\t",
    " +\t\tnext trace;   step -> +1", " -\t\tprev trace;   step -> -1",
    " dN\t\tadv N traces; step -> N\t", " %\t\tpercentiles\t\t", " r\t\tranks\t\t\t",
    " p [n1 [n2]]  \ttabplot\t\t\t", " g [tr1 tr2] [\"opts\"]\twiggle plot",
    " f\t\twig plot Fourier Transf ", " ! key=val  \tmodify field\t\t", " ?\t\tprint this file\t\t",
    " q\t\tquit\t\t\t", "\t\t\t\t\t",
];
unsafe extern "C" {
    fn access(path: *const i8, mode: i32) -> i32;
    fn sscanf(input: *const i8, format: *const i8, ...) -> i32;
    fn qsort(base: *mut std::ffi::c_void, n: usize, size: usize, cmp: unsafe extern "C" fn(*const std::ffi::c_void, *const std::ffi::c_void) -> i32);
}
static mut SORT_DATA: *const f32 = std::ptr::null();
unsafe extern "C" fn cmp_indirect(a: *const std::ffi::c_void, b: *const std::ffi::c_void) -> i32 {
    unsafe {
        let diff = *SORT_DATA.add(*(a as *const i32) as usize) - *SORT_DATA.add(*(b as *const i32) as usize);
        if diff > 0.0 { 1 } else if diff < 0.0 { -1 } else { 0 }
    }
}
fn can_access(path: &str, mode: i32) -> bool {
    let c = std::ffi::CString::new(path).unwrap();
    unsafe { access(c.as_ptr(), mode) == 0 }
}
fn help() { for s in HELP { fprintf!(stderr(), "%s\n", *s); } }
fn plot(userin: &str, fourier: bool, infp: &CFile, tty: &CFile, tr: &mut Trace) {
    let directory = ""; let filename = temporary_filename(directory); let fp = efopen(&filename, "w+"); let mut opts = String::new();
    if fourier { tr.ntr = 0; fputtr(&fp, tr); } else {
        let mut parts = userin[1..].split_whitespace();
        let n1 = parts.next().and_then(|s| s.parse::<i32>().ok());
        let n2 = parts.next().and_then(|s| s.parse::<i32>().ok());
        match (n1, n2) {
            (Some(a), Some(b)) => {
                if b < a { warn!("must specify n1=%d < n2=%d", a, b); return; }
                opts = parts.collect::<Vec<_>>().join(" ");
                for i in a-1..=b-1 { fgettra(infp, tr, i); tr.ntr = 0; fputtr(&fp, tr); }
            }
            (Some(a), None) => { fgettra(infp, tr, a-1); tr.ntr = 0; fputtr(&fp, tr); }
            _ => fputtr(&fp, tr),
        }
    }
    erewind(&fp);
    let cmd = if fourier { format!("suspecfx <{} | suxwigb", directory) } else { format!("suxwigb <{} {}", directory, opts) };
    let _ = std::process::Command::new("sh").arg("-c").arg(cmd).status();
    efclose(&fp); eremove(&directory);
    fprintf!(stderr(), "\npress return key to continue\n"); fgetc(tty); printheader(tr);
}
fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    if argv.len() > 2 { err!("only one filename argument is allowed"); }
    let tty = efopen("/dev/tty", "r");
    let (infp, write_ok) = if !isatty(0) {
        warn!("! examine only (no header editing from STDIN)\n"); (stdin(), false)
    } else if can_access(&argv[1], 6) { (efopen(&argv[1], "r+"), true) }
    else if can_access(&argv[1], 4) { warn!("! %s is readonly (no header editing)\n", argv[1]); (efopen(&argv[1], "r"), false) }
    else { err!("can't open %s for reading", argv[1]); };
    let mut tr = Trace::new(); let ntr = fgettra(&infp, &mut tr, 0); let nt = tr.ns as i32;
    let mut rank: Vec<i32> = (0..nt).collect();
    printf!("%d traces in input file\n", ntr);
    let mut itr = ntr - 1; let mut step = -1; let mut p1 = 0i32;
    fgettra(&infp, &mut tr, itr); printheader(&tr); printf!("> "); efflush(&stdout());
    while let Some(userin) = fgets(8192, &tty) {
        match userin.as_bytes()[0] {
            b'0'..=b'9' => {
                itr = eatoi(&userin) - 1;
                if itr < 0 || itr > ntr-1 { warn!("no such trace"); itr = ntr-1; }
                fgettra(&infp, &mut tr, itr); printheader(&tr);
            }
            b'q' => exit(CWP_Exit()),
            b'p' => {
                let input = std::ffi::CString::new(&userin[1..]).unwrap();
                let mut p2 = 0i32;
                let count = unsafe { sscanf(input.as_ptr(), c"%d %d".as_ptr(), &mut p1, &mut p2) };
                if count == 2 {
                    if p2 < p1 { warn!("need p1=%d < p2=%d", p1, p2); exit(0); }
                    p2 = p2.min(nt);
                } else { p2 = (p1+19).min(nt); }
                if p1 >= nt || p1 <= 0 { p1 = 1; }
                tabplot(&tr.data, p1-1, p2-1); p1 = p2;
            }
            b'g' => plot(&userin, false, &infp, &tty, &mut tr),
            b'f' => plot(&userin, true, &infp, &tty, &mut tr),
            b'+' | b'-' | b'\n' | b'd' => {
                let c = userin.as_bytes()[0];
                if c == b'+' { itr += 1; } else if c == b'-' { itr -= 1; } else { if c == b'd' { step = eatoi(&userin[1..]); } itr += step; }
                if itr > ntr-1 { step = if c == b'+' { -1 } else { -step }; itr = ntr-1; printf!("\nBounced off end of data:\n\n"); }
                if itr < 0 { step = if c == b'-' { 1 } else { -step }; itr = 0; printf!("\nBounced off end of data:\n\n"); }
                fgettra(&infp, &mut tr, itr); printheader(&tr);
                if c == b'+' { step = 1; } else if c == b'-' { step = -1; }
            }
            b'r' | b'%' => {
                unsafe { SORT_DATA = tr.data.as_ptr(); qsort(rank.as_mut_ptr().cast(), nt as usize, 4, cmp_indirect); }
                if userin.as_bytes()[0] == b'r' {
                    let a = nt/20; let b = nt/2-a;
                    for (j,i) in [0,a,b,nt-1-b,nt-1-a,nt-1].iter().enumerate() { printf!(" rank[%d] = %8.2e", i+1, tr.data[rank[*i as usize] as usize]); if j == 2 || j == 5 { printf!("\n"); } }
                    printf!("min is at sample %d,  max at %d\n", rank[0]+1, rank[nt as usize-1]+1);
                } else {
                    for (i,q) in [(1,0.01),(5,0.05),(25,0.25),(50,0.50),(75,0.75),(95,0.95),(99,0.99)] {
                        let iq = (q * nt as f64 - 0.5) as usize;
                        if i == 1 { printf!(" %dst percentile is %8.2e\n", i+1, tr.data[rank[iq] as usize]); }
                        else if i == 5 { printf!(" %dth percentile is %8.2e\n", i+1, tr.data[rank[iq] as usize]); }
                        else { printf!("%dth percentile is %8.2e\n", i+1, tr.data[rank[iq] as usize]); }
                    }
                    printf!("min at sample %d equals %8.2e\n", rank[0]+1, tr.data[rank[0] as usize]);
                    printf!("max at sample %d equals %8.2e\n", rank[nt as usize-1]+1, tr.data[rank[nt as usize-1] as usize]);
                }
            }
            b'!' => {
                if !write_ok { warn!("file not writable"); } else if let Some((key, value)) = userin[1..].trim_start().split_once('=') {
                    let key = key.trim_end(); let mut val = Value::default(); atoval(hdtype(key), value, &mut val); puthdval(&mut tr, key, &val);
                    let databytes = nt as i64 * 4; efseeko(&infp, -(databytes + HDRBYTES as i64), SEEK_CUR);
                    efwrite(tr.hdr.as_bytes(), &infp); efseeko(&infp, databytes, SEEK_CUR); printheader(&tr);
                }
            }
            b'?' => help(),
            _ => { warn!("bad key %s\n", userin); help(); }
        }
        printf!("> "); efflush(&stdout());
    }
    exit(CWP_Exit());
}
