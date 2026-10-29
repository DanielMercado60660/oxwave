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

//! Trace I/O against the original binaries: our fgettr/fputtr pass-through
//! must produce byte-identical streams to a C gettr/puttr pass-through
//! (`sushw key=nvs a=0`), including packed (supack1/supack2) traces.

use std::path::PathBuf;
use std::process::Command;

fn su_bin() -> PathBuf {
    let ws = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    std::env::var("OXWAVE_CWPROOT").map(PathBuf::from).unwrap_or(ws.join("vendor/su-build")).join("bin")
}

fn sh(cmd: &str) -> Vec<u8> {
    let bin = su_bin();
    let path = format!("{}:{}", bin.display(), std::env::var("PATH").unwrap_or_default());
    let out = Command::new("bash")
        .arg("-o").arg("pipefail").arg("-c").arg(cmd)
        .env("PATH", path)
        .env("CWPROOT", bin.parent().unwrap())
        .output()
        .expect("bash");
    assert!(out.status.success(), "{} failed: {}", cmd, String::from_utf8_lossy(&out.stderr));
    out.stdout
}

fn passthrough(input: &[u8]) -> Vec<u8> {
    use oxwave_core::cio::*;
    use oxwave_core::trace::*;
    let inp = etmpfile();
    efwrite(input, &inp);
    erewind(&inp);
    let out = etmpfile();
    let mut tr = Trace::new();
    while fgettr(&inp, &mut tr) != 0 {
        fputtr(&out, &mut tr);
    }
    erewind(&out);
    let mut v = Vec::new();
    let mut buf = [0u8; 4096];
    loop {
        let n = fread(&mut buf[..], &out);
        if n == 0 {
            break;
        }
        v.extend_from_slice(&buf[..n]);
    }
    v
}

#[test]
fn passthrough_matches_c() {
    let cases = [
        "suplane nt=64 ntr=8",
        "suplane nt=17 ntr=3 | sushw key=d1,f1,d2,f2,ungpow,unscale,ntr,mark,shortpad a=1.5,-2.25,1e20,-0.0,3,4,5,6,7",
        "suplane nt=8 ntr=4 | sushw key=trid a=1 | sushw key=offset,gx,sx a=-1,-2147483648,2147483647",
        "suplane nt=8 ntr=2 | supack1",
        "suplane nt=8 ntr=2 | supack2",
        "suspike nt=100 ntr=5 | sufilter f=10,20,30,40",
    ];
    for c in cases {
        let input = sh(c);
        let want = sh(&format!("{} | sushw key=nvs a=0", c));
        let got = passthrough(&input);
        assert_eq!(got.len(), want.len(), "{}", c);
        assert!(got == want, "pass-through differs from C for: {}", c);
    }
}
