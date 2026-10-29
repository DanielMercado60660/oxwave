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

//! Build script: link the ORIGINAL SU static libraries (libsu, libpar,
//! libcwp in that order, plus libm) from `$OXWAVE_SU_LIB` (default
//! `<workspace>/vendor/su-build/lib`), and compile `csrc/shim.c` (copies of
//! SU functions that are `static` in programs) with `/usr/bin/cc` using the
//! reference SU build's flags. std only; no `cc` crate.

use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let ws = manifest.parent().unwrap().parent().unwrap().to_path_buf();
    let lib = std::env::var("OXWAVE_SU_LIB")
        .map(PathBuf::from)
        .unwrap_or_else(|_| ws.join("vendor").join("su-build").join("lib"));
    let include = lib.parent().unwrap().join("include");
    println!("cargo:rerun-if-env-changed=OXWAVE_SU_LIB");
    println!("cargo:rerun-if-changed=csrc/shim.c");
    for l in ["libsu.a", "libpar.a", "libcwp.a"] {
        let p = lib.join(l);
        if !p.exists() {
            // No reference SU build: let `cargo build` of the workspace succeed. Only the
            // differential tests need these libraries, and they fail to link without them.
            println!(
                "cargo:warning=oxwave-cref: {} not found; differential tests against original SU \
                 will not link (build reference SU, see vendor/BUILD_NOTES.md, or set OXWAVE_SU_LIB)",
                p.display()
            );
            return;
        }
        println!("cargo:rerun-if-changed={}", p.display());
    }

    // compile the shim exactly like the reference SU build compiles SU
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let obj = out.join("shim.o");
    let st = Command::new("/usr/bin/cc")
        .args(["-c", "-g", "-std=c90", "-m64", "-ansi", "-Wno-long-long", "-ffp-contract=off"])
        .args(["-DCWP_LITTLE_ENDIAN", "-DSUXDR"])
        .arg(format!("-I{}", include.display()))
        .arg(manifest.join("csrc").join("shim.c"))
        .arg("-o")
        .arg(&obj)
        .status()
        .expect("oxwave-cref: cannot run /usr/bin/cc");
    assert!(st.success(), "oxwave-cref: compiling csrc/shim.c failed");
    let ar = out.join("liboxcrefshim.a");
    let _ = std::fs::remove_file(&ar);
    let st = Command::new("/usr/bin/ar").arg("crs").arg(&ar).arg(&obj).status().expect("ar");
    assert!(st.success(), "oxwave-cref: ar failed");

    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-search=native={}", lib.display());
    println!("cargo:rustc-link-lib=static=oxcrefshim");
    println!("cargo:rustc-link-lib=static=su");
    println!("cargo:rustc-link-lib=static=par");
    println!("cargo:rustc-link-lib=static=cwp");
    println!("cargo:rustc-link-lib=dylib=m");
}
