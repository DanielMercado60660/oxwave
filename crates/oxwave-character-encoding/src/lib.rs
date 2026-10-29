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

//! Port of cwp/lib/zasc.c.
//! Author: Stew Levin of Mobil, 1997; translated by f2c.

/// Convert EBCDIC bytes to ASCII using the original translation table.
///
/// C prototype: `int zasc(char *ainput, char *aoutput, integer nchar);`
///
/// The input slice determines the number of bytes translated. The output must
/// have at least that many bytes; any remaining output bytes are untouched.
pub fn zasc(ainput: &[u8], aoutput: &mut [u8]) -> i32 {
    const T2: [i32; 96] = [
        64, 79, 127, 123, 91, 108, 80, 125, 77, 93, 92, 78, 107, 96,
        75, 97, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 122,
        94, 76, 126, 110, 111, 124, 193, 194, 195, 196, 197, 198, 199,
        200, 201, 209, 210, 211, 212, 213, 214, 215, 216, 217, 226, 227,
        228, 229, 230, 231, 232, 233, 74, 224, 90, 95, 109, 64, 129,
        130, 131, 132, 133, 134, 135, 136, 137, 145, 146, 147, 148, 149,
        150, 151, 152, 153, 162, 163, 164, 165, 166, 167, 168, 169, 192,
        106, 208, 161, 64,
    ];
    assert!(aoutput.len() >= ainput.len());
    for i in 0..ainput.len() {
        let kk = ainput[i] as std::ffi::c_char as i32;
        let mut this = b' ';
        for k in 0..96 {
            if kk == T2[k] {
                this = (32 + k) as u8;
                break;
            }
        }
        aoutput[i] = this;
    }
    0
}
