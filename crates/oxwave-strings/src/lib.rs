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

//! Port of cwp/lib/strstuff.c.
//! Credits: John Stockwell; John Stockwell and Jack K. Cohen; reversal based
//! on public-domain strrev by Bob Stout.

/// Duplicate bytes before the first NUL, returning a vector without the NUL.
/// The slice end also acts as a terminator. C: `char *cwp_strdup(char *str)`.
pub fn cwp_strdup(str: &[u8]) -> Vec<u8> {
    let mut len = 0;
    while len < str.len() && str[len] != 0 { len += 1; }
    let mut copy = vec![0; len];
    copy.copy_from_slice(&str[..len]);
    copy
}

/// Reverse bytes before the first NUL, returning a vector without the NUL.
/// The slice end also acts as a terminator. C: `char *cwp_strrev(char *str)`.
pub fn cwp_strrev(str: &[u8]) -> Vec<u8> {
    let mut len = 0;
    while len < str.len() && str[len] != 0 { len += 1; }
    let mut revstr = vec![0; len];
    let mut i = len;
    let mut j = 0;
    while i > 0 {
        i -= 1;
        revstr[j] = str[i];
        j += 1;
    }
    revstr
}

/// Copy bytes up to a comma or NUL and terminate `t` with a NUL. The slice
/// end also terminates `s`; `t` must have space for the copy and terminator.
/// C prototype: `void strchop(char *s, char *t)`.
pub fn strchop(s: &[u8], t: &mut [u8]) {
    let mut i = 0;
    while i < s.len() && s[i] != b',' && s[i] != 0 {
        t[i] = s[i];
        i += 1;
    }
    t[i] = 0;
}
