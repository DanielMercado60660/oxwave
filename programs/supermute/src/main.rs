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
//! SUPERMUTE - permute or transpose a 3d datacube.
//! Port of su/main/operations/supermute.c (Revision 1.5).
//! Credits: VT: Matthias Imhof.

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUPERMUTE - permute or transpose a 3d datacube\t \t\t",
    " \t\t\t\t\t\t\t\t\t",
    " supermute <stdin >sdout\t \t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Required parameters:\t\t\t\t\t\t\t",
    " none\t\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " n1=ns from header\t\tnumber of samples in the fast direction\t",
    " n2=ntr from header\t\tnumber of samples in the med direction\t",
    " n3=1\t\t\t\tnumber of samples in the slow direction\t",
    " \t\t\t\t\t\t\t\t\t",
    " o1=1\t\t\t\tnew fast direction\t\t\t",
    " o2=2\t\t\t\tnew med direction\t\t\t",
    " o3=3\t\t\t\tnew slow direction\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " d1=1\t\t\t\toutput interval in new fast direction\t",
    " d2=1\t\t\t\toutput interval in new med direction\t",
    " d3=1\t\t\t\toutput interval in new slow direction\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " header fields d1 and d2 default to d1=1.0 and d2=1.0\t\t\t",
];

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut n = [0i32; 4];
    let mut m = [0i32; 4];
    let mut o = [0i32; 4];
    let mut d = [0i32; 4];
    let mut p1 = [0i32; 4];
    let mut p2 = [0i32; 4];
    let mut p3 = [0i32; 4];
    let mut verbose = 0;
    let mut check = 0;
    initargs(&argv);
    requestdoc(1, SDOC);
    getparint("verbose", &mut verbose);
    getparint("check", &mut check);
    if gettr(&mut tr) == 0 { err!("cannot get first trace"); }
    if getparint("n1", &mut n[1]) == 0 { n[1] = tr.ns as i32; }
    if getparint("n2", &mut n[2]) == 0 { n[2] = tr.ntr; }
    if getparint("n3", &mut n[3]) == 0 { n[3] = 1; }
    for i in 1..=3 {
        if n[i] == 0 { err!("invalid dimensions: %d %d %d", n[1], n[2], n[3]); }
    }
    if verbose != 0 { warn!("old dimensions: %d %d %d", n[1], n[2], n[3]); }
    if getparint("o1", &mut o[1]) == 0 { o[1] = 1; }
    if getparint("o2", &mut o[2]) == 0 { o[2] = 2; }
    if getparint("o3", &mut o[3]) == 0 { o[3] = 3; }
    for i in 1..=3 {
        if o[i] < 1 || o[i] > 3 { err!("invalid permutation %d: %d", i as i32, o[i]); }
    }
    if verbose != 0 { warn!("new order: %d %d %d", o[1], o[2], o[3]); }
    for i in 1..=3 { m[i] = n[o[i] as usize]; }
    if verbose != 0 { warn!("new dimension: %d %d %d", m[1], m[2], m[3]); }
    if getparint("d1", &mut d[1]) == 0 { d[1] = 1; }
    if getparint("d2", &mut d[2]) == 0 { d[2] = 1; }
    if getparint("d3", &mut d[3]) == 0 { d[3] = 1; }
    checkpars();
    for i in 1..=3 {
        if d[i] == 0 { err!("invalid interval %d: %d", i as i32, d[i]); }
    }
    if verbose != 0 { warn!("new interval: %d %d %d", d[1], d[2], d[3]); }
    p1[o[1] as usize] = 1;
    p2[o[2] as usize] = 1;
    p3[o[3] as usize] = 1;
    if verbose != 0 {
        warn!("perm 1: %d %d %d", p1[1], p1[2], p1[3]);
        warn!("perm 2: %d %d %d", p2[1], p2[2], p2[3]);
        warn!("perm 3: %d %d %d", p3[1], p3[2], p3[3]);
    }
    let len = (m[1] as usize).checked_mul(m[2] as usize).and_then(|x| x.checked_mul(m[3] as usize));
    let mut data = Vec::<f32>::new();
    if len.is_none() || data.try_reserve_exact(len.unwrap_or(usize::MAX)).is_err() {
        err!("not enough memory: %d %d %d", m[1], m[2], m[3]);
    }
    data.resize(len.unwrap(), 0.0);
    if verbose != 0 { warn!("allocated %d*%d*%d floats", m[1], m[2], m[3]); }
    for k in 1..=n[3] {
        for j in 1..=n[2] {
            if k != 1 || j != 1 {
                if check != 0 { warn!("read trace %d %d", k, j); }
                if gettr(&mut tr) == 0 {
                    err!("wrong amount of data: %d, %d (%d), %d (%d)", n[1], n[2], j, n[3], k);
                }
            }
            for i in 1..=n[1] {
                let mut p = p1[1]*i + p1[2]*j + p1[3]*k;
                let mut q = p2[1]*i + p2[2]*j + p2[3]*k;
                let mut r = p3[1]*i + p3[2]*j + p3[3]*k;
                if d[1] < 0 { p = m[1]-p+1; }
                if d[2] < 0 { q = m[2]-q+1; }
                if d[3] < 0 { r = m[3]-r+1; }
                if verbose == 2 { fprintf!(stderr(), "%4d %4d %4d to %4d %4d %4d\n", i,j,k,p,q,r); }
                if check != 0 && (p < 1 || p > m[1] || q < 1 || q > m[2] || r < 1 || r > m[3]) {
                    err!("access error: %d (%d) %d, %d (%d) %d, %d (%d) %d", 1,p,m[1],1,q,m[2],1,r,m[3]);
                }
                let idx = (((r-1)*m[2]+q-1)*m[1]+p-1) as usize;
                data[idx] = tr.data[(i-1) as usize];
            }
            if check != 0 { warn!("next trace"); }
        }
    }
    if verbose == 2 {
        fprintf!(stderr(), "m1 m2 m3: %d %d %d\n", m[1],m[2],m[3]);
        for r in 0..m[3] { for q in 0..m[2] {
            for p in 0..m[1] { fprintf!(stderr(), "%16.6f", data[((r*m[2]+q)*m[1]+p) as usize]); }
            fprintf!(stderr(), "\n");
        } fprintf!(stderr(), "\n"); }
    }
    tr.ns = ((m[1]+d[1].abs()-1)/d[1].abs()) as u16;
    tr.ntr = (m[2]+(d[2]-1).abs()/d[2].abs()) * (m[3]+(d[3]-1).abs()/d[3].abs());
    tr.d1 = 1.0; tr.f1 = 1.0; tr.d2 = 1.0; tr.f2 = 1.0;
    let mut r = 1;
    while r <= m[3] {
        let mut q = 1;
        while q <= m[2] {
            let mut t = 0usize;
            let mut p = 1;
            while p <= m[1] {
                let value = data[(((r-1)*m[2]+q-1)*m[1]+p-1) as usize];
                if verbose == 3 { fprintf!(stderr(), "%4d %4d %4d %f\n", p,q,r,value); }
                tr.data[t] = value;
                t += 1;
                p += d[1].abs();
            }
            puttr(&mut tr);
            q += d[2].abs();
        }
        r += d[3].abs();
    }
    exit(CWP_Exit());
}
