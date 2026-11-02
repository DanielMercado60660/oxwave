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

//! Port of par/lib/modeling.c. Credits: CWP Dave Hale, 09/17/91, Colorado School of Mines.
#![allow(non_snake_case)]
use oxwave_core::par::{atof, countparname, getnparstring};
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::csplin;
use oxwave_inverse_interp::yxtoxy;
use oxwave_sinc_interp::mksinc;
use std::sync::OnceLock;
#[derive(Clone, Debug, PartialEq)]
pub struct ReflectorInput {
    pub a: f32,
    pub x: Vec<f32>,
    pub z: Vec<f32>,
}
#[derive(Clone, Debug)]
pub struct ReflectorSegment {
    pub x: f32,
    pub z: f32,
    pub s: f32,
    pub c: f32,
}
#[derive(Clone, Debug)]
pub struct Reflector {
    pub ns: i32,
    pub ds: f32,
    pub a: f32,
    pub rs: Vec<ReflectorSegment>,
}
#[derive(Clone, Debug)]
pub struct Wavelet {
    pub lw: i32,
    pub iw: i32,
    pub wv: Vec<f32>,
}
/// C: `int decodeReflector(char *string, float *aPtr, int *nxzPtr, float **xPtr, float **zPtr)`
pub fn decodeReflector(string: &str) -> Option<ReflectorInput> {
    // strtok skips leading delimiters but only overwrites the delimiter that
    // ends a token. C later walks that modified buffer with strlen(s)+1,
    // so repeated delimiters must remain attached to the following value.
    fn strtok(buffer: &mut [u8], next: &mut usize, delimiters: &[u8]) -> Option<usize> {
        while buffer[*next] != 0 && delimiters.contains(&buffer[*next]) {
            *next += 1;
        }
        if buffer[*next] == 0 {
            return None;
        }
        let start = *next;
        while buffer[*next] != 0 && !delimiters.contains(&buffer[*next]) {
            *next += 1;
        }
        if buffer[*next] != 0 {
            buffer[*next] = 0;
            *next += 1;
        }
        Some(start)
    }
    fn value(buffer: &[u8], start: usize) -> (f32, usize) {
        let len = buffer[start..].iter().position(|&b| b == 0).unwrap();
        let token = std::str::from_utf8(&buffer[start..start + len]).unwrap();
        (atof(token) as f32, start + len + 1)
    }

    let string = string.split('\0').next().unwrap();
    let mut buffer = string.as_bytes().to_vec();
    buffer.push(0);
    let mut next = 0;
    let a = if let Some(i) = string.find(':') {
        if string.find([',', ';']).unwrap_or(string.len()) <= i {
            return None;
        }
        let start = strtok(&mut buffer, &mut next, b":")?;
        value(&buffer, start).0
    } else {
        1.0
    };
    let mut start = strtok(&mut buffer, &mut next, b",;")?;
    let mut count = 1;
    while strtok(&mut buffer, &mut next, b",;").is_some() {
        count += 1;
    }
    if count < 4 || count % 2 != 0 {
        return None;
    }
    let mut x = Vec::new();
    let mut z = Vec::new();
    for _ in 0..count / 2 {
        let (v, following) = value(&buffer, start);
        x.push(v);
        let (v, following) = value(&buffer, following);
        z.push(v);
        start = following;
    }
    Some(ReflectorInput { a, x, z })
}
/// C: `void decodeReflectors(int *nrPtr, float **aPtr, int **nxzPtr, float ***xPtr, float ***zPtr)`
pub fn decodeReflectors() -> Vec<ReflectorInput> {
    let nr = countparname("ref").max(1);
    (1..=nr)
        .map(|i| {
            let mut s = String::new();
            if getnparstring(i, "ref", &mut s) == 0 {
                s = "1:1,2;4,2".into()
            }
            decodeReflector(&s).unwrap_or_else(|| {
                oxwave_core::err!("Reflector number %d specified incorrectly!\n", i)
            })
        })
        .collect()
}
/// C: `void breakReflectors(int *nr, float **ar, int **nu, float ***xu, float ***zu)`
pub fn breakReflectors(input: &[ReflectorInput]) -> Vec<ReflectorInput> {
    let mut out = Vec::new();
    for r in input {
        for i in 0..r.x.len().saturating_sub(1) {
            out.push(ReflectorInput {
                a: r.a,
                x: vec![r.x[i], r.x[i + 1]],
                z: vec![r.z[i], r.z[i + 1]],
            })
        }
    }
    out
}
/// C: `void makeref(float dsmax, int nr, float *ar, int *nu, float **xu, float **zu, Reflector **r)`
pub fn makeref(dsmax: f32, input: &[ReflectorInput]) -> Vec<Reflector> {
    let mut out = Vec::new();
    for r in input {
        let n = r.x.len();
        let u: Vec<f32> = (0..n).map(|i| i as f32).collect();
        let mut xd = vec![[0.0; 4]; n];
        let mut zd = xd.clone();
        csplin(&u, &r.x, &mut xd);
        csplin(&u, &r.z, &mut zd);
        let nuu = 20 * n;
        let duu = (u[n - 1] - u[0]) / (nuu - 1) as f32;
        let mut lengths = vec![0.0f32; nuu];
        let (mut xl, mut zl) = (r.x[0], r.z[0]);
        let mut uu = duu;
        for i in 1..nuu {
            let (mut x, mut z) = ([0.0], [0.0]);
            intcub(0, &u, &xd, &[uu], &mut x);
            intcub(0, &u, &zd, &[uu], &mut z);
            let dx = x[0] - xl;
            let dz = z[0] - zl;
            lengths[i] = (lengths[i - 1] as f64 + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
            xl = x[0];
            zl = z[0];
            uu += duu
        }
        // Both operands of C's addition are float; round before converting
        // to the signed int passed to ealloc1float.
        let ns = (1.0 + lengths[nuu - 1] / dsmax) as i32;
        let ds = lengths[nuu - 1] / ns as f32;
        let fs = (0.5 * ds as f64) as f32;
        let mut us = ealloc1float(ns);
        yxtoxy(duu, 0.0, &lengths, ds, fs, 0.0, (n - 1) as f32, &mut us);
        let mut rs = Vec::new();
        for v in us {
            let (mut x, mut z, mut dx, mut dz) = ([0.0], [0.0], [0.0], [0.0]);
            intcub(0, &u, &xd, &[v], &mut x);
            intcub(0, &u, &zd, &[v], &mut z);
            intcub(1, &u, &xd, &[v], &mut dx);
            intcub(1, &u, &zd, &[v], &mut dz);
            let norm = ((dx[0] * dx[0] + dz[0] * dz[0]) as f64).sqrt();
            rs.push(ReflectorSegment {
                x: x[0],
                z: z[0],
                c: (dx[0] as f64 / norm) as f32,
                s: (-dz[0] as f64 / norm) as f32,
            })
        }
        out.push(Reflector { ns, ds, a: r.a, rs })
    }
    out
}
// C converts the signed count to size_t before allocation. Negative counts
// therefore request an impossible allocation and fail with ENOMEM on Darwin.
fn ealloc1float(n: i32) -> Vec<f32> {
    let mut values = Vec::new();
    values.try_reserve_exact(n as usize).unwrap_or_else(|_| {
        oxwave_core::err!(
            "ealloc.c: malloc failed (%s)",
            oxwave_core::err::strerror(12)
        )
    });
    values.resize(n as usize, 0.0);
    values
}
/// C: `void makericker(float fpeak, float dt, Wavelet **w)`
pub fn makericker(fpeak: f32, dt: f32) -> Wavelet {
    let iw = -(1.0 + 1.0 / (fpeak * dt) as f64) as i32;
    let lw = 1 - 2 * iw;
    let mut wv = ealloc1float(lw);
    let mut t = iw as f32 * dt;
    for sample in &mut wv {
        let mut x = (std::f64::consts::PI * fpeak as f64 * t as f64) as f32;
        x = x * x;
        *sample = ((-(x as f64)).exp() * (1.0 - 2.0 * x as f64)) as f32;
        t += dt
    }
    Wavelet { lw, iw, wv }
}
/// C: `void addsinc(float time, float amp, int nt, float dt, float ft, float *trace)`
pub fn addsinc(time: f32, amp: f32, dt: f32, ft: f32, trace: &mut [f32]) {
    static S: OnceLock<[[f32; 8]; 101]> = OnceLock::new();
    let sinc = S.get_or_init(|| {
        let mut s = [[0.0; 8]; 101];
        for j in 1..100 {
            mksinc(j as f32 / 100.0, &mut s[j])
        }
        s[0][3] = 1.0;
        s[100][4] = 1.0;
        s
    });
    let tn = (time - ft) / dt;
    let jt = tn as i32;
    let js = ((tn - jt as f32) * 100.0) as i32;
    let lo = jt - 3;
    let hi = jt + 4;
    if hi >= 0 && lo < trace.len() as i32 {
        for i in lo.max(0)..=hi.min(trace.len() as i32 - 1) {
            trace[i as usize] += amp * sinc[js as usize][(i - jt + 3) as usize]
        }
    }
}
/// C: `void raylv2(float v00, float dvdx, float dvdz, float x0, float z0, float x, float z, float *c, float *s, float *t, float *q)`
pub fn raylv2(
    v00: f32,
    dvdx: f32,
    dvdz: f32,
    mut x0: f32,
    mut z0: f32,
    mut x: f32,
    mut z: f32,
) -> (f32, f32, f32, f32) {
    if x == x0 && z == z0 {
        return (1.0, 0.0, 0.0, 0.0);
    }
    if dvdx == 0.0 && dvdz == 0.0 {
        x -= x0;
        z -= z0;
        let r = ((x * x + z * z) as f64).sqrt() as f32;
        let o = (1.0 / r as f64) as f32;
        return (z * o, x * o, r / v00, r * v00);
    }
    let (a, cr, sr) = if dvdx != 0.0 {
        let a = ((dvdx * dvdx + dvdz * dvdz) as f64).sqrt() as f32;
        let o = (1.0 / a as f64) as f32;
        let cr = dvdz * o;
        let sr = dvdx * o;
        let tmp = cr * x0 - sr * z0;
        z0 = cr * z0 + sr * x0;
        x0 = tmp;
        let tmp = cr * x - sr * z;
        z = cr * z + sr * x;
        x = tmp;
        (a, cr, sr)
    } else {
        (dvdz, 0.0, 0.0)
    };
    let v0 = v00 + a * z0;
    let v = v00 + a * z;
    if v0 < 0.0 || v < 0.0 {
        return (1.0, 0.0, f32::MAX, f32::MAX);
    }
    x -= x0;
    z -= z0;
    let (mut c, mut s, t, q);
    if ((x * x) as f64) < 0.000001 * z as f64 * z as f64 {
        if z > 0.0 {
            s = 0.0;
            c = 1.0;
            t = (((v / v0) as f64).ln() / a as f64) as f32;
            q = (0.5 * z as f64 * (v + v0) as f64) as f32
        } else {
            s = 0.0;
            c = -1.0;
            t = (-((v / v0) as f64).ln() / a as f64) as f32;
            q = (-0.5 * z as f64 * (v + v0) as f64) as f32
        }
    } else {
        let mx = -x;
        z0 = v0 / a;
        z += z0;
        x0 = (-(x * x + z * z - z0 * z0) as f64 / (2.0 * x as f64)) as f32;
        x += x0;
        let rr = ((x * x + z * z) as f64).sqrt();
        let r = if (a > 0.0 && mx > 0.0) || (a < 0.0 && mx < 0.0) {
            rr as f32
        } else {
            -rr as f32
        };
        let o = (1.0 / r as f64) as f32;
        let c0 = x0 * o;
        c = x * o;
        s = -z * o;
        t = (((v as f64 * (1.0 + c0 as f64)) / (v0 as f64 * (1.0 + c as f64))).ln() / a as f64)
            .abs() as f32;
        q = a * r * mx
    }
    if dvdx != 0.0 {
        let tmp = cr * s + sr * c;
        c = cr * c - sr * s;
        s = tmp
    }
    (c, s, t, q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflector_amplitude_skips_leading_colons() {
        for spec in [":0,1;2,1", "::350,300;550,300"] {
            assert_eq!(decodeReflector(spec), None, "{spec}");
        }
        assert_eq!(
            decodeReflector(":1:350,300;550,300"),
            Some(ReflectorInput {
                a: 1.0,
                x: vec![350.0, 550.0],
                z: vec![300.0, 300.0],
            })
        );
    }

    #[test]
    fn reflector_coordinates_walk_strtok_buffer() {
        for (spec, x, z) in [
            ("0.5,1;;2,1", [0.5, 0.0], [1.0, 1.0]),
            ("0,,1;2,1", [0.0, 2.0], [0.0, 1.0]),
            ("0,1;2,,1", [0.0, 2.0], [1.0, 0.0]),
            ("1:100,,200;900,200", [100.0, 900.0], [0.0, 200.0]),
            ("1:100,200;;900,200", [100.0, 0.0], [200.0, 200.0]),
            (",;0,1;2,1;", [0.0, 2.0], [1.0, 1.0]),
        ] {
            let r = decodeReflector(spec).unwrap();
            assert_eq!(r.a, 1.0, "{spec}");
            assert_eq!(r.x, x, "{spec}");
            assert_eq!(r.z, z, "{spec}");
        }
    }

    #[test]
    fn reflector_segment_count_adds_in_float() {
        let r = ReflectorInput {
            a: 1.0,
            x: vec![0.0, 0.99999994],
            z: vec![1.0, 1.0],
        };
        let out = makeref(1.0, &[r]);
        assert_eq!(out[0].ns, 2);
        assert_eq!(out[0].rs.len(), 2);
        assert_eq!(out[0].ds, 0.49999997);
    }

    #[test]
    fn negative_allocation_counts_exit_with_c_error() {
        const CHILD: &str = "OXWAVE_MODELING_ALLOCATION_TEST";
        if let Ok(mode) = std::env::var(CHILD) {
            oxwave_core::par::initargs(&["modeling-test".into()]);
            match mode.as_str() {
                "reflector" => {
                    makeref(-0.1, &[decodeReflector("0,1;1,1").unwrap()]);
                }
                "ricker" => {
                    makericker(5.0, -0.04);
                }
                _ => panic!("unknown allocation test"),
            }
            panic!("negative allocation count was accepted");
        }
        for mode in ["reflector", "ricker"] {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "tests::negative_allocation_counts_exit_with_c_error",
                    "--nocapture",
                ])
                .env(CHILD, mode)
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(1), "{mode}");
            assert_eq!(
                String::from_utf8(output.stderr).unwrap(),
                "\nmodeling-test: ealloc.c: malloc failed (Cannot allocate memory)\n",
                "{mode}"
            );
        }
    }
}
