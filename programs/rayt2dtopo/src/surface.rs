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

//! Surface helpers from par/main/ray_theory/rayt2dtopo.c.
//! Copyright (c) Colorado School of Mines, 1995. All rights reserved.
//! Trino Salinas, CSM, 1996; adapted from Dave Hale's modeling library,
//! decodeReflectors, decodeReflector and makeref, CSM, 09/17/91.

use oxwave_core::prelude::*;
use oxwave_cubic_interp::intcub;
use oxwave_f_cubicspline::csplin;
use oxwave_inverse_interp::yxtoxy;

#[derive(Clone, Copy, Default)]
pub(super) struct SurfaceSegment {
    x: f32,
    z: f32,
    s: f32,
    c: f32,
}
pub(super) struct Surface {
    ss: Vec<SurfaceSegment>,
}

pub(super) fn decode_surfaces() -> Vec<(Vec<f32>, Vec<f32>)> {
    let nr = countparname("surf").max(1);
    let mut result = Vec::new();
    for ir in 0..nr {
        let mut s = String::from("0,0;99999,0");
        getnparstring(ir + 1, "surf", &mut s);
        // strtok terminates each token, leaving skipped repeated delimiters
        // intact. The C conversion loop then walks the resulting buffer by
        // strlen+1 rather than following strtok's token positions.
        let mut buffer = s.into_bytes();
        buffer.push(0);
        let mut pos = 0;
        let mut first = None;
        let mut nvalues = 0;
        while pos < buffer.len() - 1 {
            while matches!(buffer[pos], b',' | b';') {
                pos += 1;
            }
            if buffer[pos] == 0 {
                break;
            }
            if first.is_none() {
                first = Some(pos);
            }
            nvalues += 1;
            while !matches!(buffer[pos], 0 | b',' | b';') {
                pos += 1;
            }
            if buffer[pos] == 0 {
                break;
            }
            buffer[pos] = 0;
            pos += 1;
        }
        if nvalues % 2 != 0 || nvalues < 4 {
            err!("Surface number %d specified incorrectly!\n", ir + 1);
        }
        let mut x = Vec::new();
        let mut z = Vec::new();
        let mut pos = first.unwrap();
        for _ixz in 0..nvalues / 2 {
            let end = pos + buffer[pos..].iter().position(|&b| b == 0).unwrap();
            x.push(atof(&String::from_utf8_lossy(&buffer[pos..end])) as f32);
            pos = end + 1;
            let end = pos + buffer[pos..].iter().position(|&b| b == 0).unwrap();
            z.push(atof(&String::from_utf8_lossy(&buffer[pos..end])) as f32);
            pos = end + 1;
        }
        result.push((x, z));
    }
    result
}

pub(super) fn makesurf(dsmax: f32, surfaces: &[(Vec<f32>, Vec<f32>)]) -> Vec<Surface> {
    let mut rr = Vec::new();
    for (xu, zu) in surfaces {
        let nu = xu.len();
        let mut u = vec![0.0f32; nu];
        for iu in 0..nu {
            u[iu] = iu as f32;
        }
        let mut xud = vec![[0.0f32; 4]; nu];
        csplin(&u, xu, &mut xud);
        let mut zud = vec![[0.0f32; 4]; nu];
        csplin(&u, zu, &mut zud);
        let nuu = 20 * nu;
        let duu = (u[nu - 1] - u[0]) / (nuu - 1) as f32;
        let mut s = vec![0.0f32; nuu];
        s[0] = 0.0;
        let mut xlast = xu[0];
        let mut zlast = zu[0];
        let mut uu = duu;
        for iuu in 1..nuu {
            let mut x = [0.0f32];
            let mut z = [0.0f32];
            intcub(0, &u, &xud, &[uu], &mut x);
            intcub(0, &u, &zud, &[uu], &mut z);
            let dx = x[0] - xlast;
            let dz = z[0] - zlast;
            s[iuu] = (s[iuu - 1] as f64 + ((dx * dx + dz * dz) as f64).sqrt()) as f32;
            xlast = x[0];
            zlast = z[0];
            uu += duu;
        }
        let ns = (1.0f32 + s[nuu - 1] / dsmax) as i32;
        let ds = s[nuu - 1] / ns as f32;
        let fs = (0.5 * ds as f64) as f32;
        let mut us = vec![0.0f32; ns as usize];
        yxtoxy(duu, 0.0, &s, ds, fs, 0.0, (nu - 1) as f32, &mut us);
        let mut ss = vec![SurfaceSegment::default(); ns as usize];
        for is in 0..ns as usize {
            let (mut rsx, mut rsz, mut rsxd, mut rszd) = ([0.0f32], [0.0f32], [0.0f32], [0.0f32]);
            intcub(0, &u, &xud, &[us[is]], &mut rsx);
            intcub(0, &u, &zud, &[us[is]], &mut rsz);
            intcub(1, &u, &xud, &[us[is]], &mut rsxd);
            intcub(1, &u, &zud, &[us[is]], &mut rszd);
            ss[is].x = rsx[0];
            ss[is].z = rsz[0];
            ss[is].c =
                (rsxd[0] as f64 / ((rsxd[0] * rsxd[0] + rszd[0] * rszd[0]) as f64).sqrt()) as f32;
            ss[is].s =
                (-rszd[0] as f64 / ((rsxd[0] * rsxd[0] + rszd[0] * rszd[0]) as f64).sqrt()) as f32;
        }
        rr.push(Surface { ss });
    }
    rr
}

pub(super) fn zcoor_topog(fxs: f32, dxs: f32, nxs: i32, srf: &[Surface]) -> (Vec<f32>, Vec<f32>) {
    let ss = &srf[0].ss;
    let ns = ss.len();
    let mut ssx = vec![0.0f32; ns];
    let mut ssz = vec![0.0f32; ns];
    let mut sss = vec![0.0f32; ns];
    for jx in 0..ns {
        ssx[jx] = ss[jx].x;
        ssz[jx] = ss[jx].z;
        sss[jx] = ss[jx].s;
    }
    // C extrapolates both endpoints even for nxs < 3, accessing uninitialized
    // or out-of-allocation floats. Zero padding models fresh zeroed heap pages
    // without performing those invalid accesses in Rust.
    let mut sz = vec![0.0f32; nxs.max(3) as usize + 6];
    let mut nangl = vec![0.0f32; nxs.max(3) as usize + 6];
    let mut ixi = 0;
    for ix in 0..nxs.max(0) as usize {
        let x = fxs + ix as f32 * dxs;
        for ik in ixi..ns.saturating_sub(1) {
            if ssx[ik] <= x && ssx[ik + 1] >= x {
                let mut az = ssz[ik] - ssz[ik + 1];
                let ax0 = ssx[ik + 1] - x;
                let ax = ssx[ik + 1] - ssx[ik];
                sz[ix + 3] = ax0 * az / ax + ssz[ik + 1];
                az = sss[ik] - sss[ik + 1];
                let temp = ax0 * az / ax + sss[ik + 1];
                nangl[ix + 3] = (temp as f64).asin() as f32;
                ixi = ik;
            }
        }
    }
    sz[3] = 2.0f32 * sz[4] - sz[5];
    nangl[3] = 2.0f32 * nangl[4] - nangl[5];
    let end = (nxs + 2) as usize;
    sz[end] = 2.0f32 * sz[end - 1] - sz[end - 2];
    nangl[end] = 2.0f32 * nangl[end - 1] - nangl[end - 2];
    for is in 0..nxs.max(0) as usize {
        warn!("sz[%d]=%g,nangl[%d]=%g", is, sz[is + 3], is, nangl[is + 3]);
    }
    (
        sz[3..3 + nxs.max(0) as usize].to_vec(),
        nangl[3..3 + nxs.max(0) as usize].to_vec(),
    )
}
