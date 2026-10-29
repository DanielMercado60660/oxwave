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

//! Port of comp/dct/lib/lpred.c.
//! Credits: Tong Chen, 01/26/95; modified 05/25/95 for multiple plane waves.

/// C `lpFilter` storage. The original lpred does not modify its caller's filter.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LpFilter {
    pub hlfilter: usize,
    pub lfilter: usize,
    pub nslope: usize,
    pub nblock: usize,
    pub filter: Vec<Vec<Vec<f32>>>,
}

/// C prototype: `void lpred(float **f, int n1, int n2, int lblock, int hlfilter, int nslope, lpFilter *lpf, float **res);`
///
/// `n2` is inferred from `f.len()`. Supply complete blocks in both `f` and
/// `res`: C reads/writes through `ceil(n1/lblock)*lblock` even when n1 is
/// not a multiple of lblock. The extra input samples are read for RHS and
/// residuals, but the prediction matrix is zero-padded after n1, as in C.
/// C reassigns its local lpf pointer, leaving the caller's structure unchanged;
/// this port preserves that behavior and releases the private filter storage.
pub fn lpred(
    f: &[Vec<f32>],
    n1: usize,
    lblock: usize,
    hlfilter: usize,
    nslope: usize,
    _lpf: &mut LpFilter,
    res: &mut [Vec<f32>],
) {
    assert!(n1 > 0 && lblock > 0 && nslope <= f.len());
    let n2 = f.len();
    let nblock = (n1 - 1) / lblock + 1;
    let ltrace = nblock * lblock;
    let lpadded = ltrace + 2 * hlfilter;
    for i in 0..n2 {
        assert!(f[i].len() >= ltrace && res[i].len() >= ltrace);
    }
    let lfilter = 2 * hlfilter + 1;
    let mut filter = vec![vec![vec![0.0f32; lfilter]; nslope]; nblock];
    let nrow = (n2 - nslope) * lblock;
    let ncol = lfilter * nslope;
    let niter = ncol;
    let mut rhs = vec![0.0f32; nrow];
    let mut paddedtrace = vec![vec![0.0f32; lpadded]; n2];
    let mut xmatrics = vec![vec![0.0f32; ncol]; nrow];
    let mut tmp = vec![0.0f32; lblock];
    let mut unknown = vec![0.0f32; ncol];
    for itrace in 0..n2 {
        for it in 0..hlfilter {
            paddedtrace[itrace][it] = 0.0;
        }
        for i in 0..n1 {
            paddedtrace[itrace][hlfilter + i] = f[itrace][i];
        }
        for it in n1 + hlfilter..lpadded {
            paddedtrace[itrace][it] = 0.0;
        }
    }
    for iblock in 0..nblock {
        let ifirst = iblock * lblock;
        let mut xmax = 0.0f32;
        for itrace in 0..n2 - nslope {
            for islope in 0..nslope {
                for i in 0..lblock {
                    for ifilter in 0..lfilter {
                        let irow = itrace * lblock + i;
                        let icol = islope * lfilter + ifilter;
                        let j = itrace + islope;
                        let it = ifirst + i + ifilter;
                        xmatrics[irow][icol] = paddedtrace[j][it];
                        let v = xmatrics[irow][icol];
                        let xabs = if v < 0.0 { -v } else { v };
                        xmax = if xmax > xabs { xmax } else { xabs };
                    }
                }
            }
        }
        if xmax < f32::MIN_POSITIVE {
            for islope in 0..nslope {
                for ifilter in 0..lfilter {
                    filter[iblock][islope][ifilter] = 0.0;
                }
            }
            for itrace in 0..n2 {
                for i in 0..lblock {
                    let it = ifirst + i;
                    res[itrace][it] = f[itrace][it];
                }
            }
            continue;
        } else {
            let rscale = (1.0 / xmax as f64) as f32;
            for irow in 0..nrow {
                for icol in 0..ncol {
                    xmatrics[irow][icol] *= rscale;
                }
            }
            for itrace in nslope..n2 {
                for i in 0..lblock {
                    let it = ifirst + i;
                    let irow = (itrace - nslope) * lblock + i;
                    rhs[irow] = rscale * f[itrace][it];
                }
            }
            cgsolver(&xmatrics, &rhs, nrow, ncol, niter, 0.0001, &mut unknown);
        }
        for islope in 0..nslope {
            for ifilter in 0..lfilter {
                let icol = islope * lfilter + ifilter;
                filter[iblock][islope][ifilter] = unknown[icol];
            }
        }
        for islope in 0..nslope {
            for i in 0..lblock {
                let it = ifirst + i;
                res[islope][it] = f[islope][it];
            }
        }
        for itrace in nslope..n2 {
            let j = itrace - nslope;
            for i in 0..lblock {
                tmp[i] = 0.0;
                for islope in 0..nslope {
                    for ifilter in 0..lfilter {
                        let it = i + ifirst + ifilter;
                        tmp[i] += paddedtrace[j + islope][it] * filter[iblock][islope][ifilter];
                    }
                }
            }
            for i in 0..lblock {
                res[itrace][ifirst + i] = tmp[i];
            }
            for i in 0..lblock {
                let it = ifirst + i;
                res[itrace][it] = f[itrace][it] - res[itrace][it];
            }
        }
    }
}

fn cgsolver(
    a: &[Vec<f32>],
    rhs: &[f32],
    nrow: usize,
    ncol: usize,
    niter: usize,
    error: f32,
    x: &mut [f32],
) {
    let mut s = vec![0.0f32; nrow];
    let mut r = vec![0.0f32; ncol];
    let mut p = vec![0.0f32; ncol];
    let mut q = vec![0.0f32; nrow];
    let mut rn = vec![0.0f32; ncol];
    for i in 0..ncol {
        x[i] = 0.0;
    }
    for i in 0..nrow {
        s[i] = rhs[i];
    }
    for i in 0..ncol {
        p[i] = 0.0;
    }
    for j in 0..nrow {
        for i in 0..ncol {
            p[i] += a[j][i] * s[j];
        }
    }
    for i in 0..ncol {
        r[i] = p[i];
    }
    for i in 0..nrow {
        q[i] = 0.0;
        for j in 0..ncol {
            q[i] += a[i][j] * p[j];
        }
    }
    let mut rsum = 0.0f32;
    for i in 0..ncol {
        rsum += r[i] * r[i];
    }
    let r0sum = rsum * error;
    let mut icount = 0;
    while icount < niter && rsum > r0sum {
        let mut qsum = 0.0f32;
        for i in 0..nrow {
            qsum += q[i] * q[i];
        }
        let alpha = rsum / qsum;
        for i in 0..ncol {
            x[i] += alpha * p[i];
        }
        for i in 0..nrow {
            s[i] -= alpha * q[i];
        }
        for i in 0..ncol {
            rn[i] = 0.0;
        }
        for j in 0..nrow {
            for i in 0..ncol {
                rn[i] += a[j][i] * s[j];
            }
        }
        let mut rnsum = 0.0f32;
        for i in 0..ncol {
            rnsum += rn[i] * rn[i];
        }
        let beta = rnsum / rsum;
        for i in 0..ncol {
            p[i] = rn[i] + beta * p[i];
        }
        for i in 0..nrow {
            q[i] = 0.0;
            for j in 0..ncol {
                q[i] += a[i][j] * p[j];
            }
        }
        for i in 0..ncol {
            r[i] = rn[i];
        }
        rsum = rnsum;
        icount += 1;
    }
}
