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

// Original source: Copyright (c) Colorado School of Mines, 2011.
// All rights reserved.
//! SURADON - generalized Radon transforms and multiple removal.
//! Port of su/main/transforms/suradon.c (Revision 1.21).
//! Credits: CWP: John Anderson (visitor to CSM from Mobil), Spring 1993.
//! References include Anderson CWP-137 (1993) and the SEG expanded abstracts
//! cited in the original source.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::*;
use oxwave_f_vnd::*;
use oxwave_linear_interp::intlin;
use oxwave_pfa_fft::{npfar, pfacr, pfarc};
use oxwave_sinc_interp::ints8r;
use std::f64::consts::PI;

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut xopt = 1;
    let mut cdpkey = String::new();
    let mut offkey = String::new();
    let (mut choose, mut ninterp, mut nk, mut igopt, mut nxmax) = (0, 0, 1, 1, 240);
    let (mut lagc, mut lent, mut lenx) = (400, 5, 7);
    getparint("xopt", &mut xopt);
    if getparstring("cdpkey", &mut cdpkey) == 0 {
        cdpkey = "cdp".into();
    }
    if getparstring("offkey", &mut offkey) == 0 {
        offkey = "offset".into();
    }
    getparint("choose", &mut choose);
    getparint("ninterp", &mut ninterp);
    getparint("nwin", &mut nk);
    getparint("igopt", &mut igopt);
    getparint("nxmax", &mut nxmax);
    // The second read in the C source overrides the documented default.
    if getparint("ninterp", &mut ninterp) == 0 {
        ninterp = 1;
    }
    getparint("lagc", &mut lagc);
    getparint("lent", &mut lent);
    getparint("lenx", &mut lenx);
    let (mut freq1, mut freq2, mut offref, mut f1, mut f2) =
        (4.0f32, 20.0f32, 2000.0f32, 60.0f32, 80.0f32);
    let (mut pmin, mut pmax, mut dp, mut pmula, mut pmulb) =
        (-200.0f32, 400.0f32, 16.0f32, 80.0f32, 200.0f32);
    let (mut depthref, mut intercept_off, mut prewhite) = (500.0f32, 0.0f32, 0.1f32);
    getparfloat("freq1", &mut freq1);
    getparfloat("freq2", &mut freq2);
    getparfloat("offref", &mut offref);
    getparfloat("f1", &mut f1);
    getparfloat("f2", &mut f2);
    getparfloat("pmin", &mut pmin);
    getparfloat("pmax", &mut pmax);
    getparfloat("dp", &mut dp);
    getparfloat("pmula", &mut pmula);
    getparfloat("pmulb", &mut pmulb);
    getparfloat("depthref", &mut depthref);
    getparfloat("interoff", &mut intercept_off);
    getparfloat("prewhite", &mut prewhite);
    let (mut ltaper, mut anderson) = (7, 1);
    getparint("ltaper", &mut ltaper);
    getparint("anderson", &mut anderson);
    checkpars();
    let cdpindex = getindex(&cdpkey);
    let offindex = getindex(&offkey);
    if nk < 0 {
        nk = 1;
    }
    let fac = (1000.0 * gofx(igopt, offref, intercept_off, depthref) as f64) as f32;
    pmin /= fac;
    pmax /= fac;
    dp /= fac;
    pmula /= fac;
    pmulb /= fac;
    let pminorig = pmin;
    let dporig = dp;
    let mut ipa = ((pmula - pmin) / dp) as i32;
    let ipb = ((pmulb - pmin) / dp) as i32;
    let mut np = (1.0 + (pmax - pmin) / dp) as i32;
    if np < 1 {
        err!("Range of PMIN and PMAX invalid");
    }
    if choose == 0 || choose == 3 || choose == 4 {
        ipa = 0;
    }
    if ipa < 0 {
        ipa = 0;
    }
    let mut tr = Trace::new();
    let mut tro = Trace::new();
    if gettr(&mut tr) == 0 {
        err!("Can't get first trace \n");
    }
    let nt = tr.ns as usize;
    let dt = (tr.dt as f64 / 1000000.0) as f32;
    let mut hdrwd = Value::default();
    gethval(&tr, cdpindex, &mut hdrwd);
    let mut oldcmp = hdrwd.i();
    let ntfft = npfar(nt as i32) as usize;
    // Invalid dimensions lead to out-of-allocation accesses in C. Reject them
    // before passing them to VND or converting signed sizes to usize.
    if nxmax < 1 || ninterp < -1 {
        err!(
            "invalid panel dimensions: nxmax=%d ninterp=%d",
            nxmax,
            ninterp
        );
    }
    let nxinterpmax = (1i64 + ninterp as i64) * (nxmax as i64 - 1) + 1;
    if nxinterpmax > i32::MAX as i64 {
        err!("interpolated trace count exceeds integer range");
    }
    let nxinterpmax = nxinterpmax as usize;
    let headerfile = tempname("radontmp");
    let headerfp = fopen(&headerfile, "w+")
        .unwrap_or_else(|| err!("couldn't open temp file for trace headers"));
    let mut orig = panel(1_000_000, 4, nt, nxmax as usize, "radontmp");
    let mut interp = panel(2_000_000, 4, nt, nxinterpmax, "radontmp");
    let mut result = panel(
        1_000_000,
        4,
        ntfft + 2,
        nxinterpmax.max(np as usize),
        "radontmp",
    );
    let mut rt = vec![0.0f32; ntfft + 2];
    let mut trace = vec![0.0f32; nt];
    let mut nx = 0;
    let mut offset = vec![0.0f32; nxmax as usize];
    let mut g = vec![0.0f32; nxmax as usize];
    // C indexes this with interpolated trace indices. Extending its mute
    // allocation with zeros removes the overread without changing the zero
    // gather's output (ub-2026-11, round 1).
    let mut mutetime = vec![0.0f32; (nxmax as usize).max(nxinterpmax)];
    let mut icount = 0;
    let mut ieod = false;
    loop {
        gethval(&tr, cdpindex, &mut hdrwd);
        if hdrwd.i() != oldcmp || ieod {
            if nx > 1 {
                let nxinterp = 1 + (nx - 1) * (ninterp + 1) as usize;
                if choose == 4 {
                    let mut k = 1;
                    for j in 1..nx {
                        if ((offset[j] - offset[j - 1]) as f64).abs() > 0.001 {
                            k += 1;
                        }
                    }
                    np = nx as i32;
                    nx = k;
                    for j in 0..np as usize {
                        read0(&mut orig, j, &mut rt, 1001);
                        rt[nt..].fill(0.0);
                        write0(&mut result, j, &mut rt, 1002);
                    }
                } else {
                    jea_xinterpolate(
                        &mut orig,
                        &mut interp,
                        ninterp,
                        nt,
                        nx,
                        freq1,
                        freq2,
                        lagc,
                        lent,
                        lenx,
                        xopt,
                        dt,
                    );
                    let d = (1.0 / (1.0 + ninterp as f64)) as f32;
                    let mut rt = vec![0.0; nxinterp];
                    let mut xin = vec![0.0; nx];
                    let mut gg = vec![0.0; nxinterp];
                    for j in 0..nxinterp {
                        rt[j] = j as f32 * d;
                    }
                    for j in 0..nx {
                        xin[j] = j as f32;
                    }
                    intlin(&xin, &offset[..nx], offset[0], offset[nx - 1], &rt, &mut gg);
                    for j in 0..nxinterp {
                        gg[j] = gofx(igopt, gg[j], intercept_off, depthref);
                    }
                    forward_p_transform(
                        &mut interp,
                        &mut result,
                        nxinterp,
                        nt,
                        &gg,
                        dt,
                        ntfft,
                        np as usize,
                        pmin,
                        dp,
                        &mutetime,
                        nk,
                        f1,
                        f2,
                        prewhite,
                    );
                }
                let mut nxout = np as usize;
                if choose >= 1 {
                    if choose == 1 {
                        for j in 0..ipb {
                            read0(&mut result, j as usize, &mut rt, 1003);
                            taupmute(j, ipa, ipb, ntfft, ltaper, &mut rt);
                            write0(&mut result, j as usize, &mut rt, 1004);
                        }
                    }
                    inverse_p_transform(
                        &mut result,
                        nx,
                        &g,
                        dt,
                        ntfft,
                        np as usize,
                        pmin,
                        dp,
                        ipa,
                        f1,
                        f2,
                    );
                    if choose == 1 {
                        for ix in 0..nx {
                            read0(&mut orig, ix, &mut trace, 1005);
                            read0(&mut result, ix, &mut rt, 1006);
                            for it in 0..nt {
                                rt[it] = trace[it] - rt[it];
                            }
                            write0(&mut result, ix, &mut rt, 1007);
                        }
                    }
                    nxout = nx;
                }
                erewind(&headerfp);
                for ix in 0..nxout {
                    if ix < nx {
                        efread(tro.hdr.as_bytes_mut(), &headerfp);
                    }
                    read0(&mut result, ix, &mut rt, 1008);
                    tro.data[..nt].copy_from_slice(&rt[..nt]);
                    icount += 1;
                    tro.tracl = icount;
                    tro.tracr = ix as i32 + 1;
                    if choose == 0 {
                        if anderson == 1 {
                            tro.f2 = (1000.0
                                * (pmin + ix as f32 * dp) as f64
                                * gofx(igopt, offref, intercept_off, depthref) as f64)
                                as f32;
                            tro.d2 = (1000.0
                                * dp as f64
                                * gofx(igopt, offref, intercept_off, depthref) as f64)
                                as f32;
                        } else if anderson == 0 {
                            tro.f2 = pminorig;
                            tro.d2 = dporig;
                        }
                    } else if choose == 4 {
                        tro.f2 = 0.0;
                        tro.d2 = 0.0;
                    }
                    puttr(&mut tro);
                }
            } else {
                erewind(&headerfp);
                efread(tro.hdr.as_bytes_mut(), &headerfp);
                read0(&mut orig, 0, &mut rt, 1009);
                tro.data[..nt].copy_from_slice(&rt[..nt]);
                icount += 1;
                tro.tracl = icount;
                tro.tracr = 1;
                puttr(&mut tro);
            }
            oldcmp = hdrwd.i();
            nx = 0;
            erewind(&headerfp);
            if ieod {
                break;
            }
        }
        if nx >= nxmax as usize {
            // C writes offset/g/mutetime past their allocations here.
            err!("gather exceeds nxmax=%d", nxmax);
        }
        efwrite(tr.hdr.as_bytes(), &headerfp);
        gethval(&tr, offindex, &mut hdrwd);
        offset[nx] = hdrwd.i() as f32;
        g[nx] = gofx(igopt, offset[nx], intercept_off, depthref);
        mutetime[nx] = tr.muts as f32;
        write0(&mut orig, nx, &mut tr.data[..nt], 1010);
        nx += 1;
        if gettr(&mut tr) == 0 {
            ieod = true;
        }
    }
    fclose(&headerfp);
    let _ = std::fs::remove_file(headerfile);
    VNDcl(orig, 1);
    VNDcl(interp, 1);
    VNDcl(result, 1);
    exit(CWP_Exit());
}

// These adapters only expose the native float/complex bytes expected by VND.
// Both types have no padding and accept every bit pattern; VND performs all
// panel storage, blocking, and I/O. No panel is duplicated in program memory.
trait VndValue {}
impl VndValue for f32 {}
impl VndValue for Complex {}
fn bytes<T: VndValue>(v: &mut [T]) -> &mut [u8] {
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr().cast(), std::mem::size_of_val(v)) }
}
fn read0<T: VndValue>(v: &mut VND, key: usize, out: &mut [T], msg: i32) {
    if key >= v.n[1] as usize {
        err!("panel trace index %d out of range", key);
    }
    V2Dr0(v, key as i64, bytes(out), msg);
}
fn write0<T: VndValue>(v: &mut VND, key: usize, data: &mut [T], msg: i32) {
    if key >= v.n[1] as usize {
        err!("panel trace index %d out of range", key);
    }
    V2Dw0(v, key as i64, bytes(data), msg);
}
fn read1<T: VndValue>(v: &mut VND, key: usize, out: &mut [T], msg: i32) {
    V2Dr1(v, key as i64, bytes(out), msg);
}
fn write1<T: VndValue>(v: &mut VND, key: usize, data: &mut [T], msg: i32) {
    V2Dw1(v, key as i64, bytes(data), msg);
}
fn tempname(root: &str) -> String {
    let name = VNDtempname(root);
    let filename = unsafe { std::ffi::CStr::from_ptr(name.cast()) }
        .to_string_lossy()
        .into_owned();
    unsafe {
        VNDfree(name, "suradon: fname");
    }
    filename
}
fn panel(memory: i64, nbytes: i64, n0: usize, n1: usize, root: &str) -> VND {
    V2Dop(2, memory, nbytes, &tempname(root), n0 as i64, n1 as i64)
        .unwrap_or_else(|| err!("couldn't open temp file for panel"))
}

fn gofx(igopt: i32, offset: f32, intercept_off: f32, refdepth: f32) -> f32 {
    let offset = offset - intercept_off;
    if igopt == 1 {
        return offset * offset;
    }
    if igopt == 2 {
        return ((refdepth * refdepth + offset * offset) as f64).sqrt() as f32;
    }
    if igopt == 4 {
        return (offset as f64).abs() as f32;
    }
    offset
}

fn freqweight(j: usize, df: f32, f1: f32, f2: f32) -> f32 {
    let f = j as f32 * df;
    if f <= f1 {
        return 1.0;
    }
    if f >= f2 {
        return 0.0;
    }
    (f2 - f) / (f2 - f1)
}

fn taupmute(ip: i32, ipa: i32, ipb: i32, nt: usize, ltap: i32, rt: &mut [f32]) {
    if ip >= ipb {
        return;
    }
    if ip <= ipa {
        for k in 0..nt {
            rt[k] = 0.0;
        }
        return;
    }
    let mut w = (ipb - ipa).max(1) as f32;
    w = (ipb - ip) as f32 / w;
    let j = (w * nt as f32) as usize;
    for k in 0..j {
        rt[k] = 0.0;
    }
    for k in j as i32..(nt as i32).min(j as i32 + ltap) {
        // Both operands of the taper in C are integers (arm64 division).
        rt[k as usize] *= ((k - j as i32).checked_div(ltap).unwrap_or(0)) as f32;
    }
}

fn forward_p_transform(
    vnda: &mut VND,
    vndb: &mut VND,
    nx: usize,
    nt: usize,
    g: &[f32],
    dt: f32,
    ntfft: usize,
    np: usize,
    pmin: f32,
    dp: f32,
    mutetime: &[f32],
    nk: i32,
    f1: f32,
    f2: f32,
    prewhite: f32,
) {
    let fac = (1.0 / ntfft as f64) as f32;
    let ntfftny = 1 + ntfft / 2;
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    let dw = (2.0 * PI * df as f64) as f32;
    let nmax = (vndb.n[0] as usize).max(vndb.n[1] as usize);
    // The C divides by nk before using it. arm64 integer division by zero is 0.
    let nxxinc = 1 + ((nx as i32 - 1).checked_div(nk).unwrap_or(0));
    let mut kindex = vec![0.0f32; nk as usize];
    let mut tindex = vec![0.0f32; nk as usize];
    if nk > 1 {
        for k in 0..nk as usize {
            kindex[k] = k as f32;
            let nxx = (nx as f64).min((k as f64 + 0.75) * nxxinc as f64) as usize;
            // C uses the interpolated trace count to index a mute allocation
            // of nxmax floats, without interpolating those times. Keep its
            // indices wherever allocated; an out-of-allocation read is UB.
            let Some(&mute) = mutetime.get(nxx) else {
                err!("mute-time index %d exceeds nxmax=%d", nxx, mutetime.len());
            };
            tindex[k] = (0.001 * mute as f64 / dt as f64) as f32;
        }
    }
    let nwindow = np
        .checked_mul(nk as usize)
        .filter(|&n| n <= i32::MAX as usize)
        .unwrap_or_else(|| err!("window trace count exceeds integer range"));
    let mut vndc = if nk > 1 {
        Some(panel(1_000_000, 8, ntfftny, nwindow, "radontmp"))
    } else {
        None
    };
    let mut rt = vec![0.0f32; 2 * nmax];
    let mut crt = vec![Complex::default(); nmax];
    let mut ccrt = vec![Complex::default(); ((nk as usize + 1) * np).max(vndb.n[1] as usize)];
    let mut r = vec![Complex::default(); np];
    let mut rhs = vec![Complex::default(); np];
    let mut wrk1 = vec![Complex::default(); np];
    let mut wrk2 = vec![Complex::default(); np];
    let mut wrk3 = vec![Complex::default(); np];
    let mut wrk4 = vec![Complex::default(); np];

    // Transform each input trace in place in the result VND panel.
    for ix in 0..nx {
        read0(vnda, ix, &mut rt, 201);
        for it in 0..nt {
            rt[it] *= fac;
        }
        for j in nt..ntfft {
            rt[j] = 0.0;
        }
        pfarc(1, &rt[..ntfft], &mut crt);
        write0(vndb, ix, &mut crt, 202);
    }
    VNDr2c(vndb);

    // Transform one frequency at a time for each offset window.
    for iw in 0..ntfftny {
        let wa = freqweight(iw, df, f1, f2);
        if wa > 0.0 {
            let w = iw as f32 * dw;
            read1(vndb, iw, &mut crt, 203);
            if wa < 1.0 {
                for ix in 0..nx {
                    crt[ix] = crmul(crt[ix], wa);
                }
            }
            for k in 0..nk as usize {
                let nxx = nx.min((k + 1) * nxxinc as usize);
                compute_rhs(w, nxx, g, &crt, np, pmin, dp, &mut rhs);
                compute_r(w, nxx, g, np, dp, &mut r);
                r[0].r = (r[0].r as f64 * (1.0 + prewhite as f64)) as f32;
                let mut rsum = 0.0f32;
                for j in 1..np {
                    rsum =
                        (rsum as f64 + ((r[j].r * r[j].r + r[j].i * r[j].i) as f64).sqrt()) as f32;
                }
                rsum /= r[0].r;
                if rsum as f64 > 1.0 + (np / 5) as f64 {
                    ctoephcg(
                        np / 7,
                        np,
                        &r,
                        &mut ccrt[k * np..(k + 1) * np],
                        &rhs,
                        &mut wrk1,
                        &mut wrk2,
                        &mut wrk3,
                        &mut wrk4,
                    );
                } else {
                    ctoep(
                        np,
                        &r,
                        &mut ccrt[k * np..(k + 1) * np],
                        &rhs,
                        &mut wrk1,
                        &mut wrk2,
                    );
                }
            }
        } else {
            for ip in 0..nwindow {
                ccrt[ip] = Complex::default();
            }
        }
        write1(vndc.as_mut().unwrap_or(&mut *vndb), iw, &mut ccrt, 204);
    }

    // Reuse the window panel for the inverse FFT; retain only one trace.
    let vc = vndc.as_mut().unwrap_or(&mut *vndb);
    for ip in 0..nwindow {
        read0(vc, ip, &mut crt, 205);
        pfacr(-1, &crt, &mut rt[..ntfft]);
        write0(vc, ip, &mut rt, 206);
    }
    VNDc2r(vc);

    // Merge windows by reading one time slice, rather than a full panel.
    if nk > 1 {
        VNDc2r(vndb);
        let vc = vndc.as_mut().unwrap();
        let mut rrt = vec![0.0f32; nwindow];
        let mut rk = [0.0f32; 1];
        for it in 0..ntfft {
            intlin(
                &tindex,
                &kindex,
                kindex[0],
                kindex[nk as usize - 1],
                &[it as f32],
                &mut rk,
            );
            if !rk[0].is_finite() || rk[0] < 0.0 || rk[0] >= nk as f32 {
                // Invalid mute coordinates can make C's float-to-index
                // conversion address memory outside the window panel.
                err!("invalid mute-window index");
            }
            let ik = rk[0] as usize;
            let ik2 = (ik + 1).min(nk as usize - 1);
            let wb = rk[0] - ik as f32;
            let wa = 1.0 - wb;
            read1(vc, it, &mut rrt, 207);
            // C writes the complete result row from scratch; only the first
            // np values are subsequently used.
            for ip in 0..np {
                rt[ip] = wa * rrt[ik * np + ip] + wb * rrt[ik2 * np + ip];
            }
            write1(vndb, it, &mut rt, 208);
        }
    }
    if let Some(vc) = vndc {
        VNDcl(vc, 1);
    }
}

fn inverse_p_transform(
    vnda: &mut VND,
    nx: usize,
    g: &[f32],
    dt: f32,
    ntfft: usize,
    np: usize,
    pmin: f32,
    dp: f32,
    ip1: i32,
    f1: f32,
    f2: f32,
) {
    let ntfftny = 1 + ntfft / 2;
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    let dw = (2.0 * PI * df as f64) as f32;
    let fac = (1.0 / ntfft as f64) as f32;
    let nmax = ntfftny.max(vnda.n[1] as usize).max(nx);
    let mut crt = vec![Complex::default(); nmax];
    let mut rt = vec![0.0f32; 2 * nmax];
    for ip in 0..np {
        read0(vnda, ip, &mut rt, 301);
        for it in 0..ntfft {
            rt[it] *= fac;
        }
        pfarc(1, &rt[..ntfft], &mut crt);
        write0(vnda, ip, &mut crt, 302);
    }
    VNDr2c(vnda);
    let fac = (1.0 / np as f64) as f32;
    let mut ctemp = vec![Complex::default(); np];
    for iw in 0..ntfftny {
        let wa = freqweight(iw, df, f1, f2);
        if wa > 0.0 {
            let w = iw as f32 * dw;
            read1(vnda, iw, &mut crt, 303);
            if wa < 1.0 {
                for ip in 0..np {
                    crt[ip] = crmul(crt[ip], wa);
                }
            }
            for ip in 0..np {
                ctemp[ip] = crt[ip];
            }
            for ix in 0..nx {
                let (mut rsum, mut isum) = (0.0f32, 0.0f32);
                for ip in ip1 as usize..np {
                    let p = pmin + ip as f32 * dp;
                    let tr = ((w * p * g[ix]) as f64).cos() as f32;
                    let ti = ((w * p * g[ix]) as f64).sin() as f32;
                    let dr = ctemp[ip].r;
                    let di = ctemp[ip].i;
                    rsum += tr * dr - ti * di;
                    isum += tr * di + ti * dr;
                }
                crt[ix].r = fac * rsum;
                crt[ix].i = fac * isum;
            }
        } else {
            for ix in 0..nx {
                crt[ix] = Complex::default();
            }
        }
        write1(vnda, iw, &mut crt, 304);
    }
    for ix in 0..nx {
        read0(vnda, ix, &mut crt, 305);
        pfacr(-1, &crt, &mut rt[..ntfft]);
        write0(vnda, ix, &mut rt, 306);
    }
    VNDc2r(vnda);
}

fn compute_r(w: f32, nx: usize, g: &[f32], np: usize, dp: f32, r: &mut [Complex]) {
    let fac = (1.0 / (nx * np) as f64) as f32;
    for j in 0..np {
        let (mut rsum, mut isum) = (0.0f32, 0.0f32);
        for k in 0..nx {
            rsum = (rsum as f64 + ((w * j as f32 * dp * g[k]) as f64).cos()) as f32;
            isum = (isum as f64 + ((w * j as f32 * dp * g[k]) as f64).sin()) as f32;
        }
        r[j].r = fac * rsum;
        r[j].i = fac * isum;
    }
}

fn compute_rhs(
    w: f32,
    nx: usize,
    g: &[f32],
    data: &[Complex],
    np: usize,
    pmin: f32,
    dp: f32,
    rhs: &mut [Complex],
) {
    let fac = (1.0 / nx as f64) as f32;
    for ip in 0..np {
        let p = pmin + ip as f32 * dp;
        let (mut rsum, mut isum) = (0.0f32, 0.0f32);
        for ix in 0..nx {
            let tr = ((w * p * g[ix]) as f64).cos() as f32;
            let ti = -((w * p * g[ix]) as f64).sin() as f32;
            let dr = data[ix].r;
            let di = data[ix].i;
            rsum += tr * dr - ti * di;
            isum += tr * di + ti * dr;
        }
        rhs[ip].r = fac * rsum;
        rhs[ip].i = fac * isum;
    }
}

// Complex Hermitian Toeplitz solver; these helpers are static in suradon.c.
fn ctoep(
    n: usize,
    r: &[Complex],
    a: &mut [Complex],
    b: &[Complex],
    f: &mut [Complex],
    g: &mut [Complex],
) -> usize {
    if r[0].r == 0.0 {
        return 0;
    }
    f[0].r = (1.0 / r[0].r as f64) as f32;
    f[0].i = 0.0;
    a[0].r = b[0].r / r[0].r;
    a[0].i = b[0].i / r[0].r;
    let (mut vr, mut vi, mut vsq) = (1.0f32, 0.0f32, 1.0f32);
    let mut j = 1;
    while j < n {
        f[j].r = 0.0;
        f[j].i = 0.0;
        let (mut er, mut ei) = (0.0f32, 0.0f32);
        for k in 0..j {
            let jmk = j - k;
            er += r[jmk].r * f[k].r + r[jmk].i * f[k].i;
            ei += r[jmk].r * f[k].i - r[jmk].i * f[k].r;
        }
        let cr = (er * vr - ei * vi) / vsq;
        let ci = (er * vi + ei * vr) / vsq;
        vr = vr - (cr * er + ci * ei);
        vi = vi + (cr * ei - ci * er);
        vsq = vr * vr + vi * vi;
        if vsq <= 0.0 {
            break;
        }
        for k in 0..=j {
            let jmk = j - k;
            g[k].r = f[k].r - cr * f[jmk].r - ci * f[jmk].i;
            g[k].i = f[k].i + cr * f[jmk].i - ci * f[jmk].r;
        }
        for k in 0..=j {
            f[k] = g[k];
        }
        a[j].r = 0.0;
        a[j].i = 0.0;
        er = 0.0;
        ei = 0.0;
        for k in 0..j {
            let jmk = j - k;
            er += r[jmk].r * a[k].r + r[jmk].i * a[k].i;
            ei += r[jmk].r * a[k].i - r[jmk].i * a[k].r;
        }
        er = er - b[j].r;
        ei = ei - b[j].i;
        let cr = (er * vr - ei * vi) / vsq;
        let ci = (er * vi + ei * vr) / vsq;
        for k in 0..=j {
            let jmk = j - k;
            a[k].r += -cr * f[jmk].r - ci * f[jmk].i;
            a[k].i += cr * f[jmk].i - ci * f[jmk].r;
        }
        j += 1;
    }
    vr = 0.0;
    vi = 0.0;
    for k in 0..j {
        vr += r[k].r * f[k].r - r[k].i * f[k].i;
        vi += r[k].r * f[k].i + r[k].i * f[k].r;
    }
    vsq = vr * vr + vi * vi;
    let er = vr / vsq;
    let ei = -vi / vsq;
    for k in 0..j {
        f[k].r = er * f[k].r - ei * f[k].i;
        f[k].i = er * f[k].i + ei * f[k].r;
    }
    j
}

fn ctoephcg(
    niter: usize,
    n: usize,
    a: &[Complex],
    x: &mut [Complex],
    y: &[Complex],
    s: &mut [Complex],
    ss: &mut [Complex],
    g: &mut [Complex],
    rr: &mut [Complex],
) -> i32 {
    let eps = 1.0e-6f32;
    let mut rp = rcdot(n, y, y);
    let test = n as f32 * eps * eps * rp;
    for j in 0..n {
        x[j] = Complex::default();
        rr[j] = y[j];
    }
    htmul(n, a, rr, g);
    for j in 0..n {
        s[j] = g[j];
    }
    let mut gammam = rcdot(n, g, g);
    let mut iter = 0;
    while iter < niter {
        htmul(n, a, s, ss);
        let alpha = gammam / rcdot(n, ss, ss);
        for j in 0..n {
            x[j] = cadd(x[j], crmul(s[j], alpha));
            rr[j] = csub(rr[j], crmul(ss[j], alpha));
        }
        let rsq = rcdot(n, rr, rr);
        if iter > 0 && (rsq == rp || rsq < test) {
            return iter as i32 - 1;
        }
        rp = rsq;
        htmul(n, a, rr, g);
        let gamma = rcdot(n, g, g);
        if gamma < eps {
            break;
        }
        let beta = gamma / gammam;
        gammam = gamma;
        for j in 0..n {
            s[j] = cadd(g[j], crmul(s[j], beta));
        }
        iter += 1;
    }
    iter as i32
}

fn rcdot(n: usize, a: &[Complex], b: &[Complex]) -> f32 {
    let mut sum = 0.0f32;
    for j in 0..n {
        sum += a[j].r * b[j].r + a[j].i * b[j].i;
    }
    sum
}

fn htmul(n: usize, a: &[Complex], x: &[Complex], y: &mut [Complex]) {
    for irow in 0..n {
        y[irow] = Complex::default();
        for j in 0..irow {
            y[irow] = cadd(cmul(conjg(a[irow - j]), x[j]), y[irow]);
        }
        for j in irow..n {
            y[irow] = cadd(cmul(a[j - irow], x[j]), y[irow]);
        }
    }
}

// jea_xinterpolate is only called with iopt=0 by this program.
fn jea_xinterpolate(
    vndorig: &mut VND,
    vndinterp: &mut VND,
    ninterp: i32,
    nt: usize,
    nx: usize,
    freq1: f32,
    freq2: f32,
    lagc: i32,
    lent: i32,
    lenx: i32,
    mut xopt: i32,
    dt: f32,
) {
    let lent = 1 + 2 * (lent / 2);
    let lenx = 1 + 2 * (lenx / 2);
    let lagc = (1.0 + lagc as f64 * 0.001 / dt as f64) as i32;
    let ntfft = npfar(nt as i32) as usize;
    let ntfftny = 1 + ntfft / 2;
    let nxfft = npfar(nx as i32) as usize;
    let nxfftny = 1 + nxfft / 2;
    let df = (1.0 / (ntfft as f32 * dt) as f64) as f32;
    if nx < 2 || ninterp == 0 {
        let mut rt = vec![0.0f32; nt];
        for ix in 0..nx {
            read0(vndorig, ix, &mut rt, 101);
            write0(vndinterp, ix, &mut rt, 102);
        }
        return;
    }
    let mut crt = vec![Complex::default(); ntfftny.max(nxfftny)];
    let mut ccrt = vec![Complex::default(); ntfftny];
    let mut rt = vec![0.0f32; ntfft.max(nxfft) + 2];
    let mut rrt = vec![0.0f32; ntfft];
    let mut a = vec![0.0f32; nx.max(nt)];
    let mut b = vec![0.0f32; nx.max(nt)];
    let mut p = vec![0.0f32; nt];
    let mut time = vec![0.0f32; nt];
    let mut aa = vec![0.0f32; nx.max(nt)];
    let mut bb = vec![0.0f32; nx.max(nt)];
    let mut vnda = panel(500_000, 4, nt, nx, "jea_xinterpolate");
    let mut vndb = panel(500_000, 4, nt, nx, "jea_xinterpolate");
    let mut eps = 0.0f32;
    for ix in 0..nx {
        read0(vndorig, ix, &mut rt, 203);
        for j in 0..nt {
            let val = (rt[j] as f64).abs();
            eps = if eps as f64 > val { eps } else { val as f32 };
        }
    }
    eps = (if 1e-2 * eps as f64 > 1e-30 {
        1e-2 * eps as f64
    } else {
        1e-30
    }) as f32;
    let dff = (2.0 * PI / ntfft as f64) as f32;
    for ix in 0..nx {
        read0(vndorig, ix, &mut rt, 103);
        for j in 0..nt {
            a[j] = (rt[j] as f64).abs() as f32;
        }
        runav(nt, lagc, &a, &mut b);
        runav(nt, lagc, &b, &mut a);
        for j in 0..nt {
            rt[j] = rt[j] / (a[j] + eps);
        }
        for j in nt..ntfft {
            rt[j] = 0.0;
        }
        pfarc(1, &rt[..ntfft], &mut crt);
        if freq1 > 0.0 {
            for j in 0..ntfftny {
                let f = j as f32 * df;
                let mut fcl = f / freq1;
                fcl = fcl * fcl * fcl * fcl;
                let mut fch = f / freq2;
                fch = fch * fch * fch * fch;
                let f = (fcl as f64 / ((1.0 + fcl as f64) * (1.0 + fch as f64))) as f32;
                crt[j] = crmul(crt[j], f);
                ccrt[j] = cmul(crt[j], cmplx(0.0, -(j as f32) * dff));
            }
        } else {
            for j in 0..ntfftny {
                let f = j as f32 * df;
                let fch = f / freq2;
                let f = (1.0 / (1.0 + (fch * fch * fch * fch) as f64)) as f32;
                crt[j] = crmul(crt[j], f);
                ccrt[j] = cmul(crt[j], cmplx(0.0, -(j as f32) * dff));
            }
        }
        pfacr(-1, &crt, &mut rt[..ntfft]);
        write0(&mut vnda, ix, &mut rt, 104);
        pfacr(-1, &ccrt, &mut rrt);
        write0(&mut vndb, ix, &mut rrt, 105);
    }
    let dxx = (2.0 * PI / (nxfft * nxfft) as f64) as f32;
    if nx < 16 {
        xopt = 1;
    }
    for it in 0..nt {
        read1(&mut vnda, it, &mut rt, 106);
        if xopt != 0 {
            for j in 0..nx - 1 {
                rt[j] = rt[j + 1] - rt[j];
            }
            rt[nx - 1] = rt[nx - 2];
        } else {
            for j in nx..nxfft {
                rt[j] = 0.0;
            }
            pfarc(1, &rt[..nxfft], &mut crt);
            for j in 0..nxfftny {
                crt[j] = cmul(crt[j], cmplx(0.0, -(j as f32) * dxx));
            }
            pfacr(-1, &crt, &mut rt[..nxfft]);
        }
        write1(&mut vnda, it, &mut rt, 107);
    }
    for ix in 0..nx {
        read0(&mut vnda, ix, &mut a, 108);
        read0(&mut vndb, ix, &mut b, 109);
        for it in 0..nt {
            aa[it] = a[it] * b[it];
            bb[it] = b[it] * b[it];
        }
        runav(nt, lent, &aa, &mut a);
        runav(nt, lent, &a, &mut aa);
        runav(nt, lent, &bb, &mut b);
        runav(nt, lent, &b, &mut bb);
        write0(&mut vnda, ix, &mut aa, 110);
        write0(&mut vndb, ix, &mut bb, 111);
    }
    if lenx > 1 {
        for it in 0..nt {
            read1(&mut vnda, it, &mut a, 112);
            read1(&mut vndb, it, &mut b, 113);
            runav(nx, lenx, &a, &mut aa);
            runav(nx, lenx, &aa, &mut a);
            runav(nx, lenx, &b, &mut bb);
            runav(nx, lenx, &bb, &mut b);
            write1(&mut vnda, it, &mut a, 114);
            write1(&mut vndb, it, &mut b, 115);
        }
    }
    eps = 0.0;
    for ix in 0..nx {
        read0(&mut vndb, ix, &mut bb, 218);
        for it in 0..nt {
            eps = if eps > bb[it] { eps } else { bb[it] };
        }
    }
    eps = (if 1e-5 * eps as f64 > 1e-30 {
        1e-5 * eps as f64
    } else {
        1e-30
    }) as f32;
    read0(vndorig, 0, &mut a, 116);
    for ix in 1..nx {
        let ixm = ix - 1;
        read0(&mut vnda, ixm, &mut aa, 117);
        read0(&mut vndb, ixm, &mut bb, 118);
        for it in 0..nt {
            p[it] = -aa[it] / (bb[it] + eps);
        }
        read0(vndorig, ix, &mut b, 119);
        write0(vndinterp, ixm * (ninterp + 1) as usize, &mut a, 120);
        for k in 0..ninterp {
            let wa = ((1.0 + k as f64) / (1 + ninterp) as f64) as f32;
            let wb = (1.0 - wa as f64) as f32;
            for it in 0..nt {
                time[it] = it as f32 - p[it] * wa;
            }
            zero_extrapolation_coordinates(&mut time);
            ints8r(1.0, 0.0, &a[..nt], 0.0, 0.0, &time, &mut aa[..nt]);
            for it in 0..nt {
                time[it] = it as f32 + p[it] * wb;
            }
            zero_extrapolation_coordinates(&mut time);
            ints8r(1.0, 0.0, &b[..nt], 0.0, 0.0, &time, &mut bb[..nt]);
            for it in 0..nt {
                aa[it] = wb * aa[it] + wa * bb[it];
            }
            write0(
                vndinterp,
                k as usize + 1 + ixm * (ninterp + 1) as usize,
                &mut aa,
                121,
            );
        }
        std::mem::swap(&mut a, &mut b);
    }
    write0(vndinterp, (nx - 1) * (ninterp + 1) as usize, &mut a, 122);
    VNDcl(vnda, 1);
    VNDcl(vndb, 1);
}

fn zero_extrapolation_coordinates(time: &mut [f32]) {
    // These calls to ints8r use unit spacing, zero origin and zero tails.
    // inttable8.c:177 cannot convert these extreme coordinates to int.
    // The minimally fixed C returns the zero tail directly, byte-identical
    // to the original for the lowpass-picker case (ub-2026-11, round 1).
    // Use a representable coordinate with all eight taps in that zero tail
    // so the shared interpolator can compute the same result safely.
    for x in time {
        let xoutn = (8.0f32 + *x) as f64;
        if !(xoutn < 2147483648.0 && xoutn >= -2147483637.0) {
            *x = -16.0;
        }
    }
}

fn runav(n: usize, len: i32, a: &[f32], b: &mut [f32]) {
    let mut sum = 0.0f32;
    let lenh = len / 2;
    if len <= 1 {
        for j in 0..n {
            b[j] = a[j];
        }
        return;
    }
    for j in 0..(len as usize).min(n) {
        sum += a[j];
    }
    for j in 0..((lenh + 1) as usize).min(n) {
        b[j] = sum;
    }
    for j in lenh + 1..n as i32 - lenh {
        sum = sum + a[(j + lenh) as usize] - a[(j - lenh - 1) as usize];
        b[j as usize] = sum;
    }
    for j in 0.max(n as i32 - lenh)..n as i32 {
        b[j as usize] = sum;
    }
    sum = (1.0 / len as f64) as f32;
    for j in 0..n {
        b[j] *= sum;
    }
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SURADON - compute forward or reverse Radon transform or remove multiples",
    "           by using the parabolic Radon transform to estimate multiples",
    "           and subtract.\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "     suradon <stdin >stdout [Optional Parameters]\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " choose=0    0  Forward Radon transform\t\t\t\t",
    "             1  Compute data minus multiples\t\t\t\t",
    "             2  Compute estimate of multiples\t\t\t\t",
    "             3  Compute forward and reverse transform\t\t\t",
    "             4  Compute inverse Radon transform\t\t\t",
    " igopt=1     1  parabolic transform: g(x) = offset**2\t\t\t",
    "             2  Foster/Mosher psuedo hyperbolic transform\t\t",
    "                   g(x) = sqrt(depth**2 + offset**2)\t\t\t",
    "             3  Linear tau-p: g(x) = offset\t\t\t\t",
    "             4  abs linear tau-p: g(x) = abs(offset)\t\t\t",
    " offref=2000.    reference maximum offset to which maximum and minimum\t",
    "                 moveout times are associated\t\t\t\t",
    " interoff=0.     intercept offset to which tau-p times are associated\t",
    " pmin=-200       minimum moveout in ms on reference offset\t\t",
    " pmax=400        maximum moveout in ms on reference offset\t\t",
    " anderson=1\t  p values in ms, =0 pvalues in horizontal slowness\t",
    " dp=16           moveout increment in ms on reference offset\t\t",
    " pmula=80        moveout in ms on reference offset where multiples begin",
    "                     at maximum time\t\t\t\t\t",
    " pmulb=200       moveout in ms on reference offset where multiples begin",
    "                     at zero time\t\t\t\t\t",
    " depthref=500.   Reference depth for Foster/Mosher hyperbolic transform",
    " nwin=1          number of windows to use through the mute zone\t",
    " f1=60.          High-end frequency before taper off\t\t\t",
    " f2=80.          High-end frequency\t\t\t\t\t",
    " prewhite=0.1    Prewhitening factor in percent.\t\t\t",
    " cdpkey=cdp      name of header word for defining ensemble\t\t",
    " offkey=offset   name of header word with spatial information\t\t",
    " nxmax=240       maximum number of input traces per ensemble\t\t",
    " ltaper=7\t  taper (integer) for mute tapering function\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optimizing Parameters:\t\t\t\t\t\t",
    " The following parameters are occasionally used to avoid spatial aliasing",
    " problems on the linear tau-p transform.  Not recommended for other\t",
    " transforms...\t\t\t\t\t\t\t\t",
    " ninterp=0      number of traces to interpolate between each input trace",
    "                   prior to computing transform\t\t\t",
    " freq1=4.0      low-end frequency in Hz for picking (good default: 3 Hz)",
    "                (Known bug: freq1 cannot be zero) ",
    " freq2=20.0     high-end frequency in Hz for picking (good default: 20 Hz)",
    " lagc=400       length of AGC operator for picking (good default: 400 ms)",
    " lent=5         length of time smoother in samples for picker\t\t",
    "                     (good default: 5 samples)\t\t\t\t",
    " lenx=7         length of space smoother in samples for picker\t\t",
    "                     (good default: 1 sample)\t\t\t\t",
    " xopt=1         1 = use differences for spatial derivative\t\t",
    "                        (works with irregular spacing)\t\t\t",
    "                0 = use FFT derivative for spatial derivatives\t\t",
    "                      (more accurate but requires regular spacing and\t",
    "                      at least 16 input tracs--will switch to differences",
    "                      automatically if have less than 16 input traces)\t",
    "\t\t\t\t\t\t\t\t\t",
];
