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

//! SUSEMGROUP, port of su/main/noise/susemgroup.c, revision 1.00.
//! Copyright (c) Colorado School of Mines, 2022. All rights reserved.
//! Author: Andre Latour, Feb 2023. Derived from susembase.
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                            ",
    " SUSEMGROUP - Semblance Based Structure-Oriented Smooth In Trace Groups.    ",
    "                                                                            ",
    " This program performs semblance scanning to derive (apparent) dips in the  ",
    " direction of offsets within cdps (or other key values). It then uses those ",
    " dips to align trace sample values and performs a normalized sum (which     ",
    " reduces noise). Since the dip compensation produces locally-flat horizons  ",
    " during summing, the result is less smearing of those horizons.             ",
    "      This program does not interpolate. No traces will be created.         ",
    "      This program can process 2D and 3D stacked data. Since interpolation  ",
    "      is often desirable for stacked data, see program susembase for that.  ",
    "                                                                            ",
    " susemgroup   <in.su >out.su  (other parameters)                            ",
    "                                                                            ",
    " Parameter overview:                                                        ",
    "                                                                            ",
    " lockey=offset Trace key that indicates location within groupkey group.     ",
    "               (Traces within groups do not have to be sorted by this key). ",
    "         Note: Other parameters default to values corresponding to this     ",
    "               lockey being offset. For 2D stacked data, you should SPECIFY ",
    "               lockkey=cdp and groupkey=none and approximately scanstep=0.5 ",
    "               and scanmax=12 and locmaxscan=3. For 3D stacked data, you    ",
    "               should SPECIFY lockkey=igi or igc and groupkey=igc or igi and",
    "               approximately scanstep=0.5 and scanmax=12 and locmaxscan=3.  ",
    "                                                                            ",
    " groupkey=cdp  Trace key that indicates a group of traces. Traces with the  ",
    "               same groupkey value are processed as a group. Any CHANGE of  ",
    "               this groupkey value indicates a new group. For cdps this     ",
    "               grouping is typically achieved by sorting. But, for instance,",
    "               a 3D shot gather typically contains traces in receiver line  ",
    "               groups, with NEW groups being defined herein for subsequent  ",
    "               shots even though the receiver line numbers REPEAT.          ",
    "         =none There is only 1 group and all lockey values are in the group.",
    "               Typically this option is for stacked input. For instance,    ",
    "               to input a 2D stack specify lockey=cdp (and groupkey=none).  ",
    "               For a 3D stack then lockey,groupkey are igi,igc then igc,igi ",
    "               (or vice-versa, depending on which direction you sort first).",
    "               For large stacked data, also increase the maxfold value.     ",
    "         Note: If you need to use susort, remember that the sort keys are   ",
    "               listed in reverse order of lockey and groupkey here.         ",
    "                                                                            ",
    " maxfold=999  This is maximum number of traces in a group. For instance,    ",
    "              if groupkey=cdp this is maximum number of traces in any cdp.  ",
    "                                                                            ",
    " locmaxscan=75 Maximum Relative lockey to include in scanning for dips.     ",
    "               Must be greater than or equal to 0. This program scans traces",
    "               from -locmaxscan to locmaxscan (inclusive) at each location. ",
    "           =0  Do not scan for dips. You must specify parameter locmaxsum   ",
    "               and you cannot specify scanstep, scanmax, scanmin, dipext.   ",
    "               Note: The default value corresponds to a lockey=offset,      ",
    "                     so the 75 typically means 75 metres or 75 feet.        ",
    "                                                                            ",
    " scanstep=0.05      Time Step to Scan (milliseconds per 1 lockey unit).     ",
    " scanmax =1.1       Maximum Time to Scan (milliseconds per 1 lockey unit).  ",
    " scanmin =-scanmax  Minimum Time to Scan (milliseconds per 1 lockey unit).  ",
    "              Note: These defaults correspond to a lockey=offset, so these  ",
    "                    values are typically milliseconds per metre or foot.    ",
    "                                                                            ",
    " dipext=12. Dip Extender (milliseconds).                                    ",
    "            At each location the semblance scanning produces an estimate    ",
    "            of the dip at each sample time. This parameter performs a       ",
    "            weighted-average of those estimated dips symmetrically from     ",
    "            -dipext to +dipext. The weights used are the sum of the         ",
    "            squared amplitudes of seismic values along the estimated dip    ",
    "            in the locmaxscan range of input traces. In cruder terms: the   ",
    "            power along the dip value chozen by the semblance scanning.     ",
    "            Effectively, this causes the dip values from large horizons     ",
    "            to be trusted more and extended up and down nearby times.       ",
    "                                                                            ",
    " dipuse=1  Estimate and use dips at input locations. (Same options as the   ",
    "           susembase program, but since input and output locations are same ",
    "           herein, the actual decision is how you want to use the dips).    ",
    "       =0  Estimate and use dips at output locations.                       ",
    "                                                                            ",
    " locmaxsum=locmaxscan  Maximum Relative lockey to include in normalized sum.",
    "                       Must be greater than 0. This program sums the input  ",
    "                       traces from -locmaxsum to locmaxsum (inclusive) at   ",
    "                       each location (and divides by the total number).     ",
    "                 Note: Remember that the sample values of these traces are  ",
    "                       corrected by their estimated dips before summing, so ",
    "                       the result is less smeared than uncorrected values.  ",
    "                                                                            ",
    " outdip=   Do not output this file.                                         ",
    "       =(file name) Output the dip values chosen by the semblance scanning  ",
    "           into trace samples in this file. The trace header values are     ",
    "           copied from the seismic traces.                                  ",
    "     Note: These dips are after dipext processing. But even with dipext=0   ",
    "           they will not have scanstep granularity because the scanned      ",
    "           maximum semblance dip undergoes a 3-point interpolation to get   ",
    "           a better approximation of the dip at the true peak in semblance. ",
    "     Note: It is technically more correct to call the dips made and used    ",
    "           by this program APPARANT dips. The dip values are put into the   ",
    "           trace samples and are in milliseconds per 1 location unit.       ",
    "                                                                            ",
    " outpow=   Do not output this file.                                         ",
    "       =(file name) Output the average power along the dip direction chosen ",
    "           by the semblance scanning. (Same header values as for outdip).   ",
    "                                                                            ",
    " outsem=   Do not output this file.                                         ",
    "       =(file name) Output the semblance value computed and chosen by the   ",
    "           semblance scanning. (Same header values as for outdip).          ",
    "                                                                            ",
    " indip=    Name of trace file which contains (apparant) dip values.         ",
    "           Dip scanning is disabled if you use this parameter. Thus you     ",
    "           cannot specify parameters scanstep, scanmax, scanmin, dipext,    ",
    "           and locmaxscan - and you must specify parameter locmaxsum.       ",
    "           This file is usually initially made by specifying outdip above.  ",
    "     Note: The intention here is to allow you to output dips from one run   ",
    "           of this program and use them in a second run of this program.    ",
    "           The first run can input smoothed, denoised seismic to get dips,  ",
    "           and the second run can input the un-smoothed, noisy seismic and  ",
    "           use the dips from the first run. As well, you can also actually  ",
    "           smooth and denoise the DIP VALUES in this file before input here.",
    "                                                                            ",
    " Note: This program enforces non-hierarchical parameter error checking.     ",
    "       This means if you disable a main option, all parameters related      ",
    "       to it must be removed from the command line. That is, for novice     ",
    "       users this program is sufficiently complex that it refuses to guess  ",
    "       which parameter was specified in error.                              ",
    "                                                                            ",
    " ------------------------------------------------------------------------   ",
    "                                                                            ",
];
use oxwave_core::prelude::*;
use oxwave_header_case::{fromhead, GetCase};
use oxwave_sinc_interp::ints8r;

fn scanit(
    t: &[Trace],
    l: &[i32],
    at: f64,
    range: f64,
    smin: f32,
    step: f32,
    ns: usize,
    nt: usize,
) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let mut dip = vec![0.; nt];
    let mut pow = vec![0.; nt];
    let mut sem = vec![0.; nt];
    let a: Vec<usize> = (0..t.len())
        .filter(|&i| (l[i] as f64) >= at - range && (l[i] as f64) <= at + range)
        .collect();
    let count = a.len();
    if count == 0 {
        return (dip, pow, sem);
    }
    if count == 1 {
        for n in 0..nt {
            pow[n] = t[a[0]].data[n] * t[a[0]].data[n];
            sem[n] = 1.;
        }
        return (dip, pow, sem);
    }
    let mut ii = vec![vec![0i32; ns]; count];
    let mut ff = vec![vec![0f32; ns]; count];
    for j in 0..ns {
        let s = smin + j as f32 * step;
        for i in 0..count {
            let x = s as f64 * (l[a[i]] as f64 - at);
            let mut m = x as i32;
            if x > 0. {
                m += 1
            }
            ii[i][j] = m;
            ff[i][j] = (m as f64 - x) as f32;
        }
    }
    let mut w = vec![0f32; ns];
    let mut q = vec![0f32; ns];
    for n in 0..nt {
        for j in 0..ns {
            w[j] = 0.;
            q[j] = 1e-30;
            for i in 0..count {
                let k = n as i32 + ii[i][j];
                if k > 0 && k < nt as i32 {
                    let f = ff[i][j];
                    let d = &t[a[i]].data;
                    let v = ((f * d[k as usize - 1]) as f64
                        + (1. - f as f64) * d[k as usize] as f64)
                        as f32;
                    w[j] += v;
                    q[j] += v * v;
                }
            }
        }
        let mut best = 0f32;
        let mut b = ns / 2;
        for j in 0..ns {
            let x = w[j] * w[j] / q[j];
            if x > best {
                best = x;
                b = j;
            }
        }
        if (best as f64) < 1e-30 {
        } else if b > 0 && b < ns - 1 {
            let sl = w[b - 1] * w[b - 1] / q[b - 1];
            let sr = w[b + 1] * w[b + 1] / q[b + 1];
            let p = (0.5 * (sl - sr) as f64 / (sl as f64 - 2. * best as f64 + sr as f64 + 1e-30))
                as f32;
            dip[n] = smin + b as f32 * step + p * step;
            sem[n] = ((best as f64 - 0.25 * (sl - sr) as f64 * p as f64) / count as f64) as f32;
            pow[n] = q[b] / count as f32;
        } else {
            dip[n] = smin + b as f32 * step;
            sem[n] = best / count as f32;
            pow[n] = q[b] / count as f32;
        }
    }
    (dip, pow, sem)
}
fn smooth(v: &[f32], w: &[f32], r: usize) -> Vec<f32> {
    let n = v.len();
    let mut o = vec![0.; n];
    if 2 * r >= n {
        return o;
    }
    for i in r..n - r {
        let mut tw = 1e-30f32;
        let mut tv = 0f32;
        for k in i - r..=i + r {
            tw += w[k];
            tv += w[k] * v[k];
        }
        o[i] = tv / tw;
    }
    for i in 0..r {
        o[i] = o[r];
    }
    for i in n - r..n {
        o[i] = o[n - r - 1];
    }
    o
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut lockey = "offset".to_string();
    if countparval("lockey") > 0 {
        getparstring("lockey", &mut lockey);
    }
    let lk = GetCase(&lockey);
    if lk < 1 {
        err!(
            "**** Error: Specified lockey name %s is not recognized.",
            lockey
        );
    }
    let mut groupkey = "cdp".to_string();
    if countparval("groupkey") > 0 {
        getparstring("groupkey", &mut groupkey);
    }
    let group = groupkey != "none";
    let gk = GetCase(&groupkey);
    if group && gk < 1 {
        err!(
            "**** Error: Specified groupkey name %s is not recognized.",
            groupkey
        );
    }
    let mut scan = 75f64;
    getpardouble("locmaxscan", &mut scan);
    if scan < 0. {
        err!("error: locmaxscan must be greater than or equal to 0.");
    }
    if scan == 0.
        && (countparval("scanstep") > 0
            || countparval("scanmax") > 0
            || countparval("scanmin") > 0
            || countparval("dipext") > 0
            || countparval("locmaxsum") < 1)
    {
        err!("error: locmaxscan=0 does not allow parameters scanstep, scanmax, scanmin, dipext, and REQUIRES locmaxsum");
    }
    let mut sum = scan;
    getpardouble("locmaxsum", &mut sum);
    if sum <= 0. {
        err!("error: locmaxsum must be greater than 0.");
    }
    let mut dipuse = 1;
    getparint("dipuse", &mut dipuse);
    if dipuse < 0 || dipuse > 1 {
        err!("error: dipuse parameter is out-of-range.");
    }
    let mut ext = 12f32;
    getparfloat("dipext", &mut ext);
    if ext < 0. {
        err!("error: dipext must be equal to or greater than 0.0");
    }
    let mut step = 0.05f32;
    getparfloat("scanstep", &mut step);
    if step < 0. {
        err!("error: scanstep must be equal to or greater than 0.0");
    }
    let mut max = 1.1f32;
    getparfloat("scanmax", &mut max);
    let mut min = -max;
    getparfloat("scanmin", &mut min);
    if max < min {
        err!("error: scanmax must be equal to or greater than scanmin");
    }
    let mut fold = 999i32;
    getparint("maxfold", &mut fold);
    if fold < 2 {
        err!("error: maxfold must be greater than 1");
    }
    let mut name = String::new();
    let od = if countparval("outdip") > 0 {
        getparstring("outdip", &mut name);
        Some(efopen(&name, "w"))
    } else {
        None
    };
    let op = if countparval("outpow") > 0 {
        getparstring("outpow", &mut name);
        Some(efopen(&name, "w"))
    } else {
        None
    };
    let os = if countparval("outsem") > 0 {
        getparstring("outsem", &mut name);
        Some(efopen(&name, "w"))
    } else {
        None
    };
    let inp = if countparval("indip") > 0 {
        if countparval("scanstep") > 0
            || countparval("scanmax") > 0
            || countparval("scanmin") > 0
            || countparval("dipext") > 0
            || countparval("locmaxscan") > 0
        {
            err!("error: indip file does not allow parameters scanstep, scanmax, scanmin, dipext, and locmaxscan");
        }
        getparstring("indip", &mut name);
        Some(efopen(&name, "r"))
    } else {
        None
    };
    checkpars();
    let mut tr = Trace::new();
    let mut td = Trace::new();
    let mut ts: Vec<Trace> = Vec::new();
    let mut loc: Vec<i32> = Vec::new();
    let mut dips: Vec<Vec<f32>> = Vec::new();
    let mut grp = -1e30f64;
    let mut ni = 0;
    let mut no = 0;
    let mut dt = 0f32;
    let mut nt = 0usize;
    let mut radius = 0usize;
    let mut ns = 0usize;
    loop {
        let has = gettr(&mut tr) != 0;
        if has {
            if ni == 0 {
                dt = (tr.dt as f64 / 1000.) as f32;
                step /= dt;
                min /= dt;
                max /= dt;
                radius = ((ext as f64 + 0.01) / dt as f64) as usize;
                ns = (1f32 + (max - min) / step) as usize;
                nt = tr.ns as usize;
            }
            ni += 1;
        }
        let changed = has && group && grp != fromhead(&tr, gk);
        if changed {
            grp = fromhead(&tr, gk);
        }
        if !has || changed {
            if inp.is_none() {
                dips.clear();
                for j in 0..ts.len() {
                    let (d, p, s) = scanit(&ts, &loc, loc[j] as f64, scan, min, step, ns, nt);
                    let d = smooth(&d, &p, radius);
                    let mut o = ts[j].clone();
                    if let Some(ref f) = od {
                        for i in 0..nt {
                            o.data[i] = d[i] * dt;
                        }
                        fputtr(f, &mut o);
                    }
                    if let Some(ref f) = op {
                        o.data[..nt].copy_from_slice(&p);
                        fputtr(f, &mut o);
                    }
                    if let Some(ref f) = os {
                        o.data[..nt].copy_from_slice(&s);
                        fputtr(f, &mut o);
                    }
                    dips.push(d);
                }
            }
            for j in 0..ts.len() {
                let mut sums: Vec<Vec<f32>> = Vec::new();
                for k in 0..ts.len() {
                    let diff = loc[k] as f64 - loc[j] as f64;
                    if diff >= -sum && diff <= sum {
                        let d = if dipuse > 0 { &dips[k] } else { &dips[j] };
                        let mut x = vec![0f32; nt];
                        for i in 0..nt {
                            x[i] = (i as f64 + d[i] as f64 * diff) as f32;
                        }
                        let mut y = vec![0f32; nt];
                        ints8r(1., 0., &ts[k].data[..nt], 0., 0., &x, &mut y);
                        sums.push(y);
                    }
                }
                let mut o = ts[j].clone();
                for i in 0..nt {
                    o.data[i] = 0.;
                    for s in &sums {
                        o.data[i] += s[i];
                    }
                    o.data[i] /= sums.len() as f32;
                }
                puttr(&mut o);
                no += 1;
            }
            ts.clear();
            loc.clear();
            dips.clear();
        }
        if !has {
            break;
        }
        if ts.len() >= fold as usize {
            err!(
                "There are too many traces in groupkey value %g    Parameter maxfold=%d ",
                grp,
                fold
            );
        }
        let x = fromhead(&tr, lk) as i32;
        loc.push(x);
        ts.push(tr.clone());
        if let Some(ref f) = inp {
            if fgettr(f, &mut td) == 0 {
                err!("next input seismic trace needs a trace from indip file: %s   but get failed (near location=%d)",name,x);
            }
            if x as f64 != fromhead(&td, lk) {
                err!("lockey value for dipfile trace %g   and input seismic trace %d are not the same ",fromhead(&td,lk),x);
            }
            dips.push((0..nt).map(|i| td.data[i] / dt).collect());
        }
    }
    if let Some(f) = inp {
        efclose(&f);
    }
    if let Some(f) = od {
        efclose(&f);
    }
    if let Some(f) = op {
        efclose(&f);
    }
    if let Some(f) = os {
        efclose(&f);
    }
    warn!("Number of traces input=%d.  Number output=%d.", ni, no);
    exit(CWP_Exit());
}
