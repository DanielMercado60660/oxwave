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

//! SURESCSV - Surface Consistent Residual Statics Determination.
//! Port of su/main/statics/surescsv.c (Revision 1.0, 2023/08/25).
//! Credits (from the C source):
//! Andre Latour. June 2023
//! 1. Looked at suresstat.c but no code from that program survives herein.
//! 2. Copied the unique-identifier-keys code from sustackup.c
//! 3. Copied (the needed part of) the kdtree code from sunearcsv.c

use oxwave_core::prelude::*;
use oxwave_header_case::{fromhead, GetCase};
use oxwave_sinc_interp::ints8r;
use std::ffi::c_void;

unsafe extern "C" {
    fn qsort(
        base: *mut c_void,
        count: usize,
        size: usize,
        compare: unsafe extern "C" fn(*const c_void, *const c_void) -> i32,
    );
}

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SURESCSV - Surface Consistent Residual Statics Determination.         ",
    "                                                                       ",
    " Notes:                                                                ",
    " Derives surface-consistent residual statics for sources and receivers.",
    " There is one static value estimated for each identified source and    ",
    " one static value estimated for each identified receiver.              ",
    " The method here is based on stack-power maximization as described by: ",
    " Ronen, J. and Claerbout, J., 1985, Geophysics, vol. 50, 2759-2767.    ",
    "                                                                       ",
    " ----------------------------------------------------------------------",
    " Conceptually, the basic methodology herein is:                        ",
    "    - store all (nmo-corrected) input traces.                          ",
    "    - stack those traces into cdps.                                    ",
    "    - cross-correlate each trace of a source with the stacked trace    ",
    "      that has the same cdp number.                                    ",
    "    - sum those cross-correlations for each source.                    ",
    "    - find time of maximum of sum of cross-correlation for each source.",
    "    - perform a 3-point interpolation around that time using the       ",
    "      neighbouring values (yielding a fractional sample shift).        ",
    "                                                                       ",
    "    - (by default) reduce long spatial wavelengths in source shifts.   ",
    "                                                                       ",
    "    - restack the cdps with the derived source time shifts applied by  ",
    "      the ints8r routine (which is able to perform fractional shifts). ",
    "                                                                       ",
    "    - repeat the procedure above, but for receivers. That is,          ",
    "      cross-correlate each trace of a receiver with the stacked trace  ",
    "      that has the same cdp number, and so on...                       ",
    "                                                                       ",
    "    - repeat everything above for the amount of iterations you choose. ",
    "                                                                       ",
    "    - output the source and receiver statics in csv files formatted    ",
    "      for use by multiple SU programs (and spreadsheet programs).      ",
    " ----------------------------------------------------------------------",
    "                                                                       ",
    " surescsv  <stdin                                                      ",
    "                                                                       ",
    " skeyloc=fldr    List of keys that identify source locations.          ",
    " \t\t\t\t\t\t\t\t\t",
    " rkeyloc=gaps    List of keys that identify receiver locations.        ",
    " \t\t\t\t\t\t\t\t\t",
    " ckeyloc=cdp     List of keys that identify cdp locations.             ",
    "           Note: For 3D surveys, you may wish to use igi,igc keys here ",
    "                 rather than cdp key. See next parameter (cdivider).   ",
    " \t\t\t\t\t\t\t\t\t",
    " cdivider=4.0    List of dividers for ckeyloc values. This parameter   ",
    "                 makes multiple input cdps into 1 larger cdp.          ",
    "                 For complicated reasons, typical survey geometry may  ",
    "                 cause this program to produce statics with repeating  ",
    "                 larger-smaller values (a sawtooth or zigzag pattern). ",
    "                 This is known as decoupling and is more likely to     ",
    "                 occur when neighbouring cdps contain traces with      ",
    "                 unique sets of offsets or receivers or sources.       ",
    "                 This parameter effectively bundles multiple input     ",
    "                 cdps together during static derivation.               ",
    " \t\t\t\t\t\t\t\t\t",
    " maxsources=10000 Maximum Amount of Sources. This is also the          ",
    "                  maximum unique combinations of skeyloc values.       ",
    "                                                                       ",
    " maxreceivers=20000 Maximum Amount of Receivers. This is also the      ",
    "                    maximum unique combinations of rkeyloc values.     ",
    "                                                                       ",
    " maxcdps=40000 Maximum Amount of CDPs (divided by cdivider). This is   ",
    "               also maximum unique combinations of ckeyloc / cdivider. ",
    "                                                                       ",
    " maxtraces=1000000 Maximum Amount of Traces. Note that most memory for ",
    "                   the traces is only allocated for the actual amount  ",
    "                   of traces read-in.                                  ",
    "                                                                       ",
    " maxshift=40.     Maximum Time Shift Per Iteration (milliseconds).     ",
    "                  (The output statics can be larger than this).        ",
    "                                                                       ",
    " numiter=5        Number of Iterations to Perform.                     ",
    "                                                                       ",
    " useoff=          Maximum offset to include traces while deriving the  ",
    "                  time shifts. Rapid changes in NMO velocity, dip,     ",
    "                  structure, or spread layout may result in some bad   ",
    "                  static values. These issues may be reduced by only   ",
    "                  using traces with smaller offsets.                   ",
    "             ***  If not specified, this parameter value defaults to   ",
    "                  the average offset of all traces multipled by 1.5    ",
    "                                                                       ",
    " iuseoff=0        Iteration number to start only using offsets less or ",
    "                  equal to parameter useoff while deriving shifts.     ",
    "                  All iterations less than this value use all offsets. ",
    "                  Example: If this value is 2 then iterations 0,1,2    ",
    "                           use all offsets, but iterations 3,4,5,..    ",
    "                           only use offsets less or equal to useoff.   ",
    "                           The default (0) therefor means that offsets ",
    "                           are limited to less than parameter useoff   ",
    "                           for all iterations.                         ",
    "                                                                       ",
    " longwave=        For complicated reasons, inaccurate NMO correction,  ",
    "                  changes in dip, or changes in geometry layouts may   ",
    "                  cause the statics to contain errors in long (spatial)",
    "                  wavelengths (roughly the spread length). For every   ",
    "                  iteration the time shifts for all sources within this",
    "                  radius of each individual source are averaged and    ",
    "                  subtracted from the shift of that individual source. ",
    "                  This reduces the longer spatial wavelengths in the   ",
    "                  output shifts. This is also done for the receivers.  ",
    "             ***  This parameter uses the sx,sy and gx,gy key values   ",
    "                  of the first trace input for each source location    ",
    "                  and each receiver location.                          ",
    "             ***  If not specified, this parameter value defaults to   ",
    "                  the average offset of all traces.                    ",
    "            =-1.0 Do not reduce the longer spatial wavelengths.        ",
    "            Note: A large value here acts similarly to option -1.      ",
    "                                                                       ",
    " srderive=0       Derive time shifts for both sources and receivers.   ",
    "              =-1 Sources only.                                        ",
    "              =1  Receivers only.                                      ",
    "            Note: This is not the same as just applying one set of     ",
    "                  statics since all time differences are resolved as   ",
    "                  far as possible into whichever points are chosen.    ",
    "                                                                       ",
    " sout=sstat.csv   Output source statics file name.                     ",
    "                  This is a comma-separated-values q-file useable by   ",
    "                  multiple SU programs, and spreadsheet programs.      ",
    "                  It contains skeyloc values, source coordinates,      ",
    "                  input source statics (named numb1_sstat),            ",
    "                  source fold after useoff (named numb2_sfold), and    ",
    "                  surface consistent residual source static sstat.     ",
    "            Note: The input source statics are copied to numb1_sstat   ",
    "                  but NOT APPLIED OR USED IN ANY WAY by this program.  ",
    "                  They are output because they may have some benefits  ",
    "                  (including computation of simple CDP floating datum).",
    "                                                                       ",
    " rout=rstat.csv   Output receiver statics file name.                   ",
    "                  This is a comma-separated-values q-file useable by   ",
    "                  multiple SU programs, and spreadsheet programs.      ",
    "                  It contains rkeyloc values, receiver coordinates,    ",
    "                  input receiver statics (named numb1_gstat),          ",
    "                  receiver fold after useoff (named numb2_gfold), and  ",
    "                  surface consistent residual receiver static gstat.   ",
    "            Note: The input receiver statics are copied to numb1_gstat ",
    "                  but NOT APPLIED OR USED IN ANY WAY by this program.  ",
    "                  They are output because they may have some benefits  ",
    "                  (including computation of simple CDP floating datum).",
    "                                                                       ",
];

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct OrderInfo {
    ndxa: i32,
    ndxb: i32,
    ndxc: i32,
    ndxd: i32,
}

unsafe extern "C" fn comparenn(q1: *const c_void, q2: *const c_void) -> i32 {
    // SAFETY: qsort passes aligned pointers to initialized OrderInfo elements
    // or temporary copies. The comparator only reads ndxa and cannot unwind.
    let p1 = unsafe { &*q1.cast::<OrderInfo>() };
    let p2 = unsafe { &*q2.cast::<OrderInfo>() };
    if p1.ndxa < p2.ndxa {
        return -1;
    }
    if p1.ndxa > p2.ndxa {
        return 1;
    }
    0
}

fn comparedd(p1: &[f64], p2: &[f64]) -> i32 {
    for n in 0..p1.len() {
        if p1[n] < p2[n] {
            return -1;
        }
        if p1[n] > p2[n] {
            return 1;
        }
    }
    0
}
fn bhighdd(all: &[Vec<f64>], guy: &[f64]) -> usize {
    let mut low = 0;
    let mut high = all.len();
    while low < high {
        let mid = low + (high - low) / 2;
        if comparedd(guy, &all[mid]) >= 0 {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    low
}

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    if !isatty(1) {
        err!("**** Error: this program does not output traces.");
    }

    let mut skeyloc = vec!["fldr".to_string()];
    if countparval("skeyloc") > 0 {
        getparstringarray("skeyloc", &mut skeyloc);
    }
    let mut skaseloc = Vec::new();
    let mut swriteloc = Vec::new();
    for key in &skeyloc {
        let kase = GetCase(key);
        if kase < 1 {
            err!(
                "**** Error: Specified skeyloc name %s is not recognized.",
                key
            );
        }
        skaseloc.push(kase);
        swriteloc.push(key != "sx" && key != "sy" && key != "sstat");
    }
    let mut maxsources = 10000i32;
    if countparval("maxsources") > 0 {
        getparint("maxsources", &mut maxsources);
    }
    if maxsources < 1 {
        err!("**** Error: maxsources must be greater than 0");
    }
    let mut sunqlocs: Vec<Vec<f64>> = Vec::new();
    let mut sunqx: Vec<f64> = Vec::new();
    let mut sunqy: Vec<f64> = Vec::new();
    let mut sunqs = vec![0.0f32; maxsources as usize];
    let mut sunqi = vec![0.0f32; maxsources as usize];
    let mut sunqq: Vec<f32> = Vec::new();
    let mut sunqf = vec![0i32; maxsources as usize];

    let mut rkeyloc = vec!["gaps".to_string()];
    if countparval("rkeyloc") > 0 {
        getparstringarray("rkeyloc", &mut rkeyloc);
    }
    let mut rkaseloc = Vec::new();
    let mut rwriteloc = Vec::new();
    for key in &rkeyloc {
        let kase = GetCase(key);
        if kase < 1 {
            err!(
                "**** Error: Specified rkeyloc name %s is not recognized.",
                key
            );
        }
        rkaseloc.push(kase);
        rwriteloc.push(key != "gx" && key != "gy" && key != "gstat");
    }
    let mut maxreceivers = 20000i32;
    if countparval("maxreceivers") > 0 {
        getparint("maxreceivers", &mut maxreceivers);
    }
    if maxreceivers < 1 {
        err!("**** Error: maxreceivers must be greater than 0");
    }
    let mut runqlocs: Vec<Vec<f64>> = Vec::new();
    let mut runqx: Vec<f64> = Vec::new();
    let mut runqy: Vec<f64> = Vec::new();
    let mut runqs = vec![0.0f32; maxreceivers as usize];
    let mut runqi = vec![0.0f32; maxreceivers as usize];
    let mut runqq: Vec<f32> = Vec::new();
    let mut runqf = vec![0i32; maxreceivers as usize];

    let mut ckeyloc = vec!["cdp".to_string()];
    if countparval("ckeyloc") > 0 {
        getparstringarray("ckeyloc", &mut ckeyloc);
    }
    let mut ckaseloc = Vec::new();
    for key in &ckeyloc {
        let kase = GetCase(key);
        if kase < 1 {
            err!(
                "**** Error: Specified ckeyloc name %s is not recognized.",
                key
            );
        }
        ckaseloc.push(kase);
    }
    let mut maxcdps = 40000i32;
    if countparval("maxcdps") > 0 {
        getparint("maxcdps", &mut maxcdps);
    }
    if maxcdps < 1 {
        err!("**** Error: maxcdps must be greater than 0");
    }
    let mut cdivider = vec![4.0f64; ckeyloc.len()];
    if countparval("cdivider") > 0 {
        if countparval("cdivider") != ckeyloc.len() as i32 {
            err!("**** Error: cdivider list must have same amount as ckeyloc ");
        }
        getpardouble("cdivider", &mut cdivider);
        for &value in &cdivider {
            if value <= 0.0 {
                err!("**** Error: cdivider values must be greater than 0.0 ");
            }
        }
    }
    let mut cunqlocs: Vec<Vec<f64>> = Vec::new();
    let mut cunqf = vec![0i32; maxcdps as usize];
    let mut maxtraces = 1000000i32;
    if countparval("maxtraces") > 0 {
        getparint("maxtraces", &mut maxtraces);
    }
    if maxtraces < 1 {
        err!("**** Error: maxtraces must be greater than 0");
    }
    let mut strclocs: Vec<Vec<f64>> = Vec::new();
    let mut rtrclocs: Vec<Vec<f64>> = Vec::new();
    let mut ctrclocs: Vec<Vec<f64>> = Vec::new();
    let mut sndx = Vec::new();
    let mut rndx = Vec::new();
    let mut cndx = Vec::new();
    let mut allsamps: Vec<Vec<f32>> = Vec::new();

    let mut maxshift = 40.0f32;
    if countparval("maxshift") > 0 {
        getparfloat("maxshift", &mut maxshift);
        if maxshift <= 0.0 {
            err!("**** Error: maxshift must be greater than 0.0");
        }
    }
    let mut numiter = 5i32;
    if countparval("numiter") > 0 {
        getparint("numiter", &mut numiter);
    }
    if numiter < 1 {
        err!("**** Error: numiter must be greater than 0");
    }
    let mut useoff = -10.0f64;
    if countparval("useoff") > 0 {
        getpardouble("useoff", &mut useoff);
        if useoff < 0.0 {
            err!("**** Error: useoff must be greater than or equal to 0.0");
        }
    }
    let mut iuseoff = 0i32;
    if countparval("iuseoff") > 0 {
        getparint("iuseoff", &mut iuseoff);
        if iuseoff < 0 {
            err!("**** Error: iuseoff must be greater than or equal to 0");
        }
    }
    let mut longwave = -10.0f64;
    if countparval("longwave") > 0 {
        getpardouble("longwave", &mut longwave);
        if longwave != -1.0 && longwave <= 0.0 {
            err!("**** Error: longwave must be -1.0 or greater than 0.0");
        }
    }
    let mut srderive = 0i32;
    if countparval("srderive") > 0 {
        getparint("srderive", &mut srderive);
    }
    if srderive < -1 || srderive > 1 {
        err!("**** Error: srderive option is out-of-range.");
    }
    let mut sout = "sstat.csv".to_string();
    if countparval("sout") > 0 {
        getparstring("sout", &mut sout);
    }
    let fp_s = fopen(&sout, "w").unwrap_or_else(|| {
        err!("error: output source static file did not open correctly.");
    });
    let mut rout = "rstat.csv".to_string();
    if countparval("rout") > 0 {
        getparstring("rout", &mut rout);
    }
    let fp_r = fopen(&rout, "w").unwrap_or_else(|| {
        err!("error: output receiver static file did not open correctly.");
    });
    checkpars();

    let mut tr = Trace::new();
    let mut lensam = 0usize;
    let mut sampi = 0.0f32;
    let mut lenxcor = 0usize;
    let mut avoff = 0.0f64;
    let mut nproct = 0usize;
    while gettr(&mut tr) != 0 {
        if nproct < 1 {
            lensam = tr.ns as usize;
            sampi = (tr.dt as f64 / 1000.0) as f32;
            lenxcor = ((maxshift as f64 + 0.00001) / sampi as f64) as i32 as usize;
        }
        if nproct >= maxtraces as usize {
            err!(
                "error: Number of input traces is greater than maxtraces=%d",
                maxtraces
            );
        }
        let sq = tr.sstat as f32;
        let rq = tr.gstat as f32;
        let mut sx = tr.sx as f64;
        let mut sy = tr.sy as f64;
        let mut rx = tr.gx as f64;
        let mut ry = tr.gy as f64;
        if tr.scalco > 1 {
            sx *= tr.scalco as f64;
            sy *= tr.scalco as f64;
            rx *= tr.scalco as f64;
            ry *= tr.scalco as f64;
        } else if tr.scalco < 0 {
            sx /= -(tr.scalco as f64);
            sy /= -(tr.scalco as f64);
            rx /= -(tr.scalco as f64);
            ry /= -(tr.scalco as f64);
        }
        avoff += ((sx - rx) * (sx - rx) + (sy - ry) * (sy - ry)).sqrt();
        let mut sguy = vec![0.0; skaseloc.len()];
        for i in 0..skaseloc.len() {
            sguy[i] = fromhead(&tr, skaseloc[i]);
        }
        let ihere = bhighdd(&sunqlocs, &sguy);
        if ihere > 0 && comparedd(&sguy, &sunqlocs[ihere - 1]) == 0 {
            if (sx - sunqx[ihere - 1]).abs() > 0.01 || (sy - sunqy[ihere - 1]).abs() != 0.0 {
                warn!(
                    "trace %d has different sx,sy than other traces with same skeyloc values ",
                    nproct + 1
                );
            }
        } else {
            if sunqlocs.len() >= maxsources as usize {
                err!("error: At input trace=%d  number of unique skeyloc combinations is greater than maxsources=%d", nproct, maxsources);
            }
            sunqlocs.insert(ihere, sguy.clone());
            sunqx.insert(ihere, sx);
            sunqy.insert(ihere, sy);
            sunqq.insert(ihere, sq);
        }
        let mut rguy = vec![0.0; rkaseloc.len()];
        for i in 0..rkaseloc.len() {
            rguy[i] = fromhead(&tr, rkaseloc[i]);
        }
        let jhere = bhighdd(&runqlocs, &rguy);
        if jhere > 0 && comparedd(&rguy, &runqlocs[jhere - 1]) == 0 {
            if (rx - runqx[jhere - 1]).abs() > 0.01 || (ry - runqy[jhere - 1]).abs() != 0.0 {
                warn!(
                    "trace %d has different gx,gy than other traces with same rkeyloc values ",
                    nproct + 1
                );
            }
        } else {
            if runqlocs.len() >= maxreceivers as usize {
                err!("error: At input trace=%d  number of unique rkeyloc combinations is greater than maxreceivers=%d", nproct, maxreceivers);
            }
            runqlocs.insert(jhere, rguy.clone());
            runqx.insert(jhere, rx);
            runqy.insert(jhere, ry);
            runqq.insert(jhere, rq);
        }
        let mut cguy = vec![0.0; ckaseloc.len()];
        for i in 0..ckaseloc.len() {
            cguy[i] = (fromhead(&tr, ckaseloc[i]) / cdivider[i]).floor();
        }
        let khere = bhighdd(&cunqlocs, &cguy);
        if khere > 0 && comparedd(&cguy, &cunqlocs[khere - 1]) == 0 {
        } else {
            if cunqlocs.len() >= maxcdps as usize {
                err!("error: At input trace=%d  number of unique ckeyloc combinations is greater than maxcdps=%d", nproct, maxcdps);
            }
            cunqlocs.insert(khere, cguy.clone());
        }
        let mut samples = vec![0.0f32; lensam];
        for i in 0..lensam {
            samples[i] = tr.data[i];
        }
        allsamps.push(samples);
        strclocs.push(sguy);
        rtrclocs.push(rguy);
        ctrclocs.push(cguy);
        nproct += 1;
    }
    let lastsource = sunqlocs.len();
    let lastreceiver = runqlocs.len();
    let lastcdp = cunqlocs.len();
    let mut cdpamp = vec![vec![0.0f32; lensam]; lastcdp];
    let mut samplive = vec![0i32; lensam];
    let mut onetim = vec![0.0f32; lensam];
    let mut trshifted = vec![0.0f32; lensam];
    let mut tpad = vec![0.0f32; lensam + 2 * lenxcor + 1];
    let mut txcor = vec![0.0f64; 2 * lenxcor + 1];
    for j in 0..nproct {
        let ihere = bhighdd(&sunqlocs, &strclocs[j]) - 1;
        let jhere = bhighdd(&runqlocs, &rtrclocs[j]) - 1;
        let khere = bhighdd(&cunqlocs, &ctrclocs[j]) - 1;
        sndx.push(OrderInfo {
            ndxa: ihere as i32,
            ndxb: j as i32,
            ndxc: khere as i32,
            ndxd: jhere as i32,
        });
        rndx.push(OrderInfo {
            ndxa: jhere as i32,
            ndxb: j as i32,
            ndxc: khere as i32,
            ndxd: ihere as i32,
        });
        cndx.push(OrderInfo {
            ndxa: khere as i32,
            ndxb: j as i32,
            ndxc: ihere as i32,
            ndxd: jhere as i32,
        });
    }
    // libc qsort's tie order determines the floating-point accumulation order
    // and therefore the decimal statics written to the CSV files.
    // SAFETY: each vector contains nproct initialized, Copy OrderInfo elements
    // with the C layout. qsort exclusively rearranges them with the exact
    // element size; comparenn only reads their integer keys and cannot unwind.
    unsafe {
        qsort(
            sndx.as_mut_ptr().cast(),
            nproct,
            std::mem::size_of::<OrderInfo>(),
            comparenn,
        );
        qsort(
            rndx.as_mut_ptr().cast(),
            nproct,
            std::mem::size_of::<OrderInfo>(),
            comparenn,
        );
        qsort(
            cndx.as_mut_ptr().cast(),
            nproct,
            std::mem::size_of::<OrderInfo>(),
            comparenn,
        );
    }
    avoff /= nproct as f64;
    if useoff < -9.0 {
        useoff = 1.5 * avoff;
    }
    useoff *= useoff;
    if longwave < -9.0 {
        longwave = avoff;
    }
    let tree_nodes = if longwave > 0.0 {
        connect_nodes(&[&sunqx, &sunqy])
    } else {
        Vec::new()
    };
    let tree_noder = if longwave > 0.0 {
        connect_nodes(&[&runqx, &runqy])
    } else {
        Vec::new()
    };
    cndx.push(OrderInfo {
        ndxa: -1,
        ..OrderInfo::default()
    });
    sndx.push(OrderInfo {
        ndxa: -1,
        ..OrderInfo::default()
    });
    rndx.push(OrderInfo {
        ndxa: -1,
        ..OrderInfo::default()
    });
    for iiter in 0..numiter {
        if srderive <= 0 {
            // Stack the traces into their cdps (restack after source statics).
            for i in 0..lensam {
                samplive[i] = 0;
                cdpamp[0][i] = 0.0;
            }
            let mut nowcdp = 0;
            cunqf[0] = 0;
            for j in 0..nproct + 1 {
                let mut inear = 1;
                if iiter >= iuseoff {
                    let toff = (sunqx[cndx[j].ndxc as usize] - runqx[cndx[j].ndxd as usize])
                        * (sunqx[cndx[j].ndxc as usize] - runqx[cndx[j].ndxd as usize])
                        + (sunqy[cndx[j].ndxc as usize] - runqy[cndx[j].ndxd as usize])
                            * (sunqy[cndx[j].ndxc as usize] - runqy[cndx[j].ndxd as usize]);
                    if toff > useoff {
                        inear = 0;
                    }
                }
                if inear > 0 {
                    for k in 0..lensam {
                        onetim[k] = k as f32
                            - (sunqs[cndx[j].ndxc as usize] + runqs[cndx[j].ndxd as usize]);
                    }
                    ints8r(
                        1.0,
                        0.0,
                        &allsamps[cndx[j].ndxb as usize],
                        0.0,
                        0.0,
                        &onetim,
                        &mut trshifted,
                    );
                } else {
                    for k in 0..lensam {
                        trshifted[k] = 0.0;
                    }
                }
                if cndx[j].ndxa == nowcdp as i32 {
                    for i in 0..lensam {
                        if trshifted[i] != 0.0 {
                            cdpamp[nowcdp][i] += trshifted[i];
                            samplive[i] += 1;
                        }
                    }
                    cunqf[nowcdp] += inear;
                } else {
                    for i in 0..lensam {
                        if samplive[i] > 1 {
                            cdpamp[nowcdp][i] /= samplive[i] as f32;
                        }
                    }
                    nowcdp += 1;
                    if nowcdp == lastcdp {
                        break;
                    }
                    cunqf[nowcdp] = 1;
                    for i in 0..lensam {
                        cdpamp[nowcdp][i] = trshifted[i];
                        samplive[i] = if trshifted[i] != 0.0 { 1 } else { 0 };
                    }
                }
            }
            // Derive the new source time differences from summed cross-correlations.
            for i in 0..2 * lenxcor + 1 {
                txcor[i] = 0.0;
            }
            let mut nowsource = 0usize;
            sunqf[0] = 0;
            for j in 0..nproct + 1 {
                // C accesses ndxa=-1 at the sentinel before breaking. Its shifted
                // samples are discarded; use valid coordinates for that unused work.
                let a = sndx[j].ndxa.max(0) as usize;
                let mut inear = 1;
                if iiter >= iuseoff {
                    let toff = (sunqx[a] - runqx[sndx[j].ndxd as usize])
                        * (sunqx[a] - runqx[sndx[j].ndxd as usize])
                        + (sunqy[a] - runqy[sndx[j].ndxd as usize])
                            * (sunqy[a] - runqy[sndx[j].ndxd as usize]);
                    if toff > useoff {
                        inear = 0;
                    }
                }
                if inear > 0 {
                    for k in 0..lensam {
                        onetim[k] = k as f32 - (sunqs[a] + runqs[sndx[j].ndxd as usize]);
                    }
                    ints8r(
                        1.0,
                        0.0,
                        &allsamps[sndx[j].ndxb as usize],
                        0.0,
                        0.0,
                        &onetim,
                        &mut trshifted,
                    );
                } else {
                    for k in 0..lensam {
                        trshifted[k] = 0.0;
                    }
                }
                if sndx[j].ndxa == nowsource as i32 {
                    for i in 0..lensam {
                        tpad[i + lenxcor] = cdpamp[sndx[j].ndxc as usize][i];
                    }
                    for i in 0..2 * lenxcor + 1 {
                        for k in 0..lensam {
                            txcor[i] += (trshifted[k] * tpad[k + i]) as f64;
                        }
                    }
                    sunqf[nowsource] += inear;
                } else {
                    let mut bigx = txcor[lenxcor];
                    let mut ibig = lenxcor;
                    for i in 0..2 * lenxcor + 1 {
                        if txcor[i] > bigx {
                            bigx = txcor[i];
                            ibig = i;
                        }
                    }
                    let mut tbig = ibig as f64;
                    if ibig > 0
                        && ibig < 2 * lenxcor
                        && txcor[ibig - 1] - 2.0 * bigx + txcor[ibig + 1] != 0.0
                    {
                        tbig += 0.5 * (txcor[ibig - 1] - txcor[ibig + 1])
                            / (txcor[ibig - 1] - 2.0 * bigx + txcor[ibig + 1]);
                    }
                    sunqs[nowsource] = (sunqs[nowsource] as f64 + (tbig - lenxcor as f64)) as f32;
                    nowsource += 1;
                    if nowsource == lastsource {
                        break;
                    }
                    for i in 0..2 * lenxcor + 1 {
                        txcor[i] = 0.0;
                    }
                    if inear > 0 {
                        for i in 0..lensam {
                            tpad[i + lenxcor] = cdpamp[sndx[j].ndxc as usize][i];
                        }
                        for i in 0..2 * lenxcor + 1 {
                            for k in 0..lensam {
                                txcor[i] += (trshifted[k] * tpad[k + i]) as f64;
                            }
                        }
                    }
                    sunqf[nowsource] = inear;
                }
            }
            // Remove spatial long waves in source statics.
            if longwave > 0.0 {
                for i in 0..lastsource {
                    let extent_min = [sunqx[i] - longwave, sunqy[i] - longwave];
                    let extent_max = [sunqx[i] + longwave, sunqy[i] + longwave];
                    let mut out_elem = Vec::new();
                    find_in(
                        &tree_nodes,
                        &[&sunqx, &sunqy],
                        &extent_min,
                        &extent_max,
                        &mut out_elem,
                    );
                    sunqi[i] = 0.0;
                    let mut k = 0i32;
                    for j in 0..out_elem.len() {
                        if (sunqx[i] - sunqx[out_elem[j]]) * (sunqx[i] - sunqx[out_elem[j]])
                            + (sunqy[i] - sunqy[out_elem[j]]) * (sunqy[i] - sunqy[out_elem[j]])
                            <= longwave * longwave
                        {
                            sunqi[i] += sunqs[out_elem[j]];
                            k += 1;
                        }
                    }
                    sunqi[i] /= k as f32;
                }
                for i in 0..lastsource {
                    sunqs[i] -= sunqi[i];
                }
            }
        }
        if srderive >= 0 {
            // Stack the traces into their cdps (restack after source statics).
            for i in 0..lensam {
                samplive[i] = 0;
                cdpamp[0][i] = 0.0;
            }
            let mut nowcdp = 0;
            cunqf[0] = 0;
            for j in 0..nproct + 1 {
                let mut inear = 1;
                if iiter >= iuseoff {
                    let toff = (sunqx[cndx[j].ndxc as usize] - runqx[cndx[j].ndxd as usize])
                        * (sunqx[cndx[j].ndxc as usize] - runqx[cndx[j].ndxd as usize])
                        + (sunqy[cndx[j].ndxc as usize] - runqy[cndx[j].ndxd as usize])
                            * (sunqy[cndx[j].ndxc as usize] - runqy[cndx[j].ndxd as usize]);
                    if toff > useoff {
                        inear = 0;
                    }
                }
                if inear > 0 {
                    for k in 0..lensam {
                        onetim[k] = k as f32
                            - (sunqs[cndx[j].ndxc as usize] + runqs[cndx[j].ndxd as usize]);
                    }
                    ints8r(
                        1.0,
                        0.0,
                        &allsamps[cndx[j].ndxb as usize],
                        0.0,
                        0.0,
                        &onetim,
                        &mut trshifted,
                    );
                } else {
                    for k in 0..lensam {
                        trshifted[k] = 0.0;
                    }
                }
                if cndx[j].ndxa == nowcdp as i32 {
                    for i in 0..lensam {
                        if trshifted[i] != 0.0 {
                            cdpamp[nowcdp][i] += trshifted[i];
                            samplive[i] += 1;
                        }
                    }
                    cunqf[nowcdp] += inear;
                } else {
                    for i in 0..lensam {
                        if samplive[i] > 1 {
                            cdpamp[nowcdp][i] /= samplive[i] as f32;
                        }
                    }
                    nowcdp += 1;
                    if nowcdp == lastcdp {
                        break;
                    }
                    cunqf[nowcdp] = 1;
                    for i in 0..lensam {
                        cdpamp[nowcdp][i] = trshifted[i];
                        samplive[i] = if trshifted[i] != 0.0 { 1 } else { 0 };
                    }
                }
            }
            // Derive the new receiver time differences from summed cross-correlations.
            for i in 0..2 * lenxcor + 1 {
                txcor[i] = 0.0;
            }
            let mut nowreceiver = 0usize;
            runqf[0] = 0;
            for j in 0..nproct + 1 {
                // C accesses ndxa=-1 at the sentinel before breaking. Its shifted
                // samples are discarded; use valid coordinates for that unused work.
                let a = rndx[j].ndxa.max(0) as usize;
                let mut inear = 1;
                if iiter >= iuseoff {
                    let toff = (sunqx[rndx[j].ndxd as usize] - runqx[a])
                        * (sunqx[rndx[j].ndxd as usize] - runqx[a])
                        + (sunqy[rndx[j].ndxd as usize] - runqy[a])
                            * (sunqy[rndx[j].ndxd as usize] - runqy[a]);
                    if toff > useoff {
                        inear = 0;
                    }
                }
                if inear > 0 {
                    for k in 0..lensam {
                        onetim[k] = k as f32 - (sunqs[rndx[j].ndxd as usize] + runqs[a]);
                    }
                    ints8r(
                        1.0,
                        0.0,
                        &allsamps[rndx[j].ndxb as usize],
                        0.0,
                        0.0,
                        &onetim,
                        &mut trshifted,
                    );
                } else {
                    for k in 0..lensam {
                        trshifted[k] = 0.0;
                    }
                }
                if rndx[j].ndxa == nowreceiver as i32 {
                    for i in 0..lensam {
                        tpad[i + lenxcor] = cdpamp[rndx[j].ndxc as usize][i];
                    }
                    for i in 0..2 * lenxcor + 1 {
                        for k in 0..lensam {
                            txcor[i] += (trshifted[k] * tpad[k + i]) as f64;
                        }
                    }
                    runqf[nowreceiver] += inear;
                } else {
                    let mut bigx = txcor[lenxcor];
                    let mut ibig = lenxcor;
                    for i in 0..2 * lenxcor + 1 {
                        if txcor[i] > bigx {
                            bigx = txcor[i];
                            ibig = i;
                        }
                    }
                    let mut tbig = ibig as f64;
                    if ibig > 0
                        && ibig < 2 * lenxcor
                        && txcor[ibig - 1] - 2.0 * bigx + txcor[ibig + 1] != 0.0
                    {
                        tbig += 0.5 * (txcor[ibig - 1] - txcor[ibig + 1])
                            / (txcor[ibig - 1] - 2.0 * bigx + txcor[ibig + 1]);
                    }
                    runqs[nowreceiver] =
                        (runqs[nowreceiver] as f64 + (tbig - lenxcor as f64)) as f32;
                    nowreceiver += 1;
                    if nowreceiver == lastreceiver {
                        break;
                    }
                    for i in 0..2 * lenxcor + 1 {
                        txcor[i] = 0.0;
                    }
                    if inear > 0 {
                        for i in 0..lensam {
                            tpad[i + lenxcor] = cdpamp[rndx[j].ndxc as usize][i];
                        }
                        for i in 0..2 * lenxcor + 1 {
                            for k in 0..lensam {
                                txcor[i] += (trshifted[k] * tpad[k + i]) as f64;
                            }
                        }
                    }
                    runqf[nowreceiver] = inear;
                }
            }
            // Remove spatial long waves in receiver statics.
            if longwave > 0.0 {
                for i in 0..lastreceiver {
                    let extent_min = [runqx[i] - longwave, runqy[i] - longwave];
                    let extent_max = [runqx[i] + longwave, runqy[i] + longwave];
                    let mut out_elem = Vec::new();
                    find_in(
                        &tree_noder,
                        &[&runqx, &runqy],
                        &extent_min,
                        &extent_max,
                        &mut out_elem,
                    );
                    runqi[i] = 0.0;
                    let mut k = 0i32;
                    for j in 0..out_elem.len() {
                        if (runqx[i] - runqx[out_elem[j]]) * (runqx[i] - runqx[out_elem[j]])
                            + (runqy[i] - runqy[out_elem[j]]) * (runqy[i] - runqy[out_elem[j]])
                            <= longwave * longwave
                        {
                            runqi[i] += runqs[out_elem[j]];
                            k += 1;
                        }
                    }
                    runqi[i] /= k as f32;
                }
                for i in 0..lastreceiver {
                    runqs[i] -= runqi[i];
                }
            }
        }
    }
    fprintf!(
        fp_s,
        "C_SURESCSV numiter avrgoff longwave useoff iuseoff =,%d,%.12g,%.12g,%.12g,%d\n",
        numiter,
        avoff,
        longwave,
        useoff.sqrt(),
        iuseoff
    );
    fprintf!(fp_s, "C_SU_MATCH");
    for i in 0..skeyloc.len() {
        fprintf!(fp_s, ",%s", &skeyloc[i]);
    }
    fprintf!(fp_s, "\nC_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
    for i in 0..skeyloc.len() {
        if swriteloc[i] {
            fprintf!(fp_s, "%s", ",%.12g");
        }
    }
    fprintf!(
        fp_s,
        "%s",
        ",%.12g,%.12g,%.12g,%.12g,%.12g\nC_SU_NAMES\nC_SU_ID"
    );
    for i in 0..skeyloc.len() {
        if swriteloc[i] {
            fprintf!(fp_s, ",%s", &skeyloc[i]);
        }
    }
    fprintf!(fp_s, ",sx,sy,numb1_sstat,numb2_sfold,sstat\n");
    for i in 0..lastsource {
        fprintf!(fp_s, "Q");
        for j in 0..skeyloc.len() {
            if swriteloc[j] {
                fprintf!(fp_s, ",%.12g", sunqlocs[i][j]);
            }
        }
        fprintf!(
            fp_s,
            ",%.12g,%.12g,%.12g,%d,%.12g\n",
            sunqx[i],
            sunqy[i],
            sunqq[i],
            sunqf[i],
            sunqs[i] * sampi
        );
    }
    fprintf!(
        fp_r,
        "C_SURESCSV numiter avrgoff longwave useoff iuseoff =,%d,%.12g,%.12g,%.12g,%d\n",
        numiter,
        avoff,
        longwave,
        useoff.sqrt(),
        iuseoff
    );
    fprintf!(fp_r, "C_SU_MATCH");
    for i in 0..rkeyloc.len() {
        fprintf!(fp_r, ",%s", &rkeyloc[i]);
    }
    fprintf!(fp_r, "\nC_SU_SETID,Q\nC_SU_FORMS\nC_SU_ID");
    for i in 0..rkeyloc.len() {
        if rwriteloc[i] {
            fprintf!(fp_r, "%s", ",%.12g");
        }
    }
    fprintf!(
        fp_r,
        "%s",
        ",%.12g,%.12g,%.12g,%.12g,%.12g\nC_SU_NAMES\nC_SU_ID"
    );
    for i in 0..rkeyloc.len() {
        if rwriteloc[i] {
            fprintf!(fp_r, ",%s", &rkeyloc[i]);
        }
    }
    fprintf!(fp_r, ",gx,gy,numb1_gstat,numb2_gfold,gstat\n");
    for i in 0..lastreceiver {
        fprintf!(fp_r, "Q");
        for j in 0..rkeyloc.len() {
            if rwriteloc[j] {
                fprintf!(fp_r, ",%.12g", runqlocs[i][j]);
            }
        }
        fprintf!(
            fp_r,
            ",%.12g,%.12g,%.12g,%d,%.12g\n",
            runqx[i],
            runqy[i],
            runqq[i],
            runqf[i],
            runqs[i] * sampi
        );
    }
    exit(CWP_Exit());
}

#[derive(Clone, Default)]
struct Node {
    elem: usize,
    l: Option<usize>,
    r: Option<usize>,
}

// connect_nodes, connect_all, find_in and find_in_rcur are private helpers
// from surescsv.c, rather than library kd-tree routines.
fn connect_nodes(tree_dl: &[&Vec<f64>; 2]) -> Vec<Node> {
    let tree_numc = tree_dl[0].len();
    let mut tree_nodes = vec![Node::default(); tree_numc];
    for m in 0..tree_numc {
        tree_nodes[m].elem = m;
    }
    let mut m = 0;
    let mut nrat = tree_numc;
    while nrat > 0 {
        let nd = if nrat > 6 {
            (nrat as f64).sqrt() as usize
        } else {
            0
        };
        let mut n = nd + nrat / 2;
        while n < tree_numc {
            if tree_nodes[n].l.is_none() {
                tree_nodes[n].l = Some(0); // pointer used as a flag
                tree_nodes[m].elem = n;
                m += 1;
            }
            n = n + nd + nrat;
        }
        nrat = (nrat as f64 * 0.6) as usize;
    }
    for m in 0..tree_numc {
        tree_nodes[m].l = None;
    }
    connect_all(&mut tree_nodes, tree_dl);
    tree_nodes
}
fn connect_all(tree_nodes: &mut [Node], tree_dl: &[&Vec<f64>; 2]) {
    for m in 1..tree_nodes.len() {
        let mut now_node = 0;
        let mut next_node = Some(0usize);
        let mut naxe = 0;
        while let Some(next) = next_node {
            now_node = next;
            naxe += 1;
            if naxe == 2 {
                naxe = 0;
            }
            next_node =
                if tree_dl[naxe][tree_nodes[m].elem] < tree_dl[naxe][tree_nodes[now_node].elem] {
                    tree_nodes[now_node].l
                } else {
                    tree_nodes[now_node].r
                };
        }
        if tree_dl[naxe][tree_nodes[m].elem] < tree_dl[naxe][tree_nodes[now_node].elem] {
            tree_nodes[now_node].l = Some(m);
        } else {
            tree_nodes[now_node].r = Some(m);
        }
    }
}
fn find_in(
    tree_nodes: &[Node],
    tree_dl: &[&Vec<f64>; 2],
    extent_min: &[f64; 2],
    extent_max: &[f64; 2],
    out_elem: &mut Vec<usize>,
) {
    out_elem.clear();
    find_in_rcur(
        tree_nodes,
        tree_dl,
        Some(0),
        1,
        extent_min,
        extent_max,
        out_elem,
    );
}
fn find_in_rcur(
    tree_nodes: &[Node],
    tree_dl: &[&Vec<f64>; 2],
    now_node: Option<usize>,
    mut naxe: usize,
    extent_min: &[f64; 2],
    extent_max: &[f64; 2],
    out_elem: &mut Vec<usize>,
) {
    let Some(now) = now_node else {
        return;
    };
    let node = &tree_nodes[now];
    let mut inside = true;
    for i in 0..2 {
        if tree_dl[i][node.elem] < extent_min[i] || tree_dl[i][node.elem] >= extent_max[i] {
            inside = false;
            break;
        }
    }
    if inside {
        out_elem.push(node.elem);
    }
    if naxe == 2 {
        naxe = 0;
    }
    if tree_dl[naxe][node.elem] >= extent_min[naxe] && node.l.is_some() {
        find_in_rcur(
            tree_nodes,
            tree_dl,
            node.l,
            naxe + 1,
            extent_min,
            extent_max,
            out_elem,
        );
    }
    if tree_dl[naxe][node.elem] < extent_max[naxe] && node.r.is_some() {
        find_in_rcur(
            tree_nodes,
            tree_dl,
            node.r,
            naxe + 1,
            extent_min,
            extent_max,
            out_elem,
        );
    }
}
