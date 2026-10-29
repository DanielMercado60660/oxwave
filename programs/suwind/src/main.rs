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

//! SUWIND - window traces by key word, port of su/main/windowing_sorting_muting/suwind.c (Revision 1.55).
//! Credits: SEP: Einar Kjartansson; CWP: Shuki Ronen, Jack Cohen, Chris Liner; Warnemuende: Toralf Foerster; CENPET: Werner M. Heigl; Andre Latour.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUWIND - window traces by key word\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "  suwind <stdin >stdout [options]\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    "  none \t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " verbose=0\t\t=1 for verbose\t\t\t\t\t",
    " key=tracl\t\tKey header word to window on (see segy.h)\t",
    " min=LONG_MIN\t\tmin value of key header word to pass\t\t",
    " max=LONG_MAX\t\tmax value of key header word to pass\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " abs=0\t\t\t=1 to take absolute value of key header word\t",
    " j=1\t\t\tPass every j-th trace ...\t\t\t",
    " s=0\t\t\t... based at s  (if ((key - s)%j) == 0)\t\t",
    " skip=0\t\tskip the initial N traces                       ",
    " count=ULONG_MAX\t... up to count traces\t\t\t\t",
    " reject=none\t\tSkip traces with specified key values\t\t",
    " accept=none\t\tPass traces with specified key values(see notes)",
    "\t\t\tprocessing, but do no window the data\t\t",
    " ordered=0\t\t=1 if traces sorted in increasing keyword value ",
    "\t\t\t=-1  if traces are sorted in a decreasing order ",
    "\t\t\t\t\t\t\t\t\t",
    " Options for vertical windowing (time gating):\t\t\t\t",
    " dt=tr.dt (from header) time sampling interval (sec)\t(seismic data)\t",
    " \t\t\t =tr.d1  (nonseismic)\t\t\t\t",
    " f1=tr.delrt (from header) first sample\t\t(seismic data)\t",
    " \t\t\t =tr.f1  (nonseismic)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmin=0.0\t\tmin time to pass\t\t\t\t",
    " tmax=(from header)\tmax time to pass\t\t\t\t",
    " itmin=0\t\tmin time sample to pass\t\t\t\t",
    " itmax=(from header)   max time sample to pass\t\t\t\t",
    " nt=itmax-itmin+1\tnumber of time samples to pass\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " On large data sets, the count parameter should be set if\t\t",
    " possible.  Otherwise, every trace in the data set will be\t\t",
    " examined.  However, the count parameter overrides the accept\t\t",
    " parameter, so you can't specify count if you want true\t\t",
    " unconditional acceptance.\t\t\t\t\t\t",
    "                                                                       ",
    " The skip= option allows the user to skip over traces, which helps\t",
    " for selecting traces far from the beginning of the dataset.\t\t",
    " Caveat: skip only works with disk input.                        \t",
    "                                                                       ",
    " The ordered= option will speed up the process if the data are   \t",
    " sorted in according to the key.                                 \t",
    "\t\t\t\t\t\t\t\t\t",
    " The accept option is a bit strange--it does NOT mean accept ONLY\t",
    " the traces on the accept list!  It means accept these traces,   \t",
    " even if they would otherwise be rejected (except as noted in the\t",
    " previous paragraph).  To implement accept-only, you can use the \t",
    " max=0 option (rejecting everything).  For example, to accept    \t",
    " only the tracl values 4, 5 and 6:\t\t\t\t\t",
    "\t... | suwind max=0 accept=4,5,6 | ...\t\t   \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Another example is the case of suppressing nonseismic traces in \t",
    " a seismic data set. By the SEGY standard header field trace id, \t",
    " trid=1 designates traces as being seismic traces. Other traces, \t",
    " such as calibration traces may be designated by another value.  \t",
    "\t\t\t\t\t\t\t\t\t",
    " Example:  trid=1 seismic and trid=0 is nonseismic. To reject    \t",
    "       nonseismic traces\t\t\t\t\t\t",
    "       ... | suwind key=trid reject=0 | ...\t\t\t\t",
    "      \t\t\t\t\t\t\t\t\t",
    " On most 32 bit machines, LONG_MIN, LONG_MAX and ULONG_MAX are   \t",
    " about -2E9,+2E9 and 4E9, they are defined in limits.h.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Selecting times beyond the maximum in the data induces\t\t",
    " zero padding (up to SU_NFLTS).\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The time gating here is to the nearest neighboring sample or    \t",
    " time value. Gating to the exact temporal value requires\t \t",
    " resampling if the selected times fall between samples on the    \t",
    " trace. Use suresamp to perform the time gating in this case.    \t",
    "\t\t\t\t\t\t\t\t\t",
    " It doesn't really make sense to specify both itmin and tmin,\t\t",
    " but specifying itmin takes precedence over specifying tmin.\t\t",
    " Similarly, itmax takes precedence over tmax and tmax over nt.\t\t",
    " If dt in header is not set, then dt is mandatory\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
];
fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut key = String::new();
    if getparstring("key", &mut key) == 0 {
        key = "tracl".into();
    }
    let mut min = i32::MIN as i64;
    let mut max = i32::MAX as i64;
    let mut j = 1;
    let mut s = 0;
    let mut count = u64::MAX;
    let mut ab: i16 = 0;
    let mut verbose: i16 = 0;
    let mut skip = 0;
    let mut ordered = 0;
    getparlong("min", &mut min);
    getparlong("max", &mut max);
    getparint("j", &mut j);
    getparint("s", &mut s);
    getparulong("count", &mut count);
    getparshort("abs", &mut ab);
    getparshort("verbose", &mut verbose);
    getparint("skip", &mut skip);
    getparint("ordered", &mut ordered);
    let mut itr = 0;
    let mut ntr = 0;
    if skip != 0 {
        ntr = fgettra(&stdin(), &mut tr, itr);
        if ntr == 0 {
            err!("can't get first trace");
        }
        if skip > ntr - 1 {
            err!("can't skip past last trace");
        }
        itr = skip;
        fgettra(&stdin(), &mut tr, itr);
    } else if gettr(&mut tr) == 0 {
        err!("can't get first trace");
    }
    let ntfirst = tr.ns as i32;
    let seismic = (0..=17).contains(&tr.trid) || tr.trid == 109;
    let mut dt = 0f32;
    let mut f1 = 0f32;
    if seismic {
        if verbose != 0 {
            warn!("input is seismic data, trid=%d", tr.trid);
        }
        if getparfloat("dt", &mut dt) == 0 {
            dt = (tr.dt as f64 / 1000000.0) as f32;
        }
        if dt == 0.0 {
            err!("dt not set in header or not getparred");
        }
        if getparfloat("f1", &mut f1) == 0 {
            f1 = (tr.delrt as f64 / 1000.0) as f32;
        }
        if f1 == 0.0 && verbose != 0 {
            warn!("delrt not set in header or not getparred");
        }
    } else {
        if verbose != 0 {
            warn!("input is not seismic data, trid=%d", tr.trid);
        }
        if getparfloat("dt", &mut dt) == 0 {
            if tr.d1 != 0.0 {
                dt = tr.d1;
            } else {
                if verbose != 0 {
                    warn!("d1 field not set");
                    warn!("using dt field");
                }
                dt = tr.dt as f32;
                if dt == 0.0 {
                    err!("d1 and dt not set in header or dt not getparred");
                }
            }
        }
        if getparfloat("f1", &mut f1) == 0 {
            f1 = tr.f1;
        }
        if f1 == 0.0 && verbose != 0 {
            warn!("f1 not set in header or not getparred");
        }
    }
    let mut tmin = 0f32;
    let mut tmax = 0f32;
    let mut itmin = 0i32;
    let mut itmax = 0i32;
    let mut nt = 0i32;
    if getparint("itmin", &mut itmin) != 0 {
        tmin = itmin as f32 * dt + f1;
    } else if getparfloat("tmin", &mut tmin) != 0 {
        itmin = nint(((tmin - f1) / dt) as f64);
    } else {
        itmin = 0;
        tmin = f1;
    }
    if getparint("itmax", &mut itmax) != 0 {
        tmax = itmax as f32 * dt + f1;
        nt = itmax.wrapping_sub(itmin).wrapping_add(1);
    } else if getparfloat("tmax", &mut tmax) != 0 {
        itmax = nint(((tmax - f1) / dt) as f64);
        nt = itmax.wrapping_sub(itmin).wrapping_add(1);
    } else if getparint("nt", &mut nt) != 0 {
        itmax = itmin.wrapping_add(nt).wrapping_sub(1);
        tmax = itmax as f32 * dt + f1;
    } else {
        itmax = ntfirst - 1;
        tmax = itmax as f32 * dt + f1;
        nt = itmax.wrapping_sub(itmin).wrapping_add(1);
    }
    if itmin < 0 {
        err!("itmin=%d should be positive", itmin);
    }
    if nt > SU_NFLTS as i32 {
        err!("nt=%d exceeds SU_NFLTS=%d", nt, SU_NFLTS);
    }
    if itmin > itmax {
        err!("itmin=%d, itmax=%d conflict", itmin, itmax);
    }
    if nt <= 0 {
        err!("invalid window length nt=%d", nt);
    }
    let nzeros = (nt - ntfirst).max(0) as usize;
    let type_ = hdtype(&key);
    let indx = getindex(&key);
    let nbad = countparval("reject") as usize;
    let mut bad = vec![Value::default(); nbad];
    if nbad != 0 {
        getparval("reject", type_, nbad, &mut bad);
    }
    let ngood = countparval("accept") as usize;
    let mut good = vec![Value::default(); ngood];
    if ngood != 0 {
        getparval("accept", type_, ngood, &mut good);
    }
    checkpars();
    if verbose != 0 {
        warn!("key = %s", key);
        warn!("type = %s", type_);
        warn!("min = %ld", min);
        warn!("max = %ld", max);
        warn!("j = %d", j);
        warn!("s = %d", s);
        warn!("count = %lu", count);
        warn!("abs = %d", ab);
        warn!("ordered = %d", ordered);
        for (i, v) in good.iter().enumerate() {
            fprintf!(stderr(), "accept[%d] = ", i);
            fprintfval(&stderr(), type_, *v);
            fputc(10, &stderr());
        }
        for (i, v) in bad.iter().enumerate() {
            fprintf!(stderr(), "reject[%d] = ", i);
            fprintfval(&stderr(), type_, *v);
            fputc(10, &stderr());
        }
        warn!("tmin=%f tmax=%f", tmin, tmax);
        warn!("itmin=%d itmax=%d nt=%u", itmin, itmax, nt);
        warn!("Padding %d zeroes", nzeros);
    }
    let mut minin = i32::MAX as i64;
    let mut maxin = i32::MIN as i64;
    let mut maxless = i32::MIN as i64;
    let mut minmore = i32::MAX as i64;
    let mut numin = 0i64;
    loop {
        let mut val = Value::default();
        gethval(&tr, indx, &mut val);
        let isbad = bad.iter().any(|v| valcmp(type_, val, *v) == 0);
        let isgood = good.iter().any(|v| valcmp(type_, val, *v) == 0);
        if ab != 0 {
            val = valtoabs(type_, val);
        }
        let ival = vtoi(type_, val);
        let iv = ival as i64;
        if verbose != 0 {
            if iv < min && iv > maxless {
                maxless = iv;
            }
            if iv > max && iv < minmore {
                minmore = iv;
            }
        }
        if (ordered == 1 && max < iv) || (ordered == -1 && min > iv) {
            break;
        }
        // arm64 division by zero yields zero, leaving the dividend as remainder.
        let rem = ival
            .wrapping_sub(s)
            .checked_rem(j)
            .unwrap_or(ival.wrapping_sub(s));
        if isgood || (min <= iv && iv <= max && rem == 0 && !isbad) {
            // The C program accesses outside tr.data for these windows.
            if itmax >= SU_NFLTS as i32 {
                err!("itmax=%d exceeds trace data buffer", itmax);
            }
            if nzeros != 0 && itmin > ntfirst {
                err!("padding starts before trace data buffer");
            }
            if verbose != 0 {
                if iv < minin {
                    minin = iv;
                }
                if iv > maxin {
                    maxin = iv;
                }
                numin += 1;
            }
            if itmin > 0 {
                for i in itmin..=itmax {
                    tr.data[(i - itmin) as usize] = tr.data[i as usize];
                }
                if seismic {
                    tr.delrt =
                        nint(itmin.wrapping_mul(tr.dt as i32) as f64 * 0.001 + tr.delrt as f64)
                            as i16;
                } else {
                    tr.f1 = itmin as f32 * tr.d1 + tr.f1;
                }
            }
            if nzeros != 0 {
                let start = (ntfirst - itmin) as usize;
                tr.data[start..start + nzeros].fill(0.0);
            }
            tr.ns = nt as u16;
            tr.ntr = 0;
            puttr(&mut tr);
            count = count.wrapping_sub(1);
            if count == 0 {
                break;
            }
        }
        if skip != 0 {
            itr += 1;
            if itr >= ntr || fgettra(&stdin(), &mut tr, itr) == 0 {
                break;
            }
        } else if gettr(&mut tr) == 0 {
            break;
        }
    }
    if verbose != 0 {
        warn!("numinside=%ld", numin);
        warn!("mininside=%ld", minin);
        warn!("maxinside=%ld", maxin);
        warn!("maxless=%ld", maxless);
        warn!("minmore=%ld", minmore);
    }
    exit(CWP_Exit());
}
