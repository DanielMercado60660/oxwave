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

//! Layout and duplicate checks, in their original order from su/main/headers/sutoolcsv.c.
//! Credits: Andre Latour; started from sugeom by Fernando Roxo <petro@roxo.org>.
#![allow(non_snake_case, unused_parens, unused_assignments, unused_mut)]
use super::{longt, PointView};
use oxwave_core::prelude::*;

pub(crate) fn check(
    RecInfo: &[PointView<'_>],
    mapx: &[i32; 10],
    match_: &[String],
    dtolh: f64,
    dtol: f64,
) {
    let numR = RecInfo.len();
    let num_to_sort_by = match_.len();
    // Channel ranges are grouped by every match key except the last (_cf).
    // Overlapping ranges with different increments may represent multiple
    // components at the same receiver, so both increments and residues matter.
    if (mapx[3] > -1) {
        let mut lapover: i32 = 0;
        let mut l1verr: i32 = 0;
        let mut l2verr: i32 = 0;
        let mut l3verr: i32 = 0;
        let mut l4verr: i32 = 0;
        let mut l5verr: i32 = 0;
        let mut l6verr: i32 = 0;
        let mut l7verr: i32 = 0;
        let mut l8verr: i32 = 0;
        let mut l1same: i32 = 0;
        let mut l2same: i32 = 0;
        let mut l3same: i32 = 0;
        let mut l4same: i32 = 0;
        let mut l5same: i32 = 0;
        let mut l6same: i32 = 0;
        let mut l7same: i32 = 0;
        let mut l8same: i32 = 0;
        let mut ntop_grnofr: i64 = 0;
        let mut ntop_grnlof: i64 = 0;

        let num_of_others = num_to_sort_by - 1;

        let mut ntop: usize = 0;
        let mut npchan: i32 = 0;
        let mut npsegs: i32 = 0;
        let mut rpinc: f64 = 0.;

        for n in 1..=numR {
            if (n == numR
                || RecInfo[ntop].lfield[..num_of_others] != RecInfo[n].lfield[..num_of_others])
            {
                l1same = 0;
                l2same = 0;
                l3same = 0;
                l4same = 0;
                l5same = 0;
                l6same = 0;
                l7same = 0;
                l8same = 0;

                let mut nchan: i32 = 0;
                let mut nsegs: i32 = 0;

                for m in ntop..n {
                    let mut mcf: i64 = longt(RecInfo[m].dfield[mapx[0] as usize], dtolh, dtol);
                    let mut mct: i64 = longt(RecInfo[m].dfield[mapx[1] as usize], dtolh, dtol);
                    let mut mci: i64 = longt(RecInfo[m].dfield[mapx[2] as usize], dtolh, dtol);
                    if (mci == 0) {
                        mci = 1;
                    }

                    nchan = nchan.wrapping_add(((mct - mcf) / mci + 1) as i32);
                    nsegs += 1;

                    let mut rinc: f64 = 1.;
                    if (mct != mcf) {
                        rinc = (RecInfo[m].dfield[mapx[3] as usize]
                            - RecInfo[m].dfield[mapx[4] as usize])
                            / (((mct - mcf) / mci) as f64);
                    }

                    if (ntop > 0 && l5same == 0 && (rinc - rpinc).abs() > dtolh * 2.) {
                        l5same = 1;
                        l5verr = l5verr + 1;
                        if (l5verr < 4) {
                            let mut j: i32 = (RecInfo[m].lfield[0] as f64 / dtol + 0.5) as i32;
                            warn!(
                                "Warning: Receiver points per channel changed at %s= %d ",
                                match_[0], j
                            );
                            if (l5verr == 3) {
                                warn!("Have 3 Receiver-points-per-channel changed warnings, no more will be printed.");
                            }
                        }
                    }
                    rpinc = rinc;

                    if (l1same == 0 && mcf % mci != mct % mci) {
                        l1same = 1;
                        l1verr = l1verr + 1;
                        if (l1verr < 4) {
                            let mut j: i32 = (RecInfo[m].lfield[0] as f64 / dtol + 0.5) as i32;
                            warn!(
                                "Error: Layout ends do not conform to same increment at %s= %d ",
                                match_[0], j
                            );
                            if (l1verr == 3) {
                                warn!("Have 3 Layout-ends-do-not-conform-to-same-increment errors, no more will be printed.");
                            }
                        }
                    }

                    if (mapx[7] > -1) {
                        if (m == ntop) {
                            ntop_grnofr = longt(RecInfo[m].dfield[mapx[7] as usize], dtolh, dtol);
                        } else {
                            if (l7same == 0
                                && ntop_grnofr
                                    != longt(RecInfo[m].dfield[mapx[7] as usize], dtolh, dtol))
                            {
                                l7same = 1;
                                l7verr = l7verr + 1;
                                if (l7verr < 4) {
                                    let mut j: i32 =
                                        (RecInfo[m].lfield[0] as f64 / dtol + 0.5) as i32;
                                    warn!("Error: The grnofr values are different for same shot at %s= %d ",match_[0],j);
                                    if (l7verr == 3) {
                                        warn!("Have 3 grnofr-values-are-different-for-same-shot, no more will be printed.");
                                    }
                                }
                            }
                        }
                    }

                    if (mapx[8] > -1) {
                        if (m == ntop) {
                            ntop_grnlof = longt(RecInfo[m].dfield[mapx[8] as usize], dtolh, dtol);
                        } else {
                            if (l8same == 0
                                && ntop_grnlof
                                    != longt(RecInfo[m].dfield[mapx[8] as usize], dtolh, dtol))
                            {
                                l8same = 1;
                                l8verr = l8verr + 1;
                                if (l8verr < 4) {
                                    let mut j: i32 =
                                        (RecInfo[m].lfield[0] as f64 / dtol + 0.5) as i32;
                                    warn!("Error: The grnlof values are different for same shot at %s= %d ",match_[0],j);
                                    if (l8verr == 3) {
                                        warn!("Have 3 grnlof-values-are-different-for-same-shot, no more will be printed.");
                                    }
                                }
                            }
                        }
                    }

                    for i in m + 1..n - 1 {
                        let mut icf: i64 = longt(RecInfo[i].dfield[mapx[0] as usize], dtolh, dtol);
                        let mut ici: i64 = longt(RecInfo[i].dfield[mapx[2] as usize], dtolh, dtol);
                        if (icf > mct) {
                            break;
                        }
                        if (lapover < 1) {
                            let mut j: i32 = (RecInfo[i].lfield[0] as f64 / dtol + 0.5) as i32;
                            warn!("Warning: Overlapping channel range at %s= %d  Unusual, but not always an error.",match_[0],j);
                        }
                        lapover += 1;

                        if (l2same == 0 && ici != mci) {
                            l2same = 1;
                            l2verr = l2verr + 1;
                            if (l2verr < 4) {
                                let mut j: i32 = (RecInfo[i].lfield[0] as f64 / dtol + 0.5) as i32;
                                warn!(
                                    "Error: Different increments in overlapping layout at %s= %d ",
                                    match_[0], j
                                );
                                if (l2verr == 3) {
                                    warn!("Have 3 Different-increments-in-overlapping layout errors, no more will be printed.");
                                }
                            }
                        }

                        // ARM64 integer remainder with a zero divisor yields the
                        // dividend (SDIV gives zero); avoid a Rust panic here.
                        if (l3same == 0 && mcf % mci == icf.checked_rem(ici).unwrap_or(icf)) {
                            l3same = 1;
                            l3verr = l3verr + 1;
                            if (l3verr < 4) {
                                let mut j: i32 = (RecInfo[i].lfield[0] as f64 / dtol + 0.5) as i32;
                                warn!(
                                    "Error: Overlapping layout hits same channels at %s= %d ",
                                    match_[0], j
                                );
                                if (l3verr == 3) {
                                    warn!("Have 3 Overlapping-layout-hits-same-channels errors, no more will be printed.");
                                }
                            }
                        }
                    }
                }

                if (ntop > 0 && l4same == 0 && nchan != npchan) {
                    l4same = 1;
                    l4verr = l4verr + 1;
                    if (l4verr < 4) {
                        let mut j: i32 = (RecInfo[ntop].lfield[0] as f64 / dtol + 0.5) as i32;
                        warn!(
                            "Warning: Number of channels in layout changed at %s= %d ",
                            match_[0], j
                        );
                        if (l4verr == 3) {
                            warn!("Have 3 Number-of-channels-in-layout-changed warnings, no more will be printed.");
                        }
                    }
                }
                npchan = nchan;

                if (ntop > 0 && l6same == 0 && nsegs != npsegs) {
                    l6same = 1;
                    l6verr = l6verr + 1;
                    if (l6verr < 4) {
                        let mut j: i32 = (RecInfo[ntop].lfield[0] as f64 / dtol + 0.5) as i32;
                        if (nsegs == npsegs * 2 || nsegs * 2 == npsegs) {
                            warn!("Warning: Number of segments in layout changed by ratio 2 at %s= %d (see unrepeat=).",match_[0],j);
                        } else {
                            warn!(
                                "Warning: Number of segments in layout changed at %s= %d ",
                                match_[0], j
                            );
                        }
                        if (l6verr == 3) {
                            warn!("Have 3 Number-of-segments-in-layout-changed warnings, no more will be printed.");
                        }
                    }
                }
                npsegs = nsegs;

                ntop = n;
            }
        }

        if (lapover > 0) {
            warn!(
                "Total Overlapping-channel-ranges in layouts:        %d (very unusual)",
                lapover
            );
        }
        if (l2verr > 0) {
            warn!("Total Different-increments-in-overlapping layout:   %d (will error-halt SUGEOMCSV)",l2verr);
        }
        if (l3verr > 0) {
            warn!("Total Overlapping-layout-hits-same-channels:        %d (will error-halt SUGEOMCSV)",l3verr);
        }
        if (l1verr > 0) {
            warn!("Total Layout-ends-do-not-conform-to-same-increment: %d (will error-halt SUGEOMCSV)",l1verr);
        }
        if (l4verr > 0) {
            warn!(
                "Total Number-of-channels-in-layout-changed:         %d (unusual)",
                l4verr
            );
        }
        if (l5verr > 0) {
            warn!(
                "Total Receiver-points-per-channel changed:          %d (very unusual)",
                l5verr
            );
        }
        if (l6verr > 0) {
            warn!(
                "Total Number-of-segments-in-layout-changed:         %d (unusual)",
                l6verr
            );
        }
        if (l7verr > 0) {
            warn!("Total grnofr-values-are-different-for-same-shot:    %d (will error-halt SUGEOMCSV)",l7verr);
        }
        if (l8verr > 0) {
            warn!("Total grnlof-values-are-different-for-same-shot:    %d (will error-halt SUGEOMCSV)",l8verr);
        }
    } else {
        let mut l7verr: i32 = 0;
        for n in 1..numR {
            if (RecInfo[n].lfield == RecInfo[n - 1].lfield) {
                l7verr = l7verr + 1;
                if (l7verr < 4) {
                    if (num_to_sort_by == 1) {
                        let mut j: i32 = (RecInfo[n].lfield[0] as f64 / dtol + 0.5) as i32;
                        warn!(
                            "Error: Records have duplicate match values of %s=%d",
                            match_[0], j
                        );
                    } else if (num_to_sort_by == 2) {
                        let mut j: i32 = (RecInfo[n].lfield[0] as f64 / dtol + 0.5) as i32;
                        let mut k: i32 = (RecInfo[n].lfield[1] as f64 / dtol + 0.5) as i32;
                        warn!(
                            "Error: Records have duplicate match values of %s=%d   %s=%d",
                            match_[0], j, match_[1], k
                        );
                    } else if (num_to_sort_by == 3) {
                        let mut j: i32 = (RecInfo[n].lfield[0] as f64 / dtol + 0.5) as i32;
                        let mut k: i32 = (RecInfo[n].lfield[1] as f64 / dtol + 0.5) as i32;
                        let mut l: i32 = (RecInfo[n].lfield[2] as f64 / dtol + 0.5) as i32;
                        warn!(
                            "Error: Records have duplicate match values of %s=%d   %s=%d   %s=%d",
                            match_[0], j, match_[1], k, match_[2], l
                        );
                    } else {
                        let mut j: i32 = (RecInfo[n].lfield[0] as f64 / dtol + 0.5) as i32;
                        let mut k: i32 = (RecInfo[n].lfield[1] as f64 / dtol + 0.5) as i32;
                        let mut l: i32 = (RecInfo[n].lfield[2] as f64 / dtol + 0.5) as i32;
                        let mut m: i32 = (RecInfo[n].lfield[3] as f64 / dtol + 0.5) as i32;
                        warn!("Error: Records have duplicate match values of %s=%d   %s=%d   %s=%d   %s=%d",
                  match_[0],j,match_[1],k,match_[2],l,match_[3],m);
                    }
                    if (l7verr == 3) {
                        warn!("Have 3 Records-have-duplicate-match-values errors, no more will be printed.");
                    }
                }
            }
        }

        if (l7verr > 0) {
            warn!("Total  errors  for Records-have-duplicate-match-values: %d (will error-halt SUGEOMCSV)",l7verr);
        }

        if (num_to_sort_by == 1 || num_to_sort_by == 2) {
            let mut jn: i32 = 0;
            let mut jp: i32 = 0;
            let mut jd: i32 = 0;
            let mut jc: i32 = 1;
            let mut je: i32 = -1;
            let mut kn: i32 = 0;
            let mut kp: i32 = 0;
            let mut kd: i32 = 0;
            let mut kc: i32 = -1;

            for n in 0..numR {
                if (num_to_sort_by == 2) {
                    jn = (RecInfo[n].lfield[0] as f64 / dtol + 0.5) as i32;
                    kn = (RecInfo[n].lfield[1] as f64 / dtol + 0.5) as i32;
                } else {
                    kn = (RecInfo[n].lfield[0] as f64 / dtol + 0.5) as i32;
                }

                if (n == 0) {
                    jp = jn;
                    kp = kn;
                }

                if (jn != jp) {
                    jc += 1;
                    if (jd != jn - jp) {
                        je += 1;
                    }
                    jd = jn - jp;
                } else {
                    if (kn - kp != kd) {
                        kc += 1;
                        kd = kn - kp;
                    }
                }
                jp = jn;
                kp = kn;
            }

            if (num_to_sort_by == 2) {
                warn!(
                    "Note: There are: %d sets of %s values (lines?).",
                    jc, match_[0]
                );
                if (je > 1) {
                    warn!("Warning: There are: %d irregular %s increments between the sets (missing lines?).",(je+1)/2,match_[0]);
                }
                if (kc > 1) {
                    warn!("Warning: There are: %d irregular %s increments within the lines (missing points?).",(kc+1)/2,match_[1]);
                }
            } else {
                if (kc > 1) {
                    warn!("Warning: There are: %d irregular %s increments within the line (missing points?).",(kc+1)/2,match_[0]);
                }
            }
        }
    }
}
