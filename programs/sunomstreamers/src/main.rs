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

//! SUNOMSTREAMERS - nominal marine streamer geometry, port of su/main/headers/sunomstreamers.c (Revision 1.00).
//! Credits: Andre Latour, Dec 2022.

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::nint;
use oxwave_header_case::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                            ",
    " SUNOMSTREAMERS - Nominal Geometry and 2D CDPS for Towed-Streamer Marine.   ",
    "                                                                            ",
    " sunomstreamers  >out.su (many parameters, input traces are optional).      ",
    "                                                                            ",
    " This program allows you to define a towed-streamer marine layout with      ",
    " multiple airgun arrays (or other sources) amd multiple towed-streamers.    ",
    " This basically consists of inline and crossline offset positions as well   ",
    " as airgun and streamer identifiers and streamer channel ranges and spacing.",
    " From this definition, output traces are (optionally) updated with offset,  ",
    " sx,sy,gx,gy,cdp, and airgun and streamer sequence (grnofr and grnors keys).",
    " These values are nominal geometry and nominal inline cdps since they are   ",
    " computed only from the layout definition and not any INPUT coordinates.    ",
    "                                                                            ",
    "   **********************************************************************   ",
    "   To output this documentation:  sunomstreamers 2> nomstreamers.txt        ",
    "   **********************************************************************   ",
    "   Also consult: src/demos/Geom3D/Sunomstreamers for layout diagram.        ",
    "   **********************************************************************   ",
    "                                                                            ",
    " Parameters:                                                                ",
    "                                                                            ",
    " shotkey=fldr    Shot Key. This value is multiplied by shotspacing to shift ",
    "                 the shot reference point forward in the inline direction.  ",
    "                 The reference point is the point from where all inline and ",
    "                 crossline offsets are measured. The reference point may    ",
    "                 also be called antenna position or centre-of-boat.         ",
    "                                                                            ",
    " shotspacing=    Inline distance between shotkey values. Must be specified. ",
    "                                                                            ",
    " gunkey=         Key that indicates which airgun array is firing for shots. ",
    "                 The default is the shotkey specified above.              . ",
    "                 This value and gunmod identify which gun array is firing   ",
    "                 (typical acquizition has just port and starboard arrays).  ",
    "           Note: Usually the gunkey is a sequential shot or field record    ",
    "                 number and odd shots are from the port airgun array and    ",
    "                 even shots are from starboard airgun array (or vice-versa) ",
    "                 so the gunid list should be 0,1 (or 1,0).                  ",
    "                                                                            ",
    " gunmod=    Modulus value for gunkey. The remainder is the value to         ",
    "            specify in the gunid list. The default here is the amount       ",
    "            of values specified in the gunid list (typically 2).            ",
    "                                                                            ",
    " gunid=     List of gun array identifier values. Must be specified.         ",
    "            See parameters gunkey and gunmod.                               ",
    "                                                                            ",
    " guncross=  List of gun array crossline offset distances.                   ",
    "            Must specify the same amount of values as the gunid list.       ",
    "            Usually port offsets are negative, starboard are positive.      ",
    "                                                                            ",
    " guninline= List of gun array inline offset distances.                      ",
    "            You can specify either 1 value for all gun arrays or you        ",
    "            can list the same amount of values as the gunid list.           ",
    "            Fore offsets are positive, aft offsets are negative. Since      ",
    "            the gun arrays are typically behind the boat, these values      ",
    "            are usually negative (and also usually all the same).           ",
    "                                                                            ",
    " streamercross=    List of streamer crossline offset distances.             ",
    "                   Must specify a value for each streamer.                  ",
    "                   Usually port offsets are negative, starboard positive.   ",
    "                                                                            ",
    " streamerinline=   List of inline offset distances at start of streamers.   ",
    "                   You can specify either 1 value for all streamers or you  ",
    "                   can list the same amount of values as streamercross.     ",
    "                   Fore offsets are positive, aft are negative. Since the   ",
    "                   streamers are typically behind the boat, these values    ",
    "                   are usually negative (and also usually all the same).    ",
    "                                                                            ",
    " streamerkey=      Streamer Identifier key. If not specified then also the  ",
    "                   streameridnt cannot be specified, and channelstart and   ",
    "                   channelend range must not overlap amoung any streamers.  ",
    "                                                                            ",
    " streameridnt=     List of streamer identifier numbers. If streamerkey is   ",
    "                   not specified, you cannot specify this parameter either. ",
    "                   Otherwise, list same amount of values as streamercross.  ",
    "                                                                            ",
    " channelkey=tracf  Channel Identifier key. Channel numbers also often have  ",
    "                   the same values as trace numbers.                        ",
    "                                                                            ",
    " channelstart=     List of channel numbers at start of streamers.           ",
    "                   You can specify either 1 channelstart and 1 channelend   ",
    "                   or each with the same amount of values as streamercross. ",
    "                   If you do specify only 1 channelstart and 1 channelend   ",
    "                   and streamerkey is specified, the channel range for      ",
    "                   subsequent streamers is identical to the first streamer. ",
    "                   If you do specify only 1 channelstart and 1 channelend   ",
    "                   and streamerkey is NOT specified, you must specify the   ",
    "                   lowest channel range. Other streamers will have the      ",
    "                   same channel amount and their channel numbers will be    ",
    "                   contiguous with and greater than your specified range.   ",
    "                                                                            ",
    " channelend=       List of channel numbers at end of streamers. This value  ",
    "                   can be less than the channelstart of this streamer.      ",
    "                   (Channel numbers can decrease away from streamer start). ",
    "                                                                            ",
    " channelspacing=   List of channel spacing on streamers.                    ",
    "                   You can specify either 1 value for all streamers or you  ",
    "                   can list the same amount of values as streamercross.     ",
    "                                                                            ",
    " cdpspacing=  Inline distance between cdps. Default is the minimum channel  ",
    "              spacing of any streamer divided by 2. Any positive value sets ",
    "              the cdp key to a computed nominal INLINE cdp number. The cdp  ",
    "              is computed only from the layout parameters specified here,   ",
    "              and only from the INLINE offset components of the layout.     ",
    "           =0. Do not reset the cdp key values in the traces.               ",
    "                                                                            ",
    " offpi=2  Set offset key to Pythagorean value computed from the inline and  ",
    "          crossline offsets and other values of the layout definition.      ",
    "      =1  Set offset key to the absolute difference of INLINE components.   ",
    "          (This is also known as the nominal inline offset distance).       ",
    "      =0  Do not reset offset key.                                          ",
    "    Note: This program NEVER computes offset from input trace coordinates.  ",
    "    Note: Option 1 is used in some differential NMO processing situations.  ",
    "                                                                            ",
    " linelocs=1 Set the grnofr,grnors,grnlof,gaps keys as follows:              ",
    "              grnofr to airgun sequence (from 1 to amount of airguns)       ",
    "              grnors to streamer sequence (from 1 to amount of streamers)   ",
    "              grnlof to nominal inline cdp computed at shot location        ",
    "              gaps   to nominal inline cdp computed at channel location     ",
    "         =0 Do not reset the grnofr,grnors,grnlof,gaps keys.                ",
    "      Note: grnofr multiplied by grnors forms a subsurface line number      ",
    "            which is useful in some data processing situations.             ",
    "                                                                            ",
    " scalco=10  Multiply coordinates by this power of 10 (1,10,100...)          ",
    "            before putting them in traces. Default is 10 which means        ",
    "            that sx,sy,gx,gy are multiplied by 10. The actual value of      ",
    "            scalco in the traces is therefore set to -10 (meaning these     ",
    "            values need to be divided by 10). Also counit key is set to 1.  ",
    "      Note: The x coordinates have 100000 added to them before output and   ",
    "            the y coordinates have 10000 added to them before output.       ",
    "            The 100000 and 10000 values are added just to avoid having      ",
    "            negative coordinates for the output sx,sy,gx,gy key values.     ",
    "       =0   Do not reset sx,sy,gx,gy,scalco,counit. There are situations    ",
    "            when you want to compute nominal inline cdps and other values   ",
    "            but traces already have coordinates (which you wish to retain). ",
    "                                                                            ",
    " verbose=  Default is not to print expanded gun and streamer information.   ",
    "        =1 Print expanded gun and streamer information. All airgun array    ",
    "           and streamer information is printed after expanding all defaults.",
    "                                                                            ",
    " create=     Default is not to create output traces.                        ",
    "       =n    Create n shots (all traces of each shot).                      ",
    "       Note: Trying to also input a trace file error-halts unless spikes=n. ",
    "                                                                            ",
    " firstshot=1 This is the shotkey value of the first shot. The shotkey       ",
    "             values increment by 1 for each output shot.                    ",
    "                                                                            ",
    " spikes=     This option will error-halt if create= is not specified.       ",
    "             By default, create= outputs only 1 zero-sample at 4 ms.        ",
    "             This parameter accepts a list of time,amplitude pairs.         ",
    "                                                                            ",
    "               Example: spikes=4,0.001,800,1000,1500,-2000,3000,0           ",
    "                                                                            ",
    "             The maximum time sets the trace length (3000 ms. above).       ",
    "             The first pair in the list is special, it specifies the        ",
    "             sample interval (in ms.) and base-amplitude.                   ",
    "             In this example, the sample interval is 4, and the             ",
    "             base-amplitude is 0.001. The base-amplitude is the value       ",
    "             samples are set to if they are not spiked by remainder of      ",
    "             the list. (For geometry checking and other QC and survey       ",
    "             design reasons, base-amplitude of 0.0 is not a nice value.     ",
    "             For instance, 0.0 makes it difficult to check mutes).          ",
    "             All times should be whole multiples of the sample interval     ",
    "             but will be rounded to nearest sample time if not.             ",
    "             Times do not have to increase but cannot be negative.          ",
    "                                                                            ",
    "             -or-                                                           ",
    " spikes=n    Where n is the sequential number of an input trace. The        ",
    "             sample values of this trace will be copied to all created      ",
    "             output traces. The intention is to allow input reflectivity    ",
    "             type of trace, and duplicate the samples for all output.       ",
    "                                                                            ",
    " ------------------------------------------------------------------------   ",
    "                                                                            ",
];

fn main() {
    initargs(&args());
    requestdoc(1, SDOC);
    let mut verbose = 0;
    getparint("verbose", &mut verbose);
    if verbose < 0 || verbose > 1 {
        err!("error: verbose= is out of range.");
    }
    let mut offpi = 2;
    getparint("offpi", &mut offpi);
    if offpi < 0 || offpi > 2 {
        err!("error: offpi= is out of range.");
    }
    let mut linelocs = 1;
    getparint("linelocs", &mut linelocs);
    if linelocs < 0 || linelocs > 1 {
        err!("error: linelocs= is out of range.");
    }
    let mut scalco: i32 = 10;
    getparint("scalco", &mut scalco);
    let mut scale = 0.0;
    if scalco != 0 {
        if scalco == -1 {
            scalco = 1;
        }
        if ![1, 10, 100, 1000, 10000, 100000, 1000000, 10000000].contains(&scalco.wrapping_abs()) {
            err!("**** Error: scalco= must be signed powers of 10 (1,10,100...-10,-100,...)");
        }
        scale = if scalco > 0 {
            scalco as f64
        } else {
            -1.0 / scalco as f64
        };
        scalco = -scalco;
    }
    let mut shotkey = String::from("fldr");
    getparstring("shotkey", &mut shotkey);
    let shot = GetCase(&shotkey);
    if shot < 1 {
        err!("error: shotkey= is not a recognized header key name");
    }
    let mut spacing = 0.0;
    if getpardouble("shotspacing", &mut spacing) == 0 {
        err!("error: shot spacing must be specified.");
    }
    if spacing <= 0.0 {
        err!("error: shot spacing must be greater than 0.");
    }
    let ns = countparval("guncross") as usize;
    if ns == 0 {
        err!("**** Error: at least 1 gun array crossline offset must be specified.");
    }
    if countparval("guninline") != 1 && countparval("guninline") != ns as i32 {
        err!("**** Error: must specify 1 gun array inline offset -OR- same amount as gun array crossline offsets.");
    }
    if countparval("gunid") != ns as i32 {
        err!("**** Error: different amount of gun array crossline offsets and gun identifiers.");
    }
    let mut gunmod = ns as i32;
    getparint("gunmod", &mut gunmod);
    if gunmod < 1 {
        err!("**** Error: gunmod must be greater than 0.");
    }
    let mut key = String::new();
    let gun = if getparstring("gunkey", &mut key) == 0 {
        shot
    } else {
        let k = GetCase(&key);
        if k < 1 {
            err!("error: gunkey= is not a recognized header key name");
        }
        k
    };
    let mut sx = vec![0.0; ns];
    let mut si = vec![0.0; ns];
    let mut sg = vec![0; ns + 1];
    getpardouble("guncross", &mut sx);
    getpardouble("guninline", &mut si);
    if countparval("guninline") == 1 {
        for i in 1..ns {
            si[i] = si[0];
        }
    }
    getparint("gunid", &mut sg);
    let mut back = 0.0;
    for i in 0..ns {
        if sg[i] >= gunmod || sg[i] < 0 {
            err!(
                "error: specified gunid=%d cannot occur with gunmod=%d",
                sg[i],
                gunmod
            );
        }
        for j in i + 1..ns {
            if sg[i] == sg[j] {
                err!("error: gunid=%d has been specified at least twice.", sg[i]);
            }
            if sx[i] == sx[j] {
                warn!("Unusual: some gun arrays have the same crossline offset.");
            }
            if si[i] != si[j] {
                warn!("Unusual: some gun arrays have different inline offsets.");
            }
        }
        if si[i] < back {
            back = si[i];
        }
    }
    sg[ns] = -2147483645;
    key = String::from("tracf");
    getparstring("channelkey", &mut key);
    let chan = GetCase(&key);
    if chan < 1 {
        err!("error: channelkey= is not a recognized header key name");
    }
    let nt = countparval("streamercross") as usize;
    if nt == 0 {
        err!("**** Error: at least 1 streamer crossline offset must be specified.");
    }
    if countparval("streamerinline") != 1 && countparval("streamerinline") != nt as i32 {
        err!("**** Error: must specify 1 streamer inline offset -OR- same amount as streamer crossline offsets.");
    }
    let useid = countparval("streameridnt") > 0;
    if useid && countparval("streameridnt") != nt as i32 {
        err!("**** Error: different amount of streamer crosslines and streamer idnts specified.");
    }
    let mut idkey = 0;
    key.clear();
    if useid {
        if getparstring("streamerkey", &mut key) == 0 {
            err!("**** Error: streameridnt and streamerkey must both be specified (or neither).");
        }
        idkey = GetCase(&key);
        if idkey < 1 {
            err!("error: streamerkey= is not a recognized header key name");
        }
    } else if getparstring("idntkey", &mut key) != 0 {
        err!("**** Error: streameridnt and streamerkey must both be specified (or neither).");
    }
    if countparval("channelstart") != 1 || countparval("channelend") != 1 {
        if countparval("channelstart") != nt as i32 {
            err!(
                "**** Error: different amount of streamer crosslines and channel starts specified."
            );
        }
        if countparval("channelend") != nt as i32 {
            err!("**** Error: different amount of streamer crosslines and channel ends specified.");
        }
    }
    if countparval("channelspacing") != 1 && countparval("channelspacing") != nt as i32 {
        err!("**** Error: must specify 1 streamer channel spacing -OR- same amount as streamer crossline offsets.");
    }
    let mut tx = vec![0.0; nt];
    let mut ti = vec![0.0; nt];
    let mut ids = vec![0; nt + 1];
    let mut beg = vec![0i32; nt + 1];
    let mut end = vec![0i32; nt + 1];
    let mut dist = vec![0.0; nt];
    getpardouble("streamercross", &mut tx);
    getpardouble("streamerinline", &mut ti);
    if countparval("streamerinline") == 1 {
        for i in 1..nt {
            ti[i] = ti[0];
        }
    }
    if useid {
        getparint("streameridnt", &mut ids);
        for i in 0..nt {
            for j in i + 1..nt {
                if ids[i] == ids[j] {
                    err!(
                        "error: streameridnt=%d has been specified at least twice.",
                        ids[i]
                    );
                }
            }
        }
    } else {
        for x in &mut ids[..nt] {
            *x = -9999;
        }
    }
    getpardouble("channelspacing", &mut dist);
    if countparval("channelspacing") == 1 {
        for i in 1..nt {
            dist[i] = dist[0];
        }
    }
    getparint("channelstart", &mut beg);
    getparint("channelend", &mut end);
    if countparval("channelstart") == 1 && countparval("channelend") == 1 {
        for i in 1..nt {
            if useid {
                beg[i] = beg[0];
                end[i] = end[0];
            } else {
                let d = (end[0] - beg[0]).abs() + 1;
                beg[i] = beg[0] + i as i32 * d;
                end[i] = end[0] + i as i32 * d;
            }
        }
    }
    for i in 0..nt {
        for j in i + 1..nt {
            if !useid {
                if beg[i] >= beg[j] && beg[i] <= end[j] {
                    err!("error: channelstart=%d for streamer=%d is in channel range of streamer=%d.",beg[i],i+1,j+1);
                }
                if end[i] >= beg[j] && end[i] <= end[j] {
                    err!(
                        "error: channelend=%d for streamer=%d is in channel range of streamer=%d.",
                        end[i],
                        i + 1,
                        j + 1
                    );
                }
            }
            if tx[i] == tx[j] {
                warn!("Unusual: some streamers have the same crossline offset.");
            }
            if ti[i] != ti[j] {
                warn!("Unusual: some streamers have different inline offsets.");
            }
            if dist[i] != dist[j] {
                warn!("Unusual: some streamers have different channelspacings.");
            }
        }
    }
    let mut dmin = dist[0];
    for i in 0..nt {
        if dist[i] <= 0.0 {
            err!("**** Error: every channelspacing must be positive.");
        }
        if dist[i] < dmin {
            dmin = dist[i];
        }
        if beg[i] > end[i] {
            ti[i] -= (beg[i] - end[i]) as f64 * dist[i];
            std::mem::swap(&mut beg[i], &mut end[i]);
            dist[i] = -dist[i];
        }
        if dist[i] > 0.0 {
            let b = ti[i] + dist[i] * (beg[i] - end[i]) as f64;
            if b < back {
                back = b;
            }
        } else if ti[i] < back {
            back = ti[i];
        }
    }
    ids[nt] = -2147483645;
    beg[nt] = 2147483645;
    end[nt] = -2147483645;
    let mut cdp = dmin * 0.5;
    getpardouble("cdpspacing", &mut cdp);
    if cdp < 0.0 {
        err!("error: cdp spacing parameter cannot be negative.");
    }
    if cdp != 0.0 {
        if spacing % cdp != 0.0 {
            warn!("Unusual: shot spacing is not a multiple of cdp spacing.");
        }
        for i in 0..nt {
            if dist[i] % cdp != 0.0 {
                warn!(
                    "Unusual: channel spacing is not a multiple of cdp spacing on streamer=%d",
                    i + 1
                );
            }
        }
    }
    if verbose > 0 {
        warn!(" -----------------------------------------------------------------------------");
        warn!("     gunmod=%d", gunmod);
        warn!("      gunid     guninline  guncross   output grnofr value");
        for i in 0..ns {
            warn!("%10d %10g %10g %15d ", sg[i], si[i], sx[i], i + 1);
        }
        warn!(" --- ");
        warn!("   streameridnt  streamerinline streamercross  channelstart channelend channelspacing  output grnors value");
        let mut n = 0;
        for i in 0..nt {
            warn!(
                "%10d %15g %12g %16d %10d %12g %15d ",
                ids[i],
                ti[i],
                tx[i],
                beg[i],
                end[i],
                dist[i],
                i + 1
            );
            if dist[i] < 0.0 {
                n += 1;
            }
        }
        if n > 0 {
            warn!(" There are %d negative channelspacings. Which means channel range was exchanged and streamerinline adjusted.",n);
        }
        warn!(" -----------------------------------------------------------------------------");
    }
    let mut create = 0;
    getparint("create", &mut create);
    if create < 0 {
        err!("**** Error create= cannot be negative.");
    }
    let mut fake = 1;
    getparint("firstshot", &mut fake);
    let nsp = countparval("spikes") as usize;
    let mut spikes = vec![0.0f32; nsp.max(1)];
    let mut replicate = 0;
    if nsp > 0 {
        if create < 1 {
            err!("**** Error: spikes= can only be specified if create= is also specified.");
        }
        if nsp > 1 {
            if nsp % 2 == 1 {
                err!("**** Error: spikes= must have an even number of values (time amplitude pairs).");
            }
            getparfloat("spikes", &mut spikes);
        } else {
            getparint("spikes", &mut replicate);
            if replicate < 1 {
                err!("**** Error: For spikes=n the n must be +integer (the trace to input).");
            }
        }
    }
    let mut last = 0;
    if create > 0 {
        if !isatty(0) {
            if replicate < 1 {
                err!("**** Error: Input trace file not allowed when create= (unless spikes=n).");
            }
        } else if replicate > 0 {
            err!("**** Error: Input trace file must be specified when spikes=n.");
        }
        if cdp == 0.0 {
            warn!("Caution: cdpspacing=0 so CDP key will not be set in these created traces.");
        }
        if offpi == 0 {
            warn!("Caution: offpi=0 so OFFSET key will not be set in these created traces.");
        }
        fake -= 1;
        last = create + fake;
    } else if isatty(0) {
        err!("**** Error: Input trace file must be specified when not create= ");
    }
    checkpars();
    let mut tr = Trace::new();
    if create < 1 || replicate > 0 {
        if gettr(&mut tr) == 0 {
            err!("can't read first trace");
        }
        for _ in 1..replicate {
            if gettr(&mut tr) == 0 {
                err!("can't read enough traces for spikes=n option.");
            }
        }
    } else if nsp == 0 {
        tr.ns = 1;
        tr.dt = 4000;
        tr.data[0] = 0.0;
    } else {
        let mut tmax = 0.0f32;
        for i in (2..nsp).step_by(2) {
            if spikes[i] > tmax {
                tmax = spikes[i];
            }
        }
        let n = nint((tmax / spikes[0]) as f64);
        // The C trace buffer has a fixed capacity; larger counts would write past it.
        if n > tr.data.len() as i32 {
            err!("Error: Spike sample count exceeds trace capacity.");
        }
        // Validate every index before writing either the base samples or spikes.
        for i in 0..nsp / 2 {
            let j = nint((spikes[2 * i] / spikes[0]) as f64).wrapping_sub(1);
            if j < 0 {
                err!("Error: Spike times less than first sample are not permitted.");
            }
            if j as usize >= tr.data.len() {
                err!("Error: Spike sample index exceeds trace capacity.");
            }
        }
        tr.ns = n as u16;
        tr.dt = nint(1000.0 * spikes[0] as f64) as u16;
        for i in 0..n {
            tr.data[i as usize] = spikes[1];
        }
        for i in 0..nsp / 2 {
            let j = nint((spikes[2 * i] / spikes[0]) as f64).wrapping_sub(1);
            tr.data[j as usize] = spikes[2 * i + 1];
        }
    }
    let mut source = ns;
    let mut streamer = nt;
    let mut fg = -1;
    let mut fc = 2147483645;
    let mut fs = nt - 1;
    let mut idvalue = 0;
    let mut count = 0;
    loop {
        if create > 0 {
            fc += 1;
            if fc > end[fs] {
                fs += 1;
                if fs == nt {
                    fake += 1;
                    if fake > last {
                        break;
                    }
                    fs = 0;
                    tohead(&mut tr, shot, fake as f64);
                    if gun != shot {
                        fg += 1;
                        if fg == ns as i32 {
                            fg = 0;
                        }
                        tohead(&mut tr, gun, sg[fg as usize] as f64);
                    }
                }
                fc = beg[fs];
                if useid {
                    tohead(&mut tr, idkey, ids[fs] as f64);
                }
            }
            tohead(&mut tr, chan, fc as f64);
        }
        let g = fromhead(&tr, gun) as i32 % gunmod;
        if g != sg[source] {
            source = ns;
            for i in 0..ns {
                if g == sg[i] {
                    source = i;
                    break;
                }
            }
            if source == ns {
                err!("**** Error: source gun array ident is not any specified gunid. Shot=%d Idnt=%d",tr.fldr,idvalue);
            }
        }
        if useid {
            idvalue = fromhead(&tr, idkey) as i32;
            if idvalue != ids[streamer] {
                streamer = nt;
                for i in 0..nt {
                    if idvalue == ids[i] {
                        streamer = i;
                        break;
                    }
                }
                if streamer == nt {
                    err!(
                        "**** Error: streamer idnt not in any specified streamer. Shot=%d Idnt=%d",
                        tr.fldr,
                        idvalue
                    );
                }
            }
        }
        let ch = fromhead(&tr, chan) as i32;
        if ch > end[streamer] || ch < beg[streamer] {
            streamer = nt;
            for i in 0..nt {
                if ch <= end[i] && ch >= beg[i] {
                    streamer = i;
                    break;
                }
            }
            if streamer == nt {
                err!("**** Error: channel not in range of any specified streamer. Shot=%d Channel=%d",tr.fldr,ch);
            }
        }
        let reference = fromhead(&tr, shot) * spacing - back;
        let x = si[source] + reference;
        let y = sx[source];
        let gx = ti[streamer] + reference + dist[streamer] * (beg[streamer] - ch) as f64;
        let gy = tx[streamer];
        if offpi > 1 {
            tr.offset = ((x - gx) * (x - gx) + (y - gy) * (y - gy))
                .sqrt()
                .round_ties_even() as i32;
        } else if offpi == 1 {
            tr.offset = (x - gx).abs().round_ties_even() as i32;
        }
        if cdp != 0.0 {
            tr.cdp = ((x + gx) * 0.5 / cdp).round_ties_even() as i32;
        }
        if linelocs > 0 {
            tr.grnofr = (source + 1) as i16;
            tr.grnors = (streamer + 1) as i16;
            tr.grnlof = (x / cdp).round_ties_even() as i64 as i16;
            tr.gaps = (gx / cdp).round_ties_even() as i64 as i16;
        }
        if scalco != 0 {
            tr.sx = ((x + 100000.0) * scale).round_ties_even() as i32;
            tr.sy = ((y + 10000.0) * scale).round_ties_even() as i32;
            tr.gx = ((gx + 100000.0) * scale).round_ties_even() as i32;
            tr.gy = ((gy + 10000.0) * scale).round_ties_even() as i32;
            tr.scalco = scalco as i16;
            tr.counit = 1;
        }
        puttr(&mut tr);
        count += 1;
        if create < 1 && gettr(&mut tr) == 0 {
            break;
        }
    }
    warn!("Number of traces output= %d ", count);
    exit(CWP_Exit());
}
