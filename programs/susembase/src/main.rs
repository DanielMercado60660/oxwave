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

//! SUSEMBASE - Semblance Based Structure-Oriented Interpolate and Smooth.
//! Port of su/main/noise/susembase.c (Revision 1.1).
//! Credits: Author: Andre Latour, Feb 2023.

use oxwave_core::prelude::*;
use oxwave_header_case::*;
use oxwave_sinc_interp::ints8r;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                            ",
    " SUSEMBASE - Semblance Based Structure-Oriented Interpolate and Smooth.     ",
    "                                                                            ",
    " This program performs semblance scanning to derive (apparent) dips in the  ",
    " direction of cdp numbers (or other key values). It then uses those dips to ",
    " align trace sample values while performing interpolation and/or smoothing. ",
    "                                                                            ",
    " susembase    <in.su >out.su  (other parameters)                            ",
    "                                                                            ",
    " Parameter overview:                                                        ",
    "                                                                            ",
    " lockey=cdp   Trace key that indicates trace location. Input traces must be ",
    "              in increasing-OR-EQUAL order of these values. This key value  ",
    "              is reset on output. See outputinc and outputalign parameters. ",
    "        Note: The input locations must be in increasing-OR-EQUAL order.     ",
    "              So, traces can have duplicate lockeys provided the numrolling ",
    "              parameter is set large enough. Whether these duplicates have  ",
    "              serious effects on the output traces varies with situation.   ",
    "                                                                            ",
    " groupkey=    Trace key that indicates the end of a set of lockey values.   ",
    "              If specified, traces must be ordered by lockey and groupkey.  ",
    "              If a 3D stacked dataset is input, this parameter should often ",
    "              be either igi or igc, depending on which direction is being   ",
    "              processed in this run (see program subincvs).                 ",
    "        Note: The keys for susort are listed in reverse order with respect  ",
    "              to lockey and groupkey here.                                  ",
    "                                                                            ",
    " outputinc=1.0    Increment for output locations. Example: if a cdp stack   ",
    "                  is input and cdp numbers increment by 1 then setting this ",
    "                  parameter to 0.5 will output about twice as many traces.  ",
    "            Note: Since lockey is (almost always) integer, it cannot hold   ",
    "                  decimal fractions, so if this outputinc parameter is      ",
    "                  between 0 and 1, the lockey output values are inversed:   ",
    "                   (input_values - outputalign) / outputinc + outputshift   ",
    "                  If outputinc is greater than 1 then the lockey output     ",
    "                  values are:    (input_values - outputalign) + outputshift ",
    "                                                                            ",
    " outputalign=0.0  Alignment for output locations. Example: if a cdp stack   ",
    "                  is input and cdp numbers increment by 1 and outputinc=1.0 ",
    "                  then setting this parameter to 0.25 will create output    ",
    "                  traces at input cdp numbers plus 0.25                     ",
    "                                                                            ",
    " outputshift=0.0  Add this value to output location numbers. If lockey is   ",
    "                  cdp and first input cdp number is small, traces may be    ",
    "                  extrapolated to negative cdp numbers. This may be an      ",
    "                  issue for some situations, so you can add this constant.  ",
    "                                                                            ",
    " locmaxscan=3 Maximum Relative lockey to include in scanning for dips.      ",
    "              Must be greater than or equal to 0. This program scans traces ",
    "              from -locmaxscan to locmaxscan (inclusive) at each location.  ",
    "           =0 Do not scan for dips. You must specify parameter locmaxsum    ",
    "              and you cannot specify scanstep, scanmax, scanmin, dipext.    ",
    "                                                                            ",
    " scanstep=0.5   Time Step to Scan (milliseconds per 1 lockey unit).         ",
    " scanmax =12.0  Maximum Time to Scan (milliseconds per 1 lockey unit).      ",
    " scanmin =-scanmax  Minimum Time to Scan (milliseconds per 1 lockey unit).  ",
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
    " dipuse=1  Estimate and use dips at input locations.                        ",
    "       =0  Estimate and use dips at output locations.                       ",
    "                                                                            ",
    " locmaxsum=locmaxscan  Maximum Relative lockey to include in sumopt options.",
    "                       Must be greater than 0. This program uses the input  ",
    "                       traces from -locmaxsum to locmaxsum (inclusive) at   ",
    "                       each output location.                                ",
    "                                                                            ",
    " sumopt=stack All options here only include input traces within -locmaxsum  ",
    "              to +locmaxsum range (inclusive) around each output location.  ",
    "              This default just adds them all together and divides by total.",
    "              Note that correct dips from the semblance scanning will cause ",
    "              these traces to be locally flat during stacking and therefore ",
    "              the output will not have as much smearing as traces stacked   ",
    "              without dip compensation.                                     ",
    "       =sum   Just like option stack but do not divide by total.            ",
    "              (This option is intended to help with parameter setting).     ",
    "       =one   Find the input trace with lockey value nearest the output     ",
    "              location, apply the dips, and output it. Effectively, this    ",
    "              is interpolation using just the nearest trace.                ",
    "       =two   Find the two input traces with lockey value nearest output,   ",
    "              apply the dips, sum both traces and divide by 2, and output.  ",
    "       =twolinear Find the two input traces with lockey values nearest      ",
    "                  output, apply the dips, stack them with linear weight     ",
    "                  based on their locations, and output. (Often one of the   ",
    "                  two traces will get all the weight, of course).           ",
    "       =dadstack  Detect-Amplitude-Differences. You should expect to have   ",
    "                  to do considerable testing of the other dad parameters    ",
    "                  for this option. Precise setting of the dad parameters    ",
    "                  may result in the output traces being sharper at rapid    ",
    "                  amplitude changes (such as locations where horizons start ",
    "                  or terminate. These dad options search for rapid changes  ",
    "                  in amplitude ACROSS LOCATIONS at each sample time and     ",
    "                  only stacks the sample values from -locmaxsum TO the      ",
    "                  detected location -or- FROM the detected location to the  ",
    "                  +locmaxsum location. Which side is stacked depends on     ",
    "                  WHERE detected location is relative to output location.   ",
    "             ***  If a rapid change in amplitude is not detected at this    ",
    "                  sample time then ALL samples are stacked.                 ",
    "       =dadsum  Just like option dadstack but do not divide by total.       ",
    "                (This option is intended to help with parameter setting).   ",
    "                                                                            ",
    " dadscale=1.0 Detect-Amplitude-Differences scaler.                          ",
    "             This parameter can only be set if sumopt=dadstack or dadsum.   ",
    "             When a rapid change in amplitude is detected this parameter    ",
    "             causes the range of locations stacked at this sample to        ",
    "             increase or decrease.                                          ",
    "       Note: Usual range for this value is 0.9 to 1.1                       ",
    "                                                                            ",
    " dadnoise=0.0  Detect-Amplitude-Differences Noise Absolute Amplitude.       ",
    "              This parameter can only be set if sumopt=dadstack or dadsum.  ",
    "              If the average absolute amplitude at the sample time across   ",
    "              the range of input traces is less than this value, this       ",
    "              sample time is not considered to have detected a rapid change ",
    "              in amplitude (so all locations are stacked for output sample).",
    "                                                                            ",
    " dadext=12.   Detect-Amplitude-Differences extender (milliseconds).         ",
    "              This parameter can only be set if sumopt=dadstack or dadsum.  ",
    "              At each output location the amplitude difference detector     ",
    "              estimates if-and-where a rapid difference is located in the   ",
    "              locmaxsum range of input traces. This parameter performs a    ",
    "              weighted-average of those estimated rapid difference locations",
    "              symmetrically in time from -dadext to +dadext. The weights    ",
    "              used are the sum of absolute amplitudes of the seismic values ",
    "              at the same samples in the locmaxsum range of input traces.   ",
    "              Effectively, this means the detected-amplitude-difference     ",
    "              locations at large horizons are trusted more, and thus they   ",
    "              are extended up and down in time near those larger horizons.  ",
    "                                                                            ",
    " outdip=   Do not output this file.                                         ",
    "       =(file name) Output the dip values chosen by the semblance scanning  ",
    "           into trace samples in this file. The trace header values are     ",
    "           copied from the seismic traces. If dipuse=1 then they are the    ",
    "           input headers, if dipuse=0 then they are the output headers.     ",
    "     Note: These dips are after dipext processing. But even with dipext=0   ",
    "           they will not have scanstep granularity because the scanned      ",
    "           maximum semblance dip undergoes a 3-point interpolation to get   ",
    "           a better approximation of the dip at the true peak in semblance. ",
    "     Note: It is technically more correct to call the dips made and used    ",
    "           by this program APPARANT dips. The dip values are put into the   ",
    "           trace samples and are in milliseconds per 1 location unit.       ",
    "           That is, if the lockey=cdp then the dips are ms per cdp. There   ",
    "           is no dip azimuth output. The dips are in the direction from one ",
    "           location to next (if lockey=cdp this is one cdp toward next cdp).",
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
    "           If dipuse=1 then this file must have the same order as the       ",
    "           order as the seismic traces you are inputting in this run.       ",
    "           If dipuse=0 then this file must have the same INPUT order as     ",
    "           the OUTPUT order of this run.                                    ",
    "     Note: The intention here is to allow you to output dips from one run   ",
    "           of this program and use them in a second run of this program.    ",
    "           The first run can input smoothed, denoised seismic to get dips,  ",
    "           and the second run can input the un-smoothed, noisy seismic and  ",
    "           use the dips from the first run. As well, you can also actually  ",
    "           smooth and denoise the DIP VALUES in this file before input here.",
    "                                                                            ",
    "                                                                            ",
    " bypass=0  Treat all output locations the same.                             ",
    "       =1  If an input trace exists exactly at the output location, COPY    ",
    "           the input trace without any change except to reset the lockey    ",
    "           to the output location value. Process all other output locations ",
    "           as indicated by the other options you have specified.            ",
    "       =-1 Only process and output traces where an input trace is exactly   ",
    "           on an output location. Typically you use this option when you    ",
    "           just want to smooth the input, without creating any more output  ",
    "           locations (so set outputalign=0 and usually set outputinc=1).    ",
    "  CAUTION: Options 1 and -1 only retain the first trace with duplicate      ",
    "           input location value. For instance if lockey=offset and          ",
    "           groupkey=cdp then only the first trace with the same offset      ",
    "           within cdps will be copied or processed to output.               ",
    "                                                                            ",
    " numrolling=2*locmaxscan+2*locmaxsum+5. This is maximum number of traces    ",
    "              to store internally. Traces are read-in as needed to keep the ",
    "              output locations centred in this storage (a rolling buffer).  ",
    "          *** For certain setups, this default is too large. For instance,  ",
    "              when lockey=offset and groupkey=cdp this value can be set to  ",
    "              the maximum cdp fold (but must always be set greater than 4). ",
    "          *** For certain setups, this default is too small. For instance,  ",
    "              when duplicate input lockey values exist.                     ",
    "        Note: This program expects input locations to have approximately    ",
    "              the same interval between them and to be distributed nicely.  ",
    "              The input locations must be in increasing-OR-EQUAL order.     ",
    "              So, traces can have duplicate lockey values provided this     ",
    "              parameter is set large enough. Whether these duplicates have  ",
    "              serious effects on the output traces varies with situation.   ",
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

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut trd = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut lockey = String::from("cdp");
    if countparval("lockey") > 0 {
        getparstring("lockey", &mut lockey);
    }
    let lockase = GetCase(&lockey);
    if lockase < 1 {
        err!(
            "**** Error: Specified lockey name %s is not recognized.",
            lockey
        );
    }
    let mut groupkey = String::new();
    let mut groupkase = -1;
    let mut usegroup = 0;
    let mut thisgroup = -1.0e30;
    let mut newgroup = 0;
    if countparval("groupkey") > 0 {
        getparstring("groupkey", &mut groupkey);
        groupkase = GetCase(&groupkey);
        if groupkase < 1 {
            err!(
                "**** Error: Specified groupkey name %s is not recognized.",
                groupkey
            );
        }
        usegroup = 1;
    }
    let mut locmaxscan = 3.0f64;
    if getpardouble("locmaxscan", &mut locmaxscan) == 0 {
        locmaxscan = 3.0;
    }
    if locmaxscan < 0.0 {
        err!("error: locmaxscan must be greater than or equal to 0.");
    }
    let locminscan = 0.0 - locmaxscan;
    if locmaxscan == 0.0 {
        if countparval("scanstep") > 0
            || countparval("scanmax") > 0
            || countparval("scanmin") > 0
            || countparval("dipext") > 0
            || countparval("locmaxsum") < 1
        {
            err!("error: locmaxscan=0 does not allow parameters scanstep, scanmax, scanmin, dipext, and REQUIRES locmaxsum");
        }
    }
    let mut locmaxsum = locmaxscan;
    if getpardouble("locmaxsum", &mut locmaxsum) == 0 {
        locmaxsum = locmaxscan;
    }
    if locmaxsum <= 0.0 {
        err!("error: locmaxsum must be greater than 0.");
    }
    let locminsum = 0.0 - locmaxsum;
    let mut dipuse = 1i32;
    if getparint("dipuse", &mut dipuse) == 0 {
        dipuse = 1;
    }
    if dipuse < 0 || dipuse > 1 {
        err!("error: dipuse parameter is out-of-range.");
    }
    let mut bypass = 0i32;
    if getparint("bypass", &mut bypass) == 0 {
        bypass = 0;
    }
    if bypass < -1 || bypass > 1 {
        err!("error: bypass parameter is out-of-range.");
    }
    let mut dipext = 12.0f32;
    if getparfloat("dipext", &mut dipext) == 0 {
        dipext = 12.0;
    }
    if dipext < 0.0 {
        err!("error: dipext must be equal to or greater than 0.0");
    }
    let mut sstep = 1.0f32;
    if getparfloat("scanstep", &mut sstep) == 0 {
        sstep = 1.0;
    }
    if sstep < 0.0 {
        err!("error: scanstep must be equal to or greater than 0.0");
    }
    let mut smax = 12.0f32;
    if getparfloat("scanmax", &mut smax) == 0 {
        smax = 12.0;
    }
    let mut smin = -8.0f32;
    if getparfloat("scanmin", &mut smin) == 0 {
        smin = (0.0 - smax as f64) as f32;
    }
    if smax < smin {
        err!("error: scanmax must be equal to or greater than scanmin");
    }
    let mut numrolling = 5i32;
    // Inclusive ranges and floating point tolerance account for the extra five.
    if getparint("numrolling", &mut numrolling) == 0 {
        numrolling = (locmaxscan - locminscan + locmaxsum - locminsum + 5.0) as i32;
    }
    if numrolling < 5 {
        err!("error: numrolling must be greater than 5");
    }
    let numrolling = numrolling as usize;
    let mut outputinc = 1.0f64;
    if getpardouble("outputinc", &mut outputinc) == 0 {
        outputinc = 1.0;
    }
    if outputinc <= 0.0 {
        err!("error: outputinc must be greater than 0.0");
    }
    let outputcni = if outputinc < 1.0 { outputinc } else { 1.0 };
    let mut outputalign = 0.0f64;
    if getpardouble("outputalign", &mut outputalign) == 0 {
        outputalign = 0.0;
    }
    if outputalign < 0.0 {
        err!("error: outputalign must be equal to or greater than 0.0");
    }
    outputalign %= outputinc;
    let mut outputshift = 0.0f64;
    if getpardouble("outputshift", &mut outputshift) == 0 {
        outputshift = 0.0;
    }
    let mut isumopt = 3;
    let mut ilinear = 0;
    let mut sumopt = String::new();
    if countparval("sumopt") > 0 {
        isumopt = 0;
        getparstring("sumopt", &mut sumopt);
        if sumopt == "one" {
            isumopt = 1;
        }
        if sumopt == "two" {
            isumopt = 2;
        }
        if sumopt == "twolinear" {
            isumopt = 2;
            ilinear = 1;
        }
        if sumopt == "stack" {
            isumopt = 3;
        }
        if sumopt == "sum" {
            isumopt = 4;
        }
        if sumopt == "dadstack" {
            isumopt = 5;
        }
        if sumopt == "dadsum" {
            isumopt = 6;
        }
        if isumopt < 1 {
            err!(
                "**** Error: Specified sumopt=%s option not recognized.",
                sumopt
            );
        }
    }
    if isumopt < 5
        && (countparval("dadscale") > 0 || countparval("dadext") > 0 || countparval("dadnoise") > 0)
    {
        err!("error: only sumopt=dadstack or dadsum allows specifying parameters dadscale, dadext, and dadnoise");
    }
    let mut dadscale = 1.0f32;
    if getparfloat("dadscale", &mut dadscale) == 0 {
        dadscale = 1.0;
    }
    let mut dadext = 12.0f32;
    if getparfloat("dadext", &mut dadext) == 0 {
        dadext = 12.0;
    }
    if dadext < 0.0 {
        err!("error: dadext must be equal to or greater than 0.0");
    }
    let mut dadnoise = 0.0f32;
    if getparfloat("dadnoise", &mut dadnoise) == 0 {
        dadnoise = 0.0;
    }
    if dadnoise < 0.0 {
        err!("error: dadnoise must be equal to or greater than 0.0");
    }
    let mut dname = String::new();
    let mut fpoutdip = None;
    let mut fpoutpow = None;
    let mut fpoutsem = None;
    let mut fpindip = None;
    if countparval("outdip") > 0 {
        getparstring("outdip", &mut dname);
        fpoutdip = Some(efopen(&dname, "w"));
    }
    if countparval("outpow") > 0 {
        getparstring("outpow", &mut dname);
        fpoutpow = Some(efopen(&dname, "w"));
    }
    if countparval("outsem") > 0 {
        getparstring("outsem", &mut dname);
        fpoutsem = Some(efopen(&dname, "w"));
    }
    // Do this last to retain dname in case of an error.
    if countparval("indip") > 0 {
        if countparval("scanstep") > 0
            || countparval("scanmax") > 0
            || countparval("scanmin") > 0
            || countparval("dipext") > 0
            || countparval("locmaxscan") > 0
        {
            err!("error: indip file does not allow parameters scanstep, scanmax, scanmin, dipext, and locmaxscan");
        }
        getparstring("indip", &mut dname);
        fpindip = Some(efopen(&dname, "r"));
    }
    let mut rollloc = vec![0i32; numrolling];
    let mut rolldun = vec![0i32; numrolling];
    let mut rolltrace = vec![Trace::new(); numrolling];
    let mut onetrace = Trace::new();
    let mut rolldip: Vec<Vec<f32>> = vec![Vec::new(); numrolling];
    let mut jsumdat: Vec<Vec<f32>> = vec![Vec::new(); numrolling];
    let mut locit = vec![0i32; numrolling];
    let mut arng = vec![0usize; numrolling];
    let mut ifori = vec![Vec::<i32>::new(); numrolling];
    let mut afori = vec![Vec::<f32>::new(); numrolling];
    checkpars();
    let mut nproct = 0i32;
    let mut mproct = 0i32;
    let mut icycle = 1;
    let mut iroll = 1usize;
    let mut iend = 0i32;
    let mut maxstoredloc = 0i32;
    let mut outputloc = 0.0f64;
    let mut locprev = 0.0f64;
    let mut locdups = 0i32;
    let mut dt = 0.0f32;
    let mut maxsd = 1i32;
    let mut minsd = -1i32;
    let mut maxdad = 1i32;
    let mut mindad = -1i32;
    let mut jsmax = 0.0f32;
    let mut lendip = 0usize;
    let mut dipwrk = Vec::<f32>::new();
    let mut dipwrq = Vec::<f32>::new();
    let mut onedip = Vec::<f32>::new();
    let mut onepow = Vec::<f32>::new();
    let mut onepw2 = Vec::<f32>::new();
    let mut onetmp = Vec::<f32>::new();
    let mut onetm2 = Vec::<f32>::new();
    let mut onesem = Vec::<f32>::new();
    let mut onesm2 = Vec::<f32>::new();
    let mut onetim = Vec::<f32>::new();
    let mut oneand = Vec::<f32>::new();
    let mut jsumn = 0usize;
    let mut jsumn2 = 0usize;
    while icycle == 1 {
        if gettr(&mut tr) != 0 {
            if nproct == 0 {
                iroll = 0;
                dt = (tr.dt as f64 / 1000.0) as f32;
                sstep = sstep / dt;
                smin = smin / dt;
                smax = smax / dt;
                maxsd = ((dipext as f64 + 0.01) / dt as f64) as i32;
                minsd = 0 - maxsd;
                maxdad = ((dadext as f64 + 0.01) / dt as f64) as i32;
                mindad = 0 - maxdad;
                jsmax = 1.0f32 + (smax - smin) / sstep;
                dipwrk = vec![0.0; jsmax as usize];
                dipwrq = vec![0.0; jsmax as usize];
                lendip = tr.ns as usize;
                onedip = vec![0.0; lendip];
                onepow = vec![0.0; lendip];
                onepw2 = vec![0.0; lendip];
                onetmp = vec![0.0; lendip];
                onetm2 = vec![0.0; lendip];
                onesem = vec![0.0; lendip];
                onesm2 = vec![0.0; lendip];
                onetim = vec![0.0; lendip];
                oneand = vec![0.0; lendip];
                for i in 0..numrolling {
                    rolldip[i] = vec![0.0; lendip];
                    jsumdat[i] = vec![0.0; lendip];
                    ifori[i] = vec![0; jsmax as usize];
                    afori[i] = vec![0.0; jsmax as usize];
                    rollloc[i] = (fromhead(&tr, lockase) + 3.0 * (locminscan + locminsum)
                        - i as f64
                        - 5.0) as i32;
                }
                maxstoredloc = rollloc[0];
                outputloc =
                    fromhead(&tr, lockase) * outputinc + outputalign + locminscan + locminsum - 1.0;
                locprev = fromhead(&tr, lockase) - 1.0;
            }
            nproct += 1;
            if usegroup > 0 {
                newgroup = 0;
                if thisgroup != fromhead(&tr, groupkase) {
                    thisgroup = fromhead(&tr, groupkase);
                    locprev = fromhead(&tr, lockase) - 1.0;
                    newgroup = 1;
                    iend = (2.0 * (locmaxscan + locmaxsum) + 2.0) as i32;
                }
            }
            if fromhead(&tr, lockase) <= locprev {
                if fromhead(&tr, lockase) == locprev {
                    locdups += 1;
                } else {
                    err!(
                        "input lockey %g is less than previous trace %g",
                        fromhead(&tr, lockase),
                        locprev
                    );
                }
            }
            locprev = fromhead(&tr, lockase);
        } else {
            icycle = 0;
            iend = (2.0 * (locmaxscan + locmaxsum) + 2.0) as i32;
        }

        while outputloc < maxstoredloc as f64 - locmaxscan - locmaxsum + iend as f64 {
            if dipuse < 1 {
                if fpindip.is_none() {
                    scanit(
                        &rolltrace,
                        &rollloc,
                        outputloc,
                        locminscan,
                        locmaxscan,
                        smin,
                        sstep,
                        jsmax as usize,
                        &mut onedip,
                        &mut onepw2,
                        &mut onesm2,
                        &mut dipwrk,
                        &mut dipwrq,
                        &mut arng,
                        &mut ifori,
                        &mut afori,
                        lendip,
                    );
                    wsmooth(&onedip, &onepw2, lendip, minsd, maxsd, &mut onetim);
                    for i in 0..lendip {
                        onedip[i] = onetim[i];
                    }
                } else {
                    let mut ksum = 0;
                    for j in 0..numrolling {
                        if rollloc[j] as f64 >= outputloc + locminsum
                            && rollloc[j] as f64 <= outputloc + locmaxsum
                        {
                            ksum = 1;
                            break;
                        }
                    }
                    if ksum > 0 {
                        if fgettr(fpindip.as_ref().unwrap(), &mut trd) == 0 {
                            err!("next output location needs a trace from dipfile: %s   but get failed (near location=%g)", dname, outputloc);
                        }
                        if ((outputloc - outputalign) / outputcni + outputshift
                            - fromhead(&trd, lockase))
                        .abs()
                            > 1.0e-9
                        {
                            err!("lockey value for dipfile input trace and output location are not the same %g %g", (outputloc-outputalign)/outputcni+outputshift, fromhead(&trd,lockase));
                        }
                        for i in 0..lendip {
                            onedip[i] = trd.data[i] / dt;
                        }
                    }
                }
            }
            let mut jsum = 0usize;
            let mut jnear = 0usize;
            let mut diffj = 1.0e30;
            let mut diffj2 = 1.0e30;
            let mut centloc = 0.0f64;
            let mut klow = 0usize;
            let mut diffk = 1.0e30;
            for j in 0..numrolling {
                if rollloc[j] as f64 >= outputloc + locminsum
                    && rollloc[j] as f64 <= outputloc + locmaxsum
                {
                    if dipuse > 0 {
                        if rolldun[j] < 1 {
                            scanit(
                                &rolltrace,
                                &rollloc,
                                rollloc[j] as f64,
                                locminscan,
                                locmaxscan,
                                smin,
                                sstep,
                                jsmax as usize,
                                &mut onedip,
                                &mut onepow,
                                &mut onesem,
                                &mut dipwrk,
                                &mut dipwrq,
                                &mut arng,
                                &mut ifori,
                                &mut afori,
                                lendip,
                            );
                            wsmooth(&onedip, &onepow, lendip, minsd, maxsd, &mut rolldip[j]);
                            rolldun[j] = 1;
                            if fpoutdip.is_some() || fpoutpow.is_some() || fpoutsem.is_some() {
                                onetrace.clone_from(&rolltrace[j]);
                                if let Some(fp) = &fpoutdip {
                                    for i in 0..lendip {
                                        onetrace.data[i] = rolldip[j][i] * dt;
                                    }
                                    fputtr(fp, &mut onetrace);
                                }
                                if let Some(fp) = &fpoutpow {
                                    for i in 0..lendip {
                                        onetrace.data[i] = onepow[i];
                                    }
                                    fputtr(fp, &mut onetrace);
                                }
                                if let Some(fp) = &fpoutsem {
                                    for i in 0..lendip {
                                        onetrace.data[i] = onesem[i];
                                    }
                                    fputtr(fp, &mut onetrace);
                                }
                            }
                        }
                        for i in 0..lendip {
                            onetim[i] = (i as f32 as f64
                                + rolldip[j][i] as f64 * (rollloc[j] as f64 - outputloc))
                                as f32;
                        }
                    } else {
                        for i in 0..lendip {
                            onetim[i] = (i as f32 as f64
                                + onedip[i] as f64 * (rollloc[j] as f64 - outputloc))
                                as f32;
                        }
                    }
                    // Shift by the time computed from dips and location difference.
                    ints8r(
                        1.0,
                        0.0,
                        &rolltrace[j].data[..lendip],
                        0.0,
                        0.0,
                        &onetim,
                        &mut jsumdat[jsum],
                    );
                    if (rollloc[j] as f64 - outputloc).abs() < diffj {
                        diffj2 = diffj;
                        jsumn2 = jsumn;
                        jnear = j;
                        jsumn = jsum;
                        diffj = (rollloc[j] as f64 - outputloc).abs();
                    } else if (rollloc[j] as f64 - outputloc).abs() < diffj2 {
                        diffj2 = (rollloc[j] as f64 - outputloc).abs();
                        jsumn2 = jsum;
                    }
                    if (rollloc[j] as f64) < diffk {
                        diffk = rollloc[j] as f64;
                        klow = jsum;
                    }
                    centloc += rollloc[j] as f64;
                    locit[jsum] = rollloc[j];
                    jsum += 1;
                }
            }
            if jsum > 0 {
                onetrace.clone_from(&rolltrace[jnear]);
                tohead(
                    &mut onetrace,
                    lockase,
                    (outputloc - outputalign) / outputcni + outputshift,
                );
                let mut iput = 1;
                if bypass > 0 && diffj < 1.0e-9 {
                    // Keep the input samples copied above.
                } else if bypass < 0 && diffj > 1.0e-9 {
                    iput = 0;
                } else if isumopt == 1 || jsum == 1 {
                    for i in 0..lendip {
                        onetrace.data[i] = jsumdat[jsumn][i];
                    }
                } else if isumopt == 2 || jsum == 2 {
                    if ilinear == 0
                        || (outputloc - locit[jsumn] as f64).abs()
                            == (outputloc - locit[jsumn2] as f64).abs()
                    {
                        for i in 0..lendip {
                            onetrace.data[i] =
                                (0.5 * (jsumdat[jsumn][i] + jsumdat[jsumn2][i]) as f64) as f32;
                        }
                    } else if outputloc == locit[jsumn] as f64 {
                        for i in 0..lendip {
                            onetrace.data[i] = jsumdat[jsumn][i];
                        }
                    } else if outputloc == locit[jsumn2] as f64 {
                        for i in 0..lendip {
                            onetrace.data[i] = jsumdat[jsumn2][i];
                        }
                    } else {
                        let ratio = ((outputloc - locit[jsumn] as f64).abs()
                            / ((outputloc - locit[jsumn] as f64).abs()
                                + (outputloc - locit[jsumn2] as f64).abs()))
                            as f32;
                        for i in 0..lendip {
                            onetrace.data[i] = ((1.0 - ratio as f64) * jsumdat[jsumn][i] as f64
                                + (ratio * jsumdat[jsumn2][i]) as f64)
                                as f32;
                        }
                    }
                } else if isumopt == 3 {
                    for i in 0..lendip {
                        onetrace.data[i] = 0.0;
                        for j in 0..jsum {
                            onetrace.data[i] += jsumdat[j][i];
                        }
                        onetrace.data[i] /= jsum as f32;
                    }
                } else if isumopt == 4 {
                    for i in 0..lendip {
                        onetrace.data[i] = 0.0;
                        for j in 0..jsum {
                            onetrace.data[i] += jsumdat[j][i];
                        }
                    }
                } else {
                    // Detect rapid amplitude changes across locations. Absolute-amplitude
                    // weighting extends their estimated locations around strong horizons.
                    centloc /= jsum as f64;
                    for i in 0..lendip {
                        oneand[i] = 0.0;
                        onetmp[i] = 1.0e-30;
                        for j in 0..jsum {
                            let k = (j + klow) % jsum;
                            oneand[i] = (oneand[i] as f64
                                + (jsumdat[k][i] as f64).abs()
                                    * (locit[k] as f64 - centloc) as f32 as f64)
                                as f32;
                            onetmp[i] = (onetmp[i] as f64 + (jsumdat[k][i] as f64).abs()) as f32;
                        }
                        if onetmp[i] / jsum as f32 > dadnoise {
                            oneand[i] = oneand[i] / onetmp[i];
                        } else {
                            oneand[i] = 0.0;
                            onetmp[i] = 1.0e-30;
                        }
                    }
                    wsmooth(&oneand, &onetmp, lendip, mindad, maxdad, &mut onetm2);
                    for i in 0..lendip {
                        onetm2[i] = (2.0 * onetm2[i] as f64 * dadscale as f64
                            + (locmaxsum - locminsum)) as f32;
                        if onetm2[i] as f64 > locmaxsum - locminsum + 0.1 {
                            onetm2[i] = (onetm2[i] as f64 - (locmaxsum - locminsum) - 1.0) as f32;
                        }
                        let (mut jfirst, mut jlast);
                        if onetm2[i] as f64 > (locmaxsum - locminsum) / 2.0 - 0.5 {
                            jfirst = 0i32;
                            jlast = (1.5 + onetm2[i] as f64) as i32;
                            if jlast > jsum as i32 {
                                jlast = jsum as i32;
                            }
                        } else {
                            jfirst = (1.5 + onetm2[i] as f64) as i32;
                            jlast = jsum as i32;
                            if jfirst < 0 {
                                jfirst = 0;
                            }
                        }
                        onetrace.data[i] = 0.0;
                        for j in jfirst..jlast {
                            let k = (j as usize + klow) % jsum;
                            onetrace.data[i] += jsumdat[k][i];
                        }
                        if isumopt == 5 && jlast - jfirst > 1 {
                            onetrace.data[i] /= (jlast - jfirst) as f32;
                        }
                    }
                }
                if iput > 0 {
                    puttr(&mut onetrace);
                    mproct += 1;
                    if dipuse < 1 {
                        if let Some(fp) = &fpoutdip {
                            for i in 0..lendip {
                                onetrace.data[i] = onedip[i] * dt;
                            }
                            fputtr(fp, &mut onetrace);
                        }
                        if let Some(fp) = &fpoutpow {
                            for i in 0..lendip {
                                onetrace.data[i] = onepw2[i];
                            }
                            fputtr(fp, &mut onetrace);
                        }
                        if let Some(fp) = &fpoutsem {
                            for i in 0..lendip {
                                onetrace.data[i] = onesm2[i];
                            }
                            fputtr(fp, &mut onetrace);
                        }
                    }
                }
            }
            outputloc += outputinc;
        }
        if newgroup > 0 {
            for i in 0..numrolling {
                rollloc[i] = (fromhead(&tr, lockase) + 3.0 * (locminscan + locminsum)
                    - i as f64
                    - 5.0) as i32;
            }
            outputloc =
                fromhead(&tr, lockase) * outputinc + outputalign + locminscan + locminsum - 1.0;
            iroll = 0;
            iend = 0;
        }
        iroll += 1;
        if iroll == numrolling {
            iroll = 0;
        }
        rollloc[iroll] = fromhead(&tr, lockase) as i32;
        rolldun[iroll] = 0;
        rolltrace[iroll].clone_from(&tr);
        maxstoredloc = fromhead(&tr, lockase) as i32;
        if fpindip.is_some() && dipuse > 0 && icycle > 0 {
            if fgettr(fpindip.as_ref().unwrap(), &mut trd) == 0 {
                err!("next input seismic trace needs a trace from indip file: %s   but get failed (near location=%d)", dname, rollloc[iroll]);
            }
            if rollloc[iroll] as f64 != fromhead(&trd, lockase) {
                err!("lockey value for dipfile trace %g   and input seismic trace %d are not the same ", fromhead(&trd,lockase),rollloc[iroll]);
            }
            for i in 0..lendip {
                rolldip[iroll][i] = trd.data[i];
            }
        }
    }
    if let Some(fp) = &fpindip {
        efclose(fp);
    }
    if let Some(fp) = &fpoutdip {
        efclose(fp);
    }
    warn!(
        "Number of traces input=%d.  Number output=%d.",
        nproct, mproct
    );
    if locdups > 0 {
        if usegroup < 1 {
            warn!(
                "There were %d input traces with duplicate lockey values.",
                locdups
            );
        } else {
            warn!(
                "There were %d input traces with duplicate lockey values (within their groups).",
                locdups
            );
        }
    }
    exit(CWP_Exit());
}

// Semblance scanner: sum squared amplitudes divided by power and trace count.
// Interpolate the best scan with its two neighbours to reduce scan granularity.
#[allow(clippy::too_many_arguments)]
fn scanit(
    atrace: &[Trace],
    aloc: &[i32],
    scanloc: f64,
    locminscan: f64,
    locmaxscan: f64,
    smin: f32,
    sstep: f32,
    jsmax: usize,
    adip: &mut [f32],
    apow: &mut [f32],
    asem: &mut [f32],
    awrk: &mut [f32],
    awrq: &mut [f32],
    arng: &mut [usize],
    ifori: &mut [Vec<i32>],
    afori: &mut [Vec<f32>],
    lendip: usize,
) {
    for n in 0..lendip {
        adip[n] = 0.0;
        apow[n] = 0.0;
        asem[n] = 0.0;
    }
    let mut numinrange = 0usize;
    for i in 0..aloc.len() {
        if aloc[i] as f64 >= scanloc + locminscan && aloc[i] as f64 <= scanloc + locmaxscan {
            arng[numinrange] = i;
            numinrange += 1;
        }
    }
    if numinrange < 1 {
        return;
    }
    if numinrange < 2 {
        for n in 0..lendip {
            apow[n] = atrace[arng[0]].data[n] * atrace[arng[0]].data[n];
            asem[n] = 1.0;
        }
        return;
    }
    // Precompute whole and fractional samples for linear interpolation.
    for j in 0..jsmax {
        let sval = smin + j as f32 * sstep;
        for i in 0..numinrange {
            let mut m = (sval as f64 * (aloc[arng[i]] as f64 - scanloc)) as i32;
            if sval as f64 * (aloc[arng[i]] as f64 - scanloc) > 0.0 {
                m += 1;
            }
            let alin = (m as f64 - sval as f64 * (aloc[arng[i]] as f64 - scanloc)) as f32;
            ifori[i][j] = m;
            afori[i][j] = alin;
        }
    }
    for n in 0..lendip {
        for j in 0..jsmax {
            awrk[j] = 0.0;
            awrq[j] = 1.0e-30; // Prevent divide by zero in mute zones.
            for i in 0..numinrange {
                let m = ifori[i][j];
                let nm = n as i32 + m;
                if nm > 0 && nm < lendip as i32 {
                    let alin = afori[i][j];
                    let aval = ((alin * atrace[arng[i]].data[nm as usize - 1]) as f64
                        + (1.0 - alin as f64) * atrace[arng[i]].data[nm as usize] as f64)
                        as f32;
                    awrk[j] += aval;
                    awrq[j] += aval * aval;
                }
            }
        }
        let mut thesem = 0.0f32;
        let mut thej = jsmax / 2;
        for j in 0..jsmax {
            if awrk[j] * awrk[j] / awrq[j] > thesem {
                thesem = awrk[j] * awrk[j] / awrq[j];
                thej = j;
            }
        }
        if (thesem as f64) < 1.0e-30 {
            adip[n] = 0.0;
            asem[n] = 0.0;
            apow[n] = 0.0;
        } else if thej > 0 && thej < jsmax - 1 {
            let sl = awrk[thej - 1] * awrk[thej - 1] / awrq[thej - 1];
            let sr = awrk[thej + 1] * awrk[thej + 1] / awrq[thej + 1];
            let sp = (0.5 * (sl - sr) as f64
                / (sl as f64 - 2.0 * thesem as f64 + sr as f64 + 1.0e-30))
                as f32;
            adip[n] = smin + thej as f32 * sstep + sp * sstep;
            asem[n] =
                ((thesem as f64 - 0.25 * (sl - sr) as f64 * sp as f64) / numinrange as f64) as f32;
            apow[n] = awrq[thej] / numinrange as f32;
        } else {
            adip[n] = smin + thej as f32 * sstep;
            asem[n] = thesem / numinrange as f32;
            apow[n] = awrq[thej] / numinrange as f32;
        }
    }
}

// Weighted smoothing, with constant extension of the computed endpoints.
fn wsmooth(val: &[f32], wgt: &[f32], len: usize, minsd: i32, maxsd: i32, sval: &mut [f32]) {
    if len as i32 <= maxsd || len as i32 <= -minsd {
        // With very short traces C's endpoint loops access beyond its malloc
        // block, including negative indices. The original binary returns zero
        // dips for nt=1,2,3 with its default extent (verified with outdip).
        // Supply zero guard storage and execute those same assignments rather
        // than clamping the extent or changing the smoothing window.
        let pad = maxsd.max(-minsd) as usize + 1;
        let mut guarded = vec![0.0f32; len + 2 * pad];
        guarded[pad..pad + len].copy_from_slice(&sval[..len]);
        for i in 0..0 - minsd {
            guarded[(pad as i32 + i) as usize] = guarded[(pad as i32 - minsd) as usize];
        }
        for i in len as i32 - maxsd..len as i32 {
            guarded[(pad as i32 + i) as usize] =
                guarded[(pad as i32 + len as i32 - maxsd - 1) as usize];
        }
        sval[..len].copy_from_slice(&guarded[pad..pad + len]);
        return;
    }
    for i in 0 - minsd..len as i32 - maxsd {
        let mut twgt = 1.0e-30f32;
        let mut tval = 0.0f32;
        for k in minsd..=maxsd {
            twgt += wgt[(i + k) as usize];
            tval += wgt[(i + k) as usize] * val[(i + k) as usize];
        }
        sval[i as usize] = tval / twgt;
    }
    for i in 0..0 - minsd {
        sval[i as usize] = sval[(0 - minsd) as usize];
    }
    for i in len as i32 - maxsd..len as i32 {
        sval[i as usize] = sval[(len as i32 - maxsd - 1) as usize];
    }
}
