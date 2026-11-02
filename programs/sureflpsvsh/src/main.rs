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

//! SUREFLPSVSH - reflectivity modeling of PSV or SH waves for layered earth.
//! Port of Refl/main/sureflpsvsh.c.
//! Copyright (c) Colorado School of Mines, 2011. All rights reserved.
//! Credits: Original Fortran 77 PSV version written by Subshashis Mallick (1988).
//! Original Fortran 77 SH version written by Mrinal Sen (1988), based on the
//! PSV version by Subshashis Mallick. Translated to C, expanded and reformatted
//! for SU by Gabriel Alvarez (1995).

use oxwave_core::prelude::*;
use oxwave_lib_refl::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUREFLPSVSH - REFlectivity modeling of either PSV or SH waves for layered",
    "\t\tearth model\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " surelfpsvsh  required parameters [optional parameters]\t \t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required Parameters:\t\t\t\t\t\t\t",
    " m0=\t\t\tSeismic moment\t\t\t\t\t",
    " p2w=\t\t\tmaximum ray parameter value to which the\t",
    "\t\t\tsynthetics are computed\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional Parameters:\t\t\t\t\t\t\t",
    " Main Flags:\t\t\t\t\t\t\t\t",
    " wtype=1\t\t=1 for PSV, =2 for SH\t\t\t\t",
    " stype=1\t\t=1 if the moemnt tensor components are given\t",
    "\t\t\t=2 if they are to be computed from fault\t",
    "\t\t\t plane mechanism parameters\t\t\t",
    " wfield=2\t\t=1 for displacement, =2 for particle velocity,\t",
    "\t\t\t=3 for particle acceleration\t\t\t",
    " flt=0\t\t\t=1 to apply earth flattening correction\t\t",
    "\t\t\t(important for earthquake data)\t\t\t",
    " vsp=0\t\t\t=0 for surface data, =1 for VSP data\t\t",
    " int_type=1\t\t=1 to compute the slowness integration using\t",
    "\t\t\tthe trapezoidal rule. =2 to use a first order\t",
    "\t\t\tFilon scheme (faster but maybe noisier)\t\t",
    " verbose=0\t\t=0 no processing information is output (makes\t",
    "\t\t\tthe program run a litle faster), =1 to output\t",
    "\t\t\tprocessing information to the screen, =2 to\t",
    "\t\t\toutput processing information to a user\t\t",
    "\t\t\tsupplied file, =3 output processing\t\t",
    "\t\t\tinformation to both, the screen and a file\t",
    "\t\t\t\t\t\t\t\t\t",
    " Flags for special non-standard options (used only in rare ocasions)\t",
    " rand=0\t\t=1 to include random velocity and q layers\t",
    " qopt=0\t\t=1 if a q-correction is desired\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Flags for output data\t\t\t\t\t\t\t",
    " win=0\t\t\t=1 to apply a frequency domain Hanning\t\t",
    "\t\t\twindow prior to inverse FFT to time domain\t",
    " wavelet_type\t\t=1 for a spike, =2 for a Ricker wavelet\t\t",
    "\t\t\t=3 for an akb wavelet\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Main input parameters\t\t\t\t\t\t\t",
    " tsec=2.048\t\tlength of computed traces (in seconds)\t\t",
    " dt=0.004\t\ttime sampling interval (in seconds)\t\t",
    " nt=tsec/dt\t\tnumber of samples per trace\t\t\t",
    " nx=60\t\t\tnumber of traces (ranges) per shot\t\t",
    " nw=100\t\tnumber of frequencies to process\t\t",
    " nor=1\t\t\tnumber of depths at which receivers are\t\t",
    "\t\t\tlocated (different from one for VSP's)\t\t",
    " nlayers=10\t\tnumber of horizontal layers in the model\t",
    " fref=1.0\t\treference frequency (Hz)\t\t\t",
    " bx=0.0 \t\tfirst range (first trace offset) (km)\t\t",
    " dx=0.05\t\trange increment (trace to trace distance) km\t",
    " pw1=0.0 pw2=0.1\tapply Hanning window tapering to lower end\t",
    "\t\t\tof ray parameter computations (if set to zero\t",
    "\t\t\tgood default values are used)\t\t\t",
    " pw3=6.7 pw4=7.0\tapply Hanning window tapering to higher end\t",
    "\t\t\tof ray parameter computations (if set to zero\t",
    "\t\t\tgood default values are used)\t\t\t",
    " fs=0.07\t\tFilon sampling parameter, usually between\t",
    "\t\t\t0.07 and 0.12. Sampling is finer as fs\t\t",
    "\t\t\tincreases. For short range synthetics\t\t",
    "\t\t\t(<100km) use 0.07, for medium range (<1000km)\t",
    "\t\t\tuse 0.09 and for large\tranges (>1000km) use\t",
    "\t\t\t0.12 (this parameter is ignored if int_type=1)\t",
    " np=1300\t\tnumber of ray parameters used to compute the\t",
    "\t\t\tseismograms. Automatically set if int_type=2\t",
    " bp=0.0\t\tslowest ray parameter used to compute the\t",
    "\t\t\tseismograms. Set to zero if int_type=2\t\t",
    " decay=50.0\t\tdecay factor use to avoid time series\t\t",
    "\t\t\twraparound. A value 'n' for decay means that\t",
    "\t\t\tthe wrapped around signals are diminished by\t",
    "\t\t\ta factor of 'n'. The default of 50 is the\t",
    "\t\t\trecomended value\t\t\t\t",
    " lobs=\t\t\tarray[nor] of layers on top of which the\t",
    "\t\t\treceivers are located. If any receiver is\t",
    "\t\t\twithin a layer, introduce a ficticious layer\t",
    " cl=\t\t\tarray[nlayers] of compressional velocities\t",
    "\t\t\tfor each layer in km/s\t\t\t\t",
    " ct=\t\t\tarray[nlayers] of shear velocities for each\t",
    "\t\t\tlayer in km/s. A value of -1 assumes\t\t",
    "\t\t\tct=cl/sqrt(3.\t\t\t\t\t",
    " ql=\t\t\tarray[nlayers] of compressional q values\t",
    "\t\t\tfor each layer\t\t\t\t\t",
    " qt=\t\t\tarray[nlayers] of shear q values for each\t",
    "\t\t\tlayer. A value of -1 assumes qt=ql/2.\t\t",
    " rho=\t\t\tarray[nlayers] of densities for each layer\t",
    "\t\t\tin gr/cc\t\t\t\t\t",
    " t=\t\t\tabsolute thickness of the layer in kms\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Note: any of these arrays can also be input in a file as lobsfile=,clfile=,",
    " ctfile=, etc.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Source parameters:\t\t\t\t\t\t\t",
    " lsource=1\t\tlayer on top of which the source is located\t",
    " h1=1.0 h2=0.0\t\tvertical and horizontal linear part of the\t",
    "\t\t\tsource. For unit point source set h1=0, h2=0\t",
    "  \t\t\tand all other source parameters to zero. For\t",
    "\t\t\tearthquake source both h1 and h2 are usually\t",
    "\t\t\tzero\t\t\t\t\t \t",
    " m1,m2,m3\t\t[1][1],[1][2] and [2][2] components of the\t",
    "\t\t\tmoment tensor ([1][2]=[2][1] component)\t\t",
    " delta=0.0\t\tdip in degrees (necessary only if stype=1)\t",
    " lambda=0.0\t\trake in degrees (necessary only if stype=1)\t",
    " phis=0.0\t\tfault plane azimuth in degrees (necessary\t",
    "\t\t\tonly if stype=1)\t\t\t\t",
    " phi=0.0\t\tazimuth of the receiver location in degrees\t",
    "\t\t\t(necessary only if stype=1)\t\t  \t",
    " Parameters for output data:\t\t\t\t\t\t",
    " tlag=0.0\t\ttime lag to appy to the  seismograms (sec)\t",
    " nf=nw\t\t\tnumber of frequencies in output data\t\t",
    " fpeak=25.0\t\tpeak frequency for Ricker or akb wavelet (hz)\t",
    " red_vel=0.0\t\treducing velocity (km/s). If set to zero, the\t",
    "\t\t\tmaximum compressional velocity is used\t\t",
    "\t\t\t(this parameter is garanteed not to be <1.5)\t",
    " w1=0.0\t\tA Hanning window will be applied to freqs\t",
    "\t\t\tless than w1.If set to zero, the default of\t",
    "\t\t\t0.15 of maximum frequency is used (provided\t",
    "\t\t\tthat win=1)\t\t\t\t\t",
    " w2=0.0\t\tA Hanning window will be applied to freqs\t",
    "\t\t\tgreater than w2. If set to zero, the default\t",
    " \t\t\t0.85 of maximum frequency is used (provided\t",
    "\t\t\tthat win=1)\t\t\t\t\t",
    " nfilters=0\t\tnumber of filters to apply to the synthetics\t",
    " filters_phase=\tarray[nfilters] of: 0 for zero phase filters\t",
    " \t\t\tor 1 for minimum phase filters(can also be\t",
    "\t\t\tinput in a file via fphfile=\t\t\t",
    " filters_type=\t\tarray[nfilters] of: 1 for high cut filters,\t",
    "\t\t\tinput in a file via filphasefile=\t\t",
    " \t\t\t2 for low cut filters, 3 for notch filters\t",
    "\t\t\t(can also be input in a file via filtypefile=)\t",
    " dbpo=\t\t\tarray[nfilters] of: filter slopes in db/oct\t",
    "\t\t\t(can also be input in a file via dbpofile=)\t",
    " f1=\t\t\tarray[nfilters] of: frequency to start filter\t",
    "\t\t\taction (Hz)(can also be input in a file via\t",
    "\t\t\tf1file=)\t\t\t\t\t",
    " f2=\t\t\tarray[nfilters] of: frequency to end filter\t",
    "\t\t\taction (Hz). Only for notch filters\t\t",
    "\t\t\t(can also be input in a file via f2file=)\t",
    " wfp=\t\t\tname of output pressure seismogram file\t\t",
    "\t\t\t(only if wtype=1)\t\t\t\t",
    " wfr=\t\t\tname of output radial comp seismogram file\t",
    "\t\t\t(only if wtype=1)\t\t\t\t",
    " wfz=\t\t\tname of output vertical comp seismogram file\t",
    "\t\t\t(only if wtype=1)\t\t\t\t",
    " wft=\t\t\tname of output tangential comp seismogram file\t",
    "\t\t\t(only if wtype=1)\t\t\t  \t",
    " outf=info\t\tname of output processing information file\t",
    "\t\t\t(only if verbos=2 or 3)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Interpolation parameters (required only if layers with gradients desired)",
    " nlint=0\t\tnumber of times layer interp is required\t",
    " nintlayers=\t\tarray[nlint] of number of layers to interpol\t",
    "\t\t\teach time (can be input as file nintlayfile=)\t",
    " intlayers=\t\tarray[nlint] of layers on top of which to\t",
    "\t\t\tstart each interpolation\t\t\t",
    "\t\t\t(can be input as file intlayfile=)\t\t",
    " intlayth=\t\tarray[nlint] of layer thicknesses to interp\t",
    "\t\t\t(can be input as file intlaythfile=)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Other parameters (required only under very special circumstances)\t",
    " nrand_layers=0\tmaximum number of random layers allowed\t\t",
    "\t\t\t(only if rand=1)\t\t\t\t",
    " layer=0\t\tlayer on top of which the random velocity\t",
    "\t\t\tlayers are inserted (only if rand=1)\t\t",
    " zlayer=0.0\t\tthickness of random layers (only if rand=1)\t",
    "   \t\t\t(if zlayer<t[il], then zlayer=t[il])\t\t",
    " sdcl=0.0\t\tstandard deviation for compressional vels\t",
    "\t\t\t(only if rand=1)\t\t\t\t",
    " sdct=0.0\t\tstandard deviation for shear velocities\t\t",
    "\t\t\t(only if rand=1)\t\t\t\t",
    " layern=0\t\tlayer on top of which the q-option is invoked\t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " wrefp=1.0\t\treference frequency for compressional vels\t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " wrefs=1.0\t\treference frequency for shear velocities\t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " epsp=0.001\t\treference amplitude for comporessional vels\t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " epss=0.001\t\treference amplitude for shear velocities\t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " sigp=0.1\t\txxxxxx for comporessional vels\t\t\t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " sigs=0.1\t\txxxxxx for shear velocities\t\t  \t",
    "\t\t\t(only if qopt=1)\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Gradient zones between two layers can be handled with the use of the\t",
    " layer interpolation option. The program will automatically compute and",
    " insert the required number of layers with the appropriate thicknesses\t",
    " between the layers above and below. There is no restriction as to how\t",
    " many layers the program can handle, provided enough computer power is\t",
    " available.\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of frequencies to be processed is the most critical parameter",
    " to determine how long the program will take to run. The maximum frequency",
    " that can be present in the seismogram is nw/tsec Hz, however, a lower\t",
    " frequency can be selected for the output via the nf param.\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The number of computed time samples is tsec/dt, however, the user may\t",
    " choose a smaller number of samples for the output traces by setting\t",
    " nt this can be useful when a short seismogram is required with a\t",
    " broadband frequency range (directly computing a small seismogram can\t",
    " be hazardous)\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The decay parameter should be chosen with care. A value of 50 seems\t",
    " to give good results, but depending on the data this parameter can\t",
    " boost up late wrapped around noncausal energy.\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The integration flag is important for short range data, in particular ",
    " for oil exploration, for which the first order Filon scheme, though\t",
    " faster, can produce noisier seismograms, specialy in the late parts of",
    " the record. A standard trapezoidal rule seems to work better but is\t",
    " slower by about 25% in a normal situation.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " When arrays are required as input, they can also be input as files,\t",
    " however, the program will only check that the number of parameters is ",
    " the same if the arrays are used or if the number of elements in the\t",
    " files are set. A combination of arrays and files is also permited.\t",
    "\t\t\t\t\t\t\t\t\t",
    " Examples of Source Parameters:\t\t\t\t\t",
    " For a vertical point force, set stype=1, h1=1, h2=0 and all moment\t",
    " tensor components to zero. Ignore all other source parameters.\t",
    " For an explosion, set stype=1, h1=h2=0, m1=m3=A, m2=0 where A is some\t",
    " constant (normally one). If A is negative, an implosion is generated\t",
    " instead.  Ignore all other source parameters.\t\t\t\t",
    " For a fault slip, set stype=2, set phis,phi,lambda,alpha and m0 to their",
    " corresponding values and ignore the moment tensor components\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(0, SDOC); // No input data.
    let mut m0 = 0.0f32;
    if getparfloat("m0", &mut m0) == 0 {
        err!("error: the seismic moment, m0, is a required parameter\n");
    }
    let mut rand: i32 = 0;
    if getparint("rand", &mut rand) == 0 {
        rand = 0;
    }
    let mut qopt: i32 = 0;
    if getparint("qopt", &mut qopt) == 0 {
        qopt = 0;
    }
    let _ = qopt; // The C parses this flag but never uses it.
    let mut stype: i32 = 1;
    if getparint("stype", &mut stype) == 0 {
        stype = 1;
    }
    let mut wtype: i32 = 1;
    if getparint("wtype", &mut wtype) == 0 {
        wtype = 1;
    }
    let mut wfield: i32 = 1;
    if getparint("wfield", &mut wfield) == 0 {
        wfield = 1;
    }
    let mut int_type: i32 = 1;
    if getparint("int_type", &mut int_type) == 0 {
        int_type = 1;
    }
    let mut flt: i32 = 0;
    if getparint("flt", &mut flt) == 0 {
        flt = 0;
    }
    let mut vsp: i32 = 0;
    if getparint("vsp", &mut vsp) == 0 {
        vsp = 0;
    }
    let mut win: i32 = 0;
    if getparint("win", &mut win) == 0 {
        win = 0;
    }
    let mut wavelet_type: i32 = 1;
    if getparint("wavelet_type", &mut wavelet_type) == 0 {
        wavelet_type = 1;
    }
    let mut verbose: i32 = 0;
    if getparint("verbose", &mut verbose) == 0 {
        verbose = 0;
    }
    let mut lsource: i32 = 0;
    if getparint("lsource", &mut lsource) == 0 {
        lsource = 0;
    }
    let mut fs: f32 = 0.07f32;
    if getparfloat("fs", &mut fs) == 0 {
        fs = 0.07f32;
    }
    let mut decay: f32 = 50.0f32;
    if getparfloat("decay", &mut decay) == 0 {
        decay = 50.0f32;
    }
    let mut tsec: f32 = 2.048f32;
    if getparfloat("tsec", &mut tsec) == 0 {
        tsec = 2.048f32;
    }
    let mut fref: f32 = 1.0f32;
    if getparfloat("fref", &mut fref) == 0 {
        fref = 1.0f32;
    }
    let mut nw: i32 = 100;
    if getparint("nw", &mut nw) == 0 {
        nw = 100;
    }
    let mut nor: i32 = 100;
    if getparint("nor", &mut nor) == 0 {
        nor = 100;
    }
    let mut np: i32 = 1300;
    if getparint("np", &mut np) == 0 {
        np = 1300;
    }
    let mut p2w: f32 = 5.0f32;
    if getparfloat("p2w", &mut p2w) == 0 {
        p2w = 5.0f32;
    }
    let mut bx: f32 = 0.005f32;
    if getparfloat("bx", &mut bx) == 0 {
        bx = 0.005f32;
    }
    let mut bp: f32 = 0.0f32;
    if getparfloat("bp", &mut bp) == 0 {
        bp = 0.0f32;
    }
    let mut fx: f32 = 0.1f32;
    if getparfloat("fx", &mut fx) == 0 {
        fx = 0.1f32;
    }
    let mut dx: f32 = 0.001f32;
    if getparfloat("dx", &mut dx) == 0 {
        dx = 0.001f32;
    }
    let mut pw1: f32 = 0.0f32;
    if getparfloat("pw1", &mut pw1) == 0 {
        pw1 = 0.0f32;
    }
    let mut pw2: f32 = 0.1f32;
    if getparfloat("pw2", &mut pw2) == 0 {
        pw2 = 0.1f32;
    }
    let mut pw3: f32 = 6.7f32;
    if getparfloat("pw3", &mut pw3) == 0 {
        pw3 = 6.7f32;
    }
    let mut pw4: f32 = 7.0f32;
    if getparfloat("pw4", &mut pw4) == 0 {
        pw4 = 7.0f32;
    }
    let mut h1: f32 = 1.0f32;
    if getparfloat("h1", &mut h1) == 0 {
        h1 = 1.0f32;
    }
    let mut h2: f32 = 0.0f32;
    if getparfloat("h2", &mut h2) == 0 {
        h2 = 0.0f32;
    }
    let mut nx: i32 = 100;
    if getparint("nx", &mut nx) == 0 {
        nx = 100;
    }
    let mut dt: f32 = 0.004f32;
    if getparfloat("dt", &mut dt) == 0 {
        dt = 0.004f32;
    }
    let mut nt: i32 = (tsec / dt) as i32;
    if getparint("nt", &mut nt) == 0 {
        nt = (tsec / dt) as i32;
    }
    let mut nf: i32 = 50;
    if getparint("nf", &mut nf) == 0 {
        nf = 50;
    }
    let mut red_vel: f32 = 5.0f32;
    if getparfloat("red_vel", &mut red_vel) == 0 {
        red_vel = 5.0f32;
    }
    let mut fpeak: f32 = 25.0f32;
    if getparfloat("fpeak", &mut fpeak) == 0 {
        fpeak = 25.0f32;
    }
    let mut tlag: f32 = 0.0f32;
    if getparfloat("tlag", &mut tlag) == 0 {
        tlag = 0.0f32;
    }

    let mut clfile = String::new();
    let mut ctfile = String::new();
    let mut qtfile = String::new();
    let mut qlfile = String::new();
    let mut rhofile = String::new();
    let mut tfile = String::new();
    getparstring("clfile", &mut clfile);
    getparstring("ctfile", &mut ctfile);
    getparstring("qtfile", &mut qtfile);
    getparstring("qlfile", &mut qlfile);
    getparstring("rhofile", &mut rhofile);
    getparstring("tfile", &mut tfile);

    let mut wfp = String::new();
    let mut wfr = String::new();
    let mut wfz = String::new();
    let mut wft = String::new();
    if wtype == 1 {
        getparstring("wfp", &mut wfp);
        getparstring("wfr", &mut wfr);
        getparstring("wfz", &mut wfz);
    } else if wtype == 2 {
        getparstring("wft", &mut wft);
    } else {
        err!("wtype has to be zero or one");
    }
    let (mut m1, mut m2, mut m3) = (0.0f32, 0.0f32, 0.0f32);
    if stype == 1 {
        let (mut delta, mut lambda, mut phis, mut phi) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
        if getparfloat("delta", &mut delta) == 0 {
            err!("if stype==1, delta is a required parameter\n");
        }
        if getparfloat("lambda", &mut lambda) == 0 {
            err!("if stype==1, lambda is a required parameter\n");
        }
        if getparfloat("phis", &mut phis) == 0 {
            err!("if stype==1, phis is a required parameter\n");
        }
        if getparfloat("phi", &mut phi) == 0 {
            err!("if stype==1, phi is a required parameter\n");
        }
        compute_moment_tensor(
            wtype, phi, lambda, delta, phis, m0, &mut m1, &mut m2, &mut m3,
        );
    } else if stype == 2 {
        if getparfloat("m1", &mut m1) == 0 {
            err!("if stype==2, m1 is a required parameter\n");
        }
        if getparfloat("m2", &mut m2) == 0 {
            err!("if stype==2, m2 is a required parameter\n");
        }
        if getparfloat("m3", &mut m3) == 0 {
            err!("if stype==2, m3 is a required parameter\n");
        }
    } else {
        err!("error, stype flag has to be one or two\n");
    }
    let mut layern = 0i32;
    getparint("layern", &mut layern);
    let mut wrefp = 0.0f32;
    getparfloat("wrefp", &mut wrefp);
    let mut wrefs = 0.0f32;
    getparfloat("wrefs", &mut wrefs);
    let mut epsp = 0.0f32;
    getparfloat("epsp", &mut epsp);
    let mut epss = 0.0f32;
    getparfloat("epss", &mut epss);
    let mut sigp = 0.0f32;
    getparfloat("sigp", &mut sigp);
    let mut sigs = 0.0f32;
    getparfloat("sigs", &mut sigs);

    // The C only requests six of the documented model filenames.
    // Receiver, interpolation and filter filenames remain empty and their
    // file branches are therefore unreachable.
    let mut nlayers = 0i32;
    if clfile.is_empty() {
        nlayers = countparval("cl");
    } else {
        getparint("nlayers", &mut nlayers);
    }
    if ctfile.is_empty() && nlayers != countparval("cl") {
        err!("number of p-wave and s-wave velocitieshas to be the same");
    }
    if qlfile.is_empty() && nlayers != countparval("ql") {
        err!("number of p-wave velocities and q-valueshas to be the same");
    }
    if qtfile.is_empty() && nlayers != countparval("qt") {
        err!("number of p-wave velocities and shear q-valueshas to be the same");
    }
    if rhofile.is_empty() && nlayers != countparval("rho") {
        err!("number of p-wave velocities and densitieshas to be the same");
    }
    if tfile.is_empty() && nlayers != countparval("t") {
        err!("number of p-wave velocities and thicknesseshas to be the same");
    }
    if int_type != 1 && int_type != 2 {
        err!("int_type flag has to be one or two");
    }

    let mut nlint = 0i32;
    let mut intlayers = Vec::<i32>::new();
    let mut nintlayers = Vec::<i32>::new();
    let mut intlayth = Vec::<f32>::new();
    if countparval("intlayers") != 0 {
        nlint = countparval("intlayers");
        intlayers = vec![0; nlint as usize];
        getparint("intlayers", &mut intlayers);
    }
    if countparval("nintlayers") != 0 {
        if nlint != countparval("nintlayers") {
            err!("number of values in intlay and nintlay not equal");
        }
        nintlayers = vec![0; nlint as usize];
        getparint("nintlayers", &mut nintlayers);
    }
    if countparval("intlayth") != 0 {
        if nlint != countparval("intlayth") {
            err!("# of values in intlay and intlayth not equal");
        }
        intlayth = vec![0.0; nlint as usize];
        getparfloat("intlayth", &mut intlayth);
    }
    if nlint != 0 {
        for i in 0..nlint as usize {
            nlayers += intlayers[i] - 1;
        }
    }
    // C leaves these indeterminate when rand=1 and the arguments are absent.
    let (mut layer, mut nrand_layers) = (0i32, 0i32);
    let (mut zlayer, mut sdcl, mut sdct) = (0.0f32, 0.0f32, 0.0f32);
    if rand == 1 {
        getparint("layer", &mut layer);
        getparint("nrand_layers", &mut nrand_layers);
        getparfloat("zlayer", &mut zlayer);
        getparfloat("sdcl", &mut sdcl);
        getparfloat("sdct", &mut sdct);
    }
    getparint("nlayers", &mut nlayers);
    fprintf!(stderr(), " Number of Layers %i \n", nlayers);
    let nalloc = (nlayers + nrand_layers) as usize;
    let mut cl = vec![0.0f32; nalloc];
    let mut ct = vec![0.0f32; nalloc];
    let mut ql = vec![0.0f32; nalloc];
    let mut qt = vec![0.0f32; nalloc];
    let mut rho = vec![0.0f32; nalloc];
    let mut t = vec![0.0f32; nalloc];
    let mut lobs = vec![0i32; (nor + 1) as usize];
    lobs[nor as usize] = 0;

    let outfp = if verbose == 2 || verbose == 3 {
        let mut outf = String::from("info");
        getparstring("outf", &mut outf);
        Some(efopen(&outf, "w"))
    } else {
        None
    };
    fprintf!(stderr(), "clfile= %s\n", clfile);
    fprintf!(stderr(), "ctfile= %s\n", ctfile);
    fprintf!(stderr(), "rhofile= %s\n", rhofile);
    fprintf!(stderr(), "tfile= %s\n", tfile);
    fprintf!(stderr(), "qlfile= %s\n", qlfile);
    fprintf!(stderr(), "qtfile= %s\n\n", qtfile);

    // The original unconditionally calls fprintf(outfp,...) here, with a
    // null FILE* unless verbose is 2 or 3. Preserve its SIGSEGV exit.
    if outfp.is_none() {
        null_file_fault();
    }
    let outfp = outfp.unwrap();
    fprintf!(outfp, "clfile= %s\n", clfile);
    fprintf!(outfp, "ctfile= %s\n", ctfile);
    fprintf!(outfp, "rhofile= %s\n", rhofile);
    fprintf!(outfp, "tfile= %s\n", tfile);
    fprintf!(outfp, "qlfile= %s\n", qlfile);
    fprintf!(outfp, "qtfile= %s\n\n", qtfile);
    if !clfile.is_empty() {
        fprintf!(stderr(), "In loop clfile= %s\n", clfile);
    }
    read_model("cl", &clfile, nlayers, &mut cl);
    read_model("ql", &qlfile, nlayers, &mut ql);
    read_model("ct", &ctfile, nlayers, &mut ct);
    read_model("qt", &qtfile, nlayers, &mut qt);
    read_model("rho", &rhofile, nlayers, &mut rho);
    read_model("t", &tfile, nlayers, &mut t);
    getparint("lobs", &mut lobs);

    if nlint != 0 {
        // auxil.c:906-918 needs both endpoints inside the model. The
        // minimally checked C reports this error before reading either one.
        // Follow the C scan: a request beyond the model is never consumed.
        let (mut i, mut ii) = (0i32, 0usize);
        while i < nlayers {
            if intlayers.get(ii) == Some(&i) {
                let n = nintlayers[ii];
                if i < 1 || n < 0 || n >= nlayers - i {
                    err!("interpolation endpoints outside model");
                }
                i += n;
                ii += 1;
            } else {
                i += 1;
            }
        }
        // C reads the next interpolation index after consuming the final
        // entry. Supply an explicit end sentinel for that allocation overread.
        intlayers.push(i32::MAX);
        parameter_interpolation(
            &intlayers,
            &nintlayers,
            &intlayth,
            &mut cl[..nlayers as usize],
            &mut ql[..nlayers as usize],
            &mut ct[..nlayers as usize],
            &mut qt[..nlayers as usize],
            &mut rho[..nlayers as usize],
            &mut t[..nlayers as usize],
        );
    }
    if rand == 1 {
        // auxil.c:1003-1008 allocates only nlayers scratch elements. Walk
        // the generation counters before entering the checked library API.
        // Include the write at nlmax: C tests that limit after the write.
        let mut i = 0;
        let nlmax = (nrand_layers as f32 + nlayers as f32 - layer as f32 + 1.0) as i32;
        'layers: for il in (layer - 1).max(0) as usize..(nlayers - 1).max(0) as usize {
            let icnt = if t[il] <= zlayer {
                1
            } else {
                (t[il] / zlayer) as i32
            };
            for _ in 0..icnt {
                if i >= nlayers {
                    fprintf!(stderr(), "before rand: nlayers=%d ", nlayers);
                    err!("random velocity layer scratch capacity exceeded");
                }
                if i == nlmax {
                    break 'layers;
                }
                i += 1;
            }
        }
        random_velocity_layers(
            &mut nlayers,
            &mut lsource,
            nrand_layers,
            sdcl,
            sdct,
            layer as f32,
            zlayer,
            &mut cl,
            &mut ql,
            &mut ct,
            &mut qt,
            &mut rho,
            &mut t,
        );
    }
    if flt == 1 {
        apply_earth_flattening(
            0.0,
            &mut cl[..nlayers as usize],
            &mut ct[..nlayers as usize],
            &mut rho[..nlayers as usize],
            &mut t[..nlayers as usize],
        );
    }

    let nfilters = countparval("filters_type");
    let mut filters_type = vec![0i32; nfilters as usize];
    getparint("filters_type", &mut filters_type);
    if nfilters != countparval("filters_phase") {
        err!("number of elements infilterstype and phase must be equal");
    }
    let mut filters_phase = vec![0i32; nfilters as usize];
    getparint("filters_phase", &mut filters_phase);
    if nfilters != countparval("dbpo") {
        err!("number of elements in filters_type and dbpo must be equal");
    }
    let mut dbpo = vec![0.0f32; nfilters as usize];
    getparfloat("dbpo", &mut dbpo);
    if nfilters != countparval("f1") {
        err!("number of elements in filters_type and f1 must be equal");
    }
    let mut f1 = vec![0.0f32; nfilters as usize];
    getparfloat("f1", &mut f1);
    if nfilters != countparval("f2") {
        err!("number of elements in filters_type and f2 must be equal");
    }
    let mut f2 = vec![0.0f32; nfilters as usize];
    getparfloat("f2", &mut f2);

    let mut wavefield1 = vec![vec![0.0f32; nt as usize]; nx as usize];
    // The checked library API requires all three matrices even for SH.
    let mut wavefield2 = vec![vec![0.0f32; nt as usize]; nx as usize];
    let mut wavefield3 = vec![vec![0.0f32; nt as usize]; nx as usize];
    let ntc = (tsec / dt) as i32;
    if int_type == 2 {
        bp = 0.0;
    }
    if verbose == 2 || verbose == 3 {
        fprintf!(outfp, "lay # Vp  Vs rho thickness Qp Qs \n");
    }
    fprintf!(stderr(), "lay # Vp  Vs rho thickness Qp Qs \n");
    for i in 0..nlayers as usize {
        fprintf!(
            stderr(),
            "%i    %5.2f %5.2f %5.2f %5.2f   %5.2f  %5.2f \n",
            i as i32,
            cl[i],
            ct[i],
            rho[i],
            t[i],
            ql[i],
            qt[i]
        );
        if verbose == 2 || verbose == 3 {
            fprintf!(
                outfp,
                "%i    %5.2f %5.2f %5.2f %5.2f   %5.2f  %5.2f \n",
                i as i32,
                cl[i],
                ct[i],
                rho[i],
                t[i],
                ql[i],
                qt[i]
            );
        }
    }

    // These checks turn the confirmed C workspace errors into diagnostics
    // before the frozen library's assertions can panic. They preserve the
    // fixed C's exit status and error message; detection happens earlier.
    if nt > fft_length(ntc) {
        err!("nt exceeds FFT trace length");
    }
    if wtype == 1 && lobs[..nor as usize].windows(2).any(|w| w[1] <= w[0]) {
        // psvreflect.c:1053-1112 advances only once at each layer, leaving
        // later repeated or unsorted receivers' workspaces uninitialized.
        err!("receiver workspace not initialized for every receiver");
    }
    if wtype == 1 && nor > 3 {
        // psvreflect.c:180-187 has room for three 100-slowness blocks.
        err!("receiver workspace capacity exceeded");
    }

    // Keep the C's source/receiver indices, including its default lsource=0.
    // The frozen API independently rejects lsource=0 or 1 (interfaces.rs:56),
    // lobs=0 (:61), and unequal SH nx/ntc (:401). The lower-level public
    // reflectivity interfaces reject the same model indices. C's auxil.c
    // accesses lsource-2 and lobs-1, so its boundary behavior also requires
    // library workspace support, beyond removing the interface assertions.
    // cases.txt preserves isolated regressions, including component output.
    // Renumbering these indices would change the original calculation.
    compute_reflectivities(
        int_type,
        verbose,
        wtype,
        wfield,
        vsp,
        flt,
        win,
        ntc,
        nf,
        lsource,
        layern,
        &filters_phase,
        nw,
        np,
        bp,
        tlag,
        red_vel,
        0.0,
        0.0,
        fx,
        dx,
        bx,
        fs,
        decay,
        p2w,
        tsec,
        fref,
        wrefp,
        wrefs,
        epsp,
        epss,
        sigp,
        sigs,
        pw1,
        pw2,
        pw3,
        pw4,
        h1,
        h2,
        m1,
        m2,
        m3,
        fref,
        &lobs[..nor as usize],
        &filters_type,
        &dbpo,
        &f1,
        &f2,
        &cl[..nlayers as usize],
        &ct[..nlayers as usize],
        &ql[..nlayers as usize],
        &qt[..nlayers as usize],
        &rho[..nlayers as usize],
        &t[..nlayers as usize],
        &mut wavefield1,
        &mut wavefield2,
        &mut wavefield3,
        &outfp,
    );
    if verbose == 2 || verbose == 3 {
        efclose(&outfp);
    }

    if wtype == 1 {
        convolve_wavelet(wavelet_type, dt, fpeak, &mut wavefield1);
        convolve_wavelet(wavelet_type, dt, fpeak, &mut wavefield2);
        convolve_wavelet(wavelet_type, dt, fpeak, &mut wavefield3);
        write_wavefield(&wfp, nt, nx, dt, bx, dx, &wavefield1);
        write_wavefield(&wfr, nt, nx, dt, bx, dx, &wavefield2);
        write_wavefield(&wfz, nt, nx, dt, bx, dx, &wavefield3);
    } else if wtype == 2 {
        convolve_wavelet(wavelet_type, dt, fpeak, &mut wavefield1);
        write_wavefield(&wft, nt, nx, dt, bx, dx, &wavefield1);
    }
    // There is no checkpars() call in the original.
    exit(CWP_Exit());
}

fn read_model(name: &str, filename: &str, nlayers: i32, values: &mut Vec<f32>) {
    if !filename.is_empty() {
        let infp = efopen(filename, "r");
        // C uses fread rather than efread: a short file is not an error.
        fread(&mut values[..nlayers as usize], &infp);
        efclose(&infp);
    } else {
        getparfloat(name, values);
    }
}

fn write_wavefield(
    filename: &str,
    nt: i32,
    nx: i32,
    dt: f32,
    bx: f32,
    dx: f32,
    wavefield: &[Vec<f32>],
) {
    if filename.is_empty() {
        return;
    }
    let fp = efopen(filename, "w");
    let mut tr = Trace::new();
    for ix in 0..nx as usize {
        for it in 0..nt as usize {
            tr.data[it] = wavefield[ix][it];
        }
        tr.ns = nt as u16;
        tr.dt = (1000i32.wrapping_mul((1000.0f32 * dt) as i32)) as u16;
        tr.offset = ((bx + ix as f32 * dx) * 1000.0f32) as i32;
        fputtr(&fp, &mut tr);
    }
    efclose(&fp);
}

fn null_file_fault() -> ! {
    unsafe extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
        fn raise(sig: i32) -> i32;
    }
    // Rust installs a SIGSEGV stack-overflow handler. Restore SIG_DFL so
    // raising the original null-FILE fault terminates with the C signal.
    unsafe {
        signal(11, 0);
        raise(11);
    }
    exit(139);
}

// cwp/lib/pfafft.c:npfa uses the first table length >= nmin, with
// saturation at 720720. Keep its lengths here because this program has no
// direct FFT crate dependency and checks the allocation before calling Refl.
fn fft_length(nmin: i32) -> i32 {
    const LENGTHS: [i32; 240] = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 21, 22, 24, 26, 28, 30, 33,
        35, 36, 39, 40, 42, 44, 45, 48, 52, 55, 56, 60, 63, 65, 66, 70, 72, 77, 78, 80, 84, 88, 90,
        91, 99, 104, 105, 110, 112, 117, 120, 126, 130, 132, 140, 143, 144, 154, 156, 165, 168,
        176, 180, 182, 195, 198, 208, 210, 220, 231, 234, 240, 252, 260, 264, 273, 280, 286, 308,
        312, 315, 330, 336, 360, 364, 385, 390, 396, 420, 429, 440, 455, 462, 468, 495, 504, 520,
        528, 546, 560, 572, 585, 616, 624, 630, 660, 693, 715, 720, 728, 770, 780, 792, 819, 840,
        858, 880, 910, 924, 936, 990, 1001, 1008, 1040, 1092, 1144, 1155, 1170, 1232, 1260, 1287,
        1320, 1365, 1386, 1430, 1456, 1540, 1560, 1584, 1638, 1680, 1716, 1820, 1848, 1872, 1980,
        2002, 2145, 2184, 2288, 2310, 2340, 2520, 2574, 2640, 2730, 2772, 2860, 3003, 3080, 3120,
        3276, 3432, 3465, 3640, 3696, 3960, 4004, 4095, 4290, 4368, 4620, 4680, 5005, 5040, 5148,
        5460, 5544, 5720, 6006, 6160, 6435, 6552, 6864, 6930, 7280, 7920, 8008, 8190, 8580, 9009,
        9240, 9360, 10010, 10296, 10920, 11088, 11440, 12012, 12870, 13104, 13860, 15015, 16016,
        16380, 17160, 18018, 18480, 20020, 20592, 21840, 24024, 25740, 27720, 30030, 32760, 34320,
        36036, 40040, 45045, 48048, 51480, 55440, 60060, 65520, 72072, 80080, 90090, 102960,
        120120, 144144, 180180, 240240, 360360, 720720,
    ];
    LENGTHS.into_iter().find(|&n| n >= nmin).unwrap_or(720720)
}
