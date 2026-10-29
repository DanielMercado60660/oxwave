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
//! SUK1K2FILTER - symmetric box-like wavenumber filter.
//! Port of su/main/filters/suk1k2filter.c (Revision 1.13).
//! Credits: CWP: John Stockwell, November 1995.

use oxwave_core::prelude::*;
use oxwave_cwp::complex::{Complex, crmul};
use oxwave_cwp::cwp::nint;
use oxwave_pfa_fft::*;
use oxwave_tempfile::temporary_filename;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    " \t\t\t\t\t\t\t\t\t",
    " SUK1K2FILTER - symmetric box-like K-domain filter defined by the\t",
    "\t\t  cartesian product of two sin^2-tapered polygonal\t",
    "\t\t  filters defined in k1 and k2\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    "     suk1k2filter <infile >outfile [optional parameters]\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " k1=val1,val2,...\tarray of K1 filter wavenumbers\t\t\t",
    " k2=val1,val2,...\tarray of K2 filter wavenumbers\t\t\t",
    " amps1=a1,a2,...\tarray of K1 filter amplitudes\t\t\t",
    " amps2=a1,a2,...\tarray of K2 filter amplitudes\t\t\t",
    " d1=tr.d1 or 1.0\tsampling interval in first (fast) dimension\t",
    " d2=tr.d1 or 1.0\tsampling interval in second (slow) dimension\t",
    " quad=0\t\t=0 all four quandrants\t\t\t\t",
    "\t\t\t=1 (quadrants 1 and 4) \t\t\t\t",
    "\t\t\t=2 (quadrants 2 and 3) \t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Defaults:\t\t\t\t\t\t\t\t",
    " k1=.10*(nyq1),.15*(nyq1),.45*(nyq1),.50*(nyq1)\t\t\t",
    " k2=.10*(nyq2),.15*(nyq2),.45*(nyq2),.50*(nyq2)\t\t\t",
    " amps1=0.,1.,...,1.,0.  trapezoid-like bandpass filter\t\t\t",
    " amps2=0.,1.,...,1.,0.  trapezoid-like bandpass filter\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " The nyquist wavenumbers, nyq1 and nyq2, are computed internally.\t",
    "\t\t\t\t\t\t\t\t\t",
    " verbose=0\tverbose = 1 echoes information\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tmpdir= \t if non-empty, use the value as a directory path\t",
    "\t\t prefix for storing temporary files; else if the\t",
    "\t         the CWP_TMPDIR environment variable is set use\t\t",
    "\t         its value for the path; else use tmpfile()\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " The filter is assumed to be symmetric, to yield real output\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Because the data are assumed to be purely spatial (i.e. non-seismic), ",
    " the data are assumed to have trace id (30), corresponding to (z,x) data",
    "\t\t\t\t\t\t\t\t\t",
    " The relation: w = 2 pi F is well known for frequency, but there\t",
    " doesn't seem to be a commonly used letter corresponding to F for the\t",
    " spatial conjugate transform variables.  We use K1 and K2 for this.\t",
    " More specifically we assume a phase:\t\t\t\t\t",
    "\t\t-i(k1 x1 + k2 x2) = -2 pi i(K1 x1 + K2 x2).\t\t",
    " and K1, K2 define our respective wavenumbers.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
];

unsafe extern "C" { fn access(path: *const i8, mode: i32) -> i32; }

fn polygonal_filter(f: &[f32], a: &[f32], nfft: i32, dt: f32, filter: &mut [f32]) {
    let nf = nfft/2+1;
    let onfft = (1.0/nfft as f64) as f32;
    let df = onfft/dt;
    let mut bins = vec![0i32; f.len()];
    for i in 0..f.len() { bins[i] = nint((f[i]/df) as f64); if bins[i] > nf-1 { bins[i] = nf-1; } }
    for j in 0..bins[0] { filter[j as usize] = a[0]*onfft; }
    let mut taper = 0;
    for i in 0..f.len()-1 {
        if a[i] != a[i+1] {
            taper += 1;
            for j in bins[i]..=bins[i+1] {
                let c = (1.57079632679490 / (bins[i+1]-bins[i]+2) as f64) as f32;
                let rise = a[i] < a[i+1];
                let t = if rise { j-bins[i]+1 } else { bins[i+1]-j+1 };
                let s = ((c*t as f32) as f64).sin() as f32;
                let diff = if rise { a[i+1]-a[i] } else { a[i]-a[i+1] };
                let base = if rise { a[i] } else { a[i+1] };
                filter[j as usize] = (base+diff*s*s)*onfft;
            }
        } else {
            for j in bins[i]+if taper == 0 {0} else {1}..=bins[i+1] { filter[j as usize] = a[i]*onfft; }
        }
    }
    for j in bins[f.len()-1]+1..nf { filter[j as usize] = a[f.len()-1]*onfft; }
}

fn parameters(axis: i32, nyq: f32) -> (Vec<f32>, Vec<f32>) {
    let key = if axis == 1 { "k1" } else { "k2" };
    let akey = if axis == 1 { "amps1" } else { "amps2" };
    let n = countparval(key);
    let mut k = if n != 0 { vec![0.0; n as usize] } else { vec![(0.10*nyq as f64) as f32, (0.15*nyq as f64) as f32, (0.45*nyq as f64) as f32, (0.50*nyq as f64) as f32] };
    if n != 0 { getparfloat(key, &mut k); }
    if k.len() < 2 { warn!("Only %d value defining filter", k.len() as i32); }
    for i in 0..k.len()-1 { if k[i] < 0.0 || k[i] > k[i+1] { err!("Bad filter parameters"); } }
    let na = countparval(akey);
    let mut a = if na != 0 { vec![0.0; na as usize] } else { let mut x = vec![1.0; k.len()]; x[0]=0.0; let last=x.len()-1; x[last]=0.0; x };
    if na != 0 { getparfloat(akey, &mut a); }
    if a.len() != k.len() { if axis == 1 { err!("number of k1 values must = number of amps1 values"); } else { err!("number of k2 values must = number of amps2 values"); } }
    let mut count = 0i32;
    for &x in &a { count = (count as f32+x) as i32; if x < 0.0 { if axis == 1 { err!("amps1 values must be positive"); } else { err!("amps2 values must be positive"); } } }
    if count == 0 { if axis == 1 { err!("All amps1 values are zero"); } else { err!("All amps2 values are zero"); } }
    let mut differences=0; for i in 0..a.len()-1 { if a[i] != a[i+1] { differences+=1; } }
    if differences == 0 { if axis == 1 { warn!("All amps1 values are the same"); } else { warn!("All amps2 values are the same"); } }
    (k,a)
}

fn main() {
    let argv=args(); initargs(&argv); requestdoc(1, SDOC);
    let mut tr=Trace::new(); if gettr(&mut tr)==0 { err!("can't get first trace"); }
    if tr.trid != 30 { warn!("tr.trid = %d", tr.trid as i32); }
    let nx1=tr.ns as usize;
    let mut dx1=0.0; let mut dx2=0.0;
    if getparfloat("d1", &mut dx1)==0 { if tr.d1 != 0.0 { dx1=tr.d1; } else { dx1=1.0; warn!("tr.d1 not set, assuming d1=1.0"); } }
    if getparfloat("d2", &mut dx2)==0 { if tr.d2 != 0.0 { dx2=tr.d2; } else { dx2=1.0; warn!("tr.d2 not set, assuming d2=1.0"); } }
    let mut verbose=0; getparint("verbose", &mut verbose);
    let mut tmpdir=String::new(); if getparstring("tmpdir", &mut tmpdir)==0 { tmpdir=std::env::var("CWP_TMPDIR").unwrap_or_default(); }
    if !tmpdir.is_empty() { let p=std::ffi::CString::new(tmpdir.as_str()).unwrap_or_default(); if unsafe { access(p.as_ptr(),2) } != 0 { err!("you can't write in %s (or it doesn't exist)", tmpdir); } }
    let nyq1=(0.5/dx1 as f64) as f32; let nyq2=(0.5/dx2 as f64) as f32;
    let mut names=Vec::new();
    let (tracefp,headerfp)=if tmpdir.is_empty() { if verbose != 0 { warn!("using tmpfile() call"); } (etmpfile(),etmpfile()) } else { let a=temporary_filename(&tmpdir); let mut b=temporary_filename(&tmpdir); if a == b { b.push_str("_header"); } let af=efopen(&a,"w+"); let bf=efopen(&b,"w+"); names.push(a); names.push(b); if verbose != 0 { warn!("putting temporary files in %s",tmpdir); } (af,bf) };
    let mut headers=Vec::new(); loop { headers.push(tr.hdr); efwrite(&tr.data[..nx1], &tracefp); if gettr(&mut tr)==0 { break; } }
    let nx2=headers.len(); let n1=npfaro(nx1 as i32,2*nx1 as i32); let n2=npfa(nx2 as i32);
    if n1 >= SU_NFLTS as i32 || n1 >= 720720 { err!("Padded nx1=%d--too big",n1); }
    if n2 >= SU_NFLTS as i32 || n2 >= 720720 { err!("Padded nx2=%d--too big",n2); }
    let (k1,a1)=parameters(1,nyq1); let (k2,a2)=parameters(2,nyq2);
    let mut quad=0; getparint("quad",&mut quad);
    let n1u=n1 as usize; let n2u=n2 as usize; let nk1=n1u/2+1; let nk2=n2u/2+1;
    let mut rt=vec![0.0f32;n1u*n2u]; let mut ct=vec![Complex::default();nk1*n2u];
    let mut kfilter=vec![0.0f32;n1u*n2u]; let mut f1=vec![0.0;nk1]; let mut f2=vec![0.0;nk2];
    polygonal_filter(&k1,&a1,n1,dx1,&mut f1); polygonal_filter(&k2,&a2,n2,dx2,&mut f2);
    if quad==0 || quad==1 { for i2 in 0..nk2 { for i1 in 0..nk1 { kfilter[i2*n1u+i1]=f1[i1]*f2[i2]; } } }
    if quad==0 || quad==2 { for i2 in nk2..n2u { for i1 in 0..nk1 { kfilter[i2*n1u+i1]=f1[i1]*f2[n2u-i2]; } } }
    erewind(&tracefp); for i in 0..nx2 { efread(&mut rt[i*n1u..i*n1u+nx1],&tracefp); }
    pfa2rc(-1,1,n1,nx2 as i32,&rt,&mut ct); pfa2cc(-1,2,nk1 as i32,n2,&mut ct);
    for i2 in 0..n2u { for i1 in 0..nk1 { let j=i2*nk1+i1; ct[j]=crmul(ct[j],kfilter[i2*n1u+i1]); } }
    pfa2cc(1,2,nk1 as i32,n2,&mut ct); pfa2cr(1,1,n1,nx2 as i32,&ct,&mut rt);
    for i in 0..nx2 { tr.hdr=headers[i]; tr.data[..nx1].copy_from_slice(&rt[i*n1u..i*n1u+nx1]); puttr(&mut tr); }
    efclose(&headerfp); efclose(&tracefp); for name in names { eremove(&name); } exit(CWP_Exit());
}
