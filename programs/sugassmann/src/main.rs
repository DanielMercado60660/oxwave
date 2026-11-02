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

//! SUGASSMANN - port of su/main/well_logs/sugassmann.c (Revision 1.3).
//! Credits: UHouston: Chris Liner, 9/23/2009.
use oxwave_core::prelude::*;
/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SUGASSMAN - Model reflectivity change with rock/fluid properties\t",
    "\t\t\t\t\t\t\t\t\t",
    "\tsugassman [optional parameters] > data.su\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " nt=500 \tnumber of time samples\t\t\t\t\t",
    " ntr=200\tnumber of traces\t\t\t\t\t",
    " dt=0.004 \ttime sampling interval in seconds\t\t\t",
    " mode=0\tmodel isolated gassmann refl coefficient\t\t",
    "\t\t=1 embed gassmann RC in random RC series\t\t",
    "\t\t=2 R0 parameter sensitivity output\t\t\t",
    " p=0.15 \tparameter sensitivity test range (if mode=2)\t\t",
    " .... Environment variables ...\t\t\t\t\t",
    " temp=140 \tTemperature in degrees C\t\t\t\t",
    " pres=20 \tPressure in megaPascals\t\t\t\t\t",
    " .... Caprock variables ....\t\t\t\t\t\t",
    " v1=37900 \tcaprock P-wave speed (m/s)\t\t\t\t",
    " r1=44300 \tcaprock mass density (g/cc)\t\t\t\t",
    " .... Reservoir fluid variables ....\t\t\t\t\t",
    " g=0.56\tGas specific gravity 0.56 (methane)-1.8 (condensate)\t",
    " api=50 \tGas specific gravity 10 (heavy)-50 (ultra light)\t",
    " s=35\t\tBrine salinity in ppm/(1000 000\t\t\t\t",
    " so=.7 \tOil saturation (0-1)\t\t\t\t\t",
    " sg=.2 \tGas saturation (0-1)\t\t\t\t\t",
    " .... Reservoir rock frame variables ....\t\t\t\t",
    " kmin=37900 \tBulk modulus (MPa) of rock frame mineral(s) [default=quartz]",
    " mumin=44300 \tShear modulus (MPa) of rock frame mineral(s) [default=quartz]",
    " rmin=2.67 \tMass density (g/cc) of rock frame mineral(s) [default=quartz]",
    " phi=.24 \tRock frame porosity (0-1)\t\t\t\t",
    " a=1 \t\tFitting parameters: Mdry/Mmineral ~ 1/(a + b phi^c)\t",
    " b=15 \t\t... where M is P-wave modulus and defaults are for\t",
    " c=1 \t\t... Glenn sandstone [see Liner (2nd ed, table 26.2)]\t",
    "\th=20 \t\t\tReservoir thickness (m)\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Creates a reflection coefficient series based on Gassmann\t\t",
    " theory of velocity and density for porous elastic media\t\t",
];

#[allow(non_snake_case)] fn Power(x:f64,y:f64)->f64{x.powf(y)}
#[allow(non_snake_case)] fn Sqrt(x:f64)->f64{x.sqrt()}
#[allow(non_snake_case)] fn Abs(x:f64)->f64{x.abs()}
const E:f64=2.71828182845904523536029;

#[allow(non_snake_case)] fn brineModulus(s:f32, t:f32, p:f32)->f32 {
let s=s as f64;
let t=t as f64;
let p=p as f64;
((Power(1402.85 + 1.524*p + 0.003437*Power(p,2.0) - 
	0.000011970000000000002*Power(p,3.0) + 
	(780.0 - 10.0*p + 0.16*Power(p,2.0))*Power(s,1.5) - 
	1820.0*Power(s,2.0) + 4.871*t - 0.0111*p*t + 
	0.00017390000000000003*Power(p,2.0)*t - 
	1.628e-6*Power(p,3.0)*t - 0.04783*Power(t,2.0) + 
	0.0002747*p*Power(t,2.0) - 
	2.135e-6*Power(p,2.0)*Power(t,2.0) + 
	1.2370000000000001e-8*Power(p,3.0)*Power(t,2.0) + 
	0.0001487*Power(t,3.0) - 
	6.502999999999999e-7*p*Power(t,3.0) - 
	1.4550000000000001e-8*Power(p,2.0)*Power(t,3.0) + 
	1.327e-10*Power(p,3.0)*Power(t,3.0) - 
	2.197e-7*Power(t,4.0) + 7.987e-10*p*Power(t,4.0) + 
	5.23e-11*Power(p,2.0)*Power(t,4.0) - 
	4.614e-13*Power(p,3.0)*Power(t,4.0) + 
	s*(1170.0 + 2.6*p - 0.0476*Power(p,2.0) - 9.6*t - 
		0.0029*p*t + 0.055*Power(t,2.0) - 
		0.000085*Power(t,3.0)),2.0)*
	(1.0 + (489.0*p - 0.333*Power(p,2.0) - 80.0*t - 2.0*p*t - 
		0.002*Power(p,2.0)*t - 3.3*Power(t,2.0) + 
		0.016*p*Power(t,2.0) + 0.00175*Power(t,3.0) - 
		0.000013000000000000001*p*Power(t,3.0))/1.0e6 + 
	s*(0.668 + 0.44*s + 
		(300.0*p - 2400.0*p*s + 
		t*(80.0 - 13.0*p - 3300.0*s + 47.0*p*s + 3.0*t))/1.0e6
		)))/
	(1000.0*(1.0 + 0.0494*Power(10.0,
	 -4.0 - (7.786*s)/Power(17.78 + t,0.306))*
	(3676.0*Power(p,0.64) + 
		0.712*p*Power(Abs(-76.71 + t),1.5))))) as f32
}

#[allow(non_snake_case)] fn brineDensity(s:f32, t:f32, p:f32)->f32 {
let s=s as f64;
let t=t as f64;
let p=p as f64;
(1.0 + (489.0*p - 0.333*Power(p,2.0) - 80.0*t - 2.0*p*t - 
	0.002*Power(p,2.0)*t - 3.3*Power(t,2.0) + 
	0.016*p*Power(t,2.0) + 0.00175*Power(t,3.0) - 
	0.000013000000000000001*p*Power(t,3.0))/1.0e6 + 
	s*(0.668 + 0.44*s + 
	(300.0*p - 2400.0*p*s + 
	 t*(80.0 - 13.0*p - 3300.0*s + 47.0*p*s + 3.0*t))/1.0e6)) as f32
}

#[allow(non_snake_case)] fn oilVelocity(api:f32, g:f32, t:f32, p:f32)->f32 {
let api=api as f64;
let g=g as f64;
let t=t as f64;
let p=p as f64;
(4.64*p - 3.7*t + 24932.71072306419*
	Power(1.0/
	((131.5 + api)*(1.0 + 
		0.0020299999999999997*g*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205))*(0.972 + 
		0.00038*Power(17.8 + 
		0.409570869105443*g*Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205) + t,1.175))*
	(2.6 - 141.5/
		((131.5 + api)*
		(1.0 + 0.0020299999999999997*g*
		Power(Power(E,0.02878*api - 0.00377*t)*
			p,1.205))*
		(0.972 + 0.00038*
		Power(17.8 + 
			0.409570869105443*g*
			Sqrt((131.5 + api)*g)*
			Power(Power(E,
			0.02878*api - 0.00377*t)*p,1.205) + 
			t,1.175))))),0.5) + 
	0.0115*p*t*(-1.0 + 4.12*
	Sqrt(-1.0 + 0.007632508833922262*(131.5 + api)*
		(1.0 + 0.0020299999999999997*g*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205))*(0.972 + 
		0.00038*Power(17.8 + 
		0.409570869105443*g*
		Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*
			p,1.205) + t,1.175))))) as f32
}

#[allow(non_snake_case)] fn oilDensity(api:f32, g:f32, t:f32, p:f32)->f32 {
let api=api as f64;
let g=g as f64;
let t=t as f64;
let p=p as f64;
(0.000349*p + (141.5/(131.5 + api) + 
	0.0024359999999999994*Power(g,2.0)*
	Power(Power(E,0.02878*api - 0.00377*t)*p,1.205))/
	(0.972 + 0.00038*Power(17.8 + 
	 0.409570869105443*g*Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205) + t,1.175)) + 
	(0.00277*p - 1.7099999999999998e-7*Power(p,3.0))*
	Power(-1.15 + (141.5/(131.5 + api) + 
	 0.0024359999999999994*Power(g,2.0)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205))/
	(0.972 + 0.00038*
		Power(17.8 + 
		0.409570869105443*g*Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205) + t,1.175)),2.0)) as f32
}

#[allow(non_snake_case)] fn gasModulus(g:f32, t:f32, p:f32)->f32 {
let g=g as f64;
let t=t as f64;
let p=p as f64;
((p*(0.85 - 8.7/
	Power(E,0.65*(1.0 + p/(4.892 - 0.4048*g))) + 
	5.6/(2.0 + p/(4.892 - 0.4048*g)) + 
	27.1/Power(3.5 + p/(4.892 - 0.4048*g),2.0)))/
	(1.0 - (p*(0.03 + 0.00527*
		Power(3.5 - (273.15 + t)/(94.72 + 170.75*g),
		3.0) + (0.0228*
		Power(E,((94.72 + 170.75*g)*
		 Power(p/(4.892 - 0.4048*g),1.2)*
		 (-0.45 - 
			8.0*Power(0.56 - 
			(94.72 + 170.75*g)/(273.15 + t),2.0)))
		/(273.15 + t))*(94.72 + 170.75*g)*
		Power(p/(4.892 - 0.4048*g),0.2)*
		Power(3.85 - 
		(273.15 + t)/(94.72 + 170.75*g),2.0)*
		(-0.45 - 8.0*
		Power(0.56 - 
			(94.72 + 170.75*g)/(273.15 + t),2.0)))/
		(273.15 + t)))/
	((4.892 - 0.4048*g)*
	(-0.52 + (0.642*(273.15 + t))/
		(94.72 + 170.75*g) - 
		(0.007*Power(273.15 + t,4.0))/
		Power(94.72 + 170.75*g,4.0) + 
		0.019*Power(3.85 - 
		(273.15 + t)/(94.72 + 170.75*g),2.0) + 
		(p*(0.03 + 0.00527*
		Power(3.5 - 
			(273.15 + t)/(94.72 + 170.75*g),3.0)))/
		(4.892 - 0.4048*g))))) as f32
}

#[allow(non_snake_case)] fn gasDensity(g:f32, t:f32, p:f32)->f32 {
let g=g as f64;
let t=t as f64;
let p=p as f64;
((3.4638699124410666*g*p)/
	((273.15 + t)*(-0.52 + 
	(0.642*(273.15 + t))/(94.72 + 170.75*g) - 
	(0.007*Power(273.15 + t,4.0))/
	Power(94.72 + 170.75*g,4.0) + 
	0.109*Power(E,((94.72 + 170.75*g)*
		Power(p/(4.892 - 0.4048*g),1.2)*
		(-0.45 - 8.0*
		Power(0.56 - 
		 (94.72 + 170.75*g)/(273.15 + t),2.0)))/
		(273.15 + t))*
	Power(3.85 - (273.15 + t)/(94.72 + 170.75*g),2.0)
	+ (p*(0.03 + 0.00527*
		Power(3.5 - 
		(273.15 + t)/(94.72 + 170.75*g),3.0)))/
	(4.892 - 0.4048*g)))) as f32
}

#[allow(non_snake_case)] fn oilModulus(api:f32,g:f32,t:f32,p:f32)->f32 { (oilDensity(api,g,t,p) as f64*Power(oilVelocity(api,g,t,p) as f64,2.0)/1000.0) as f32 }
#[derive(Clone,Copy)] struct Params {v1:f32,r1:f32,kmin:f32,mumin:f32,rmin:f32,phi:f32,a:f32,b:f32,c:f32,api:f32,g:f32,s:f32,so:f32,sg:f32,temp:f32,pres:f32}
fn rock_r0(q:Params)->f32 {
 let Params{v1,r1,kmin,mumin,rmin,phi,a,b,c,api,g,s,so,sg,temp:t,pres:p}=q;
 if so+sg>1.0 {err!("so + sg > 1...impossible");}
 let brine=brineModulus(s,t,p); let gas=gasModulus(g,t,p); let oil=oilModulus(api,g,t,p);
 let kf=(1.0f64/((1.0f32-so-sg)/brine+so/oil+sg/gas) as f64) as f32;
 let mmin=kmin+4.0*mumin/3.0; let mdry=(mmin as f64/(a as f64+b as f64*Power(phi as f64,c as f64))) as f32;
 let rhs=mdry/(mmin-mdry)+kf/(phi*(mmin-kf)); let msat=mmin*rhs/(1.0+rhs);
 let brine=brineDensity(s,t,p); let gas=gasDensity(g,t,p); let oil=oilDensity(api,g,t,p);
 let fluid=(1.0-so-sg)*brine+so*oil+sg*gas; let rsat=rmin*(1.0-phi)+phi*fluid;
 let vsat=((1000.0f32*msat/rsat) as f64).sqrt() as f32; let ip1=v1*r1; let ip2=rsat*vsat; (ip2-ip1)/(ip2+ip1)
}
fn main(){
 let argv=args(); initargs(&argv); requestdoc(0,SDOC); let mut tr=Trace::new();
 let mut nt=500; getparint("nt",&mut nt); CHECK_NT("nt",nt); tr.ns=nt as u16;
 let mut ntr=200; getparint("ntr",&mut ntr); let mut dt=0.004f32; getparfloat("dt",&mut dt); tr.dt=(dt*1000000.0) as u16;
 let mut mode=0; getparint("mode",&mut mode); if mode!=0 && mode!=1 {mode=0;} let mut it1=nt/2; getparint("it1",&mut it1);
 let mut q=Params{temp:140.0,pres:200.0,v1:3300.0,r1:2.6,g:0.56,api:50.0,s:0.35,so:0.7,sg:0.2,kmin:37900.0,mumin:44300.0,rmin:2.67,phi:0.24,a:1.0,b:15.0,c:1.0};
 getparfloat("temp",&mut q.temp); getparfloat("pres",&mut q.pres); getparfloat("v1",&mut q.v1); getparfloat("r1",&mut q.r1); getparfloat("g",&mut q.g); getparfloat("api",&mut q.api); getparfloat("s",&mut q.s); getparfloat("so",&mut q.so); getparfloat("sg",&mut q.sg); getparfloat("kmin",&mut q.kmin); getparfloat("mumin",&mut q.mumin); getparfloat("rmin",&mut q.rmin); getparfloat("phi",&mut q.phi); getparfloat("a",&mut q.a); getparfloat("b",&mut q.b); getparfloat("c",&mut q.c);
 let mut h=20.0f32; getparfloat("h",&mut h); let mut p=0.15f32; getparfloat("p",&mut p);
 if mode==0 {tr.ns=200; let x=q.sg; for it in 0..200 {let y=(1.0f32-p)*x+2.0*(it as f32)*p*x/200.0; let mut varied=q; varied.sg=y; tr.data[it]=100.0*rock_r0(varied)/rock_r0(q);} tr.f1=(1.0-p)*x; tr.d1=2.0*p*x/200.0; puttr(&mut tr);}
 exit(CWP_Exit());
}
