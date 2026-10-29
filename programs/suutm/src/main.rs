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
//! SUUTM - UTM projection; port of su/main/headers/suutm.c (revision 1.5).
//! Author: Nils Maercklin. Projection adapted from Geo-Coordinates-UTM-0.06
//! by G. Crookham (CPAN, March 2007).

use oxwave_core::prelude::*;
use oxwave_cwp::cwp::{nint, PI};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " SUUTM - UTM projection of longitude and latitude in SU trace headers  ",
    "                                                                       ",
    " suutm <stdin >stdout [optional parameters]                            ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "    counit=(from header) input coordinate units code:                  ",
    "                    =1: length (meters or feet; no UTM projection)     ",
    "                    =2: seconds of arc                                 ",
    "                    =3: decimal degrees                                ",
    "                    =4: degrees, minutes, seconds                      ",
    "    idx=23          reference ellipsoid index (default is WGS 1984)    ",
    "    a=(from idx)    user-specified semimajor axis of ellipsoid         ",
    "    f=(from idx)    user-specified flattening of ellipsoid             ",
    "    zkey=           if set, header key to store UTM zone number        ",
    "    verbose=0       =1: echo ellipsoid parameters                      ",
    "                                                                       ",
    "    lon0=           central meridian for TM projection in degrees      ",
    "                    (default uses the 60 standard UTM longitude zones) ",
    "    xoff=500000     false Easting (default: UTM)                       ",
    "    ysoff=10000000  false Northing, southern hemisphere (default: UTM) ",
    "    ynoff=0         false Northing, northern hemisphere (default: UTM) ",
    "                                                                       ",
    " Notes:                                                                ",
    "    Universal Transverse Mercator (UTM) coordinates are defined between",
    "    latitudes 80S (-80) and 84N (84). Longitude values must be between ",
    "    -180 degrees (west) and 179.999... degrees (east).                 ",
    "                                                                       ",
    "    Latitudes are read from sy and gy (N positive), and longitudes     ",
    "    are read from sx and gx (E positive).                              ",
    "    The UTM zone is determined from the receiver coordinates gy and gx.",
    "                                                                       ",
    "    Use suazimuth to calculate shot-receiver azimuths and offsets.     ",
    "                                                                       ",
    " Reference ellipsoids:                                                 ",
    "    An ellipsoid may be specified by its semimajor axis a and its      ",
    "    flattening f, or one of the following ellipsoids may be selected   ",
    "    by its index idx (semimajor axes in meters):                       ",
    "     0  Sphere with radius of 6371000 m                                ",
    "     1  Airy 1830                                                      ",
    "     2  Australian National 1965                                       ",
    "     3  Bessel 1841 (Ethiopia, Indonesia, Japan, Korea)                ",
    "     4  Bessel 1841 (Namibia)                                          ",
    "     5  Clarke 1866                                                    ",
    "     6  Clarke 1880                                                    ",
    "     7  Everest (Brunei, E. Malaysia)                                  ",
    "     8  Everest (India 1830)                                           ",
    "     9  Everest (India 1956)                                           ",
    "    10  Everest (Pakistan)                                             ",
    "    11  Everest (W. Malaysia, Singapore 1948)                          ",
    "    12  Everest (W. Malaysia 1969)                                     ",
    "    13  Geodetic Reference System 1980 (GRS 1980)                      ",
    "    14  Helmert 1906                                                   ",
    "    15  Hough 1960                                                     ",
    "    16  Indonesian 1974                                                ",
    "    17  International 1924 / Hayford 1909                              ",
    "    18  Krassovsky 1940                                                ",
    "    19  Modified Airy                                                  ",
    "    20  Modified Fischer 1960                                          ",
    "    21  South American 1969                                            ",
    "    22  World Geodetic System 1972 (WGS 1972)                          ",
    "    23  World Geodetic System 1984 (WGS 1984) / NAD 1983               ",
    "                                                                       ",
];

fn ellipsoid(idx: usize) -> (f64, f64, &'static str) {
    let e = [
        (6371000.,0.,"Sphere"), (6377563.396,1./299.3249646,"Airy 1830"),
        (6378160.,1./298.25,"Australian National 1965"),
        (6377397.155,1./299.1528128,"Bessel 1841 (Ethiopia, Indonesia, Japan, Korea)"),
        (6377483.865,1./299.1528128,"Bessel 1841 (Namibia)"),
        (6378206.4,1./294.9786982,"Clarke 1866"), (6378249.145,1./293.465,"Clarke 1880"),
        (6377298.556,1./300.80178,"Everest (Brunei, E. Malaysia)"),
        (6377276.345,1./300.80178,"Everest (India 1830)"),
        (6377301.243,1./300.80178,"Everest (India 1956)"),
        (6377309.613,1./300.80178,"Everest (Pakistan)"),
        (6377304.063,1./300.80178,"Everest (W. Malaysia, Singapore 1948)"),
        (6377295.664,1./300.80178,"Everest (W. Malaysia 1969)"),
        (6378137.,1./298.257222101,"Geodetic Reference System 1980 (GRS 1980)"),
        (6378200.,1./298.3,"Helmert 1906"), (6378270.,1./297.,"Hough 1960"),
        (6378160.,1./298.247,"Indonesian 1974"),
        (6378388.,1./297.,"International 1924 / Hayford 1909"),
        (6378245.,1./298.3,"Krassovsky 1940"), (6377340.189,1./299.3249646,"Modified Airy"),
        (6378155.,1./298.3,"Modified Fischer 1960"), (6378160.,1./298.25,"South American 1969"),
        (6378135.,1./298.26,"WGS 1972"), (6378137.,1./298.257223563,"WGS 1984 / NAD 1983")
    ];
    e[idx]
}

fn zone(lat: f64, mut lon: f64) -> i16 {
    lon = (lon+180.)-((lon+180.)/360.).floor()*360.-180.;
    let mut z = ((lon+180.)/6.) as i16 + 1;
    if lat >= 56. && lat < 64. && lon >= 3. && lon < 12. { z = 32; }
    if lat >= 72. && lat < 84. {
        if lon >= 0. && lon < 9. { z = 31; }
        else if lon >= 9. && lon < 21. { z = 33; }
        else if lon >= 21. && lon < 33. { z = 35; }
        else if lon >= 33. && lon < 42. { z = 37; }
    }
    z
}

fn decode(x: f64, y: f64, scalco: i16, counit: i16) -> (f64,f64) {
    let mut factor = if scalco == 0 { 1. } else { scalco as f64 };
    if factor < 0. { factor = -1./factor; }
    match counit {
        2 => (factor*x/3600., factor*y/3600.),
        4 => {
            fn dms(v: f64) -> f64 {
                let d = (v/10000.).trunc();
                let m0 = 10000.*(v/10000.-d);
                let m = (m0/100.).trunc();
                let s = 100.*(m0/100.-m);
                d+m/60.+s/3600.
            }
            (dms(x*factor),dms(y*factor))
        },
        _ => (factor*x, factor*y)
    }
}

fn project(a: f64, f: f64, mut lat: f64, mut lon: f64, mut lon0: f64) -> (f64,f64) {
    lat *= PI/180.; lon *= PI/180.; lon0 *= PI/180.;
    let e2 = f*(2.-f);
    let ep2 = e2/(1.-e2);
    let k0 = 0.9996;
    let cn = a/(1.-e2*lat.sin()*lat.sin()).sqrt();
    let ct = lat.tan()*lat.tan();
    let cc = ep2*lat.cos()*lat.cos();
    let ca = lat.cos()*(lon-lon0);
    let cm = a*((1.-e2/4.-3.*e2*e2/64.-5.*e2*e2*e2/256.)*lat
        -(3.*e2/8.+3.*e2*e2/32.+45.*e2*e2*e2/1024.)*(2.*lat).sin()
        +(15.*e2*e2/256.+45.*e2*e2*e2/1024.)*(4.*lat).sin()
        -(35.*e2*e2*e2/3072.)*(6.*lat).sin());
    let x = k0*cn*(ca+(1.-ct+cc)*ca*ca*ca/6.
        +(5.-18.*ct+ct*ct+72.*cc-58.*ep2)*ca*ca*ca*ca*ca/120.);
    let y = k0*(cm+cn*lat.tan()*(ca*ca/2.+(5.-ct+9.*cc+4.*cc*cc)*ca*ca*ca*ca/24.
        +(61.-58.*ct+ct*ct+600.*cc-330.*ep2)*ca*ca*ca*ca*ca*ca/720.));
    (x,y)
}

fn setval(t: &str, d: f64) -> Value {
    let mut v = Value::default();
    match t {
        "s" => err!("can't change char header word"),
        "h" => v.set_h(nint(d) as i16), "u" => v.set_u(nint(d) as u16),
        "l" => v.set_l(nint(d) as i64), "v" => v.set_v(nint(d) as u64),
        "i" => v.set_i(nint(d)), "p" => v.set_p(nint(d) as u32),
        "f" => v.set_f(d as f32), "d" => v.set_d(d),
        _ => err!("unknown type %s", t),
    }
    v
}

fn main() {
    let argv = args(); initargs(&argv); requestdoc(1, SDOC);
    let mut idx=23; let mut verbose=0; let mut counit=0; let mut zkey=String::new();
    let mut xoff=500000.; let mut ynoff=0.; let mut ysoff=10000000.; let mut lon0=0.;
    getparint("idx", &mut idx); getparint("verbose", &mut verbose); getparint("counit", &mut counit);
    let has_zkey=getparstring("zkey", &mut zkey)!=0;
    getpardouble("xoff", &mut xoff); getpardouble("ynoff", &mut ynoff); getpardouble("ysoff", &mut ysoff);
    let userlon0=getpardouble("lon0", &mut lon0)!=0;
    if idx<0 || idx>=24 { err!("unknown ellipsoid idx=%d", idx); }
    if userlon0 && (lon0 < -180. || lon0 > 180.) { err!("lon0=%g must be between -180 and 180", lon0); }
    let mut ztype=""; let mut zindex=0;
    if has_zkey { ztype=hdtype(&zkey); zindex=getindex(&zkey); }
    let (mut a,mut f,mut name)=ellipsoid(idx as usize);
    if getpardouble("a", &mut a)!=0 {
        if a<=0. { err!("semimajor axis a=%g must be positive", a); } name="User-specified";
    }
    if getpardouble("f", &mut f)!=0 {
        if f<0. { err!("flattening f=%g must be non-negative", f); } name="User-specified";
    }
    checkpars();
    if verbose!=0 {
        warn!("ellipsoid: %s", name);
        if userlon0 { warn!("a=%.12g f=%.12g lon0=%g", a,f,lon0); }
        else { warn!("a=%.12g f=%.12g", a,f); }
    }
    let mut tr=Trace::new();
    while gettr(&mut tr)!=0 {
        if counit!=0 { tr.counit=counit as i16; }
        let (mut sx,sy)=decode(tr.sx as f64,tr.sy as f64,tr.scalco,tr.counit);
        let (mut gx,gy)=decode(tr.gx as f64,tr.gy as f64,tr.scalco,tr.counit);
        if tr.counit>=2 && tr.counit<=4 {
            if gx==180. { gx=-180.; } if sx==180. { sx=-180.; }
            if gy < -89.999 || gy > 89.999 || gx < -180. || gx >= 180. {
                err!("invalid coordinates: tracl=%d gx=%g gy=%g", tr.tracl,gx,gy);
            }
            if sy < -89.999 || sy > 89.999 || sx < -180. || sx >= 180. {
                err!("invalid coordinates: tracl=%d sx=%g sy=%g", tr.tracl,sx,sy);
            }
            if !userlon0 {
                let mut z=zone(gy,gx); lon0=((z as i32-1)*6-180+3) as f64;
                if has_zkey { if gy<0. { z *= -1; } let v=setval(ztype,z as f64); puthval(&mut tr,zindex,&v); }
            }
            let (mut gx,mut gy)=project(a,f,gy,gx,lon0);
            let (mut sx,mut sy)=project(a,f,sy,sx,lon0);
            gx+=xoff; sx+=xoff;
            gy+=if gy<0. {ysoff} else {ynoff};
            sy+=if gy<0. {ysoff} else {ynoff};
            tr.counit=1; tr.scalco=1;
            tr.gx=nint(gx); tr.gy=nint(gy); tr.sx=nint(sx); tr.sy=nint(sy);
        }
        puttr(&mut tr);
    }
    exit(CWP_Exit());
}
