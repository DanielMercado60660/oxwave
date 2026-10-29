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

//! UTMCONV - port of par/main/data_conversion/utmconv.c, revision 1.4.
//! Author: Nils Maercklin. Projection routines adopted from Geo-Coordinates-UTM-0.06 by G. Crookham.
//! References: NIMA (2000); J. P. Snyder (1987).

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "                                                                       ",
    " UTMCONV - CONVert longitude and latitude to UTM, and vice versa       ",
    "                                                                       ",
    " utmconv <stdin >stdout [optional parameters]                          ",
    "                                                                       ",
    " Optional parameters:                                                  ",
    "    idx=23          reference ellipsoid index (default is WGS 1984)    ",
    "    format=%.3f     output number format (printf style for one float)  ",
    "    a=(from idx)    user-specified semimajor axis of ellipsoid         ",
    "    f=(from idx)    user-specified flattening of ellipsoid             ",
    "    letter=0        =1: use UTM letter designator for latitude/Northing",
    "    invert=0        =0: convert latitude and longitude to UTM          ",
    "                    =1: convert UTM to latitude and longitude          ",
    "    verbose=0       =1: echo parameters and number of converted coords ",
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
    "    Input and output is in ASCII format. For a conversion from lon/lat ",
    "    to UTM (invert=0), input is a two-column table of longitude and    ",
    "    latitude in decimal degrees. Output is a three-column table of     ",
    "    UTM Easting, UTM Northing, and UTM zone. The zone is given either  ",
    "    by the zone number only (default, negative on southern hemisphere) ",
    "    or by the positive zone number plus a letter designator (letter=1).",
    "                                                                       ",
    " Example:                                                              ",
    "    Convert 40.822N, 14.125E to UTM with zone number and letter,       ",
    "    output values rounded to nearest integer:                          ",
    "        echo 14.125 40.822 | utmconv letter=1 format=%.0f              ",
    "    The output is \"426213 4519366 33T\" (Easting, Northing, UTM zone).",
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
const E: [(f64, f64, &str); 24] = [
    (6371000., 0., "Sphere"),
    (6377563.396, 1. / 299.3249646, "Airy 1830"),
    (6378160., 1. / 298.25, "Australian National 1965"),
    (
        6377397.155,
        1. / 299.1528128,
        "Bessel 1841 (Ethiopia, Indonesia, Japan, Korea)",
    ),
    (6377483.865, 1. / 299.1528128, "Bessel 1841 (Namibia)"),
    (6378206.4, 1. / 294.9786982, "Clarke 1866"),
    (6378249.145, 1. / 293.465, "Clarke 1880"),
    (6377298.556, 1. / 300.80178, "Everest (Brunei, E. Malaysia)"),
    (6377276.345, 1. / 300.80178, "Everest (India 1830)"),
    (6377301.243, 1. / 300.80178, "Everest (India 1956)"),
    (6377309.613, 1. / 300.80178, "Everest (Pakistan)"),
    (
        6377304.063,
        1. / 300.80178,
        "Everest (W. Malaysia, Singapore 1948)",
    ),
    (6377295.664, 1. / 300.80178, "Everest (W. Malaysia 1969)"),
    (
        6378137.,
        1. / 298.257222101,
        "Geodetic Reference System 1980 (GRS 1980)",
    ),
    (6378200., 1. / 298.3, "Helmert 1906"),
    (6378270., 1. / 297., "Hough 1960"),
    (6378160., 1. / 298.247, "Indonesian 1974"),
    (6378388., 1. / 297., "International 1924 / Hayford 1909"),
    (6378245., 1. / 298.3, "Krassovsky 1940"),
    (6377340.189, 1. / 299.3249646, "Modified Airy"),
    (6378155., 1. / 298.3, "Modified Fischer 1960"),
    (6378160., 1. / 298.25, "South American 1969"),
    (6378135., 1. / 298.26, "WGS 1972"),
    (6378137., 1. / 298.257223563, "WGS 1984 / NAD 1983"),
];
fn letter(lat: f64) -> i32 {
    if lat > 84. {
        return 90;
    }
    if lat >= 72. {
        return 88;
    }
    let lim = [
        64., 56., 48., 40., 32., 24., 16., 8., 0., -8., -16., -24., -32., -40., -48., -56., -64.,
        -72., -80.,
    ];
    let ls = b"WVUTSRQPNMLKJHGFEDC";
    for i in 0..lim.len() {
        if lat >= lim[i] {
            return ls[i] as i32;
        }
    }
    65
}
fn zone(lat: f64, mut lon: f64) -> i16 {
    lon = (lon + 180.) - ((lon + 180.) / 360.).floor() * 360. - 180.;
    let mut z = ((lon + 180.) / 6.) as i16 + 1;
    if lat >= 56. && lat < 64. && lon >= 3. && lon < 12. {
        z = 32
    }
    if lat >= 72. && lat < 84. {
        if lon >= 0. && lon < 9. {
            z = 31
        } else if lon >= 9. && lon < 21. {
            z = 33
        } else if lon >= 21. && lon < 33. {
            z = 35
        } else if lon >= 33. && lon < 42. {
            z = 37
        }
    }
    z
}
fn forward(a: f64, f: f64, mut lat: f64, mut lon: f64, mut lon0: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    lat *= pi / 180.;
    lon *= pi / 180.;
    lon0 *= pi / 180.;
    let e2 = f * (2. - f);
    let ep2 = e2 / (1. - e2);
    let k0 = 0.9996;
    let cn = a / (1. - e2 * lat.sin() * lat.sin()).sqrt();
    let ct = lat.tan() * lat.tan();
    let cc = ep2 * lat.cos() * lat.cos();
    let ca = lat.cos() * (lon - lon0);
    let cm = a
        * ((1. - e2 / 4. - 3. * e2 * e2 / 64. - 5. * e2 * e2 * e2 / 256.) * lat
            - (3. * e2 / 8. + 3. * e2 * e2 / 32. + 45. * e2 * e2 * e2 / 1024.) * (2. * lat).sin()
            + (15. * e2 * e2 / 256. + 45. * e2 * e2 * e2 / 1024.) * (4. * lat).sin()
            - (35. * e2 * e2 * e2 / 3072.) * (6. * lat).sin());
    let x = k0
        * cn
        * (ca
            + (1. - ct + cc) * ca * ca * ca / 6.
            + (5. - 18. * ct + ct * ct + 72. * cc - 58. * ep2) * ca * ca * ca * ca * ca / 120.);
    let y = k0
        * (cm
            + cn * lat.tan()
                * (ca * ca / 2.
                    + (5. - ct + 9. * cc + 4. * cc * cc) * ca * ca * ca * ca / 24.
                    + (61. - 58. * ct + ct * ct + 600. * cc - 330. * ep2)
                        * ca
                        * ca
                        * ca
                        * ca
                        * ca
                        * ca
                        / 720.));
    (x, y)
}
fn inverse(a: f64, f: f64, x: f64, y: f64, lon0: f64) -> (f64, f64) {
    let e2 = f * (2. - f);
    let ep2 = e2 / (1. - e2);
    let e1 = (1. - (1. - e2).sqrt()) / (1. + (1. - e2).sqrt());
    let k0 = 0.9996;
    let cm = y / k0;
    let mu = cm / (a * (1. - e2 / 4. - 3. * e2 * e2 / 64. - 5. * e2 * e2 * e2 / 256.));
    let phi = mu
        + (3. * e1 / 2. - 27. * e1 * e1 * e1 / 32.) * (2. * mu).sin()
        + (21. * e1 * e1 / 16. - 55. * e1 * e1 * e1 * e1 / 32.) * (4. * mu).sin()
        + (151. * e1 * e1 * e1 / 96.) * (6. * mu).sin();
    let cn = a / (1. - e2 * phi.sin() * phi.sin()).sqrt();
    let ct = phi.tan() * phi.tan();
    let cc = ep2 * phi.cos() * phi.cos();
    let cr = a * (1. - e2) / (1. - e2 * phi.sin() * phi.sin()).powf(1.5);
    let cd = x / (cn * k0);
    let mut lat = phi
        - (cn * phi.tan() / cr)
            * (cd * cd / 2.
                - (5. + 3. * ct + 10. * cc - 4. * cc * cc - 9. * ep2) * cd * cd * cd * cd / 24.
                + (61. + 90. * ct + 298. * cc + 45. * ct * ct - 252. * ep2 - 3. * cc * cc)
                    * cd
                    * cd
                    * cd
                    * cd
                    * cd
                    * cd
                    / 720.);
    lat *= 180. / std::f64::consts::PI;
    let mut lon = (cd - (1. + 2. * ct + cc) * cd * cd * cd / 6.
        + (5. - 2. * cc + 28. * ct - 3. * cc * cc + 8. * ep2 + 24. * ct * ct)
            * cd
            * cd
            * cd
            * cd
            * cd
            / 120.)
        / phi.cos();
    lon = lon0 + (lon * 180. / std::f64::consts::PI);
    (lat, lon)
}
fn main() {
    let argv = args();
    initargs(&argv);
    requestdoc(1, SDOC);
    let mut format = String::new();
    if getparstring("format", &mut format) == 0 {
        format = "%.3f".into()
    }
    let mut idx = 23;
    let mut verbose = 0;
    let mut invert = 0;
    let mut useletter = 0;
    getparint("idx", &mut idx);
    getparint("verbose", &mut verbose);
    getparint("invert", &mut invert);
    getparint("letter", &mut useletter);
    let mut xoff = 500000.;
    let mut ynoff = 0.;
    let mut ysoff = 10000000.;
    let mut lon0 = 0.;
    getpardouble("xoff", &mut xoff);
    getpardouble("ynoff", &mut ynoff);
    getpardouble("ysoff", &mut ysoff);
    let userlon0 = getpardouble("lon0", &mut lon0) != 0;
    if idx < 0 || idx >= 24 {
        err!("unknown ellipsoid idx=%d", idx)
    }
    if userlon0 && (lon0 < -180. || lon0 > 180.) {
        err!("lon0=%g must be between -180 and 180", lon0)
    }
    if format.len() > 20 {
        err!("format string too long")
    }
    let (mut a, mut f, mut name) = E[idx as usize];
    if getpardouble("a", &mut a) != 0 {
        if a <= 0. {
            err!("semimajor axis a=%g must be positive", a)
        }
        name = "User-specified";
    }
    if getpardouble("f", &mut f) != 0 {
        if f < 0. {
            err!("flattening f=%g must be non-negative", f)
        }
        name = "User-specified";
    }
    checkpars();
    if verbose != 0 {
        warn!("ellipsoid: %s", name);
        if userlon0 {
            warn!("a=%.12g f=%.12g lon0=%g", a, f, lon0)
        } else {
            warn!("a=%.12g f=%.12g", a, f)
        }
    }
    let input = stdin();
    let mut n = 0;
    loop {
        let Some(first) = fscanf_f64(&input) else {
            break;
        };
        let Some(second) = fscanf_f64(&input) else {
            break;
        };
        n += 1;
        if invert != 0 {
            let Some(mut tmp) = fscanf_f32(&input) else {
                if useletter != 0 {
                    err!("missing zone number/letter in record %d", n)
                } else {
                    err!("missing zone number for record %d", n)
                }
            };
            if useletter != 0 {
                let c = fgetc(&input);
                if c < 0 {
                    err!("missing zone number/letter in record %d", n)
                }
                if (65..=90).contains(&c) {
                    if c < 78 {
                        tmp = -tmp
                    }
                } else if (97..=122).contains(&c) {
                    if c < 110 {
                        tmp = -tmp
                    }
                } else {
                    err!("missing or invalid zone letter in record %d", n)
                }
            }
            let z = tmp as i16;
            if !userlon0 {
                if (z as i32).abs() < 1 || (z as i32).abs() > 60 {
                    err!("invalid zone=%d in record %d", z as i32, n)
                }
                lon0 = (((z as i32).abs() - 1) * 6 - 180 + 3) as f64;
            }
            let (lat, lon) = inverse(
                a,
                f,
                first - xoff,
                second - if z < 0 { ysoff } else { ynoff },
                lon0,
            );
            let fmt = format!("{} {}\n", format, format);
            printf!(&fmt, lon, lat);
        } else {
            let mut lon = first;
            let lat = second;
            if lon == 180. {
                lon = -180.
            }
            if lat < -89.999 || lat > 89.999 || lon < -180. || lon >= 180. {
                err!(
                    "invalid coordinates in record %d: lat=%g lon=%g",
                    n,
                    lat,
                    lon
                )
            }
            let mut z = zone(lat, lon);
            let l = letter(lat);
            if !userlon0 {
                lon0 = ((z as i32 - 1) * 6 - 180 + 3) as f64;
            }
            let (mut x, mut y) = forward(a, f, lat, lon, lon0);
            x += xoff;
            y += if lat < 0. { ysoff } else { ynoff };
            if useletter != 0 {
                let fmt = format!("{} {} %2d%c\n", format, format);
                printf!(&fmt, x, y, z as i32, l);
            } else {
                if lat < 0. {
                    z = -z
                }
                let fmt = format!("{} {} %3d\n", format, format);
                printf!(&fmt, x, y, z as i32);
            }
        }
    }
    if verbose != 0 {
        warn!(
            "converted %d coordinate pair%s",
            n,
            if n == 1 { "" } else { "s" }
        )
    }
    exit(CWP_Exit());
}
