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

//! Differential tests for psplot/lib/basic.c; credits: Dave Hale, Craig Artley, Robert Krueger.
use oxwave_cref::XorShift;
use oxwave_f_basic as r;
use std::{
    ffi::{c_char, CString},
    fs::OpenOptions,
    io::{Read, Seek, SeekFrom},
    os::fd::AsRawFd,
    sync::Mutex,
};
static LOCK: Mutex<()> = Mutex::new(());
#[link(name = "psplot", kind = "static")]
extern "C" {
    fn beginps();
    fn begineps();
    fn endps();
    fn endeps();
    fn showpage();
    fn gsave();
    fn grestore();
    fn newpath();
    fn closepath();
    fn clip();
    fn stroke();
    fn fill();
    fn translate(x: f32, y: f32);
    fn scale(x: f32, y: f32);
    fn moveto(x: f32, y: f32);
    fn rmoveto(x: f32, y: f32);
    fn lineto(x: f32, y: f32);
    fn rlineto(x: f32, y: f32);
    fn rotate(x: f32);
    fn setgray(x: f32);
    fn setlinewidth(x: f32);
    fn fontheight(name: *const c_char, size: f32) -> f32;
    fn fontwidth(name: *const c_char, size: f32) -> f32;
    fn fontcapheight(name: *const c_char, size: f32) -> f32;
    fn fontxheight(name: *const c_char, size: f32) -> f32;
    fn fontdescender(name: *const c_char, size: f32) -> f32;
    fn fontascender(name: *const c_char, size: f32) -> f32;
    fn boundingbox(a: i32, b: i32, c: i32, d: i32);
    fn newpage(name: *const c_char, n: i32);
    fn concat(m: *const f32);
    fn setrgbcolor(r: f32, g: f32, b: f32);
    fn setcymkcolor(c: f32, m: f32, y: f32, k: f32);
    fn setcolor(name: *const c_char);
    fn setlinejoin(n: i32);
    fn setdash(d: *const f32, n: i32, o: f32);
    fn arc(x: f32, y: f32, r: f32, a: f32, b: f32);
    fn show(s: *const c_char);
    fn justshow(j: f32, s: *const c_char);
    fn setfont(s: *const c_char, size: f32);
    fn fontbbox(s: *const c_char, size: f32, b: *mut f32);
    fn image(w: i32, h: i32, b: i32, m: *const f32, s: *const u8);
    fn rgbimage(w: i32, h: i32, b: i32, m: *const f32, s: *const u8);
    fn polyline(x: *const f32, y: *const f32, n: i32);
    fn markto(x: f32, y: f32, i: i32, size: f32);
    fn rectclip(x: f32, y: f32, w: f32, h: f32);
    fn rectfill(x: f32, y: f32, w: f32, h: f32);
    fn rectstroke(x: f32, y: f32, w: f32, h: f32);
    fn dup(fd: i32) -> i32;
    fn dup2(a: i32, b: i32) -> i32;
    fn close(fd: i32) -> i32;
    fn fflush(f: *mut std::ffi::c_void) -> i32;
}
fn compare<F: Fn(bool)>(run: F) {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let path = std::env::var_os("CARGO_TARGET_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../.targets/lib-f_basic")
        })
        .join(format!("oxwave-basic-{}", std::process::id()));
    let mut file = OpenOptions::new()
        .create(true)
        .truncate(true)
        .read(true)
        .write(true)
        .open(&path)
        .unwrap();
    unsafe {
        fflush(std::ptr::null_mut());
    }
    oxwave_core::cio::flush_all();
    let saved = unsafe { dup(1) };
    let mut outputs = Vec::new();
    for c in [false, true] {
        file.set_len(0).unwrap();
        file.seek(SeekFrom::Start(0)).unwrap();
        unsafe {
            assert_eq!(dup2(file.as_raw_fd(), 1), 1);
        }
        run(c);
        oxwave_core::cio::flush_all();
        unsafe {
            fflush(std::ptr::null_mut());
            dup2(saved, 1);
        }
        file.seek(SeekFrom::Start(0)).unwrap();
        let mut v = Vec::new();
        file.read_to_end(&mut v).unwrap();
        outputs.push(v);
    }
    unsafe {
        dup2(saved, 1);
        close(saved);
    }
    std::fs::remove_file(path).unwrap();
    let clean = |v: &[u8]| {
        String::from_utf8_lossy(v)
            .lines()
            .filter(|line| !line.starts_with("test test_"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let a = clean(&outputs[0]);
    let b = clean(&outputs[1]);
    if a != b {
        let pos = a
            .bytes()
            .zip(b.bytes())
            .position(|(x, y)| x != y)
            .unwrap_or(a.len().min(b.len()));
        panic!(
            "output differs at byte {pos}: {:?} vs {:?}",
            a.get(pos..pos.saturating_add(100).min(a.len())),
            b.get(pos..pos.saturating_add(100).min(b.len()))
        );
    }
}
fn val(r: &mut XorShift) -> f32 {
    match r.next_u32() % 32 {
        0 => 0.0,
        1 => -0.0,
        2 => 1e-20,
        3 => -1e20,
        _ => r.uniform_f32(-10000.0, 10000.0),
    }
}
#[test]
fn test_beginps() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { beginps() }
            } else {
                r::beginps()
            }
        }
    });
}
#[test]
fn test_begineps() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { begineps() }
            } else {
                r::begineps()
            }
        }
    });
}
#[test]
fn test_endps() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { endps() }
            } else {
                r::endps()
            }
        }
    });
}
#[test]
fn test_endeps() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { endeps() }
            } else {
                r::endeps()
            }
        }
    });
}
#[test]
fn test_showpage() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { showpage() }
            } else {
                r::showpage()
            }
        }
    });
}
#[test]
fn test_gsave() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { gsave() }
            } else {
                r::gsave()
            }
        }
    });
}
#[test]
fn test_grestore() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { grestore() }
            } else {
                r::grestore()
            }
        }
    });
}
#[test]
fn test_newpath() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { newpath() }
            } else {
                r::newpath()
            }
        }
    });
}
#[test]
fn test_closepath() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { closepath() }
            } else {
                r::closepath()
            }
        }
    });
}
#[test]
fn test_clip() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { clip() }
            } else {
                r::clip()
            }
        }
    });
}
#[test]
fn test_stroke() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { stroke() }
            } else {
                r::stroke()
            }
        }
    });
}
#[test]
fn test_fill() {
    compare(|c| {
        for _ in 0..2000 {
            if c {
                unsafe { fill() }
            } else {
                r::fill()
            }
        }
    });
}
#[test]
fn test_translate() {
    compare(|c| {
        let mut rng = XorShift::new(47);
        for _ in 0..2000 {
            let x = val(&mut rng);
            let y = val(&mut rng);
            if c {
                unsafe { translate(x, y) }
            } else {
                r::translate(x, y)
            }
        }
    });
}
#[test]
fn test_scale() {
    compare(|c| {
        let mut rng = XorShift::new(47);
        for _ in 0..2000 {
            let x = val(&mut rng);
            let y = val(&mut rng);
            if c {
                unsafe { scale(x, y) }
            } else {
                r::scale(x, y)
            }
        }
    });
}
#[test]
fn test_moveto() {
    compare(|c| {
        let mut rng = XorShift::new(47);
        for _ in 0..2000 {
            let x = val(&mut rng);
            let y = val(&mut rng);
            if c {
                unsafe { moveto(x, y) }
            } else {
                r::moveto(x, y)
            }
        }
    });
}
#[test]
fn test_rmoveto() {
    compare(|c| {
        let mut rng = XorShift::new(47);
        for _ in 0..2000 {
            let x = val(&mut rng);
            let y = val(&mut rng);
            if c {
                unsafe { rmoveto(x, y) }
            } else {
                r::rmoveto(x, y)
            }
        }
    });
}
#[test]
fn test_lineto() {
    compare(|c| {
        let mut rng = XorShift::new(47);
        for _ in 0..2000 {
            let x = val(&mut rng);
            let y = val(&mut rng);
            if c {
                unsafe { lineto(x, y) }
            } else {
                r::lineto(x, y)
            }
        }
    });
}
#[test]
fn test_rlineto() {
    compare(|c| {
        let mut rng = XorShift::new(47);
        for _ in 0..2000 {
            let x = val(&mut rng);
            let y = val(&mut rng);
            if c {
                unsafe { rlineto(x, y) }
            } else {
                r::rlineto(x, y)
            }
        }
    });
}
#[test]
fn test_rotate() {
    compare(|c| {
        let mut rng = XorShift::new(48);
        for _ in 0..2000 {
            let x = val(&mut rng);
            if c {
                unsafe { rotate(x) }
            } else {
                r::rotate(x)
            }
        }
    });
}
#[test]
fn test_setgray() {
    compare(|c| {
        let mut rng = XorShift::new(48);
        for _ in 0..2000 {
            let x = val(&mut rng);
            if c {
                unsafe { setgray(x) }
            } else {
                r::setgray(x)
            }
        }
    });
}
#[test]
fn test_setlinewidth() {
    compare(|c| {
        let mut rng = XorShift::new(48);
        for _ in 0..2000 {
            let x = val(&mut rng);
            if c {
                unsafe { setlinewidth(x) }
            } else {
                r::setlinewidth(x)
            }
        }
    });
}

#[test]
fn test_boundingbox() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let a = q.next_u32() as i32;
            let b = q.next_u32() as i32;
            let d = q.next_u32() as i32;
            let e = q.next_u32() as i32;
            if c {
                unsafe {
                    boundingbox(a, b, d, e);
                    begineps()
                }
            } else {
                r::boundingbox(a, b, d, e);
                r::begineps()
            }
        }
    });
}
#[test]
fn test_newpage() {
    compare(|c| {
        let s = CString::new("page").unwrap();
        for i in 0..2000 {
            if c {
                unsafe { newpage(s.as_ptr(), i) }
            } else {
                r::newpage("page", i)
            }
        }
    });
}
#[test]
fn test_concat() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let m = [
                val(&mut q),
                val(&mut q),
                val(&mut q),
                val(&mut q),
                val(&mut q),
                val(&mut q),
            ];
            if c {
                unsafe { concat(m.as_ptr()) }
            } else {
                r::concat(&m)
            }
        }
    });
}
#[test]
fn test_setrgbcolor() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let (a, b, d) = (val(&mut q), val(&mut q), val(&mut q));
            if c {
                unsafe { setrgbcolor(a, b, d) }
            } else {
                r::setrgbcolor(a, b, d)
            }
        }
    });
}
#[test]
fn test_setcymkcolor() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let (a, b, d, e) = (val(&mut q), val(&mut q), val(&mut q), val(&mut q));
            if c {
                unsafe { setcymkcolor(a, b, d, e) }
            } else {
                r::setcymkcolor(a, b, d, e)
            }
        }
    });
}
#[test]
fn test_setcolor() {
    compare(|c| {
        let mut q = XorShift::new(3);
        for i in 0..2000 {
            let s = if i % 5 == 0 {
                "GREEN".to_string()
            } else if i % 5 == 1 {
                "unknown".to_string()
            } else {
                format!("#{:06x}", q.next_u32() & 0xffffff)
            };
            if c {
                unsafe { setcolor(CString::new(s).unwrap().as_ptr()) }
            } else {
                r::setcolor(&s)
            }
        }
    });
}
#[test]
fn test_setlinejoin() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let v = q.next_u32() as i32;
            if c {
                unsafe { setlinejoin(v) }
            } else {
                r::setlinejoin(v)
            }
        }
    });
}
#[test]
fn test_setdash() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for i in 0..2000 {
            let v: Vec<_> = (0..i % 17).map(|_| val(&mut q)).collect();
            let o = val(&mut q);
            if c {
                unsafe { setdash(v.as_ptr(), v.len() as i32, o) }
            } else {
                r::setdash(&v, o)
            }
        }
    });
}
#[test]
fn test_arc() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let (a, b, d, e, f) = (
                val(&mut q),
                val(&mut q),
                val(&mut q),
                val(&mut q),
                val(&mut q),
            );
            if c {
                unsafe { arc(a, b, d, e, f) }
            } else {
                r::arc(a, b, d, e, f)
            }
        }
    });
}
#[test]
fn test_show() {
    compare(|c| {
        for i in 0..2000 {
            let s = format!("string ({i})");
            if c {
                unsafe { show(CString::new(s).unwrap().as_ptr()) }
            } else {
                r::show(&s)
            }
        }
    });
}
#[test]
fn test_justshow() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for i in 0..2000 {
            let s = format!("string ({i})");
            let j = val(&mut q);
            if c {
                unsafe { justshow(j, CString::new(s).unwrap().as_ptr()) }
            } else {
                r::justshow(j, &s)
            }
        }
    });
}
#[test]
fn test_setfont() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let x = val(&mut q);
            if c {
                unsafe { setfont(c"Helvetica".as_ptr(), x) }
            } else {
                r::setfont("Helvetica", x)
            }
        }
    });
}
const FONTS: &[&str] = &[
    "AvantGarde-Book",
    "Courier",
    "Helvetica",
    "Helvetica-Narrow-Oblique",
    "NewCentrySchlbk-Bold",
    "Palatino-Roman",
    "Times-Roman",
    "Times-Italic",
    "Unknown",
];
#[test]
fn test_fontheight() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let a = r::fontheight(name, s);
        let b = unsafe { fontheight(CString::new(name).unwrap().as_ptr(), s) };
        assert_eq!(a.to_bits(), b.to_bits(), "{name} {s}");
    }
}
#[test]
fn test_fontwidth() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let a = r::fontwidth(name, s);
        let b = unsafe { fontwidth(CString::new(name).unwrap().as_ptr(), s) };
        assert_eq!(a.to_bits(), b.to_bits(), "{name} {s}");
    }
}
#[test]
fn test_fontcapheight() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let a = r::fontcapheight(name, s);
        let b = unsafe { fontcapheight(CString::new(name).unwrap().as_ptr(), s) };
        assert_eq!(a.to_bits(), b.to_bits(), "{name} {s}");
    }
}
#[test]
fn test_fontxheight() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let a = r::fontxheight(name, s);
        let b = unsafe { fontxheight(CString::new(name).unwrap().as_ptr(), s) };
        assert_eq!(a.to_bits(), b.to_bits(), "{name} {s}");
    }
}
#[test]
fn test_fontdescender() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let a = r::fontdescender(name, s);
        let b = unsafe { fontdescender(CString::new(name).unwrap().as_ptr(), s) };
        assert_eq!(a.to_bits(), b.to_bits(), "{name} {s}");
    }
}
#[test]
fn test_fontascender() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let a = r::fontascender(name, s);
        let b = unsafe { fontascender(CString::new(name).unwrap().as_ptr(), s) };
        assert_eq!(a.to_bits(), b.to_bits(), "{name} {s}");
    }
}

#[test]
fn test_fontbbox() {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let mut q = XorShift::new(1);
    for i in 0..5000 {
        let name = FONTS[i % FONTS.len()];
        let s = val(&mut q);
        let mut a = [0.0; 4];
        let mut b = [0.0; 4];
        r::fontbbox(name, s, &mut a);
        unsafe { fontbbox(CString::new(name).unwrap().as_ptr(), s, b.as_mut_ptr()) };
        assert_eq!(a, b);
    }
}
#[test]
fn test_image() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for i in 0..2000 {
            let w = [1, 2, 7, 16, 39, 67][i % 6];
            let h = 1 + i % 3;
            let b = [1, 2, 4, 8, 12][i % 5];
            let count = (1 + (w * b * 1 - 1) / 8) * h;
            let data: Vec<u8> = (0..count).map(|_| q.next_u32() as u8).collect();
            let m = [val(&mut q), 0.0, 0.0, val(&mut q), 1.0, -1.0];
            if c {
                unsafe { image(w as i32, h as i32, b as i32, m.as_ptr(), data.as_ptr()) }
            } else {
                r::image(w as i32, h as i32, b as i32, &m, &data)
            }
        }
    });
}
#[test]
fn test_rgbimage() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for i in 0..2000 {
            let w = [1, 2, 7, 16, 39, 67][i % 6];
            let h = 1 + i % 3;
            let b = [1, 2, 4, 8, 12][i % 5];
            let count = (1 + (w * b * 3 - 1) / 8) * h;
            let data: Vec<u8> = (0..count).map(|_| q.next_u32() as u8).collect();
            let m = [val(&mut q), 0.0, 0.0, val(&mut q), 1.0, -1.0];
            if c {
                unsafe { rgbimage(w as i32, h as i32, b as i32, m.as_ptr(), data.as_ptr()) }
            } else {
                r::rgbimage(w as i32, h as i32, b as i32, &m, &data)
            }
        }
    });
}
#[test]
fn test_polyline() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for i in 0..2000 {
            let n = [1, 2, 7, 16, 199, 200, 201, 257][i % 8];
            let x: Vec<_> = (0..n).map(|_| val(&mut q)).collect();
            let y: Vec<_> = (0..n).map(|_| val(&mut q)).collect();
            if c {
                unsafe { polyline(x.as_ptr(), y.as_ptr(), n as i32) }
            } else {
                r::polyline(&x, &y)
            }
        }
    });
}
#[test]
fn test_markto() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for i in 0..2000 {
            let (x, y, z) = (val(&mut q), val(&mut q), val(&mut q));
            let n = i % 31 - 10;
            if c {
                unsafe { markto(x, y, n, z) }
            } else {
                r::markto(x, y, n, z)
            }
        }
    });
}
#[test]
fn test_rectclip() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let (a, b, d, e) = (val(&mut q), val(&mut q), val(&mut q), val(&mut q));
            if c {
                unsafe { rectclip(a, b, d, e) }
            } else {
                r::rectclip(a, b, d, e)
            }
        }
    });
}
#[test]
fn test_rectfill() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let (a, b, d, e) = (val(&mut q), val(&mut q), val(&mut q), val(&mut q));
            if c {
                unsafe { rectfill(a, b, d, e) }
            } else {
                r::rectfill(a, b, d, e)
            }
        }
    });
}
#[test]
fn test_rectstroke() {
    compare(|c| {
        let mut q = XorShift::new(1);
        for _ in 0..2000 {
            let (a, b, d, e) = (val(&mut q), val(&mut q), val(&mut q), val(&mut q));
            if c {
                unsafe { rectstroke(a, b, d, e) }
            } else {
                r::rectstroke(a, b, d, e)
            }
        }
    });
}
