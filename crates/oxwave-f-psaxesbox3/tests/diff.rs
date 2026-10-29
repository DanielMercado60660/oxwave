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

//! Differential tests for psplot/lib/psAxesBox3.c.
//! Credits: Dave Hale and Zhiming Li, Colorado School of Mines.
use oxwave_cref::XorShift;
use oxwave_f_psaxesbox3::{psAxesBBox3 as rust_bbox, psAxesBox3 as rust_box};
use std::ffi::{c_char, CString};
#[link(name="psplot", kind="static")]
extern "C" {
    fn psAxesBBox3(x:f32,y:f32,w:f32,h:f32,lf:*const c_char,ls:f32,tf:*const c_char,ts:f32,style:i32,b:*mut i32);
    fn psAxesBox3(x:f32,y:f32,w:f32,h:f32,a:f32,b:f32,c:f32,d:f32,e:f32,f:f32,n:i32,g:i32,l:*const c_char,a2:f32,b2:f32,c2:f32,d2:f32,e2:f32,f2:f32,n2:i32,g2:i32,l2:*const c_char,lf:*const c_char,ls:f32,t:*const c_char,tf:*const c_char,ts:f32,style:i32,t2:*const c_char);
    fn dup(fd:i32)->i32; fn dup2(a:i32,b:i32)->i32; fn close(fd:i32)->i32; fn fflush(stream:*mut std::ffi::c_void)->i32;
}
#[test]
fn bbox_matches() {
    let mut rng=XorShift::new(123); let font=CString::new("Helvetica").unwrap();
    for i in 0..10000 { let mut v=[0.0f32;6]; for x in &mut v {*x=((rng.next_f64()*2.0-1.0)*1000.0) as f32;} let style=(i%3) as i32;
        let mut c=[0;4]; let mut r=[0;4]; unsafe {psAxesBBox3(v[0],v[1],v[2],v[3],font.as_ptr(),v[4],font.as_ptr(),v[5],style,c.as_mut_ptr())};
        rust_bbox(v[0],v[1],v[2],v[3],"Helvetica",v[4],"Helvetica",v[5],style,&mut r); assert_eq!(r,c);
    }
}
#[test]
fn box_matches() {
    if std::env::var_os("OXWAVE_BOX_CHILD").is_none() {
        let result = std::process::Command::new(std::env::current_exe().unwrap()).args(["--exact", "box_matches"]).env("OXWAVE_BOX_CHILD", "1").output().unwrap();
        assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stdout));
        return;
    }
    use std::os::fd::AsRawFd;
    use std::io::{Read,Seek,SeekFrom};
    let path=std::env::var("CARGO_TARGET_DIR").unwrap_or_else(|_| ".".into());
    let path=std::path::Path::new(&path).join("psaxesbox3-diff-output");
    let mut file=std::fs::OpenOptions::new().create(true).truncate(true).read(true).write(true).open(&path).unwrap();
    let mut rng=XorShift::new(57); let font=CString::new("Helvetica").unwrap(); let label=CString::new("axis").unwrap();
    for i in 0..2000 {
        let x=(rng.next_f64()*100.0) as f32; let style=i%2; let grid=i%4; let size=50.0+(rng.next_f64()*300.0) as f32; let step=if i%3==0 {0.0} else {1.0}; let ntic=[1,2,3,8,17][i as usize%5];
        let saved=unsafe {fflush(std::ptr::null_mut()); dup(1)}; unsafe {dup2(file.as_raw_fd(),1);}
        unsafe {psAxesBox3(x,-x,size,100.0,-3.0,5.0,0.0,0.0,step,-3.0,ntic,grid,label.as_ptr(),-2.0,4.0,0.0,0.0,step,-2.0,ntic,grid,label.as_ptr(),font.as_ptr(),12.0,label.as_ptr(),font.as_ptr(),14.0,style,label.as_ptr()); fflush(std::ptr::null_mut());}
        let split=file.stream_position().unwrap();
        rust_box(x,-x,size,100.0,-3.0,5.0,0.0,0.0,step,-3.0,ntic,grid,"axis",-2.0,4.0,0.0,0.0,step,-2.0,ntic,grid,"axis","Helvetica",12.0,"axis","Helvetica",14.0,style,"axis"); oxwave_core::cio::flush_all();
        unsafe {dup2(saved,1);close(saved);} let end=file.stream_position().unwrap(); file.seek(SeekFrom::Start(0)).unwrap(); let mut data=Vec::new(); file.read_to_end(&mut data).unwrap(); assert_eq!(&data[..split as usize],&data[split as usize..end as usize],"case {i}"); file.set_len(0).unwrap(); file.seek(SeekFrom::Start(0)).unwrap();
    }
    std::fs::remove_file(path).unwrap();
}

