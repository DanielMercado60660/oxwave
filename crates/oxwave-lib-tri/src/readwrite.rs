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

//! Port of tri/lib/readwrite.c.
//! Credits: Jack K. Cohen, Center for Wave Phenomena, September 1990.
//! Modified by Dave Hale, November 1990 (binary format and attributes),
//! and Craig Artley, April 1994 (exterior edge and vertex uses).

use super::*;
use std::io::Write;

/// C address relocation tree, from tri/include/triP.h.
#[repr(C)]
pub struct Address {
    pub oaddress: *mut c_void,
    pub naddress: *mut c_void,
    pub aLeft: *mut Address,
    pub aRight: *mut Address,
}

extern "C" {
    fn fread(p: *mut c_void, size: usize, count: usize, fp: *mut c_void) -> usize;
    fn fwrite(p: *const c_void, size: usize, count: usize, fp: *mut c_void) -> usize;
}

fn fail(message: &str) -> ! {
    let _ = std::io::stderr().write_all(message.as_bytes());
    unsafe { exit(1) }
}

/// C: Address *updateAddressTree(Address *p, address_t oa, int *found);
/// # Safety
/// The tree must contain valid Address nodes, uniquely mutable for the call.
/// The caller owns the returned tree. Address keys are compared, never dereferenced.
pub unsafe fn updateAddressTree(
    mut p: *mut Address,
    oa: *mut c_void,
    found: &mut i32,
) -> *mut Address {
    unsafe {
        if p.is_null() {
            p = allocate();
            (*p).oaddress = oa;
            (*p).aLeft = std::ptr::null_mut();
            (*p).aRight = std::ptr::null_mut();
            *found = 0;
        } else {
            // C off_t subtraction on this 64-bit target.
            let sign = (oa as i64).wrapping_sub((*p).oaddress as i64);
            if sign < 0 {
                (*p).aLeft = updateAddressTree((*p).aLeft, oa, found);
            } else if sign > 0 {
                (*p).aRight = updateAddressTree((*p).aRight, oa, found);
            } else {
                *found = 1;
            }
        }
        p
    }
}

unsafe fn address(mut p: *mut Address, oa: *mut c_void) -> *mut Address {
    unsafe {
        while !p.is_null() {
            let sign = (oa as i64).wrapping_sub((*p).oaddress as i64);
            if sign < 0 {
                p = (*p).aLeft;
            } else if sign > 0 {
                p = (*p).aRight;
            } else {
                return p;
            }
        }
        p
    }
}

// DONODE and DOFIELD: allocate forward references only once, then relocate.
unsafe fn node<T>(head: &mut *mut Address, old: *mut T) -> *mut T {
    unsafe {
        let mut found = 0;
        *head = updateAddressTree(*head, old.cast(), &mut found);
        let a = address(*head, old.cast());
        if found == 0 {
            (*a).naddress = allocate::<T>().cast();
        }
        (*a).naddress.cast()
    }
}

// SU leaves these temporary trees allocated. They are no longer referenced by
// the returned graph, so release them without changing the model or stream.
unsafe fn release_tree(p: *mut Address) {
    unsafe {
        if !p.is_null() {
            release_tree((*p).aLeft);
            release_tree((*p).aRight);
            free(p.cast());
        }
    }
}

unsafe fn read_record<T>(p: *mut T, fp: *mut c_void, message: &str) {
    unsafe {
        if fread(p.cast(), std::mem::size_of::<T>(), 1, fp) != 1 {
            fail(message);
        }
    }
}
unsafe fn read_attributes(p: &mut *mut c_void, size: usize, fp: *mut c_void, message: &str) {
    unsafe {
        if !p.is_null() && size != 0 {
            *p = calloc(1, size);
            if fread(*p, size, 1, fp) != 1 {
                fail(message);
            }
        } else {
            *p = std::ptr::null_mut();
        }
    }
}

/// C: Model *readModel(FILE *fp);
/// Reads SU's native binary model format, including attributes and forward links.
/// # Safety
/// fp must be a readable C FILE stream containing a valid model, or EOF.
/// The caller owns the returned graph. Callbacks are reset to NULL as in C.
pub unsafe fn readModel(fp: *mut c_void) -> *mut Model {
    unsafe {
        let mut m: *mut Model = std::ptr::null_mut();
        let mut fh = std::ptr::null_mut();
        let mut euh = std::ptr::null_mut();
        let mut vuh = std::ptr::null_mut();
        let mut eh = std::ptr::null_mut();
        let mut vh = std::ptr::null_mut();
        let mut kind = 0i32;
        while fread((&raw mut kind).cast(), std::mem::size_of::<i32>(), 1, fp) != 0 {
            match kind {
                0 => {
                    m = allocate();
                    // The first address is present but unused in the C reader.
                    fread(m.cast(), std::mem::size_of::<*mut c_void>(), 1, fp);
                    read_record(m, fp, "error reading Model\n");
                    let mut found = 0;
                    fh = updateAddressTree(fh, (*m).f.cast(), &mut found);
                    if found != 0 {
                        fail("Model's first face already exists\n");
                    }
                    let a = address(fh, (*m).f.cast());
                    (*a).naddress = allocate::<Face>().cast();
                    (*m).f = (*a).naddress.cast();
                    (*m).vAdd = None;
                    (*m).vDel = None;
                    (*m).tAdd = None;
                    (*m).tDel = None;
                    read_attributes(
                        &mut (*m).ma,
                        (*m).sma,
                        fp,
                        "error reading model attributes\n",
                    );
                }
                1 => {
                    let mut old = std::ptr::null_mut::<Face>();
                    fread(
                        (&raw mut old).cast(),
                        std::mem::size_of::<*mut Face>(),
                        1,
                        fp,
                    );
                    let f = node(&mut fh, old);
                    read_record(f, fp, "error reading Face\n");
                    (*f).m = m;
                    (*f).fPrev = node(&mut fh, (*f).fPrev);
                    (*f).fNext = node(&mut fh, (*f).fNext);
                    (*f).eu = node(&mut euh, (*f).eu);
                    read_attributes(&mut (*f).fa, (*m).sfa, fp, "error reading FaceAttributes\n");
                }
                2 => {
                    let mut old = std::ptr::null_mut::<EdgeUse>();
                    fread(
                        (&raw mut old).cast(),
                        std::mem::size_of::<*mut EdgeUse>(),
                        1,
                        fp,
                    );
                    let eu = node(&mut euh, old);
                    read_record(eu, fp, "error reading EdgeUse\n");
                    if !(*eu).f.is_null() {
                        (*eu).f = node(&mut fh, (*eu).f);
                    }
                    (*eu).vu = node(&mut vuh, (*eu).vu);
                    (*eu).euMate = node(&mut euh, (*eu).euMate);
                    (*eu).euCW = node(&mut euh, (*eu).euCW);
                    (*eu).euCCW = node(&mut euh, (*eu).euCCW);
                    (*eu).e = node(&mut eh, (*eu).e);
                    read_attributes(
                        &mut (*eu).eua,
                        (*m).seua,
                        fp,
                        "error reading edge-use attributes\n",
                    );
                }
                3 => {
                    let mut old = std::ptr::null_mut::<Edge>();
                    fread(
                        (&raw mut old).cast(),
                        std::mem::size_of::<*mut Edge>(),
                        1,
                        fp,
                    );
                    let e = node(&mut eh, old);
                    read_record(e, fp, "error reading Edge\n");
                    (*e).eu = node(&mut euh, (*e).eu);
                    read_attributes(
                        &mut (*e).ea,
                        (*m).sea,
                        fp,
                        "error reading edge attributes\n",
                    );
                }
                4 => {
                    let mut old = std::ptr::null_mut::<VertexUse>();
                    fread(
                        (&raw mut old).cast(),
                        std::mem::size_of::<*mut VertexUse>(),
                        1,
                        fp,
                    );
                    let vu = node(&mut vuh, old);
                    read_record(vu, fp, "error reading VertexUse\n");
                    (*vu).eu = node(&mut euh, (*vu).eu);
                    (*vu).vuPrev = node(&mut vuh, (*vu).vuPrev);
                    (*vu).vuNext = node(&mut vuh, (*vu).vuNext);
                    (*vu).v = node(&mut vh, (*vu).v);
                    read_attributes(
                        &mut (*vu).vua,
                        (*m).svua,
                        fp,
                        "error reading vertex-use attributes\n",
                    );
                }
                5 => {
                    let mut old = std::ptr::null_mut::<Vertex>();
                    fread(
                        (&raw mut old).cast(),
                        std::mem::size_of::<*mut Vertex>(),
                        1,
                        fp,
                    );
                    let v = node(&mut vh, old);
                    read_record(v, fp, "error reading Vertex");
                    (*v).vu = node(&mut vuh, (*v).vu);
                    read_attributes(
                        &mut (*v).va,
                        (*m).sva,
                        fp,
                        "error reading vertex attributes\n",
                    );
                }
                _ => fail("invalid type in model\n"),
            }
        }
        for h in [fh, euh, vuh, eh, vh] {
            release_tree(h);
        }
        m
    }
}

unsafe fn write_record<T>(fp: *mut c_void, kind: i32, p: *mut T, attr: *mut c_void, size: usize) {
    unsafe {
        fwrite((&raw const kind).cast(), std::mem::size_of::<i32>(), 1, fp);
        fwrite((&raw const p).cast(), std::mem::size_of::<*mut T>(), 1, fp);
        fwrite(p.cast(), std::mem::size_of::<T>(), 1, fp);
        if !attr.is_null() && size != 0 {
            fwrite(attr, size, 1, fp);
        }
    }
}

/// C: void writeModel(Model *m, FILE *fp);
/// Writes the original native binary format, including pointer IDs and attributes.
/// # Safety
/// All model links and non-null attribute buffers must be valid, including
/// structure padding bytes. fp must be a writable C FILE stream.
pub unsafe fn writeModel(m: *mut Model, fp: *mut c_void) {
    unsafe {
        let mut vh = std::ptr::null_mut();
        let mut eh = std::ptr::null_mut();
        let mut repeat = 0;
        write_record(fp, 0, m, (*m).ma, (*m).sma);
        let mut f = (*m).f;
        let mut eu;
        loop {
            write_record(fp, 1, f, (*f).fa, (*m).sfa);
            eu = (*f).eu;
            loop {
                write_record(fp, 2, eu, (*eu).eua, (*m).seua);
                let e = (*eu).e;
                eh = updateAddressTree(eh, e.cast(), &mut repeat);
                if repeat == 0 {
                    write_record(fp, 3, e, (*e).ea, (*m).sea);
                }
                let vu = (*eu).vu;
                write_record(fp, 4, vu, (*vu).vua, (*m).svua);
                let v = (*vu).v;
                vh = updateAddressTree(vh, v.cast(), &mut repeat);
                if repeat == 0 {
                    write_record(fp, 5, v, (*v).va, (*m).sva);
                }
                eu = (*eu).euCW;
                if eu == (*f).eu {
                    break;
                }
            }
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
        f = (*m).f;
        loop {
            if (*(*(*f).eu).euMate).f.is_null() {
                eu = (*(*f).eu).euMate;
                break;
            } else if (*(*(*(*f).eu).euCW).euMate).f.is_null() {
                eu = (*(*(*f).eu).euCW).euMate;
                break;
            } else if (*(*(*(*f).eu).euCCW).euMate).f.is_null() {
                eu = (*(*(*f).eu).euCCW).euMate;
                break;
            }
            f = (*f).fNext;
            if f == (*m).f {
                break;
            }
        }
        if !(*eu).f.is_null() {
            fail("writeModel: can't find exterior face\n");
        }
        let first = eu;
        loop {
            write_record(fp, 2, eu, (*eu).eua, (*m).seua);
            let vu = (*eu).vu;
            write_record(fp, 4, vu, (*vu).vua, (*m).svua);
            eu = (*eu).euCW;
            if eu == first {
                break;
            }
        }
        release_tree(vh);
        release_tree(eh);
    }
}
