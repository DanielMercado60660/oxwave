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

//! Differential tests for tri/lib/{check,circum,colinear,create,delete,
//! distance,fixedges,inside,nearest,project,readwrite}.c.
//! Credits: Dave Hale, Chris Elmer, Jack K. Cohen, and Craig Artley,
//! Colorado School of Mines / Center for Wave Phenomena, 1990–1994.

use oxwave_cref::XorShift;
use oxwave_lib_tri::{self as tri, Vertex};
use std::ffi::c_void;
use std::ptr::null_mut;

#[test]
fn unfix_edge_swaps_restore_c_delaunay_mesh() {
    let mut rng = XorShift::new(47);
    let mut swaps = 0;
    for i in 0..2048 {
        let p = unsafe { Pair::bounds(-16.0, -16.0, 16.0, 16.0) };
        let (x, y) = p.point(&mut rng, 8);
        p.add(x, y);
        unsafe {
            let triangles = |m| {
                let mut out = Vec::new();
                for (k, p) in graph(m) {
                    if k == 1 {
                        let f = p.cast::<tri::Face>();
                        let mut eu = (*f).eu;
                        let mut vs = Vec::new();
                        loop {
                            let v = &*(*(*eu).vu).v;
                            vs.push((v.x.to_bits(), v.y.to_bits()));
                            eu = (*eu).euCW;
                            if eu == (*f).eu {
                                break;
                            }
                        }
                        vs.sort();
                        out.push(vs);
                    }
                }
                out.sort();
                out
            };
            // The C pre-constraint model supplies the Delaunay oracle. C's
            // unfixEdge would reread freed edge uses when this constraint
            // flips, so its undefined read is not used as a numeric oracle.
            let expected = triangles(p.c);
            let (a, b) = constraint_pair(p.r, 1, i, 0);
            tri::fixEdgeBetweenVertices(a, b);
            let constrained = triangles(p.r);
            if constrained != expected {
                swaps += 1;
            }
            loop {
                let e = graph(p.r)
                    .into_iter()
                    .find(|(k, p)| {
                        *k == 3 && {
                            let e = &*p.cast::<tri::Edge>();
                            e.fixed != 0 && !(*e.eu).f.is_null() && !(*(*e.eu).euMate).f.is_null()
                        }
                    })
                    .map(|(_, p)| p.cast::<tri::Edge>());
                let Some(e) = e else {
                    break;
                };
                assert_eq!(tri::unfixEdge(e), 1);
            }
            assert_eq!(triangles(p.r), expected);
            tri::checkModel(p.r);
        }
    }
    assert!(swaps > 500, "edge correction must exercise actual flips");
}

#[test]
fn geometry_edge_ranges_match_c() {
    let mut rng = XorShift::new(45);
    let mut f = Fixture::new();
    let mut g = Fixture::new();
    for i in 0..12000 {
        let scale = [1.0e-12f32, 1.0e-6, 1.0, 1.0e6, 1.0e12][i % 5];
        for v in &mut *f.vertices {
            v.x = value(&mut rng) * scale;
            v.y = value(&mut rng) * scale;
        }
        for v in &mut *g.vertices {
            v.x = value(&mut rng) * scale;
            v.y = value(&mut rng) * scale;
        }
        if i % 7 == 0 {
            f.vertices[2].x = f.vertices[0].x;
            f.vertices[2].y = f.vertices[0].y;
        }
        if i % 11 == 0 {
            f.vertices[2].x = f.vertices[0].x + (f.vertices[1].x - f.vertices[0].x) * 0.5;
            f.vertices[2].y = f.vertices[0].y + (f.vertices[1].y - f.vertices[0].y) * 0.5;
        }
        if i % 13 == 0 {
            for v in &mut *f.vertices {
                v.x = 0.0;
                v.y = -0.0;
            }
        }
        let (a, b, c) = (&f.vertices[0], &f.vertices[1], &f.vertices[2]);
        let (mut x, mut y, mut rs) = (0.0, 0.0, 0.0);
        let (mut cx, mut cy, mut crs) = (0.0, 0.0, 0.0);
        tri::circum(a.x, a.y, b.x, b.y, c.x, c.y, &mut x, &mut y, &mut rs);
        unsafe {
            circum(a.x, a.y, b.x, b.y, c.x, c.y, &mut cx, &mut cy, &mut crs);
        }
        eq_float(x, cx);
        eq_float(y, cy);
        eq_float(rs, crs);
        let v = vertex(a.x + (b.x - a.x) * 0.5, a.y + (b.y - a.y) * 0.5);
        assert_eq!(tri::vertexBetweenVertices(&v, a, b), unsafe {
            vertexBetweenVertices(&v, a, b)
        });
        let (qx, qy) = if i % 3 == 0 { (a.x, a.y) } else { (v.x, v.y) };
        assert_eq!(tri::in3Vertices(qx, qy, a, b, c), unsafe {
            in3Vertices(qx, qy, a, b, c)
        });
        assert_eq!(tri::inCircum(qx, qy, x, y, rs), unsafe {
            inCircum(qx, qy, cx, cy, crs)
        });
        unsafe {
            assert_eq!(
                tri::edgesColinear(&f.edge, &g.edge),
                edgesColinear(&*f.edge, &*g.edge)
            );
            assert_eq!(
                tri::edgesColinear(&f.edge, &f.edge),
                edgesColinear(&*f.edge, &*f.edge)
            );
            eq_float(
                tri::distanceToEdge(&f.edge, qx, qy),
                distanceToEdge(&*f.edge, qx, qy),
            );
            let (mut x, mut y, mut cx, mut cy) = (qx, qy, qx, qy);
            tri::projectToEdge(&f.edge, &mut x, &mut y);
            projectToEdge(&*f.edge, &mut cx, &mut cy);
            eq_float(x, cx);
            eq_float(y, cy);
        }
    }
}

#[test]
fn read_model_errors_match_c() {
    if let Ok(s) = std::env::var("OXWAVE_TRI_CHILD") {
        let args: Vec<_> = s.split(':').collect();
        let rust = args[0] == "rust";
        let case: usize = args[1].parse().unwrap();
        let stream = Stream::new();
        let mut rng = XorShift::new(46);
        let p = Pair::new(&mut rng, 0);
        unsafe {
            if case == 0 {
                let kind = -1i32;
                fwrite(
                    (&raw const kind).cast(),
                    std::mem::size_of::<i32>(),
                    1,
                    stream.fp,
                );
            } else {
                let full = Stream::new();
                writeModel(p.c, full.fp);
                let bytes = full.bytes();
                // Truncate inside each kind of structure, after complete
                // preceding records, so every fread error is compared.
                let mut starts = [None; 6];
                let mut offset = 0usize;
                while offset < bytes.len() {
                    let k =
                        i32::from_ne_bytes(bytes[offset..offset + 4].try_into().unwrap()) as usize;
                    let size = [
                        std::mem::size_of::<tri::Model>(),
                        std::mem::size_of::<tri::Face>(),
                        std::mem::size_of::<tri::EdgeUse>(),
                        std::mem::size_of::<tri::Edge>(),
                        std::mem::size_of::<tri::VertexUse>(),
                        std::mem::size_of::<Vertex>(),
                    ][k];
                    starts[k].get_or_insert(offset + 4 + std::mem::size_of::<usize>() + size - 1);
                    offset += 4 + std::mem::size_of::<usize>() + size;
                }
                let n = starts[case - 1].unwrap();
                fwrite(bytes.as_ptr().cast(), 1, n, stream.fp);
            }
            rewind(stream.fp);
            if rust {
                tri::readModel(stream.fp);
            } else {
                readModel(stream.fp);
            }
        }
        panic!("invalid model stream was accepted");
    }
    for case in 0..7 {
        let r = child("read_model_errors_match_c", &format!("rust:{case}"));
        let c = child("read_model_errors_match_c", &format!("c:{case}"));
        assert_eq!(r.status.code(), Some(1));
        assert_eq!(r.status.code(), c.status.code());
        assert_eq!(r.stderr, c.stderr);
    }
}

// The SU constraint walker can fail to terminate for arbitrary lines passing
// almost through a third vertex (its collinearity tolerance is about 1%).
// Exercise existing chains in larger random meshes, and edge replacement in
// random five-vertex meshes, with no discarded inputs or timeout-based skips.
fn constraint_pair(
    m: *mut tri::Model,
    i: usize,
    a: usize,
    _b: usize,
) -> (*mut Vertex, *mut Vertex) {
    unsafe {
        if i % 2 == 0 {
            let e = choose::<tri::Edge>(m, 3, a);
            let eu = (*e).eu;
            ((*(*eu).vu).v, (*(*(*eu).euCW).vu).v)
        } else {
            let find = |x, y| {
                graph(m)
                    .into_iter()
                    .filter(|(k, _)| *k == 5)
                    .map(|(_, p)| p.cast::<Vertex>())
                    .find(|v| (**v).x == x && (**v).y == y)
                    .unwrap()
            };
            (find((*m).xmin, (*m).ymax), find((*m).xmax, (*m).ymin))
        }
    }
}

thread_local! {
    static EVENTS:std::cell::RefCell<Vec<Vec<u32>>>=const {std::cell::RefCell::new(Vec::new())};
}
unsafe extern "C" fn on_v_add(_: *mut tri::Model, v: *mut Vertex) {
    EVENTS.with(|e| unsafe {
        e.borrow_mut()
            .push(vec![0, (*v).x.to_bits(), (*v).y.to_bits()]);
    });
}
unsafe extern "C" fn on_v_del(_: *mut tri::Model, v: *mut Vertex) {
    EVENTS.with(|e| unsafe {
        e.borrow_mut()
            .push(vec![1, (*v).x.to_bits(), (*v).y.to_bits()]);
    });
}
unsafe fn on_tri(kind: u32, t: *mut tri::Tri) {
    unsafe {
        let mut event = vec![kind];
        let mut eu = (*t).eu;
        loop {
            let v = (*(*eu).vu).v;
            event.extend([(*v).x.to_bits(), (*v).y.to_bits()]);
            eu = (*eu).euCW;
            if eu == (*t).eu {
                break;
            }
        }
        EVENTS.with(|e| e.borrow_mut().push(event));
    }
}
unsafe extern "C" fn on_t_add(_: *mut tri::Model, t: *mut tri::Tri) {
    unsafe {
        on_tri(2, t);
    }
}
unsafe extern "C" fn on_t_del(_: *mut tri::Model, t: *mut tri::Tri) {
    unsafe {
        on_tri(3, t);
    }
}
fn events(f: impl FnOnce()) -> Vec<Vec<u32>> {
    EVENTS.with(|e| e.borrow_mut().clear());
    f();
    EVENTS.with(|e| std::mem::take(&mut *e.borrow_mut()))
}
#[test]
fn topology_callbacks_match_c() {
    let mut rng = XorShift::new(32);
    for i in 0..512 {
        let p = unsafe { Pair::bounds(-16.0, -16.0, 16.0, 16.0) };
        unsafe {
            for m in [p.r, p.c] {
                (*m).vAdd = Some(on_v_add);
                (*m).vDel = Some(on_v_del);
                (*m).tAdd = Some(on_t_add);
                (*m).tDel = Some(on_t_del);
            }
            for _ in 0..8 {
                let (x, y) = p.point(&mut rng, 8);
                assert_eq!(
                    events(|| {
                        tri::addVertexToModel(p.r, x, y);
                    }),
                    events(|| {
                        addVertexToModel(p.c, x, y);
                    })
                );
                p.compare(true);
            }
            let (rv, cv) = (choose(p.r, 5, i), choose(p.c, 5, i));
            assert_eq!(
                events(|| {
                    tri::deleteVertexFromModel(p.r, rv);
                }),
                events(|| {
                    deleteVertexFromModel(p.c, cv);
                })
            );
            p.compare(true);
            // Existing Delaunay edges do not trigger the original C unfixEdge
            // use-after-free; flips are exercised by insertion and deletion.
            let (ra, rb) = constraint_pair(p.r, 2 * i, i, i + 1);
            let (ca, cb) = constraint_pair(p.c, 2 * i, i, i + 1);
            assert_eq!(
                events(|| {
                    tri::fixEdgeBetweenVertices(ra, rb);
                }),
                events(|| {
                    fixEdgeBetweenVertices(ca, cb);
                })
            );
            p.compare(true);
            let (re, ce) = (choose(p.r, 3, i), choose(p.c, 3, i));
            assert_eq!(
                events(|| {
                    tri::unfixEdge(re);
                }),
                events(|| {
                    unfixEdge(ce);
                })
            );
            p.compare(true);
        }
    }
}

fn query_suite(kind: usize) {
    let mut rng = XorShift::new(40 + kind as u64);
    for i in 0..72 {
        let p = Pair::new(&mut rng, i);
        p.cloud(&mut rng, sizes(i));
        for j in 0..256 {
            unsafe {
                let m = p.r;
                let (x, y) = if j % 5 == 0 {
                    p.point(&mut rng, j)
                } else {
                    (
                        (*m).xmin + ((rng.next_f64() * 2.0 - 0.5) as f32) * ((*m).xmax - (*m).xmin),
                        (*m).ymin + ((rng.next_f64() * 2.0 - 0.5) as f32) * ((*m).ymax - (*m).ymin),
                    )
                };
                match kind {
                    0 => {
                        let start = if j % 2 == 0 {
                            null_mut()
                        } else {
                            choose(m, 1, j)
                        };
                        assert_eq!(
                            tri::insideTriInModel(m, start, x, y),
                            insideTriInModel(m, start, x, y)
                        );
                    }
                    1 => {
                        let start = if j % 2 == 0 {
                            null_mut()
                        } else {
                            choose(m, 3, j)
                        };
                        assert_eq!(
                            tri::nearestEdgeInModel(m, start, x, y),
                            nearestEdgeInModel(m, start, x, y)
                        );
                    }
                    2 => {
                        let start = if j % 2 == 0 {
                            null_mut()
                        } else {
                            choose(m, 5, j)
                        };
                        assert_eq!(
                            tri::nearestVertexInModel(m, start, x, y),
                            nearestVertexInModel(m, start, x, y)
                        );
                    }
                    _ => unreachable!(),
                }
            }
        }
    }
}

unsafe fn run_check(kind: usize, rust: bool, m: *mut tri::Model, index: usize) {
    unsafe {
        match (kind, rust) {
            (0, true) => tri::checkVertexUse(choose(m, 4, index)),
            (0, false) => checkVertexUse(choose(m, 4, index)),
            (1, true) => tri::checkEdgeUse(choose(m, 2, index)),
            (1, false) => checkEdgeUse(choose(m, 2, index)),
            (2, true) => tri::checkFace(choose(m, 1, index)),
            (2, false) => checkFace(choose(m, 1, index)),
            (3, true) => tri::checkModel(m),
            (3, false) => checkModel(m),
            _ => unreachable!(),
        }
    }
}
fn child(test: &str, case: &str) -> std::process::Output {
    std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", test, "--nocapture"])
        .env("OXWAVE_TRI_CHILD", case)
        // Cargo's default target directory does not require this variable.
        // Exercise that environment in every subprocess, including I/O errors.
        .env_remove("CARGO_TARGET_DIR")
        .output()
        .unwrap()
}
fn check_cases(kind: usize, test: &str) {
    if let Ok(s) = std::env::var("OXWAVE_TRI_CHILD") {
        let args: Vec<_> = s.split(':').collect();
        let rust = args[0] == "rust";
        let case: usize = args[1].parse().unwrap();
        let mut rng = XorShift::new(24);
        let p = Pair::new(&mut rng, 0);
        unsafe {
            let m = if rust { p.r } else { p.c };
            let eu = (*(*m).f).eu;
            let vu = (*eu).vu;
            match (kind, case) {
                (0, 0) => (*(*vu).vuNext).v = (*(*(*eu).euCW).vu).v,
                (0, 1) => (*(*vu).vuNext).vuPrev = (*vu).vuNext,
                (0, 2) => (*(*vu).eu).vu = (*(*vu).eu).euCW.as_ref().unwrap().vu,
                (1, 0) => (*(*eu).euMate).euMate = (*eu).euMate,
                (1, 1) => (*(*eu).euCW).euCCW = (*eu).euCW,
                (1, 2) => (*eu).euCCW = (*eu).euCW,
                (2, 0) => (*(*(*m).f).fNext).fPrev = (*(*m).f).fNext,
                (2, 1) => (*(*m).f).eu = (*eu).euMate,
                (2, 2) => (*(*eu).euMate).euMate = (*eu).euMate,
                (3, 0) => (*(*m).f).m = null_mut(),
                (3, 1) => (*(*(*m).f).fNext).m = null_mut(),
                (3, 2) => (*(*eu).euMate).euMate = (*eu).euMate,
                _ => unreachable!(),
            }
            match (kind, rust) {
                (0, true) => tri::checkVertexUse(vu),
                (0, false) => checkVertexUse(vu),
                (1, true) => tri::checkEdgeUse(eu),
                (1, false) => checkEdgeUse(eu),
                (2, true) => tri::checkFace((*m).f),
                (2, false) => checkFace((*m).f),
                (3, true) => tri::checkModel(m),
                (3, false) => checkModel(m),
                _ => unreachable!(),
            }
        }
        panic!("invalid model was accepted");
    }
    let mut rng = XorShift::new(24 + kind as u64);
    for i in 0..2048 {
        let p = Pair::new(&mut rng, i);
        p.cloud(&mut rng, i % 7);
        let before = unsafe { snapshot(p.r, true) };
        unsafe {
            run_check(kind, true, p.r, i);
            run_check(kind, false, p.c, i);
        }
        p.compare(true);
        assert_eq!(before, unsafe { snapshot(p.r, true) });
    }
    for case in 0..3 {
        let r = child(test, &format!("rust:{case}"));
        let c = child(test, &format!("c:{case}"));
        assert_eq!(r.status.code(), Some(255));
        assert_eq!(r.status.code(), c.status.code());
        assert_eq!(r.stderr, c.stderr);
    }
}
#[test]
fn check_vertex_use_matches_c() {
    check_cases(0, "check_vertex_use_matches_c");
}
#[test]
fn check_edge_use_matches_c() {
    check_cases(1, "check_edge_use_matches_c");
}
#[test]
fn check_face_matches_c() {
    check_cases(2, "check_face_matches_c");
}
#[test]
fn check_model_matches_c() {
    check_cases(3, "check_model_matches_c");
}

#[test]
fn bad_model_matches_c() {
    if let Ok(s) = std::env::var("OXWAVE_TRI_CHILD") {
        unsafe {
            if s == "rust" {
                tri::badModel();
            } else {
                badModel();
            }
        }
        panic!("badModel returned");
    }
    // No arguments or random domain: compare the only possible behavior.
    let r = child("bad_model_matches_c", "rust");
    let c = child("bad_model_matches_c", "c");
    assert_eq!(r.status.code(), Some(255));
    assert_eq!(r.status.code(), c.status.code());
    assert_eq!(r.stderr, c.stderr);
}

#[test]
fn kill_model_matches_c() {
    let mut rng = XorShift::new(28);
    for i in 0..2048 {
        let p = std::mem::ManuallyDrop::new(Pair::new(&mut rng, i));
        p.cloud(&mut rng, i % 7);
        p.compare(true);
        unsafe {
            // Both implementations preserve SU's leak of the exterior uses,
            // boundary edges and boundary vertices. Reclaim those in the test.
            let exterior = |m| {
                let mut set = std::collections::HashSet::new();
                for (k, q) in graph(m) {
                    if k == 2 {
                        let eu = q.cast::<tri::EdgeUse>();
                        if (*eu).f.is_null() {
                            set.insert((2, eu as usize));
                            set.insert((4, (*eu).vu as usize));
                            set.insert((3, (*eu).e as usize));
                            set.insert((5, (*(*eu).vu).v as usize));
                        }
                    }
                }
                set
            };
            let (r, c) = (exterior(p.r), exterior(p.c));
            tri::killModel(p.r);
            killModel(p.c);
            // killModel nulls mates of surviving exterior uses.
            for nodes in [r, c] {
                for &(k, q) in &nodes {
                    if k == 2 {
                        assert!((*(q as *mut tri::EdgeUse)).euMate.is_null());
                    }
                }
                for (_, q) in nodes {
                    free(q as *mut c_void);
                }
            }
        }
    }
}
unsafe fn tree(p: *mut tri::Address) -> Vec<Option<usize>> {
    unsafe {
        if p.is_null() {
            vec![None]
        } else {
            let mut v = vec![Some((*p).oaddress as usize)];
            v.extend(tree((*p).aLeft));
            v.extend(tree((*p).aRight));
            v
        }
    }
}
unsafe fn free_tree(p: *mut tri::Address) {
    unsafe {
        if !p.is_null() {
            free_tree((*p).aLeft);
            free_tree((*p).aRight);
            free(p.cast());
        }
    }
}
#[test]
fn update_address_tree_matches_c() {
    let mut rng = XorShift::new(29);
    for i in 0..128 {
        let (mut r, mut c) = (null_mut(), null_mut());
        let mut keys = Vec::new();
        unsafe {
            for j in 0..sizes(i).max(32) {
                let key = match j % 11 {
                    0 => 0,
                    1 => usize::MAX,
                    2 => i64::MAX as usize,
                    3 if !keys.is_empty() => keys[j % keys.len()],
                    _ => rng.next_u64() as usize,
                };
                keys.push(key);
                let (mut rf, mut cf) = (-1, -1);
                r = tri::updateAddressTree(r, key as *mut c_void, &mut rf);
                c = updateAddressTree(c, key as *mut c_void, &mut cf);
                assert_eq!(rf, cf);
                assert_eq!(tree(r), tree(c));
            }
            free_tree(r);
            free_tree(c);
        }
    }
}

static STREAM_ID: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
struct Stream {
    fp: *mut c_void,
}
impl Stream {
    fn new() -> Self {
        use std::os::fd::{AsRawFd, IntoRawFd};
        let id = STREAM_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        // Keep files in the build directory without assuming CARGO_TARGET_DIR
        // is exported (it is optional even when Cargo launches this test).
        let executable = std::env::current_exe().unwrap();
        let path = executable
            .parent()
            .unwrap()
            .join(format!("tri-diff-{}-{id}", std::process::id()));
        let file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        // An unlinked Unix file stays open until fclose or process exit. The
        // malformed-input subprocesses exit inside C/Rust without running Drop.
        std::fs::remove_file(&path).unwrap();
        let fp = unsafe { fdopen(file.as_raw_fd(), b"w+b\0".as_ptr().cast()) };
        assert!(!fp.is_null());
        // Transfer ownership only after fdopen succeeds.
        let _ = file.into_raw_fd();
        Self { fp }
    }
    fn bytes(&self) -> Vec<u8> {
        unsafe {
            rewind(self.fp);
            let mut out = Vec::new();
            let mut buf = [0u8; 4096];
            loop {
                let n = fread(buf.as_mut_ptr().cast(), 1, buf.len(), self.fp);
                if n == 0 {
                    break;
                }
                out.extend_from_slice(&buf[..n]);
            }
            out
        }
    }
}
impl Drop for Stream {
    fn drop(&mut self) {
        unsafe {
            assert_eq!(fclose(self.fp), 0);
        }
    }
}
extern "C" {
    fn fdopen(fd: i32, mode: *const std::ffi::c_char) -> *mut c_void;
}

unsafe fn attributes(m: *mut tri::Model, rng: &mut XorShift, i: usize) {
    unsafe {
        let ss = [0, 1, 2, 3, 7, 16, 31, 64, 257];
        (*m).sma = ss[i % 9];
        (*m).sfa = ss[(i + 1) % 9];
        (*m).seua = ss[(i + 2) % 9];
        (*m).sea = ss[(i + 3) % 9];
        (*m).svua = ss[(i + 4) % 9];
        (*m).sva = ss[(i + 5) % 9];
        for (k, p) in graph(m) {
            let (slot, n) = match k {
                0 => (&mut (*p.cast::<tri::Model>()).ma, (*m).sma),
                1 => (&mut (*p.cast::<tri::Face>()).fa, (*m).sfa),
                2 => (&mut (*p.cast::<tri::EdgeUse>()).eua, (*m).seua),
                3 => (&mut (*p.cast::<tri::Edge>()).ea, (*m).sea),
                4 => (&mut (*p.cast::<tri::VertexUse>()).vua, (*m).svua),
                5 => (&mut (*p.cast::<Vertex>()).va, (*m).sva),
                _ => unreachable!(),
            };
            if n != 0 {
                *slot = calloc(1, n);
                for b in std::slice::from_raw_parts_mut((*slot).cast::<u8>(), n) {
                    *b = rng.next_u32() as u8;
                }
            } else if i % 3 == 0 {
                // A non-null attribute pointer with zero size is not written
                // and is reset to NULL by both readers.
                *slot = calloc(1, 1);
            }
        }
    }
}
#[test]
fn write_model_matches_c() {
    let mut rng = XorShift::new(30);
    for i in 0..2048 {
        let p = Pair::new(&mut rng, i);
        p.cloud(&mut rng, i % 7);
        unsafe {
            attributes(p.r, &mut rng, i);
            let (r, c) = (Stream::new(), Stream::new());
            tri::writeModel(p.r, r.fp);
            writeModel(p.r, c.fp);
            assert_eq!(r.bytes(), c.bytes(), "write case {i}");
            rewind(r.fp);
            let decoded = readModel(r.fp);
            assert_eq!(snapshot(p.r, true), snapshot(decoded, true));
            destroy(decoded);
        }
    }
}
#[test]
fn read_model_matches_c() {
    let mut rng = XorShift::new(31);
    for i in 0..2048 {
        let p = Pair::new(&mut rng, i);
        p.cloud(&mut rng, i % 7);
        unsafe {
            attributes(p.c, &mut rng, i);
            let stream = Stream::new();
            (*p.c).vAdd = Some(on_v_add);
            (*p.c).vDel = Some(on_v_del);
            (*p.c).tAdd = Some(on_t_add);
            (*p.c).tDel = Some(on_t_del);
            writeModel(p.c, stream.fp);
            // Readers must clear callback addresses found in the binary data.
            (*p.c).vAdd = None;
            (*p.c).vDel = None;
            (*p.c).tAdd = None;
            (*p.c).tDel = None;
            rewind(stream.fp);
            let r = tri::readModel(stream.fp);
            rewind(stream.fp);
            let c = readModel(stream.fp);
            assert_eq!(snapshot(r, true), snapshot(c, true));
            assert_eq!(snapshot(r, true), snapshot(p.c, true));
            tri::checkModel(r);
            checkModel(c);
            destroy(r);
            destroy(c);
        }
    }
    let stream = Stream::new();
    unsafe {
        assert!(tri::readModel(stream.fp).is_null());
        rewind(stream.fp);
        assert!(readModel(stream.fp).is_null());
    }
}

// Enumerate the complete graph, including the exterior ring. IDs are assigned
// by link traversal, independently of allocator addresses. Comparing every row
// checks aliasing, representative pointers, list order, geometry and attributes.
type Node = (u8, *mut c_void);
unsafe fn graph(m: *mut tri::Model) -> Vec<Node> {
    let mut out: Vec<Node> = vec![(0, m.cast())];
    let mut seen = std::collections::HashMap::new();
    seen.insert((0, m as usize), 0usize);
    let mut i = 0;
    while i < out.len() {
        let (kind, p) = out[i];
        let links: Vec<Node> = unsafe {
            match kind {
                0 => vec![(1, (*p.cast::<tri::Model>()).f.cast())],
                1 => {
                    let f = &*p.cast::<tri::Face>();
                    vec![
                        (0, f.m.cast()),
                        (1, f.fPrev.cast()),
                        (1, f.fNext.cast()),
                        (2, f.eu.cast()),
                    ]
                }
                2 => {
                    let eu = &*p.cast::<tri::EdgeUse>();
                    vec![
                        (1, eu.f.cast()),
                        (4, eu.vu.cast()),
                        (2, eu.euMate.cast()),
                        (2, eu.euCW.cast()),
                        (2, eu.euCCW.cast()),
                        (3, eu.e.cast()),
                    ]
                }
                3 => vec![(2, (*p.cast::<tri::Edge>()).eu.cast())],
                4 => {
                    let vu = &*p.cast::<tri::VertexUse>();
                    vec![
                        (2, vu.eu.cast()),
                        (4, vu.vuPrev.cast()),
                        (4, vu.vuNext.cast()),
                        (5, vu.v.cast()),
                    ]
                }
                5 => vec![(4, (*p.cast::<Vertex>()).vu.cast())],
                _ => unreachable!(),
            }
        };
        for (k, q) in links {
            if !q.is_null() && !seen.contains_key(&(k, q as usize)) {
                seen.insert((k, q as usize), out.len());
                out.push((k, q));
            }
        }
        i += 1;
        assert!(i < 100_000, "invalid or unexpectedly large model");
    }
    out
}
unsafe fn snapshot(m: *mut tri::Model, circles: bool) -> Vec<Vec<u64>> {
    unsafe {
        let nodes = graph(m);
        let ids: std::collections::HashMap<_, _> = nodes
            .iter()
            .enumerate()
            .map(|(i, (k, p))| ((*k, *p as usize), i as u64 + 1))
            .collect();
        let id = |k, p: *mut c_void| {
            if p.is_null() {
                0
            } else {
                ids[&(k, p as usize)]
            }
        };
        let attr = |p: *mut c_void, n| {
            if p.is_null() || n == 0 {
                Vec::new()
            } else {
                std::slice::from_raw_parts(p.cast::<u8>(), n)
                    .iter()
                    .map(|x| u64::from(*x))
                    .collect()
            }
        };
        nodes
            .iter()
            .map(|&(k, p)| {
                let (mut row, a, n) = match k {
                    0 => {
                        let a = &*p.cast::<tri::Model>();
                        (
                            vec![
                                id(1, a.f.cast()),
                                a.xmin.to_bits() as u64,
                                a.ymin.to_bits() as u64,
                                a.xmax.to_bits() as u64,
                                a.ymax.to_bits() as u64,
                                a.eps.to_bits() as u64,
                                a.sma as u64,
                                a.sfa as u64,
                                a.seua as u64,
                                a.sea as u64,
                                a.svua as u64,
                                a.sva as u64,
                                a.vAdd.is_some() as u64,
                                a.vDel.is_some() as u64,
                                a.tAdd.is_some() as u64,
                                a.tDel.is_some() as u64,
                            ],
                            a.ma,
                            a.sma,
                        )
                    }
                    1 => {
                        let a = &*p.cast::<tri::Face>();
                        let mut row = vec![
                            id(0, a.m.cast()),
                            id(1, a.fPrev.cast()),
                            id(1, a.fNext.cast()),
                            id(2, a.eu.cast()),
                        ];
                        if circles {
                            row.extend([
                                a.xc.to_bits() as u64,
                                a.yc.to_bits() as u64,
                                a.rs.to_bits() as u64,
                            ]);
                        }
                        (row, a.fa, (*m).sfa)
                    }
                    2 => {
                        let a = &*p.cast::<tri::EdgeUse>();
                        (
                            vec![
                                id(1, a.f.cast()),
                                id(4, a.vu.cast()),
                                id(2, a.euMate.cast()),
                                id(2, a.euCW.cast()),
                                id(2, a.euCCW.cast()),
                                id(3, a.e.cast()),
                            ],
                            a.eua,
                            (*m).seua,
                        )
                    }
                    3 => {
                        let a = &*p.cast::<tri::Edge>();
                        (vec![id(2, a.eu.cast()), a.fixed as u64], a.ea, (*m).sea)
                    }
                    4 => {
                        let a = &*p.cast::<tri::VertexUse>();
                        (
                            vec![
                                id(2, a.eu.cast()),
                                id(4, a.vuPrev.cast()),
                                id(4, a.vuNext.cast()),
                                id(5, a.v.cast()),
                            ],
                            a.vua,
                            (*m).svua,
                        )
                    }
                    5 => {
                        let a = &*p.cast::<Vertex>();
                        (
                            vec![
                                id(4, a.vu.cast()),
                                a.fixed as u64,
                                a.x.to_bits() as u64,
                                a.y.to_bits() as u64,
                            ],
                            a.va,
                            (*m).sva,
                        )
                    }
                    _ => unreachable!(),
                };
                row.insert(0, k as u64);
                row.push((!a.is_null() && n != 0) as u64);
                row.extend(attr(a, n));
                row
            })
            .collect()
    }
}
unsafe fn destroy(m: *mut tri::Model) {
    unsafe {
        let nodes = graph(m);
        for &(k, p) in &nodes {
            let a = match k {
                0 => (*p.cast::<tri::Model>()).ma,
                1 => (*p.cast::<tri::Face>()).fa,
                2 => (*p.cast::<tri::EdgeUse>()).eua,
                3 => (*p.cast::<tri::Edge>()).ea,
                4 => (*p.cast::<tri::VertexUse>()).vua,
                5 => (*p.cast::<Vertex>()).va,
                _ => unreachable!(),
            };
            if !a.is_null() {
                free(a);
            }
        }
        for (_, p) in nodes {
            free(p);
        }
    }
}
struct Pair {
    r: *mut tri::Model,
    c: *mut tri::Model,
}
impl Pair {
    fn new(rng: &mut XorShift, i: usize) -> Self {
        let scale = [0.001, 1.0, 1000.0][i % 3];
        let x = value(rng) * scale;
        let y = value(rng) * scale;
        let w = (value(rng).abs() + 10.0) * scale;
        let h = (value(rng).abs() + 10.0) * scale;
        unsafe { Self::bounds(x, y, x + w, y + h) }
    }
    unsafe fn bounds(x: f32, y: f32, xx: f32, yy: f32) -> Self {
        unsafe {
            let p = Self {
                r: tri::makeModel(x, y, xx, yy),
                c: makeModel(x, y, xx, yy),
            };
            // C makeModel leaves ma uninitialized. No attributes are attached.
            (*p.c).ma = null_mut();
            p
        }
    }
    fn compare(&self, circles: bool) {
        unsafe {
            assert_eq!(snapshot(self.r, circles), snapshot(self.c, circles));
        }
    }
    fn point(&self, rng: &mut XorShift, mode: usize) -> (f32, f32) {
        unsafe {
            let m = &*self.r;
            let mut u = 0.03 + 0.94 * rng.next_f64() as f32;
            let mut v = 0.03 + 0.94 * rng.next_f64() as f32;
            match mode % 13 {
                0 => u = 0.0,
                1 => v = 0.0,
                2 => u = 1.0,
                3 => v = 1.0,
                4 => {
                    u = 0.0;
                    v = 0.0;
                }
                _ => (),
            }
            (
                m.xmin + u * (m.xmax - m.xmin),
                m.ymin + v * (m.ymax - m.ymin),
            )
        }
    }
    fn add(&self, x: f32, y: f32) -> (*mut Vertex, *mut Vertex) {
        unsafe {
            let r = tri::addVertexToModel(self.r, x, y);
            let c = addVertexToModel(self.c, x, y);
            assert_eq!(r.is_null(), c.is_null());
            if !r.is_null() {
                eq_float((*r).x, (*c).x);
                eq_float((*r).y, (*c).y);
            }
            self.compare(true);
            (r, c)
        }
    }
    fn cloud(&self, rng: &mut XorShift, n: usize) {
        for i in 0..n {
            let (x, y) = self.point(rng, i);
            self.add(x, y);
        }
    }
}
impl Drop for Pair {
    fn drop(&mut self) {
        unsafe {
            destroy(self.r);
            destroy(self.c);
        }
    }
}
fn choose<T>(m: *mut tri::Model, kind: u8, index: usize) -> *mut T {
    let nodes = unsafe { graph(m) };
    let selected: Vec<_> = nodes.into_iter().filter(|(k, _)| *k == kind).collect();
    selected[index % selected.len()].1.cast()
}
fn sizes(i: usize) -> usize {
    [1, 2, 3, 7, 16, 31, 64, 127, 256][i % 9]
}

#[test]
fn make_model_matches_c() {
    let mut rng = XorShift::new(14);
    for i in 0..3000 {
        let p = if i % 7 == 0 {
            // Degenerate and reversed rectangles are accepted by makeModel.
            let x = value(&mut rng);
            let y = value(&mut rng);
            unsafe { Pair::bounds(x, y, if i % 2 == 0 { x } else { -x }, y) }
        } else {
            Pair::new(&mut rng, i)
        };
        p.compare(true);
    }
}
#[test]
fn add_vertex_to_model_matches_c() {
    let mut rng = XorShift::new(15);
    for i in 0..80 {
        let p = Pair::new(&mut rng, i);
        p.cloud(&mut rng, sizes(i));
        unsafe {
            tri::checkModel(p.r);
            checkModel(p.c);
        }
    }
}

#[test]
fn delete_vertex_from_model_matches_c() {
    let mut rng = XorShift::new(16);
    for i in 0..2048 {
        // Dyadic rectangles avoid the original C deletion loop's numerical
        // nontermination for some almost-collinear boundary configurations.
        let x = (i % 31) as f32 - 16.0;
        let y = (i % 17) as f32 - 8.0;
        let p = unsafe { Pair::bounds(x, y, x + 32.0, y + 32.0) };
        for _ in 0..1 + i % 7 {
            let (x, y) = p.point(&mut rng, 8);
            p.add(x, y);
        }
        let (r, c) = if i % 4 == 0 {
            p.add(x, y + 16.0)
        } else {
            (choose::<Vertex>(p.r, 5, i), choose::<Vertex>(p.c, 5, i))
        };
        unsafe {
            deleteVertexFromModel(p.c, c);
            tri::deleteVertexFromModel(p.r, r);
        }
        p.compare(true);
    }
}

#[test]
fn make_edge_vertex_matches_c() {
    let mut rng = XorShift::new(17);
    for i in 0..2048 {
        let p = Pair::new(&mut rng, i);
        let (x, y) = p.point(&mut rng, 8);
        unsafe {
            let (rf, cf) = ((*p.r).f, (*p.c).f);
            let (rv, cv) = ((*(*rf).eu).vu, (*(*cf).eu).vu);
            let (mut re, mut ce, mut rn, mut cn) = (null_mut(), null_mut(), null_mut(), null_mut());
            tri::makeEdgeVertex((*rv).v, x, y, rf, &mut re, &mut rn);
            makeEdgeVertex((*cv).v, x, y, cf, &mut ce, &mut cn);
            eq_float((*rn).x, (*cn).x);
            eq_float((*rn).y, (*cn).y);
            assert_eq!(re, (*(*rn).vu).eu.as_ref().unwrap().e);
            assert_eq!(ce, (*(*cn).vu).eu.as_ref().unwrap().e);
            p.compare(false);
        }
    }
}

#[test]
fn kill_edge_matches_c() {
    let mut rng = XorShift::new(18);
    for i in 0..2048 {
        let p = Pair::new(&mut rng, i);
        unsafe {
            let interior = |m| {
                graph(m)
                    .into_iter()
                    .find(|(k, p)| {
                        *k == 3 && {
                            let eu = (*p.cast::<tri::Edge>()).eu;
                            !(*eu).f.is_null() && !(*(*eu).euMate).f.is_null()
                        }
                    })
                    .unwrap()
                    .1
                    .cast::<tri::Edge>()
            };
            let (re, ce) = (interior(p.r), interior(p.c));
            let (mut rf, mut cf) = (null_mut(), null_mut());
            tri::killEdge(re, &mut rf);
            killEdge(ce, &mut cf);
            assert_eq!(rf, (*p.r).f);
            assert_eq!(cf, (*p.c).f);
            p.compare(false);
            // C killEdge intentionally leaves Edge allocation to the caller.
            free(re.cast());
            free(ce.cast());
        }
    }
}

#[test]
fn make_edge_face_matches_c() {
    let mut rng = XorShift::new(19);
    for i in 0..2048 {
        let p = Pair::new(&mut rng, i);
        unsafe {
            let (rf, cf) = ((*p.r).f, (*p.c).f);
            let (re, ce) = ((*(*(*rf).eu).euCW).e, (*(*(*cf).eu).euCW).e);
            let (a, b) = ((*(*(*re).eu).vu).v, (*(*(*(*re).eu).euCW).vu).v);
            let (ca, cb) = ((*(*(*ce).eu).vu).v, (*(*(*(*ce).eu).euCW).vu).v);
            let (mut rs, mut cs) = (null_mut(), null_mut());
            tri::killEdge(re, &mut rs);
            killEdge(ce, &mut cs);
            let (mut en, mut cen, mut fn_, mut cfn) =
                (null_mut(), null_mut(), null_mut(), null_mut());
            tri::makeEdgeFace(a, b, &mut en, &mut fn_);
            makeEdgeFace(ca, cb, &mut cen, &mut cfn);
            assert!(!en.is_null() && !fn_.is_null() && !cen.is_null() && !cfn.is_null());
            p.compare(false);
            free(re.cast());
            free(ce.cast());
        }
    }
}

#[test]
fn make_boundary_edge_tri_matches_c() {
    let mut rng = XorShift::new(20);
    for i in 0..2048 {
        let x = (i % 31) as f32 - 16.0;
        let p = unsafe { Pair::bounds(x, -16.0, x + 32.0, 16.0) };
        let (x, y) = p.point(&mut rng, 0);
        let (r, c) = p.add(x, y);
        assert!(!r.is_null());
        unsafe {
            let (mut re, mut ce, mut rt, mut ct) = (null_mut(), null_mut(), null_mut(), null_mut());
            tri::makeBoundaryEdgeTri(r, &mut re, &mut rt);
            makeBoundaryEdgeTri(c, &mut ce, &mut ct);
            eq_float((*rt).rs, (*ct).rs);
            assert_eq!((*re).fixed, (*ce).fixed);
        }
        p.compare(true);
    }
}

#[test]
fn kill_boundary_edge_matches_c() {
    let mut rng = XorShift::new(21);
    for i in 0..2048 {
        let x = (i % 31) as f32 - 16.0;
        let p = unsafe { Pair::bounds(x, -16.0, x + 32.0, 16.0) };
        let (x, y) = p.point(&mut rng, 0);
        let (r, c) = p.add(x, y);
        unsafe {
            let (mut re, mut ce, mut rt, mut ct) = (null_mut(), null_mut(), null_mut(), null_mut());
            tri::makeBoundaryEdgeTri(r, &mut re, &mut rt);
            makeBoundaryEdgeTri(c, &mut ce, &mut ct);
            tri::killBoundaryEdge(re);
            killBoundaryEdge(ce);
        }
        p.compare(true);
    }
}

#[test]
fn fix_edge_between_vertices_matches_c() {
    let mut rng = XorShift::new(22);
    for i in 0..2048 {
        let x = (i % 31) as f32 - 16.0;
        let p = unsafe { Pair::bounds(x, -16.0, x + 32.0, 16.0) };
        for _ in 0..if i % 2 == 0 { 1 + i % 11 } else { 1 } {
            let (x, y) = p.point(&mut rng, 8);
            p.add(x, y);
        }
        let n = unsafe { graph(p.r) }
            .iter()
            .filter(|(k, _)| *k == 5)
            .count();
        let a = i % n;
        let b = (a + 1 + i % (n - 1)) % n;
        let (ra, rb) = constraint_pair(p.r, i, a, b);
        let (ca, cb) = constraint_pair(p.c, i, a, b);
        let cv = unsafe { fixEdgeBetweenVertices(ca, cb) };
        unsafe {
            assert_eq!(tri::fixEdgeBetweenVertices(ra, rb), cv);
        }
        p.compare(true);
    }
}
#[test]
fn unfix_edge_matches_c() {
    let mut rng = XorShift::new(23);
    for i in 0..2048 {
        let x = (i % 31) as f32 - 16.0;
        let p = unsafe { Pair::bounds(x, -16.0, x + 32.0, 16.0) };
        for _ in 0..if i % 2 == 0 { 1 + i % 11 } else { 1 } {
            let (x, y) = p.point(&mut rng, 8);
            p.add(x, y);
        }
        let (ra, rb) = constraint_pair(p.r, 2 * i, i, i + 1);
        let (ca, cb) = constraint_pair(p.c, 2 * i, i, i + 1);
        unsafe {
            assert_eq!(
                tri::fixEdgeBetweenVertices(ra, rb),
                fixEdgeBetweenVertices(ca, cb)
            );
            p.compare(true);
            let (re, ce) = (choose(p.r, 3, i), choose(p.c, 3, i));
            assert_eq!(tri::unfixEdge(re), unfixEdge(ce));
        }
        p.compare(true);
    }
}

#[link(name = "triang", kind = "static")]
extern "C" {
    fn addVertexToModel(m: *mut tri::Model, x: f32, y: f32) -> *mut Vertex;
    fn badModel();
    fn checkEdgeUse(eu: *mut tri::EdgeUse);
    fn checkFace(f: *mut tri::Face);
    fn checkModel(m: *mut tri::Model);
    fn checkVertexUse(vu: *mut tri::VertexUse);
    fn deleteVertexFromModel(m: *mut tri::Model, v: *mut Vertex);
    fn fixEdgeBetweenVertices(v1: *mut Vertex, v2: *mut Vertex) -> i32;
    fn killBoundaryEdge(e: *mut tri::Edge);
    fn killEdge(e: *mut tri::Edge, fs: *mut *mut tri::Face);
    fn makeBoundaryEdgeTri(v: *mut Vertex, enew: *mut *mut tri::Edge, tnew: *mut *mut tri::Tri);
    fn makeEdgeFace(
        v1: *mut Vertex,
        v2: *mut Vertex,
        enew: *mut *mut tri::Edge,
        fnew: *mut *mut tri::Face,
    );
    fn makeEdgeVertex(
        v1: *mut Vertex,
        x: f32,
        y: f32,
        f: *mut tri::Face,
        enew: *mut *mut tri::Edge,
        vnew: *mut *mut Vertex,
    );
    fn readModel(fp: *mut c_void) -> *mut tri::Model;
    fn unfixEdge(e: *mut tri::Edge) -> i32;
    fn updateAddressTree(
        p: *mut tri::Address,
        oa: *mut c_void,
        found: *mut i32,
    ) -> *mut tri::Address;
    fn writeModel(m: *mut tri::Model, fp: *mut c_void);
    fn circum(
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        x3: f32,
        y3: f32,
        xc: *mut f32,
        yc: *mut f32,
        rs: *mut f32,
    );
    fn inCircum(x: f32, y: f32, xc: f32, yc: f32, rs: f32) -> i32;
    fn in3Vertices(x: f32, y: f32, v1: *const Vertex, v2: *const Vertex, v3: *const Vertex) -> i32;
    fn vertexBetweenVertices(v: *const Vertex, v1: *const Vertex, v2: *const Vertex) -> i32;
    fn circumTri(t: *mut tri::Tri);
    fn inCircumTri(x: f32, y: f32, t: *const tri::Tri) -> i32;
    fn inTri(x: f32, y: f32, t: *const tri::Tri) -> i32;
    fn edgesColinear(e1: *const tri::Edge, e2: *const tri::Edge) -> i32;
    fn distanceToEdge(e: *const tri::Edge, x: f32, y: f32) -> f32;
    fn projectToEdge(e: *const tri::Edge, x: *mut f32, y: *mut f32);
    fn makeModel(xmin: f32, ymin: f32, xmax: f32, ymax: f32) -> *mut tri::Model;
    fn killModel(m: *mut tri::Model);
    fn insideTriInModel(m: *mut tri::Model, start: *mut tri::Tri, x: f32, y: f32) -> *mut tri::Tri;
    fn nearestEdgeInModel(
        m: *mut tri::Model,
        start: *mut tri::Edge,
        x: f32,
        y: f32,
    ) -> *mut tri::Edge;
    fn nearestVertexInModel(m: *mut tri::Model, start: *mut Vertex, x: f32, y: f32) -> *mut Vertex;
}

extern "C" {
    fn calloc(count: usize, size: usize) -> *mut c_void;
    fn free(p: *mut c_void);
    fn fclose(fp: *mut c_void) -> i32;
    fn rewind(fp: *mut c_void);
    fn fread(p: *mut c_void, size: usize, count: usize, fp: *mut c_void) -> usize;
    fn fwrite(p: *const c_void, size: usize, count: usize, fp: *mut c_void) -> usize;
}
fn value(r: &mut XorShift) -> f32 {
    ((r.next_f64() * 2.0 - 1.0) * 1000.0) as f32
}
fn vertex(x: f32, y: f32) -> Vertex {
    Vertex {
        vu: std::ptr::null_mut(),
        fixed: 0,
        x,
        y,
        va: std::ptr::null_mut(),
    }
}
#[test]
fn circum_matches_c() {
    let mut r = XorShift::new(1);
    for i in 0..20000 {
        let a = [
            value(&mut r),
            value(&mut r),
            value(&mut r),
            value(&mut r),
            value(&mut r),
            value(&mut r),
        ];
        let a = if i % 13 == 0 { [0.0; 6] } else { a };
        let (mut x, mut y, mut s) = (0.0, 0.0, 0.0);
        let (mut cx, mut cy, mut cs) = (0.0, 0.0, 0.0);
        tri::circum(a[0], a[1], a[2], a[3], a[4], a[5], &mut x, &mut y, &mut s);
        unsafe {
            circum(
                a[0], a[1], a[2], a[3], a[4], a[5], &mut cx, &mut cy, &mut cs,
            );
        }
        for (u, v) in [(x, cx), (y, cy), (s, cs)] {
            assert_eq!(u.to_bits(), v.to_bits());
        }
    }
}
#[test]
fn in_circum_matches_c() {
    let mut r = XorShift::new(2);
    for i in 0..20000 {
        let (x, y, cx, cy) = (value(&mut r), value(&mut r), value(&mut r), value(&mut r));
        let rs = if i % 7 == 0 {
            f32::MAX
        } else {
            value(&mut r) * 1000.0
        };
        assert_eq!(tri::inCircum(x, y, cx, cy, rs), unsafe {
            inCircum(x, y, cx, cy, rs)
        });
    }
}
#[test]
fn in_3_vertices_matches_c() {
    let mut r = XorShift::new(3);
    for _ in 0..20000 {
        let (x, y) = (value(&mut r), value(&mut r));
        let a = vertex(value(&mut r), value(&mut r));
        let b = vertex(value(&mut r), value(&mut r));
        let c = vertex(value(&mut r), value(&mut r));
        assert_eq!(tri::in3Vertices(x, y, &a, &b, &c), unsafe {
            in3Vertices(x, y, &a, &b, &c)
        });
    }
}
#[test]
fn vertex_between_vertices_matches_c() {
    let mut r = XorShift::new(4);
    for i in 0..20000 {
        let a = vertex(value(&mut r), value(&mut r));
        let b = vertex(value(&mut r), value(&mut r));
        let v = if i % 11 == 0 {
            vertex(a.x, a.y)
        } else {
            vertex(value(&mut r), value(&mut r))
        };
        assert_eq!(tri::vertexBetweenVertices(&v, &a, &b), unsafe {
            vertexBetweenVertices(&v, &a, &b)
        });
    }
}

struct Fixture {
    vertices: Box<[Vertex; 3]>,
    uses: Box<[tri::VertexUse; 3]>,
    edges: Box<[tri::EdgeUse; 3]>,
    edge: Box<tri::Edge>,
    face: Box<tri::Face>,
}
impl Fixture {
    fn new() -> Self {
        let mut f = Self {
            vertices: Box::new(std::array::from_fn(|_| vertex(0.0, 0.0))),
            uses: Box::new(std::array::from_fn(|_| tri::VertexUse {
                eu: std::ptr::null_mut(),
                vuPrev: std::ptr::null_mut(),
                vuNext: std::ptr::null_mut(),
                v: std::ptr::null_mut(),
                vua: std::ptr::null_mut(),
            })),
            edges: Box::new(std::array::from_fn(|_| tri::EdgeUse {
                f: std::ptr::null_mut(),
                vu: std::ptr::null_mut(),
                euMate: std::ptr::null_mut(),
                euCW: std::ptr::null_mut(),
                euCCW: std::ptr::null_mut(),
                e: std::ptr::null_mut(),
                eua: std::ptr::null_mut(),
            })),
            edge: Box::new(tri::Edge {
                eu: std::ptr::null_mut(),
                fixed: 0,
                ea: std::ptr::null_mut(),
            }),
            face: Box::new(tri::Face {
                m: std::ptr::null_mut(),
                fPrev: std::ptr::null_mut(),
                fNext: std::ptr::null_mut(),
                eu: std::ptr::null_mut(),
                xc: 0.0,
                yc: 0.0,
                rs: 0.0,
                fa: std::ptr::null_mut(),
            }),
        };
        for i in 0..3 {
            f.uses[i].v = &mut f.vertices[i];
            f.edges[i].vu = &mut f.uses[i];
        }
        let p = f.edges.as_mut_ptr();
        for i in 0..3 {
            unsafe {
                f.edges[i].euCW = p.add((i + 1) % 3);
                f.edges[i].euCCW = p.add((i + 2) % 3);
            }
        }
        f.edge.eu = p;
        f.face.eu = p;
        f
    }
    fn randomize(&mut self, r: &mut XorShift, i: usize) {
        for v in self.vertices.iter_mut() {
            v.x = value(r);
            v.y = value(r);
        }
        if i % 13 == 0 {
            self.vertices[1].x = self.vertices[0].x;
            self.vertices[1].y = self.vertices[0].y;
        }
    }
}
fn eq_float(a: f32, b: f32) {
    assert!(
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
        "{a} != {b}"
    );
}
#[test]
fn circum_tri_matches_c() {
    let mut r = XorShift::new(5);
    let mut f = Fixture::new();
    for i in 0..20000 {
        f.randomize(&mut r, i);
        unsafe { tri::circumTri(&mut f.face) };
        let a = (f.face.xc, f.face.yc, f.face.rs);
        unsafe { circumTri(&mut *f.face) };
        eq_float(a.0, f.face.xc);
        eq_float(a.1, f.face.yc);
        eq_float(a.2, f.face.rs);
    }
}
#[test]
fn in_circum_tri_matches_c() {
    let mut r = XorShift::new(6);
    let mut f = Fixture::new();
    for i in 0..20000 {
        f.face.xc = value(&mut r);
        f.face.yc = value(&mut r);
        f.face.rs = if i % 7 == 0 {
            f32::MAX
        } else {
            value(&mut r) * 1000.0
        };
        let x = value(&mut r);
        let y = value(&mut r);
        assert_eq!(tri::inCircumTri(x, y, &f.face), unsafe {
            inCircumTri(x, y, &*f.face)
        });
    }
}
#[test]
fn in_tri_matches_c() {
    let mut r = XorShift::new(7);
    let mut f = Fixture::new();
    for i in 0..20000 {
        f.randomize(&mut r, i);
        let x = value(&mut r);
        let y = value(&mut r);
        assert_eq!(unsafe { tri::inTri(x, y, &f.face) }, unsafe {
            inTri(x, y, &*f.face)
        });
    }
}
#[test]
fn edges_colinear_matches_c() {
    let mut r = XorShift::new(8);
    let mut f = Fixture::new();
    let mut g = Fixture::new();
    for i in 0..20000 {
        f.randomize(&mut r, i);
        g.randomize(&mut r, i);
        assert_eq!(unsafe { tri::edgesColinear(&f.edge, &g.edge) }, unsafe {
            edgesColinear(&*f.edge, &*g.edge)
        });
    }
}
#[test]
fn distance_to_edge_matches_c() {
    let mut r = XorShift::new(9);
    let mut f = Fixture::new();
    for i in 0..20000 {
        f.randomize(&mut r, i);
        let x = value(&mut r);
        let y = value(&mut r);
        eq_float(unsafe { tri::distanceToEdge(&f.edge, x, y) }, unsafe {
            distanceToEdge(&*f.edge, x, y)
        });
    }
}
#[test]
fn project_to_edge_matches_c() {
    let mut r = XorShift::new(10);
    let mut f = Fixture::new();
    for i in 0..20000 {
        f.randomize(&mut r, i);
        let mut x = value(&mut r);
        let mut y = value(&mut r);
        let (mut cx, mut cy) = (x, y);
        unsafe {
            tri::projectToEdge(&f.edge, &mut x, &mut y);
            projectToEdge(&*f.edge, &mut cx, &mut cy)
        };
        eq_float(x, cx);
        eq_float(y, cy);
    }
}
#[test]
fn inside_tri_in_model_matches_c() {
    query_suite(0);
}
#[test]
fn nearest_edge_in_model_matches_c() {
    query_suite(1);
}
#[test]
fn nearest_vertex_in_model_matches_c() {
    query_suite(2);
}
