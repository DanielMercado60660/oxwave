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

//! Port of cwp/lib/sort.c. Author: Dave Hale, Colorado School of Mines,
//! 12/26/89 and 01/13/89; quicksort adapted from C. A. R. Hoare.

use std::sync::Mutex;
static VALUE_SEED: Mutex<i64> = Mutex::new(0);
static INDEX_SEED: Mutex<i64> = Mutex::new(0);

/// Sort floats using a heap. C prototype: `void hpsort(int n, float a[]);`
pub fn hpsort(a: &mut [f32]) {
    let mut root = a.len() as isize / 2 - 1;
    let mut last = a.len() as isize - 1;
    while last > 0 {
        let mut parent;
        let value;
        if root >= 0 { parent = root; value = a[root as usize]; root -= 1; }
        else { parent = 0; value = a[last as usize]; a[last as usize] = a[0]; last -= 1; }
        let mut child = parent * 2 + 1;
        while child <= last {
            if child < last && a[(child + 1) as usize] > a[child as usize] { child += 1; }
            if a[child as usize] > value { a[parent as usize] = a[child as usize]; parent = child; }
            else { break; }
            child = parent * 2 + 1;
        }
        a[parent as usize] = value;
    }
}

// The two C partition helpers have independent static seeds.
trait Items {
    type Item: Copy;
    fn value(&self, pos: usize) -> f32;
    fn get(&self, pos: usize) -> Self::Item;
    fn put(&mut self, pos: usize, item: Self::Item);
    fn item_value(&self, item: Self::Item) -> f32;
    fn swap(&mut self, a: usize, b: usize) {
        let t = self.get(a); self.put(a, self.get(b)); self.put(b, t);
    }
}
struct Values<'a>(&'a mut [f32]);
impl Items for Values<'_> {
    type Item = f32;
    fn value(&self, p: usize) -> f32 { self.0[p] }
    fn get(&self, p: usize) -> f32 { self.0[p] }
    fn put(&mut self, p: usize, x: f32) { self.0[p] = x; }
    fn item_value(&self, x: f32) -> f32 { x }
}
struct Indices<'a> { a: &'a [f32], i: &'a mut [i32] }
impl Items for Indices<'_> {
    type Item = i32;
    fn value(&self, p: usize) -> f32 { self.a[self.i[p] as usize] }
    fn get(&self, p: usize) -> i32 { self.i[p] }
    fn put(&mut self, p: usize, x: i32) { self.i[p] = x; }
    fn item_value(&self, x: i32) -> f32 { self.a[x as usize] }
}
fn partition<T: Items>(a: &mut T, p: isize, q: isize, seed: &mut i64) -> (isize, isize) {
    *seed = (*seed * 211 + 1663) % 7875;
    // C evaluates the entire pivot expression in float, including p.
    let pivot = ((p as f32) + ((q - p) as f32 * *seed as f32) / 7875.0) as isize;
    let pivot = pivot.clamp(p, q);
    let value = a.value(pivot as usize);
    let (mut left, mut right) = (p, q);
    loop {
        while a.value(left as usize) <= value && left < q { left += 1; }
        while a.value(right as usize) >= value && right > p { right -= 1; }
        if left < right { a.swap(left as usize, right as usize); left += 1; right -= 1; }
        else { break; }
    }
    if left < pivot { a.swap(left as usize, pivot as usize); left += 1; }
    else if pivot < right { a.swap(right as usize, pivot as usize); right -= 1; }
    (right, left)
}
fn insertion<T: Items>(a: &mut T, p: isize, q: isize) {
    for j in p + 1..=q {
        let item = a.get(j as usize);
        let value = a.item_value(item);
        let mut k = j;
        while k > p && a.value((k - 1) as usize) > value {
            a.put(k as usize, a.get((k - 1) as usize)); k -= 1;
        }
        a.put(k as usize, item);
    }
}
fn quicksort<T: Items>(a: &mut T, n: usize, seed: &mut i64) {
    let mut stack = vec![(0isize, n as isize - 1)];
    while let Some((mut p, mut q)) = stack.pop() {
        while q - p > 7 {
            let (j, k) = partition(a, p, q, seed);
            if j - p < q - k { stack.push((k, q)); q = j; }
            else { stack.push((p, j)); p = k; }
        }
        insertion(a, p, q);
    }
}
fn quickfind<T: Items>(a: &mut T, n: usize, m: usize, seed: &mut i64) {
    let (mut p, mut q) = (0isize, n as isize - 1);
    while q - p > 7 {
        let (j, k) = partition(a, p, q, seed);
        if m as isize <= j { q = j; }
        else if m as isize >= k { p = k; }
        else { return; }
    }
    insertion(a, p, q);
}

/// Sort floats using quicksort. C prototype: `void qksort(int n, float a[]);`
pub fn qksort(a: &mut [f32]) {
    let n = a.len(); quicksort(&mut Values(a), n, &mut VALUE_SEED.lock().unwrap_or_else(|e| e.into_inner()));
}
/// Partially sort floats to select `m`. C prototype: `void qkfind(int m, int n, float a[]);`
pub fn qkfind(m: usize, a: &mut [f32]) {
    assert!(m < a.len()); let n = a.len();
    quickfind(&mut Values(a), n, m, &mut VALUE_SEED.lock().unwrap_or_else(|e| e.into_inner()));
}
/// Sort indices by their corresponding values. C prototype: `void qkisort(int n, float a[], int i[]);`
/// Each index must address `a`; the number sorted is `i.len()`.
pub fn qkisort(a: &[f32], i: &mut [i32]) {
    let n = i.len(); quicksort(&mut Indices { a, i }, n, &mut INDEX_SEED.lock().unwrap_or_else(|e| e.into_inner()));
}
/// Partially sort indices to select `m`. C prototype: `void qkifind(int m, int n, float a[], int i[]);`
/// Each index must address `a`; the number considered is `i.len()`.
pub fn qkifind(m: usize, a: &[f32], i: &mut [i32]) {
    assert!(m < i.len()); let n = i.len();
    quickfind(&mut Indices { a, i }, n, m, &mut INDEX_SEED.lock().unwrap_or_else(|e| e.into_inner()));
}

