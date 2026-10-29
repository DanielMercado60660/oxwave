# Audit: triangulated/tetrahedral model family

Scope: crates `oxwave-lib-tri`, `oxwave-lib-trielas`, `oxwave-lib-tetra`; programs `tetramod`, `sutetraray`.

## How this was checked

- **Golden tests:** `tools/check_program.py` passes for both programs (tetramod 54/54, sutetraray 64/64).
- **Reading:** every file in scope was compared expression by expression against the C (`vendor/su-original/src/{tri,Trielas,tetra}/...`).
- **Differential fuzzing:** all runs were outside the repo, in `/private/tmp/claude-501/{tri-audit,tetramod-audit,sutetraray-audit,trielas-audit}`. Nothing in the repo was edited and cargo was never run in the repo.

| Target | Runs | Result |
|---|---|---|
| tetramod | about 750 random cases. Ranges: nxhz 2..40, nyhz 3..41, nhz 1..6, crossing or zero-thickness horizons, irregular binary grids, ficth, verbose. | Compared stdout, stderr, exit code and `hz.bin` byte for byte. |
| sutetraray | about 500 random cases. Ranges: models with nhz 1..4, dipping corners, velocity gradients, irefseq, shootupward, every aux file. | About 300 runs where C exits cleanly matched byte for byte. |
| lib-tri, through `uni2tri \| tri2uni` | 160 random smooth and stepped fields (n1, n2 up to 30), methods 1–3 and various errmax/tol. | 0 diffs, f32 compare. Exercises `makeModel`, `addVertexToModel`, `nearest*`, `insideTriInModel`, `readModel`/`writeModel`. |
| lib-tri, through `trimodel \| tri2uni/triray` | 160 random interface sets. | The only divergence on input where C does not crash is M4. About 9 cases had differing stderr. I checked two of them individually, and both crash on both sides through recursion in the `trimodel` program (outside this scope). The other seven share the same self-crossing-interface pattern but were not checked one by one. |
| lib-trielas | A bit-exact C-vs-Rust harness linking `libtrielas.a`. 200k physical-range cases per function, 300k edge-value cases, and `rt_ani_real` with the original elaray.c helpers compiled `-ffp-contract=off`. | 11 of 12 functions are bitwise identical. The exception is M3. |
| lib-tetra | Read by hand. `tetra_volume` and `area3d` are faithful, and their diff test calls the real C. | No issues. |

- **Hygiene:**
  - Every `.rs` in scope starts with `tools/cwp_header.txt` verbatim.
  - There is no `#[ignore]` and no golden special-casing.
  - Every library `tests/diff.rs` really calls the C library. The one exception is `unfixEdge` (see Low 1). lib-tri compares bit-exact topology snapshots of whole models built by C and Rust.

---

## HIGH

### H1. sutetraray: when a ray reaches `nseg == maxnsegs`, Rust always segfaults, so `err("maxnsegs too small")` can never be reached
- **Rust:** `programs/sutetraray/src/main.rs:7303-7316` (`array2` leaves `rows[n2]` null on purpose: "A null guard row preserves its SIGSEGV"). It is used for `xseg` at main.rs:1466.
- **C:** `sutetraray.c:840` does `ealloc2float(3,maxnsegs)`, but the checks at `sutetraray.c:2736-2740` and `2796-2800` are `if (ray->nseg>maxnsegs) err(...)`. So C writes `xseg[maxnsegs]` (one row past the end) before it errors.
  - In the reference build that row pointer is a valid neighbouring heap pointer. C therefore either completes or reaches its own error message.
  - Rust dereferences NULL instead.
- **Repro** (verified):
  ```
  tetramod nxhz=3 nyhz=3 nhz=1 hzfile=hz.bin > model.txt
  sutetraray maxnsegs=16 < model.txt
  ```
  - C: rc=0 with 6804 bytes of traces, byte-identical to the default `maxnsegs=100` run.
  - Rust: rc=139 with 0 bytes.
- **Second repro** (sutetraray sub-audit, 8/8 runs): `tetramod nxhz=3 nyhz=7 nhz=1 z00=0.0 v00=1.203 z01=0.013 v01=1.109 dvdz01=0.02 z10=0.0 v10=1.333 dvdz10=0.395 z11=0.0 v11=1.323 dvdz11=0.245 hzfile=hz.bin`, then `sutetraray nxgd=5 nygd=3 takeoff=45 nazimuth=24 maxnsegs=20`.
  - C: rc=1 with `maxnsegs too small`.
  - Rust: rc=139.
  - So a user who sets maxnsegs too low gets a crash instead of the documented error.
- **Why the goldens miss it:** golden `maxnsegs=1` crashes in both. The guard was evidently chosen to match that case.
- **Fix:** set `rows[n2] = data.add(n2*n1)`. The data block and `tseg` already reserve the `+1` slot.
  - After the fix, `nseg == maxnsegs` succeeds and `nseg > maxnsegs` prints the C error.
  - Remaining caveat: C's stray row aliases another ray's segment 0, so `rayfile=` output can still differ in that corner.

### H2. tetramod: when `ficth` marks every layer at or below `ihz`, `iiihz = -1` and Rust panics (high, borderline)
- **Rust:** `programs/tetramod/src/main.rs:754` stores `iihz[ihz] = -1`. Line 813 then indexes `v01[iiihz as usize]`, which panics with index 18446744073709551615 and exit 101.
- **C:** `tetramod.c:601-606` leaves `iiihz = -1`. Lines 638-650 then read `v01[-1]`, `v0hztem[-1]`, `zflat[-1]` out of bounds.
- **Repro** (verified):
  ```
  tetramod nxhz=2 nyhz=41 hzfile=hz.bin nhz=4 ficth=2,1,0 blt=0
  ```
  - C: rc=0, a full model built from heap garbage, deterministic over 5/5 runs.
  - Rust: rc=101, `index out of bounds: the len is 4 but the index is 18446744073709551615` plus the Rust panic note.
- **Variant:** `tetramod nxhz=10 nyhz=31 hzfile=hz.bin nhz=4 ficth=0,1,3 xmin=1 xmax=1.5` gives C SIGSEGV (139) vs Rust 101 with panic text. More than 10 fuzz cases behave like this.
- **Why high:** the input is degenerate, and when C "succeeds" its model is built from garbage. I still rate it high because `ficth=` is a documented user parameter, and Rust gives a Rust panic, which C never does, on inputs the harness accepts. Downgrading to medium would also be reasonable.
- **Fix:** never panic. When the effective `iiihz < 0`, call the existing `c_segv()`, which matches the common C outcome. Document that the occasional C "success on garbage" cannot be reproduced.

### H3. tetramod: C's getpar string aliasing is not reproduced, so Rust reads grid files that C ignores
- **C:** `tetramod.c:478-482` (also 488, 536, 561).
  - `getparstring(xfile,&xfile)` repoints `xfile` into getpar's argument buffer.
  - The next horizon's `sprintf(xfile,"x%dfile",ihz)` then writes into that value, spilling over the next argv entry's name.
  - x/y/z/v names are 6 chars, so this triggers when the previous value is ≤5 chars. dvdz names are 9 chars, so it triggers when the previous value is ≤8 chars.
- **Rust:** `programs/tetramod/src/main.rs:520-528, 624, 658` builds a fresh String, so there is no aliasing.
- **Repro** (verified; `z.bin` = 9×0.1f, `w.bin` = 9×0.5f):
  ```
  tetramod nxhz=3 nyhz=3 hzfile=hz.bin nhz=2 z0file=z.bin z1file=w.bin
  ```
  - C: `Using uniform cartesian grids for z-horizon1`, then `1.600000=zmax`.
  - Rust: `read in w.bin`, then `2.100000=zmax`. Points, facets and hz.bin all differ.
- **dvdz variant:** `dvdz0file=dvdz.bin dvdz1file=w.bin` gives 33 differing stdout lines.
- **Why the goldens miss it:** they use `xgrid.bin` and `dvdzgrid.bin`, which are long enough not to spill.
- **Disposition:** C is the buggy side here, since it corrupts its own argument table. Rust's behaviour is the sane one. But plausible short file names give different models.
- **Recommendation:** record this as an accepted deviation (UB exclusion) rather than emulate argv corruption, which would need getpar internals in the core crate.

## MEDIUM

### M1. sutetraray: `ifwf2dump`/`nwf2dump` are zeroed in Rust but read uninitialized in C
- **Rust:** `programs/sutetraray/src/main.rs:304-308`. **C:** `sutetraray.c:297-298`. They are set only when `wffile` is open (398-401), but read at 404 (verbose print) and 487 (ntwf check).
- **Repro** (verified):
  ```
  tetramod nxhz=3 nyhz=3 nhz=1 hzfile=hz.bin > model.txt
  sutetraray verbose=1 < model.txt
  ```
  stderr line 14 differs:
  - C: `ifwf2dump=1, nwf2dump=-193501640`, stable over runs.
  - Rust: `ifwf2dump=0, nwf2dump=0`.
- **ntwf variant** (sub-audit): with `ntwf=-1` and no wffile, Rust prints an extra `ntwf < ifwf2dump+nwf2dump, reset=0` and C does not.
- **Why the goldens miss it:** every golden `verbose=1` case passes `wffile=`.
- **Fix:** either hard-code the observed reference stack values for parity, or list this as a UB exclusion.

### M2. sutetraray: C error paths for `maxnsegs <= 0` are missing
- **maxnsegs=0:**
  - C: rc=1, `maxnsegs too small`. Rust: rc=139.
  - Cause: CWP `alloc2` always sets `p[0]` (cwp/lib/alloc.c:255), but Rust `array2` with n2=0 leaves row 0 null.
- **maxnsegs=-1:**
  - C: rc=1, `ealloc.c: malloc failed (Cannot allocate memory)`. Rust: SIGBUS (138).
  - Cause: `-1 as usize` wraps, and `allocate`/`array2` do not reproduce the ealloc failure message.
  - The same file at main.rs:7295 `unwrap()`s a `LayoutError`, which gives a Rust panic (exit 101) for `maxnrays=-1`, where C gives SIGSEGV.
- **Verified** on the golden model (`tetramod nxhz=3 nyhz=3 nhz=1`):
  - `maxnsegs=0`: C rc=1 with `maxnsegs too small`; Rust rc=139.
  - `maxnsegs=-1`: C rc=1 with `ealloc.c: malloc failed (Cannot allocate memory)`; Rust rc=138.
  - `maxnrays=-1`: C rc=139; Rust rc=101 with `panicked at programs/sutetraray/src/main.rs:7295:85` (`unwrap` on a `LayoutError`).
- **Fix:** the H1 fix, plus a valid row 0 when n2=0, plus emulation of ealloc's error for negative or overflowing sizes.

### M3. lib-trielas: `rt_iso_cmplx` coefficient is 1 ulp off in about 0.7% of cases, and the diff-test tolerance hides it
- **Rust:** `crates/oxwave-lib-trielas/src/lib.rs:608, 652-664` call the core `rcabs` (`crates/oxwave-cwp/src/complex.rs:111-125`), which squares the ratio in f64.
- **C:** `cwp/lib/complex.c` `rcabs` computes `temp*temp` in float, then promotes.
- **Inconsistency:** the port has a local C-faithful `reciprocal` (lib.rs:559-565) because core `cinv` is unfaithful, but no such local `rcabs`.
- **Repro** (bit-exact harness against `libtrielas.a`): 2840/400000 mismatches. For example vp1=3.5854886, vp2=2.0740638, vs1=1.0834757, vs2=1.3912363, rho1=2.14403, rho2=2.0636044, pl=0.47742823, modei=1, modet=0, rort=1:
  - C coeff=0.11859905
  - Rust coeff=0.11859904
  - Swapping in a C-faithful local `rcabs` gives 0/400000 mismatches.
- **Effect:** the `rcabs(dd) < FLT_EPSILON` pole guard could in principle flip, though this was not observed.
- **Test tolerance:** `tests/diff.rs:212-227` `same()` uses `1e-5 * max|v|` over [coeff, phase]. A phase near 2π then sets the scale, so a small coeff can carry about 0.6% slack. The other 11 functions are bitwise identical, so `exact()` could be used everywhere.
- **Fix:** add a local f32-squaring `rcabs` closure next to `reciprocal`, and switch the tests to `exact()`.
- **Note:** the same core `rcabs` imprecision is documented in DESIGN.md as at most 1.3e-7 relative, so this is within harness tolerance. It is still an avoidable precision change.

### M4. lib-tri/trimodel: an inherited C use-after-free makes Rust hang nondeterministically where C deterministically exits 0
- **Repro** (verified):
  ```
  trimodel xedge=0,1 xedge=1,0 zedge=0,0 zedge=1,1 sedge=1,1 sedge=1,1 \
    xedge=0.3585,0.1314,0.2146,0.3341,0.9570 zedge=0.5518,0.3053,0.1453,0.4095,0.3232 \
    sedge=2.766,1.869,1.608,2.677,2.913 maxangle=1
  ```
  - C: rc=0, 454888 bytes, 6/6 runs, also 3/3 under `MallocScribble=1`.
  - Rust: hangs in 4/6 runs and completes in 2/6. It did not hang in the few lldb runs I tried, where ASLR is off.
- **Where it hangs:** `sample` and an lldb attach place the loop in `oxwave_lib_tri::create::addVertexToModel` (create.rs:642, the inlined `nearestVertexInModel` walk). The query point is (0.1445, 0.2200).
- **Root cause is shared C UB:** under Guard Malloc (`DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib`), both binaries fault at the same place.
  - C faults at `nearest.c:124`, `vmin = m->f->eu->vu->v`, called from create.c:530.
  - The model's head face `m->f` is live, but its `eu` points to a freed EdgeUse.
  - In C the freed block is promptly reused by a same-size allocation, so the walk survives. Rust's allocation pattern differs (trimodel's attribute `Box`es, `calloc`), so the stale pointer lands on garbage and the vertex-use ring never closes.
  - The likely origin is `killEdge`'s `f->eu = eu->euCW->euCW` (delete.c:303, delete.rs:249) on a degenerate face, but this is not proven.
- **Caveat:** the C output for this input is itself corrupt. `tri2uni n1=50 n2=50` on C's model also loops forever. So this is an edge case. The audit of the tri programs (claude-wave3-tri-progs.md) separately notes trimodel hangs on curved edges.
- **Fix:** not a translation error. Either document it as UB, or validate `m->f->eu` (or make `killEdge` keep `f->eu` pointing at a surviving edge-use) on both sides of the comparison.

## LOW (brief)

1. **lib-tri `unfixEdge` deliberately departs from C** (`fixedges.rs:289-303` vs `fixedges.c:280-281`). It saves both endpoints before C's use-after-free re-read. Its test (`tests/diff.rs:45`) checks against a Delaunay oracle rather than the C function. No SU program calls `unfixEdge` or `deleteVertexFromModel`, so this is harmless.
2. **lib-tri exits through libc `exit()`** (`lib.rs:58`; `readwrite.rs:56-58` `fail`; `check.rs:44` `badModel`). This bypasses `oxwave_core::flush_all`, so any oxwave-core-buffered stdout written before a `readModel` or `badModel` error is lost, whereas C flushes it. Malformed-model tests of tri2uni gave identical results today, because programs read the model before writing output.
3. **tetramod emulation of C undefined behaviour:**
   - `a[2][7]` is never set in C (tetramod.c:1967 assigns `a[2][6]` twice); Rust uses 0 (main.rs:2253).
   - `region()` (main.rs:2514-2560) emulates the macOS malloc layout to reproduce C's `tetra[-1]`/`tetra[-2]` reads.
   - The crash paths for nhz≤0 and missing files are hand-pinned.
   - All of these are honest, but platform-pinned.
4. **lib-trielas:** private `findqPqSV` drops C's `err(" wrong mode in findqP \n ")` (ani.rs:137-141). It is unreachable, since only modes 3 and 4 are passed. No diff test exercises any `err!` path (invalid rort/mode), so the verbatim messages are untested.
5. **Aux files on crash (core, seen via sutetraray):** C has flushed 4 KB blocks of `wf.txt`/`ray.txt` before a crash; Rust's 64 KB `cio` buffer leaves them empty. Goldens never compare aux files on crashing runs.

## Unconfirmed / out of scope
- C nondeterminism: the vendor C `tetramod` and `sutetraray` sometimes die with SIGBUS when the `tetra[-1]` read lands on an unmapped page (seen with `tetramod nxhz=7 nyhz=5 nhz=2`, 3/12 runs). Goldens on larger models will be flaky on the C side. This is not a Rust divergence.
- `trimodel` with crossing interfaces: both sides crash, but Rust aborts with "has overflowed its stack" (SIGABRT, 134) where C gets SIGSEGV (139). The recursion is in the trimodel and cubic-interp code, not in this scope.
- sutetraray: the C code leaves the last wavefront slice and the ray array uninitialized (malloc); Rust zeroes them. No observable difference was found.

## Overall judgement
The arithmetic is very faithful across the whole family:
- lib-tri reproduces C topology bit for bit (verified by strong whole-model snapshot tests and 160 clean uni2tri/tri2uni fuzz runs).
- lib-trielas is bitwise identical in 11 of 12 functions.
- lib-tetra is exact.
- tetramod and sutetraray matched the C byte for byte in roughly 1000 clean differential runs, including multi-layer, dipping and reflecting models.

No test-specific code, golden special-casing or license problems were found. Every real divergence sits where the C itself has undefined behaviour: out-of-bounds rows, uninitialized locals, getpar buffer aliasing, negative layer indices, use-after-free. That is where the ports' hand-written "emulate the crash" choices sometimes go wrong:
- a null guard row that turns the documented `maxnsegs too small` error into a segfault (H1);
- a Rust index panic in place of C's behaviour (H2);
- zeroed uninitialized variables that change verbose output (M1).

These are cheap to fix or to document as UB exclusions. The numerical core needs no rework.
