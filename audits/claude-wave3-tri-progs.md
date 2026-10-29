# Audit: triangulated-model programs

Scope: elacheck, elamodel, elaps, elaray, elatriuni, gbbeam, normray, tri2uni,
trimodel, triray, triseis, uni2tri (plus the crates they pull in:
oxwave-lib-tri, oxwave-lib-trielas, oxwave-f-basic/psaxesbox, the interp/spline crates).

Method: the goldens for these programs almost all use trivial models (a unit
square, homogeneous sloth, or a 2x2 `uni2tri` grid). So I ran the **stock SU
demo pipelines** (`demos/Synthetic/Tri/Seismograms/PSdemo1`,
`Tri/NormalRays/PSnormrays`, `Trielas/Elasyn`, `Trielas/Waves`, `Trielas/DisplModel`)
and parameter batteries on those heterogeneous models through both the C
binaries (`vendor/su-build/bin` first on PATH) and the Rust ones (a shim dir with
only the 12 in-scope Rust binaries ahead of `vendor/su-build/bin`). I compared every
output file (raw floats within 1e-5 rel, text/PS byte-exact). I also ran an
edge-input battery (negative or zero sizes, unknown strings, too few values).
Experiments are in `/private/tmp/claude-501/triaudit` (`cmp.sh`, `edge.sh`, `fcmp.py`).
The existing goldens pass today (`check_program.py`: gbbeam 110/110, tri2uni 13/13,
elatriuni 15/15, elaps 10/10).

Counts: **high 0, medium 2, low 2.**

---

## M1 (medium): gbbeam (and triseis) receiver geometry depends on C heap garbage, so the stock demo gives different seismograms

- Program: gbbeam. triseis has the same C code path.
- Rust: `programs/gbbeam/src/main.rs:265-280` (`makexzs` uses a zero guard for `s[-1]`); `programs/triseis/src/main.rs:299-312` (same guard).
- C: `vendor/su-original/src/tri/main/gbbeam.c:182`, `vendor/su-original/src/tri/main/triseis.c:275`:
  `s[iuu] = s[iuu-1]+sqrt(dx*dx+dz*dz);` with `iuu==0` reads `s[-1]` (heap memory before the allocation).
- What differs: C sets `s[0] = s[-1]`, and that value is whatever sits in the heap.
  Usually it is ~0 (1.4e-45), but once the rayend input is large enough (>= 94 records
  here, i.e. >8 KB) it is a real float. Then `ds=(s[nuu-1]-s[0])/(ns-1)` is wrong and every
  receiver is mislocated. Rust always uses 0, the intended value. Rust's gbbeam output
  equals the sum of single-ray runs (to 1e-7) and does not depend on ray order. The C output
  is not linear in the rays, depends on ray order, and even changes when `infofile=` is
  added (one more FILE allocation shifts the heap). An instrumented rebuild of `gbbeam.c`
  (same flags, same libs; it reproduces the original binary's md5) shows the cause directly.
  With all 130 rays: `s[-1]=1.2124` and the receivers come out at
  `xg = 0, 0, 0.1209, 0.7876, 1.454, 2.121, 2.788` instead of
  `0, 0.6667, 1.333, 2, 2.667, 3.333, 4`.
- Reproduction (stock PSdemo1 pipeline):
  ```
  M='trimodel xmin=-1 zmin=0 xmax=5.0 zmax=2.0 xedge=-1,0,1,2,3,4,5 zedge=0,0,0,0,0,0,0 sedge=0,0,0,0,0,0,0 xedge=-1,0,0.5,1.0,1.5,2,2.5,3,4,5 zedge=1,1.1,1.2,1.5,1.7,1.5,1.2,1.1,1,1 sedge=0,0,0,0,0,0,0,0,0,0 xedge=-1,0,1,2,3,4,5 zedge=2,2,2,2,2,2,2 sedge=0,0,0,0,0,0,0 sfill=1,0.5,0,0,0.25,0,0 sfill=1,1.8,0,0,0.1,0,0 kedge=1,2,3'
  $M > model.data
  triray <model.data >rayends.data xs=1.5 zs=0 rayfile=rays.data krecord=1 nxz=201 nangle=130 fangle=-65 langle=65 refseq=2,1 refseq=1,-1
  gbbeam <rayends.data bw=0 xg=0,4 zg=0,0 ng=100 lscale=0.5 krecord=1 fpeak=8 nt=101 dt=0.025 ft=0 > seis.bin
  ```
  model.data, rayends.data and rays.data are identical between C and Rust (rayends byte-identical).
  seis.bin: 9855 of 10100 samples differ, worst |diff| 1.389 at sample 4919
  (C 1.1386, Rust -0.2507; max |amp| 1.275). With `head -c $((93*88))` of the rayends,
  C and Rust agree to 1e-7. Adding `infofile=/dev/null` to the 94-ray C run also makes it agree.
  A second repro: `gbbeam <rayends.data xg=0,4 zg=0,0 ng=7 nt=101 dt=0.02 krecord=1 atten=1 fpeak=10 aperture=20 prim=1`
  gives 597 of 707 samples differing, worst |diff| 0.633.
- Judgement: this is a bug in the C, not a porting error, and the Rust gives the intended
  answer. **Do not "fix" Rust toward C.** It still matters because (a) the canonical demo
  pipeline produces visibly different seismograms from the reference, and (b) no stable golden
  is possible for large inputs: `make_goldens.py` would record heap garbage that depends on input
  size and other open files. At audit time the Rust comments pointed elsewhere for the explanation.
  It should be recorded next to the goldens.
- Suggested action: add a comment to `goldens/gbbeam/cases.txt` and
  `goldens/triseis/cases.txt` explaining that the C reads `s[-1]`, so goldens must keep
  rayend input small (here < 94 records / 8 KB), since above that the C reads heap garbage.
  Pin the Rust behavior with a unit test instead of a golden: `makexzs(2, [0,4], [0,0], ng=7)`
  must return `0, 0.667, ..., 4`. Adding `infofile=` to the C run is not a reliable way to
  get a clean reference. It fixed the 94-ray case but not the 130-ray `atten=1` case.

## M2 (medium): negative grid sizes panic (or silently write data) instead of SU's malloc error, in tri2uni and elatriuni

- tri2uni. Rust `programs/tri2uni/src/main.rs:98`
  (`vec![0.0f32; (n1 as usize).wrapping_mul(n2 as usize)]`). C `tri2uni.c:72` (`ealloc2float(n1,n2)`).
- elatriuni. Rust `programs/elatriuni/src/main.rs:109-110` (same `wrapping_mul` pattern).
  C `Trielas/main/elatriuni.c:136-145` (`ealloc2float(nz,nx)`).
- What differs: C's `ealloc2float` fails in `malloc` for any negative dimension. It prints
  `<prog>: ealloc.c: malloc failed (Cannot allocate memory)` and exits 1. Rust:
  - one negative size: `capacity overflow` panic, exit 101, Rust backtrace hint on stderr;
  - **tri2uni with both negative (`n1=-1 n2=-1`)**: the product wraps to 1. Rust exits **0 and
    writes 4 bytes** (one 0.0 float). That is silently wrong output where C errors.
- Reproduction (any tri model on stdin, e.g. `model.data` from M1; `emodel.bin` = the DisplModel elamodel output):
  ```
  tri2uni <model.data n1=5 n2=-1     # C: "tri2uni: ealloc.c: malloc failed (Cannot allocate memory)" rc=1 | Rust: panicked ... capacity overflow, rc=101
  tri2uni <model.data n1=-1 n2=5     # same
  tri2uni <model.data n1=-1 n2=-1    # C: same malloc error, rc=1, 0 bytes | Rust: rc=0, 4 bytes on stdout
  elatriuni <emodel.bin nx=-1 nz=5   # C: "elatriuni: ealloc.c: malloc failed (Cannot allocate memory)" rc=1 | Rust: capacity overflow panic, rc=101
  elatriuni <emodel.bin nx=5 nz=-1   # same
  ```
- Suggested fix: emulate `alloc2`. C mallocs the row-pointer array (`n2` entries) first, then
  `n1*n2` floats, so error when `n2 < 0 || n1 < 0` (or when the product overflows) with
  `err!("ealloc.c: malloc failed (Cannot allocate memory)")`. This is the pattern gbbeam already
  uses at `programs/gbbeam/src/main.rs:224`. Checking only the wrapped product is not enough
  because (-1)*(-1)=1.
- Golden lines (the harness compares exit code and first stderr line):
  ```
  [stderr] python3 -c 'import sys,struct; sys.stdout.buffer.write(struct.pack("9f",*[1]*9))' | uni2tri n1=3 n2=3 verbose=0 | tri2uni n1=5 n2=-1
  [stderr] python3 -c 'import sys,struct; sys.stdout.buffer.write(struct.pack("9f",*[1]*9))' | uni2tri n1=3 n2=3 verbose=0 | tri2uni n1=-1 n2=-1
  [stderr] elamodel fill=0.5,0.5,2,1,2 > m && elatriuni < m nx=-1 nz=5
  ```

---

## L1 (low): elaps (and any psplot user) with an unknown `labelfont`/`titlefont`

- Rust `crates/oxwave-f-basic/src/lib.rs:466-470` (`font()` does `.expect("unknown font...")`).
  C `vendor/su-original/src/psplot/lib/basic.c:390-403`: `fontindex` returns `NFONTS` (38) when
  the name is not found, which indexes past the end of `FontMetrics[NFONTS]` (UB).
- Repro: `elaps <emodel.bin labelfont=Bogus`. C: SIGSEGV (rc 139), 0 bytes out. Rust: panic, rc 101,
  454 bytes of partial PostScript on stdout.
- Fix: raise SIGSEGV the same way the tri programs emulate C's null-model fault, or
  error cleanly. Exit code is the only thing a golden can pin:
  `[exact] elamodel fill=0.5,0.5,2,1,2 > m; elaps labelfont=Bogus < m`.

## L2 (low): gbbeam `nt<0`

- Rust `programs/gbbeam/src/main.rs:542` (`&syn[..nt as usize]` panics, rc 101).
  C `gbbeam.c:506` (`efwrite(syn,sizeof(float),nt,stdout)` with `nt<0` turns into a huge `size_t`).
  C dumps about 4 MB of process memory to stdout and then dies (rc 1 into a pipe, SIGBUS rc 138 into
  /dev/null).
- Repro: `gbbeam <rayends.data xg=0,4 zg=0,0 ng=5 nt=-3 krecord=1`.
- C behavior is UB either way, so matching it is not desirable. A clean `err!` would be better than a panic. No golden is possible.

---

## Checked and found faithful (negative results)

- License header: every `.rs` in the 12 programs and in oxwave-lib-tri, oxwave-lib-trielas,
  f-lincoeff, cubic-interp, f-cubicspline, f-basic, f-psaxesbox, sort,
  inverse-interp, linear-interp and pfa-fft starts with `tools/cwp_header.txt` verbatim.
- Test-specific code: none found. No `Command::new`, no references to goldens or vendor paths.
  The `raise(11)` sites (elaray, elatriuni, tri2uni, triray, trimodel, normray, uni2tri)
  are documented emulations of the C null-model or overflow SIGSEGV. Other C segfaults on
  real models (e.g. `triray ... fangle=-170 langle=170 refseq=2,0 ...`, `normray nrays=2 nangle=3`
  on the syncline model) happen naturally in the mechanically translated unsafe code, at the
  same point as in C.
- Byte-identical or within 1e-7, on heterogeneous models:
  - triray: multi-layer syncline with refseq reflect/transmit/stop, krecord, prim, infofile,
    wavefile, rayfile, fresnelfile, outparfile.
  - normray: 42-edge NormalRays demo with caustic=0/1, surface=0/1, nonsurface=0/1, nangle>1,
    ashift, dangle, wavefile, infofile, fresnel, refseq, krecord, prim.
  - triseis: ns>1, ds/fs/dg/fg/dgds, kreflect, krecord=3, reftrans=1, atten=1/2, bw, lscale,
    aperture, prim.
  - gbbeam: < 94 rays, including atten=1/2, reftrans, bw, aperture, ang/amp.
  - uni2tri: method=1/2/3, tol, errmax, verbose=2/efile, ifile; tri2uni round trips.
  - elamodel, elaray, elaps, elatriuni, elacheck on the Trielas demos:
    - elaray mode 0/1/2 on isotropic and anisotropic (e, d, g, azimuth) layers, with reftrans=1,
      mode-converting refseq (2/3), fdip/fazi, tw, infofile;
    - elaps p=0..3, style=normal, explicit axes/ticks/grids/colors/fonts;
    - all ten elatriuni output grids;
    - elacheck interactive probes.
- The gbbeam `makesyn` float/double promotion order was checked line by line against C and matches.
- Inherited C behavior (both sides behave the same, not a divergence): `elacheck` loops forever on
  stdin EOF; `elaps dxnum<0` / `nxtic=0` hang; `trimodel maxangle=0.5` on a curved edge hangs;
  `elamodel` calls `checkpars()` before reading `fill=`, so `fill` is reported as unknown.

## Overall judgement

This group of programs is a high-fidelity port. Most of the tri/Trielas code is a mechanical, pointer-level
translation (raw `*mut` model graphs, the C float/double promotion kept, the same crash
points). On realistic heterogeneous models, well beyond the trivial geometries in the goldens,
every program reproduced the original outputs byte-for-byte or to float rounding: ray paths,
rayends, info and Fresnel files, wavefronts, synthetics, resampled grids, stiffness grids and
PostScript. I found no wrong-number bug in the Rust and no test-specific code. The one large divergence
on a stock demo (M1) comes from undefined behavior in the original gbbeam/triseis. There the
Rust's deliberate zero guard gives the physically correct answer, but it should be documented
so nobody chases it or "fixes" it. The remaining issues are robustness on invalid input:
Rust panics with exit 101 instead of SU's `ealloc` error, and in one case (tri2uni with both
sizes negative) silently writes a float. These are cheap to fix and worth one `[stderr]`
golden each.
