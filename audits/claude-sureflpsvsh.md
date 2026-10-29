# sureflpsvsh: diagnosis of the 8 failing golden cases

Auditor: Claude, 2026-10-29. The repo was not modified except for this file.
Workdir: `/private/tmp/claude-501/refl_audit` (builds, shims, and all run directories are kept there).

## Verdict

All 8 failing cases are **class (a)**: the C original has undefined behaviour (heap out-of-bounds
reads, and in one case an out-of-bounds **write**), and that UB decides the golden. None of them is a
Rust bug. The Rust panics (exit 101) come from deliberate guard asserts in the library API.
Those asserts reject exactly the inputs that put C out of bounds:
`crates/oxwave-lib-refl/src/interfaces.rs:56` (lsource >= 2), `:61` (1 <= lobs <= nlayers) and
`:401` (SH needs nx == ntc).

The earlier claim that "the C original writes all-NaN output and SH runs segfault, so it is not a usable
oracle" is **true of the current case set but false of the program**:
* All SU-trace goldens in the suite are NaN (or all zero in 35 and 71). This comes from the suite's
  receiver choice, not from SU in general. Every case uses `lobs=1` (or the UB value `lobs=0`). For a
  receiver at `lobs=1`, psvreflect never fills its free-surface arrays (details below), so C computes
  from **uninitialized heap memory**: zeros give 0/0 = NaN, and random garbage gives finite garbage.
* With a receiver below the top layer (`lobs>=2`, `2<=lsource`), the C writes finite, non-zero
  seismograms. They are ASan/UBSan clean, stable across runs, and robust to heap contents once
  one known uninit read is accounted for (see the trace-0 note below). On these cases the current Rust
  binary is **bit-exact** to the original for every sample, header and stderr byte.
* SH segfaults only when `nx != tsec/dt`, because of the argument swap at compute.c:230. With `nx == ntc`,
  SH runs and gives finite output. However, that output depends on one never-written element
  (shreflect.c:647), so I do not recommend SH numeric goldens.

So the program **can pass meaningful goldens**. Exclude the 8 UB cases and add the finite cases below,
which the Rust already passes bit-exactly.

## Method

* ASan+UBSan copy of the C: `cc -O0 -g -std=c90 -m64 -ansi -ffp-contract=off -DCWP_LITTLE_ENDIAN -DSUXDR
  -fsanitize=address,undefined` on `Refl/main/sureflpsvsh.c` + `Refl/lib/{compute,convolve,auxil,psvreflect,shreflect,synthetic}.c`,
  linked with `-lsu -lpar -lcwp -lm` (`$W/build.sh`). I also built a `-fsanitize-recover=address` copy to
  enumerate every OOB site, a plain copy with SU's flags (bit-identical to `vendor/su-build/bin/sureflpsvsh` on clean inputs), and
  instrumented variants:
  * `pad`: model arrays get guard words with a chosen value before them.
  * `rnd`: alloc1complex in one file returns zeroed, random or NaN memory.
  * `fix`: `refw` is zeroed, and/or shreflect.c:647 has the loop bound fixed.
* I ran all 73 manifest cases with orig, ASan, plain and Rust (a copy of `target/release/sureflpsvsh`), in a fresh dir each.

## Per-case classification of the failing cases

| k | key inputs | golden | ASan (first site / all sites) | class |
|---|---|---|---|---|
| 1 (inventory) | lsource=0 (default), lobs=0 | rc 0 | auxil.c:258 `cl[lmn]` with lmn=-1 (source_receiver_type); also auxil.c:262/264, compute_slowness auxil.c:338-341 (`cl/t[lsource-2]` = index -2), compute_prs auxil.c:498-504, compute_alpha_betha auxil.c:709-710 | (a) |
| 2 (inventory) | lobs=0 **and** SH nx=2 != ntc=64 | rc -11 | same lobs=0 reads, then SH overrun (see 72); the C segfaults | (a) |
| 68 | = k1 + `cat p.su r.su z.su` | rc 0, all-NaN traces | as k1. **Pad test:** guard value 1e30 before cl/ct/... changes the output from NaN to finite zeros, so the golden depends on the OOB value | (a) |
| 69 | lsource=0, lobs=1 | rc 0, NaN | auxil.c:262 `ct[jkl]` jkl=-1; compute_slowness auxil.c:338-341; compute_sigmas auxil.c:764-768; compute_p_aux_arrays auxil.c:184. Pad 1e30 changes output | (a) |
| 70 | lsource=2, lobs=0 | rc 0, NaN | auxil.c:258/278; compute_prs auxil.c:498-504; psvreflect.c:1412-1422 (`gl/gt[ijk2]`, ijk2=(lrec-1)*block_size < 0) | (a) |
| 71 | lsource=1 | rc 0, all-zero traces | compute_slowness auxil.c:344-346 (`ct/t[lsource-2]` = index -1), plus an **OOB WRITE** of 8 bytes before an 800-byte heap block at psvreflect.c:1646 (heap corruption) | (a) |
| 72 | SH, nx=2, ntc=tsec/dt=64 | rc -11 | compute.c:230 calls `compute_synthetics(verbose, nt, nx, ntc, ...)` but the prototype order is `(verbose, nt, ntc, nx, ...)`. The trace loop runs 64 times over a 2-row wavefield and refw is sized npfa(2): synthetic.c:352 / construct_tx_trace synthetic.c:435 overflow, then SIGSEGV | (a) |
| 73 | SH nx=ntc=64 (well-defined swap) but lobs=0 | rc 0, NaN | same lobs=0 reads as k1. Pad 1e30 changes output | (a) |

None of the 8 is (b), and none is (c) on its own. For 70 and 71 my pad values did not change the bytes,
because the OOB reads land in other heap blocks or are masked by the `p2w` cap. Still, the reads, and
for 71 the write, are UB on any input, so the golden is not defined behaviour.

Rust side (RUST_BACKTRACE=1): the panics are `oxwave_lib_refl::interfaces::compute_reflectivities`, which calls
`model()` at interfaces.rs:56 / :61 (k 1, 2, 68, 69, 70, 71, 73) and the SH guard at interfaces.rs:401 (k 72).
These are intentional and correct guards on the domain where C is defined. **No Rust fix is required.**
An optional nicety is to turn these asserts into an `err(...)` message instead of a panic. That is not
golden-testable, since C has no defined behaviour there.

### Exact lines to add to goldens/sureflpsvsh/EXCLUDE.txt

These are the bare commands as they appear in manifest.json. EXCLUDE matching is done after tag
stripping, so the `[stderr] ` prefix from cases.txt must **not** be included. k1 and k2 come from
`inventory/w3-progs-1.json`, so they can only be removed through EXCLUDE.txt.

```
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 outf=info nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lobs=0 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wfp=p.su wfr=r.su wfz=z.su
sureflpsvsh m0=1 stype=1 delta=45 lambda=90 phis=0 phi=30 wtype=2 win=1 verbose=2 outf=info nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lobs=0 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wft=t.su
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 outf=info nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lobs=0 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lsource=0 lobs=1 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lsource=2 lobs=0 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lsource=1 lobs=1 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lsource=2 lobs=1 nor=1 nx=2 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wtype=2 wft=t.su
sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 verbose=2 nlayers=2 cl=2,3 ct=1,1.5 ql=100,100 qt=50,50 rho=2,2.2 t=1,1 lobs=0 nor=1 nx=64 nt=64 tsec=0.256 dt=0.004 nw=16 nf=16 np=32 wtype=2 wft=t.su && cat t.su
```

Suggested header (same rule style as wpc1comp2):
`# The original reads (k71: also writes) outside the heap model/work arrays when
lsource<2 or lobs<1 (auxil.c:258-264, 338-346; psvreflect.c:1412-1422, 1646), and in SH when
nx != tsec/dt (compute.c:230 argument swap -> synthetic.c:352/435). ASan: heap-buffer-overflow.`
After excluding them, also delete the 6 matching lines from cases.txt (lines 94, 97-99, 102, 104 of
the current file); they only re-test UB.

## Are the passing goldens meaningful?

**No passing golden compares a finite seismogram.** Of the SU-output cases (k 19, 21-43, 45, 56, 57, 59, 64, 66,
67), every sample is NaN, except k35 (45 samples, all 0.0). The cause is `lobs=1` throughout:

* psvreflect.c:1053 fills the free-surface arrays (`rurf*/turf*/rdfr*/tdfr*`, allocated uninitialized at
  psvreflect.c:180-187) only `if (il==lobs[ijk]-2)`. The layer loop is `il=0..nlayers-2` (psvreflect.c:302),
  so for lobs=1 they are **never written**. The "receiver above source" branch then reads them
  (psvreflect.c:1342-1347 `det = turf11*turf22 - turf21*turf12`, with rtb = turf/det).
* With zeroed (fresh) memory, and with any uniform fill byte, det = 0, so 0/0 gives NaN everywhere. That
  is what the golden contains. With the `rnd` build (random garbage in those arrays), case 66 becomes **finite**
  (max 0.87). So these goldens record indeterminate memory. They keep the same bytes under the zero-init
  "one-line fix" rule, so they need not be excluded, but they test control flow, headers, stderr and the
  info file, **not numerics**.
* The text goldens (k 20, 44, 58: `cat info`/`report`) are meaningful for the processing-info printout.
  Error-path goldens (rc 1) are meaningful.
* k16-18 and k65 (rc -11) are a deterministic NULL `FILE*` deref: sureflpsvsh.c:612 `fprintf(outfp, ...)`
  with outfp=NULL when verbose<2. This is also UB, but a reliable crash; the Rust emulates it. I would keep them
  (a judgement call). Practical consequence: **the program cannot run at all unless verbose=2 or 3.**
* k56: ASan heap-buffer-overflow READ at auxil.c:893 (`intlayers[ii]` with ii == nlint after the last
  interpolation, allocated at sureflpsvsh.c:538). The golden is NaN (lobs=1) anyway and unchanged under a bounds fix. Keep it.

## Other latent C defects found (relevant to choosing new goldens)

1. **synthetic.c:215** `refw=alloc1complex(ntfft)` is never zeroed before the first
   `construct_tx_trace`. Only bins 1..nw1 are written (synthetic.c:322-330), and the array is zeroed only after
   each inverse FFT (synthetic.c ~440). So **trace 0 of every component** includes uninit bins 0 and
   nw1+1..ntfft-1. The orig binary reliably sees zeros here (fresh pages), and the Rust zero-initialises. Using
   the "fix" build (memset after line 215) with ASan malloc_fill_byte 0x07/0xbe/0xff, the outputs of all
   proposed PSV cases are byte-identical to the orig golden. Under the EXCLUDE.txt rule ("excluded only if
   output changes when that one line is fixed"), they therefore stay.
2. **shreflect.c:647** `for (ip=0; ip<block_size-1; ip++) texp[ip]=...`, but the integration reads
   `texp[jj]`, jj = block_size-1 (shreflect.c:659, 664-666). PSV uses `ip<block_size` (psvreflect.c:1964).
   The SH output depends on that uninit element: zero gives the orig bytes, and fixing the bound changes the output
   by about 2e-3 of the trace max. **SH numeric goldens are therefore not robust.**
3. **wavelet_type=3**: convolve.c:44 calls `ricker2_wavelet(nt/2, ...)` into `wavelet[nt]`, which writes only
   2*(nt/2)-1 samples (cwp/lib/waveforms.c:76-82). wavelet[nt-1] (odd nt: last 2) is uninit. Avoid
   wavelet_type=3 in numeric goldens.
4. libcwp `convolve_cwp` reads `x[jlow-1]` (cwp/lib/convolution.c:181, value unused). In my -O0 plain
   rebuild this intermittently SIGBUSes (about 10% of runs) when the scratch array sits at a page start. The vendored
   orig binary did not crash in 60 runs. This is harmless for goldens, but it explains flaky crashes in custom C builds.

## Proposed new golden cases (verified)

All of these pass every check below. Orig rc 0; stable over 3 runs; 100% finite, mostly non-zero samples.
ASan+UBSan clean. Byte-identical between orig and the refw-zeroed ASan build under fill bytes 0xbe and 0xff,
so there is no dependence on indeterminate memory. **Current Rust binary bit-exact on stdout and stderr.**
I recommend the `[stderr]` tag, as in the existing cases.

```
[stderr] sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 h1=0 h2=0 verbose=2 nlayers=4 cl=3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.9 ql=50,50,50,1000 qt=30,30,30,600 rho=2.65,2.65,2.65,2.8 t=0.001,0.343,0.655,2 lsource=3 lobs=2 nor=1 nx=4 bx=0.025 dx=0.25 tsec=1.024 dt=0.004 nt=256 nw=64 nf=64 p2w=1 pw1=0 pw2=0 pw3=0 pw4=0 decay=5 red_vel=3.5 tlag=0.2 wfield=2 wavelet_type=2 fpeak=25 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
[stderr] sureflpsvsh m0=1 stype=2 m1=1 m2=0.2 m3=-0.3 h1=0 h2=0 verbose=2 nlayers=3 cl=3,4,5 ct=1.5,2,2.5 ql=100,120,150 qt=50,60,80 rho=2,2.2,2.5 t=0.01,0.2,1 lsource=3 lobs=2 nor=1 nx=3 bx=0.1 dx=0.2 tsec=0.512 dt=0.004 nt=128 nw=32 nf=32 p2w=0.8 pw1=0 pw2=0 pw3=0 pw4=0 decay=10 red_vel=0 tlag=0.1 wfield=3 wavelet_type=4 fpeak=20 vsp=0 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
[stderr] sureflpsvsh m0=1 stype=2 m1=1 m2=0 m3=1 h1=0 h2=0 verbose=2 nlayers=4 cl=3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.9 ql=50,50,50,1000 qt=30,30,30,600 rho=2.65,2.65,2.65,2.8 t=0.001,0.343,0.655,2 lsource=3 lobs=2 nor=1 nx=4 bx=0.025 dx=0.25 tsec=1.024 dt=0.004 nt=256 nw=64 nf=64 p2w=1 decay=5 red_vel=3.5 tlag=0.2 wfield=2 wavelet_type=2 fpeak=25 filters_type=1,2 filters_phase=0,0 dbpo=12,12 f1=5,40 f2=10,50 wfp=p.su && cat p.su
[stderr] sureflpsvsh stype=1 delta=60 lambda=30 phis=10 phi=45 m0=1 h1=0 h2=0 verbose=2 nlayers=4 cl=3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.9 ql=50,50,50,1000 qt=30,30,30,600 rho=2.65,2.65,2.65,2.8 t=0.001,0.343,0.655,2 lsource=3 lobs=2 nor=1 nx=4 bx=0.025 dx=0.25 tsec=1.024 dt=0.004 nt=256 nw=64 nf=64 p2w=1 pw1=0 pw2=0 pw3=0 pw4=0 decay=5 red_vel=3.5 tlag=0.2 wfield=1 wavelet_type=1 int_type=2 win=1 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
[stderr] sureflpsvsh stype=2 m1=1 m2=0 m3=1 m0=1 h1=0 h2=0 verbose=2 nlayers=4 cl=3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.9 ql=50,50,50,1000 qt=30,30,30,600 rho=2.65,2.65,2.65,2.8 t=0.001,0.343,0.655,2 lsource=3 lobs=2 nor=1 nx=4 bx=0.025 dx=0.25 tsec=1.024 dt=0.004 nt=256 nw=64 nf=64 p2w=1 pw1=0 pw2=0 pw3=0 pw4=0 decay=5 red_vel=3.5 tlag=0.2 wfield=2 wavelet_type=2 vsp=1 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
[stderr] sureflpsvsh stype=2 m1=1 m2=0 m3=1 m0=1 h1=0 h2=0 verbose=2 nlayers=4 cl=3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.9 ql=50,50,50,1000 qt=30,30,30,600 rho=2.65,2.65,2.65,2.8 t=0.001,0.343,0.655,2 lsource=3 lobs=2 nor=1 nx=4 bx=0.025 dx=0.25 tsec=1.024 dt=0.004 nt=256 nw=64 nf=64 p2w=1 pw1=0 pw2=0 pw3=0 pw4=0 decay=5 red_vel=3.5 tlag=0.2 wfield=2 wavelet_type=2 flt=1 wfz=z.su && cat z.su
[stderr] sureflpsvsh stype=2 m1=1 m2=0.5 m3=1 nlayers=5 cl=3.5,3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.2,2.9 ql=50,50,50,50,1000 qt=30,30,30,30,600 rho=2.65,2.65,2.65,2.65,2.8 t=0.001,0.001,0.343,0.655,2 lsource=4 lobs=2,3 nor=2 m0=1 h1=0 h2=0 verbose=2 nx=3 bx=0.025 dx=0.25 tsec=1.024 dt=0.004 nt=256 nw=64 nf=64 p2w=1 pw1=0 pw2=0 pw3=0 pw4=0 decay=5 red_vel=3.5 tlag=0.2 wfield=2 wavelet_type=2 wfp=p.su && cat p.su
```

What they cover: (1) PSV, stype=2, velocity, ricker1, reducing velocity and time lag, all 3 components;
(2) 3-layer model, red_vel=0 (computed), acceleration, akb wavelet; (3) output filters (lo+hi cut);
(4) fault-plane source (stype=1), displacement, spike, Filon integration, frequency window;
(5) vsp=1; (6) earth flattening; (7) nor=2 receivers.
Results (samples / max |x|): 3072/1840, 1152/7.7e4, 1024/141, 3072/2.56, 3072/13.1, 1024/656, 768/1988.

Optional SH case. It is finite and the Rust is bit-exact, but it depends on the shreflect.c:647 uninit element,
so add it only if you accept zero-init as the reference:
```
[stderr] sureflpsvsh m0=1 stype=1 delta=90 lambda=0 phis=0 phi=30 h1=0 h2=0 wtype=2 verbose=2 nlayers=4 cl=3.5,3.5,3.5,5.9 ct=2.2,2.2,2.2,2.9 ql=50,50,50,1000 qt=30,30,30,600 rho=2.65,2.65,2.65,2.8 t=0.001,0.343,0.655,2 lsource=3 lobs=2 nor=1 nx=64 bx=0.025 dx=0.05 tsec=0.256 dt=0.004 nt=64 nw=32 nf=32 p2w=1 pw1=0 pw2=0 pw3=0 pw4=0 decay=5 red_vel=3.5 tlag=0.05 wavelet_type=2 fpeak=25 wft=t.su && cat t.su
```

Cheap extra win: the existing parameter-coverage block (cases.txt lines 14-82) uses
`cl=3,4,5 ... lobs=1 lsource=2`. Changing it to `lobs=2 lsource=3` makes it produce finite output: I checked
the base with wfp; it gives 45 finite samples, max 72, and is ASan-clean. Each derived case would need its own
ASan check before being kept.

## Recommendation

1. Add the 8 lines above to `goldens/sureflpsvsh/EXCLUDE.txt`, and drop cases.txt lines 94,
   97, 98, 99, 102, 104, since they duplicate the UB cases.
2. Add the 7 verified finite PSV cases. Optionally re-base the coverage block to `lobs=2 lsource=3`, and keep
   only variants that are ASan-clean and robust to heap fill.
3. Rust: no correctness fix needed. The guards at interfaces.rs:56/61/401 are correct. Keep zero-initialised
   scratch for refw, texp and the psv free-surface arrays: that is what reproduces the orig bytes.
4. Expected result: 66 passing + 7 new finite cases, all passing now. This is a meaningful pass: the Rust is bit-exact to
   the C on real seismograms, including a 20 x 1024 x 3 demo-sized run (Refl demo Rpsv3 parameters).
