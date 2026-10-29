# OxWave golden-case coverage audit

Snapshot: 2026-10-29 00:15 CDT, 398 program crates. The repo was being changed while this ran: `supolar` and `sureflpsvsh` were still being ported. Everything was produced in `/private/tmp/claude-501/oxcov/`; nothing in the repo was changed except this report and its JSON:

- `coverage.py`: the metric script; writes `metrics_final.json` and `table_final.txt`.
- `runcases.py`: runs proposed cases against the C and Rust binaries the same way `check_program.py` does.
- `cases/<prog>.txt`: the proposed cases.
- `results/<prog>.txt`: the run logs.

## Summary

- **3 divergences found** (details at the end):
  1. **suprofcsv: the Rust port is an unimplemented stub.** It ends with `err!("suprofcsv: profile generation is not yet implemented")`. It passes only because none of its 12 golden cases ever produces a profile: every one exits 1 in C.
  2. **surescsv: a sort-order difference.** The Rust uses a stable `sort_by` where the C uses `qsort`. That changes the order of float summation, and the `%.12g` CSV text differs in the 8th–9th digit. The relative error is about 1e-7, but the text compare in the harness fails.
  3. **sufctanismod: the C reads uninitialised memory** in the flux-corrected-transport step (dofct=1). This makes long `dofct=1 isurf=1` runs differ by 0.8% relative. It is not a port bug: the port zero-fills, and C run with `MallocPreScribble=1` matches the Rust to 4e-39. Four existing goldens rest on this undefined behaviour.
- 35 programs got proposals: the top 25 plus 10 that were in the top 25 under earlier versions of the metric. **141 new cases were run against both binaries.** Apart from the cases above, everything matched within 1e-5 (sureflpsvsh was run in C only).
- **197 golden cases in 84 programs check nothing but the exit code.** Their stdout goes to a file (`> out.su`) that is never `cat`ed. 169 of them come from the inventory `test_cases`. Examples:
  - all 4 suwaveform inventory cases;
  - 8 segyread cases;
  - both raydata inventory cases, where `elamodel … && elaray … && raydata` never even runs raydata, because elamodel and elaray exit 1 on success.

## Method

### Per-program metrics (`coverage.py`)

| Metric | How it is computed |
|---|---|
| **cases** | `oxharness.collect_cases(prog)`: @selfdoc + inventory test_cases + cases.txt, after EXCLUDE.txt and de-duplication. These are the cases the harness actually runs. |
| **data / err / empty** | Uses each case's original-C outcome from `goldens/<prog>/manifest.json` and the `.out`/`.err` sizes. For the programs without a manifest at snapshot time (sufilter, sureflpsvsh) the C binary is run. |
| **params** | Every name read through `getpar*` / `getnpar*` / `MUSTGETPAR*` / `countparval` / `countparname` with a literal name in the C source. If there are none, sdoc `name=` rows (that also appear as a quoted literal in the C) plus inventory params. `verbose` and `par` are excluded. |
| **never named** | Params that appear in no case at all, error cases included. |
| **weak** (<2 values) | Params that take fewer than 2 distinct values across the data cases, so at most one setting of the branch they control is ever checked. File-name and temp-dir params are excluded, since they naturally take one value. |
| **string modes** | Literals compared with `STREQ`/`strcmp`/`strncmp` in the C (mode names such as `type=berlage` or `op=diff`), minus SU header key names. A mode is untested if it never appears in a data case. |
| **C lines** | `wc -l` of the C main file (inventory `path`). |

How a case is classified:
- **data:** the original succeeded (rc 0) with non-empty stdout, or wrote at least 200 bytes of stdout even with rc ≠ 0 (several Trielas programs always exit 1).
- **data, for programs that report on stderr:** a `[stderr]` case with rc 0 and at least 64 bytes of stderr.
- **err:** the original exited non-zero with no real output.
- **empty:** rc 0 and nothing compared but the exit code.

How parameters are matched to cases: each case line is split with a shell-aware tokenizer, so `surf="-10,0;40,0"` stays one word. A stage counts as the program's when the program's bare name appears anywhere in it, not only first: the csv suite runs behind `python3 -c '…pty…' prog …` or `script -q /dev/null prog …`. Parameters are the `name=value` words after the program name. Quoted command lines (`subprocess.call(" prog …")`) are parsed recursively, and `par=` files written inside the case are scanned too.

### Weakness score

```
score = W * log2(1 + Clines/50) * ( 0.5  * fw * sqrt(Nw)    # params with <2 values in data cases
                                  + 0.5  * fv * sqrt(Uv)    # string modes never used in a data case
                                  + 0.25 * fu * sqrt(U)     # params never named in any case
                                  + 1 / sqrt(1 + data) )    # depth: how many cases check real output
```

`f*` is the untested fraction of each set and `N*`/`U*` is its count. Why this form:

- **Fraction × √count:** a program with 20 params, half of them weak, scores higher than one with 2 params, half weak. Very large parameter lists (psplot, csv) still don't swamp everything.
- **The depth term** keeps big programs with few data cases visible. Parameter-name coverage saturates in well-built ports: nearly every program names every parameter in some case.
- **log(C size)** stands for the number of branches.
- **W (importance):**

| W | Programs |
|---|---|
| 3 | ~100 widely used processing programs: sufilter, sugain, sunmo, sustack, sumute, suxcor, sufft, suresamp, the Kirchhoff/phase-shift migrations and datumers, decon, velan, statics, segyread/segywrite, the main generators, … |
| 2 | Other su/main numerical families |
| 1.5 | Headers, par/main numerics, Trielas/tri/tetra/Refl, compression, test-pattern generators |
| 1 | Data conversion, well logs, supromax, cwp/main, epidemic/ODE toys |
| 0.75 | The csv geometry suite: one contributor's add-on, niche |
| 0.5 | PostScript plotting (their weak params are mostly cosmetic: colours, fonts, labels) |

Metric limitations, found while validating it against the forks' close reading:
- Programs whose result is an exit code (sucmp) or a file of tiny records (raydata) score as weak while actually being well tested. Both are in the top 25 and are flagged in their sections.
- A param counted as "weak" may still be well covered when its two settings are separate branches tested through other params.
- The table is a triage list. The per-program sections below are the actual judgement.

Where the main programs of interest rank:

| Program | Rank |
|---|---|
| sugain | #3 |
| sumute | #12 |
| sukdmig2d | #21 |
| sunmo | #30 |
| sustack | #79 |
| suresamp | #91 |
| sufilter | #112 |
| suxcor | #147 |
| sufft | #159 |
| sumigtk | #215 |

These are comparatively well tested. Every parameter is named, and most take several values in cases with real output.

### Running proposed cases (`runcases.py`)

1. Each case runs with `oxharness.run_case` in a fresh temp dir. The C side uses `vendor/su-build/bin` first on PATH. The Rust side puts a temp dir holding only a `<prog>` symlink to `target/release/<prog>` in front of it, like `check_program.py`.
2. The C side runs twice; a case is marked UNSTABLE if the two runs differ.
3. The comparison follows `check_program.py`:
   - exit code;
   - then stdout through `oxharness.compare_stdout` (SU traces at rtol 1e-5 of the trace max; text exact or `[ws]`; `[f32]`; bytes);
   - the first stderr line when C failed;
   - the whole stderr for `[stderr]` cases.

## Ranked table (top 40)

Columns:
- **never named:** params that appear in no case.
- **<2 values:** params that take fewer than 2 distinct values in data cases.
- **string modes:** STREQ literals never used in a data case.
- **data / err:** cases with real compared output / cases that are error exits.

Ranks 1, 4, 18 and 20 carry notes:
- **sureflpsvsh (#1):** still being ported. It has no goldens dir, no release binary, and only 2 inventory cases, both of which redirect to files.
- **suprofcsv (#4):** a stub (see Divergences).
- **sucmp (#18) and raydata (#20):** metric false positives, explained in their sections.

| # | program | score | W | C lines | cases | data | err | never named | <2 values in data cases | untested string modes | weakest params (sample) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | sureflpsvsh | 40.0 | 1.5 | 904 | 3 | 0 | 1 | 53/85 | 79/85 | 0/0 | bp, bx, cl, ct, dbpo, decay, delta, dt |
| 2 | sudatumk2dr | 20.2 | 3 | 1242 | 77 | 33 | 43 | 0/38 | 21/38 | 0/0 | angmax, antiali, ds, dxgo, dxi, dxso, dxt, dzt |
| 3 | sugain | 14.5 | 3 | 877 | 33 | 22 | 10 | 0/25 | 13/25 | 0/0 | agc, dt, gagc, jon, mark, maxbal, mbal, norm |
| 4 | suprofcsv | 13.8 | 0.75 | 1029 | 13 | 0 | 12 | 8/18 | 18/18 | 2/2 | chordf, chordi, dextra, ename, formxy, keyp, keyx, keyy |
| 5 | segywrite | 13.4 | 3 | 735 | 20 | 7 | 10 | 0/12 | 7/12 | 0/0 | buff, conv, endian, errmax, tape, trmax, trmin |
| 6 | sufctanismod | 13.3 | 1.5 | 2589 | 69 | 55 | 11 | 0/64 | 32/64 | 0/0 | aa00, cc00, daadx, daadz, dccdx, dccdz, deta0dx, deta0dz |
| 7 | suremel2dan | 12.5 | 1.5 | 2158 | 89 | 13 | 74 | 0/44 | 23/44 | 0/0 | abso, dtsnap, dx, dz, fmax, fx, fz, irx |
| 8 | unif2aniso | 11.0 | 1.5 | 764 | 14 | 7 | 4 | 0/45 | 26/45 | 0/0 | dddx, dddz, dedx, dedz, delta00, dgdx, dgdz, dqdx |
| 9 | suharlan | 10.8 | 3 | 1857 | 71 | 43 | 24 | 0/34 | 11/34 | 0/0 | depthref, dt, dx, interoff, ninterp, norm, offref, pmula |
| 10 | sugeomcsv | 9.8 | 0.75 | 3390 | 119 | 62 | 56 | 0/17 | 1/17 | 34/50 | scalco |
| 11 | sumatrix | 8.7 | 2 | 484 | 27 | 7 | 17 | 0/10 | 7/10 | 0/0 | maxreceivers, maxsources, maxtraces, rdivider, rkeyloc, sdivider, skeyloc |
| 12 | sumute | 8.2 | 3 | 345 | 29 | 15 | 13 | 0/14 | 7/14 | 0/0 | absolute, below, hmute, key, linvel, nmute, twindow |
| 13 | sukdmig3d | 8.1 | 3 | 1334 | 96 | 56 | 39 | 0/48 | 12/48 | 0/0 | angmax, dxgd, dxm, dxo, dygd, dyo, ft, fxgd |
| 14 | sutxtaper | 7.7 | 3 | 511 | 23 | 17 | 5 | 0/12 | 0/12 | 1/1 |  |
| 15 | sukdmdcr | 7.5 | 2 | 931 | 62 | 22 | 39 | 0/35 | 13/35 | 0/0 | ds, dxgo, dxso, dxt, dzt, fs, fxt, fzt |
| 16 | sugassmann | 7.5 | 1 | 760 | 17 | 8 | 2 | 0/23 | 17/23 | 0/0 | a, api, b, c, g, h, it1, kmin |
| 17 | sunomstreamers | 6.8 | 1.5 | 717 | 32 | 7 | 22 | 0/23 | 11/23 | 0/0 | channelkey, firstshot, gunkey, gunmod, idntkey, linelocs, offpi, shotkey |
| 18 | sucmp | 6.7 | 2 | 184 | 19 | 0 | 11 | 0/1 | 1/1 | 0/0 | limit |
| 19 | suvelan_nsel | 6.6 | 3 | 346 | 14 | 8 | 4 | 0/10 | 4/10 | 0/0 | dv, dx, fv, pwr |
| 20 | raydata | 6.5 | 1.5 | 257 | 54 | 5 | 48 | 0/9 | 8/9 | 0/0 | ascci, pol, px, pxvx, pz, t, vgx, vgz |
| 21 | sukdmig2d | 6.1 | 3 | 919 | 79 | 45 | 33 | 0/41 | 9/41 | 0/0 | ds, dxt, dzt, fs, fxt, fzt, ns, nxt |
| 22 | pswigp | 6.0 | 0.5 | 499 | 15 | 8 | 6 | 0/51 | 46/51 | 1/5 | axescolor, axeswidth, backcolor, bias, clip, curve, curvecolor, curvedash |
| 23 | sumedian | 6.0 | 3 | 488 | 36 | 15 | 20 | 0/12 | 4/12 | 0/0 | key, median, nshift, subtract |
| 24 | sudatumk2ds | 5.9 | 3 | 1166 | 70 | 27 | 41 | 0/39 | 7/39 | 0/0 | ds, dxi, dxt, exgo, fxgo, fxso, fzo |
| 25 | pscubecontour | 5.9 | 0.5 | 607 | 16 | 6 | 9 | 0/60 | 48/60 | 0/4 | angle, axescolor, cgray, d1num, d2num, d3num, dc, f1 |
| 26 | psgraph | 5.8 | 0.5 | 528 | 17 | 10 | 6 | 0/47 | 40/47 | 2/5 | axescolor, axeswidth, d1num, d2num, f1num, f2num, gridcolor, gridwidth |
| 27 | sutetraray | 5.6 | 1.5 | 3373 | 64 | 32 | 28 | 0/36 | 10/36 | 0/0 | dsz, dtwf, dxgd, dygd, fxgd, fygd, maxnsegs, nsy |
| 28 | suwindpoly | 5.6 | 2 | 741 | 24 | 8 | 12 | 0/11 | 4/11 | 0/0 | ilkey, xkey, xlkey, ykey |
| 29 | unif2ti2 | 5.5 | 1.5 | 688 | 18 | 10 | 6 | 1/41 | 14/41 | 0/0 | dddx, dddz, dedx, dedz, drdx, drdz, dvpdx, dvpdz |
| 30 | sunmo | 5.5 | 3 | 422 | 22 | 12 | 9 | 0/9 | 3/9 | 0/0 | invert, sscale, upward |
| 31 | subincsv | 5.3 | 0.75 | 1057 | 23 | 9 | 11 | 0/7 | 4/7 | 2/2 | check, offset, rpkey, spkey |
| 32 | sugeom | 5.3 | 1.5 | 736 | 18 | 9 | 8 | 0/7 | 4/7 | 0/0 | dt, rps, sps, xps |
| 33 | sumax | 5.3 | 3 | 458 | 35 | 26 | 6 | 0/8 | 2/8 | 2/9 | key1, key2 |
| 34 | sustackup | 5.1 | 2 | 432 | 29 | 12 | 16 | 0/11 | 5/11 | 0/0 | keyabs, keyabs2, keysign, keysign2, maxstack |
| 35 | sutoolcsv | 5.1 | 0.75 | 3358 | 219 | 148 | 70 | 0/18 | 0/18 | 23/53 |  |
| 36 | suazimuth | 4.8 | 1.5 | 370 | 15 | 9 | 5 | 0/10 | 6/10 | 0/0 | cmp, mxkey, mykey, offkey, offset, signedflag |
| 37 | makevel | 4.8 | 1.5 | 214 | 19 | 8 | 4 | 0/29 | 15/29 | 0/0 | dlens, dvdy, dy, exc, l1c, l2c, lambda, tlens |
| 38 | sunmo_a | 4.7 | 2 | 478 | 21 | 11 | 9 | 0/10 | 4/10 | 0/0 | cdp, invert, sscale, upward |
| 39 | pswigb | 4.7 | 0.5 | 513 | 14 | 8 | 5 | 0/50 | 38/50 | 0/4 | axescolor, axeswidth, bias, clip, d1, d1num, d2, d2num |
| 40 | suvelan_usel | 4.5 | 2 | 325 | 13 | 6 | 4 | 0/10 | 4/10 | 0/0 | dv, fv, nx, pwr |

## Proposed cases: top 25

Each program's section ends with its case file. Each case line has a comment line above it saying which branch it targets, and the lines are in `cases.txt` syntax. Unless a section says otherwise, every case succeeds in C with real output. The same lines (without comments) are in `audits/claude-coverage-cases.json`.

### #1 sureflpsvsh (score 40.0, 904 C lines, 3 cases / 0 data)

**Status:** still being ported at the time. There is no `goldens/sureflpsvsh` directory and no `target/release/sureflpsvsh` binary. The only cases are 2 inventory cases, and both write `wf*.su` files without `cat`, so they check only the exit code.

**The original C is broken in this build:**
- Every PSV run writes all-NaN samples; stderr shows `ntfft=2` for nt=64 and nt=128.
- SH (`wtype=2`) segfaults with rc 139, including the second inventory case.
- Leaving out `verbose=2 outf=info` also segfaults.
- A likely cause is in `Refl/lib/synthetic.c`: the `construct_tx_trace` prototype at line 24 lacks the `ntc` argument.

**What the cases below can and cannot do:** they succeed in C (rc 0), but the samples are NaN. So they pin only headers, trace counts and the NaN pattern, and are run in C only. Fix or EXCLUDE the C before using it as an oracle.

```
# WARNING: in this C build every PSV run (including the inventory case) writes all-NaN samples (stderr shows ntfft=2 for nt=64/128), and SH (wtype=2) segfaults (rc 139). These cases only pin headers/trace counts/NaN pattern until the C behaviour is understood.
# PSV, moment tensor given (stype=2), particle velocity, 3-layer model with receivers at the surface; cat the radial/vertical seismograms (the inventory cases never cat their outputs, so they only check the exit code)
sureflpsvsh verbose=2 outf=info m0=1 stype=2 m1=1 m2=0 m3=1 wtype=1 wfield=2 nlayers=3 cl=2,3,4 ct=1,1.6,2.2 ql=100,200,300 qt=50,100,150 rho=2,2.2,2.4 t=0.5,1,1 lsource=2 lobs=0 nor=1 bx=0.1 dx=0.2 nx=3 tsec=0.512 dt=0.004 nt=128 nw=32 nf=32 np=64 wavelet_type=2 fpeak=10 wfp=p.su wfr=r.su wfz=z.su && cat p.su r.su z.su
# PSV from fault-plane mechanism (stype=1: dip/rake/strike/azimuth), displacement output (wfield=1), Hanning frequency window, akb wavelet. (SH, wtype=2, segfaults in this C build - rc 139 - like the 2nd inventory case, so it cannot carry a data golden)
sureflpsvsh verbose=2 outf=info m0=1 stype=1 delta=45 lambda=90 phis=10 phi=30 wtype=1 wfield=1 win=1 wavelet_type=3 fpeak=12 nlayers=3 cl=2,3,4 ct=1,1.6,2.2 ql=100,200,300 qt=50,100,150 rho=2,2.2,2.4 t=0.5,1,1 lsource=2 lobs=0 nor=1 bx=0.1 dx=0.2 nx=3 tsec=0.512 dt=0.004 nt=128 nw=32 nf=32 np=64 wfp=p.su wfz=z.su && cat p.su z.su
# Filon slowness integration (int_type=2), earth-flattening (flt=1), acceleration output, two output filters (zero-phase high-cut + min-phase low-cut) and a time lag / reducing velocity
sureflpsvsh verbose=2 outf=info m0=1 stype=2 m1=1 m2=0.5 m3=-1 wtype=1 wfield=3 int_type=2 fs=0.09 flt=1 nlayers=3 cl=2,3,4 ct=1,1.6,2.2 ql=100,200,300 qt=50,100,150 rho=2,2.2,2.4 t=0.5,1,1 lsource=1 lobs=0 nor=1 bx=0.1 dx=0.2 nx=3 tsec=0.512 dt=0.004 nt=128 nw=32 nf=32 tlag=0.05 red_vel=3 nfilters=2 filters_phase=0,1 filters_type=1,2 dbpo=12,12 f1=30,3 f2=0,0 wavelet_type=2 fpeak=10 wfz=z.su && cat z.su
# VSP geometry (vsp=1) with two receiver depths (nor=2, lobs=1,2), ct=-1/qt=-1 defaults (ct=cl/sqrt3, qt=ql/2), spike wavelet
sureflpsvsh verbose=2 outf=info m0=1 stype=2 m1=1 m2=0 m3=1 wtype=1 vsp=1 nlayers=3 cl=2,3,4 ct=-1,-1,-1 ql=100,200,300 qt=-1,-1,-1 rho=2,2.2,2.4 t=0.5,1,1 lsource=1 lobs=1,2 nor=2 bx=0.1 dx=0.2 nx=3 tsec=0.512 dt=0.004 nt=128 nw=32 nf=32 np=64 wavelet_type=1 wfz=z.su && cat z.su
```

### #2 sudatumk2dr (score 20.2, 1242 C lines, 77 cases / 33 data)

**Gaps:**
- All cases use one traveltime-table geometry.
- No case has sources between table nodes.
- No case skips traces while keeping others.
- No case combines topography with `dvz`.

**C quirk pinned by case 2:** when the first trace of a gather is skipped, the C writes the kept headers with all-zero samples.

All 4 cases: SAME.

```
# Off-grid source positions (sx=3,17 between tt shots): exercises the res interpolation of ttab/tb/pb/sigb/angb between shots and the new-gather fxin=gx branch, on a 5x6x5 table.
[stderr] suspike nt=65 ntr=6 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=3,3 b=0,5 c=14,14 j=3,3 > in.su; python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",max(5,math.hypot(x*10-s*10,z*8))/1500) for s in range(5) for x in range(5) for z in range(6)))'; sudatumk2dr infile=in.su outfile=out.su ttfile=tt nxt=5 fxt=0 dxt=10 nzt=6 fzt=0 dzt=8 ns=5 fs=0 ds=10 nxi=9 fxi=0 dxi=5 sgn=1 surf="-10,0;50,0" surf="-10,12;50,12" && cat out.su
# Sources outside [fs,es] are skipped (first gather sx=-20) and the next gather (sx=10) is read: C outputs 3 headers with all-zero samples here (gather index/fxin quirk) -- pins that quirk.
[stderr] suspike nt=65 ntr=6 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=-20,-20 b=0,5 c=30.05,30 j=3,3 > in.su; python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",max(5,math.hypot(x*10-s*10,z*8))/1500) for s in range(5) for x in range(5) for z in range(6)))'; sudatumk2dr infile=in.su outfile=out.su ttfile=tt nxt=5 fxt=0 dxt=10 nzt=6 fzt=0 dzt=8 ns=5 fs=0 ds=10 nxi=9 fxi=0 dxi=5 sgn=1 surf="-10,0;50,0" surf="-10,12;50,12" && cat out.su
# Topographic source surface combined with a vertical velocity gradient (v00 = v0 + sz*dvz with sz != 0) and a small aperture/angle so dat2d's nxf/nxe edge clipping is hit.
[stderr] suspike nt=65 ntr=6 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=10,10 b=0,5 c=10,10 j=3,3 > in.su; python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",max(5,math.hypot(x*10-s*10,z*8))/1500) for s in range(5) for x in range(5) for z in range(6)))'; sudatumk2dr infile=in.su outfile=out.su ttfile=tt nxt=5 fxt=0 dxt=10 nzt=6 fzt=0 dzt=8 ns=5 fs=0 ds=10 nxi=9 fxi=0 dxi=5 sgn=1 surf="-10,2;10,5;25,1;50,3" surf="-10,14;15,18;50,13" v0=1700 dvz=0.8 aperx=6 angmax=35 && cat out.su
# Noisy multi-sample input instead of a single spike, non-default output grid (dzo, nzo) and fmax below Nyquist.
[stderr] suplane nt=65 ntr=6 | suaddnoise seed=11 sn=4 | sushw key=sx,gx a=0,0 b=0,5 c=20,20 j=3,3 > in.su; python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",max(5,math.hypot(x*10-s*10,z*8))/1500) for s in range(5) for x in range(5) for z in range(6)))'; sudatumk2dr infile=in.su outfile=out.su ttfile=tt nxt=5 fxt=0 dxt=10 nzt=6 fzt=0 dzt=8 ns=5 fs=0 ds=10 nxi=9 fxi=0 dxi=5 sgn=1 surf="-10,0;50,0" surf="-10,12;50,12" nzo=13 dzo=3 fmax=50 && cat out.su
```

### #3 sugain (score 14.5, 877 C lines, 33 cases / 22 data)

Every existing data case uses sparse `suplane` data, so these use dense noisy data instead. Branches newly exercised:
- agc and gagc on dense data;
- `panel=1` for tpow, qclip, agc and qbal;
- `gpow=0.5` (the sqrt branch) and the general `pow` branch;
- `pbal`, `mbal`, `maxbal` and `trap`;
- `mark` with `jon=1`;
- `vred` with negative offsets.

All 7 cases: SAME. Three further probes also matched and are not listed: panel+vred+pbal, gpow=-1 with clip/pclip/nclip/scale, and epow=-2 with bias/mbal/maxbal.

```
# agc on noisy (non-sparse) data, with norm rescaling (scale/=norm); all existing agc cases use sparse suplane data
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sugain agc=1 wagc=0.05 norm=2
# gagc on noisy data, then trap and maxbal (max subtraction) branches
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sugain gagc=1 wagc=0.03 trap=1.5 maxbal=1
# panel=1 (whole data set as one nt*ntr vector): tpow over the concatenated panel + quantile clip on the panel
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sugain panel=1 tpow=1 qclip=0.9
# gpow=0.5 sqrt branch, rms balance (pbal) and mean removal (mbal) on non-trivial data
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sugain gpow=0.5 pbal=1 mbal=1
# general pow branch (gpow=1.5) with quantile balance (qbal, qclip<1) and per-trace mark gating with jon=1
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sushw key=mark a=0 b=1 | sugain mark=1 jon=1
# epow/etpow with vred>0 and negative offsets (tred sign removal) on noisy data
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sushw key=offset a=-300 b=200 | sugain vred=2000 tpow=2 epow=1.5 etpow=1 gpow=1.5 qbal=1 qclip=0.8
# panel=1 with agc (window over the concatenated panel) and qbal at qclip<1 (quantile index from nt*ntr)
suplane nt=64 ntr=4 | suaddnoise seed=3 sn=4 | sugain panel=1 agc=1 wagc=0.1 qbal=1 qclip=0.95
```

### #4 suprofcsv (score 13.8, 1029 C lines, 13 cases / 0 data)

**Divergence 1 — the Rust is a stub.**
- No existing case ever produces a profile: all 12 exit 1 in C.
- 8 params are never named in any case: ename, formxy, keyx, keyy, nmins, point_crz, shiftdist, shiftrec.

These are the first successful cases. All 5 succeed in C (rc 0, 297–1153 bytes of Q-file) and all 5 fail in Rust.

```
# No existing case ever produces a profile (all 12 are validation errors). Parameter-mode input, initial chording, averaging and smoothing on a bent receiver line.
[stderr] python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' suprofcsv gaps=1,2,3,4,5,6,7,8 gx=0,25,50,75,100,124,147,168 gy=0,0,1,3,6,12,20,30 qout=p.q nmaxa=2 nmaxs=1 chordi=12.5 formxy=%.6f && cat p.q
# Default %.20g output, nmina/nmins > 0 (annulus averaging/smoothing), chordf different from chordi, explicit dextra.
[stderr] python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' suprofcsv gaps=1,2,3,4,5,6,7,8 gx=0,25,50,75,100,124,147,168 gy=0,0,1,3,6,12,20,30 qout=p.q nmina=1 nmaxa=3 nmins=1 nmaxs=2 chordi=10 chordf=20 dextra=15 && cat p.q
# qin= path: unsorted Q-file with an extra tuple-free name, sorted by keyp=gaps; ename= drops a name; keyy=none uses only keyx; chordf=0.
[stderr] printf 'C_SU_NAMES\nC_SU_ID,gaps,gx,gy,numb1\nQ,3,50,1,7\nQ,1,0,0,5\nQ,4,75,3,8\nQ,2,25,0,6\nQ,5,100,6,9\n' > in.q && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' suprofcsv qin=in.q qout=p.q keyy=none ename=numb1 nmaxa=1 nmaxs=0 chordi=12.5 chordf=0 formxy=%.4f && cat p.q
# shiftrec/shiftdist along-profile shifts on an existing Q-file, no chording/averaging (keyp=asis order).
[stderr] printf 'C_SU_NAMES\nC_SU_ID,gaps,gx,gy\nQ,1,0,0\nQ,2,25,0\nQ,3,50,1\nQ,4,75,3\nQ,5,100,6\nQ,6,124,12\n' > in.q && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' suprofcsv qin=in.q qout=p.q keyp=asis chordi=0 nmaxa=0 nmaxs=0 chordf=0 shiftrec=2,5 shiftdist=4,-6 formxy=%.5f && cat p.q
# point_cru/point_crz cdp numbering with the subincsv-compatible setup (keyx=keyp=gaps, keyy=none, chordi=1/point_cru).
[stderr] python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' suprofcsv gaps=1,2,3,4,5,6,7,8 gx=0,25,50,75,100,124,147,168 gy=0,0,1,3,6,12,20,30 qout=p.q keyp=gaps keyx=gaps keyy=none chordi=0.5 nmaxa=0 nmaxs=0 point_cru=2 point_crz=-3 formxy=%.3f && cat p.q
```

### #5 segywrite (score 13.4, 735 C lines, 20 cases / 7 data)

The existing cases only ever convert the values 0 and 1 to IBM floats. These cases add:
- IBM conversion of noisy, signed, varied-exponent values;
- tiny and subnormal values;
- `endian=1`;
- `conv=0` (IEEE output) together with `buff=0` (the stdio path) and `trmax`;
- a round trip through segyread with large values and set headers.

All 5 cases: SAME.

```
# IBM float conversion (conv=1, default little-endian byte swap) on noisy, signed, multi-exponent samples; existing cases only write 0/1 spike/plane values
[exact] suplane nt=17 ntr=3 | suaddnoise seed=5 sn=2 | sugain scale=1000 | segyhdrs | segywrite tape=out && cat out
# same with endian=1 (no byte swap of IBM words) and tiny magnitudes (IBM exponent near the low end)
[exact] suplane nt=17 ntr=3 | suaddnoise seed=5 sn=2 | sugain scale=1e-30 | segyhdrs | segywrite tape=out endian=1 && cat out
# conv=0 (IEEE through xdr_float) with default endian on noisy data; buff=0 (stdio fwrite path) with trmax only
[exact] suplane nt=17 ntr=5 | suaddnoise seed=6 sn=2 | segyhdrs format=5 | segywrite tape=out conv=0 buff=0 trmax=4 format=5 && cat out
# round trip through original segyread: IBM encode/decode of large values with headers incl. offsets/coords
suplane nt=33 ntr=4 | suaddnoise seed=7 sn=1 | sugain scale=12345.678 | sushw key=offset,sx,gx a=-50,100,200 b=25,10,-10 | segyhdrs | segywrite tape=out trmin=2 && segyread tape=out verbose=0 | sushw key=tracr a=0
# IEEE subnormal samples (~1e-40) through float_to_ibm: exponent underflow path of the bit manipulation
[exact] suplane nt=17 ntr=2 | suaddnoise seed=5 | sugain scale=1e-20 | sugain scale=1e-20 | segyhdrs | segywrite tape=out && cat out
```

### #6 sufctanismod (score 13.3, 2589 C lines, 69 cases / 55 data)

**Gap:** the long-run absorbing-surface and FCT paths.

**Undefined behaviour in the C (Divergence 3):**
- Any `dofct=1` run with non-zero eta reads uninitialised heap cells in the C.
- Note the defaults: the C defaults are dofct=0 and isurf=1, while the sdoc says dofct=1 and isurf=2.
- The divergence case is listed under Divergences and is **not** in the JSON.
- Cases 2 and 4 below use `dofct=1`. They pass today, but the C is heap-sensitive there.
- `movebc=1` with nx=41, nz=45 and sz=3 makes the C die with SIGBUS about 2 runs in 6.

4 cases listed: all SAME.

```
# Same long run, gradients and absorbing surface WITHOUT FCT (heap-insensitive): absorb1 isurf=1 + anisotropic gradients over 60 steps
[f32] sufctanismod nx=31 nz=33 nt=60 sx=15 sz=10 source=1 indexux=1 indexuy=1 indexuz=1 wavelet=1 dofct=0 isurf=1 aa00=3 cc00=2.2 ff00=1.1 ll00=0.8 nn00=1.2 daadz=0.5 dccdx=0.3 dffdz=-0.2 dlldx=0.1 dnndz=0.2 drhodz=0.4 > final; cat snapshotx.data snapshoty.data snapshotz.data final
# FCT on heterogeneous file-backed layered model (all six material files), source on reflectors (source=2) with Ricker, SU-headed reflection + VSP seismograms (passes today but C is heap-sensitive with dofct=1)
[f32] python3 -c 'import struct;nx,nz=21,25;f=lambda n,g: open(n,"wb").write(b"".join(struct.pack("=f",g(ix,iz)) for ix in range(nx) for iz in range(nz)));f("d",lambda x,z:2.0+(0.4 if z>12 else 0));f("a",lambda x,z:2.0+0.05*z+0.02*x);f("c",lambda x,z:2.4+(0.8 if z>12 else 0));f("ff",lambda x,z:1.0);f("l",lambda x,z:0.7+0.01*z);f("n",lambda x,z:0.9)' && sufctanismod nx=21 nz=25 nt=40 source=2 indexux=1 indexuz=1 wavelet=2 fpeak=15 dofct=1 eta0=0.02 eta=0.03 dfile=d afile=a cfile=c ffile=ff lfile=l nfile=n receiverdepth=4 vspnx=10 reflxfile=rx reflzfile=rz vspxfile=vx vspzfile=vz > final; cat rx rz vx vz | sustrip outpar=/dev/null; cat final
# FCT restricted window params + spatially varying eta/eta0, zero-value surface (isurf=3), shallow source; snapshots every 10 steps (movebc=1 omitted: the C sporadically dies with SIGBUS on this grid)
[f32] sufctanismod nx=41 nz=45 nt=50 sx=20 sz=3 source=1 indexux=1 indexuz=1 wavelet=2 fpeak=25 dofct=1 isurf=3 fctxbeg=5 fctxend=35 fctzbeg=2 fctzend=30 eta0=0.03 eta=0.04 deta0dx=0.0005 detadz=-0.0004 mt=10 > final; cat snapshotx.data snapshotz.data final
# No FCT, long run with strong density/stiffness gradients, source=3 (dipping reflector) + impulse, zero-gradient surface: accumulated FD error
[f32] sufctanismod nx=41 nz=31 nt=80 source=3 indexux=1 indexuy=1 indexuz=1 wavelet=3 isurf=2 dofct=0 rho00=1.5 drhodx=0.5 drhodz=1.0 aa00=4 daadx=1 cc00=3 dccdz=2 ff00=1.5 ll00=1 nn00=1.3 dnndx=0.4 > final; cat snapshoty.data final
```

### #7 suremel2dan (score 12.5, 2158 C lines, 89 cases / 13 data)

**Gaps:**
- The C traps (rc 133, fortified `sprintf`) whenever receiver or snapshot files are requested. So no golden ever compares a seismogram; the wavefield is visible only through the `verbose=2` per-term log.
- Every existing velocity grid is constant (exp(0)·K).

These cases use Python-generated graded and layered grids, with:
- dx ≠ dz;
- two source types;
- an anisotropic mode;
- a larger grid;
- a user wavelet with resampling.

All 5 cases: SAME.

```
# NOTE: the C build traps (rc 133, fortified sprintf) whenever receiver/snapshot files are requested, so the wavefield is only observable through verbose=2 per-term amplitudes.
# Heterogeneous isotropic model (every existing case uses constant exp(0)*scale grids): z-gradient + layer + lateral gradient; Ricker (default sflag=2); dx!=dz; absorbing iabso=1; longer nt
python3 -c 'import struct;nx,nz=33,35;f=lambda n,g: open(n,"wb").write(b"".join(struct.pack("=f",g(ix,iz)) for ix in range(nx) for iz in range(nz)));f("vp",lambda x,z:2000+15*z+(400 if z>17 else 0)+3*x);f("vs",lambda x,z:(2000+15*z+(400 if z>17 else 0))/1.8);f("dens",lambda x,z:1.0+0.02*z+(0.3 if z>17 else 0))' && python3 -c 'import subprocess,sys; p=subprocess.run(sys.argv[1:]); sys.exit(128-p.returncode if p.returncode<0 else p.returncode)' suremel2dan amode=0 nx=33 nz=35 nt=20 dx=8 dz=12 dt=0.002 isx=16 isz=14 styp=0 samp=1 fmax=40 nbwx=6 nbwz=6 verbose=2 2>log && sed "/CPU time/d" log
# Two simultaneous sources of different type, Gaussian source width w=2, non-default abso, heterogeneous model, Ricker
python3 -c 'import struct;nx,nz=33,35;f=lambda n,g: open(n,"wb").write(b"".join(struct.pack("=f",g(ix,iz)) for ix in range(nx) for iz in range(nz)));f("vp",lambda x,z:2000+15*z+(400 if z>17 else 0)+3*x);f("vs",lambda x,z:(2000+15*z+(400 if z>17 else 0))/1.8);f("dens",lambda x,z:1.0+0.02*z+(0.3 if z>17 else 0))' && python3 -c 'import subprocess,sys; p=subprocess.run(sys.argv[1:]); sys.exit(128-p.returncode if p.returncode<0 else p.returncode)' suremel2dan amode=0 nx=33 nz=35 nt=15 dx=10 dz=10 dt=0.003 isx=12,20 isz=12,16 styp=1,2 samp=1,-0.5 fmax=30 w=2 iabso=1 nbwx=5 nbwz=5 abso=0.05 verbose=2 2>log && sed "/CPU time/d" log
# Anisotropic (amode=1) with spatially varying c11/c33/c15/c55 and density, Ricker, pressure + vertical-force sources, vmaxu override
python3 -c 'import struct;nx,nz=33,35;f=lambda n,g: open(n,"wb").write(b"".join(struct.pack("=f",g(ix,iz)) for ix in range(nx) for iz in range(nz)));f("dens",lambda x,z:1.0+0.01*z);f("c11",lambda x,z:4.0e6+2e4*z);f("c13",lambda x,z:1.5e6);f("c15",lambda x,z:1e4*(x-16)/16);f("c33",lambda x,z:3.6e6+3e4*z);f("c35",lambda x,z:-5e3);f("c55",lambda x,z:1.0e6+5e3*x)' && python3 -c 'import subprocess,sys; p=subprocess.run(sys.argv[1:]); sys.exit(128-p.returncode if p.returncode<0 else p.returncode)' suremel2dan amode=1 vmax=2400 vmin=900 vmaxu=2600 nx=33 nz=35 nt=12 dx=10 dz=10 dt=0.003 isx=12,20 isz=12,12 styp=0,3 samp=1,1 fmax=30 iabso=1 nbwx=4 nbwz=4 verbose=2 2>log && sed "/CPU time/d" log
# Larger valid FFT sizes (45x63) with the default 20-point absorbing border, horizontal single force, Ricker
python3 -c 'import struct;nx,nz=45,63;f=lambda n,g: open(n,"wb").write(b"".join(struct.pack("=f",g(ix,iz)) for ix in range(nx) for iz in range(nz)));f("vp",lambda x,z:1800+10*z);f("vs",lambda x,z:900+5*z);f("dens",lambda x,z:2.0)' && python3 -c 'import subprocess,sys; p=subprocess.run(sys.argv[1:]); sys.exit(128-p.returncode if p.returncode<0 else p.returncode)' suremel2dan amode=0 nx=45 nz=63 nt=8 dx=10 dz=10 dt=0.004 isx=22 isz=25 styp=2 samp=1 fmax=25 verbose=2 2>log && sed "/CPU time/d" log
# User wavelet (sflag=0) with a wavelet dt different from the modeling dt (resampling path) on the heterogeneous model
python3 -c 'import struct;nx,nz=33,35;f=lambda n,g: open(n,"wb").write(b"".join(struct.pack("=f",g(ix,iz)) for ix in range(nx) for iz in range(nz)));f("vp",lambda x,z:2000+15*z+(400 if z>17 else 0)+3*x);f("vs",lambda x,z:(2000+15*z+(400 if z>17 else 0))/1.8);f("dens",lambda x,z:1.0+0.02*z+(0.3 if z>17 else 0))' && suspike nt=20 ntr=1 dt=0.001 ix1=1 it1=4 nspk=1 | sufilter f=10,20,60,80 > wavelet.su && python3 -c 'import subprocess,sys; p=subprocess.run(sys.argv[1:]); sys.exit(128-p.returncode if p.returncode<0 else p.returncode)' suremel2dan amode=0 nx=33 nz=35 nt=10 dx=10 dz=10 dt=0.003 isx=16 isz=16 styp=0 samp=1 sflag=0 fmax=40 iabso=0 verbose=2 2>log && sed "/CPU time/d" log
```

### #8 unif2aniso (score 11.0, 764 C lines, 14 cases / 7 data)

Branches newly exercised:
- Thomsen parameters with rotation, which existing cases test only with Sayers parameters;
- per-layer gradients with non-zero x0/z0;
- a grid that extends beyond the model.

**C quirk pinned:** line 314 reads `getparint("verbose",&paramtype)`, so `verbose=0` silently selects Sayers.

All 4 cases: SAME.

```
# Thomsen parameters (default paramtype=1) with NON-ZERO eps/delta/gamma and per-layer rotation phi00 (rotation is only tested with Sayers today), custom 2-interface dipping/curved model, default linear interpolation
[f32] printf '0 0\n1000 0\n1 -99999\n0 300\n1000 500\n1 -99999\n0 700\n300 650\n600 760\n1000 720\n1 -99999\n' > m.txt && unif2aniso ninf=2 nx=6 nz=10 dx=200 dz=100 vp00=1800,2400,3000 vs00=1000,1300,1700 rho00=2000,2200,2400 q00=50,80,120 eps00=.1,.2,.05 delta00=.05,-.05,.1 gamma00=.03,.1,.2 phi00=0,.3,-.5 < m.txt && cat c11_file c13_file c15_file c33_file c35_file c44_file c55_file c66_file rho_file q_file
# C quirk: getparint("verbose",&paramtype) - passing verbose=0 silently switches to Sayers even with paramtype=1; spline interfaces
[f32] printf '0 0\n1000 0\n1 -99999\n0 300\n1000 500\n1 -99999\n0 700\n300 650\n600 760\n1000 720\n1 -99999\n' > m.txt && unif2aniso ninf=2 nx=6 nz=10 dx=200 dz=100 method=spline paramtype=1 verbose=0 vp00=1800,2400,3000 vs00=1000,1300,1700 rho00=2000,2200,2400 eps00=.1,.2,.05 delta00=.05,-.05,.1 gamma00=.03,.1,.2 < m.txt && cat c11_file c13_file c33_file c55_file c44_file c66_file
# Per-layer gradients referenced to non-zero x0/z0 inside each layer (differing per layer, incl. negative), Akima interfaces, Thomsen
[f32] printf '0 0\n1000 0\n1 -99999\n0 300\n1000 500\n1 -99999\n0 700\n300 650\n600 760\n1000 720\n1 -99999\n' > m.txt && unif2aniso ninf=2 nx=6 nz=10 dx=200 dz=100 method=akima x0=0,500,900 z0=0,400,800 vp00=1800,2400,3000 vs00=1000,1300,1700 rho00=2000,2200,2400 q00=50,80,120 dvpdx=.2,-.3,0 dvpdz=.5,1,.8 dvsdx=0,.1,-.1 dvsdz=.3,.4,.5 drdx=.1,0,-.1 drdz=.2,.3,0 dqdx=0,.01,.02 dqdz=.05,0,.01 eps00=.1,.1,.1 dedz=.0001,-.0001,0 dddx=.00005,0,0 delta00=.02,.02,.02 dgdz=.0002,0,.0001 gamma00=.05,.05,.05 < m.txt && cat c11_file c13_file c33_file c44_file c55_file c66_file rho_file q_file
# Grid finer/larger than the model (z beyond the deepest interface, fx/fz offsets) with the shipped tfile model and monotonic interpolation
[f32] unif2aniso tfile=in.txt && unif2aniso nx=9 nz=13 dx=125 dz=100 fx=-50 fz=-100 method=mono eps00=.1,.1,.1,.1,.1,.1 phi00=.1,.1,.1,.1,.1,.1 < in.txt && cat c11_file c15_file c35_file rho_file
```

### #9 suharlan (score 10.8, 1857 C lines, 71 cases / 43 data)

Every existing data case uses linear or flat `suplane`/`suspike` events. These use hyperbolic `susynlv` gathers with seeded noise, under `gopt` 1, 2 and 4, with:
- `ninterp`;
- `pmula`/`pmulb`;
- `norm=0`;
- `prewhite`;
- asymmetric p ranges.

All 3 cases: SAME.

```
# Hyperbolic (non-linear-moveout) CMP-like events from susynlv + seeded noise: every existing data case uses suplane/suspike (linear/flat events); parabolic gopt=1 with denser np, non-default r1/r2 smoothing, two iterations
susynlv nt=128 dt=0.004 nxo=12 fxo=0 dxo=40 nxs=1 fxs=0 v00=2000 dvdz=0.5 ref="1:0,300;4000,300" ref="1:0,600;4000,640" | suaddnoise seed=11 sn=4 | suharlan seed=5 plot=0 verbose=0 dx=40 offref=440 np=21 nintlh=25 sditer=3 r1=4 r2=1 niter=2 pmin1=-100 pmax1=300 && cat out_signal out_noise
# Trace interpolation (ninterp=2) with Foster/Mosher (gopt=2) and an explicit mute window pmula/pmulb on the hyperbolic gather; norm=0
susynlv nt=128 dt=0.004 nxo=12 fxo=0 dxo=40 nxs=1 fxs=0 v00=2000 dvdz=0.5 ref="1:0,300;4000,300" ref="1:0,600;4000,640" | suaddnoise seed=12 sn=3 | suharlan seed=6 plot=0 verbose=0 dx=40 offref=440 depthref=400 gopt=2 ninterp=2 np=15 nintlh=21 sditer=2 pmula=150 pmulb=200 norm=0 && cat out_signal out_noise
# Absolute-linear transform (gopt=4) with asymmetric forward/inverse p ranges and prewhite=0.1 on noisy data; rel thresholds near 1
susynlv nt=128 dt=0.004 nxo=12 fxo=0 dxo=40 nxs=1 fxs=0 v00=2000 dvdz=0.5 ref="1:0,300;4000,300" | suaddnoise seed=13 sn=2 | suharlan seed=7 plot=0 verbose=0 dx=40 gopt=4 np=13 nintlh=21 sditer=2 pmin1=-200 pmax1=400 pmin2=-50 pmax2=300 prewhite=0.1 rel1=0.9 rel2=0.95 && cat out_signal out_noise
```

### #10 sugeomcsv (score 9.8, 3390 C lines, 119 cases / 62 data)

Branches newly exercised:
- the `match1`/`matche1` aliases;
- fixed-column records with real data;
- `create=all`;
- `action=add` with scalco.

All 4 cases: SAME. The "string modes" column overstates this program's gap: its strcmp literals are mostly C_SU_* record tags.

```
# names= aliases match1/matche1 resolved to the match= keys (two-key sorted lookup)
[stderr] printf '1,2,10\n1,3,20\n2,2,30\n' > geom.csv && sunull nt=3 ntr=3 | sushw key=fldr,tracf a=1,2 b=0,1 c=1,0 j=2,2 | sugeomcsv rfile=geom.csv ufile=command setid=any match=fldr,tracf names=match1,matche1,selev
# Fixed-column records with real data (rtype=fixed, start_name_end column ranges, C_su_id column)
[stderr] printf '%s\n' 'S  1  10' 'S  2 -25' 'R  3  99' 'S  3 300' > geom.txt && sunull nt=2 ntr=3 | sugeomcsv rfile=geom.txt ufile=command rtype=fixed setid=S match=tracl names=C_su_id,2_tracl_4,5_selev_8
# X-record csv with two records, channel increments (ci=2) and reversed receiver range: create=all trace generation
[stderr] printf '0,1,1,5,2,101,103\n0,2,9,5,2,203,201\n' > geom.csv && python3 -c 'import pty,subprocess,sys; m,s=pty.openpty(); sys.exit(subprocess.call(" sugeomcsv rfile=geom.csv ufile=command setid=any match=fldr,tracf names=null,fldr,tracf_cf,tracf_ct,tracf_ci,gaps_rf,gaps_rt create=all",shell=True,stdin=s))'
# Coordinate update with action=add on existing headers, scalco=100 and fractional values (rounding into int keys)
[stderr] printf '1,10.126,-3.004\n2,-0.555,7.5\n' > geom.csv && sunull nt=2 ntr=2 | sushw key=sx,gx,scalco a=1000,2000,-100 | sugeomcsv rfile=geom.csv ufile=command setid=any match=tracl names=tracl,sx,gx action=add scalco=100
```

### #11 sumatrix (score 8.7, 484 C lines, 27 cases / 7 data)

Branches newly exercised:
- fractional and negative dividers;
- float header keys;
- input longer than the zone;
- two-key sources with `rfill=2`.

All 4 cases: SAME.

```
# Fractional and negative dividers on coordinate keys: floor() of negative quotients, tohead of int*double products back into int headers
[stderr] suspike ntr=8 nt=6 dt=0.004 | sushw key=sx,gx a=-130,45 b=37,-23 j=4,8 c=0,0 | sumatrix skeyloc=sx sdivider=2.5 rkeyloc=gx rdivider=-3 lenzone=8 numzone=2
# Float header keys (f2,d2) as locations: fromhead/tohead float path with fractional divider
[stderr] suspike ntr=6 nt=5 dt=0.004 | sushw key=f2,d2 a=1.25,0.5 b=0.75,-0.25 | sumatrix skeyloc=f2 sdivider=0.5 rkeyloc=d2 rdivider=0.25 lenzone=12 numzone=3 rfill=0
# Input traces longer than the zone (lensam truncated) and lenzone rounded by NINT (10ms/4ms -> 3 samples)
[stderr] suplane ntr=6 nt=16 dt=0.004 | sushw key=fldr,gaps a=1,1 b=0,1 c=1,0 j=3,3 | sumatrix lenzone=10 numzone=2
# Duplicate source/receiver cells (REPLACED warning), multi-key locations, numzone larger than source count, rfill=2
[stderr] suspike ntr=8 nt=4 dt=0.002 | sushw key=fldr,ep,gaps a=1,3,1 b=0,0,1 c=1,0,0 j=4,8,2 | sumatrix skeyloc=fldr,ep rkeyloc=gaps lenzone=6 numzone=5 rfill=2
```

### #12 sumute (score 8.2, 345 C lines, 29 cases / 15 data)

Branches newly exercised:
- non-zero `delrt`;
- offsets beyond the table;
- bottom mute on a depth section;
- `hmute`;
- hyperbolic and air-wave mutes with tapers.

All 5 cases: SAME.

**Do not add as goldens:**
- `mode=4` reads an uninitialised `twfval` in the C; the port uses 0.
- `mode=2` with `absolute=0` and negative offsets can produce a negative `memset` length in the C.

```
# mode=0 with nonzero delrt (tmin), negative offsets, extrapolation outside xmute range, 3-point interpolation, taper on noisy data
suplane nt=64 ntr=5 | suaddnoise seed=2 sn=3 | sushw key=delrt,offset a=40,-200 b=0,100 | sumute xmute=-150,0,150 tmute=0.1,0.05,0.1 ntaper=5
# depth section (trid=130, d1) with mode=1 bottom mute and taper; mute header written as NINT(t) in depth units
suplane nt=64 ntr=4 | suaddnoise seed=2 sn=3 | sushw key=trid,d1 a=130,5 | sumute key=tracl xmute=1,4 tmute=120,250 mode=1 ntaper=3
# hmute (per-trace mute time from header) with taper and delrt
suplane nt=64 ntr=4 | suaddnoise seed=2 sn=3 | sushw key=muts,delrt a=40,20 b=30,0 | sumute hmute=muts ntaper=4
# mode=3 hyperbolic mute with delrt, negative offsets (absolute=1 default) and taper on both sides
suplane nt=128 ntr=5 | suaddnoise seed=2 sn=3 | sushw key=offset,delrt a=-200,8 b=100,0 | sumute xmute=-200,200 tmute=0.02,0.03 mode=3 linvel=2000 tm0=0.1 ntaper=4
# key of a different header type (sx, int) plus mode=2 air-wave mute with tm0 and taper on noisy data
suplane nt=128 ntr=5 | suaddnoise seed=2 sn=3 | sushw key=sx a=0 b=150 | sumute key=sx xmute=0,600 tmute=0.02,0.06 mode=2 linvel=3000 tm0=0.05 ntaper=3
```

### #13 sukdmig3d (score 8.1, 1334 C lines, 96 cases / 56 data)

**Gap:** every existing traveltime table is a constant (or random) 3x3 table.

**Gotcha for case writers:** the C build has its DEBUG block active, which forces sy=gy=3. Tables therefore need fys ≤ 3.

All 4 cases: SAME.

```
# Physically varying traveltime tables (hypot distance / 2500 m/s, ms) on a 5x4 grid with 2x2x3 sources, instead of the constant/random 3x3 tables every existing case uses; 4 traces with distinct offsets all migrate.
[stderr] python3 -c 'import struct,math; open("t.bin","wb").write(b"".join(struct.pack("=f",1000*math.sqrt((x-xs)**2+(y-ys)**2+zs**2)/2500) for ys in (0,10) for xs in (10,30) for zs in (10,20,30) for y in (0,10,20,30) for x in (0,10,20,30,40)))'; suaddhead ns=5 <t.bin | sushw key=f1,d1,f2,d2,ntr,selev,gelev,nhs,nvs,duse,sdel,gdel,sdepth,swdep,gwdep,ep a=0,10,0,10,4,5,4,2,2,3,1,0,10000,2,1,10000 >tt; suspike nt=65 ntr=4 dt=0.001 nspk=1 ix1=1 it1=30 | suaddnoise seed=3 sn=20 | sushw key=d2,sx,gx,sy,gy a=10,10,20,10,20 b=0,5,-2,2,1 >data; sukdmig3d datain=data ttfile=tt >out && cat out
# Real x-offset binning: offsets 10,3,-4,-11 into nxoffset=2 panels (fxoffset=-10 dxoffset=10) and xoffsetmax=10 rejects the -11 trace.
[stderr] python3 -c 'import struct,math; open("t.bin","wb").write(b"".join(struct.pack("=f",1000*math.sqrt((x-xs)**2+(y-ys)**2+zs**2)/2500) for ys in (0,10) for xs in (10,30) for zs in (10,20,30) for y in (0,10,20,30) for x in (0,10,20,30,40)))'; suaddhead ns=5 <t.bin | sushw key=f1,d1,f2,d2,ntr,selev,gelev,nhs,nvs,duse,sdel,gdel,sdepth,swdep,gwdep,ep a=0,10,0,10,4,5,4,2,2,3,1,0,10000,2,1,10000 >tt; suspike nt=65 ntr=4 dt=0.001 nspk=1 ix1=1 it1=30 | suaddnoise seed=3 sn=20 | sushw key=d2,sx,gx,sy,gy a=10,10,20,10,20 b=0,5,-2,2,1 >data; sukdmig3d datain=data ttfile=tt nxoffset=2 fxoffset=-10 dxoffset=10 xoffsetmax=10 >out && cat out
# Amplitude-preserving branch ls=1 with a narrow aperture/angle and a finer custom output grid.
[stderr] python3 -c 'import struct,math; open("t.bin","wb").write(b"".join(struct.pack("=f",1000*math.sqrt((x-xs)**2+(y-ys)**2+zs**2)/2500) for ys in (0,10) for xs in (10,30) for zs in (10,20,30) for y in (0,10,20,30) for x in (0,10,20,30,40)))'; suaddhead ns=5 <t.bin | sushw key=f1,d1,f2,d2,ntr,selev,gelev,nhs,nvs,duse,sdel,gdel,sdepth,swdep,gwdep,ep a=0,10,0,10,4,5,4,2,2,3,1,0,10000,2,1,10000 >tt; suspike nt=65 ntr=4 dt=0.001 nspk=1 ix1=1 it1=30 | suaddnoise seed=3 sn=20 | sushw key=d2,sx,gx,sy,gy a=10,10,20,10,20 b=0,5,-2,2,1 >data; sukdmig3d datain=data ttfile=tt ls=1 xaper=8 angmax=40 nzo=21 dzo=1 nxo=17 dxo=2.5 >out && cat out
# Midpoints outside [fxs,exs] (xm<10) are skipped while in-range traces migrate (sx=0..15, gx=0..6).
[stderr] python3 -c 'import struct,math; open("t.bin","wb").write(b"".join(struct.pack("=f",1000*math.sqrt((x-xs)**2+(y-ys)**2+zs**2)/2500) for ys in (0,10) for xs in (10,30) for zs in (10,20,30) for y in (0,10,20,30) for x in (0,10,20,30,40)))'; suaddhead ns=5 <t.bin | sushw key=f1,d1,f2,d2,ntr,selev,gelev,nhs,nvs,duse,sdel,gdel,sdepth,swdep,gwdep,ep a=0,10,0,10,4,5,4,2,2,3,1,0,10000,2,1,10000 >tt; suspike nt=65 ntr=4 dt=0.001 nspk=1 ix1=1 it1=30 | suaddnoise seed=3 sn=20 | sushw key=d2,sx,gx,sy,gy a=10,0,0,10,20 b=0,5,2,2,1 >data; sukdmig3d datain=data ttfile=tt >out && cat out
```

### #14 sutxtaper (score 7.7, 511 C lines, 23 cases / 17 data)

Branches newly exercised: the `scalco` and `scalel` coordinate scaling, taper type 5, and gaussian tapers.

All 3 cases: SAME.

```
# coordinate key (isCoord): x = sx * 10^scalco (C uses pow(10,scalco), not the SEG-Y divide convention), gaussian x-taper
suplane nt=31 ntr=7 | suaddnoise seed=7 sn=2 | sushw key=sx,scalco a=1000,-1 b=500,0 | sutxtaper key=sx min=150 max=350 dx=100 taper=4
# elevation key (isDepth, scalel) with taper=5 gaussian, low>0 and an end-only time taper (tbeg=0, tend>0)
suplane nt=31 ntr=7 | suaddnoise seed=7 sn=2 | sushw key=gelev,scalel a=100,1 b=50,0 | sutxtaper key=gelev min=1500 max=3500 dx=1500 taper=5 low=0.2 tend=40
# trace-count taper with tr1 != tr2 (separate end-taper branch) and gaussian shape plus begin time taper on noisy data
suplane nt=40 ntr=9 | suaddnoise seed=7 sn=2 | sutxtaper tr1=3 tr2=2 taper=4 tbeg=20 low=0.25
```

### #15 sukdmdcr (score 7.5, 931 C lines, 62 cases / 22 data)

Branches newly exercised:
- off-grid receivers;
- a finer table with fzt/fxso ≠ 0 and `offmax`;
- upward continuation.

All 3 cases: SAME.

**C bugs to stay clear of:**
- nxso>1 makes the C err "dxso must = dxgo".
- A negative `fxgo` breaks the `nxi` formula.

```
# Off-grid receivers (gx=5,15,25 between table nodes -> res=0.5 in sum2 interpolation of ttab).
[stderr] python3 -c 'import struct,math; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(10*x-10*s,10*z)/1500) for s in range(5) for x in range(5) for z in range(7)))'; suspike nt=65 ntr=6 nspk=4 ix1=1 it1=20 ix2=2 it2=30 ix3=4 it3=40 ix4=6 it4=50 | sushw key=sx,gx a=0,5 b=0,10 c=10,10 j=3,3 > in; sukdmdcr ttfile=tt nxt=5 fxt=0 dxt=10 nzt=7 fzt=0 dzt=10 ns=5 fs=0 ds=10 nxso=2 dxso=10 nxgo=3 dxgo=10 zrec=0 zdat=15 datain=in dataout=out jpfile=job; cat out
# Different table geometry: finer spacing (ds=dxt=5, ns=nxt=9), nonzero fzt and fxso, v0=1700, aperx=8; noise input; offmax=12 drops the far (15 m) trace of each gather. (Negative fxgo is avoided: nxi=(int)((fxgo+...)/dxi)+1 shrinks the topography grid and errs.)
[stderr] python3 -c 'import struct,math; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(5*x-5*s,4+6*z)/1700) for s in range(9) for x in range(9) for z in range(6)))'; suplane nt=65 ntr=8 | suaddnoise seed=17 sn=3 | sushw key=sx,gx a=5,5 b=0,5 c=10,10 j=4,4 > in; sukdmdcr ttfile=tt nxt=9 fxt=0 dxt=5 nzt=6 fzt=4 dzt=6 ns=9 fs=0 ds=5 nxso=4 fxso=5 dxso=5 nxgo=4 fxgo=0 dxgo=5 zrec=4 zdat=20 v0=1700 offmax=12 aperx=8 datain=in dataout=out jpfile=job; cat out
# Upward continuation (zdat < zrec) with a narrow aperture and non-default freq/scale.
[stderr] python3 -c 'import struct,math; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(10*x-10*s,10*z)/1500) for s in range(5) for x in range(5) for z in range(7)))'; suspike nt=65 ntr=6 nspk=4 ix1=1 it1=20 ix2=2 it2=30 ix3=4 it3=40 ix4=6 it4=50 | sushw key=sx,gx a=0,0 b=0,10 c=10,10 j=3,3 > in; sukdmdcr ttfile=tt nxt=5 fxt=0 dxt=10 nzt=7 fzt=0 dzt=10 ns=5 fs=0 ds=10 nxso=2 dxso=10 nxgo=3 dxgo=10 zrec=40 zdat=10 aperx=8 freq=30 scale=2 datain=in dataout=out jpfile=job; cat out
```

### #16 sugassmann (score 7.5, 760 C lines, 17 cases / 8 data)

Only mode 0 ever produces output: mode 1 writes nothing, and the mode=2 block is dead code. These cases vary the physics inputs instead.

All 4 cases: SAME.

```
# Heavy oil / condensate gas / fresh brine at shallow, cold conditions: fluid-property polynomials far from defaults
sugassmann api=10 g=1.8 s=.02 temp=15 pres=5 so=.3 sg=.4
# Stiff low-porosity frame with non-unit exponent c: mdry=mmin/(a+b*phi^c) path
sugassmann phi=.05 c=2.5 a=.8 b=30 kmin=70000 mumin=30000 rmin=2.71 v1=4500 r1=2.5
# Full sensitivity sweep p=1 (y from 0 to 2*sg, crossing sg=0 at first sample) with tiny gas saturation and no oil
sugassmann sg=.01 so=0 p=1 g=.6 api=30
# High-temperature/high-pressure brine with salt-saturated water and soft caprock
sugassmann temp=250 pres=300 s=.8 so=.1 sg=.05 v1=2000 r1=2.1 p=.5
```

### #17 sunomstreamers (score 6.8, 717 C lines, 32 cases / 7 data)

Branches newly exercised:
- `offpi=1`;
- reversed channel ranges;
- streamer ids with `channelkey`;
- `cdpspacing=0`;
- `create=` with `gunkey`/`firstshot`.

All 5 cases: SAME.

```
# offpi=1 (|inline difference|), gunid order 1,0 with fldr modulus, one channel range expanded contiguously over two streamers (no streamerkey)
[stderr] sunull nt=3 ntr=12 | sushw key=fldr,tracf a=1,1 b=0,1 c=1,0 j=6,6 | sunomstreamers shotspacing=25 gunid=1,0 guncross=-25,25 guninline=-15 streamercross=-50,50 streamerinline=-120,-130 channelstart=1 channelend=3 channelspacing=12.5 offpi=1
# Reversed single channel range expanded to two streamers (stcbeg>stcend contiguous branch, stioff shift), default offpi/linelocs/scalco
[stderr] sunull nt=2 ntr=8 | sushw key=fldr,tracf a=3,1 b=0,1 | sunomstreamers shotspacing=37.5 gunid=0 guncross=10 guninline=-8 streamercross=-30,30 streamerinline=-90 channelstart=4 channelend=1 channelspacing=15
# streamerkey/streameridnt with a single channel range shared by all streamers, non-default channelkey (tracr), scalco=1
[stderr] sunull nt=2 ntr=8 | sushw key=fldr,ep,tracr a=5,7,3 b=0,1,0 c=1,0,0 j=4,4,4 | sunomstreamers shotspacing=50 gunid=0 guncross=0 guninline=-20 streamercross=-75,-25,25,75 streamerinline=-100 streamerkey=ep streameridnt=7,8,9,10 channelkey=tracr channelstart=3 channelend=3 channelspacing=12.5 scalco=1
# cdpspacing=0 while linelocs=1 stays on: grnlof/gaps = lrint(x/0) (inf/nan to int conversion)
[stderr] sunull nt=2 ntr=3 | sushw key=fldr,tracf a=0,1 b=0,1 | sunomstreamers shotspacing=25 gunid=0 guncross=0 guninline=0 streamercross=0 streamerinline=-100 channelstart=1 channelend=3 channelspacing=50 cdpspacing=0
# create=3 shots from firstshot=10 with separate gunkey (cycled gun idents) and spike times not on the sample grid
[stderr] python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s); sys.exit(p.returncode)' sunomstreamers shotspacing=25 gunkey=ep gunid=0,1 guncross=-20,20 guninline=-10 streamercross=-40,40 streamerinline=-100 channelstart=1 channelend=2 channelspacing=12.5 create=3 firstshot=10 spikes=4,0.5,10,1,22,-1
```

### #18 sucmp (score 6.7, 184 C lines, 19 cases / 0 data)

**Metric false positive.** sucmp reports its result through the exit code plus a short stdout line, so it scores "0 data". The proposals add `[stderr]` comparisons on noisy data, the sign-flip ratio branch, and a header difference on a later trace.

All 4 cases: SAME.

```
# Dense noisy data, B = A with independent small noise: first out-of-band sample found mid-trace; the %15g print of A/B values (existing cases only compare constants 1/-1).
[stderr] suplane nt=64 ntr=4 | suaddnoise seed=3 sn=2 > a && suaddnoise < a seed=4 sn=2000 > b && sucmp a b limit=0.001
# Same pair within a loose limit (all samples pass, exit 0) and after a sign flip (ratio negative -> fails at the first nonzero sample).
[stderr] suplane nt=64 ntr=4 | suaddnoise seed=3 sn=2 > a && suaddnoise < a seed=4 sn=2000 > b && sucmp a b limit=0.5
[stderr] suplane nt=64 ntr=4 | suaddnoise seed=3 sn=2 > a && suop < a op=neg > b && sucmp a b
# Difference only in the last trace's header (later-trace header branch, not trace 1).
[stderr] suplane nt=16 ntr=5 > a && sushw < a key=offset a=0 b=0 c=0 j=4 | sushw key=cdp a=0 b=0 j=5 > b && sushw < a key=cdp a=0 b=0 j=5 > c && sucmp c b; echo rc=$?
```

### #19 suvelan_nsel (score 6.6, 346 C lines, 14 cases / 8 data)

Every existing data case uses linear `suplane` events. These use hyperbolic `susynlv` CMP gathers, so NMO and the stretch mute actually act.

All 4 cases: SAME.

**C notes:**
- `pwr` is checked for validity but never used.
- Input with more than one CDP is undefined behaviour in the C: `ix` is never reset, so the next trace writes past the end of `aoff`/`datm`. Avoid it.

```
# Hyperbolic CMP gather (susynlv, offsets 50..400 m) instead of suplane's linear events: the NMO interpolation and stretch mute (smute=1.5) actually act; non-default fv/dv/dx/tau/dtratio.
susynlv nt=126 dt=0.004 nxo=8 fxo=0.05 dxo=0.05 nxs=1 ref="0,0.25;2,0.25" ref="0,0.5;2,0.5" v00=2 fpeak=30 | suvelan_nsel nx=8 dx=50 nv=6 fv=1600 dv=200 tau=0.3 smute=1.5 dtratio=2
# Same gather with seeded noise, negative offsets (|offset| path for ol), delrt=20 ms so ft>0 changes itmute and the interpolation origin; tau=0.7, nsmooth=9, pwr=3.
susynlv nt=126 dt=0.004 nxo=8 fxo=0.05 dxo=0.05 nxs=1 ref="0,0.25;2,0.25" ref="0,0.5;2,0.5" v00=2 fpeak=30 | suaddnoise seed=31 sn=6 | sushw key=offset,delrt a=-400,20 b=50,0 | suvelan_nsel nx=8 dx=50 nv=5 fv=1800 dv=150 tau=0.7 nsmooth=9 dtratio=3 pwr=3
# nx and dx taken from cdpt/d2 headers with default nv=50 fv=1500 dv=50 dtratio=5 on the hyperbolic gather (only flat suplane data used the defaults before).
susynlv nt=126 dt=0.004 nxo=8 fxo=0.05 dxo=0.05 nxs=1 ref="0,0.25;2,0.25" ref="0,0.5;2,0.5" v00=2 fpeak=30 | sushw key=cdpt,d2 a=8,50 | suvelan_nsel
# Fractional dx=12.5 with large minimum offset (ol=m/dx=4 => different j0t/beta), tau=0.1 includes most crosscorrelations.
susynlv nt=126 dt=0.004 nxo=8 fxo=0.05 dxo=0.05 nxs=1 ref="0,0.25;2,0.25" ref="0,0.5;2,0.5" v00=2 fpeak=30 | sushw key=offset a=50 b=12 | suvelan_nsel nx=8 dx=12.5 nv=4 fv=1900 dv=40 tau=0.1 dtratio=1 nsmooth=3 smute=3
```

### #20 raydata (score 6.5, 257 C lines, 54 cases / 5 data)

**Metric false positive.** The existing cases are thorough, but raydata exits 1 even on success and its records are tiny. This still found that both inventory cases (`elamodel … && elaray … && raydata`) never run raydata at all.

The proposals use real elaray endpoints with 21 rays.

All 3 cases: SAME.

```
# The two inventory cases chain `elamodel ... && elaray ... && raydata`, but elamodel/elaray exit 1 on success, so raydata never runs in them. Real elaray endpoints for a 2-layer model with many rays, all binary outputs.
[f32] elamodel fill=0.5,0.5,2,1,1 > model.bin; elaray nangle=21 fangle=-40 langle=40 < model.bin > rays.bin; { raydata t=1 px=1 pz=1 vgx=1 vgz=1 pxvx=1 pol=1 < rays.bin; s=$?; cat x_t.data x_px.data x_pz.data x_vgx.data x_vgz.data px_vgx.data x_pol.data; exit $s; }
# kend selection on real endpoints (only rays ending at interface 0) plus ASCII output and outpar count.
[exact] elamodel fill=0.5,0.5,2,1,1 > model.bin; elaray nangle=21 fangle=-40 langle=40 < model.bin > rays.bin; { raydata kend=0 t=1 pol=1 ascci=1 < rays.bin; s=$?; cat x_t.data x_pol.data; exit $s; }
[exact] elamodel fill=0.5,0.5,2,1,1 > model.bin; elaray nangle=21 fangle=-40 langle=40 < model.bin > rays.bin; { raydata t=1 < rays.bin; s=$?; cat outpar; exit $s; }
```

### #21 sukdmig2d (score 6.1, 919 C lines, 79 cases / 45 data)

**Gap:** all 79 existing cases use one 5x5x3 traveltime table at fs=0. These use a new 9x6x5 table at fs=-25 with noisy input, adding:
- off-node snapping and interpolation;
- offset panels with `absoff`/`limoff`;
- `offmax` together with `mtr`.

All 4 cases: SAME.

```
# New table geometry (nxt=9 fxt=-25 dxt=25, nzt=6 fzt=10 dzt=30, ns=5 fs=-25 ds=50: every existing case uses one 5x5x3 table at fs=0), dense noisy multi-sample input, and a first source outside [fs,es] that is skipped (7 of 8 traces migrate).
[stderr] python3 -c 'import math,struct;open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(-25+x*25-(-25+s*50),10+z*30)/1700) for s in range(5) for x in range(9) for z in range(6)))' && suplane nt=101 ntr=8 dt=0.002 | suaddnoise seed=21 sn=4 | sushw key=sx,gx,d2 a=-50,175,25 b=30,-25,0 > in && sukdmig2d infile=in outfile=out ttfile=tt nxt=9 fxt=-25 dxt=25 nzt=6 fzt=10 dzt=30 ns=5 fs=-25 ds=50 v0=1700 && cat out
# Sources/receivers 0.5 m and 24.6 m off table nodes: hits the res<=0.01 / res>=0.99 snapping and true interpolation in sum2 on a non-unit geometry; vertical gradient + narrow aperture.
[stderr] python3 -c 'import math,struct;open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(-25+x*25-(-25+s*50),10+z*30)/1700) for s in range(5) for x in range(9) for z in range(6)))' && suplane nt=101 ntr=4 dt=0.002 | suaddnoise seed=22 sn=4 | sushw key=sx,gx,d2 a=0.4,74.6,25 b=49.2,0.3,0 > in && sukdmig2d infile=in outfile=out ttfile=tt nxt=9 fxt=-25 dxt=25 nzt=6 fzt=10 dzt=30 ns=5 fs=-25 ds=50 v0=1700 dvz=0.4 aperx=40 angmax=45 && cat out
# Offset panels on the new geometry: noff=3 with absoff=1 limoff=1 drops traces outside the bins, ls=0, finer output grid.
[stderr] python3 -c 'import math,struct;open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(-25+x*25-(-25+s*50),10+z*30)/1700) for s in range(5) for x in range(9) for z in range(6)))' && suplane nt=101 ntr=8 dt=0.002 | suaddnoise seed=23 sn=4 | sushw key=sx,gx,d2 a=-25,175,25 b=25,-25,0 > in && sukdmig2d infile=in outfile=out ttfile=tt nxt=9 fxt=-25 dxt=25 nzt=6 fzt=10 dzt=30 ns=5 fs=-25 ds=50 v0=1700 noff=3 off0=0 doff=60 absoff=1 limoff=1 ls=0 nzo=31 dzo=5 nxo=21 dxo=10 && cat out
# offmax=120 rejects far traces and mtr/ntr progress limits on a longer gather (ntr=6 of 8).
[stderr] python3 -c 'import math,struct;open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(-25+x*25-(-25+s*50),10+z*30)/1700) for s in range(5) for x in range(9) for z in range(6)))' && suplane nt=101 ntr=8 dt=0.002 | suaddnoise seed=24 sn=4 | sushw key=sx,gx,d2 a=-25,175,25 b=25,-25,0 > in && sukdmig2d infile=in outfile=out ttfile=tt nxt=9 fxt=-25 dxt=25 nzt=6 fzt=10 dzt=30 ns=5 fs=-25 ds=50 v0=1700 offmax=120 ntr=6 mtr=2 && cat out
```

### #22 pswigp (score 6.0, 499 C lines, 15 cases / 8 data)

Branches newly exercised:
- two curves with list propagation;
- `perc` clipping on dense data;
- an axis-2 window with reversed limits and a negative d1;
- irregular `x2` with `style=normal`.

All 4 cases: SAME.

```
# Two curves with list propagation: curvecolor/curvewidth/curvedash given once (or omitted) are copied to the second curve; npair per curve.
printf '0 0\n4 1\n8 2\n' > c1; printf '1 2\n3 1.5\n5 0.5\n7 0\n' > c2; suplane ntr=3 nt=9 | sustrip outpar=/dev/null | pswigp n1=9 n2=3 curve=c1,c2 npair=3,4 curvecolor=red curvewidth=1.5 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# Clip from percentile on dense noisy data (qkfind on real |z| values) with bias and default style=seismic.
suplane ntr=6 nt=33 | suaddnoise seed=11 sn=2 | sustrip outpar=/dev/null | pswigp n1=33 n2=6 perc=90 bias=0.1 fill=1 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# Axis-2 window keeps only a subset of traces (n2in < n2, xcur scaled by window), reversed x2beg>x2end, and x1beg/x1end snapped (NINT, clamped) on a negative d1 axis.
suplane ntr=8 nt=21 | suaddnoise seed=4 sn=3 | sustrip outpar=/dev/null | pswigp n1=21 n2=8 d2=10 f2=100 x2beg=160 x2end=120 d1=-0.01 f1=0.2 x1beg=0.137 x1end=-5 xcur=2 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# Irregular x2 positions (x2= list, non-monotonic) with style=normal and verbose clip echo on stderr.
[stderr] suplane ntr=4 nt=17 | suaddnoise seed=21 sn=1 | sustrip outpar=/dev/null | pswigp n1=17 x2=3,0.5,7,2 style=normal perc=50 verbose=1 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
```

### #23 sumedian (score 6.0, 488 C lines, 36 cases / 15 data)

Branches newly exercised: varying, non-zero shifts along dipping events on noisy data, the mix path, a 3-knot shift curve with `subtract=0`, and `nmed=9`.

All 4 cases: SAME.

**C crashes, not port bugs:**
- `sign=1` with any non-zero shift writes outside the buffer, so the C segfaults (rc 139). The port skips those writes.
- At the end of the gather the C reads `inshift[ntr]`, one past the end of the array.

Keep both out of goldens.

```
# Median filter (median=1, default nmed=5) along a dipping event: varying non-zero tshift (0 -> 56 ms over 15 traces) on noisy data, default subtract=1 and sign=-1 -- every existing data case uses flat suplane events with tshift=0 or a small constant shift.
suplane npl=1 nt=64 ntr=15 dt=0.004 ct1=8 dip1=4 len1=15 | suaddnoise seed=41 sn=5 | sumedian xshift=1,15 tshift=0,0.056 median=1
# Weighted-mix (non-median) path with 5 unequal weights along the same varying shift curve, subtract=1.
suplane npl=1 nt=64 ntr=15 dt=0.004 ct1=8 dip1=4 len1=15 | suaddnoise seed=42 sn=5 | sumedian xshift=1,15 tshift=0,0.056 mix=0.5,1,2,1,0.5
# Piecewise-linear 3-knot shift curve on key=offset, median=1 nmed=7, subtract=0 (filtered estimate output), tmpdir=. disk path.
suplane npl=1 nt=64 ntr=15 dt=0.004 ct1=10 dip1=3 len1=15 | suaddnoise seed=43 sn=4 | sushw key=offset a=100 b=50 | sumedian key=offset xshift=100,450,800 tshift=0,0.016,0.044 median=1 nmed=7 subtract=0 tmpdir=.
# Larger gather (ntr=31) with two crossing events, median along one of them with nmed=9 (wider ring buffer, larger shiftmax=30).
suplane npl=2 nt=48 ntr=31 dt=0.004 ct1=6 dip1=1 len1=31 ct2=40 dip2=-1 len2=31 | suaddnoise seed=44 sn=8 | sumedian xshift=1,31 tshift=0,0.12 median=1 nmed=9
```

### #24 sudatumk2ds (score 5.9, 1166 C lines, 70 cases / 27 data)

Branches newly exercised:
- off-grid sources;
- out-of-order receiver gathers;
- the "Topography definition is out of range" error;
- `offmax` skipping;
- topography with `dvz`.

All 5 cases: SAME. **C quirk:** if the first trace of a gather is skipped, the C outputs all zeros.

```
# Several receiver gathers with off-grid sources (sx=7,13,27 -> res interpolation between tt shots and between surface nodes) and the gx!=s new-gather branch (fxin=sx).
suspike nt=65 ntr=6 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=7,10 b=6,0 c=0,20 j=3,3 >in; python3 -c 'import struct,math; open("tt","wb").write(struct.pack("=75f",*[math.hypot(x-s,z)/1500 for s in (0,20,40) for x in (0,10,20,30,40) for z in (0,5,10,15,20)]))'; sudatumk2ds infile=in outfile=out ttfile=tt nxt=5 fxt=0 dxt=10 nzt=5 fzt=0 dzt=5 ns=3 fs=0 ds=20 nxi=9 fxi=0 dxi=5 sgn=1 'surf=-10,0;50,0' 'surf=-10,5;50,5' && cat out
# Receiver gathers arrive out of order (gx=30 then gx=10): ng[io] counts per io but headers are written back in input order -- checks the C header/data pairing quirk.
suspike nt=65 ntr=4 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=0,30 b=20,0 c=0,-20 j=2,2 >in; python3 -c 'import struct,math; open("tt","wb").write(struct.pack("=75f",*[math.hypot(x-s,z)/1500 for s in (0,20,40) for x in (0,10,20,30,40) for z in (0,5,10,15,20)]))'; sudatumk2ds infile=in outfile=out ttfile=tt nxt=5 fxt=0 dxt=10 nzt=5 fzt=0 dzt=5 ns=3 fs=0 ds=20 nxi=9 fxi=0 dxi=5 sgn=1 'surf=-10,0;50,0' 'surf=-10,5;50,5' && cat out
# Error path: a gather whose first source is outside the table (fxin=sx=-10) trips "Topography definition is out of range!".
[stderr] suspike nt=65 ntr=4 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=-10,20 b=25,0 >in; python3 -c 'import struct,math; open("tt","wb").write(struct.pack("=75f",*[math.hypot(x-s,z)/1500 for s in (0,20,40) for x in (0,10,20,30,40) for z in (0,5,10,15,20)]))'; sudatumk2ds infile=in outfile=out ttfile=tt nxt=5 fxt=0 dxt=10 nzt=5 fzt=0 dzt=5 ns=3 fs=0 ds=20 nxi=9 fxi=0 dxi=5 sgn=1 'surf=-10,0;50,0' 'surf=-10,5;50,5' && cat out
# offmax filter: of sx=10,20,30,40 at gx=20 the |offset|=20 trace is skipped (3 of 4 output traces). Note: when the FIRST trace of a gather is skipped the C output is all zeros (fxin quirk), so keep the first trace in range.
suspike nt=65 ntr=4 nspk=1 ix1=1 it1=17 | sushw key=sx,gx a=10,20 b=10,0 >in; python3 -c 'import struct,math; open("tt","wb").write(struct.pack("=75f",*[math.hypot(x-s,z)/1500 for s in (0,20,40) for x in (0,10,20,30,40) for z in (0,5,10,15,20)]))'; sudatumk2ds infile=in outfile=out ttfile=tt nxt=5 fxt=0 dxt=10 nzt=5 fzt=0 dzt=5 ns=3 fs=0 ds=20 nxi=9 fxi=0 dxi=5 sgn=1 'surf=-10,0;50,0' 'surf=-10,5;50,5' offmax=15 && cat out
# Topographic surfaces with a velocity gradient (sz != 0 in timeb) and narrow aperture.
suplane nt=65 ntr=3 dt=0.004 | suaddnoise seed=5 sn=3 | sushw key=sx,gx a=5,20 b=15,0 >in; python3 -c 'import struct,math; open("tt","wb").write(struct.pack("=75f",*[math.hypot(x-s,z)/1500 for s in (0,20,40) for x in (0,10,20,30,40) for z in (0,5,10,15,20)]))'; sudatumk2ds infile=in outfile=out ttfile=tt nxt=5 fxt=0 dxt=10 nzt=5 fzt=0 dzt=5 ns=3 fs=0 ds=20 nxi=9 fxi=0 dxi=5 sgn=1 'surf=-10,1;12,3;30,0;50,2' 'surf=-10,6;20,9;50,7' v0=1700 dvz=0.6 aperx=8 angmax=40 && cat out
```

### #25 pscubecontour (score 5.9, 607 C lines, 16 cases / 6 data)

Branches newly exercised:
- a smooth 3-D field instead of sparse spikes or a ramp;
- an interior window on all three axes;
- face files with contour values outside the data range.

All 3 cases: SAME.

```
# Smooth 3-D field (every existing data case is sparse suplane spikes or a ramp): default nc/fc/dc picked from real min/max, contour tracing across all three faces, faces=1, steeper angle.
python3 -c 'import math,struct,sys; sys.stdout.buffer.write(b"".join(struct.pack("f",math.sin(.5*i1)*math.cos(.4*i2)+.3*i3) for i3 in range(6) for i2 in range(7) for i1 in range(9)))' | pscubecontour n1=9 n2=7 n3=6 d1=.1 d2=.2 d3=.3 faces=1 angle=60 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# Window inside the cube on all three axes (x1beg/x1end, x2beg/x2end, x3beg/x3end interior), so the faces are interpolated off-grid; nc=4 with labelled contours every 2nd line.
python3 -c 'import math,struct,sys; sys.stdout.buffer.write(b"".join(struct.pack("f",math.sin(.5*i1)*math.cos(.4*i2)+.3*i3) for i3 in range(6) for i2 in range(7) for i1 in range(9)))' | pscubecontour n1=9 n2=7 n3=6 x1beg=1 x1end=6 x2beg=1.5 x2end=5 x3beg=0.5 x3end=4 nc=4 labelcper=2 nlabelc=3 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# Explicit front/side/top face files with real (smooth) values and a contour list containing values outside the data range; grid styles on each axis.
python3 -c 'import math,struct; open("front","wb").write(struct.pack("20f",*[math.sin(i*.7) for i in range(20)])); open("side","wb").write(struct.pack("15f",*[math.cos(i*.5) for i in range(15)])); open("top","wb").write(struct.pack("12f",*[i/6-1 for i in range(12)]))' && pscubecontour n1=5 n2=4 n3=3 front=front side=side top=top c=-2,-.5,0,.25,.5,3 grid1=dash grid2=dot grid3=none | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
```

## Additional proposals (ranked below 25 after the metric fixes; all run against both binaries)

### #26 psgraph (score 5.8, 528 C lines, 17 cases / 10 data)

Branches newly exercised: three curves of different lengths with clipping, `style=normal`, and `pairs=0` with large magnitudes. All 3 cases: SAME.

```
# Three curves from stdin with different lengths (n=4,5,6), per-curve marks, and an explicit window that clips lines leaving the box.
python3 -c 'import math,struct,sys; sys.stdout.buffer.write(b"".join(struct.pack("2f",i*0.5,math.sin(i)*3+k) for k,n in enumerate((4,5,6)) for i in range(n)))' | psgraph nplot=3 n=4,5,6 mark=0,2,5 marksize=4 linecolor=red,green,blue x1beg=0.2 x1end=2 x2beg=-2 x2end=3 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# style=normal (never used in a case) with reverse=1 and negative data, automatic axis limits from the data.
python3 -c 'import struct,sys; sys.stdout.buffer.write(struct.pack("10f",-3,-120,-1,40,0.5,-7,2,300,4.25,-55))' | psgraph n=5 style=normal reverse=1 grid1=dot grid2=none | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
# Two curves without pairs (pairs=0: y only, uniform x from d1/f1 per curve) with large magnitudes for the auto tic spacing.
python3 -c 'import struct,sys; sys.stdout.buffer.write(struct.pack("7f",12345,23000,9000.5,-4000,15000,250000,-125000))' | psgraph nplot=2 n=4,3 pairs=0 d1=0.25,2 f1=100,90 linewidth=0,1 mark=3,6 | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'
```

### #27 sutetraray (score 5.6, 3373 C lines, 64 cases / 32 data)

**Coverage hole that cannot be closed:** no golden ever traces through a deeper interface. The existing 2-layer downward case (k=33) segfaults in C (rc 139), and the port copies that crash on purpose with `raise(11)`. Every downward-shooting 2-layer configuration crashed the C.

The cases below avoid those configurations. All 4: SAME.

```
# NOTE: on the 2-layer model, every downward shot that reaches the lower interface (shootupward=-1 with irefseq=1,..; or tmax/takeoff wide enough) segfaults in C (use-after-free of spoiled triangles; the port raises SIGSEGV deliberately, see programs/sutetraray/src/main.rs retire_triangle). Existing golden k=33 pins that crash (rc 139) — so no golden checks traced data through a sub-surface interface. These cases are the non-crashing configurations.
# Multiple sources in x and y (nsx=nsy=2) on the laterally varying 2-layer model: 4 shots written, per-shot table origin logic (ixsf/iysf vs nxt/2).
tetramod nxhz=4 nyhz=3 nhz=2 hzfile=hz.bin v00=1,2 v01=1.1,2.1 v10=1.2,2.2 v11=1.3,2.3 dvdz00=0.1,0.2 dvdz01=0.1,0.2 dvdz10=0.1,0.2 dvdz11=0.1,0.2 z00=0,0.5 z01=0,0.6 z10=0,0.7 z11=0,0.8 > model.txt && sutetraray nxgd=5 nygd=5 nxt=3 nyt=3 ntakeoff=3 nazimuth=8 ntwf=30 nsx=2 nsy=2 fsx=0.5 fsy=0.5 dsx=0.5 dsy=0.5 fsz=0.3 < model.txt
# Source buried in the lower layer (fsz=1.0) shooting up through interface 1 with reflection sequence irefseq=1,0 (transmission at hz1, then surface); ray segments and surface traveltimes dumped.
tetramod nxhz=4 nyhz=3 nhz=2 hzfile=hz.bin v00=1,2 v01=1.1,2.1 v10=1.2,2.2 v11=1.3,2.3 dvdz00=0.1,0.2 dvdz01=0.1,0.2 dvdz10=0.1,0.2 dvdz11=0.1,0.2 z00=0,0.5 z01=0,0.6 z10=0,0.7 z11=0,0.8 > model.txt && sutetraray nxgd=5 nygd=5 nxt=3 nyt=3 ntakeoff=3 nazimuth=8 ntwf=30 fsz=1.0 shootupward=1 irefseq=1,0 rayfile=ray.txt sttfile=s.txt < model.txt > tt.su && cat ray.txt s.txt
# Two sources in depth straddling the 0.5-0.8 interface (nsz=2 dsz=0.6), one per layer, with rayfile dump.
tetramod nxhz=4 nyhz=3 nhz=2 hzfile=hz.bin v00=1,2 v01=1.1,2.1 v10=1.2,2.2 v11=1.3,2.3 dvdz00=0.1,0.2 dvdz01=0.1,0.2 dvdz10=0.1,0.2 dvdz11=0.1,0.2 z00=0,0.5 z01=0,0.6 z10=0,0.7 z11=0,0.8 > model.txt && sutetraray nxgd=5 nygd=5 nxt=3 nyt=3 ntakeoff=3 nazimuth=8 ntwf=30 nsz=2 fsz=0.3 dsz=0.6 rayfile=ray.txt < model.txt > tt.su && cat ray.txt
# tmax/dtwf set instead of ntwf (ntwf derived), non-default traveltime grid (nxgd/nygd/fxgd/fygd/dxgd/dygd) and crfile amplitudes, 2-layer model.
[f32] tetramod nxhz=4 nyhz=3 nhz=2 hzfile=hz.bin v00=1,2 v01=1.1,2.1 v10=1.2,2.2 v11=1.3,2.3 dvdz00=0.1,0.2 dvdz01=0.1,0.2 dvdz10=0.1,0.2 dvdz11=0.1,0.2 z00=0,0.5 z01=0,0.6 z10=0,0.7 z11=0,0.8 > model.txt && sutetraray nxgd=5 nygd=5 nxt=3 nyt=3 ntakeoff=3 nazimuth=8 ntwf=30 nxgd=7 nygd=4 fxgd=0.1 fygd=0.2 dxgd=0.25 dygd=0.4 tmax=0.8 dtwf=0.04 takeoff=50 fsz=0.3 crfile=cr.bin < model.txt > tt.su && cat cr.bin
```

### #28 suwindpoly (score 5.6, 741 C lines, 24 cases / 8 data)

Every existing case feeds all-zero `sunull` traces. These use real samples, polyline mode, multiple polygons and a non-convex polygon. All 3 cases: SAME.

```
# polyline mode on real (non-zero) samples, coordinate keys sx/gx, dx != dy scaling, short-typed ilkey/xlkey (changeval to short), traces output once per near segment
suplane nt=8 ntr=16 | suaddnoise seed=9 | sushw key=sx,gx a=0,0 b=1,0 c=0,1 j=4,4 | suwindpoly xkey=sx ykey=gx x=0,3,3 y=0,1,3 dx=4 dy=2 dw=3 ilkey=nvs xlkey=nhs file=o; cat o
# pass=2 (exterior to all of two polygons) on real data with offset/cdp keys
suplane nt=8 ntr=16 | suaddnoise seed=9 | sushw key=offset,cdp a=0,0 b=1,0 c=0,1 j=4,4 | suwindpoly pass=2 xkey=offset ykey=cdp x=-1,1,1,-1 y=-1,-1,1,1 x=2,4,4 y=2,2,4 file=o; cat o
# pass=1 non-convex polygon: points on a vertex (InPolygon=1), on an edge (=2) and strictly inside (=3)
suplane nt=8 ntr=16 | suaddnoise seed=9 | sushw key=fldr,ep a=0,0 b=1,0 c=0,1 j=4,4 | suwindpoly pass=1 x=0,3,3,1.5,0 y=0,0,3,1,3 file=o; cat o
```

### #35 sutoolcsv (score 5.1, 3358 C lines, 219 cases / 148 data)

These add the SPS1/SPS2 presets, which are the program's main use. All 4 cases: SAME.

**C bug:** with SPS X records the C copies 16-character names into 16-byte buffers, a heap overflow. The port reproduces this on purpose. This case is an EXCLUDE candidate and is not proposed.

```
# SPS2 standard expansion (match=sps2 names=sps2 forms=sps2) on fixed-format S records -> csv: never exercised (no case uses SPS options).
[stderr] printf 'S%10.2f%10.2f  1E1%4d%4.1f%4d%2d%6.1f%9.1f%10.1f%6.1f%3d%06d\n' 1001 2001 5 1.5 10 2 3.0 500100.5 6100200.2 345.6 123 101010 1001 2002 -7 2.0 12 3 0 500125.0 6100225.8 344.1 123 101530 1001 2003 0 0 0 0 0 500150.0 6100250.0 343.0 124 000001 > s.sps && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' sutoolcsv rfile=s.sps wfile=out.csv setid=S match=sps2 names=sps2 forms=sps2 && cat out.csv
# Round trip csv -> fixed with all options read back from C_SU_* records (no setid/match/names/forms on the command line); exercises C_SU_MORE continuation on fixed output.
[stderr] printf 'S%10.2f%10.2f  1E1%4d%4.1f%4d%2d%6.1f%9.1f%10.1f%6.1f%3d%06d\n' 1001 2001 5 1.5 10 2 3.0 500100.5 6100200.2 345.6 123 101010 1001 2002 -7 2.0 12 3 0 500125.0 6100225.8 344.1 123 101530 1001 2003 0 0 0 0 0 500150.0 6100250.0 343.0 124 000001 > s.sps && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' sutoolcsv rfile=s.sps wfile=out.csv setid=S match=sps2 names=sps2 forms=sps2 && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' sutoolcsv rfile=out.csv wfile=back.txt && cat out.csv back.txt
# names=sps2all: the numbN (non-key passthrough) fields of S records, fixed -> fixed.
[stderr] printf 'S%10.2f%10.2f  1E1%4d%4.1f%4d%2d%6.1f%9.1f%10.1f%6.1f%3d%06d\n' 1001 2001 5 1.5 10 2 3.0 500100.5 6100200.2 345.6 123 101010 1001 2002 -7 2.0 12 3 0 500125.0 6100225.8 344.1 123 101530 1001 2003 0 0 0 0 0 500150.0 6100250.0 343.0 124 000001 > s.sps && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' sutoolcsv rfile=s.sps wfile=out.txt setid=S match=sps2 names=sps2all forms=sps2 && cat out.txt
# SPS1 R-record expansion (grnors/gaps match, sps1 formats) with duplicate-free multi-line receivers.
[stderr] printf 'R%10.2f%10.2f  1G1%4d%4.1f%4d%2d%6.1f%9.1f%10.1f%6.1f%3d%06d\n' 3001 4001 0 0 0 0 0 600000.0 6200000.0 100.0 200 120000 3001 4002 4 0 0 0 0 600025.0 6200000.0 101.5 200 120100 3002 4001 -3 0 0 0 0 600000.0 6200050.0 99.0 200 120200 > r.sps && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' sutoolcsv rfile=r.sps wfile=out.csv setid=R match=sps1 names=sps1 forms=sps1 && cat out.csv
# NOT PROPOSED (C heap overflow): SPS2/SPS1 X presets strcpy 16-char names into ealloc1(16) buffers; on macOS the C then fails "names= list has an entry that parses incorrectly". The port reproduces this deliberately (names[8..=11].concat()). Kept here for the record only:
#   [stderr] printf 'X%6d%8d11%10.2f%10.2f1%5d%5d1%10.2f%10.2f%10.2f1\n' 1 77 1001 2001 1 3 3001 4001 4003 1 78 1001 2002 1 2 3002 4001 4002 > x.sps && python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' sutoolcsv rfile=x.sps wfile=out.csv setid=X match=sps2 names=sps2 forms=sps2 && cat out.csv
```

### #71 sufarldcsv (score 3.2, 2083 C lines, 141 cases / 46 data)

Branches newly exercised: `cseq=gx`, `cseq=sx` with 2-D locations, duplicate cseq values, and `asis` with two-valued lists. All 4 cases: SAME.

```
# cseq=<key> other than cdp/asis: Q-records re-sorted by gx (compSort1 on a different column), never exercised.
[stderr] python3 -c 'import os,pty,subprocess,sys,pathlib; q="C_SU_NAMES\nC_SU_ID,cdp,gx,gy,sx,sy,sstat,gstat\n"+"".join("Q,%d,%d,%d,%d,%d,%d,%d\n"%(i+1,i*7,(i%3)*2,i*7+3,(i%2)*3,(i*i)%17-6,i*3-9) for i in [4,0,6,2,5,1,3]); [pathlib.Path(n).write_text(q) for n in ["c.csv","sstat.csv","rstat.csv"]]; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); out=pathlib.Path(next((a[5:] for a in sys.argv[1:] if a.startswith("dout=")),"tstat.csv")); sys.stdout.buffer.write(out.read_bytes() if out.is_file() else b""); sys.exit(p.returncode)' sufarldcsv cin=c.csv cseq=gx navrg=2 nsmth=1 nback=1
# cseq=sx with 2-D cloc (gx,gy): sequencing column differs from both location columns.
[stderr] python3 -c 'import os,pty,subprocess,sys,pathlib; q="C_SU_NAMES\nC_SU_ID,cdp,gx,gy,sx,sy,sstat,gstat\n"+"".join("Q,%d,%d,%d,%d,%d,%d,%d\n"%(i+1,i*7,(i%3)*2,i*7+3,(i%2)*3,(i*i)%17-6,i*3-9) for i in [4,0,6,2,5,1,3]); [pathlib.Path(n).write_text(q) for n in ["c.csv","sstat.csv","rstat.csv"]]; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); out=pathlib.Path(next((a[5:] for a in sys.argv[1:] if a.startswith("dout=")),"tstat.csv")); sys.stdout.buffer.write(out.read_bytes() if out.is_file() else b""); sys.exit(p.returncode)' sufarldcsv cin=c.csv cseq=sx cloc=gx,gy sloc=sx,sy rloc=gx,gy navrg=1 nsmth=0 nback=0
# Duplicate cseq value -> "Two sets of values input for cin cseq" error branch.
[stderr] python3 -c 'import os,pty,subprocess,sys,pathlib; q="C_SU_NAMES\nC_SU_ID,cdp,gx,gy,sx,sy,sstat,gstat\n"+"".join("Q,%d,%d,%d,%d,%d,%d,%d\n"%(i+1,i*7,(i%3)*2,i*7+3,(i%2)*3,(i*i)%17-6,i*3-9) for i in [4,0,6,2,5,1,3,2]); [pathlib.Path(n).write_text(q) for n in ["c.csv","sstat.csv","rstat.csv"]]; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); out=pathlib.Path(next((a[5:] for a in sys.argv[1:] if a.startswith("dout=")),"tstat.csv")); sys.stdout.buffer.write(out.read_bytes() if out.is_file() else b""); sys.exit(p.returncode)' sufarldcsv cin=c.csv navrg=0 nsmth=0 nback=0
# cseq=asis with nback/nsmth/navrg all 2-valued lists (separate first/second values) on the default 2-D grid-free path.
[stderr] python3 -c 'import os,pty,subprocess,sys,pathlib; q="C_SU_NAMES\nC_SU_ID,cdp,gx,gy,sx,sy,sstat,gstat\n"+"".join("Q,%d,%d,%d,%d,%d,%d,%d\n"%(i+1,i*7,(i%3)*2,i*7+3,(i%2)*3,(i*i)%17-6,i*3-9) for i in [4,0,6,2,5,1,3]); [pathlib.Path(n).write_text(q) for n in ["c.csv","sstat.csv","rstat.csv"]]; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); out=pathlib.Path(next((a[5:] for a in sys.argv[1:] if a.startswith("dout=")),"tstat.csv")); sys.stdout.buffer.write(out.read_bytes() if out.is_file() else b""); sys.exit(p.returncode)' sufarldcsv cin=c.csv cseq=asis navrg=1,2 nsmth=2,1 nback=1,2 kdist=3 kmult=1.5
```

### #166 sumigtopo2d (score 1.7, 1149 C lines, 89 cases / 27 data)

Branches newly exercised: realistic shot gathers binned into offset panels, `offmax` with dvz, off-node receivers, and a fine output grid with `ls=0`. All 4 cases: SAME.

```
# Realistic shot gathers (2 shots x 4 receivers, split-spread) binned into 3 offset panels; receivers outside [fs,es] (gx=-6, gx=40) are skipped while the rest migrate.
[stderr] python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(x-s,z)/1600.0) for s in (0,10,20,30) for x in (0,10,20,30) for z in (5,12,19,26)))'; suspike nt=81 ntr=8 dt=0.001 nspk=2 ix1=2 it1=20 ix2=7 it2=35 | sushw key=sx,gx,d2 a=4,-6,10 b=0,10,0 c=16,16,0 j=4,4,4 >in; sumigtopo2d ttfile=tt infile=in outfile=out nxt=4 fxt=0 dxt=10 nzt=4 fzt=5 dzt=7 ns=4 fs=0 ds=10 nxi=7 fxi=0 dxi=5 v0=1600 'surf=-10,0;40,0' noff=3 off0=-10 doff=10 && cat out
# offmax rejects the far traces (|offset|>12) while near ones still migrate; topography + positive dvz.
[stderr] python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(x-s,z)/1600.0) for s in (0,10,20,30) for x in (0,10,20,30) for z in (5,12,19,26)))'; suplane nt=81 ntr=6 dt=0.001 | suaddnoise seed=33 sn=5 | sushw key=sx,gx,d2 a=3,27,10 b=4,-4,0 >in; sumigtopo2d ttfile=tt infile=in outfile=out nxt=4 fxt=0 dxt=10 nzt=4 fzt=5 dzt=7 ns=4 fs=0 ds=10 nxi=7 fxi=0 dxi=5 v0=1600 'surf=-10,1;8,-1;22,2;40,0' offmax=12 dvz=0.4 noff=2 off0=-8 doff=16 && cat out
# Receivers 1 m off table nodes (sx=10, gx=19,20,21) give res=0.9/0.1 in both the traveltime-table (ds=10) and surface (dxi=5) interpolation, plus an on-node trace.
[stderr] python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(x-s,z)/1600.0) for s in (0,10,20,30) for x in (0,10,20,30) for z in (5,12,19,26)))'; suspike nt=81 ntr=3 dt=0.001 nspk=1 ix1=2 it1=25 | sushw key=sx,gx,d2 a=10,19,10 b=0,1,0 >in; sumigtopo2d ttfile=tt infile=in outfile=out nxt=4 fxt=0 dxt=10 nzt=4 fzt=5 dzt=7 ns=4 fs=0 ds=10 nxi=7 fxi=0 dxi=5 v0=1600 'surf=-10,0;40,0' && cat out
# Large output depth grid finer than the table (nzo=25, dzo=0.875) with ls=0 and a narrow angle aperture.
[stderr] python3 -c 'import math,struct; open("tt","wb").write(b"".join(struct.pack("=f",math.hypot(x-s,z)/1600.0) for s in (0,10,20,30) for x in (0,10,20,30) for z in (5,12,19,26)))'; suspike nt=81 ntr=3 dt=0.001 nspk=1 ix1=2 it1=25 | sushw key=sx,gx,d2 a=5,15,10 b=5,0,0 >in; sumigtopo2d ttfile=tt infile=in outfile=out nxt=4 fxt=0 dxt=10 nzt=4 fzt=5 dzt=7 ns=4 fs=0 ds=10 nxi=7 fxi=0 dxi=5 v0=1600 'surf=-10,0;40,0' nzo=25 fzo=5 dzo=0.875 nxo=13 fxo=0 dxo=2.5 ls=0 angmax=25 && cat out
```

### #226 susmoqcsv (score 1.3, 1444 C lines, 153 cases / 67 data)

These add longer unordered lines and 3-D grids, using the default navrg/nsmth/nback, `keyloc=aux`, and `outlist`/`extrapt`. All 3 cases: SAME.

```
# 2D line of 25 unordered records with default navrg=10/nsmth=10/nback=3 (full-width averaging + smoothing + end extrapolation), tuples smoothed, aux copied
[exact] [stderr] python3 -c 'import math; print("C_SU_NAMES"); print("C_SU_ID,cdp,aux,vel,eta,vel,eta"); print("C_SU_NDIMS,tims,null,0,0,1500,1500"); [print("Q,%d,%g,%g,%g,%g,%g" % (c, (c*7)%11, 1500+40*math.sin(c), 0.1+0.01*c, 2500+3*c*c, 0.3-0.02*(c%5))) for c in range(25,0,-1)]' > in.q; python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); sys.exit(p.returncode)' susmoqcsv qin=in.q qout=out.q formtv=%.4f && cat out.q
# Same line: keyloc=aux (sorted by a non-cdp name with duplicates), sname=aux, nback=0 (duplicate end padding), outind resampling of tuples with extrapt=1
[exact] [stderr] python3 -c 'import math; print("C_SU_NAMES"); print("C_SU_ID,cdp,aux,vel,eta,vel,eta"); print("C_SU_NDIMS,tims,null,0,0,1500,1500"); [print("Q,%d,%g,%g,%g,%g,%g" % (c, (c*7)%11, 1500+40*math.sin(c), 0.1+0.01*c, 2500+3*c*c, 0.3-0.02*(c%5))) for c in range(25,0,-1)]' > in.q; python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); sys.exit(p.returncode)' susmoqcsv qin=in.q qout=out.q keyloc=aux sname=aux navrg=3 nsmth=2 nback=0 outind=0,3000,750 extrapt=1 formtv=%.4f && cat out.q
# 3D 7x5 grid, cdp-located records: navrg=-1 inline (whole-line mean) with crossline smoothing, cdpt=1 transposed output, cname=eta copy, outlist with extrapt=3
[exact] [stderr] printf %b 'C_SU_SETID,K\nC_SU_FORMS\nC_SU_ID,%g,%g,%g,%g,%g,%g,%g,%g\nC_SU_NAMES\nC_SU_ID,grid_xa,grid_ya,grid_xb,grid_yb,grid_xc,grid_yc,grid_wb,grid_wc\nK,1000,2000,1090,2120,920,2060,25,25\n' > grid.k; python3 -c 'print("C_SU_NAMES"); print("C_SU_ID,cdp,aux,vel,eta,vel,eta"); print("C_SU_NDIMS,tims,null,0,0,1000,1000"); [print("Q,%d,%d,%d,%g,%d,%g" % (c, c*c%13, 1500+11*c, 0.05*c, 2400-7*c, 0.2+0.03*(c%4))) for c in range(1,36)]' > in.q; python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(p.stderr); sys.exit(p.returncode)' susmoqcsv qin=in.q qout=out.q rfile=grid.k navrg=-1,1 nsmth=0,1 nback=1,2 cdpt=1 cname=eta outlist=1500,0,2500 extrapt=3 formtv=%.5g formxy=%.1f && cat out.q
```

### #229 subinqcsv (score 1.2, 1574 C lines, 163 cases / 82 data)

These add 2-D `keyloc=gaps`, re-input through `qin=`, and `ename`/`invp2`. All 3 cases: SAME.

```
# 2D keyloc other than cdp/fldr (gaps) with three locations and per-location tuples of different time samples; extrapolation beyond input times (extrapt) and regular output grid outind.
[stderr] python3 -c 'import os,subprocess,sys;a,b=os.openpty();sys.exit(subprocess.call(sys.argv[1:],stdin=b,stdout=b))' subinqcsv qout=out.csv keyloc=gaps gaps=10,20,40 tims=0,800,2000 vnmo=1500,1900,2600 tims=0,1000,2500 vnmo=1550,2000,2800 tims=0,1200,3000 vnmo=1600,2100,3000 outind=0,3500,500 extrapt=1 && cat out.csv
# Re-input of a Q-file written by this program (qin= path, C_SU records read back) with an irregular output time list (outlist) and formtv=.
[stderr] python3 -c 'import os,subprocess,sys;a,b=os.openpty();sys.exit(subprocess.call(sys.argv[1:],stdin=b,stdout=b))' subinqcsv qout=out.csv cdp=1,5,9 tims=0,1000,3000 vels=1500,2000,3000 vels=1600,2100,3100 vels=1700,2200,3200 outind=0,3000,1000 && python3 -c 'import os,subprocess,sys;a,b=os.openpty();sys.exit(subprocess.call(sys.argv[1:],stdin=b,stdout=b))' subinqcsv qin=out.csv qout=out2.csv outlist=250,750,1500,2500 formtv=%.5g && cat out2.csv
# invp2=tuple with ename removing one tuple name on qin input.
[stderr] python3 -c 'import os,subprocess,sys;a,b=os.openpty();sys.exit(subprocess.call(sys.argv[1:],stdin=b,stdout=b))' subinqcsv qout=out.csv cdp=1,5 tims=0,1000,3000 vels=1500,2000,3000 vnmo=1400,1900,2900 vels=1600,2100,3100 vnmo=1450,1950,2950 outind=0,3000,500 && python3 -c 'import os,subprocess,sys;a,b=os.openpty();sys.exit(subprocess.call(sys.argv[1:],stdin=b,stdout=b))' subinqcsv qin=out.csv qout=out2.csv ename=vnmo invp2=tuple outind=0,3000,1000 && cat out2.csv
```

### #313 suquantile (score 0.8, 236 C lines, 56 cases / 38 data)

**Metric false positive.** Its report goes to stderr and all 56 existing cases are `[stderr]` cases; the metric was fixed to count them. The cases below cover large n, heavy ties, and NaN/inf input. All 6: SAME.

```
# Whole-panel quantiles on 10000 dense noisy samples (large-n qsort path, many percentile indices > 0).
[stderr] suplane nt=500 ntr=20 | suaddnoise seed=3 | suquantile
# Panel rank mode where ~95% of samples tie at zero: reported min/max sample indices depend on libc qsort's ordering of equal keys.
[stderr] suplane nt=301 ntr=12 | suquantile quantiles=0
# Trace mode on 1000-sample traces of only -1/0/+1 values (massive ties) in both rank and quantile mode.
[stderr] suplane nt=1000 ntr=2 | suaddnoise seed=7 | suop op=sgn | suquantile panel=0 quantiles=0
[stderr] suplane nt=1000 ntr=2 | suaddnoise seed=7 | suop op=sgn | suquantile panel=0
# NaN and +-inf samples (built natively, suaddhead): cmp_indirect returns 0 for NaN (and inf-inf), an inconsistent qsort comparator; %8.2e of nan/inf.
[stderr] python3 -c 'import struct,sys; v=[float(i%7)-3 for i in range(48)]; v[5]=float("nan"); v[17]=float("inf"); v[30]=float("-inf"); v[40]=float("nan"); sys.stdout.buffer.write(struct.pack("=48f",*v))' | suaddhead ns=16 | suquantile
[stderr] python3 -c 'import struct,sys; v=[float(i%7)-3 for i in range(48)]; v[5]=float("nan"); v[17]=float("inf"); v[30]=float("-inf"); v[40]=float("nan"); sys.stdout.buffer.write(struct.pack("=48f",*v))' | suaddhead ns=16 | suquantile panel=0 quantiles=0
```

### #350 surescsv (score 0.6, 1462 C lines, 88 cases / 44 data)

Most existing data cases have zero static shifts. These apply known source and receiver statics with `sustatic`, so the correlation-peak, interpolation and restack code really runs.

5 cases: SAME. The last case is the text-form divergence (Divergence 2).

```
# Known source (-8..8 ms) and receiver (4..-6 ms) statics applied by sustatic hdrs=1 to a band-limited flat event: exercises real cross-correlation peaks, 3-point interpolation and ints8r fractional restack (existing data cases are mostly sunull/suplane with zero shifts). No long-wave reduction, all offsets.
[f32] [stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=3 longwave=-1 maxshift=12 cdivider=1 useoff=1000 && python3 -c 'import csv,struct,sys; values=[float(row[-1]) for name in ("sstat.csv","rstat.csv") for row in csv.reader(open(name)) if row and row[0]=="Q"]; sys.stdout.buffer.write(struct.pack("<%df"%len(values),*values))'
# Same data with defaults: long-wavelength reduction (longwave = average offset), cdivider=4 bundling and default useoff (1.5 x average offset excludes far traces).
[f32] [stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=3 maxshift=12 && python3 -c 'import csv,struct,sys; values=[float(row[-1]) for name in ("sstat.csv","rstat.csv") for row in csv.reader(open(name)) if row and row[0]=="Q"]; sys.stdout.buffer.write(struct.pack("<%df"%len(values),*values))'
# Offset limiting switched on only from iteration 2 (iuseoff=1, useoff=150) with maxshift=5 ms (just above the 4 ms sample) limiting the correlation lag window; result only partly converged, so sensitive to iteration order.
[f32] [stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=4 iuseoff=1 useoff=150 maxshift=5 longwave=-1 cdivider=1 && python3 -c 'import csv,struct,sys; values=[float(row[-1]) for name in ("sstat.csv","rstat.csv") for row in csv.reader(open(name)) if row and row[0]=="Q"]; sys.stdout.buffer.write(struct.pack("<%df"%len(values),*values))'
# Sources-only and receivers-only derivation on data with real shifts (srderive=-1 / 1).
[f32] [stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=3 srderive=-1 longwave=-1 cdivider=1 useoff=1000 && python3 -c 'import csv,struct,sys; values=[float(row[-1]) for name in ("sstat.csv","rstat.csv") for row in csv.reader(open(name)) if row and row[0]=="Q"]; sys.stdout.buffer.write(struct.pack("<%df"%len(values),*values))'
[f32] [stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=3 srderive=1 longwave=300 cdivider=2 useoff=1000 && python3 -c 'import csv,struct,sys; values=[float(row[-1]) for name in ("sstat.csv","rstat.csv") for row in csv.reader(open(name)) if row and row[0]=="Q"]; sys.stdout.buffer.write(struct.pack("<%df"%len(values),*values))'
# Text form of case 1 (golden style `cat sstat.csv rstat.csv`, %.12g): exposes the qsort-vs-stable-sort ordering difference (last-digit f32 differences => exact text mismatch).
[stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=3 longwave=-1 maxshift=12 cdivider=1 useoff=1000 && cat sstat.csv rstat.csv
```

## Divergences found

All three were reproduced with `runcases.py`. The C output was stable across two runs in each case.

### 1. suprofcsv: the Rust port is an unimplemented stub (real bug, high severity)

- `programs/suprofcsv/src/main.rs` is 285 lines against 1029 in the C. After validating its parameters it ends with:
  ```rust
  efclose(&fp);
  err!("suprofcsv: profile generation is not yet implemented");
  ```
- It passes its golden tests only because every golden case is an error exit. The 2 inventory cases fail in C too: they run without a pty, and the C refuses that.
- All 5 cases in `cases/suprofcsv.txt` give: exit code C 0, Rust 1.

Command (needs the pty wrapper because the C errs when stdin is not a tty):
```
[stderr] python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); r=subprocess.run(sys.argv[1:],stdin=s,stdout=s,stderr=subprocess.PIPE); os.close(s); os.close(m); sys.stderr.buffer.write(r.stderr); sys.exit(r.returncode)' suprofcsv gaps=1,2,3,4,5,6,7,8 gx=0,25,50,75,100,124,147,168 gy=0,0,1,3,6,12,20,30 qout=p.q nmaxa=2 nmaxs=1 chordi=12.5 formxy=%.6f && cat p.q
```
| | C | Rust |
|---|---|---|
| Exit code | 0 | 1 |
| stdout | 941 bytes of Q-file | empty |
| stderr | empty | `suprofcsv: suprofcsv: profile generation is not yet implemented` |

A grep of all programs for "not yet implemented", `unimplemented!` and `todo!` found no other stub. The other hits (sucwt, sudmovz, sudmotivz) are faithful copies of C error messages.

### 2. surescsv: `qsort` order versus stable sort (harness-visible text divergence)

- The C sorts its source, receiver and CDP index arrays with `qsort(..., comparenn)`, which compares only the key. macOS `qsort` reorders equal keys.
- The Rust (`programs/surescsv/src/main.rs:550-552`) uses `sndx/rndx/cndx.sort_by(comparenn)`, which is stable. Its comment says "Preserve input order for ties".
- So traces within one source, receiver or CDP are summed in a different order in the CDP stack and in the cross-correlation sums. The derived statics then differ in the last float digit, and `%.12g` prints the difference.

Command (last line of `cases/surescsv.txt`):
```
[stderr] suplane npl=1 dip1=0 ct1=50 len1=30 nt=101 ntr=30 | sufilter f=5,10,30,40 | sushw key=fldr,gaps,cdp,offset,sx,gx,sstat,gstat a=1,1,2,0,0,0,-8,4 b=0,1,1,100,0,100,0,-2 c=1,0,1,-100,100,0,4,0 j=6,6,6,6,6,6,6,6 | sustatic hdrs=1 | sushw key=sstat,gstat a=0,0 | python3 -c 'import os,pty,subprocess,sys; m,s=pty.openpty(); p=subprocess.run(sys.argv[1:],stdin=sys.stdin.buffer,stdout=s); os.close(s); os.close(m); sys.exit(p.returncode)' surescsv numiter=3 longwave=-1 maxshift=12 cdivider=1 useoff=1000 && cat sstat.csv rstat.csv
```
Output: rc 0 on both sides. `text line 9: golden 'Q,2,100,0,0,6,-4.08975839615' actual 'Q,2,100,0,0,6,-4.08975887299'`. The same happens with `numiter=1` and with `srderive=-1`.

- **Size:** about 1.2e-7 relative, so the `[f32]` versions of the same runs pass.
- **Why it still matters:** any golden that `cat`s the CSV on data with non-zero statics fails. Existing goldens avoid this only because their statics are zero or extracted as f32.
- **Fix:** call libc `qsort` through FFI, as susort, suquantile, suedit and sutoolcsv already do.

### 3. sufctanismod: uninitialised memory in the C (not a port bug; affects existing goldens)

Command:
```
[f32] sufctanismod nx=31 nz=33 nt=60 sx=15 sz=10 source=1 indexux=1 indexuy=1 indexuz=1 wavelet=1 dofct=1 isurf=1 aa00=3 cc00=2.2 ff00=1.1 ll00=0.8 nn00=1.2 daadz=0.5 dccdx=0.3 dffdz=-0.2 dlldx=0.1 dnndz=0.2 drhodz=0.4 > final; cat snapshotx.data snapshoty.data snapshotz.data final
```
Output: rc 0 on both sides. `float32 sample 60854: golden -0.00468435045 actual -0.00512528187 |diff| 0.000441 > allowed 0.000221`.

**Smaller repro:** `sufctanismod nx=31 nz=33 nt=60 sx=15 sz=10 source=1 indexux=1 wavelet=4 dofct=1 isurf=1`
- 0.84% maximum relative error; 70 of 1023 samples out of tolerance.
- The error grows with the number of steps: about 1e-14 at nt=20, 4e-10 at nt=30, 1e-6 at nt=40, and 8e-3 at nt=60.
- `dofct=0`, `isurf=2` and `isurf=3` each give zero difference.

**Cause:**
- `anis_solver2` allocates `u1`/`f0`/`f1` with `alloc2float` on every call and never zeroes them.
- `difference()` writes only interior cells.
- `fct2d1o` then reads boundary cells of `u1` that nothing has written, so it sees whatever the previous step left in that heap block.
- The port zero-fills these arrays (`programs/sufctanismod/src/solver.rs:89`, `Grid::new`).

**Evidence:** C run with `MallocPreScribble=1` matches the Rust to 3.9e-39, while the plain C run differs by 8.4e-3.

**Recommendation:** treat this like the existing wpc EXCLUDE entries. Four existing `dofct=1` goldens pass only because the stale values stay within tolerance:
- the `order=4 dofct=1` case;
- the `deta*`/`fct*` case;
- the `movebc=1 impulse=1 dofct=1` case;
- the eta-gradient case.

The port's zero-initialised result is the defensible one.

## Other findings

- **Exit-code-only cases:** 197 cases in 84 programs. The ones that matter most are listed below; all are listed in `metrics_final.json` under `redirect_only`. 169 come from the inventory `test_cases`. The fix is to drop the redirect, or append `&& cat <file>`.

  | Program | Cases | Program | Cases |
  |---|---|---|---|
  | segyread | 8 | unisam | 4 |
  | suwellrf | 6 | velconv | 4 |
  | pdfhistogram | 4 | suwaveform | 4 |
  | suvibro | 4 | sugetgthr | 3 |
  | susynvxzcs | 3 | supack1/2, suunpack1/2 | 3 each |
  | sutihaledmo | 3 | supswigb/supswigp | 3 each |
  | sudmofk, sudmotx, sudmovz, sudmotivz | 2 each | sumigfd | 2 |
  | sucdpbin | 2 | sustrip | 2 |
  | susort | 2 | | |
- **suwaveform:** all 4 of its inventory cases redirect to `out.su`. Its cases.txt cases do check output.
- **sureflpsvsh and supolar** were still being ported during the audit. Re-run `coverage.py` once they land.
- **sufilter** has no recorded goldens (`manifest.json`) in this checkout, although its cases.txt exists.
- **C undefined behaviour or crashes to keep out of goldens** (the port behaves sanely in each):

  | Program | Problem in the C |
  |---|---|
  | sumute | `mode=4` uses an uninitialised `twfval` |
  | sumedian | `sign=1` with non-zero shifts segfaults; reads `inshift[ntr]` at the end of the gather |
  | suvelan_nsel | multi-CDP input writes past `aoff`/`datm` |
  | sutoolcsv | SPS X presets overflow a heap buffer |
  | sutetraray | downward 2-layer runs segfault |
  | suremel2dan | receiver/snapshot files trap (rc 133) |
  | sureflpsvsh | NaN output and SH segfault |
  | sufctanismod | `movebc=1` sporadic SIGBUS |
  | sukdmig3d | DEBUG block forces sy=gy=3 |

## Top 10 weakest (final ranking)

1. sureflpsvsh (still being ported: no binary, no goldens; broken C oracle)
2. sudatumk2dr
3. sugain
4. **suprofcsv (stub port: Divergence 1)**
5. segywrite
6. sufctanismod (C undefined behaviour: Divergence 3)
7. suremel2dan
8. unif2aniso
9. suharlan
10. sugeomcsv

**Run against both binaries:**
- Ranks 2–10, plus sumatrix (#11) in place of sureflpsvsh, which has no binary. sureflpsvsh was run in C only.
- All other top-25 programs.
- 10 more programs: sumigtopo2d, sutoolcsv, sufarldcsv, surescsv (Divergence 2), susmoqcsv, subinqcsv, suquantile, sutetraray, suwindpoly, psgraph. They ranked in the top 25 under earlier versions of the metric, before the quote, wrapper and stderr parsing fixes.
