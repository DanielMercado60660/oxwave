# Audit (misc): spsplot, succepstrum, sufxdecon, sumigpspi, sust

Auditor: Claude (read-only on the repo; experiments in /private/tmp/claude-501/audit3).
Method: I read each C source and its Rust port side by side, including the library routines they
call (phase_unwrap, LU_decomposition, cwp_cexp, pfacc). I ran tools/check_program.py for all five
programs. I then ran differential tests (original SU vs target/release binaries, harness tolerance)
on new synthetic inputs, and built ASan copies of the C programs (with the relevant cwp/lib
sources instrumented) to check the claims of undefined behaviour.

**Summary:** 1 high, 5 medium, 5 low.

All the numerics I checked match: float/double promotion, NINT, literal types, loop order and
FFT conventions in all five programs. Every new realistic run that was not on an undefined-behaviour
path matched the original within tolerance:
- sumigpspi: a dipping reflector with constant velocity, a laterally varying gradient model, a
  makevel low-velocity lens, and random blocky velocities (PSPI branch).
- succepstrum: about 130 trace lengths, both phase modes, up to nt=3001.
- sust: about 50 trace lengths, up to nt=1500.
- sufxdecon: window and filter grid.
- spsplot: a three-layer trimodel model, 394 shaded strips, three parameter sets.

All five files start with tools/cwp_header.txt verbatim. No `rtol`/`atol` loosening appears in
these programs' cases.txt.

## Parameter coverage

| program | selfdoc params | exercised by goldens | important untested branch (now tested) |
|---|---|---|---|
| spsplot | gedge gtri gmin gmax sgmin sgmax labelx labelz title style labelfont titlefont labelsize titlesize titlecolor axescolor gridcolor xbox ybox wbox hbox xbeg xend zbeg zend dxnum fxnum nxtic gridx dznum fznum nztic gridz axeswidth ticwidth gridwidth | all (on one base64 2-face model) | real trimodel multi-layer model: matches (L4) |
| succepstrum | sign1 sign2 mode unwrap trend zeromean | all | many trace lengths: match except the pfacc-UB lengths |
| sufxdecon | taper fmin fmax twlen ntrw ntrf verbose tmpdir (+dt) | all | many time windows (ntw ≳ 25): **panic (H1)**; ntrf>2·ntrw: panic (M2); negative sizes: panic (M1) |
| sumigpspi | nz dz vfile dt dx tmpdir verbose | all | tmpdir root-owned/non-writable: wrong message (M3) |
| sust | dt fmin fmax verbose | all | frequency bins outside [0,nfft]: invented error (M4); lo=hi=0: C segfaults (L1) |

## Excluded golden cases (succepstrum, sust): the undefined-behaviour claim checks out

- **The UB is real, and it is in pfacc.** I rebuilt sust.c and succepstrum.c with `-fsanitize=address`,
  compiling `cwp/lib/pfafft.c` instrumented. Every excluded length ends in a heap-buffer-overflow
  READ. Frame #0 falls inside `_pfacc` (between it and the next symbol, `_pfacr`) for:
  - sust nt=95, 96, 212, 282, 2000;
  - succepstrum nt=46, 48, 121, 126, 318.

  The cause is in pfacc's mu loop: when `m` has no inverse mod `ifac`, the loop runs to the end,
  leaving `mm=n` and `j01=2n`. The original also *crashes* on some of these lengths (succepstrum
  nt=47: SIGTRAP, rc 133).
- **Excluding these cases did not hide a Rust bug.** I swept succepstrum (ouphase nt=2..200 step 3,
  suphase 3..330 step 5, plus 501–3001) and sust (2..300 step 7, plus 500–1500) with nonzero data.
  The only differences were the rejected lengths. Every rejection is a clean `err!`
  (`pfafft.c: pfacc: invalid transform length`, rc 1), never a panic.
- **Impact.** Of trace lengths 2..600 plus a few large ones, 30 of
  608 are rejected by each program:
  - sust: 91–96, 211–224, 281–288, **2000, 2001**;
  - succepstrum: 46–48, 121–126, 316–336.

  nt=2000/2001 is a common record length, so OxWave sust refuses some everyday inputs.
- **Note on the "control" cases.** Even sust nt=64 (default) trips ASan. That is a different overflow:
  a WRITE 0 bytes past the 9216-byte `vgr` block (36 rows × 64 × 4). It happens because `st()`
  writes `nfft` values per row while `ealloc2float(nt,nf)` provides `nt`. It occurs on every
  run with nfft>nt. It does not affect the output, and the Rust port emulates it safely
  (programs/sust/src/main.rs:156-168). So every retained sust golden also comes from a C run with
  heap corruption.
- **The goldens are stale.** See M5.

---

## H1 (high): sufxdecon panics on long traces with short time windows

- Rust: programs/sufxdecon/src/main.rs:335-336 (reads `tdatain`); writes at :525, :535, :553 (`tdataout`).
- C: sufxdecon.c:329 (read), :506/:518 (writes).
- **What differs.** When `ntw = NINT((nsp-1)*dt/twlen)` rounds up and there are many windows, the
  last *intermediate* window (iw = ntw-2, `nspw = nspwi`) extends past `nsp`. The index is
  `it + iw*nspws - nstaper/2`, which reaches `(ntw-1)*nspws + nstaper/2 - 1 ≥ nsp`.
  - For traces other than the last, both C and Rust read the next trace's leading samples. They
    agree because both use a flat array.
  - For the **last trace**, C reads (and later writes) a few floats past the `ealloc2float` block
    and still produces normal output. Rust indexes past the Vec and panics (rc 101, Rust panic
    text instead of output).
  - The defaults (dt=4 ms, taper=.1) with twlen=0.3 trigger this for nt in 2064, 2139–2140,
    2214–2216, …, 2964–2976 (91 of the lengths 1000–3000). These are 8–12 s records with
    0.3 s windows, a normal f-x decon setup. twlen=0.5 or 1.0 does not trigger it below nt=3000.
- **Repro.**
  ```
  suplane nt=2064 ntr=12 | suaddnoise seed=3 sn=4 | sufxdecon twlen=0.3
  SU: rc 0, 12 traces
  OX: rc 101  thread 'main' panicked at programs/sufxdecon/src/main.rs:336:32:
      index out of bounds: the len is 24768 but the index is 24768
  ASan C build: heap-buffer-overflow READ 0 bytes after 99072-byte region (= 12*2064*4)
  ```
  The C output is deterministic: md5 is identical under environment-size perturbation. A patched
  C build that allocates one extra **zeroed** row for `tdatain`/`tdataout` produces output within
  the harness tolerance of the original for nt=2064 and nt=2976 (`compare_stdout` → None).
- **Suggested fix.** Allocate `tdatain` and `tdataout` with `nspwi` (or `nsp`) zeroed slack floats
  past `ntrc*nsp`, output only the first `ntrc*nsp`, and add a comment citing C's over-read
  (DESIGN rule 5). If this is treated as UB instead, at minimum replace the panic with a
  clean `err!`.
- **Golden case.**
  `suplane nt=2064 ntr=12 | suaddnoise seed=3 sn=4 | sufxdecon twlen=0.3`
  and `suplane nt=2976 ntr=12 | suaddnoise seed=3 sn=4 | sufxdecon twlen=0.3`

## M1 (medium): sufxdecon panics on negative ntrw/ntrf instead of C's allocation error

- Rust: programs/sufxdecon/src/main.rs:279-301 (`vec![..; (2*ntrw) as usize]`, `(2*ntrf) as usize`, …).
- C: sufxdecon.c:266-291 (`ealloc*`).
- **What differs.** C fails in `ealloc` with `ealloc.c: malloc failed (Cannot allocate memory)`, rc 1.
  Rust panics with `capacity overflow`, rc 101.
- **Repro.**
  ```
  suplane nt=251 ntr=20 | sufxdecon ntrw=-3   -> SU rc1 "sufxdecon: ealloc.c: malloc failed (Cannot allocate memory)"; OX rc101 panic "capacity overflow"
  suplane nt=251 ntr=20 | sufxdecon ntrf=-1   -> same
  ```
- **Suggested fix.** Check the sizes before allocating (as sust does at programs/sust/src/main.rs:141-143) and emit
  the same `err!`.
- **Golden cases.** `[stderr] suplane nt=101 ntr=2 | sufxdecon ntrw=-3` and
  `[stderr] suplane nt=101 ntr=2 | sufxdecon ntrf=-1`

## M2 (medium): sufxdecon with ntrf > 2·ntrw: C produces output, Rust panics in LU_decomposition

- Rust: programs/sufxdecon/src/main.rs:300 (`ipvt` = 4·ntrw), :469; crates/oxwave-f-axb/src/lib.rs:58 (`assert!(index.len() >= n …)`).
- C: sufxdecon.c:290, :453.
- **What differs.** `ipvt` has 4·ntrw entries, but `LU_decomposition(2*ntrf, …)` writes 2·ntrf
  entries. With a single spatial window, the negative-index path of L2 does not fire, so C writes
  past `ipvt` (ASan with axb.c instrumented: WRITE 0 bytes after a 32-byte region) and exits 0
  with output. Rust trips the library assertion (rc 101). The selfdoc says "ntrf … smaller than
  ntrw", so this input is documented as invalid, but it is still a crash.
- **Repro.**
  ```
  suplane nt=151 ntr=3 | suaddnoise seed=7 sn=3 | sufxdecon ntrw=2 ntrf=5   -> SU rc0; OX rc101 "assertion failed: index.len() >= n …"
  suplane nt=151 ntr=1 | suaddnoise seed=7 sn=3 | sufxdecon ntrw=1 ntrf=3   -> same
  ```
- **Suggested fix.** Size `ipvt`/`info` as `max(4*ntrw, 2*ntrf)`; the extra entries emulate C's
  scribble into the adjacent heap. Or rule the input invalid and use a clean `err!`.
- **Golden case.** `suplane nt=151 ntr=3 | suaddnoise seed=7 sn=3 | sufxdecon ntrw=2 ntrf=5`
  (only if accepted as a golden; C writes out of bounds).

## M3 (medium): sumigpspi checks tmpdir writability with permission bits instead of access(W_OK)

- Rust: programs/sumigpspi/src/main.rs:151-155 (`metadata().permissions().readonly()`).
- C: sumigpspi.c:138 (`access(tmpdir, WRITE_OK)`).
- **What differs.** `readonly()` is true only when *no* write bit is set. For a directory that is
  not writable by the user but has an owner write bit (root-owned 755: `/`, `/usr`, `/System`), Rust
  passes the check and fails later with a different message.
- **Repro** (both exit 1):
  ```
  suspike nt=7 ntr=2 | sumigpspi nz=6 dz=10 vfile=v12 tmpdir=/
  SU: sumigpspi: you can't write in / (or it doesn't exist)
  OX: sumigpspi: subcalls.c: efopen: fopen failed
  ```
- **Suggested fix.** Use `access(path, 2)` through the same `extern "C"` that sufxdecon already
  uses (programs/sufxdecon/src/main.rs:82, 212).
- **Golden case.** `[stderr] sunull nt=7 ntr=1 | sushw key=d2 a=10 | sumigpspi nz=3 dz=10 vfile=missing tmpdir=/`

## M4 (medium, needs a ruling): sust invents an error for frequency bins outside [0, nfft]

- Rust: programs/sust/src/main.rs:148-155.
- C: sust.c:268-274 (`k=n+i; if(k>=len) k-=len; G[i]=h[k]*s`).
- **What differs.** For `fmax` above about 2×Nyquist, or `fmin` below about −d1, C reads `h` out of
  range. ASan, `fmax=300`: READ 0 bytes after the 560-byte `h` region. The C still exits 0 with
  deterministic output (md5 stable under environment perturbation). Rust exits 1 with a message
  that does not exist in SU:
  `frequency-bin range [1, 85] outside supported range [0, 70]`.
  Unlike the pfacc exclusions, this rejection was never ruled on and is not documented in
  cases.txt. "fmax=300 with dt=4 ms" is a plausible typo-level input.
- **Repro.**
  ```
  suplane nt=64 ntr=1 | sust fmax=300        -> SU rc0 (output); OX rc1 "sust: frequency-bin range [1, 85] outside supported range [0, 70]"
  suplane nt=64 ntr=1 | sust fmin=-20 fmax=50 -> SU rc0; OX rc1 "... [-4, 15] ..."
  ```
- **Suggested fix.** Decide, consistently with the pfacc policy, and record the decision as an
  `# UB-removed:` comment with these lines. If it is kept, the message should at least be prefixed
  so users know it is an OxWave-only check.
- **Golden case.** None possible (C reads uninitialised heap); document instead.

## M5 (medium, process): the UB exclusions are not applied, so check_program fails these programs

- Files: goldens/succepstrum/manifest.json and goldens/sust/manifest.json were last regenerated
  at 16:40/16:45. The `# UB-removed:` edits to cases.txt came later (20:27).
- **What differs.** `python3 tools/check_program.py succepstrum` reports 43 pass / **4 FAIL**
  (cases 43–46: nt=46, 48, 121, 126). `sust` reports 41 pass / **3 FAIL** (cases 6, 14, 15:
  nt=96, 96, 95). These are exactly the removed lines. Anyone running check_program sees both
  programs fail.
- **Fix.** Rerun `tools/make_goldens.py succepstrum` and `tools/make_goldens.py sust`.

---

## Low

- **L1. Rust quietly produces output where the C segfaults.** In both cases below C is UB and
  crashes, while Rust prints traces with rc 0.
  - sust with lo=hi=0 (`suplane nt=64 ntr=1 | sust fmin=-3 fmax=-2`, or
    `sunull nt=8 ntr=1 | suaddnoise seed=29 | sust dt=0.015625 fmin=-4 fmax=-6`): C rc 139, Rust
    rc 0 with 1 trace. `st()` resets hi=len/2 and writes len/2+1 rows into a 1-row `vgr`.
  - sumigpspi nz=-1 with an existing vfile: C rc 139, Rust rc 0 with ns=0 traces. Rust clamps
    nz<0 to 0 at programs/sumigpspi/src/main.rs:210.

  Consider a clean `err!`, or just document it.
- **L2. sufxdecon with ntrw < ntrf ≤ 2·ntrw and several spatial windows crashes in both.**
  For example `suplane nt=251 ntr=20 | sufxdecon ntrw=3 ntrf=4`: C SIGSEGV (negative
  `fdata[jw + jx*ntrw - ntrf]` row), Rust index-out-of-bounds panic at main.rs:391.
  Documented-invalid input; both crash.
- **L3. Stale TODO comments.** programs/sumigpspi/src/main.rs:39-44 and
  goldens/sumigpspi/cases.txt:81-84 still say cwp_cexp keeps a double intermediate and that the
  "fix belongs in oxwave-cwp". It is already fixed (crates/oxwave-cwp/src/complex.rs:157-161
  rounds `exp` to `$T`), and the PSdemo cases pass.
- **L4. spsplot goldens 1 and 2 are dead.** They are marked "unstable" only because `trimodel …
  sfill=0.5,0.5,0,0,1,0.1,0.2` itself segfaults (rc 139, empty model). So the real
  trimodel→spsplot pipeline is covered only by a base64 fixture. I verified a working replacement
  that matches the original:
  `trimodel xmin=0 xmax=2 zmin=0 zmax=1 xedge=0,2 zedge=0,0 sedge=0,0 xedge=0,0.5,1,1.5,2 zedge=0.3,0.35,0.4,0.3,0.25 sedge=0,0,0,0,0 xedge=0,2 zedge=1,1 sedge=0,0 sfill=0.1,0.1,0,0,0.44,0,0 sfill=0.1,0.8,0,0,0.25,0.01,0.02 kedge=1,2,3 > m.bin && spsplot < m.bin | sed -e '/^%%CreationDate/d' -e '/^%%Creator/d'`
  (also with `style=normal gridx=dash … title=Model`).

  Separately, Rust spsplot on empty stdin also dies with SIGSEGV (rc 139), the same as C: a
  segfault inside the unsafe `readModel` port rather than an error.
- **L5. Minor fidelity.** Rust sumigpspi zero-initialises buffers that C leaves uninitialised
  (`v` on a short vfile read, PSPI `V[]` slots that are never assigned), and sufxdecon zero-fills
  out-of-band `ffreq` bins. These only matter on paths that are already UB in C. Also,
  sumigpspi's named temp files are unlinked immediately and no signal handler is installed.
  Neither is observable in normal runs.

## Judgement

The numerical cores of all five programs are faithful line-by-line translations, and I found no
precision regressions:
- C float stays f32 throughout;
- double promotions from unsuffixed literals and math calls are reproduced;
- the cwp_cexp rounding fix is in place;
- the phase-unwrap and LU library ports match the C.

Differential runs on new realistic synthetics agreed everywhere except on paths where the C itself
has undefined behaviour. The succepstrum/sust exclusions are justified: ASan places
every rejected length's overflow inside pfacc, and the sweeps show no Rust bug hiding behind them.
The real gaps are where the golden suite never went:
- sufxdecon panics on a plausible long-trace/short-window configuration (H1); the fix is cheap and
  validated (zero slack matches the C);
- sufxdecon has two more panic-instead-of-error paths (M1, M2);
- sumigpspi's tmpdir check is not `access()` (M3);
- sust adds an error message that has not been ruled on (M4);
- the stale manifests make two programs fail check_program today (M5).

I found no tests weakened or special-cased to make the port pass.
