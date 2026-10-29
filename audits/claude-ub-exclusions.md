# UB-dependent goldens: exclusion lists

Auditor: Claude. Read-only on the repo apart from this file and `audits/claude-ub-exclusions.json`.
Experiments were run in `/private/tmp/claude-501/ubx`, which is ephemeral; the recipe below rebuilds it.

Follow-up to `audits/claude-wave3-compression.md` (H1, H2, M3) and `audits/claude-coverage.md`
(Divergence 3). It covers 6 programs and lists 24 cases. The verbatim case strings and the EXCLUDE.txt
header text are in `audits/claude-ub-exclusions.json`; they are generated from
`oxharness.collect_cases` and checked against each `manifest.json`.

## Summary

| program | exclude (both criteria) | of which the minimal fix changes the output | of which only scribble changes it |
|---|---|---|---|
| wpc1uncomp2 | 6: k 3,4,5,6,7,8 | 3,4,8 (compressor realloc) | 5,6,7 (4 is in both columns) |
| wpcuncompress | 3: k 4,10,11 | 4,11 (compressor realloc); 10 (decoder zeroflag) | none |
| sufctanismod | 8: k 16,17,22,25,28,38,39,40 | 25 only | 16,17,22,28,38,39,40 |
| dctcomp | 2: k 20,21 | none | 20,21 |
| wptcomp | 4: k 35,36,55,56 | none | all 4 |
| wtcomp | 1: k 10 | none | 10 |

`k` is the index in the current manifest. The indices shift once an EXCLUDE.txt exists, so use the strings,
not the numbers.

**About the two criteria.** The exclusion rule has two parts. A case is excluded if its output changes
under a minimal fix of the C. For uninitialised memory, it is also excluded if its output changes when
fresh malloc memory is scribbled. The JSON `exclude` lists apply both parts. With the fix criterion
alone, drop the "scribble only" cases: keep wpc1uncomp2 5,6,7; sufctanismod 16,17,22,28,38,39,40;
all of dctcomp, wptcomp and wtcomp. Every scribble-only case has these properties:
- the zero-initialising fix reproduces the golden byte for byte;
- the Rust port zero-initialises and also matches the golden.

So the zero-heap result is what these goldens record, but no C code guarantees it.

**Consequence for coverage.** Under both criteria, wpc1uncomp2 keeps only its error-path cases (1, 2, 9, 10;
8 is already unstable). Even 1 and 2 never reach the decoder: their `sustrip` has no `outpar=`,
so it fails in the harness (`efopen: fopen failed`). The same is true of wpcuncompress 1 and 2.
wpcuncompress keeps 12 round-trip cases: 3, 5 to 9, and 12 to 17. If you apply both, consider adding wpc1uncomp2
round-trip cases whose average block is not ~0 (`maxave > FLT_MIN`). In that case `wpc1Decoder` sets `ave`.

## Method

- **Sources:** copied from `vendor/su-build/src` into `ubx/src-<variant>`. Libraries `libwpc1`, `libwpc` and
  `libcomp` were rebuilt as static archives from the copied lib sources and linked first. The vendor
  `libpar`, `libcwp` and `libsu` were used as they are.
- **Flags:** exactly the ones in `vendor/su-build/build.log`: `-g -std=c90 -m64 -Wall -ansi -Wno-long-long
  -ffp-contract=off -DCWP_LITTLE_ENDIAN`.
  - The comp programs and libraries have no `-DSUXDR`, because their Makefiles don't pass it. They don't
    use SU trace I/O.
  - sufctanismod gets `-DSUXDR` and links `-lsu -lpar -lcwp -lm`, as in build.log:1600.
- **Variants:**
  - `orig`: unmodified.
  - `fixcomp`: the one-line realloc fixes, `code = wpc1->code = realloc(...)` at wpc1.c:230 and
    `code = wpc->code = realloc(...)` at wpc.c:209/261/280.
  - `fixdec`, the decoder and sufctanismod fixes:
    - `wpccoding.c` `codeDecoder`: `quant->zeroflag = zeroflag; quant->ave = ave;` on the all-zero path,
      and `quant->zeroflag = zeroflag;` on the normal path.
    - `wpc1coding.c` `wpc1Decoder`: `else for(i...) quant[i]->ave = 0.;` when `maxave <= FLT_MIN`, in
      both the DC and the normal branch.
    - sufctanismod `anis_solver2`: `memset` of `u1`, `f0` and `f1` right after `alloc2float`.
  - `fixz`: dctcomp, wptcomp and wtcomp with `memset(f[0],0,..)` and `memset(g[0],0,..)` after allocation
    (calloc semantics).
  - `zall`: sufctanismod with every `alloc1float`/`alloc2float` zero-filled.
  - `asan-orig`: `orig` built with `-O0 -fsanitize=address`.
- **Runner:** `ubx/runv.py <prog> <config>` reruns every manifest case through `oxharness.run_case`, the
  harness's own environment. A shim directory goes in front of `vendor/su-build/bin` and replaces only
  the named stages, e.g. `wpccompress=fixcomp+wpcuncompress=rust`. The runner records:
  - rc, length, md5 and NaN count of stdout;
  - whether stdout is byte-equal to the golden;
  - whether `oxharness.compare_stdout` still passes ("passes cmp").
- **Scribble ("pre"):** a `DYLD_INSERT_LIBRARIES` interposer (`ubx/prescribble.c`). It fills every fresh
  `malloc` block, and the grown tail of `realloc`, with 0xAA, and leaves `free` alone. So it detects
  reads of uninitialised memory only, not freed memory.
  - **Caution:** on this macOS (Darwin 27), `MallocPreScribble=1` on its own does nothing. A test
    program sees 0x00 in fresh blocks. Only `MallocScribble=1` turns on 0xAA pre-fill, and it also fills
    freed blocks with 0x55.
  - This may be why the earlier audit's `MallocPreScribble=1` sufctanismod run looked like "matches
    Rust". Use `MallocScribble=1` or the interposer.
  - The two methods gave identical results on every case here, except wpcuncompress 4 and
    wpc1uncomp2 8. Both decode freed-memory streams from the upstream compressor.
- **Baseline:** rerunning the vendor binaries reproduces every golden byte for byte, with two exceptions.
  Both come from the vendor compressor's freed-memory output, which depends on the environment:
  - wpc1uncomp2 8 is already marked unstable;
  - wpcuncompress 4: today's vendor run decodes to 0 NaN, the golden holds 16384 NaN.

  A rerun with exactly the harness PATH, with no shim, gives the same result. The compressor stream
  itself changes with the length of PATH (md5 d3de7ffe vs 7eee0876 for the same command).

In the tables, a cell reads `md5[:8]/bytes [nan=N]`. `=G` means byte-equal to the golden. "passes cmp" means
the bytes differ but the harness comparison still passes.

---

## 1+4. wpc1uncomp2 (compressor realloc fix, and the M3 decoder `ave`)

| k | golden | vendor rerun | fixcomp upstream | fixdec decoder | fixcomp + fixdec (all fixed) | scribble decoder | Rust (harness) | Rust + fixcomp upstream |
|---|---|---|---|---|---|---|---|---|
| 3 | d41d8cd9/0 | =G | 4ac07bc4/1536 | =G | 4ac07bc4/1536 | =G | =G | 4ac07bc4/1536 |
| 4 | a75d7d42/400 | =G | 15ffd782/1200 | =G | 15ffd782/1200 | db5abefd/400 | =G | 15ffd782/1200 |
| 5 | f1d3ff84/4 | =G | =G | =G | =G | c87ad371/4 | =G | =G |
| 6 | 2cb67822/1028 | =G | =G | =G | =G | d8109d7e/1028 | =G | =G |
| 7 | 522fd2be/3084 | =G | =G | =G | =G | 54750a13/3084 | =G | =G |
| 8 (unstable) | 538e7ffd/8192 | d41d8cd9/0 | 1ec6b9e9/8192 | e7007624/8192 | 1ec6b9e9/8192 | 67258aeb/8192 | c4d4c95e/8192 | 1ec6b9e9/8192 |

**Evidence**
- **ASan on the upstream `wpc1comp2` (asan-orig):** heap-use-after-free in memcpy on 3, 4 and 8 (rc 134).
  It is clean on 5, 6 and 7.
- **ASan on the decoder:** heap-buffer-overflow in `bookDecoder` wpc1coding.c:598 on 3, which is the
  garbage stream.
- **3:** the golden is 0 bytes because decoding the garbage stream fails. With the fixed upstream the
  real 1536 bytes come out.
- **4:** the golden holds only trace 0 (400 bytes); fixed, it is 1200 bytes.
- **4 to 7, scribble:** sample 0 or 256 becomes −3.03e-13. That is the 0xAAAAAAAA float, added as `ave`,
  and the harness comparison fails.
- **Zero-ave fix (fixdec):** gives exactly the golden on 4 to 7.
- **8:** the vendor rerun no longer reproduces the golden. Every variant gives a different answer.

**Rust:** matches the fixed C byte for byte on every excluded case. Run with the fixed upstream, it gives
the same hashes as fixcomp+fixdec, including 8.

## 1+2. wpcuncompress (compressor realloc fix, and the H1 decoder zeroflag)

| k | golden | vendor rerun | fixcomp upstream | fixdec decoder | fixcomp + fixdec (all fixed) | scribble decoder | Rust (harness) | Rust + fixcomp upstream |
|---|---|---|---|---|---|---|---|---|
| 4 | 3591e7e3/65536 nan=16384 | fcd6bcb5 nan=0 | 1c036714 nan=0 | fcd6bcb5 nan=0 | 1c036714 nan=0 | 9802287a nan=0 | =G | 1c036714 nan=0 |
| 10 | 8ba905c9/266252 nan=1027 | =G | =G | 0bb626dd nan=0 | 0bb626dd nan=0 | 1e6642df nan=0 | =G | =G (nan=1027) |
| 11 | f0aac468/660 nan=0 | =G | f413294c nan=165 | =G | =G | 2988ba82 nan=0 | =G | f413294c nan=165 |

**Per case**
- **4:** an H2 case.
  - With the fixed compressor the 16384 all-NaN samples become real data (max ≈ 3).
  - The decoder fix alone doesn't change it.
  - The golden is not reproducible today, even by the vendor binaries; see Method.
- **10:** an H1 case, and the only case the decoder fix changes. The 1027 NaN in the zero edge tiles
  become 0. The fixed compressor doesn't change it.
  - **Correction to the audit:** it says the H1 fix changes "4 and 10". Only 10 changes. Case 4 is H2
    and does not depend on the decoder fix.
- **11:** the subtle one.
  - With the vendor compressor the stream is garbage. The decode fails, and C writes out its output
    buffer without filling it. The golden zeros are fresh-heap zeros: scribble gives −3.03e-13.
  - With the fixed compressor alone it decodes the real stream, but H1 turns the all-zero tile into
    165 NaN.
  - With both fixes it gives zeros again, equal to the golden.
  - It qualifies under the item-1 test (the compressor fix changes it) and under the scribble test.
    But the all-fixed output equals the golden. If you would rather keep it, it is the one case where
    that is arguable.

**Evidence**
- The wpccompress ASan build is not usable here. It exits 139 on every case, including cases that
  should exit 1, with no ASan report, so it is not cited.
- The existing `goldens/wpccompress/EXCLUDE.txt` header already documents the use-after-free.
- Decoder ASan finds nothing on 4, 10 and 11. Its 0xBE fill changes all three; uninitialised reads are
  not ASan errors.

**Rust does NOT match the fixed C on 10 and 11.** The Rust decoder deliberately reproduces H1
(`crates/oxwave-lib-wpc2/src/codec.rs:125-140` keeps the zero-heap quant state):
- **10:** Rust gives the golden's NaN, with either upstream. The fixed C gives 0.
- **11:** with the fixed upstream, Rust gives 165 NaN. That equals fixcomp + orig decoder, not the fully
  fixed C, which gives zeros.

On 4, Rust with the fixed upstream equals the fixed C (1c036714). Rust passes all three today only because
it is fed the vendor garbage stream and repeats the zero-heap behaviour. After exclusion nothing tests
the H1 path. If the fixed behaviour is wanted, the Rust `codeDecoder` needs the same two
assignments.

## 3. sufctanismod

| k | golden | vendor rerun | orig rebuilt | fixdec (u1/f0/f1 memset) | zall (all allocs zeroed) | scribble orig | scribble fixdec | Rust |
|---|---|---|---|---|---|---|---|---|
| 16 | 7bf118e0/1020 | =G | =G | =G | =G | e30bb6fc (passes cmp) | =G | =G |
| 17 | 80dd9299/1020 | =G | =G | =G | =G | 0a6eaf91 (passes cmp) | =G | =G |
| 22 | 748dc6e3/572 | =G | =G | =G | =G | eeb74851 (passes cmp) | =G | =G |
| 25 | 195b652d/20604 | =G | =G | **d15234e1** (passes cmp) | d15234e1 | eb4f2cd3 (passes cmp) | d15234e1 | **d15234e1** (passes cmp) |
| 28 | d6c8a8db/396 | =G | =G | =G | =G | 08a851c1 (passes cmp) | =G | =G |
| 38 | 7acc682c/616 | =G | =G | =G | =G | e513beca (passes cmp) | e513beca | =G |
| 39 | 2bce3a6f/396 nan=99 | =G | =G | =G | =G | 0a15f509 nan=0 (**FAILS cmp**) | 0a15f509 | dbc7c104 nan=99 (passes cmp) |
| 40 | dbc7c104/396 nan=99 | =G | =G | =G | =G | be0c6ea7 nan=99 (passes cmp) | be0c6ea7 | =G |

Other cases:
- 2 (`dofct=1`) writes to a file and its stdout is empty.
- 18 (`dofct=1 eta0=0 eta=0`) is unaffected by both the fix and the scribble.
- No other case of the 68 changes under scribble.

**Mechanisms**
- **16, 17, 25, 28 (dofct=1), as in the audit, plus 22:**
  - `anis_solver2` allocates `u1`, `f0` and `f1` on every call and never zeroes them.
  - `fct2d1o` reads `u1` boundary cells that nothing wrote. On step 1 these are fresh zeros. Later they
    are recycled blocks that still hold the previous step's values. That is why only 25 changes under
    the memset fix.
  - **22** (`isurf=0`, dofct=0) is new. With `isurf=0`, `absorb1` never writes the top row
    `u1[i][mbz1]`, and the update loop copies that row into `u`.
  - With the memset, the scribble no longer changes any of the five. The scribble-fixdec column shows this.
- **38:** `vsp_x` is allocated with `alloc2float(nt,nz)`, 220 bytes, but never filled when `vspnx=-1`.
  It is written to `vx` anyway. A single-allocation bisect with an interposer that scribbles only the
  220-byte block changes the output.
- **39, 40:** the case's `sustrip > grid` has no `outpar=`, so the harness's sustrip fails
  (`efopen: fopen failed`; see 39.err) and `grid` is empty.
  - `read_parameter`'s `fread` leaves `rho` uninitialised. Bisect: only the 15th 396-byte `alloc2float`,
    `rho`, matters.
  - The golden's NaN come from 1/0 on zero-heap `rho`.
  - Alternative to excluding them: add `outpar=/dev/null` to the case and
    regenerate the goldens.

**Tolerance:** every changed output except scribbled 39 still passes the harness comparison. The bytes differ but stay
within rtol 1e-5. So under the fix criterion, 25 is a "changes but would still pass" exclusion. Rust
already gives the fixed C's exact bytes on 25.

**Rust:** matches the fixed C byte for byte on 16, 17, 22, 25, 28, 38 and 40. On 39 it matches neither
golden nor zall byte for byte. All 99 values are NaN in both, and the only difference is the NaN sign
bit: Rust writes 0x7fc00000, C writes 0xffc00000. It passes the comparison. Unrelated: the current
`target/release/sufctanismod` fails cases 66, 67 and 68 (`mkdir snapshot?.data; …`): Rust exits 1, C exits
0 with 252 bytes.

## 4. M3 compressor cases (dctcomp, wptcomp, wtcomp)

| prog k | golden | vendor rerun / orig rebuilt | zero-fix (fixz) | scribble | Rust |
|---|---|---|---|---|---|
| dctcomp 20 `dctcomp n1=1 n2=1` | e40efe79/40 | =G / =G | =G | ee711824/284 | =G |
| dctcomp 21 | d2e3d3e7/40 | =G / =G | =G | 56612ffb/286 | =G |
| wptcomp 35 `wptcomp n1=8 n2=4` | 43cf67f1/44 | =G / =G | =G | 3bf33e2a/48 | =G |
| wptcomp 36 | 43cf67f1/44 | =G / =G | =G | 6f8a7feb/68 | =G |
| wptcomp 55 `nstage1=3` | b95fce53/44 | =G / =G | =G | 8daa865f/44 | =G |
| wptcomp 56 `nstage2=3` | 7bc46b07/44 | =G / =G | =G | 4ba67e24/44 | =G |
| wtcomp 10 `</dev/null` | b82f0ffc/33 | =G / =G | =G | 6d028fcf/33 | =G |

The M3 claims hold exactly. No other case in these three programs changes under scribble. ASan's 0xBE
fill also changes all seven; uninitialised reads are not ASan errors. Rust matches the zero-filled C on
all seven. The M3 wpc1uncomp2 claim (4 to 7) also holds; see the wpc1uncomp2 section above.

Unrelated: the rebuilt `orig` dctcomp gives rc −9 instead of −11 on dctcomp 38 (`n1=1 n2=-1`), a known UB
case. The vendor binary matches the golden.

## Rerun recipe

```sh
U=/private/tmp/claude-501/ubx                       # ./build.sh <variant>  (orig fixcomp fixdec fixz), build-asan.sh orig
python3 $U/runv.py wpcuncompress wpccompress=fixcomp+wpcuncompress=fixdec   # any '+'-joined stage=variant list;
# variants: <dir suffix of bin-*>, rust, pre (vendor + 0xAA interposer), pre:<variant>, scribble (MallocScribble=1)
```
