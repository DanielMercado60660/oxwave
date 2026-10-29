# Audit: seismic compression family

Scope: `crates/oxwave-lib-comp-dct`, `crates/oxwave-lib-wpc1`, `crates/oxwave-lib-wpc2`; programs
`dctcomp dctuncomp wpc1comp2 wpc1uncomp2 wpccompress wpcuncompress wptcomp wptuncomp wtcomp wtuncomp`.
Auditor: Claude (read-only on the repo; experiments in `/private/tmp/claude-501/oxaudit`).

All 10 programs pass `tools/check_program.py` (271 cases; wpc1uncomp2 has 2 skipped as unstable).
License headers match `tools/cwp_header.txt` on every Rust file in scope. I found no test-specific
code: no code reads the environment, no golden-specific constants, and the EXCLUDE handling sits in the
harness, not in the ports.

**Summary: 3 high, 4 medium, 2 low.**

## How the evidence was produced (so it can be rebuilt)

Scratch dir: `/private/tmp/claude-501/oxaudit` (ephemeral). Sources were copied from
`vendor/su-original/src/comp/{dct,dwpt}` and built four ways with SU's own flags
(`-g -std=c90 -m64 -ansi -Wno-long-long -ffp-contract=off -DCWP_LITTLE_ENDIAN`, no `-O`):

```sh
SB=vendor/su-build
CF="-g -std=c90 -m64 -w -ansi -Wno-long-long -ffp-contract=off -DCWP_LITTLE_ENDIAN -I$SB/include"
L="-L$SB/lib -lpar -lcwp -lm"
gcc $CF -I$S/dwpt/1d/include $S/dwpt/1d/main/wpc1comp2.c $S/dwpt/1d/lib/*.c $L -o wpc1comp2   # same pattern for the other 9
# variants: orig (as-is), fix (one-line realloc fix), asan (+ -O0 -fsanitize=address), asanfix
# the one-line fix, exactly as the UB rule requires:
sed -i '' 's/code = (unsigned char \*)realloc(code, size\*sizeof(char));/code = wpc1->code = (unsigned char *)realloc(code, size*sizeof(char));/' dwpt/1d/lib/wpc1.c
sed -i '' 's/code = (unsigned char \*)realloc(code, mbound\*sizeof(char));/code = wpc->code = (unsigned char *)realloc(code, mbound*sizeof(char));/; s/code = (unsigned char \*)realloc(code, size\*sizeof(char));/code = wpc->code = (unsigned char *)realloc(code, size*sizeof(char));/' dwpt/2d/lib/wpc.c
```

A small runner (`runv.py`) reran every manifest case plus the EXCLUDE.txt lines through the harness's
own `oxharness.run_case` environment. A shim directory let me replace any subset of pipeline stages
(target and/or upstream compressor) with a chosen build. It recorded stdout, exit code and ASan
reports. "Scribble" runs wrap only the target vendor binary with `MallocPreScribble=1
MallocScribble=1` (macOS libmalloc fills new blocks with 0xAA and freed blocks with 0x55). Any
output change in those runs means the program read uninitialised or freed heap memory.

Sanity check: the rebuilt `orig` binaries reproduce every non-excluded golden byte for byte. The
only exceptions are dctcomp 38 (`n1=1 n2=-1`, a UB case) and two decoder cases (wpc1uncomp2 8,
wpcuncompress 4) discussed in H2. Those two differ between two builds of identical source.

---

## Verdict on the EXCLUDE.txt decisions

**Compressors (wpc1comp2: 6, wpccompress: 8): yes, correct and minimal.**
- All 14 excluded cases change under the one-line fix. Their hashes differ between the vendor
  binary, the rebuilt original and the fixed build. Example: wpc1comp2 X0 gives `48753941` (vendor)
  and `dbe635c4` (fixed); wpccompress X3 gives `9105d51f` (vendor) and `6d9b4cc9` (fixed).
- Every non-excluded compressor golden is byte-identical with the fix. That includes cases where
  ASan does report the heap-use-after-free but the freed block still held the right bytes:
  wpc1comp2 14/32/34 and wpccompress 8–27. Keeping those as tests is consistent with the rule.
- **(c) Rust on the excluded cases:** the Rust compressors are byte-identical to the fixed C build
  on all 14 excluded cases and on every other case. This is the intended, non-UB behaviour.
- I checked the other realloc sites for the same pattern. `buffRealloc1` in `dct/lib/buffalloc.c:82`,
  `dwpt/1d/lib/wbuffalloc.c:80` and `dwpt/2d/lib/wpcbuffal.c:88` all assign through
  `buff->code = realloc(...)`, and no caller caches the old pointer. The stale-pointer pattern exists
  only at `wpc1.c:230` and `wpc.c:209/261/280`. After the fix, the ASan run (`bin-asanfix`) reports
  only the deliberate error-path cases: negative or zero dimensions and empty input.

**Decoders: no, incomplete.** The same rule, applied one stage downstream, is not met. Five decoder
goldens consume the vendor compressor's freed-memory output, and they change when only the
compressor line is fixed (finding H2).

---

## High

### H1 — wpcuncompress turns every all-zero 256×256 tile into NaN; golden 10 enshrines 1027 NaNs

- **C:** `vendor/su-original/src/comp/dwpt/2d/lib/wpccoding.c:373-383`. `codeDecoder` reads
  `zeroflag` and, for an all-zero tile, the tile average `ave` into locals. It returns 1 without
  storing either into `quant`.
- **C:** `dwpt/2d/lib/wpcquant.c:223-240`. `quantDe` then reads the uninitialised
  `quant->zeroflag/step/ave/qx` from a fresh `malloc`. With the zero-filled pages this allocator
  hands out, `zeroflag==0` sends it down the non-zero path, giving `qx(0) * (1/step(0) = inf) = NaN`.
- **Rust:** `crates/oxwave-lib-wpc2/src/codec.rs:125-140` (`wpcDecoder` builds `WpcQuant::new`,
  all zero, and never sets `zeroflag`/`ave`); `quant.rs:142-160`. The comment at codec.rs:125 says
  this is intentional ("initializes quantization state to zero"). The result is NaN, not zero.
- **Golden that depends on UB:** goldens/wpcuncompress case 10,
  `suplane nt=257 ntr=259 | … | wpccompress n1=257 n2=259 | wpcuncompress`. The golden output
  contains **1027 NaN samples**: sample 256 of traces 0–255, then all of traces 256–258. The input
  is exactly 0 in that region. Rust passes because the comparator treats NaN as equal to NaN.
  Scribbling only the vendor decoder changes cases 4, 10 and 11.
- **Plausible-input repro:** a 512×256 section whose first 256 samples are zero (a top mute), made
  from a susynlv gather:
  ```
  wpccompress n1=512 n2=256 error=0.1 < mute.bin | wpcuncompress
  vendor C decoder : 131072 samples, NaN = 65536 (the whole muted tile)
  Rust decoder     : 131072 samples, NaN = 65536
  C with fix below : NaN = 0, relRMS = 0.115
  ```
  Any section with n1 or n2 not a multiple of 256 and a zero edge tile, any muted region of at
  least 256×256, or any dead-trace block is affected.
- **Suggested fix:** store the decoded flag and average (`q.zeroflag = zeroflag; q.ave = ave`) in
  Rust `codeDecoder`, and fill the last sub-block with `ave` on the zero path. That is what the
  encoder wrote (`wpcquant.c:126,146`). I patched the C this way (`bin-fix3`). Only wpcuncompress
  goldens 4 and 10 change, and both are UB-dependent anyway. Regenerate them from a fixed decoder,
  or exclude them.

### H2 — Five decoder goldens are built on the compressor's freed-memory output and were not excluded

The check harness substitutes only the program under test, so the upstream `wpc1comp2` /
`wpccompress` stage is the vendor binary with the realloc bug. I ran the vendor decoders with
`wpc1comp2=bin-fix` / `wpccompress=bin-fix` upstream:

| golden | what the golden contains | with fixed upstream |
|---|---|---|
| wpc1uncomp2 3 `suplane nt=128 ntr=3 … wpc1comp2 n1=128 \| wpc1uncomp2` | **0 bytes** (decode fails on the garbage stream; ASan: heap-buffer-overflow in `bookDecoder` wpc1coding.c:598) | 1536 bytes, the real data |
| wpc1uncomp2 4 (`n1=100 error=0.05`, 3 traces) | 400 bytes (only trace 0) | 1200 bytes |
| wpc1uncomp2 8 (`suplane nt=1024 ntr=2`) | decoded garbage: relRMS vs input 0.34 | relRMS 0.0019 |
| wpcuncompress 4 (`suplane 128×128 error=0.05`) | **16384 / 16384 samples NaN** | real data, max 2.997, no NaN |
| wpcuncompress 11 (`sunull 33×5`) | zeros, but only because the garbage stream failed the decode and C wrote its fresh zero buffer | 165 NaN (this is H1 again) |

- wpc1uncomp2 8 is also environment-dependent. The vendor `wpc1comp2 n1=1024` output in my shell
  differs from the one recorded in the harness, so the freed memory holds environment-dependent
  bytes. Decoding it gives a third answer (relRMS 0.998). The golden is stable only because the
  harness environment is fixed.
- Separately, a rebuilt-from-identical-source decoder already differs from the vendor decoder on
  cases 8 and 4. The decoders themselves read uninitialised memory (see M4 and H1).
- **(c) Rust behaviour is right:** with fixed upstream streams, the Rust decoders are byte-identical
  to the fixed C decoders on all 32 wpc1uncomp2 and wpcuncompress cases.
- **Suggested fix:** regenerate these five decoder goldens from a fixed compressor, or check in
  `.wpc` fixtures produced by the Rust or fixed compressor and `cat` them. Otherwise add them to the
  decoders' EXCLUDE.txt. As they stand they test "decode freed memory like libmalloc did in that
  run".

### H3 — wpcCompress's "almost never happens" retry re-encodes an already-transformed block (inherited; Rust is faithful)

- **C:** `dwpt/2d/lib/wpc.c:203-216` and `:255-268`. When `wpcEncoder` returns `WPC_EOB`, the code
  grows the buffer and re-encodes the same `block`. But `wavePack2` writes its result back into the
  block in place (`wpcpack2.c:106`, `wavePack_row(z,x,…)`). The retry therefore transforms twice.
- **Rust:** `crates/oxwave-lib-wpc2/src/codec.rs:144-188`. The doc comment at line 144 states the
  behaviour is retained.
- **Trigger:** the encoded size exceeds `mbound = n1*n2` bytes, i.e. compression ratio below 4.
  This happens with the default `error=0.01` on ordinary data.
- **Repro** (512×128 susynlv gather + sufilter + suaddnoise sn=30; relRMS = rms(out−in)/rms(in)):

  | | error=0.01 | error=0.1 |
  |---|---|---|
  | wpccompress/wpcuncompress, vendor C | 1.000 | 0.0807 |
  | fixed C (realloc only) | **1.029** | 0.0807 |
  | Rust | **1.029** | 0.0807 |
  | C with retry restoring the saved block (`bin-fix2`) | 0.0081 | — (0.0405 at error=0.05) |

  The Rust compressor is byte-identical to the fixed C. No golden reaches the EOB path with real
  data, because the full 256×256 cases compress well.
- **Suggested fix:** this is a project decision, not a port defect. Either keep the fidelity and
  document that wpccompress is unusable below ratio 4, or copy the block before `wpcEncoder` and
  restore it on `-1`. The second changes only the excluded cases. Either way, add a pinning case,
  e.g. `sunull nt=256 ntr=256 | suaddnoise seed=1 | … | wpccompress n1=256 n2=256 error=0.001`
  through the round trip.

---

## Medium

### M1 — wtuncomp/wptuncomp panic (exit 101) on streams their own compressors produce; C gets SIGSEGV (139)

- **Rust:** `programs/wtuncomp/src/main.rs:77` and `programs/wptuncomp/src/main.rs:83` call
  `buffAlloc1(nsize as usize)`. With `nsize = -1` this panics with `capacity overflow`.
  wtuncomp main.rs:69-70 also panics in `(n1 - 1).ilog2()` when n1 ≤ 0, where C loops forever.
  The `waveGetfilter(..).unwrap()` calls at wtuncomp:79 and wptuncomp:80-81 panic on a bad nfilter.
- **C:** `dct/main/wtuncomp.c` and `wptuncomp.c`: `buffAlloc1(nsize)`, then `fread` into NULL, then
  SIGSEGV.
- The sibling port `programs/dctuncomp/src/main.rs:106-115` handles the same `nsize<0` stream with
  `raise(11)`. The two newer ports are inconsistent with it.
- **Root cause (inherited):** dctcomp/wptcomp/wtcomp write `nsize=-1` with exit 0 when `pEncode`'s
  bounded buffer cannot hold the Huffman output. That happens on small sections at any error, and on
  noise at error ≤ 0.001. On the realistic 512×128 gather it happens at error=1e-4. C and Rust agree
  on this, so the file is silently undecodable in both.
- **Repro:** a differential sweep over n1 ∈ {1…513}, n2 ∈ {1…257}, error ∈ {0.001, 0.05, 0.3},
  plane and noise inputs: 3180 compress+decompress pairs. All compressed streams and decoded outputs
  are byte-identical between Rust and C (fixed C for wpc), except **188 decoder runs where C exits
  139 and Rust exits 101 with a panic message**. Example:
  ```
  sunull nt=128 ntr=1 | suaddnoise seed=129 | sugain scale=1000 | sustrip outpar=/dev/null \
    | wtcomp n1=128 n2=1 error=0.001 | wtuncomp     # C: rc 139   Rust: rc 101 "capacity overflow"
  ```
- **Suggested fix:** mirror dctuncomp. Route a negative `nsize`, n1/n2 ≤ 0 and a failed
  `waveGetfilter` through the same `signal(11,0); raise(11)` path. Add one `[stderr]` case per
  decoder fed by a `nsize=-1` stream.

### M2 — Negative `nstage1/nstage2` in wptcomp/wtcomp: Rust silently compresses zeros

- **Rust:** `programs/wptcomp/src/main.rs:181-182` and `programs/wtcomp/src/main.rs:87` pass
  `nstage as usize`, so −1 becomes huge. `stage > pow` then makes `wavePack_2`/`waveTrans_2` return
  without transforming, and `g` stays zero.
- **C:** `wavePack_1` takes a partial path with `1<<stage` for a negative shift (undefined).
- **Repro** (`suplane nt=16 ntr=8`, n1=16 n2=8, md5 of stdout):
  - wptcomp `nstage1=-1`: C `7fe25adb…` exit 0 (deterministic; the rebuilt original agrees);
    Rust `9f72bc57…` exit 0, an all-zero payload.
  - wtcomp `nstage1=-1`: C SIGSEGV (139, no output); Rust exit 0, writes a 5-byte "compressed"
    zero section ("compression ratio = 102.4").
- **Why medium:** valid-looking parameters cause silent data loss in Rust. The C behaviour is UB,
  so do not try to reproduce it.
- **Suggested fix:** treat `nstage < 0` like C's crash (raise SIGSEGV), or at minimum not as "no
  transform". Exclude or document; no golden covers it.

### M3 — The decoders read uninitialised averages (wpc1) and the C compressors quantise uninitialised samples; goldens pass by allocator coincidence

Scribbling only the vendor target changes these passing goldens:
- dctcomp 20, 21 (`dctcomp n1=1 n2=1` on empty stdin; `printf '\000\000' | dctcomp n1=3 n2=2`):
  the `fread` short read leaves `alloc2float` memory uninitialised.
- wptcomp 35, 36 (short reads), and wptcomp 55, 56 (`nstage1/2=3 > npow`). In 55/56,
  `wavePack_2` returns 0 without writing `g`, and C then quantises the uninitialised `g`
  (`wptcomp.c:116-120`).
- wtcomp 10 (`wtcomp n1=1 n2=1 </dev/null`).
- wpc1uncomp2 4, 5, 6, 7. `wpc1Decoder` (`wpc1coding.c:250-300`) never sets `quant[i]->ave` when
  `maxave <= FLT_MIN`, and `wpc1Dequant` adds it.

The Rust ports zero-initialise in every one of these cases (e.g. `crates/oxwave-lib-wpc1/src/lib.rs:481`).
That gives the deterministic, intended result, which equals the golden only because this macOS
heap hands out zeroed fresh blocks. No change to the ports is needed. The goldens should be marked
fragile: a different libc, allocator or link order breaks them, as the rebuilt binaries already show
for wpc1uncomp2 8 and wpcuncompress 4.

### M4 — wpccompress goldens 8–27 read freed memory but are kept (consistent with the rule, but fragile)

ASan reports a heap-use-after-free in `wpcWrite` (wpc.c:506) on every non-excluded full-tile
wpccompress case (8–27), and target-only scribbling changes all of them. The fixed build gives the
same bytes, so under the UB rule they correctly stay. The Rust output equals the fixed C
output, which is what matters. I note it so nobody later treats a scribble or ASan failure on these
cases as a Rust regression.

---

## Low

- **L1** Stray build debris inside crate dirs: `programs/wpccompress/.tmp/actual.wpc` and
  `*/.tmp/xcrun_db` in every crate in scope (`programs/dctcomp/.tmp/xcrun_db` is non-empty). Delete
  them, or add `.tmp/` to `.gitignore`.
- **L2** dctcomp 38 (`dctcomp n1=1 n2=-1`) asserts the vendor's SIGSEGV (golden rc −11) even though
  the comment in `goldens/dctcomp/cases.txt` says "do not assert a fixed C signal". The
  rebuilt-original C exits 0, and the ASan build times out. Rust passes only through its `raise(11)`
  emulation. Exclude the case or drop it.

---

## (d) Round-trip sanity

Input: `susynlv nt=512 nxm=1 nxo=128 … | sufilter f=5,10,40,60 | suaddnoise sn=30 seed=7`
(512×128, rms 0.26). Values are relRMS of the round trip against the requested `error=`, with the
compression ratio for the Rust tools.

| pair | error=0.01 | error=0.1 | error=0.001 |
|---|---|---|---|
| wpc1comp2/wpc1uncomp2, vendor C | **no output** (decode fails on the UAF stream) | **no output** | — |
| wpc1, fixed C / Rust (identical) | 0.0111 (ratio 4.3) | 0.111 (7.6) | — |
| wpccompress/wpcuncompress, vendor C | 1.000 | 0.0807 | — |
| wpc2, fixed C / Rust (identical) | **1.029** (H3) | 0.0807 (8.4) | — |
| dctcomp/dctuncomp C / Rust (identical) | 0.00997 (5.4) | 0.1002 (11.9) | 0.00100 (3.1) |
| wptcomp/wptuncomp C / Rust (identical) | 0.00999 | 0.0999 | 0.00100 |
| wtcomp/wtuncomp C / Rust (identical) | 0.00999 | 0.1000 | 0.00100 |

At error=1e-4, dct/wpt/wt fail in both C and Rust (`nsize=-1`, M1). The Rust tools always match the
fixed C round trip. The vendor binaries are worse than the Rust port for wpc1, and both are broken
for wpc2 at the default error.

## Precision

The 3180-case sweep compared compressed streams and decoded floats byte-exactly between Rust and C.
It found no f32/f64 divergence. Spot checks of `uniQuant`, `pEncode`/`pDecode`, `wpc1Quant` and the
compression-ratio `printf` confirm C float promotion is kept (literal `3.4641` as double, `1./n`
narrowed to float, `(int)(atmp+.5)` in double).

## Judgement

The Rust ports in this family are careful, line-faithful translations. They are byte-identical to
the C over a broad differential sweep and to the fixed C on every UB-tainted case. Wherever the C
reads freed or uninitialised memory, the Rust gives the deterministic answer, with one exception:
the 2D WPC decoder's dead-tile path (H1). There the port deliberately reproduced the zero-heap
outcome, so it outputs NaN for any all-zero 256×256 tile, and a passing golden (wpcuncompress 10)
records those NaNs. The compressor EXCLUDE lists hold up: they are correct and minimal,
and Rust equals the fixed C on them. However, five decoder goldens downstream of the same realloc
bug still encode freed-memory garbage (empty output, all-NaN output, or environment-dependent
garbage) and need regenerating or excluding. Two inherited SU defects make wpccompress unusable at
its default error on real data (H3), and the wt/wpt compressors emit undecodable `nsize=-1` files.
Neither is a port bug, but goldens cannot catch either, so how to handle them is a project
decision. The remaining issues are robustness gaps in the two newest decoders (panic instead of
SIGSEGV) and a silent zero-output on negative `nstage`.
