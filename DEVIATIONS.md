# Where OxWave deliberately differs from SU

OxWave reproduces CWP/SU's behaviour exactly (same bytes, exit codes, messages), including its
quirks. The exceptions are listed here. Each is either a bug fix the maintainers chose to make, or
a place where SU's own behaviour is undefined (it reads freed or uninitialised memory, or writes
out of bounds), so there is nothing well-defined to copy. Golden cases that depended on undefined
behaviour are listed in `goldens/<prog>/EXCLUDE.txt`; the audits behind them are in `audits/`.

## Bug fixes (SU's behaviour is well-defined but wrong)

| Program / library | SU behaviour | OxWave behaviour | Since |
|---|---|---|---|
| `wpccompress` (`oxwave-lib-wpc2`, `wpcCompress`) | When a block's encoding overflows the output buffer (compression ratio below about 4, e.g. noisy data at the default `error=0.01`), SU grows the buffer and re-encodes the block it has already transformed in place (`comp/dwpt/2d/lib/wpc.c:203-216, 255-268`). The round trip then has ~100% relative error. | Re-encodes from the original block. On a 512x128 noisy gather at `error=0.01`: SU 100% round-trip error, OxWave 0.81% (and a smaller file). Output is byte-identical to SU whenever no overflow happens. Tests: `crates/oxwave-lib-wpc2/tests/roundtrip.rs`. | 2026-10-29 |

## Undefined behaviour in SU (OxWave does the intended thing)

| Program / library | SU behaviour | OxWave behaviour |
|---|---|---|
| `wpc1comp2`, `wpccompress` | After `realloc`, the struct keeps the old pointer, so the output contains freed memory (`wpc1.c:230`, `wpc.c:209/261/280`). | Writes the encoded bytes. |
| `wpcuncompress` (`oxwave-lib-wpc2`) | `codeDecoder` never stores an all-zero tile's flag and average; `quantDe` reads them uninitialised and writes NaN. | All-zero tiles decode to 0. |
| `sureflpsvsh` | With `lobs<=1` or `lsource<=1` it reads unfilled free-surface arrays (NaN output); SH mode overruns when `nx != tsec/dt` (`compute.c:230` argument swap). | Rejects those inputs with an error. |
| `gbbeam`, `triseis` | Read `s[-1]` (heap garbage) for the first receiver position. | Treats it as 0 (linear receiver spacing). |
| `tetramod` | Short grid-file names corrupt getpar's argument buffer, so later file parameters are silently ignored. | Reads every named file. |
| `sutetraray` | `ifwf2dump`/`nwf2dump` uninitialised when no `wffile` is given (visible with `verbose=1`). | Zero. |
| Several decoders and compressors (`dctcomp`, `wptcomp`, `wtcomp`, `wpc1uncomp2`) | Short reads or skipped transforms leave arrays uninitialised. | Zero-initialised. On macOS SU happens to get zeroed memory too, so outputs match today. |

## Not changed (documented SU quirks OxWave copies)

- `elamodel` (and other Trielas programs) exit with status 1 on success, because SU calls
  `checkpars()` before reading `fill`.
- Programs whose C crashes on some inputs (e.g. `sufctanismod` in some modes, `trimodel` on a
  curved edge with `maxangle=0.5`) are matched where the crash is deterministic, e.g. by raising
  the same signal.
