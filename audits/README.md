# OxWave audits (2026-10-29)

Independent Claude audits of the port, run after every program and library crate had passed its tests.
They looked for what golden tests cannot catch. Each finding below was reproduced against the
original SU binaries before anything was changed.

| Report | Scope | Outcome |
|---|---|---|
| `claude-wave3-misc.md` | spsplot, succepstrum, sufxdecon, sumigpspi, sust | sufxdecon long-trace panic and 2 bad-input panics fixed; sumigpspi tmpdir check fixed; stale goldens regenerated |
| `claude-wave3-compression.md` | dct/wpc/wpt/wt libraries and programs | decoder goldens built on SU UB excluded; zero-tile NaN fixed in lib_wpc2; wt/wpt decoder and negative-nstage panics fixed |
| `claude-wave3-tri-progs.md` | ela*, gbbeam, normray, tri*, uni2tri | negative grid sizes fixed; gbbeam/triseis s[-1] behaviour documented and pinned |
| `claude-wave3-tri-libs.md` | lib-tri, lib-trielas, lib-tetra, tetramod, sutetraray | sutetraray maxnsegs crash and error paths fixed; tetramod ficth panic fixed; lib-trielas 1-ulp rcabs fixed |
| `claude-coverage.md` (+ `-cases.json`) | all 393 programs, ranked by test weakness | **suprofcsv was a stub** (now fully ported); surescsv sort order fixed; 124 new cases added to 29 programs, all passing; ela*/raydata goldens rewritten to run the real computation |
| `claude-ub-exclusions.md` (+ `.json`) | goldens that depend on SU undefined behaviour | exclusion lists for wpc1uncomp2, wpcuncompress; sufctanismod grid cases repaired |
| `claude-sureflpsvsh.md` | the one program still failing its goldens | 8 failing goldens were SU UB (excluded); 7 finite PSV seismogram cases added; now passes |

## Rules applied

- **Undefined behaviour in SU.** A golden is excluded (`goldens/<prog>/EXCLUDE.txt`) only if its
  output changes when the UB is removed by a minimal fix of the C, built with SU's own flags. Where the
  UB happens but the output is unaffected, the case stays. Where the C reads garbage, the Rust does the
  intended thing (e.g. zeros) and says so in a comment.
- **Coverage rule.** Every program needs at least 3 golden cases where the original exits 0 and its
  output is compared. Exemptions with reasons: `goldens/DATA_EXEMPT.txt` (psepsi: needs Ghostscript).

## SU bugs found (see ../DEVIATIONS.md for how OxWave handles each)

- `wpccompress` at default `error=0.01` on noisy data: SU's overflow retry re-encodes an already
  transformed block (wpc.c:203-216/255-268), so the round-trip error reaches ~100%. **Fixed in OxWave
  2026-10-29**: it retries from the original block; see `../DEVIATIONS.md`.
- `wpc1comp2`/`wpccompress` write freed memory after realloc (wpc1.c:230, wpc.c:209/261/280); OxWave
  writes the intended bytes.
- `sureflpsvsh` with `lobs<=1` or `lsource<=1` reads unfilled arrays in SU (NaN output); OxWave rejects
  those inputs. SH mode segfaults in SU when `nx != tsec/dt` (compute.c:230 argument swap).
- `elamodel` exits 1 on success (checkpars before getpar of `fill`); OxWave matches.
