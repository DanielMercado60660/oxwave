# OxWave license audit (CWP/SU source)

Source: github.com/JohnWStockwellJr/SeisUnix, commit 4db4181 (cloned 2026-10-28) into `vendor/su-original/`.
Scope audited: phase 1 directories `src/su/main`, `src/su/lib`, `src/su/include`, `src/cwp/lib`, `src/cwp/include`, `src/par/lib`, `src/par/include`.

## Main license

BSD-3-style (Copyright 2008 Colorado School of Mines). Porting is allowed. Requirements we must keep:
- Every ported file keeps the CWP copyright notice and disclaimer in its header.
- No use of the Colorado School of Mines or contributor names to endorse or promote OxWave. Describe it as "a Rust port of CWP/SU."
- README cites SU in the approved reference format: Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x Release No. __: an open source software package for seismic research and processing, Center for Wave Phenomena, Colorado School of Mines.

## Files needing special handling

| File(s) | Issue | Handling |
|---|---|---|
| `su/main/filters/susvdfilt.c` | Copyright Durham University 2013, "All rights reserved", no grant | Not ported |
| `su/main/synthetics_waveforms_testpatterns/suea2df.c` | Copyright Uppsala University 1998, "All rights reserved", no grant | Not ported |
| `cwp/lib/svd.c` | "Adapted from svdecomp.c in XLISP-STAT, which is code from Numerical Recipes" | Clean-room from the public EISPACK/Golub-Reinsch algorithm; still tested against the C output |
| `cwp/lib/axb.c`, `par/lib/smooth.c` | "Adapted from Numerical Recipes in C" | Clean-room from the cited algorithm, tested against the C output |
| `par/lib/minfunc.c` | "Based on a similar routine in Numerical Recipes", also Netlib | Clean-room from the Netlib/textbook algorithm |
| `cwp/lib/complex.c`, `cwp/lib/complexd.c` | "Adapted from Numerical Recipes (Appendix E)": complex arithmetic | Low risk (textbook operations); write from scratch, no copying |

Clean-room ports were written from the algorithm description and the C function's signature and goldens, without the C body.

## SEG code (decided 2026-10-28: port with attribution)

SEG's terms of use (https://software.seg.org/disclaimer.txt, copy in `LICENSES/SEG-disclaimer.txt`; the live site is behind a bot check, so the
Internet Archive copy of 2023-10-03 was read) let anyone use, modify and redistribute the code, on four conditions:
1. Any publication or software distribution that uses it references the GEOPHYSICS paper and the code's software.seg.org URL.
2. Source files keep a one-line SEG copyright notice, a link to the disclaimer and a link to the code URL.
3. Imported code, original or modified, is clearly marked and attributed (paper, code URL, disclaimer; patents if any - none are mentioned).
4. Modified, redistributed code carries a notice of the nature, extent, purpose, date, author and affiliation of the change.

What is actually SEG's (Werner M. Heigl, Apache Corporation; code at http://software.seg.org/2007/0004):

| C source | SEG part | OxWave home |
|---|---|---|
| `cwp/lib/orthpolynomials.c` | `hermite_n_polynomial` | crate `oxwave-orthogonal-polynomials` |
| `cwp/lib/waveforms.c` | `deriv_n_gauss` only; the other wavelets (Ricker, AKB, Berlage, Gaussian...) are CWP code | crate `oxwave-waveforms` |
| `su/main/synthetics_waveforms_testpatterns/sudgwaveform.c` | the program (Heigl, 2007) | `programs/sudgwaveform` |

Paper: Heigl, W. M., 2007, Computing Gaussian derivative waveforms of any order: Geophysics, 72(4), H39-H42, doi:10.1190/1.2716624.

How OxWave meets the terms: each Rust file holding SEG-derived code carries `tools/seg_header.txt` after the CWP header
(SEG copyright line, code URL, disclaimer link, paper, modification notice), and the SEG code sits between
`// SEG-BEGIN <file>:<function>` and `// SEG-END` lines.
The OxWave README must cite the paper and code URL (condition 1).

`suwaveform.c`, `suattributes.c` and `supolar.c` were held by mistake: they only cite papers published by SEG in their
reference lists, and are plain CWP code (Colorado School of Mines copyright). `suwaveform` calls `deriv_n_gauss` through
`oxwave-waveforms`, so it needs no notice of its own.

## Citation-only (no action)

`cwp/lib/pfafft.c` and `dpfafft.c` (states it "differs significantly" from Numerical Recipes), `cwp/lib/sort.c`, `cwp/lib/linfit.c` ("inspired by"), `su/main/transforms/suslowft.c`, `suslowift.c`, `su/main/multicomponent/sueipofi.c`. Virginia Tech mentions in `segyread.c`, `sucommand.c` and `suhtmath.c` are author credits only.

## Out of scope for phase 1 (separate terms or not CLI processing)

`BHP_SU_V7.0` (own BHP license), `tksu-0.1.27` (GPL v2), `Third_Party/*`, `Sfio`, `Xtcwp`, `Xmcwp`, `xplot`, `xtri`, `Mesa`, `PVM`, `Fortran`, `Mathematica`, `CWPMatlab`. `PotashCorp_SU` and `SLT_SU` use the same CSM license and can be added later.

## Size

286 programs in `su/main` across 27 families, plus 72 files in `cwp/lib`, 19 in `par/lib`, 19 in `su/lib`.
