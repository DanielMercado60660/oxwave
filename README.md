# OxWave

A Rust port of CWP/SU (Seismic Un*x) that reproduces the original programs' output.

## What it is

OxWave is a line-by-line Rust translation of the CWP/SU programs and the C libraries behind
them. Each binary has the same name, parameters, selfdoc, error messages and exit codes as its
SU counterpart, so pipelines that use ported programs run without changes to parameters or
file formats.

It ports CWP/SU release 45R00 (30 April 2023) as of SeisUnix commit `4db4181`
([github.com/JohnWStockwellJr/SeisUnix](https://github.com/JohnWStockwellJr/SeisUnix), September 2024): 398 programs, which is every compiled program of a standard SU
`make install` except two held back for licence reasons. The Rust workspace uses no external
crates.
It is not free of `unsafe`: 91 programs and 8 library crates use it, mostly for libc calls
(`signal`, `raise`, `access`) and for viewing sample buffers as bytes.

## Status and fidelity

- **493 crates pass their tests**: 398 programs and 95 library crates. Library crates are tested
  differentially against the original C functions.
- Each program was compared with the original SU binary on golden cases (selfdoc plus pipelines
  built from SU generators such as `suplane` and `suspike`). Rules, from
  [`tools/check_program.py`](tools/check_program.py): trace headers bit-identical, samples
  within `1e-5 × max|sample|` of each trace; text exact; other bytes exact; selfdoc exact;
  exit code always; first stderr line when the original fails.
- **SU's own demos**, re-run on the final binaries (122 directories with non-X11 scripts):
  96 MATCH, 0 DIFF, 14 NO_OUTPUT, 9 ORIG_FAILS (original SU fails too), 2 MATCH_ORIG_ONLY,
  1 TIMEOUT. The matching flows ran 129 distinct OxWave programs; SU shell scripts they called
  (`susum`, `merge2`, `merge4`, `sudiff`) came from original SU.
  See [`integrations/demo_flows/REPORT.md`](integrations/demo_flows/REPORT.md).
- After every program and crate had passed, independent audits fixed what golden tests missed (for example,
  `suprofcsv` had been a stub): [`audits/README.md`](audits/README.md).
- **Differential testing round (November 2026):** 10,971 more golden cases, so that every parameter a
  program reads is exercised (boundary sizes, header-driven paths, error exits). They found 291 places where
  OxWave still differed from SU in 70 programs (most often a panic where SU exits cleanly, a different exit
  code or diagnostic, or edge-case numerics); all but one were fixed in the ports. Shared-library fixes
  followed (C-style scanf/strtod/printf, stdio buffering and flush order, mixed packed trace output, and a
  heap use-after-free in the triangulation library that SU's C also has). Cases whose SU output depends on
  undefined behaviour are listed with evidence in `goldens/<prog>/UB_CONFIRMED.txt`; the 14 cases still
  known to differ are in `goldens/<prog>/PENDING_LIB.txt`.
- **Intentional differences** are all in [`DEVIATIONS.md`](DEVIATIONS.md): one SU bug fixed on
  purpose (`wpccompress`); the rest are places where SU has undefined behaviour.

Tested on macOS (Apple silicon) only. If you build and test it on Linux or another platform, please [open an issue](https://github.com/DanielMercado60660/oxwave/issues) with the result.
<!-- TODO(Daniel): test on a Linux VM. -->

## Build

Nothing is on crates.io or PyPI yet. Build from source with Rust 1.88 or later:

```sh
git clone https://github.com/DanielMercado60660/oxwave
cd oxwave
cargo build --release
export PATH="$PWD/target/release:$PATH"
```

`oxwave-cref`, the test-only link to the original C libraries, builds as an empty crate (with a
warning) until a reference SU build is present; see [Testing](#testing-against-reference-su).

- **The binaries use SU's names** (`suplane`, `sugain`, ...). Whichever of OxWave and SU comes
  first on your `PATH` shadows the other.
- **Byte order:** like SU built with `-DSUXDR`, OxWave reads and writes big-endian (XDR) SU
  files only.

## Example

```sh
suplane npl=2 nt=128 ntr=24 | sufilter f=8,12,40,60 | sugain agc=1 wagc=0.2 > out.su
surange < out.su
```

On the development machine (macOS arm64), `out.su` is byte-identical to what the original SU
binaries produce for the same pipeline.

## Python bindings and MCP server

[`integrations/`](integrations/) contains a Python package, `oxwave` (Python 3.10+), and an
MCP server, `oxwave-mcp`.

- `oxwave` reads and writes SU and SEG-Y files as NumPy arrays, gives header access by SU key
  name, and runs any OxWave program as a function (`oxwave.run`, `oxwave.pipe`) by piping SU
  bytes through the binaries. A small PyO3 extension handles IBM floats, EBCDIC and SU decoding.
- `oxwave-mcp` lets AI assistants list and describe programs, run flows on files in one
  workspace directory, inspect files and plot sections. It never runs `sucommand`, `suedit`
  or `suxedit`.

```sh
pip install maturin
maturin build --release -m integrations/py/Cargo.toml --out dist
pip install "$(ls dist/oxwave-*.whl)[mcp]"
export OXWAVE_BIN="$PWD/target/release"
oxwave-mcp --workspace ./ox-workspace
```

```python
import oxwave as ox
tr = ox.run("suplane", npl=2, nt=128, ntr=24)
out = ox.pipe([("sufilter", {"f": [8, 12, 40, 60]}), ("sugain", {"agc": 1, "wagc": 0.2})], tr)
print(ox.surange_text(out))
```

Binaries are found only through `OXWAVE_BIN` (or `oxwave.set_bin_dir`); `PATH` is never
searched and original SU is never used. Details: [`integrations/DESIGN.md`](integrations/DESIGN.md).

## Testing against reference SU

1. Build the original SU into `vendor/su-build` as in [`vendor/BUILD_NOTES.md`](vendor/BUILD_NOTES.md).
   `cargo test` links it through `oxwave-cref` (or set `OXWAVE_SU_LIB`).
2. Record goldens: `python3 tools/make_goldens.py <program>`. Only the case lists
   (`goldens/<program>/cases.txt`) and the program catalogue (`inventory/`) are in the repository.
3. Compare: `cargo build --release -p <program> && python3 tools/check_program.py <program>`
   (exit 0 means every case passed).

## Not ported

- X11, Motif and OpenGL viewers (`suxwigb`, `suximage`, `xwigb`, ...). The PostScript plotting
  programs (`psimage`, `supswigb`, `pswigb`, ...) are ported.
- SU's shell scripts (55 in the reference build, e.g. `suagc`, `susum`, `sudiff`, `merge2`).
- `susvdfilt` and `suea2df` (licence, below).
- BHP_SU, tksu, Third_Party, Sfio, PVM, Fortran, Mathematica, CWPMatlab, PotashCorp_SU and SLT_SU.

## Licence and credits

OxWave is distributed under CWP/SU's BSD-3-Clause licence (Copyright (c) 2008, Colorado School
of Mines); see [LICENSE](LICENSE). Ported files keep the CWP copyright notice and disclaimer.
OxWave is not affiliated with or endorsed by the Colorado School of Mines or the SU contributors.

CWP/SU was written at the Center for Wave Phenomena, Colorado School of Mines, by Jack K. Cohen,
John W. Stockwell Jr. and many contributors. Please cite SU in its approved format:

> Cohen, J. K. and Stockwell, Jr. J. W., (2023), CWP/SU: Seismic Un*x Release No. 45R00: an open
> source software package for seismic research and processing, Center for Wave Phenomena,
> Colorado School of Mines.

**SEG code.** `sudgwaveform`, `deriv_n_gauss` (crate `oxwave-waveforms`) and
`hermite_n_polynomial` (crate `oxwave-orthogonal-polynomials`) translate code Copyright (c) 2007
Society of Exploration Geophysicists, used under SEG's terms
([`LICENSES/SEG-disclaimer.txt`](LICENSES/SEG-disclaimer.txt)) and marked in the source.
Heigl, W. M., 2007, Computing Gaussian derivative waveforms of any order: Geophysics, 72(4),
H39-H42, doi:10.1190/1.2716624; code: http://software.seg.org/2007/0004.

**Excluded and rewritten.** `susvdfilt` (Durham University) and `suea2df` (Uppsala University)
are "All rights reserved" with no licence grant, so they are not ported. Routines SU adapted
from Numerical Recipes (`svd`, `axb`, `smooth`, `minfunc`, complex arithmetic) were rewritten
from the published algorithms and tested against SU. See [`LICENSE_AUDIT.md`](LICENSE_AUDIT.md).
