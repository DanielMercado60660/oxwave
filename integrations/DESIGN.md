# OxWave integrations: Python bindings and MCP server

OxWave is a Rust port of CWP/SU (Seismic Un*x): 398 command-line programs that read
and write SU trace streams (240-byte header + float32 samples, **big-endian/XDR** in this
build) and are chained with pipes. The program catalogue bundled with the Python package
lists 400 SU programs, 398 of them available. The integrations (this root,
`oxwave/integrations/`) add:

* the Python package `oxwave` (pyo3 + maturin, abi3): SU/SEG-Y <-> NumPy, header access
  by SU key name, every OxWave program as a function (`oxwave.run`, `oxwave.pipe`) by
  piping SU bytes through the real binaries, selfdoc lookup, PNG section plots;
* `oxwave.mcp_server` (`oxwave-mcp`): an MCP server (official `mcp` SDK 2.3) so agents
  can list/describe programs, run flows on files of a sandboxed workspace, inspect and plot.

## Layout (paths relative to this root)

```
DESIGN.md                                            this file
pyproject.toml  uv.lock  .venv/                      uv env (Python 3.12): numpy matplotlib mcp pytest maturin ruff
Cargo.toml  Cargo.lock                               cargo workspace (member: py); pyo3 =0.30.1 abi3-py310
py/Cargo.toml  py/pyproject.toml                     extension crate `oxwave-py` (lib `_native`), maturin config
py/src/lib.rs                                        the `oxwave._native` extension
py/python/oxwave/
  __init__.py  _headers.py  traces.py                re-exports, header dtypes, Traces + SU codec
  _native.pyi                                        contract of the Rust extension
  _programs.json                                     program catalogue (from ../inventory)
  su.py segy.py runner.py plot.py                    SU I/O, SEG-Y, program runner, plotting
  mcp_server/__init__.py __main__.py config.py server.py   server, tool schemas, limits
  mcp_server/workspace.py tools_info.py              workspace sandbox, flow validation, info tools
  mcp_server/tools_flow.py                           run_flow, plot_section
judge/tests/  judge/fixtures/                        pytest tests; fixtures made by ORIGINAL SU
tools/                                               helpers (build_ext, stage_bins, make_fixtures, run_demo_flows)
```

References outside this root: `../crates/oxwave-core` (Rust core:
header tables, IBM floats, EBCDIC tables), `../vendor/su-original/src` (original C sources;
`su/main/<family>/<prog>.c`), `../tools/oxharness.py` (comparison rules).

## The contract

The **signatures and docstrings in the stub files are the specification**. The shared
modules (`__init__.py`, `_headers.py`, `traces.py`, `_native.pyi`, and the MCP server's
`__init__.py`, `config.py` and `server.py`) define the common types; the other modules build
on them rather than re-implementing them:

* `oxwave.HEADER_DTYPE`: 240-byte structured dtype, big-endian fields at the segy.h
  offsets plus `unass` (14 shorts), so headers round-trip bit-exactly. `oxwave.KEYS` are
  the 80 SU keys in key-index order. `BHEADER_DTYPE`/`BKEYS`: SEG-Y binary header.
* `oxwave.Traces(data, headers=None, *, dt=None)`: `headers` (ntr,) HEADER_DTYPE,
  `data` float32 (ntr, ns). `tr["offset"]` / `tr["offset"] = ...`, `tr.header(i)`,
  `tr[i:j]`, `Traces.from_bytes(b)` / `tr.to_bytes()` = the reference SU codec with
  fgettr's error messages (`SUFormatError`).
* Every module keeps its `__all__`. Signatures, defaults and docstrings are the public API:
  changing one, or adding a required parameter, is an API change. Private helpers (leading
  `_`) are fine.

## Environment and commands

* Use the venv binaries: `.venv/bin/python`, `.venv/bin/maturin`, `.venv/bin/ruff`.
  The build is offline: `uv run`/`uv sync`/`pip install` would try to re-lock and need
  the network.
* Rust 1.88 or later, `cargo ... --offline`. `uv.lock` and
  `Cargo.lock` pin the dependencies.
* `.venv/bin/python tools/build_ext.py` builds the wheel offline (target/wheels) and puts
  `_native*.so` in target/native/; the tests' conftest loads it as `oxwave._native` and
  imports the Python package from `py/python` (sources, not the wheel).
* OxWave binaries: `$OXWAVE_BIN` = `target/oxbin` (a snapshot of the built programs,
  staged by `tools/stage_bins.py`). Original SU binaries are never run by tests or by `oxwave`.
* Tests need no network:
  `.venv/bin/python tools/build_ext.py && .venv/bin/python tools/stage_bins.py && .venv/bin/python -m pytest -p no:cacheprovider judge/tests`
  (after `cargo build --release` at the repository root)
  (or a single `judge/tests/test_<module>.py`).

## Conventions

1. Faithful to SU: where the contract says "like `<prog>`", the result must be
   byte-identical to the original program on this machine (arm64 macOS, SUXDR build).
   Read the C source when in doubt; the goldens in judge/fixtures decide.
2. Subprocesses: argv lists only, never `shell=True`, never `os.system`; programs are
   launched by absolute path with `argv[0]` = bare program name.
3. Pure Python modules depend only on the stdlib, NumPy and (lazily, plot only)
   matplotlib; only `segy.py` requires `oxwave._native`.
4. MCP tools: every anticipated failure raises `config.ToolError(message)` (other
   exceptions reach the agent as an opaque "Error executing tool"); results are
   JSON-serialisable dicts (no NumPy scalars, no NaN/inf). The server (`server.py`) calls
   the tool functions with a `Workspace` and passes `ToolError` through.
5. File header: Python files keep the first comment line
   (`# Part of OxWave, a Rust port of CWP/SU.`); `py/src/lib.rs` keeps the CWP license
   header (`../tools/cwp_header.txt`). Never use the Colorado School of Mines or
   contributor names to promote OxWave.

## Tolerances

Byte-exact everywhere (codec, SEG-Y, sushw, surange text, run/pipe vs the same binaries
in bash). Where a test compares OxWave output with an ORIGINAL golden (pipes/*.su) it
uses check_program's rule: headers bit-exact, samples within 1e-5 of the trace max.
