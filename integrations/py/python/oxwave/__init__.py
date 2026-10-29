# Part of OxWave, a Rust port of CWP/SU.
"""OxWave for Python: SU/SEG-Y data in NumPy, and every OxWave program as a function.

OxWave is a Rust port of CWP/SU (Seismic Un*x). This package reads and writes SU
(XDR, big-endian) and SEG-Y files into NumPy, gives header access by SU key name, and
runs the OxWave command-line programs over arrays or files by piping SU streams through
the real binaries (exactly the bytes the command line would produce).

Quick tour::

    import oxwave as ox
    tr = ox.run("suplane", nt=64, ntr=8)             # -> Traces (headers + float32 data)
    tr = ox.run("sugain", tr, scale=2.5)              # array in, array out
    tr["offset"] = -200 + 50 * np.arange(tr.ntr)      # header access by SU key
    out = ox.pipe([("sufilter", {"f": [5, 10, 40, 50]}), ("sugain", {"agc": 1})], tr)
    ox.write_su("out.su", out)
    print(ox.surange_text(out))
    print(ox.selfdoc("sunmo"))

Binaries are found in, in order: :func:`set_bin_dir`, ``$OXWAVE_BIN`` (os.pathsep-separated
directories), then ``oxwave/bin`` inside the installed package. Original SU binaries on
``$PATH`` are never used.

Module map (each documents its contract; this file only re-exports):
``traces`` (Traces, codec), ``su`` (SU files, header tools, surange), ``segy`` (SEG-Y),
``runner`` (binaries, run, pipe, selfdoc, catalogue), ``plot`` (section plots),
``mcp_server`` (MCP server for agents), ``_native`` (Rust extension: IBM floats, EBCDIC,
fast SU codec; optional for everything except ``segy``).
"""
from ._headers import (BHEADER_DTYPE, BHEADER_KEYS, BKEYS, HEADER_DTYPE, HEADER_KEYS, KEYS,
                       key_type)
from .traces import SUFormatError, Traces
from .su import concat, gethw, iter_su, read_su, sethw, surange, surange_text, write_su
from .segy import SegyFile, read_segy, write_segy
from .runner import (OutputLimitExceeded, OxWaveError, ProgramError, ProgramNotFound,
                     ProgramTimeout, Repeat, RunResult, UnknownParameterError, available_programs,
                     bin_dirs, catalog, format_argv, help, pipe, program_info, run, selfdoc,
                     set_bin_dir, validate_params, which)

__version__ = "0.1.0"

__all__ = [
    "BHEADER_DTYPE", "BHEADER_KEYS", "BKEYS", "HEADER_DTYPE", "HEADER_KEYS", "KEYS", "key_type",
    "SUFormatError", "Traces",
    "concat", "gethw", "iter_su", "read_su", "sethw", "surange", "surange_text", "write_su",
    "SegyFile", "read_segy", "write_segy",
    "OutputLimitExceeded", "OxWaveError", "ProgramError", "ProgramNotFound", "ProgramTimeout",
    "Repeat", "RunResult", "UnknownParameterError", "available_programs", "bin_dirs", "catalog",
    "format_argv", "help", "pipe", "program_info", "run", "selfdoc", "set_bin_dir",
    "validate_params", "which",
]
