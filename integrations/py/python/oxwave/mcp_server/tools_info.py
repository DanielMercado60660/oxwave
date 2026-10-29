# Part of OxWave, a Rust port of CWP/SU.
"""Read-only MCP tools: program catalogue, file inspection, workspace listing.

Contract (stable signatures/docstrings). Each function takes the server's
:class:`~oxwave.mcp_server.workspace.Workspace` first and returns a JSON-serialisable
dict (no NumPy types, no NaN/inf: use None). Anticipated failures raise
``config.ToolError``.
"""
import math
import os

import numpy as np

from .. import runner, su
from ..traces import SUFormatError
from .config import BLOCKED_PROGRAMS, ToolError

__all__ = ["list_programs", "describe_program", "inspect_file", "list_files"]


def list_programs(ws, filter=None, available_only=True):
    """Catalogue listing.

    Returns ``{"count": n, "programs": [{"name", "purpose", "family", "stdin", "stdout",
    "status", "available"}...]}`` sorted by name. ``available`` is True when the OxWave
    binary is found and the program is not in BLOCKED_PROGRAMS. ``available_only``
    keeps only those. ``filter`` (case-insensitive substring) keeps programs whose name,
    purpose or family contains it.
    """

    available = set(runner.available_programs()) - BLOCKED_PROGRAMS
    query = filter.lower() if filter is not None else None
    programs = []
    for name, entry in sorted(runner.catalog().items()):
        found = name in available
        if available_only and not found:
            continue
        if query is not None and not any(query in (entry.get(k) or "").lower()
                                         for k in ("name", "purpose", "family")):
            continue
        item = {key: entry.get(key) for key in ("name", "purpose", "family", "stdin", "stdout", "status")}
        item["available"] = found
        programs.append(item)
    return {"count": len(programs), "programs": programs}


def describe_program(ws, name):
    """Everything an agent needs to call a program.

    Returns ``{"name", "purpose", "family", "status", "available", "stdin", "stdout",
    "params": [{"name", "type", "default", "required", "description"}...],
    "params_known": bool (catalogue params non-empty), "related_programs": [...],
    "selfdoc": str or None}``; ``selfdoc`` is ``oxwave.selfdoc(name)`` when available
    (exact text), else None. Unknown name -> ToolError ``"unknown program"``.
    """

    entry = runner.catalog().get(name)
    if entry is None:
        raise ToolError("unknown program: %s" % name)
    try:
        runner.which(name)
        available = name not in BLOCKED_PROGRAMS
    except runner.ProgramNotFound:
        available = False
    result = {key: entry.get(key) for key in
              ("name", "purpose", "family", "status", "stdin", "stdout", "params", "related_programs")}
    result["available"] = available
    result["params_known"] = bool(entry["params"])
    try:
        result["selfdoc"] = runner.selfdoc(name) if available else None
    except (runner.OxWaveError, OSError) as exc:
        raise ToolError("unable to read selfdoc: %s" % exc) from exc
    return result


def inspect_file(ws, path, max_traces=None):
    """Summary of an SU file in the workspace (like surange + basic stats).

    ``path`` is resolved with ``ws.resolve(path, must_exist=True)``. Reads at most
    ``min(max_traces or inf, ws.limits.max_inspect_traces)`` traces. Returns
    ``{"path": relative posix path, "format": "su", "size_bytes", "ntr", "ns",
    "dt": seconds or None, "truncated": bool (stopped early),
    "ranges": {key: [min, max]} (exactly the keys/values oxwave.surange gives),
    "first_header": {key: value} (nonzero keys of trace 0),
    "stats": {"min", "max", "mean", "rms", "nan_count"} over all samples read (float64
    arithmetic; None when not finite), "surange": oxwave.surange_text(traces)}``.
    A file that is not a valid SU stream -> ToolError containing ``"not a valid SU"``
    and the SUFormatError message.
    """

    resolved = ws.resolve(path, must_exist=True)
    if max_traces is not None and (isinstance(max_traces, bool) or not isinstance(max_traces, int)
                                   or max_traces < 0):
        raise ToolError("max_traces must be a nonnegative integer")
    limit = min(max_traces or ws.limits.max_inspect_traces, ws.limits.max_inspect_traces)
    if limit <= 0:
        raise ToolError("inspection trace limit must be positive")
    try:
        with resolved.open("rb") as f:
            traces = su.read_su(f, max_traces=limit)
            truncated = traces.ntr == limit and len(f.read(240)) >= 240
        size = resolved.stat().st_size
        ranges = su.surange(traces)["keys"]
        samples = traces.data.astype(np.float64)
        with np.errstate(invalid="ignore", over="ignore"):
            stats = {"min": float(np.min(samples)), "max": float(np.max(samples)),
                     "mean": float(np.mean(samples)), "rms": float(np.sqrt(np.mean(samples * samples))),
                     "nan_count": int(np.isnan(samples).sum())}
        stats = {k: (None if isinstance(v, float) and not math.isfinite(v) else v)
                 for k, v in stats.items()}
        def clean(value):
            if isinstance(value, (tuple, list)):
                return [clean(v) for v in value]
            if isinstance(value, float) and not math.isfinite(value):
                return None
            return value
        header = {k: clean(v) for k, v in traces.header(0).items() if v != 0}
        return {"path": ws.relative(resolved), "format": "su", "size_bytes": size,
                "ntr": traces.ntr, "ns": traces.ns, "dt": traces.dt, "truncated": truncated,
                "ranges": {k: clean(v) for k, v in ranges.items()}, "first_header": header,
                "stats": stats, "surange": su.surange_text(traces)}
    except SUFormatError as exc:
        raise ToolError("not a valid SU: %s" % exc) from exc
    except (OSError, ValueError) as exc:
        raise ToolError("cannot inspect file: %s" % exc) from exc


def list_files(ws, subdir=""):
    """Files under ``subdir`` (recursive, resolved inside the workspace; must exist).

    Returns ``{"root": ".", "files": [{"path": relative posix path, "size": bytes,
    "kind": "su" | "segy" | "png" | "text" | "other"}...]}`` sorted by path; ``kind`` by
    extension (``.su``; ``.sgy``/``.segy``; ``.png``; ``.txt``/``.par``/``.json``/``.csv``).
    Hidden entries (name starting with ``.``) are skipped; at most 1000 entries.
    """

    directory = ws.resolve(subdir or ".", must_exist=True)
    if not directory.is_dir():
        raise ToolError("subdir is not a directory")
    kinds = {".su": "su", ".sgy": "segy", ".segy": "segy", ".png": "png",
             ".txt": "text", ".par": "text", ".json": "text", ".csv": "text"}
    files = []
    try:
        for root, dirs, names in os.walk(directory, followlinks=False):
            dirs[:] = sorted(d for d in dirs if not d.startswith(".") and not (os.path.islink(os.path.join(root, d))))
            for name in sorted(names):
                if name.startswith("."):
                    continue
                p = ws.resolve(os.path.join(root, name), must_exist=True)
                if p.is_file():
                    files.append({"path": ws.relative(p), "size": p.stat().st_size,
                                  "kind": kinds.get(p.suffix.lower(), "other")})
        files.sort(key=lambda item: item["path"])
        return {"root": ".", "files": files[:1000]}
    except OSError as exc:
        raise ToolError("cannot list files: %s" % exc) from exc
