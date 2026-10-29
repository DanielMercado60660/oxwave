# Part of OxWave, a Rust port of CWP/SU.
"""MCP tools that execute: run a processing flow, plot a section.

Contract (stable signatures/docstrings). JSON-serialisable dict results; anticipated
failures raise ``config.ToolError``. Execution goes through ``oxwave.pipe`` (argv lists,
never a shell) with ``cwd = ws.root`` and the workspace's binary lookup.
"""
import math

import oxwave

from .config import ToolError
from .tools_info import inspect_file

__all__ = ["run_flow", "plot_section"]


def run_flow(ws, stages, input=None, output=None, timeout_s=None):
    """Run ``stages`` (see ``Workspace.stage_argv``) as one pipeline.

    ``input``: workspace path of an existing file fed to stage 0's stdin, or None (no
    stdin). ``output``: workspace path for the final stdout (parents created, overwritten
    if it exists); None -> ``"flows/flow-<n><ext>"`` where ``<ext>`` is ``.su`` when the
    last program's catalogue stdout is ``su_traces`` and ``.out`` otherwise, and n is the
    smallest integer >= 1 for which that file does not exist yet. ``timeout_s``: None -> ``ws.limits.default_timeout_s``; must be
    > 0 and <= ``max_timeout_s`` (ToolError otherwise).

    Validation (all before running, ToolError): ``ws.flow_argvs(stages)``, input exists
    inside the workspace, output inside the workspace.
    Execution failures are NOT ToolErrors: a stage exiting non-zero returns
    ``ok: false``. Timeout -> ToolError containing ``"timed out"``; final stdout larger
    than ``ws.limits.max_output_bytes`` -> ToolError containing ``"output limit"``; in
    both cases the partial output file is removed.

    Returns ``{"ok": bool (every stage exited 0), "output": relative path,
    "returncodes": [int...], "argv": [[str...]...] (bare program names),
    "stderr": combined stderr, at most ws.limits.max_stderr_chars (tail kept),
    "elapsed_s": float, "size_bytes": int,
    "summary": inspect_file(ws, output) when ok and the output parses as SU, else None}``.
    """

    argvs = ws.flow_argvs(stages)
    timeout = ws.limits.default_timeout_s if timeout_s is None else timeout_s
    if (isinstance(timeout, bool) or not isinstance(timeout, (int, float))
            or not math.isfinite(timeout) or not 0 < timeout <= ws.limits.max_timeout_s):
        raise ToolError("invalid timeout")
    source = ws.resolve(input, must_exist=True) if input is not None else None
    if source is not None and not source.is_file():
        raise ToolError("input is not a file")
    if output is None:
        ext = ".su" if oxwave.catalog()[argvs[-1][0]].get("stdout") == "su_traces" else ".out"
        n = 1
        while ws.resolve("flows/flow-%d%s" % (n, ext)).exists():
            n += 1
        output = "flows/flow-%d%s" % (n, ext)
    dest = ws.resolve(output, for_write=True)
    try:
        result = oxwave.pipe([(argv[0], None, argv[1:]) for argv in argvs], input=source,
                             output=dest, timeout=timeout, check=False, raw=True,
                             cwd=ws.root, max_output_bytes=ws.limits.max_output_bytes)
    except oxwave.ProgramTimeout as exc:
        dest.unlink(missing_ok=True)
        raise ToolError("flow timed out: %s" % exc) from exc
    except oxwave.OutputLimitExceeded as exc:
        dest.unlink(missing_ok=True)
        raise ToolError("flow output limit exceeded: %s" % exc) from exc
    except (oxwave.OxWaveError, OSError, ValueError) as exc:
        raise ToolError("cannot run flow: %s" % exc) from exc
    ok = all(code == 0 or (code == -13 and i < len(result.returncodes) - 1)
             for i, code in enumerate(result.returncodes))
    summary = None
    if ok:
        try:
            summary = inspect_file(ws, dest)
        except ToolError:
            pass
    limit = ws.limits.max_stderr_chars
    try:
        size = dest.stat().st_size
    except OSError as exc:
        raise ToolError("cannot read flow output: %s" % exc) from exc
    return {"ok": ok, "output": ws.relative(dest), "returncodes": result.returncodes,
            "argv": argvs, "stderr": result.stderr[-limit:] if limit > 0 else "",
            "elapsed_s": result.elapsed, "size_bytes": size, "summary": summary}


def plot_section(ws, path, output=None, kind="image", perc=99.0, clip=None, title=None):
    """Plot an SU file of the workspace to PNG with ``oxwave.plot.plot_section``.

    ``path`` must exist in the workspace and be a valid SU stream (ToolError
    ``"not a valid SU"`` otherwise). ``output``: workspace path (``.png`` appended if it
    lacks that suffix); None -> ``path`` with its suffix replaced by ``.png``.
    ``kind``: ``"image"`` or ``"wiggle"``. Reads at most ``ws.limits.max_inspect_traces``
    traces. Returns ``{"png": relative path, "width_px": int, "height_px": int,
    "ntr": int, "ns": int}`` (800x600 px: the plot module's defaults).
    """

    from .. import plot, su
    from ..traces import SUFormatError

    source = ws.resolve(path, must_exist=True)
    if output is None:
        destination = source.with_suffix(".png")
    else:
        destination = ws.resolve(output)
        if destination.suffix.lower() != ".png":
            destination = destination.with_name(destination.name + ".png")
    destination = ws.resolve(destination)
    if kind not in ("image", "wiggle"):
        raise ToolError("unknown plot kind: %s" % kind)
    try:
        traces = su.read_su(source, max_traces=ws.limits.max_inspect_traces)
        if not traces.ntr:
            raise SUFormatError("no traces")
        plot.clip_value(traces, perc, clip)
        destination = ws.resolve(destination, for_write=True)
        plot.plot_section(traces, destination, kind=kind, perc=perc, clip=clip, title=title)
    except SUFormatError as exc:
        raise ToolError("not a valid SU: %s" % exc) from exc
    except (OSError, ValueError, TypeError) as exc:
        raise ToolError("cannot plot section: %s" % exc) from exc
    return {"png": ws.relative(destination), "width_px": 800, "height_px": 600,
            "ntr": traces.ntr, "ns": traces.ns}
