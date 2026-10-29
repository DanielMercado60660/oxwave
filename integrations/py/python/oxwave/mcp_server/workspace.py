# Part of OxWave, a Rust port of CWP/SU.
"""Workspace sandbox and flow validation for the MCP server.

Contract (stable signatures/docstrings). Every anticipated failure raises
``config.ToolError`` with a message an agent can act on; the substrings quoted below
are checked by the judges.
"""
import math
import os
import re
from collections.abc import Mapping
from pathlib import Path

from .. import runner

from .config import (BLOCKED_PROGRAMS, PARAM_KEY, SAFE_SHELL_VALUE, SHELL_SENSITIVE_PROGRAMS,  # noqa: F401
                     Limits, ToolError)


__all__ = ["Workspace"]


class Workspace:
    """A directory the server may read and write, plus limits and the binary lookup.

    ``Workspace(root, *, bin_dir=None, limits=None)``: ``root`` must be an existing
    directory (else ValueError); it is stored resolved (``Path.resolve()``) as
    ``self.root``. ``bin_dir`` (directory or list of directories, or None) is stored as
    ``self.bin_dir`` and passed to ``oxwave.set_bin_dir`` by the server before running
    anything (None = keep ``$OXWAVE_BIN`` lookup). ``limits`` defaults to ``Limits()``.
    """

    def __init__(self, root, *, bin_dir=None, limits=None):
        try:
            root = Path(root).resolve()
            if not root.is_dir():
                raise ValueError("workspace root must be an existing directory")
        except (OSError, TypeError, RuntimeError) as exc:
            raise ValueError("workspace root must be an existing directory") from exc
        self.root = root
        self.bin_dir = bin_dir
        self.limits = limits if limits is not None else Limits()

    def resolve(self, path, *, must_exist=False, for_write=False):
        """Map an agent-supplied path to an absolute ``Path`` inside ``root``.

        Relative paths are taken from ``root``; absolute paths are accepted only if they
        lie inside ``root``. The result is checked after resolving symlinks and ``..``:
        anything outside ``root`` -> ToolError containing ``"outside the workspace"``.
        Empty path, NUL bytes, or a leading ``~`` -> ToolError. ``must_exist`` and missing
        -> ToolError containing ``"not found"``. ``for_write``: parent directories are
        created (inside root); an existing directory at that path -> ToolError.
        """

        try:
            raw = os.fspath(path)
            if not isinstance(raw, str) or not raw or "\0" in raw or raw.startswith("~"):
                raise ToolError("invalid workspace path")
            resolved = (self.root / raw).resolve()
            if not resolved.is_relative_to(self.root):
                raise ToolError("path outside the workspace: %s" % raw)
            if must_exist and not resolved.exists():
                raise ToolError("path not found: %s" % raw)
            if for_write:
                if resolved.is_dir():
                    raise ToolError("write path is a directory: %s" % raw)
                resolved.parent.mkdir(parents=True, exist_ok=True)
            return resolved
        except ToolError:
            raise
        except (OSError, ValueError, TypeError, RuntimeError) as exc:
            raise ToolError("invalid workspace path: %s" % exc) from exc

    def relative(self, path):
        """POSIX-style path of ``path`` (absolute, inside root) relative to ``root``."""

        try:
            return Path(path).relative_to(self.root).as_posix()
        except ValueError as exc:
            raise ToolError("path outside the workspace") from exc

    def check_program(self, name):
        """Validate a program name for running: must be a catalogue program
        (ToolError ``"unknown program"``), not in BLOCKED_PROGRAMS (ToolError
        ``"not allowed"``), and available as an OxWave binary (ToolError
        ``"not available"``, mentioning its catalogue status such as BLOCKED)."""

        catalog = runner.catalog()
        if not isinstance(name, str) or name not in catalog:
            raise ToolError("unknown program: %s" % name)
        if name in BLOCKED_PROGRAMS:
            raise ToolError("program not allowed: %s" % name)
        try:
            runner.which(name)
        except runner.ProgramNotFound as exc:
            raise ToolError("program not available: %s (status %s)" %
                            (name, catalog[name].get("status"))) from exc

    def stage_argv(self, stage):
        """Validate one flow stage and return its argv list (``[program, *args, *params]``,
        program as a bare name; encoding = ``oxwave.format_argv``).

        ``stage`` is a mapping with ``program`` (str, required), optional ``params``
        (mapping key -> str | int | float | list of those; lists are comma-joined),
        optional ``repeat`` (mapping key -> list of values; emitted after ``params`` as
        repeated ``key=value``), optional ``args`` (list of str). Other stage keys ->
        ToolError. Rules (ToolError on violation):

        * ``check_program``;
        * keys match PARAM_KEY; a key both in ``params`` and ``repeat`` -> error;
        * if the catalogue entry has a non-empty ``params`` list, every key must be in it
          (``par`` always allowed): message contains ``"unknown parameter"`` and the key;
        * ``bool`` or ``None`` values, nested lists, non-finite floats -> error;
        * string values (and list elements) containing NUL or newline -> error;
        * path-like string values -- containing ``/``, starting with ``.`` or ``~``, or the
          value of ``par`` -- must ``resolve`` inside the workspace (``"outside the
          workspace"``); they are passed through unchanged (stages run with cwd = root);
        * for SHELL_SENSITIVE_PROGRAMS every string value must fully match
          SAFE_SHELL_VALUE (``"unsafe value"``);
        * each ``args`` word must not start with ``-`` and must ``resolve`` inside the
          workspace (it names a file).
        """

        if not isinstance(stage, Mapping) or set(stage) - {"program", "params", "repeat", "args"}:
            raise ToolError("invalid stage")
        name = stage.get("program")
        self.check_program(name)
        params = stage.get("params", {})
        repeat = stage.get("repeat", {})
        args = stage.get("args", [])
        if not isinstance(params, Mapping) or not isinstance(repeat, Mapping) or not isinstance(args, list):
            raise ToolError("invalid stage parameters or arguments")
        known = runner.catalog()[name].get("params", [])
        allowed = {p["name"] for p in known} | {"par"}
        if set(params) & set(repeat):
            raise ToolError("parameter appears in both params and repeat")

        def validate_value(value, key):
            if isinstance(value, list):
                for item in value:
                    if isinstance(item, list):
                        raise ToolError("nested parameter list")
                    validate_value(item, key)
            elif isinstance(value, bool) or not isinstance(value, (str, int, float)):
                raise ToolError("invalid parameter value for %s" % key)
            elif isinstance(value, float) and not math.isfinite(value):
                raise ToolError("non-finite parameter value for %s" % key)
            elif isinstance(value, str):
                if "\0" in value or "\n" in value or "\r" in value:
                    raise ToolError("invalid parameter value for %s" % key)
                if name in SHELL_SENSITIVE_PROGRAMS and re.fullmatch(SAFE_SHELL_VALUE, value) is None:
                    raise ToolError("unsafe value for %s" % key)
                if "/" in value or value.startswith((".", "~")) or key == "par":
                    self.resolve(value)

        combined = dict(params)
        for key, values in repeat.items():
            if not isinstance(values, list):
                raise ToolError("repeat values must be a list")
            combined[key] = runner.Repeat(values)
        for key, value in combined.items():
            if not isinstance(key, str) or re.fullmatch(PARAM_KEY, key) is None:
                raise ToolError("invalid parameter name: %s" % key)
            if known and key not in allowed:
                raise ToolError("unknown parameter: %s" % key)
            if isinstance(value, runner.Repeat):
                for item in value:
                    validate_value(item, key)
            else:
                validate_value(value, key)
        for arg in args:
            if not isinstance(arg, str) or arg.startswith("-") or "\n" in arg or "\r" in arg:
                raise ToolError("invalid positional argument")
            self.resolve(arg)
        try:
            return runner.format_argv(name, combined, args)
        except (ValueError, TypeError) as exc:
            raise ToolError(str(exc)) from exc

    def flow_argvs(self, stages):
        """Validate a whole flow: ``stages`` must be a non-empty list of at most
        ``limits.max_stages`` stages (ToolError otherwise); returns the list of argv lists.
        Nothing is executed."""

        if not isinstance(stages, list) or not stages or len(stages) > self.limits.max_stages:
            raise ToolError("flow must contain between 1 and %d stages" % self.limits.max_stages)
        return [self.stage_argv(stage) for stage in stages]
