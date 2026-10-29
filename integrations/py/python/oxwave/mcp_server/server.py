# Part of OxWave, a Rust port of CWP/SU.
"""MCP protocol layer: registers the six tools on an ``mcp`` ``MCPServer``.

This file defines the contract: tool names, argument schemas and descriptions.
The tool bodies live in ``tools_info`` and ``tools_flow``;
they raise ``ToolError`` for anticipated failures (the message reaches the agent) and
return JSON-serialisable dicts (sent as JSON text content).
"""
import argparse
import os
import sys
from typing import Any, Literal, Optional, Union

from typing_extensions import NotRequired, TypedDict

from mcp.server.mcpserver import MCPServer

from .. import runner
from . import tools_flow, tools_info
from .config import Limits, ToolError
from .workspace import Workspace

__all__ = ["make_server", "main", "Stage"]

Scalar = Union[str, int, float]
ParamValue = Union[str, int, float, list[Scalar]]


class Stage(TypedDict):
    """One program of a flow."""
    program: str
    params: NotRequired[dict[str, ParamValue]]
    repeat: NotRequired[dict[str, list[ParamValue]]]
    args: NotRequired[list[str]]


INSTRUCTIONS = """OxWave: a Rust port of CWP/SU (Seismic Un*x). Programs read and write SU trace
streams (240-byte headers + float32 samples) and are chained like a Unix pipeline.
Use list_programs / describe_program to find programs and their key=value parameters,
run_flow to run a pipeline on files of the workspace (all paths are relative to it),
inspect_file for header ranges and statistics, plot_section for a PNG.
Synthetic data: e.g. stages [{"program": "suplane", "params": {"nt": 64, "ntr": 32}}]."""


def make_server(workspace, *, bin_dir=None, limits=None, name="oxwave"):
    """Build an MCPServer bound to ``workspace`` (existing directory). ``bin_dir``: OxWave
    binary directory/directories (None = ``$OXWAVE_BIN``). ``limits``: :class:`Limits`."""
    ws = Workspace(workspace, bin_dir=bin_dir, limits=limits or Limits())
    srv = MCPServer(name, instructions=INSTRUCTIONS)

    def _bins():
        if ws.bin_dir is not None:
            runner.set_bin_dir(ws.bin_dir)

    @srv.tool()
    def list_programs(filter: Optional[str] = None, available_only: bool = True) -> dict[str, Any]:
        """List OxWave programs with a one-line purpose. `filter`: case-insensitive substring of
        name, purpose or family. `available_only`: only programs that can run here."""
        _bins()
        return tools_info.list_programs(ws, filter=filter, available_only=available_only)

    @srv.tool()
    def describe_program(name: str) -> dict[str, Any]:
        """Parameters (name, type, default, description), I/O kinds and the program's own
        self-documentation text."""
        _bins()
        return tools_info.describe_program(ws, name)

    @srv.tool()
    def run_flow(stages: list[Stage], input: Optional[str] = None, output: Optional[str] = None,
                 timeout_s: Optional[float] = None) -> dict[str, Any]:
        """Run a pipeline stage0 | stage1 | ... inside the workspace. Each stage:
        {"program": "sugain", "params": {"agc": 1, "wagc": 0.2}, "repeat": {"key": [..]},
        "args": ["file.su"]}; list values are comma-joined (f=5,10,40,50). `input`: workspace
        file for stdin; `output`: workspace file for stdout (default flows/flow-N.su).
        Returns ok, output path, return codes, stderr and a summary of the output."""
        _bins()
        return tools_flow.run_flow(ws, stages, input=input, output=output, timeout_s=timeout_s)

    @srv.tool()
    def inspect_file(path: str, max_traces: Optional[int] = None) -> dict[str, Any]:
        """Header ranges (like surange), first header and sample statistics of an SU file."""
        _bins()
        return tools_info.inspect_file(ws, path, max_traces=max_traces)

    @srv.tool()
    def plot_section(path: str, output: Optional[str] = None,
                     kind: Literal["image", "wiggle"] = "image", perc: float = 99.0,
                     clip: Optional[float] = None, title: Optional[str] = None) -> dict[str, Any]:
        """Render an SU file as a PNG (grey image or wiggles) in the workspace; returns its path."""
        _bins()
        return tools_flow.plot_section(ws, path, output=output, kind=kind, perc=perc, clip=clip,
                                       title=title)

    @srv.tool()
    def list_files(subdir: str = "") -> dict[str, Any]:
        """Files in the workspace (or a subdirectory) with size and kind (su, segy, png, ...)."""
        _bins()
        return tools_info.list_files(ws, subdir=subdir)

    srv.oxwave_workspace = ws  # for tests and embedding
    return srv


def main(argv=None):
    """``oxwave-mcp --workspace DIR [--bin DIR] [--timeout S] [--max-output-mb N]`` (stdio)."""
    ap = argparse.ArgumentParser(prog="oxwave-mcp", description="OxWave MCP server (stdio)")
    ap.add_argument("--workspace", default=os.environ.get("OXWAVE_WORKSPACE"),
                    help="directory the server may read/write (default $OXWAVE_WORKSPACE)")
    ap.add_argument("--bin", default=None, help="OxWave binary directory (default $OXWAVE_BIN)")
    ap.add_argument("--timeout", type=float, default=Limits.default_timeout_s,
                    help="default run_flow timeout in seconds")
    ap.add_argument("--max-output-mb", type=float, default=Limits.max_output_bytes / 2**20)
    a = ap.parse_args(argv)
    if not a.workspace:
        ap.error("--workspace (or $OXWAVE_WORKSPACE) is required")
    os.makedirs(a.workspace, exist_ok=True)
    limits = Limits(max_output_bytes=int(a.max_output_mb * 2**20), default_timeout_s=a.timeout)
    srv = make_server(a.workspace, bin_dir=a.bin, limits=limits)
    srv.run()  # stdio transport
    return 0


if __name__ == "__main__":
    sys.exit(main())


_ = ToolError  # re-export for type checkers
