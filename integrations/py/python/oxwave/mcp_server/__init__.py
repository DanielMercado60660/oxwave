# Part of OxWave, a Rust port of CWP/SU.
"""MCP server that lets agents (Claude Code, Codex, ...) build and run OxWave processing
flows inside a sandboxed workspace directory. Requires the ``mcp`` package (``oxwave[mcp]``).

Entry points: ``oxwave-mcp --workspace DIR [--bin DIR]`` (stdio) or
``python -m oxwave.mcp_server ...``; in-process: :func:`make_server`.
"""
from .config import Limits

__all__ = ["Limits", "make_server", "main"]


def make_server(*args, **kwargs):
    from .server import make_server as _make
    return _make(*args, **kwargs)


def main(argv=None):
    from .server import main as _main
    return _main(argv)
