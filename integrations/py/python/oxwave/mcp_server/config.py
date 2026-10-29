# Part of OxWave, a Rust port of CWP/SU.
"""Fixed policy of the MCP server: limits, program blacklist, error type."""
from dataclasses import dataclass

from mcp.server.mcpserver.exceptions import ToolError  # re-exported: the ONLY anticipated-error type

__all__ = ["Limits", "ToolError", "BLOCKED_PROGRAMS", "SHELL_SENSITIVE_PROGRAMS", "SAFE_SHELL_VALUE",
           "PARAM_KEY", "TOOL_NAMES"]


@dataclass(frozen=True)
class Limits:
    """Resource limits of one server instance."""
    max_output_bytes: int = 256 * 1024 * 1024   # final stdout of a flow / any written file
    default_timeout_s: float = 60.0             # run_flow timeout when the call gives none
    max_timeout_s: float = 600.0                # larger requested timeouts -> ToolError
    max_stages: int = 32                        # stages per flow
    max_stderr_chars: int = 8000                # stderr returned per flow (tail kept)
    max_inspect_traces: int = 100000            # inspect_file reads at most this many traces


#: Programs the server never runs: they spawn a shell with user-controlled text or are
#: interactive (sucommand: popen of an arbitrary command; suedit/suxedit: interactive, system()).
BLOCKED_PROGRAMS = frozenset({"sucommand", "suedit", "suxedit"})

#: Programs that build a shell command line from parameter values (epopen of `dd ...`).
#: For these every string parameter value must fully match SAFE_SHELL_VALUE.
SHELL_SENSITIVE_PROGRAMS = frozenset({"segyread", "segywrite", "suaddnoise"})
SAFE_SHELL_VALUE = r"[A-Za-z0-9_.,:+/=-]*"

#: Valid parameter key.
PARAM_KEY = r"[A-Za-z_][A-Za-z0-9_]*"

#: The server's tools, in registration order (frozen names).
TOOL_NAMES = ("list_programs", "describe_program", "run_flow", "inspect_file", "plot_section",
              "list_files")
