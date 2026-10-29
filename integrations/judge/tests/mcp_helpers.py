# In-process MCP client helpers for the mcp-core / mcp-flow judges.
import json
import shutil

import anyio

from conftest import FIX, OXBIN


def make_ws(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    for n in ("plane.su", "spike.su", "gather.su", "one.su"):
        shutil.copy(FIX / n, ws / n)
    (ws / "notes.txt").write_text("not seismic\n")
    return ws


def call(server, tool, args):
    """Call one tool through a real MCP client session; returns (is_error, payload) where
    payload is the parsed JSON dict on success and the error text on failure."""
    from mcp.client import Client

    async def go():
        async with Client(server) as c:
            return await c.call_tool(tool, args)

    r = anyio.run(go)
    text = r.content[0].text if r.content else ""
    if r.is_error:
        return True, text
    return False, json.loads(text)


def list_tools(server):
    from mcp.client import Client

    async def go():
        async with Client(server) as c:
            return (await c.list_tools()).tools

    return anyio.run(go)


def server_for(ws, **kw):
    from oxwave.mcp_server import make_server
    kw.setdefault("bin_dir", OXBIN)
    return make_server(ws, **kw)
