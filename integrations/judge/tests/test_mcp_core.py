# Tests for (mcp_server/workspace.py, mcp_server/tools_info.py).
import os

import pytest

from conftest import FIX, OXBIN, UNAVAILABLE, load_su
from mcp_helpers import call, list_tools, make_ws, server_for

import oxwave
from oxwave.mcp_server.config import TOOL_NAMES, Limits, ToolError
from oxwave.mcp_server.workspace import Workspace


@pytest.fixture
def ws(tmp_path, oxbin):
    oxwave.set_bin_dir(oxbin)
    yield Workspace(make_ws(tmp_path), bin_dir=oxbin)
    oxwave.set_bin_dir(None)


# ------------------------------------------------------------------ sandbox --
def test_workspace_basics(tmp_path, ws):
    assert ws.root == (tmp_path / "ws").resolve() and ws.limits == Limits()
    with pytest.raises(ValueError):
        Workspace(tmp_path / "missing")
    p = ws.resolve("plane.su", must_exist=True)
    assert p == ws.root / "plane.su" and ws.relative(p) == "plane.su"
    assert ws.resolve(str(ws.root / "plane.su")) == ws.root / "plane.su"
    assert ws.resolve("sub/dir/../x.su") == ws.root / "sub" / "x.su"
    q = ws.resolve("a/b/c.su", for_write=True)
    assert q.parent.is_dir() and ws.relative(q) == "a/b/c.su"


@pytest.mark.parametrize("bad", ["../outside.su", "/etc/passwd", "a/../../x", "sub/../../../etc/hosts"])
def test_workspace_rejects_outside(ws, bad):
    with pytest.raises(ToolError, match="outside the workspace"):
        ws.resolve(bad)


def test_workspace_rejects_symlink_escape_and_junk(tmp_path, ws):
    secret = tmp_path / "secret.su"
    secret.write_bytes(b"x")
    os.symlink(secret, ws.root / "link.su")
    os.symlink(tmp_path, ws.root / "linkdir")
    for p in ("link.su", "linkdir/secret.su"):
        with pytest.raises(ToolError, match="outside the workspace"):
            ws.resolve(p)
    for p in ("", "~/x.su", "a\0b"):
        with pytest.raises(ToolError):
            ws.resolve(p)
    with pytest.raises(ToolError, match="not found"):
        ws.resolve("nope.su", must_exist=True)
    (ws.root / "d").mkdir()
    with pytest.raises(ToolError):
        ws.resolve("d", for_write=True)


# ------------------------------------------------------------------ flow validation --
def test_check_program(ws):
    ws.check_program("sugain")
    with pytest.raises(ToolError, match="unknown program"):
        ws.check_program("rm")
    with pytest.raises(ToolError, match="not allowed"):
        ws.check_program("sucommand")
    with pytest.raises(ToolError, match="not available") as e:
        ws.check_program(UNAVAILABLE)
    assert "EXCLUDED" in str(e.value)


def test_stage_argv_encoding(ws):
    argv = ws.stage_argv({"program": "sufilter", "params": {"f": [5, 10, 40, 50], "amps": "0,1,1,0"}})
    assert argv == ["sufilter", "f=5,10,40,50", "amps=0,1,1,0"]
    argv = ws.stage_argv({"program": "suop2", "params": {"op": "diff"}, "args": ["plane.su", "plane.su"]})
    assert argv == ["suop2", "plane.su", "plane.su", "op=diff"]
    argv = ws.stage_argv({"program": "sugain", "params": {"scale": 2.5}, "repeat": {"tpow": [1, 2]}})
    assert argv == ["sugain", "scale=2.5", "tpow=1", "tpow=2"]
    assert ws.stage_argv({"program": "suhilb", "params": {"whatever": 3}}) == ["suhilb", "whatever=3"]
    assert ws.stage_argv({"program": "sugain", "params": {"par": "gain.par"}}) == ["sugain", "par=gain.par"]


@pytest.mark.parametrize("stage,msg", [
    ({"program": "sugain", "params": {"scael": 2}}, "unknown parameter"),
    ({"program": "sugain", "params": {"agc": True}}, ""),
    ({"program": "sugain", "params": {"scale": None}}, ""),
    ({"program": "sugain", "params": {"bad key": 1}}, ""),
    ({"program": "sufilter", "params": {"f": [[1, 2], 3]}}, ""),
    ({"program": "sugain", "params": {"scale": float("nan")}}, ""),
    ({"program": "sugain", "params": {"tmpdir": "a\nb"}}, ""),
    ({"program": "sugain", "params": {"tmpdir": "/tmp"}}, "outside the workspace"),
    ({"program": "sugain", "params": {"par": "../x.par"}}, "outside the workspace"),
    ({"program": "sumax", "params": {"outpar": "/dev/tty"}}, "outside the workspace"),
    ({"program": "segyread", "params": {"tape": "x.sgy", "hfile": "h;rm -rf ."}}, "unsafe value"),
    ({"program": "suop2", "args": ["../plane.su", "plane.su"]}, "outside the workspace"),
    ({"program": "suop2", "args": ["-x", "plane.su"]}, ""),
    ({"program": "sugain", "params": {"scale": 1}, "repeat": {"scale": [2]}}, ""),
    ({"program": "sugain", "flags": ["-v"]}, ""),
    ({"program": "sucommand"}, "not allowed"),
])
def test_stage_argv_rejects(ws, stage, msg):
    with pytest.raises(ToolError, match=msg or None):
        ws.stage_argv(stage)


def test_flow_argvs_limits(tmp_path, oxbin):
    w = Workspace(make_ws(tmp_path), bin_dir=oxbin, limits=Limits(max_stages=2))
    oxwave.set_bin_dir(oxbin)
    try:
        assert w.flow_argvs([{"program": "suplane"}, {"program": "sugain"}]) == [["suplane"], ["sugain"]]
        with pytest.raises(ToolError):
            w.flow_argvs([])
        with pytest.raises(ToolError):
            w.flow_argvs([{"program": "suplane"}] + [{"program": "sugain"}] * 2)
    finally:
        oxwave.set_bin_dir(None)


# ------------------------------------------------------------------ tools through MCP --
def test_tool_list_and_schemas(tmp_path):
    tools = {t.name: t for t in list_tools(server_for(make_ws(tmp_path)))}
    assert tuple(tools) == TOOL_NAMES
    rf = tools["run_flow"].input_schema
    assert rf["required"] == ["stages"]
    assert all(t.description for t in tools.values())


def test_list_programs(tmp_path):
    srv = server_for(make_ws(tmp_path))
    err, r = call(srv, "list_programs", {})
    assert not err
    names = [p["name"] for p in r["programs"]]
    assert names == sorted(names) and r["count"] == len(names)
    assert "sugain" in names and UNAVAILABLE not in names and "sucommand" not in names
    g = next(p for p in r["programs"] if p["name"] == "sugain")
    assert set(g) == {"name", "purpose", "family", "stdin", "stdout", "status", "available"}
    assert g["available"] is True and g["status"] == "PASS" and g["purpose"]
    err, r = call(srv, "list_programs", {"available_only": False})
    allp = {p["name"]: p for p in r["programs"]}
    assert len(allp) == 400 and allp[UNAVAILABLE]["available"] is False and allp[UNAVAILABLE]["status"] == "EXCLUDED"
    assert allp["susynlv"]["status"] == "PASS"  # catalogue snapshot (refreshed 2026-10-29)
    assert allp["sucommand"]["available"] is False
    err, r = call(srv, "list_programs", {"filter": "NMO"})
    assert "sunmo" in [p["name"] for p in r["programs"]]
    assert all("nmo" in (p["name"] + p["purpose"] + p["family"]).lower() for p in r["programs"])


def test_describe_program(tmp_path):
    srv = server_for(make_ws(tmp_path))
    err, r = call(srv, "describe_program", {"name": "sugain"})
    assert not err
    assert r["selfdoc"] == (FIX / "selfdoc" / "sugain.txt").read_bytes().decode("latin-1")
    assert r["params_known"] is True and "agc" in [p["name"] for p in r["params"]]
    assert r["stdin"] == "su_traces" and r["available"] is True
    err, r = call(srv, "describe_program", {"name": UNAVAILABLE})
    assert not err and r["selfdoc"] is None and r["status"] == "EXCLUDED" and r["available"] is False
    err, r = call(srv, "describe_program", {"name": "suhilb"})
    assert r["params_known"] is False
    err, msg = call(srv, "describe_program", {"name": "nosuch"})
    assert err and "unknown program" in msg


def test_inspect_file(tmp_path):
    srv = server_for(make_ws(tmp_path))
    err, r = call(srv, "inspect_file", {"path": "plane.su"})
    assert not err
    assert (r["path"], r["format"], r["ntr"], r["ns"], r["truncated"]) == ("plane.su", "su", 8, 64, False)
    assert r["size_bytes"] == 3968 and abs(r["dt"] - 0.004) < 1e-12
    assert r["ranges"]["offset"] == [-200, 150] and list(r["ranges"])[:3] == ["tracl", "tracr", "cdp"]
    assert r["surange"] == (FIX / "plane.surange.txt").read_text()
    assert r["first_header"]["offset"] == -200 and "fldr" not in r["first_header"]
    tr = load_su("plane.su")
    assert r["stats"]["max"] == pytest.approx(float(tr.data.max()))
    assert r["stats"]["rms"] == pytest.approx(float((tr.data.astype(float) ** 2).mean() ** 0.5))
    assert r["stats"]["nan_count"] == 0
    err, r = call(srv, "inspect_file", {"path": "gather.su", "max_traces": 10})
    assert r["ntr"] == 10 and r["truncated"] is True


@pytest.mark.parametrize("path,msg", [("../plane.su", "outside the workspace"),
                                      ("/etc/hosts", "outside the workspace"),
                                      ("missing.su", "not found"),
                                      ("notes.txt", "not a valid SU")])
def test_inspect_file_rejects(tmp_path, path, msg):
    err, text = call(server_for(make_ws(tmp_path)), "inspect_file", {"path": path})
    assert err and msg in text


def test_list_files(tmp_path):
    ws = make_ws(tmp_path)
    (ws / "sub").mkdir()
    (ws / "sub" / "x.png").write_bytes(b"\x89PNG")
    (ws / ".hidden").write_text("x")
    err, r = call(server_for(ws), "list_files", {})
    assert not err and r["root"] == "."
    got = {f["path"]: f for f in r["files"]}
    assert list(got) == sorted(got)
    assert got["plane.su"] == {"path": "plane.su", "size": 3968, "kind": "su"}
    assert got["sub/x.png"]["kind"] == "png" and got["notes.txt"]["kind"] == "text"
    assert ".hidden" not in got
    err, r = call(server_for(ws), "list_files", {"subdir": "sub"})
    assert [f["path"] for f in r["files"]] == ["sub/x.png"]
    err, text = call(server_for(ws), "list_files", {"subdir": ".."})
    assert err and "outside the workspace" in text


def test_bin_dir_is_respected(tmp_path):
    empty = tmp_path / "nobin"
    empty.mkdir()
    err, r = call(server_for(make_ws(tmp_path), bin_dir=empty), "list_programs", {})
    assert not err and r["count"] == 0
    oxwave.set_bin_dir(None)
    assert OXBIN.exists()
