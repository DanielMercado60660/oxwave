# Tests for (mcp_server/tools_flow.py): run_flow, plot_section.
import struct

import pytest

from conftest import FIX, UNAVAILABLE, bash, fixture_bytes, su_mismatch
from mcp_helpers import call, make_ws, server_for

import oxwave
from oxwave.mcp_server.config import Limits

PLANE_GAIN = [{"program": "suplane", "params": {"nt": 64, "ntr": 8}},
              {"program": "sugain", "params": {"scale": 2.5}}]


def png_size(path):
    b = open(path, "rb").read(24)
    assert b[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", b[16:24])


@pytest.fixture
def ws(tmp_path, oxbin):
    yield make_ws(tmp_path)
    oxwave.set_bin_dir(None)


def test_run_flow_generator_pipeline(ws):
    srv = server_for(ws)
    err, r = call(srv, "run_flow", {"stages": PLANE_GAIN, "output": "out/gain.su"})
    assert not err, r
    assert r["ok"] is True and r["output"] == "out/gain.su" and r["returncodes"] == [0, 0]
    assert r["argv"] == [["suplane", "nt=64", "ntr=8"], ["sugain", "scale=2.5"]]
    data = (ws / "out" / "gain.su").read_bytes()
    assert data == bash("suplane nt=64 ntr=8 | sugain scale=2.5").stdout
    assert su_mismatch(fixture_bytes("pipes/gain.su"), data) is None
    assert r["size_bytes"] == len(data) and r["elapsed_s"] >= 0
    assert r["summary"]["ntr"] == 8 and r["summary"]["ns"] == 64 and r["summary"]["path"] == "out/gain.su"


def test_run_flow_default_output_names(ws):
    srv = server_for(ws)
    _, a = call(srv, "run_flow", {"stages": PLANE_GAIN})
    _, b = call(srv, "run_flow", {"stages": PLANE_GAIN})
    assert (a["output"], b["output"]) == ("flows/flow-1.su", "flows/flow-2.su")
    _, c = call(srv, "run_flow", {"stages": [{"program": "sugethw", "params": {"key": "tracl", "output": "geom"}}],
                                  "input": "plane.su"})
    assert c["ok"] and c["output"] == "flows/flow-1.out" and c["summary"] is None
    assert (ws / "flows" / "flow-1.out").read_text().split() == [str(i) for i in range(1, 9)]


def test_run_flow_with_input_file(ws):
    stages = [{"program": "sushw", "params": {"key": "offset", "a": -200, "b": 50}},
              {"program": "suwind", "params": {"key": "offset", "min": -100, "max": 100}}]
    err, r = call(server_for(ws), "run_flow", {"stages": stages, "input": "plane.su", "output": "w.su"})
    assert not err and r["ok"]
    assert (ws / "w.su").read_bytes() == fixture_bytes("pipes/wind_hw.su")
    assert r["summary"]["ntr"] == 5


def test_run_flow_repeat_and_args(ws):
    err, r = call(server_for(ws), "run_flow", {"stages": [
        {"program": "suop2", "params": {"op": "diff"}, "args": ["plane.su", "plane.su"]},
        {"program": "sugain", "repeat": {"scale": [1, 3]}}], "output": "d.su"})
    assert not err and r["ok"], r
    assert r["argv"][1] == ["sugain", "scale=1", "scale=3"]
    assert (ws / "d.su").read_bytes() == bash("suop2 plane.su plane.su op=diff | sugain scale=1 scale=3", cwd=ws).stdout


def test_run_flow_failure_is_reported_not_raised(ws):
    err, r = call(server_for(ws), "run_flow", {"stages": [
        {"program": "suplane", "params": {"nt": 32, "ntr": 4}},
        {"program": "sufilter", "params": {"f": [1, 2, 3], "amps": [1, 0]}}], "output": "bad.su"})
    assert not err
    assert r["ok"] is False and r["returncodes"][1] == 1
    assert "number of f values must = number of amps values" in r["stderr"]
    assert r["summary"] is None


@pytest.mark.parametrize("args,msg", [
    ({"stages": [{"program": "sugain", "params": {"scael": 2}}], "input": "plane.su"}, "unknown parameter"),
    ({"stages": [{"program": "sugain"}], "input": "../plane.su"}, "outside the workspace"),
    ({"stages": [{"program": "sugain"}], "input": "/etc/hosts"}, "outside the workspace"),
    ({"stages": [{"program": "sugain"}], "input": "plane.su", "output": "../evil.su"}, "outside the workspace"),
    ({"stages": [{"program": "sugain"}], "input": "missing.su"}, "not found"),
    ({"stages": [{"program": UNAVAILABLE}]}, "not available"),
    ({"stages": [{"program": "sucommand"}]}, "not allowed"),
    ({"stages": []}, ""),
    ({"stages": PLANE_GAIN, "timeout_s": 100000}, ""),
    ({"stages": PLANE_GAIN, "timeout_s": 0}, ""),
])
def test_run_flow_rejections(ws, args, msg):
    err, text = call(server_for(ws), "run_flow", args)
    assert err and msg in text
    assert not (ws.parent / "evil.su").exists()


def test_run_flow_no_shell(ws):
    err, r = call(server_for(ws), "run_flow", {"stages": [
        {"program": "sugain", "params": {"tmpdir": "$(touch pwned)", "scale": "2;touch pwned2"}}],
        "input": "plane.su", "output": "x.su"})
    assert not err
    assert not (ws / "pwned").exists() and not (ws / "pwned2").exists()


def test_run_flow_timeout(ws):
    srv = server_for(ws)
    err, text = call(srv, "run_flow", {"stages": [
        {"program": "suplane", "params": {"nt": 4000, "ntr": 60000}},
        {"program": "sufilter", "params": {"f": [5, 10, 40, 50]}}], "output": "slow.su", "timeout_s": 0.3})
    assert err and "timed out" in text
    assert not (ws / "slow.su").exists()


def test_run_flow_output_limit(ws):
    srv = server_for(ws, limits=Limits(max_output_bytes=10000))
    err, text = call(srv, "run_flow", {"stages": [{"program": "suplane", "params": {"nt": 64, "ntr": 5000000}}],
                                       "output": "big.su"})
    assert err and "output limit" in text
    assert not (ws / "big.su").exists()
    err, r = call(srv, "run_flow", {"stages": [{"program": "suplane", "params": {"nt": 64, "ntr": 8}}],
                                    "output": "small.su"})
    assert not err and r["ok"]


def test_run_flow_stderr_truncated(ws):
    srv = server_for(ws, limits=Limits(max_stderr_chars=20))
    err, r = call(srv, "run_flow", {"stages": [{"program": "sugain"}], "output": "e.su"})
    assert not err and r["ok"] is False and len(r["stderr"]) <= 20
    assert r["stderr"].rstrip().endswith("bad first header")


@pytest.mark.parametrize("kind", ["image", "wiggle"])
def test_plot_section_tool(ws, kind):
    srv = server_for(ws)
    err, r = call(srv, "plot_section", {"path": "gather.su", "kind": kind, "title": "g"})
    assert not err, r
    assert r["png"] == "gather.png" and (r["width_px"], r["height_px"]) == (800, 600)
    assert (r["ntr"], r["ns"]) == (120, 250)
    assert png_size(ws / "gather.png") == (800, 600)
    err, r = call(srv, "plot_section", {"path": "plane.su", "output": "figs/p", "kind": kind})
    assert not err and r["png"] == "figs/p.png" and (ws / "figs" / "p.png").exists()


@pytest.mark.parametrize("args,msg", [
    ({"path": "notes.txt"}, "not a valid SU"),
    ({"path": "plane.su", "output": "../x.png"}, "outside the workspace"),
    ({"path": "../plane.su"}, "outside the workspace"),
    ({"path": "plane.su", "kind": "contour"}, ""),
])
def test_plot_section_rejections(ws, args, msg):
    err, text = call(server_for(ws), "plot_section", args)
    assert err and msg in text


def test_flow_then_inspect_then_plot(ws):
    srv = server_for(ws)
    _, r = call(srv, "run_flow", {"stages": [{"program": "sunmo", "params": {"vnmo": [1500, 2500], "tnmo": [0, 0.25]}}],
                                  "input": "gather.su", "output": "nmo.su"})
    assert r["ok"] and su_mismatch(fixture_bytes("pipes/nmo.su"), (ws / "nmo.su").read_bytes()) is None
    _, i = call(srv, "inspect_file", {"path": "nmo.su"})
    assert i["ntr"] == 120 and i["surange"] == (FIX / "gather.surange.txt").read_text()
    _, p = call(srv, "plot_section", {"path": "nmo.su"})
    assert p["png"] == "nmo.png"
