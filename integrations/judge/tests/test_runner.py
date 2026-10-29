# Tests for (py/python/oxwave/runner.py).
# Uses the staged OxWave binaries (target/oxbin); goldens from the original SU in judge/fixtures.
import json
import os

import numpy as np
import pytest

from conftest import FIX, OXBIN, UNAVAILABLE, bash, fixture_bytes, load_su, su_mismatch

import oxwave
from oxwave import (OutputLimitExceeded, ProgramError, ProgramNotFound, ProgramTimeout, Repeat,
                    RunResult, Traces, UnknownParameterError)

PIPES = json.loads((FIX / "pipes" / "PIPES.json").read_text())


@pytest.fixture(autouse=True)
def _bins(oxbin, monkeypatch):
    monkeypatch.setenv("OXWAVE_BIN", str(oxbin))
    oxwave.set_bin_dir(None)
    yield
    oxwave.set_bin_dir(None)


# ------------------------------------------------------------------ argv encoding --
@pytest.mark.parametrize("params,args,want", [
    ({"nt": 64, "ntr": 8}, (), ["suplane", "nt=64", "ntr=8"]),
    ({"scale": 2.5, "tpow": 1e-05, "x": 0.1}, (), ["suplane", "scale=2.5", "tpow=1e-05", "x=0.1"]),
    ({"f": [5, 10, 40.5, 50]}, (), ["suplane", "f=5,10,40.5,50"]),
    ({"f": np.array([1, 2], dtype=np.int16), "g": (np.float32(0.5), "a")}, (), ["suplane", "f=1,2", "g=0.5,a"]),
    ({"ref": Repeat(["1:1,2;4,2", "1:0,3;4,3"]), "v": None}, (), ["suplane", "ref=1:1,2;4,2", "ref=1:0,3;4,3"]),
    ({"r": Repeat([[1, 2], 3])}, (), ["suplane", "r=1,2", "r=3"]),
    ({"title": "My Data: 64ms", "n": np.int64(3)}, ("a.su", "b.su"), ["suplane", "a.su", "b.su", "title=My Data: 64ms", "n=3"]),
    (None, (), ["suplane"]),
])
def test_format_argv(params, args, want):
    assert oxwave.format_argv("suplane", params, args) == want


@pytest.mark.parametrize("params,args,exc", [
    ({"agc": True}, (), TypeError),
    ({"bad-key": 1}, (), ValueError),
    ({"1x": 1}, (), ValueError),
    ({"f": [[1, 2], [3]]}, (), TypeError),
    ({"x": float("nan")}, (), ValueError),
    ({"x": float("inf")}, (), ValueError),
    ({"s": "a\0b"}, (), ValueError),
    ({}, ("-x",), ValueError),
])
def test_format_argv_errors(params, args, exc):
    with pytest.raises(exc):
        oxwave.format_argv("suplane", params, args)


# ------------------------------------------------------------------ binaries, catalogue --
def test_bin_lookup(tmp_path, monkeypatch):
    assert oxwave.bin_dirs() == [OXBIN]
    assert oxwave.which("sugain") in (OXBIN / "sugain", (OXBIN / "sugain").resolve())
    assert oxwave.which("sugain").is_absolute()
    empty = tmp_path / "empty"
    empty.mkdir()
    oxwave.set_bin_dir(empty)
    assert oxwave.bin_dirs() == [empty]
    with pytest.raises(ProgramNotFound) as e:
        oxwave.which("sugain")
    assert "sugain" in str(e.value) and "OXWAVE_BIN" in str(e.value)
    with pytest.raises(ProgramNotFound):
        oxwave.which("ls")  # $PATH is never searched
    oxwave.set_bin_dir([empty, OXBIN])
    assert oxwave.which("sugain").name == "sugain"
    oxwave.set_bin_dir(None)
    monkeypatch.setenv("OXWAVE_BIN", os.pathsep.join([str(empty), str(OXBIN)]))
    assert oxwave.bin_dirs() == [empty, OXBIN]
    with pytest.raises(ValueError):
        oxwave.which("../oxbin/sugain")


def test_catalog():
    c = oxwave.catalog()
    assert len(c) == 400
    assert c["sugain"]["purpose"].startswith("Apply selectable gain")
    assert c["susynlv"]["status"] == "PASS"
    c["sugain"]["purpose"] = "changed"
    assert oxwave.program_info("sugain")["purpose"] != "changed"
    with pytest.raises(KeyError):
        oxwave.program_info("nosuchprog")


def test_available_programs():
    man = json.loads((OXBIN / "MANIFEST.json").read_text())
    av = oxwave.available_programs()
    assert av == sorted(set(man) & set(oxwave.catalog()))
    assert "sugain" in av and UNAVAILABLE not in av


# ------------------------------------------------------------------ run / pipe --
def test_run_generator_and_filter_bit_exact():
    tr = oxwave.run("suplane", nt=64, ntr=8)
    assert isinstance(tr, Traces) and (tr.ntr, tr.ns) == (8, 64)
    out = oxwave.run("sugain", tr, scale=2.5)
    want = bash(PIPES["gain"]).stdout
    assert out.to_bytes() == want
    assert su_mismatch(fixture_bytes("pipes/gain.su"), out.to_bytes()) is None


def test_pipe_matches_bash_and_original():
    out = oxwave.pipe([("suplane", {"nt": 128, "ntr": 12, "dt": 0.002}),
                       ("sufilter", {"f": [5, 10, 40, 50]}),
                       ("sugain", {"agc": 1, "wagc": 0.1})])
    want = bash(PIPES["filter_agc"]).stdout
    assert out.to_bytes() == want
    assert su_mismatch(fixture_bytes("pipes/filter_agc.su"), want) is None


def test_pipe_input_path_cwd_and_stage_forms(workdir):
    out = oxwave.pipe([("sushw", {"key": "offset", "a": -200, "b": 50}),
                       ("suwind", {"key": "offset", "min": -100, "max": 100}, ())],
                      "plane.su", cwd=workdir)
    want = bash(PIPES["wind_hw"], cwd=workdir).stdout
    assert out.to_bytes() == want and out.ntr == 5
    out2 = oxwave.pipe(["sugain", ("sugain",)], workdir / "plane.su")
    assert out2.to_bytes() == bash("sugain < plane.su | sugain", cwd=workdir).stdout


def test_input_bytes_traces_and_nmo(workdir):
    g = load_su(workdir / "gather.su")
    a = oxwave.run("sunmo", g, vnmo=[1500, 2500], tnmo=[0, 0.25])
    b = oxwave.run("sunmo", fixture_bytes("gather.su"), params={"vnmo": [1500, 2500], "tnmo": [0, 0.25]})
    want = bash(PIPES["nmo"], cwd=workdir).stdout
    assert a.to_bytes() == want and b.to_bytes() == want
    assert su_mismatch(fixture_bytes("pipes/nmo.su"), want) is None


def test_output_file_and_raw(workdir):
    r = oxwave.run("suplane", output=workdir / "o.su", nt=64, ntr=8)
    assert isinstance(r, RunResult) and r.returncode == 0 and r.stdout == b""
    assert (workdir / "o.su").read_bytes() == bash("suplane nt=64 ntr=8").stdout
    r = oxwave.run("suplane", output="o2.su", cwd=workdir, nt=16, ntr=2)
    assert (workdir / "o2.su").read_bytes() == bash("suplane nt=16 ntr=2").stdout
    r = oxwave.run("suplane", raw=True, nt=16, ntr=2)
    assert isinstance(r, RunResult) and r.stdout == bash("suplane nt=16 ntr=2").stdout
    assert r.argv == [["suplane", "nt=16", "ntr=2"]] and r.returncodes == [0]


def test_text_output_program():
    # `output` is run()'s own keyword, so the SU parameter output= goes through params=
    r = oxwave.run("sugethw", load_su(FIX / "plane.su"),
                   params={"key": "tracl,offset,sx,gx,cdp,d2,dt,ns", "output": "geom"})
    assert isinstance(r, RunResult)
    assert r.stdout == fixture_bytes("plane.gethw.txt")


def test_repeat_last_wins_like_getpar():
    tr = load_su(FIX / "plane.su")
    a = oxwave.run("sugain", tr, scale=Repeat([1, 2.5]))
    b = oxwave.run("sugain", tr, scale=2.5)
    assert a.to_bytes() == b.to_bytes()


def test_empty_su_output():
    out = oxwave.run("suwind", load_su(FIX / "plane.su"), key="offset", min=10000)
    assert isinstance(out, Traces) and out.ntr == 0


def test_program_error_and_check_false():
    with pytest.raises(ProgramError) as e:
        oxwave.run("sugain", b"")
    err = e.value
    assert err.program == "sugain" and err.returncode == 1 and err.stage == 0
    assert "fgettr.c: bad first header" in err.stderr and "fgettr.c: bad first header" in str(err)
    assert err.argv == ["sugain"]
    r = oxwave.run("sugain", b"", check=False)
    assert isinstance(r, RunResult) and r.returncodes == [1] and "bad first header" in r.stderr


def test_pipe_error_names_failing_stage():
    with pytest.raises(ProgramError) as e:
        oxwave.pipe([("suplane", {"nt": 32, "ntr": 4}), ("sufilter", {"f": "1,2,3", "amps": "1,0"}),
                     ("sugain", {"scale": 2})])
    assert e.value.stage == 1 and e.value.program == "sufilter" and e.value.returncode != 0


def test_sigpipe_upstream_is_not_failure():
    out = oxwave.pipe([("suplane", {"nt": 1000, "ntr": 200000}), ("suwind", {"count": 1})], timeout=60)
    assert out.ntr == 1


def test_timeout_kills(tmp_path):
    import time
    t0 = time.monotonic()
    with pytest.raises(ProgramTimeout):
        oxwave.pipe([("suplane", {"nt": 4000, "ntr": 60000}), ("sufilter", {"f": [5, 10, 40, 50]})],
                    output=tmp_path / "big.su", timeout=0.3)
    assert time.monotonic() - t0 < 5.0, "stages must be killed at the deadline"


def test_output_limit(tmp_path):
    import time
    t0 = time.monotonic()
    with pytest.raises(OutputLimitExceeded):  # ~1.5 GB producer: must be killed while streaming
        oxwave.run("suplane", output=tmp_path / "lim.su", max_output_bytes=10000, nt=64, ntr=5000000)
    assert not (tmp_path / "lim.su").exists()
    with pytest.raises(OutputLimitExceeded):
        oxwave.run("suplane", max_output_bytes=10000, nt=64, ntr=5000000)
    assert time.monotonic() - t0 < 10.0
    assert oxwave.run("suplane", max_output_bytes=10000, nt=64, ntr=8).ntr == 8


def test_validate_params():
    with pytest.raises(UnknownParameterError) as e:
        oxwave.run("sugain", load_su(FIX / "plane.su"), validate=True, scael=2)
    assert "scael" in str(e.value) and "sugain" in str(e.value)
    oxwave.validate_params("sugain", {"scale": 1, "par": "x"})
    oxwave.validate_params("suhilb", {"anything": 1})  # catalogue has no params: accept
    with pytest.raises(UnknownParameterError):
        oxwave.validate_params("sugain", {"bogus": 1})


def test_run_param_conflicts():
    with pytest.raises(TypeError):
        oxwave.run("suplane", params={"nt": 8}, nt=9)
    with pytest.raises(ProgramNotFound):
        oxwave.run(UNAVAILABLE)  # EXCLUDED from the port: no OxWave binary


@pytest.mark.parametrize("prog", ["sugain", "suplane", "sufilter", "surange", "sushw"])
def test_selfdoc_exact(prog):
    assert oxwave.selfdoc(prog) == (FIX / "selfdoc" / ("%s.txt" % prog)).read_bytes().decode("latin-1")


def test_help():
    assert oxwave.help("sugain") == oxwave.selfdoc("sugain")
    h = oxwave.help(UNAVAILABLE)
    assert UNAVAILABLE in h and "EXCLUDED" in h and "bendian" in h
    with pytest.raises(KeyError):
        oxwave.help("nosuchprog")


def test_never_uses_shell(tmp_path):
    marker = tmp_path / "pwned"
    r = oxwave.run("sugain", load_su(FIX / "plane.su"), cwd=tmp_path,
                   params={"tmpdir": "$(touch pwned)", "scale": "2;touch pwned"}, check=False, raw=True)
    assert isinstance(r, RunResult)
    assert not marker.exists()
    assert r.argv == [["sugain", "tmpdir=$(touch pwned)", "scale=2;touch pwned"]]
