# Shared setup for the OxWave integrations judges.
"""Imports `oxwave` from py/python (sources), loads the native module built by
tools/build_ext.py from target/native/ if present, points OXWAVE_BIN at the staged OxWave
binaries (target/oxbin, see tools/stage_bins.py) and keeps matplotlib headless and inside
this root. Fixtures in judge/fixtures were produced by the ORIGINAL SU binaries
(tools/make_fixtures.py); tests never run original SU."""
import glob
import importlib.machinery
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "judge" / "fixtures"
OX = ROOT.parent

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "target" / "mpl"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
OXBIN = Path(os.environ.get("OXWAVE_BIN") or (ROOT / "target" / "oxbin"))
os.environ["OXWAVE_BIN"] = str(OXBIN)

sys.path.insert(0, str(ROOT / "py" / "python"))
sys.path.insert(0, str(OX / "tools"))  # oxharness (read-only): SU comparison with tolerance


def _load_native():
    sos = sorted(glob.glob(str(ROOT / "target" / "native" / "_native*.so")))
    if not sos:
        return None
    loader = importlib.machinery.ExtensionFileLoader("oxwave._native", sos[-1])
    spec = importlib.util.spec_from_file_location("oxwave._native", sos[-1], loader=loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    sys.modules["oxwave._native"] = mod
    return mod


NATIVE = _load_native()
import oxwave  # noqa: E402

if NATIVE is not None:
    oxwave._native = NATIVE


def fixture_bytes(name):
    return (FIX / name).read_bytes()


def load_su(path):
    """Fixture traces through the core codec (not read_su): `path` is a fixture
    name or a Path."""
    from oxwave.traces import Traces
    p = path if isinstance(path, Path) else FIX / path
    return Traces.from_bytes(p.read_bytes())


# A catalogue program that is never built (EXCLUDED from the port), used wherever a judge
# needs an "unavailable" program; BLOCKED programs may PASS later and be re-staged.
UNAVAILABLE = "su2voxet"


def oxbin_env():
    """Environment whose PATH holds ONLY the staged OxWave binaries plus the system dirs."""
    env = dict(os.environ)
    env["PATH"] = str(OXBIN) + os.pathsep + "/usr/bin:/bin"
    env["LC_ALL"] = "C"
    return env


def bash(cmd, cwd=None, input=None):
    """Run `cmd` with bash -o pipefail using only OxWave binaries; returns CompletedProcess."""
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], cwd=cwd, env=oxbin_env(),
                          input=input, stdin=None if input is not None else subprocess.DEVNULL,
                          capture_output=True, timeout=300)


def su_mismatch(golden, actual, rtol=1e-5):
    """None if `actual` matches `golden` with check_program's rules (headers bit-exact,
    samples within rtol of the trace max), else a message."""
    import oxharness as H
    return H.compare_stdout(golden, actual, {"rtol": rtol}, "su")


@pytest.fixture
def native():
    if NATIVE is None:
        pytest.fail("oxwave._native not built: run .venv/bin/python tools/build_ext.py")
    return NATIVE


@pytest.fixture(scope="session")
def oxbin():
    if not (OXBIN / "suplane").exists():
        pytest.fail("OxWave binaries not staged in %s: run tools/stage_bins.py first" % OXBIN)
    return OXBIN


@pytest.fixture
def workdir(tmp_path):
    """A temp dir holding copies of the SU fixtures (plane.su, spike.su, gather.su, one.su)."""
    for n in ("plane.su", "spike.su", "gather.su", "one.su"):
        shutil.copy(FIX / n, tmp_path / n)
    return tmp_path
