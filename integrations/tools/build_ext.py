#!/usr/bin/env python3
"""Build the `oxwave` wheel offline with maturin and extract the
compiled `oxwave/_native*.so` to target/native/ (never into the source tree). The judges'
conftest imports the package from py/python and loads that file as `oxwave._native`.

  .venv/bin/python tools/build_ext.py        # run from the integrations root

Wheels land in target/wheels/ (newest is used). Exit status non-zero on build failure.
"""
import glob
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUSTUP = "/opt/homebrew/opt/rustup/bin"


def main():
    env = dict(os.environ)
    env["PATH"] = RUSTUP + os.pathsep + env.get("PATH", "")
    env.setdefault("CARGO_TARGET_DIR", str(ROOT / "target"))
    out = ROOT / "target" / "wheels"
    out.mkdir(parents=True, exist_ok=True)
    py = ROOT / ".venv" / "bin" / "python"
    maturin = ROOT / ".venv" / "bin" / "maturin"
    cmd = [str(maturin), "build", "--release", "--offline", "-m", str(ROOT / "py" / "Cargo.toml"),
           "-i", str(py), "--out", str(out)]
    r = subprocess.run(cmd, cwd=ROOT, env=env)
    if r.returncode != 0:
        sys.exit("maturin build failed (exit %d)" % r.returncode)
    wheels = sorted(glob.glob(str(out / "oxwave-*.whl")), key=os.path.getmtime)
    if not wheels:
        sys.exit("no wheel produced in %s" % out)
    whl = wheels[-1]
    dest = ROOT / "target" / "native"
    dest.mkdir(parents=True, exist_ok=True)
    for old in dest.glob("_native*.so"):
        old.unlink()
    with zipfile.ZipFile(whl) as z:
        names = [n for n in z.namelist() if n.startswith("oxwave/_native") and n.endswith(".so")]
        if len(names) != 1:
            sys.exit("wheel %s: expected one oxwave/_native*.so, found %r" % (whl, names))
        with z.open(names[0]) as src, open(dest / Path(names[0]).name, "wb") as dst:
            shutil.copyfileobj(src, dst)
    print("built %s; native module -> %s" % (os.path.basename(whl), dest / Path(names[0]).name))


if __name__ == "__main__":
    main()
