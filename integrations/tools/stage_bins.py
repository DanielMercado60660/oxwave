#!/usr/bin/env python3
"""Copy the OxWave program binaries from ../target/release
into target/oxbin/ (inside this root, so sandboxed tools and tests can read a stable
snapshot) and write target/oxbin/MANIFEST.json {program: sha256}.

  python3 tools/stage_bins.py [--src DIR]

Run after `cargo build --release --offline` of the main OxWave workspace, and again
after rebuilding. Tests read binaries from $OXWAVE_BIN (default target/oxbin).
"""
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OX = ROOT.parent


def main():
    src = Path(sys.argv[sys.argv.index("--src") + 1]) if "--src" in sys.argv else OX / "target" / "release"
    status = OX / "jobs_status.json"
    if status.exists():  # build records present: stage only programs that passed
        st = json.loads(status.read_text())
        passed = sorted(v["key"] for k, v in st.items() if v.get("kind") == "prog" and v.get("state") == "PASS")
    else:  # public repo: every program in programs/
        passed = sorted(p.name for p in (OX / "programs").iterdir() if (p / "Cargo.toml").is_file())
    dest = ROOT / "target" / "oxbin"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    man, missing = {}, []
    for p in passed:
        s = src / p
        if not (s.is_file() and os.access(s, os.X_OK)):
            missing.append(p)
            continue
        shutil.copy2(s, dest / p)
        man[p] = hashlib.sha256((dest / p).read_bytes()).hexdigest()
    (dest / "MANIFEST.json").write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
    print("staged %d binaries into %s" % (len(man), dest))
    if missing:
        print("PASS but not built (%d): %s" % (len(missing), " ".join(missing)))
        sys.exit(1)


if __name__ == "__main__":
    main()
