#!/usr/bin/env python3
"""Integration check: SU's own demo flows, original SU vs OxWave.

  python3 tools/run_demo_flows.py [--only SUBSTR] [--jobs N] [--timeout S] [--out DIR]

For every demo directory under ../vendor/su-original/src/demos that holds non-X11 shell
scripts, the directory is copied three times into a scratch area and its scripts are run
in order (Make*/Model* first, then the rest alphabetically; X11 scripts, Clean.sh and RCS
files skipped) under three PATHs:

  orig   original SU only (vendor/su-build/bin)
  orig2  the same again (stability: files that differ between orig and orig2 are UNSTABLE)
  ox     every PASSed OxWave binary (target/oxbin, see tools/stage_bins.py) first, original
         SU behind it for programs not ported/passed yet (logged as fallbacks)

Every SU program is reached through a logging wrapper (so we know which programs each run
used); programs with a `seed` parameter get `seed=17` appended when the script gives none
(identical in all three runs, so unseeded noise does not make flows unstable).
Then every file created or modified by the scripts, plus each script's stdout, is compared
orig vs ox with tools/check_program.py's rules (oxharness.compare_stdout: SU streams
header-exact and samples within rtol 1e-5 of the trace max; text exact; EPS/PS compared
as text without %%CreationDate/%%Creator/%%Title lines; other bytes exact).

Verdict per directory: MATCH (all stable outputs agree, OxWave binaries actually used),
MATCH_ORIG_ONLY (agree but no PASSed OxWave program ran), DIFF, RC (exit codes differ),
ORIG_FAILS (original scripts fail here: missing data, compilers...), NO_OUTPUT, TIMEOUT.
Writes demo_flows/REPORT.md and demo_flows/results.json.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OX = ROOT.parent
DEMOS = OX / "vendor" / "su-original" / "src" / "demos"
SUBIN = OX / "vendor" / "su-build" / "bin"
CWPROOT = OX / "vendor" / "su-build"
OXBIN = Path(os.environ.get("OXWAVE_BIN") or ROOT / "target" / "oxbin")
sys.path.insert(0, str(OX / "tools"))
import oxharness as H  # noqa: E402

X11 = re.compile(r"\b(x(wigb|image|graph|movie|contour|picker|box|pswigb)|sux(wigb|image|graph|movie|contour|"
                 r"picker|edit|ximage)|suxwigb|xpicker)\b")
# Demo dirs excluded on licence grounds (LICENSE_AUDIT.md: susvdfilt is "All rights reserved").
# Suwaveform: SEG-copyright program, held (LICENSE_AUDIT.md).
EXCLUDE_DIRS = {"Filtering/Karhunen_Loeve", "Synthetic/Suwaveform"}
SKIP_NAMES = {"README", "Clean", "Clean.sh"}
VIEWERS = re.compile(r"^\s*(gv|ghostview|open|xdg-open|evince|display)\s", re.M)
FIRST = re.compile(r"^(Make|Model|model|Data|Generate|Synth|Create)")
PS_STRIP = re.compile(rb"^%%(CreationDate|Creator|Title|For):.*$", re.M)


def catalogue():
    d = json.loads((ROOT / "py" / "python" / "oxwave" / "_programs.json").read_text())
    return {p["name"]: p for p in d["programs"]}


def demo_dirs(only=None):
    out = {}
    for f in sorted(DEMOS.rglob("*")):
        if not f.is_file() or "RCS" in f.parts or f.name == "Clean.sh" or f.name.endswith(",v"):
            continue
        try:
            head = f.read_bytes()[:64]
        except OSError:
            continue
        first = head.split(b"\n")[0]
        if not (head.startswith(b"#!") and re.search(rb"/(sh|bash)\b", first)):
            continue  # csh/perl/... scripts are not run
        rel = f.parent.relative_to(DEMOS).as_posix()
        if rel in EXCLUDE_DIRS or (only and only not in rel):
            continue
        txt = f.read_text(errors="replace")
        if X11.search(txt) or f.name.startswith("X") or f.name in SKIP_NAMES or VIEWERS.search(txt):
            continue
        out.setdefault(rel, []).append(f.name)
    for rel in out:
        out[rel].sort(key=lambda n: (0 if FIRST.match(n) else 1, n))
    return out


def make_wrappers(dest, mode, cat):
    """One bash wrapper per SU program name: log, add seed=17 if applicable, exec real binary."""
    dest.mkdir(parents=True)
    seeded = {n for n, p in cat.items() if any(q["name"] == "seed" for q in p["params"])}
    passed = set(json.loads((OXBIN / "MANIFEST.json").read_text())) if mode == "ox" else set()
    for b in sorted(SUBIN.iterdir()):
        name = b.name
        if b.suffix == ".dSYM" or not os.access(b, os.X_OK) or b.is_dir():
            continue
        real, tag = (OXBIN / name, "ox") if name in passed else (b, "orig")
        seed = ('case " $* " in *" seed="*) ;; *) set -- "$@" seed=17 ;; esac\n' if name in seeded else "")
        (dest / name).write_text('#!/bin/bash\necho "%s %s" >> "$OXW_LOG"\n%sexec -a %s "%s" "$@"\n'
                                 % (name, tag, seed, name, real))
        (dest / name).chmod(0o755)


def snapshot(d):
    out = {}
    for p in d.rglob("*"):
        if p.is_file() and not p.is_symlink():
            out[p.relative_to(d).as_posix()] = hashlib.sha1(p.read_bytes()).hexdigest()
    return out


def run_dir(rel, scripts, cat, work, timeout, wrappers):
    res = {"dir": rel, "scripts": scripts, "runs": {}}
    src = DEMOS / rel
    for mode in ("orig", "orig2", "ox"):
        base = work / mode / rel.replace("/", "__")
        d = base / "d"
        shutil.copytree(src, d, symlinks=True)
        before = snapshot(d)
        log = base / "programs.log"
        wrap = wrappers["ox" if mode == "ox" else "orig"]
        env = H.base_env(str(wrap))
        env["OXW_LOG"] = str(log)
        env["PATH"] = str(wrap) + os.pathsep + "/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin"
        env["CWPROOT"] = str(CWPROOT)
        env.pop("DISPLAY", None)
        outs = []
        t0 = time.time()
        for s in scripts:
            try:
                shell = "bash" if b"bash" in (d / s).read_bytes()[:64].split(b"\n")[0] else "sh"
                p = subprocess.run([shell, s], cwd=d, env=env, stdin=subprocess.DEVNULL, capture_output=True,
                                   timeout=timeout, start_new_session=True)
                outs.append({"script": s, "rc": p.returncode, "stdout": p.stdout, "stderr": p.stderr[-2000:]})
            except subprocess.TimeoutExpired as e:
                outs.append({"script": s, "rc": "timeout", "stdout": e.stdout or b"", "stderr": b""})
        after = snapshot(d)
        changed = sorted(k for k, v in after.items() if before.get(k) != v)
        progs = {}
        if log.exists():
            for line in log.read_text().split("\n"):
                if line.strip():
                    n, tag = line.split()
                    progs.setdefault(n, tag)
        res["runs"][mode] = {"dir": d, "outs": outs, "changed": changed, "programs": progs,
                             "elapsed": round(time.time() - t0, 1)}
    return judge(res, cat)


def norm_ps(b):
    return PS_STRIP.sub(b"", b)


def compare_file(a, b, name):
    ga, gb = a.read_bytes(), b.read_bytes()
    if name.endswith((".eps", ".ps")) or ga.startswith(b"%!PS"):
        ga, gb = norm_ps(ga), norm_ps(gb)
        return None if ga == gb else "PostScript differs"
    return H.compare_stdout(ga, gb, {})


def judge(res, cat):
    o, o2, x = res["runs"]["orig"], res["runs"]["orig2"], res["runs"]["ox"]
    rcs_o = [r["rc"] for r in o["outs"]]
    rcs_x = [r["rc"] for r in x["outs"]]
    files = sorted(set(o["changed"]) | set(x["changed"]))
    unstable, diffs, matched = [], [], []
    for f in files:
        fo, f2, fx = o["dir"] / f, o2["dir"] / f, x["dir"] / f
        if not fo.exists() or not fx.exists():
            diffs.append((f, "only in %s" % ("orig" if fo.exists() else "ox")))
            continue
        if not f2.exists() or compare_file(fo, f2, f) is not None:
            unstable.append(f)
            continue
        why = compare_file(fo, fx, f)
        (diffs if why else matched).append((f, why) if why else f)
    for so, s2, sx in zip(o["outs"], o2["outs"], x["outs"]):
        if so["stdout"] != s2["stdout"]:
            unstable.append("<stdout of %s>" % so["script"])
        elif so["stdout"] != sx["stdout"]:
            why = H.compare_stdout(so["stdout"], sx["stdout"], {}) if so["stdout"] else "stdout differs"
            if why:
                diffs.append(("<stdout of %s>" % so["script"], why))
    used_ox = sorted(n for n, t in x["programs"].items() if t == "ox")
    fallback = sorted(n for n, t in x["programs"].items() if t == "orig")
    status = {n: cat.get(n, {}).get("status", "not in catalogue") for n in fallback}
    if any(r == "timeout" for r in rcs_o + rcs_x):
        verdict = "TIMEOUT"
    elif all(r != 0 for r in rcs_o) and not o["changed"]:
        verdict = "ORIG_FAILS"
    elif rcs_o != rcs_x:
        verdict = "RC"
    elif diffs:
        verdict = "DIFF"
    elif not matched:
        verdict = "NO_OUTPUT" if not unstable else "UNSTABLE_ONLY"
    elif not used_ox:
        verdict = "MATCH_ORIG_ONLY"
    else:
        verdict = "MATCH"
    return {"dir": res["dir"], "scripts": res["scripts"], "verdict": verdict, "rc_orig": rcs_o, "rc_ox": rcs_x,
            "matched": matched, "diffs": diffs, "unstable": unstable, "oxwave_programs": used_ox,
            "fallback_programs": status, "elapsed_s": [o["elapsed"], x["elapsed"]],
            "orig_stderr_tail": [r["stderr"].decode("latin-1")[-300:] for r in o["outs"] if r["rc"] not in (0,)],
            "ox_stderr_tail": [r["stderr"].decode("latin-1")[-300:] for r in x["outs"] if r["rc"] not in (0,)]}


def report(results, outdir, elapsed):
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "results.json").write_text(json.dumps(results, indent=1, default=str) + "\n")
    from collections import Counter
    c = Counter(r["verdict"] for r in results)
    blockers = Counter()
    for r in results:
        if r["verdict"] in ("MATCH", "MATCH_ORIG_ONLY", "DIFF", "RC"):
            for n, st in r["fallback_programs"].items():
                if st != "PASS":
                    blockers[n] += 1
    L = ["# OxWave demo-flow integration check", "",
         "Generated by `tools/run_demo_flows.py` on %s (%.0f s). OxWave binaries: `%s`."
         % (time.strftime("%Y-%m-%d %H:%M"), elapsed, OXBIN), "",
         "Each demo directory of `vendor/su-original/src/demos` with non-X11 shell scripts was run three",
         "times: original SU twice (stability) and OxWave-first (PASSed programs; original SU behind",
         "it for programs not ported yet). Outputs compared with check_program.py rules.", "",
         "## Summary", "", "| verdict | dirs |", "|---|---|"]
    L += ["| %s | %d |" % (k, v) for k, v in sorted(c.items())]
    L += ["", "Unported programs most often needed by otherwise-runnable flows: " +
          (", ".join("%s (%d)" % kv for kv in blockers.most_common(15)) or "none"), "",
          "## Per directory", "",
          "| dir | verdict | OxWave programs used | fallbacks to original (status) | outputs matched | diffs | unstable |",
          "|---|---|---|---|---|---|---|"]
    for r in results:
        fb = ", ".join("%s (%s)" % kv for kv in sorted(r["fallback_programs"].items()))
        df = "; ".join("%s: %s" % (f, (w or "")[:90]) for f, w in r["diffs"][:3])
        L.append("| %s | %s | %s | %s | %d | %s | %d |" % (
            r["dir"], r["verdict"], " ".join(r["oxwave_programs"]), fb, len(r["matched"]), df or "", len(r["unstable"])))
    (outdir / "REPORT.md").write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=300)
    ap.add_argument("--out", default=str(ROOT / "demo_flows"))
    ap.add_argument("--keep", action="store_true", help="keep the scratch directories")
    a = ap.parse_args()
    if not (OXBIN / "MANIFEST.json").exists():
        sys.exit("no staged OxWave binaries in %s: run tools/stage_bins.py" % OXBIN)
    cat = catalogue()
    dirs = demo_dirs(a.only)
    work = Path(tempfile.mkdtemp(prefix="oxw_demos_", dir=str(ROOT / "target")))
    wrappers = {}
    for m in ("orig", "ox"):
        wrappers[m] = work / ("bin_" + m)
        make_wrappers(wrappers[m], m, cat)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_dir, rel, s, cat, work, a.timeout, wrappers) for rel, s in dirs.items()]
        results = []
        for f in futs:
            r = f.result()
            print("%-55s %s" % (r["dir"], r["verdict"]), flush=True)
            results.append(r)
    report(results, Path(a.out), time.time() - t0)
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("report: %s/REPORT.md" % a.out)


if __name__ == "__main__":
    main()
