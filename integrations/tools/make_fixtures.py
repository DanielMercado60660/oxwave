#!/usr/bin/env python3
"""Generate judge/fixtures/ with the ORIGINAL SU binaries
(../vendor/su-build/bin). The fixtures are committed; tests never call the
original binaries themselves.

  python3 tools/make_fixtures.py

Fixtures (all produced by original SU in a scratch dir, then copied):
  plane.su, spike.su, gather.su, one.su       SU streams (XDR big-endian)
  plane.gethw.txt, gather.gethw.txt           sugethw ... output=geom
  plane.sushw1.su, gather.sushw2.su           sushw goldens
  *.surange.txt                               surange stdout
  segy/header, segy/binary                    segyhdrs < plane.su (ASCII text header, binary header)
  segy/plane.sgy, plane5.sgy, wild.sgy        segywrite (format 1 IBM / format 5 IEEE / edge values)
  segy/<x>.sgy.su, .hfile, .bfile             segyread tape=<x>.sgy hfile= bfile=
  (endian= is ignored by segyread/segywrite in this SUXDR build: endian=0 gives identical bytes)
  selfdoc/<prog>.txt                          selfdoc (stderr of the bare program on a pty)
  pipes/<name>.su + pipes/PIPES.json          original outputs of the pipelines run()/pipe() judges use
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OX = ROOT.parent
sys.path.insert(0, str(OX / "tools"))
import oxharness as H  # noqa: E402  (read-only use of the harness)

FIX = ROOT / "judge" / "fixtures"
SUBIN = OX / "vendor" / "su-build" / "bin"

GEN = {
    "plane.su": "suplane nt=64 ntr=8 | sushw key=offset,sx,gx,d2,cdp a=-200,0,100,0.025,1 b=50,10,20,0,1",
    "spike.su": "suspike nt=100 ntr=6 dt=0.002 | sushw key=f1,scalco,gelev a=0.1,-10,-35 b=0,0,7",
    "gather.su": ("suplane nt=250 ntr=120 dt=0.002 | sushw key=fldr,tracf,offset,sx,gx a=1,1,100,1000,1100 "
                  "b=0,1,25,0,25 c=1,0,0,50,50 j=24,24,24,24,24 | suchw key1=cdp key2=gx key3=sx b=1 c=1 d=2"),
    "one.su": "suspike nt=1 ntr=1 nspk=1 ix1=1 it1=1",
}
DERIVED = [
    ("plane.gethw.txt", "sugethw key=tracl,offset,sx,gx,cdp,d2,dt,ns output=geom < plane.su"),
    ("gather.gethw.txt", "sugethw key=fldr,tracf,offset,sx,gx,cdp,trid output=geom < gather.su"),
    ("plane.sushw1.su", "sushw key=cdp,fldr,f2 a=10,1,0.5 b=2,0,0.25 c=0,1,0 j=3,3,3 < plane.su"),
    ("gather.sushw2.su", "sushw key=ep,offset,d1 a=5,-100,0.004 b=1,10,0 c=100,0,0 d=2,0,0 j=10,7,0 < gather.su"),
    ("plane.surange.txt", "surange < plane.su"),
    ("spike.surange.txt", "surange < spike.su"),
    ("gather.surange.txt", "surange < gather.su"),
    ("one.surange.txt", "surange < one.su"),
]
SEGY = [
    "segyhdrs < plane.su > /dev/null",  # writes ./header (ASCII) and ./binary
    "segywrite tape=plane.sgy < plane.su",
    "segywrite tape=plane5.sgy format=5 < plane.su",
    "segyread tape=plane.sgy hfile=plane.sgy.hfile bfile=plane.sgy.bfile > plane.sgy.su",
    "segyread tape=plane5.sgy hfile=plane5.sgy.hfile bfile=plane5.sgy.bfile > plane5.sgy.su",
    "suaddhead ns=64 < wild.f32 | sushw key=dt a=4000 > wild.su",
    "segywrite tape=wild.sgy hfile=header bfile=binary < wild.su",
    "segyread tape=wild.sgy hfile=wild.sgy.hfile bfile=wild.sgy.bfile > wild.sgy.su",
]
SELFDOC = ["sugain", "suplane", "sufilter", "surange", "sushw"]
# Pipelines the runner judges run through OxWave (bit-exact vs bash) and compare to these originals.
PIPES = {
    "gain": "suplane nt=64 ntr=8 | sugain scale=2.5",
    "filter_agc": "suplane nt=128 ntr=12 dt=0.002 | sufilter f=5,10,40,50 | sugain agc=1 wagc=0.1",
    "wind_hw": "sushw key=offset a=-200 b=50 < plane.su | suwind key=offset min=-100 max=100",
    "nmo": "sunmo vnmo=1500,2500 tnmo=0,0.25 < gather.su",
}


def wild_floats():
    """6 traces x 64 native float32 samples covering IBM/IEEE edge cases: zeros (+/-),
    denormals, tiny/huge normals, powers of 16 and 2, random mantissas, inf and NaN."""
    import random
    import struct
    rnd = random.Random(20261028)
    words = [0x00000000, 0x80000000, 0x00000001, 0x007fffff, 0x80000001, 0x00800000, 0x7f7fffff,
             0xff7fffff, 0x7f800000, 0xff800000, 0x7fc00000, 0x7fa00001, 0x3f800000, 0xbf800000,
             0x41800000, 0x43800000, 0x3d800000, 0x33800000, 0x0c800000, 0x01000000]
    while len(words) < 6 * 64:
        e = rnd.randrange(1, 255)
        words.append((rnd.getrandbits(1) << 31) | (e << 23) | rnd.getrandbits(23))
    return b"".join(struct.pack("<I", w) for w in words)


def sh(cmd, cwd):
    env = H.base_env(str(SUBIN))
    r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], cwd=cwd, env=env,
                       stdin=subprocess.DEVNULL, capture_output=True)
    if r.returncode != 0:
        sys.exit("FAILED (%d): %s\n%s" % (r.returncode, cmd, r.stderr.decode("latin-1")))
    return r.stdout


def main():
    if FIX.exists():
        shutil.rmtree(FIX)
    (FIX / "segy").mkdir(parents=True)
    (FIX / "selfdoc").mkdir()
    (FIX / "pipes").mkdir()
    work = Path(tempfile.mkdtemp(prefix="oxw_fix_"))
    for name, cmd in GEN.items():
        (work / name).write_bytes(sh(cmd, work))
        shutil.copy(work / name, FIX / name)
    for name, cmd in DERIVED:
        (FIX / name).write_bytes(sh(cmd, work))
    (work / "wild.f32").write_bytes(wild_floats())
    for cmd in SEGY:
        sh(cmd, work)
    for f in ["header", "binary", "plane.sgy", "plane5.sgy", "wild.su", "wild.sgy"] + \
            ["%s.sgy.%s" % (b, e) for b in ("plane", "plane5", "wild") for e in ("su", "hfile", "bfile")]:
        shutil.copy(work / f, FIX / "segy" / f)
    env = H.base_env(str(SUBIN))
    for p in SELFDOC:
        r = H.run_selfdoc(p, env, str(work))
        assert r["rc"] == 1 and r["err"], (p, r["rc"])
        (FIX / "selfdoc" / ("%s.txt" % p)).write_bytes(r["err"])
    for name, cmd in PIPES.items():
        (FIX / "pipes" / ("%s.su" % name)).write_bytes(sh(cmd, work))
    (FIX / "pipes" / "PIPES.json").write_text(json.dumps(PIPES, indent=1) + "\n")
    (FIX / "GENERATORS.json").write_text(json.dumps({"gen": GEN, "derived": dict(DERIVED), "segy": SEGY,
                                                     "selfdoc": SELFDOC}, indent=1) + "\n")
    shutil.rmtree(work)
    n = sum(1 for p in FIX.rglob("*") if p.is_file())
    print("wrote %d fixture files to %s" % (n, FIX))


if __name__ == "__main__":
    main()
