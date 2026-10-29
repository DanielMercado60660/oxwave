#!/usr/bin/env python3
"""check_program.py <program> [--bin PATH] [-v]

Re-run every stable golden case of <program> with the OxWave binary
substituted for the original (only that program: a temp dir holding a
symlink named <program> is put in front of vendor/su-build/bin on PATH;
every other pipeline stage is still original SU) and compare with the
goldens written by make_goldens.py.  Default binary: target/release/<program>.

Comparison: SU trace streams trace by trace (headers bitwise identical,
samples within rtol 1e-5 of the trace's max |sample| plus a tiny floor);
text exactly (or whitespace-normalised for [ws] cases); raw float32 output
within tolerance (only [f32] cases or programs whose inventory stdout is
binary_floats, as last stage); other bytes exactly; selfdoc (stderr) exactly; exit code
always; first non-empty stderr line when the original failed.
Exit status 0 only if every case passes.
"""
import argparse
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import oxharness as H  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('program')
    ap.add_argument('--bin', default=None)
    ap.add_argument('-v', action='store_true', help='print stderr on failures')
    a = ap.parse_args()
    prog = a.program
    gdir = os.path.join(H.GOLDENS, prog)
    mpath = os.path.join(gdir, 'manifest.json')
    if not os.path.exists(mpath):
        sys.exit('no goldens for %s: run tools/make_goldens.py %s' % (prog, prog))
    binpath = os.path.abspath(a.bin or os.path.join(H.ROOT, 'target', 'release', prog))
    if not os.access(binpath, os.X_OK):
        sys.exit('OxWave binary not found: %s (cargo build --release -p %s)' % (binpath, prog))
    with open(mpath) as f:
        manifest = json.load(f)
    shim = tempfile.mkdtemp(prefix='oxw_bin_')
    os.symlink(binpath, os.path.join(shim, prog))
    path = shim + os.pathsep + H.SUBIN
    nfail = npass = nskip = 0
    for c in manifest['cases']:
        k, cmd, tags = c['k'], c['cmd'], c.get('tags', {})
        if not c.get('stable', True):
            nskip += 1
            print('SKIP %3d (unstable in original): %s' % (k, cmd))
            continue
        with open(os.path.join(gdir, '%d.out' % k), 'rb') as f:
            gout = f.read()
        with open(os.path.join(gdir, '%d.err' % k), 'rb') as f:
            gerr = f.read()
        r = H.run_case(prog, cmd, path)
        why = None
        if r['rc'] != c['rc']:
            why = 'exit code: golden %s actual %s' % (c['rc'], r['rc'])
        elif cmd == '@selfdoc':
            if r['err'] != gerr:
                why = 'selfdoc differs: ' + (H.compare_stdout(gerr, r['err'], {}) or 'bytes differ')
        else:
            why = H.compare_stdout(gout, r['out'], tags, c.get('kind'))
            if why is None and c['rc'] not in (0, 'timeout'):
                if H.first_line(gerr) != H.first_line(r['err']):
                    why = 'stderr first line: golden %r actual %r' % (H.first_line(gerr), H.first_line(r['err']))
            if why is None and tags.get('stderr') and gerr != r['err']:
                why = 'stderr differs: ' + (H.compare_stdout(gerr, r['err'], {'ws': False}) or '')
        if why:
            nfail += 1
            print('FAIL %3d: %s\n         %s' % (k, cmd, why))
            if a.v:
                print('         actual stderr: %r' % r['err'][:500])
        else:
            npass += 1
            print('PASS %3d: %s' % (k, cmd))
    os.unlink(os.path.join(shim, prog))
    os.rmdir(shim)
    print('%s: %d passed, %d failed, %d skipped' % (prog, npass, nfail, nskip))
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
