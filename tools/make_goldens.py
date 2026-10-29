#!/usr/bin/env python3
"""make_goldens.py <program> [--cases N]

Run every case of <program> (selfdoc, the inventory's test pipelines, and
goldens/<program>/cases.txt) with the ORIGINAL SU binaries (vendor/su-build,
override with OXWAVE_CWPROOT) and store the results in goldens/<program>/:
  <k>.out  final stdout      <k>.err  stderr      manifest.json  cases + exit codes
Each case runs twice; cases whose results differ (randomness, time) are
flagged unstable and skipped by check_program.py.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import oxharness as H  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('program')
    ap.add_argument('--cases', type=int, default=None, help='max inventory cases to use')
    a = ap.parse_args()
    prog = a.program
    if not os.path.exists(os.path.join(H.SUBIN, prog)):
        sys.exit('original binary %s not found in %s' % (prog, H.SUBIN))
    gdir = os.path.join(H.GOLDENS, prog)
    os.makedirs(gdir, exist_ok=True)
    cases = H.collect_cases(prog, a.cases)
    manifest = {'program': prog, 'cases': []}
    for k, (tags, cmd) in enumerate(cases):
        r1 = H.run_case(prog, cmd, H.SUBIN)
        r2 = H.run_case(prog, cmd, H.SUBIN)
        stable = (r1['rc'] == r2['rc'] and r1['out'] == r2['out']
                  and H.first_line(r1['err']) == H.first_line(r2['err']))
        with open(os.path.join(gdir, '%d.out' % k), 'wb') as f:
            f.write(r1['out'])
        with open(os.path.join(gdir, '%d.err' % k), 'wb') as f:
            f.write(r1['err'])
        kind = H.classify(r1['out'], tags, prog, cmd)
        manifest['cases'].append({'k': k, 'cmd': cmd, 'tags': tags, 'rc': r1['rc'],
                                  'stable': stable, 'kind': kind})
        flag = '' if stable else '  UNSTABLE (skipped by check)'
        if r1['rc'] == 'timeout':
            flag += '  TIMEOUT'
        print('%3d rc=%-4s %-5s %6d bytes  %s%s' % (k, r1['rc'], kind, len(r1['out']), cmd, flag))
    with open(os.path.join(gdir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1)
    print('wrote %d cases to %s' % (len(cases), gdir))


if __name__ == '__main__':
    main()
