"""Shared code for make_goldens.py and check_program.py (Python 3 stdlib only).

A *case* is a shell pipeline (run with `bash -o pipefail -c`) whose final
stdout, stderr and exit code are recorded.  Inputs come from SU's own
generators inside the pipeline.  The special case `@selfdoc` runs the
program with no arguments on a pseudo-terminal (stdin+stdout) with
PAGER=cat, so the selfdoc arrives on stderr.

Case tags (in goldens/<prog>/cases.txt, before the command):
  [ws]        compare text after collapsing runs of blanks/tabs and trailing blanks
  [exact]     compare stdout bytes exactly (no float tolerance)
  [rtol=X]    relative tolerance for floats (default 1e-5, relative to the
              max |sample| of the golden trace or array)
  [atol=X]    absolute floor (default 1e-30)
  [stderr]    also compare the whole stderr exactly
  [f32]       force stdout to be compared as raw native float32 samples
"""
import json
import math
import os
import pty
import re
import select
import shutil
import struct
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CWPROOT = os.environ.get('OXWAVE_CWPROOT', os.path.join(ROOT, 'vendor', 'su-build'))
SUBIN = os.path.join(CWPROOT, 'bin')
GOLDENS = os.environ.get('OXWAVE_GOLDENS', os.path.join(ROOT, 'goldens'))
INVENTORY = os.path.join(ROOT, 'inventory')
TIMEOUT = float(os.environ.get('OXWAVE_CASE_TIMEOUT', '120'))


def inventory_cases(prog):
    """All test_cases commands for prog across every inventory/*.json file."""
    cmds = []
    for fn in sorted(os.listdir(INVENTORY)):
        if not fn.endswith('.json'):
            continue
        try:
            with open(os.path.join(INVENTORY, fn)) as f:
                d = json.load(f)
        except (OSError, ValueError):
            continue
        for p in d.get('programs', []):
            if p.get('name') == prog:
                for tc in p.get('test_cases', []):
                    cmd = (tc.get('command') or '').strip()
                    if cmd:
                        cmds.append(cmd)
    return cmds


TAG_RE = re.compile(r'^\s*\[([a-z0-9_.=+-]+)\]\s*')


def parse_tags(line):
    tags = {}
    while True:
        m = TAG_RE.match(line)
        if not m:
            break
        t = m.group(1)
        if '=' in t:
            k, v = t.split('=', 1)
            tags[k] = v
        else:
            tags[t] = True
        line = line[m.end():]
    return tags, line.strip()


def invokes(prog, cmd):
    """True if the pipeline runs `prog` by bare name (so the PATH shim applies)."""
    pat = r'(^|[\s|;&(`])' + re.escape(prog) + r'($|[\s|;&)<>`])'
    return re.search(pat, cmd) is not None


def last_stage_is(prog, cmd):
    """Crude: is `prog` the command of the last pipeline stage?"""
    last = re.split(r'\|\||&&|;|\|', cmd)[-1].split()
    words = [w for w in last if '=' not in w or w.startswith(('<', '>'))]
    return bool(words) and words[0] == prog


def inventory_stdout(prog):
    for fn in sorted(os.listdir(INVENTORY)):
        if fn.endswith('.json'):
            try:
                with open(os.path.join(INVENTORY, fn)) as f:
                    d = json.load(f)
            except (OSError, ValueError):
                continue
            for p in d.get('programs', []):
                if p.get('name') == prog:
                    return p.get('stdout')
    return None


def collect_cases(prog, max_cases=None):
    """[(tags, command)] from the inventory, then goldens/<prog>/cases.txt."""
    cases = [({}, '@selfdoc')]
    inv = [({}, cmd) for cmd in inventory_cases(prog)]
    if max_cases is not None:
        inv = inv[:max_cases]
    cases += inv
    extra = os.path.join(GOLDENS, prog, 'cases.txt')
    if os.path.exists(extra):
        with open(extra) as f:
            for line in f:
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                cases.append(parse_tags(line.rstrip('\n')))
    # Excluded cases (e.g. cases that exercise undefined behaviour in SU's C), read from the
    # repository's goldens folder, not the (possibly private) golden output folder.
    excl_path = os.path.join(ROOT, 'goldens', prog, 'EXCLUDE.txt')
    excluded = set()
    if os.path.exists(excl_path):
        with open(excl_path) as f:
            excluded = {l.rstrip('\n') for l in f if l.strip() and not l.lstrip().startswith('#')}
    seen, out = set(), []
    for t, c in cases:
        if c in excluded:
            sys.stderr.write('EXCLUDED case: %s\n' % c)
            continue
        if c != '@selfdoc' and not invokes(prog, c):
            sys.stderr.write('IGNORED case (does not invoke %s by bare name): %s\n' % (prog, c))
            continue
        if c not in seen:
            seen.add(c)
            out.append((t, c))
    return out


def base_env(path_prefix):
    env = dict(os.environ)
    env['CWPROOT'] = CWPROOT
    env['PATH'] = path_prefix + os.pathsep + env.get('PATH', '/usr/bin:/bin')
    env['PAGER'] = 'cat'
    env['LC_ALL'] = 'C'
    return env


def run_case(prog, cmd, path_prefix):
    """Run one case in a fresh temp dir; returns dict(rc, out, err)."""
    env = base_env(path_prefix)
    work = tempfile.mkdtemp(prefix='oxw_case_')
    try:
        if cmd == '@selfdoc':
            return run_selfdoc(prog, env, work)
        try:
            p = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], cwd=work, env=env,
                               stdin=subprocess.DEVNULL, capture_output=True,
                               timeout=TIMEOUT, start_new_session=True)
            return {'rc': p.returncode, 'out': p.stdout, 'err': p.stderr}
        except subprocess.TimeoutExpired as e:
            return {'rc': 'timeout', 'out': e.stdout or b'', 'err': e.stderr or b''}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run_selfdoc(prog, env, work):
    """Program with no args, stdin+stdout on a pty, stderr captured."""
    master, slave = pty.openpty()
    errr, errw = os.pipe()
    try:
        p = subprocess.Popen([prog], cwd=work, env=env, stdin=slave, stdout=slave,
                             stderr=errw, start_new_session=True)
    except FileNotFoundError:
        os.close(master), os.close(slave), os.close(errr), os.close(errw)
        return {'rc': 'notfound', 'out': b'', 'err': b''}
    os.close(slave)
    os.close(errw)
    err, out = bytearray(), bytearray()
    deadline = time.time() + 30
    fds = [errr, master]
    while fds and time.time() < deadline:
        r, _, _ = select.select(fds, [], [], 0.2)
        for fd in r:
            try:
                b = os.read(fd, 65536)
            except OSError:
                b = b''
            if not b:
                fds.remove(fd)
            elif fd == errr:
                err += b
            else:
                out += b
        if p.poll() is not None and errr not in fds:
            break
    if p.poll() is None:
        p.kill()
        p.wait()
        rc = 'timeout'
    else:
        rc = p.returncode
    os.close(errr)
    os.close(master)
    # pty output is not part of the comparison (it is a terminal)
    return {'rc': rc, 'out': b'', 'err': bytes(err)}


# ------------------------------------------------------------ comparison --

TEXT_BYTES = set(range(32, 127)) | {9, 10, 13, 12, 11}


def is_text(b):
    return all(c in TEXT_BYTES for c in b)


def su_records(b):
    """Split an XDR SU stream like fgettr would; None if it doesn't parse."""
    if len(b) < 240:
        return None
    ns = struct.unpack('>H', b[114:116])[0]
    trid = struct.unpack('>h', b[28:30])[0]
    per = 1 if trid == 201 else 2 if trid == 202 else 4
    rec = 240 + per * ns
    if len(b) % rec:
        return None
    recs = [b[i:i + rec] for i in range(0, len(b), rec)]
    for r in recs:
        if struct.unpack('>H', r[114:116])[0] != ns:
            return None
    return recs, ns, per


def classify(b, tags, prog=None, cmd=None):
    """Comparison kind of a golden stdout: text, su, f32 (raw native floats
    with tolerance: only with [f32] or when the inventory says the program
    writes binary_floats and it is the last stage) or bytes (exact)."""
    if tags.get('f32'):
        return 'f32'
    if tags.get('exact'):
        return 'bytes'
    if len(b) == 0 or is_text(b):
        return 'text'
    if su_records(b) is not None:
        return 'su'
    if (prog and cmd and len(b) % 4 == 0 and inventory_stdout(prog) == 'binary_floats'
            and last_stage_is(prog, cmd)):
        return 'f32'
    return 'bytes'


def norm_ws(s):
    lines = [re.sub(r'[ \t]+', ' ', l).rstrip() for l in s.split('\n')]
    return '\n'.join(lines)


def close(g, a, scale, rtol, atol):
    if math.isnan(g) or math.isnan(a):
        return math.isnan(g) and math.isnan(a), 0.0
    if math.isinf(g) or math.isinf(a):
        return g == a, 0.0
    allowed = rtol * scale + atol
    return abs(g - a) <= allowed, allowed


def compare_floats(g, a, rtol, atol, label):
    """g, a: lists of floats (same length). Returns None or (excess, message)."""
    finite = [abs(x) for x in g if not (math.isnan(x) or math.isinf(x))]
    scale = max(finite) if finite else 0.0
    worst = None
    for i, (x, y) in enumerate(zip(g, a)):
        ok, allowed = close(x, y, scale, rtol, atol)
        if not ok:
            excess = abs(x - y) - allowed if not (math.isnan(x) or math.isnan(y) or math.isinf(x) or math.isinf(y)) else math.inf
            if worst is None or excess > worst[0]:
                worst = (excess, i, x, y, allowed)
    if worst:
        _, i, x, y, allowed = worst
        return worst[0], '%s sample %d: golden %.9g actual %.9g |diff| %.3g > allowed %.3g' % (
            label, i, x, y, abs(x - y), allowed)
    return None


def compare_stdout(gold, act, tags, kind=None):
    kind = kind or classify(gold, tags)
    rtol = float(tags.get('rtol', 1e-5))
    atol = float(tags.get('atol', 1e-30))
    if kind == 'text':
        if tags.get('ws'):
            gs, as_ = norm_ws(gold.decode('latin-1')), norm_ws(act.decode('latin-1'))
        else:
            gs, as_ = gold.decode('latin-1'), act.decode('latin-1')
        if gs == as_:
            return None
        gl, al = gs.split('\n'), as_.split('\n')
        for i in range(max(len(gl), len(al))):
            x = gl[i] if i < len(gl) else '<EOF>'
            y = al[i] if i < len(al) else '<EOF>'
            if x != y:
                return 'text line %d: golden %r actual %r' % (i + 1, x[:120], y[:120])
        return 'text differs'
    if kind == 'bytes':
        if gold == act:
            return None
        n = min(len(gold), len(act))
        i = next((k for k in range(n) if gold[k] != act[k]), n)
        return 'bytes differ at offset %d (golden %d bytes, actual %d bytes)' % (i, len(gold), len(act))
    if kind == 'f32':
        if len(gold) != len(act):
            return 'float32 output length: golden %d bytes, actual %d bytes' % (len(gold), len(act))
        n = len(gold) // 4
        g = struct.unpack('<%df' % n, gold)
        a = struct.unpack('<%df' % n, act)
        m = compare_floats(g, a, rtol, atol, 'float32')
        return m[1] if m else None
    # SU traces
    gp = su_records(gold)
    ap = su_records(act)
    if ap is None:
        return 'actual stdout is not a valid SU (XDR) trace stream (%d bytes)' % len(act)
    grecs, gns, per = gp
    arecs, ans, _ = ap
    if len(grecs) != len(arecs) or gns != ans:
        return 'trace count/ns: golden %d traces ns=%d, actual %d traces ns=%d' % (
            len(grecs), gns, len(arecs), ans)
    worst = None
    nbad = 0
    for t, (gr, ar) in enumerate(zip(grecs, arecs)):
        if gr[:240] != ar[:240]:
            for off in range(240):
                if gr[off] != ar[off]:
                    break
            return 'trace %d header differs: %s' % (t + 1, header_field_diff(gr, ar, off))
        if per != 4:
            if gr != ar:
                return 'trace %d packed samples differ' % (t + 1)
            continue
        g = struct.unpack('>%df' % gns, gr[240:])
        a = struct.unpack('>%df' % gns, ar[240:])
        m = compare_floats(g, a, rtol, atol, 'trace %d' % (t + 1))
        if m and (worst is None or m[0] > worst[0]):
            worst = m
        if m:
            nbad += 1
    if worst:
        return 'worst: %s (%d of %d traces out of tolerance)' % (worst[1], nbad, len(grecs))
    return None


_KEYS = ('tracl:i tracr:i fldr:i tracf:i ep:i cdp:i cdpt:i trid:h nvs:h nhs:h duse:h offset:i '
         'gelev:i selev:i sdepth:i gdel:i sdel:i swdep:i gwdep:i scalel:h scalco:h sx:i sy:i gx:i '
         'gy:i counit:h wevel:h swevel:h sut:h gut:h sstat:h gstat:h tstat:h laga:h lagb:h delrt:h '
         'muts:h mute:h ns:u dt:u gain:h igc:h igi:h corr:h sfs:h sfe:h slen:h styp:h stas:h '
         'stae:h tatyp:h afilf:h afils:h nofilf:h nofils:h lcf:h hcf:h lcs:h hcs:h year:h day:h '
         'hour:h minute:h sec:h timbas:h trwf:h grnors:h grnofr:h grnlof:h gaps:h otrav:h d1:f '
         'f1:f d2:f f2:f ungpow:f unscale:f ntr:i mark:h shortpad:h').split()
HDR_FIELDS = []
_o = 0
for _kt in _KEYS:
    _k, _t = _kt.split(':')
    _n = 2 if _t in 'hu' else 4
    HDR_FIELDS.append((_k, _t, _o, _n))
    _o += _n
for _j in range(14):
    HDR_FIELDS.append(('unass[%d]' % _j, 'h', 212 + 2 * _j, 2))


def header_field_diff(g, a, off):
    """Describe the header field containing byte `off` (XDR big-endian)."""
    for k, t, o, n in HDR_FIELDS:
        if o <= off < o + n:
            fmt = {'i': '>i', 'h': '>h', 'u': '>H', 'f': '>f'}[t]
            gv = struct.unpack(fmt, g[o:o + n])[0]
            av = struct.unpack(fmt, a[o:o + n])[0]
            return '%s (bytes %d-%d): golden %r actual %r' % (k, o, o + n - 1, gv, av)
    return 'byte %d' % off


def first_line(b):
    for l in b.decode('latin-1').split('\n'):
        if l.strip():
            return l
    return ''
