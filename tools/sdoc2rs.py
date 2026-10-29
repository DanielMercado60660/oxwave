#!/usr/bin/env python3
"""Extract the `char *sdoc[] = {...};` selfdoc array from an SU C source file
and print it as a Rust `const SDOC: &[&str] = &[...];` with identical bytes.

usage: tools/sdoc2rs.py path/to/prog.c > sdoc.rs   (or paste the output)
Python 3 stdlib only.
"""
import re
import sys


def c_unescape(body: str) -> bytes:
    out = bytearray()
    i = 0
    simple = {'n': 10, 't': 9, 'r': 13, 'a': 7, 'b': 8, 'f': 12, 'v': 11,
              '\\': 92, "'": 39, '"': 34, '?': 63}
    while i < len(body):
        c = body[i]
        if c != '\\':
            out += c.encode('latin-1')
            i += 1
            continue
        i += 1
        c = body[i]
        if c in simple:
            out.append(simple[c])
            i += 1
        elif c in '01234567':
            j = i
            while j < len(body) and j < i + 3 and body[j] in '01234567':
                j += 1
            out.append(int(body[i:j], 8) & 0xFF)
            i = j
        elif c == 'x':
            j = i + 1
            while j < len(body) and body[j] in '0123456789abcdefABCDEF':
                j += 1
            out.append(int(body[i + 1:j], 16) & 0xFF)
            i = j
        elif c == '\n':          # line continuation
            i += 1
        else:                    # unknown escape: C keeps the char
            out += c.encode('latin-1')
            i += 1
    return bytes(out)


def rust_lit(b: bytes) -> str:
    s = []
    for x in b:
        if x == 0x22:
            s.append('\\"')
        elif x == 0x5C:
            s.append('\\\\')
        elif x == 9:
            s.append('\\t')
        elif x == 10:
            s.append('\\n')
        elif x == 13:
            s.append('\\r')
        elif 32 <= x < 127:
            s.append(chr(x))
        else:
            raise SystemExit('sdoc2rs: non-ASCII byte 0x%02x in sdoc; convert it to ASCII first' % x)
    return '"' + ''.join(s) + '"'


def strip_comments(src: str) -> str:
    # remove /* */ comments outside string literals
    out, i, n = [], 0, len(src)
    while i < n:
        if src.startswith('/*', i):
            j = src.find('*/', i + 2)
            j = n if j < 0 else j + 2
            out.append(' ' + '\n' * src.count('\n', i, j))
            i = j
        elif src[i] == '"':
            j = i + 1
            while j < n and src[j] != '"':
                j += 2 if src[j] == '\\' else 1
            out.append(src[i:j + 1])
            i = j + 1
        elif src[i] == "'":
            j = i + 1
            while j < n and src[j] != "'":
                j += 2 if src[j] == '\\' else 1
            out.append(src[i:j + 1])
            i = j + 1
        else:
            out.append(src[i])
            i += 1
    return ''.join(out)


# Macros defined by the reference SU build (Makefile.config CWP_FLAGS); everything
# else (e.g. SU_LINE_HEADER, SLTSU_SEGY_H) is undefined.
DEFINED = {'SUXDR', 'CWP_LITTLE_ENDIAN', '__APPLE__'}


def eval_cond(expr: str) -> bool:
    e = expr.strip()
    e = re.sub(r'defined\s*\(\s*(\w+)\s*\)', lambda m: '1' if m.group(1) in DEFINED else '0', e)
    e = re.sub(r'defined\s+(\w+)', lambda m: '1' if m.group(1) in DEFINED else '0', e)
    e = e.replace('&&', ' and ').replace('||', ' or ')
    e = re.sub(r'!(?!=)', ' not ', e)
    e = re.sub(r'\b[A-Za-z_]\w*\b', lambda m: m.group(0) if m.group(0) in ('and', 'or', 'not') else '0', e)
    try:
        return bool(eval(e, {'__builtins__': {}}, {}))
    except Exception:
        sys.stderr.write('sdoc2rs: cannot evaluate #if %r; assuming false\n' % expr)
        return False


def preprocess(src: str) -> str:
    """Apply #if/#ifdef/#ifndef/#elif/#else/#endif (lines kept/dropped)."""
    out, stack = [], []   # stack of (active, taken, parent_active)
    active = True
    for line in src.split('\n'):
        m = re.match(r'\s*#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)', line)
        if not m:
            out.append(line if active else '')
            continue
        d, rest = m.group(1), m.group(2)
        if d in ('ifdef', 'ifndef', 'if'):
            if d == 'ifdef':
                c = rest.split()[0] in DEFINED
            elif d == 'ifndef':
                c = rest.split()[0] not in DEFINED
            else:
                c = eval_cond(rest)
            stack.append((active, c))
            active = active and c
        elif d == 'elif':
            parent, taken = stack[-1]
            c = (not taken) and eval_cond(rest)
            stack[-1] = (parent, taken or c)
            active = parent and c
        elif d == 'else':
            parent, taken = stack[-1]
            active = parent and not taken
            stack[-1] = (parent, True)
        else:
            parent, _ = stack.pop()
            active = parent
        out.append('')
    return '\n'.join(out)


def extract(src: str):
    src = preprocess(strip_comments(src))
    m = re.search(r'char\s*\*\s*sdoc\s*\[\s*\]\s*=\s*\{', src)
    if not m:
        raise SystemExit('no sdoc array found')
    i = m.end()
    entries, cur = [], None
    n = len(src)
    while i < n:
        c = src[i]
        if c == '"':
            j = i + 1
            while src[j] != '"':
                j += 2 if src[j] == '\\' else 1
            piece = c_unescape(src[i + 1:j])
            cur = piece if cur is None else cur + piece   # adjacent literals concatenate
            i = j + 1
        elif c == ',':
            if cur is not None:
                entries.append(cur)
            cur = None
            i += 1
        elif src.startswith('NULL', i):
            if cur is not None:
                entries.append(cur)
            return entries
        elif c == '}':
            if cur is not None:
                entries.append(cur)
            return entries
        else:
            i += 1
    raise SystemExit('unterminated sdoc array')


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    src = open(sys.argv[1], encoding='latin-1').read()
    entries = extract(src)
    print('/// Self-documentation, byte-identical to the C `sdoc[]` array.')
    print('#[rustfmt::skip]')
    print('const SDOC: &[&str] = &[')
    for e in entries:
        print('    ' + rust_lit(e) + ',')
    print('];')


if __name__ == '__main__':
    main()
