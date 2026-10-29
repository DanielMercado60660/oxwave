# Part of OxWave, a Rust port of CWP/SU.
"""SU files and streams, header tools, and ``surange``.

Contract (signatures and docstrings are stable). Pure NumPy:
this module must not require ``oxwave._native`` (it may use it opportunistically only
if results stay bit-identical). Byte format and errors: see :class:`oxwave.Traces`
(``Traces.from_bytes`` / ``Traces.to_bytes`` are the reference codec).
"""
import io
import math
import os

import numpy as np

from ._headers import HEADER_DTYPE, KEYS, KEY_TYPE

from .traces import SUFormatError, Traces  # noqa: F401


__all__ = ["read_su", "write_su", "iter_su", "concat", "gethw", "sethw", "surange", "surange_text"]


def read_su(source, *, max_traces=None):
    """Read a whole SU stream.

    ``source``: a path (``str`` / ``os.PathLike``), a bytes-like object, or a binary
    file object (read to EOF from its current position). ``max_traces``: if given
    (>= 0), stop after that many traces (the rest of a file is not read).
    Returns :class:`Traces`. Errors exactly as ``Traces.from_bytes`` (SUFormatError);
    a missing path raises FileNotFoundError.
    """
    if max_traces is not None:
        if max_traces < 0:
            raise ValueError("max_traces must be nonnegative")
        if max_traces == 0:
            return Traces(np.empty((0, 0), dtype=np.float32))
        it = iter_su(source, chunk=1)
        try:
            parts = []
            for _ in range(max_traces):
                try:
                    parts.append(next(it))
                except StopIteration:
                    break
            return concat(parts)
        finally:
            it.close()
    if isinstance(source, (bytes, bytearray, memoryview)):
        return Traces.from_bytes(source)
    if isinstance(source, (str, os.PathLike)):
        with open(source, "rb") as f:
            return Traces.from_bytes(f.read())
    return Traces.from_bytes(source.read())


def write_su(dest, traces):
    """Write ``traces`` as an SU stream to ``dest`` (path, created/truncated; or a
    binary file object, written at its position). Bytes are exactly
    ``traces.to_bytes()``. Returns None."""
    data = traces.to_bytes()
    if isinstance(dest, (str, os.PathLike)):
        with open(dest, "wb") as f:
            f.write(data)
    else:
        dest.write(data)


def iter_su(source, chunk=1000):
    """Iterate over an SU stream in chunks without loading it all.

    ``source``: path or binary file object (bytes-like also accepted). Yields
    :class:`Traces` of up to ``chunk`` (>= 1) traces each, in order; concatenating them
    gives exactly ``read_su(source)``. Errors as ``Traces.from_bytes``, raised when the
    offending trace is reached (traces before it have already been yielded). An empty
    source raises ``SUFormatError("fgettr.c: bad first header")`` on the first ``next``.
    """
    if chunk < 1:
        raise ValueError("chunk must be at least 1")
    if isinstance(source, (bytes, bytearray, memoryview)):
        f, owned = io.BytesIO(bytes(source)), True
    elif isinstance(source, (str, os.PathLike)):
        f, owned = open(source, "rb"), True
    else:
        f, owned = source, False
    try:
        def read_n(n):
            result = bytearray()
            while len(result) < n:
                block = f.read(n - len(result))
                if not block:
                    break
                result.extend(block)
            return bytes(result)

        header = read_n(240)
        if len(header) < 240:
            Traces.from_bytes(header)
        h = np.frombuffer(header, dtype=HEADER_DTYPE)[0]
        ns = int(h["ns"])
        if int(h["trid"]) in (201, 202) or ns == 0:
            Traces.from_bytes(header)
        first = header + read_n(4 * ns)
        Traces.from_bytes(first)
        pending, index = [first], 1
        def decode(records):
            if records[0] is first:
                return Traces.from_bytes(b"".join(records))
            return Traces.from_bytes(first + b"".join(records))[1:]

        while True:
            if len(pending) == chunk:
                yield decode(pending)
                pending = []
            header = read_n(240)
            if len(header) < 240:
                break
            data = read_n(4 * ns)
            error = None
            if len(data) < 4 * ns:
                error = SUFormatError("truncated SU stream: trace %d is incomplete" % index)
            else:
                current_ns = int(np.frombuffer(header, dtype=HEADER_DTYPE)[0]["ns"])
                if current_ns != ns:
                    error = SUFormatError("fgettr.c: on trace #%d, number of samples in header (%d) differs "
                                          "from number for first trace (%d)" % (index, current_ns, ns))
            if error is not None:
                if pending:
                    yield decode(pending)
                raise error
            pending.append(header + data)
            index += 1
        if pending:
            yield decode(pending)
    finally:
        if owned:
            f.close()


def concat(parts):
    """Concatenate an iterable of :class:`Traces` (same ``ns``; ValueError otherwise)
    into one, headers and samples copied in order. Empty iterable -> ValueError."""
    parts = list(parts)
    if not parts:
        raise ValueError("cannot concatenate an empty iterable")
    if any(part.ns != parts[0].ns for part in parts):
        raise ValueError("traces must have the same ns")
    return Traces(np.concatenate([p.data for p in parts]), np.concatenate([p.headers for p in parts]))


def gethw(traces, keys):
    """Header values by SU key name, like ``sugethw key=... output=geom``.

    ``keys``: a comma-separated string (``"sx,gx"``) or a sequence of names.
    Returns ``{key: numpy array of length ntr}`` in the given key order, each array in
    the key's native type (int32/int16/uint16/float32). Unknown key -> KeyError.
    """
    return {k: traces[k] for k in _keys(keys)}


def sethw(traces, key, a=0, b=0, c=0, d=0, j=None):
    """Set header words like ``sushw key=... a=... b=... c=... d=... j=...``.

    For each key k (``key``: one name, a comma string or a sequence; ``a``..``j``:
    scalars or sequences with one entry per key, scalars broadcast) and trace
    ``itr`` (0-based), in C doubles: ``i = itr + d``;
    ``val = a + b*mod(i, j) + c*(int)(i/j)`` where ``mod`` is cwp's ``x - y*floor(x/y)``
    and ``(int)`` truncates toward zero (NOT Python's floor ``//``); ``j`` None or 0 means
    ``(double)ULONG_MAX``. Conversion to the field type exactly as sushw.c on this arm64
    build: int fields truncate and saturate to int32 (NaN -> 0); short fields saturate to
    int32 then keep the low 16 bits; unsigned short fields saturate to uint32 then keep
    the low 16 bits; float fields round to nearest float32. Returns a new
    :class:`Traces`; the result must be byte-identical to the original ``sushw``.
    Unknown key -> KeyError; mismatched list lengths -> ValueError.
    """
    keys = _keys(key)
    for k in keys:
        if k not in KEY_TYPE:
            raise KeyError("unknown SU header key %r" % (k,))
    params = [_params(v, len(keys)) for v in (a, b, c, d, j)]
    out = traces.copy()
    for idx, k in enumerate(keys):
        aa, bb, cc, dd, jj = (p[idx] for p in params)
        aa, bb, cc, dd = map(float, (aa, bb, cc, dd))
        jj = float(jj) if jj is not None and jj != 0 else float(2**64 - 1)
        values = []
        for itr in range(traces.ntr):
            i = float(itr) + dd
            q = _saturate(i / jj, -2147483648, 2147483647)
            val = (aa + bb * (i % jj)) + cc * q
            kind = KEY_TYPE[k]
            if kind == "f":
                values.append(val)
            elif kind == "u":
                values.append(_saturate(val, 0, 4294967295) & 65535)
            else:
                v = _saturate(val, -2147483648, 2147483647)
                values.append(((v + 32768) % 65536 - 32768) if kind == "h" else v)
        with np.errstate(over="ignore", invalid="ignore"):
            out.headers[k] = values
    return out


def surange(traces):
    """Header ranges like ``surange`` (no parameters), as data.

    Returns a dict: ``{"ntr": int, "keys": {key: (min, max)}, "first": {key: value},
    "last": {key: value}}`` where ``keys`` lists, in SU key order, exactly the keys that
    ``surange`` prints (a key is printed when any trace has it nonzero), values in the
    key's native Python type (int or float). For a key printed with a single value
    (constant), min == max. ``first``/``last`` are the ``(...)`` values surange shows
    for that key (first and last trace values; equal to min for constant keys).
    Also includes ``"coords"`` -> ``{"shot": {...}, "receiver": {...}, "midpoint": {...}}``
    with ``North/South/East/West`` as ``(x, y)`` tuples when surange prints those
    blocks, else ``{}``.
    """
    if not traces.ntr:
        raise ValueError("surange requires at least one trace")
    ranges, first, last = {}, {}, {}
    for k in KEYS:
        vals = traces[k]
        lo = hi = vals[0].item()
        for v in vals[1:]:
            v = v.item()
            if v < lo:
                lo = v
            if v > hi:
                hi = v
        if lo != 0 or hi != 0:
            ranges[k] = (lo, hi)
            first[k] = vals[0].item() if lo < hi else lo
            last[k] = vals[-1].item() if lo < hi else lo
    coords = {}
    points = {"shot": [], "receiver": [], "midpoint": []}
    for h in traces.headers:
        scale = int(h["scalco"])
        scale = float(scale) if scale > 0 else (-1.0 / scale if scale < 0 else 1.0)
        sx, sy, gx, gy = (int(h[k]) * scale for k in ("sx", "sy", "gx", "gy"))
        points["shot"].append((sx, sy))
        points["receiver"].append((gx, gy))
        points["midpoint"].append((0.5 * (sx + gx), 0.5 * (sy + gy)))
    for name, pts in points.items():
        north = south = east = west = pts[0]
        for pt in pts[1:]:
            if pt[1] > north[1]:
                north = pt
            if pt[1] < south[1]:
                south = pt
            if pt[0] > east[0]:
                east = pt
            if pt[0] < west[0]:
                west = pt
        if north[1] != 0 or south[1] != 0 or east[0] != 0 or west[0] != 0:
            coords[name] = {"North": north, "South": south, "East": east, "West": west}
    return {"ntr": traces.ntr, "keys": ranges, "first": first, "last": last, "coords": coords}


def surange_text(traces):
    """Exactly the stdout bytes (as ``str``) the original ``surange`` prints for this
    stream with no parameters, e.g. ``"8 traces:\\ntracl    1 8 (1 - 8)\\n..."``."""
    result = surange(traces)
    text = "%d traces:\n" % result["ntr"]
    for k, (lo, hi) in result["keys"].items():
        fmt = "%f" if KEY_TYPE[k] == "f" else "%d"
        text += "%-8s " % k + (fmt % lo)
        if lo < hi:
            text += " " + (fmt % hi) + " (" + (fmt % result["first"][k])
            text += " - " + (fmt % result["last"][k]) + ")"
        text += "\n"
    for name in ("shot", "receiver", "midpoint"):
        if name in result["coords"]:
            c = result["coords"][name]
            text += "\n%s coordinate limits:\n" % name.capitalize()
            text += "\tNorth(%g,%g) South(%g,%g) East(%g,%g) West(%g,%g)\n" % (
                *c["North"], *c["South"], *c["East"], *c["West"])
    return text


def _keys(keys):
    return keys.split(",") if isinstance(keys, str) else list(keys)


def _params(value, n):
    if value is None or np.isscalar(value):
        return [value] * n
    result = list(value)
    if len(result) != n:
        raise ValueError("parameter lengths must match the number of keys")
    return result


def _saturate(value, lo, hi):
    if math.isnan(value):
        return 0
    if value <= lo:
        return lo
    if value >= hi:
        return hi
    return math.trunc(value)
