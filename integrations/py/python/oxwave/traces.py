# Part of OxWave, a Rust port of CWP/SU.
"""`Traces`: an in-memory SU gather (headers + float32 samples) and the SU stream codec.

The codec here is the NumPy reference implementation of the XDR (big-endian) SU trace
format used by this build (``fgettr.c`` / ``fputtr.c``): each trace is a 240-byte
big-endian header followed by ``ns`` big-endian IEEE float32 samples, ``ns`` fixed by
the first trace. It is bit-exact both ways.
"""
import numpy as np

from ._headers import HEADER_DTYPE, HDRBYTES, KEYS, KEY_TYPE, native_type

__all__ = ["Traces", "SUFormatError"]


class SUFormatError(ValueError):
    """Raised when bytes are not a valid SU (XDR) trace stream. Messages follow
    fgettr.c where it has one (e.g. ``fgettr.c: bad first header``)."""


def _record_dtype(ns):
    return np.dtype([("h", HEADER_DTYPE), ("d", ">f4", (ns,))]) if ns > 0 else None


class Traces:
    """A gather of SU traces held in memory.

    Attributes
    ----------
    headers : numpy.ndarray
        Structured array, shape ``(ntr,)``, dtype :data:`oxwave.HEADER_DTYPE`
        (big-endian fields at the segy.h offsets, including ``unass``).
    data : numpy.ndarray
        ``float32`` samples, shape ``(ntr, ns)``, C-contiguous, native byte order.

    Construction
    ------------
    ``Traces(data, headers=None, *, dt=None)``: ``data`` is array-like of shape
    ``(ntr, ns)`` (a 1-D array is one trace). If ``headers`` is None, fresh headers are
    made like ``suaddhead``: zeros except ``tracl = tracr = 1..ntr``, ``trid = 1``,
    ``ns``, and ``dt = round(dt * 1e6)`` microseconds when ``dt`` (seconds) is given.
    If ``headers`` is given it must have ``ntr`` entries; it is converted to
    ``HEADER_DTYPE`` (copied) and ``dt`` is ignored.

    Header access by SU key name
    ----------------------------
    ``tr["offset"]`` -> native-typed copy of that header field for all traces;
    ``tr["offset"] = values`` broadcasts and casts into the field;
    ``tr.header(i)`` -> ``{key: value}`` for trace ``i`` (the 80 SU keys);
    ``tr[i]``, ``tr[a:b]``, ``tr[mask_or_index_array]`` -> a new ``Traces`` (copies).
    """

    __slots__ = ("headers", "data")

    def __init__(self, data, headers=None, *, dt=None):
        d = np.asarray(data, dtype=np.float32)
        if d.ndim == 1:
            d = d.reshape(1, -1)
        if d.ndim != 2:
            raise ValueError("Traces data must be 1-D or 2-D (ntr, ns), got shape %r" % (d.shape,))
        ntr, ns = d.shape
        if ns > 65535:
            raise ValueError("ns=%d exceeds SU_NFLTS=65535" % ns)
        self.data = np.ascontiguousarray(d).copy() if d is data else np.ascontiguousarray(d)
        if headers is None:
            h = np.zeros(ntr, dtype=HEADER_DTYPE)
            h["tracl"] = np.arange(1, ntr + 1)
            h["tracr"] = np.arange(1, ntr + 1)
            h["trid"] = 1
            h["ns"] = ns
            if dt is not None:
                h["dt"] = int(round(float(dt) * 1e6))
        else:
            h = np.asarray(headers)
            if h.dtype != HEADER_DTYPE:
                h = h.astype(HEADER_DTYPE)
            else:
                h = h.copy()
            h = h.reshape(-1)
            if h.shape[0] != ntr:
                raise ValueError("headers has %d entries but data has %d traces" % (h.shape[0], ntr))
        self.headers = h

    # ---------------------------------------------------------------- shape --
    @property
    def ntr(self):
        """Number of traces."""
        return int(self.data.shape[0])

    @property
    def ns(self):
        """Samples per trace (columns of ``data``)."""
        return int(self.data.shape[1])

    @property
    def dt(self):
        """Sample interval in seconds from the first header's ``dt`` (microseconds);
        None when there are no traces or ``dt`` is 0."""
        if self.ntr == 0:
            return None
        v = int(self.headers["dt"][0])
        return v * 1e-6 if v else None

    def __len__(self):
        return self.ntr

    def keys(self):
        """The 80 SU header key names in key-index order."""
        return KEYS

    # --------------------------------------------------------------- access --
    def __getitem__(self, item):
        if isinstance(item, str):
            if item not in KEY_TYPE:
                raise KeyError("unknown SU header key %r" % (item,))
            return self.headers[item].astype(native_type(item))
        if isinstance(item, (int, np.integer)):
            i = int(item)
            if i < -self.ntr or i >= self.ntr:
                raise IndexError("trace index %d out of range (ntr=%d)" % (i, self.ntr))
            item = slice(i, i + 1) if i != -1 else slice(i, None)
        return Traces(self.data[item].copy(), self.headers[item].copy())

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError("assign header fields by SU key name, e.g. tr['offset'] = ...")
        if key not in KEY_TYPE:
            raise KeyError("unknown SU header key %r" % (key,))
        v = np.asarray(value)
        if KEY_TYPE[key] != "f" and v.dtype.kind == "f":
            v = np.trunc(v)  # C assignment of a double to an integer field truncates
        self.headers[key] = np.broadcast_to(v, (self.ntr,)).astype(native_type(key), casting="unsafe")

    def header(self, i):
        """Header of trace ``i`` as ``{key: int|float}`` over the 80 SU keys."""
        h = self.headers[i]
        return {k: (float(h[k]) if KEY_TYPE[k] == "f" else int(h[k])) for k in KEYS}

    def copy(self):
        return Traces(self.data.copy(), self.headers.copy())

    def __repr__(self):
        dt = self.dt
        return "Traces(ntr=%d, ns=%d, dt=%s)" % (self.ntr, self.ns, "%g" % dt if dt else "unset")

    def identical(self, other):
        """True if both encode to the same SU bytes (headers and samples bit-exact)."""
        return isinstance(other, Traces) and self.to_bytes() == other.to_bytes()

    # ---------------------------------------------------------------- codec --
    def to_bytes(self):
        """Encode as an SU (XDR, big-endian) stream, like ``puttr``: every header
        verbatim, then ``ns`` samples. ValueError if the first header's ``ns`` differs
        from ``data.shape[1]`` (puttr would write a stream fgettr cannot read back)."""
        if self.ntr == 0:
            return b""
        ns = self.ns
        if int(self.headers["ns"][0]) != ns:
            raise ValueError("first header ns=%d differs from data.shape[1]=%d" % (int(self.headers["ns"][0]), ns))
        if ns == 0:
            raise ValueError("cannot encode traces with ns=0")
        rec = np.empty(self.ntr, dtype=_record_dtype(ns))
        rec["h"] = self.headers
        rec["d"] = self.data
        return rec.tobytes()

    @classmethod
    def from_bytes(cls, buf):
        """Decode an SU (XDR) stream like ``fgettr`` (fixed length from the first trace).

        Errors (SUFormatError): empty or < 240 bytes -> ``fgettr.c: bad first header``;
        ``trid`` 201/202 (packed) -> ``packed SU traces (trid=201/202) are not supported``;
        no complete first-trace data -> ``fgettr.c: no data on first trace``; a later
        header with a different ``ns`` -> ``fgettr.c: on trace #<k>, number of samples in
        header (<ns>) differs from number for first trace (<nsfirst>)`` (k = 0-based index
        of that trace, as fgettr prints it); a truncated later trace ->
        ``truncated SU stream: trace <k> is incomplete``. Trailing bytes shorter than a
        header are ignored, as fgettr does.
        """
        b = memoryview(buf).cast("B") if not isinstance(buf, (bytes, bytearray)) else buf
        n = len(b)
        if n < HDRBYTES:
            raise SUFormatError("fgettr.c: bad first header")
        h0 = np.frombuffer(bytes(b[:HDRBYTES]), dtype=HEADER_DTYPE)[0]
        if int(h0["trid"]) in (201, 202):
            raise SUFormatError("packed SU traces (trid=201/202) are not supported")
        ns = int(h0["ns"])
        nb = HDRBYTES + 4 * ns
        if ns == 0 or n < nb:
            if ns > 0 and n > HDRBYTES:
                raise SUFormatError("fgettr.c: first trace: tried to read %d bytes read %d bytes"
                                    % (4 * ns, n - HDRBYTES))
            raise SUFormatError("fgettr.c: no data on first trace")
        ntr = n // nb
        rem = n - ntr * nb
        if rem >= HDRBYTES:
            raise SUFormatError("truncated SU stream: trace %d is incomplete" % ntr)
        rec = np.frombuffer(bytes(b[:ntr * nb]) if rem else b, dtype=_record_dtype(ns), count=ntr)
        hs = rec["h"]
        bad = np.nonzero(hs["ns"] != ns)[0]
        if bad.size:
            k = int(bad[0])
            raise SUFormatError("fgettr.c: on trace #%d, number of samples in header (%d) differs "
                                "from number for first trace (%d)" % (k, int(hs["ns"][k]), ns))
        t = cls.__new__(cls)
        t.headers = hs.copy()
        t.data = np.ascontiguousarray(rec["d"].astype(np.float32))
        return t
