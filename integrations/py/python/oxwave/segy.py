# Part of OxWave, a Rust port of CWP/SU.
"""SEG-Y files (rev 0/1 layout as SU's ``segyread`` / ``segywrite`` see them).

Contract (stable signatures/docstrings). Semantics are those of the original SU
programs with their default parameters, and results must be byte-identical:

* ``read_segy(p).traces.to_bytes()`` == stdout of ``segyread tape=p``, and
  ``read_segy(p).text_header_ascii()`` / ``binary_header``
  bytes == the ``hfile`` / ``bfile`` files segyread writes.
* ``write_segy(p, traces, text_header=H, binary_header=B)`` writes exactly the file
  ``segywrite tape=p [format=5]`` writes when ``hfile`` holds ``H`` and
  ``bfile`` holds ``B``.

SEG-Y here is always big-endian: this SUXDR build ignores segyread/segywrite ``endian=``
(``endian=0`` produces identical bytes), so there is no byte-order option.

IBM float conversion and the text tables MUST come from ``oxwave._native``
(``ibm_to_ieee``, ``ieee_to_ibm``; ``ebcdic_to_ascii`` = segyread's ``dd conv=ascii``;
``zebc`` = what segywrite of this SUXDR build applies after reading the ASCII header
with ``fgets`` in 80-byte lines): they are ports on the frozen oxwave-core.
The C sources are readable at ``../vendor/su-original/src/su/main/data_conversion/``
(segyread.c, segywrite.c; the SUXDR branches are the ones this build uses).
"""
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._headers import BHEADER_DTYPE, BHEADER_KEYS, BKEYS, HEADER_DTYPE
from .traces import Traces


__all__ = ["SegyFile", "read_segy", "write_segy"]


@dataclass
class SegyFile:
    """A SEG-Y file in memory.

    ``text_header``: the raw 3200-byte card-image header exactly as stored (EBCDIC).
    ``binary_header``: 400 bytes identical to the ``bfile`` segyread writes (on this
    little-endian machine that is the binary header in **native** byte order).
    ``traces``: the traces converted to SU exactly as ``segyread`` does (header mapping,
    sample conversion to float32).
    """
    text_header: bytes
    binary_header: bytes
    traces: Traces

    def text_header_ascii(self):
        """The text header converted to ASCII (bytes, 3200 long), exactly as segyread's
        ``hfile`` (``dd conv=ascii`` table)."""
        from . import _native

        return _native.ebcdic_to_ascii(self.text_header)

    def binary_header_dict(self):
        """``{key: int}`` for the 27 bhdr.h keys (``oxwave.BKEYS`` order) decoded from
        ``binary_header``."""
        bh = np.frombuffer(
            self.binary_header, dtype=BHEADER_DTYPE.newbyteorder("="), count=1,
        )[0]
        return {key: int(bh[key]) for key in BKEYS}


def read_segy(source):
    """Read a SEG-Y file (path or bytes-like) like ``segyread tape=...``.

    Sample formats 1 (IBM float), 2 (int32), 3 (int16), 5 (IEEE float) and
    8 (int8) as segyread handles them, ``ns`` from the binary header. Returns
    :class:`SegyFile`. A file shorter than 3600 bytes, or a format code segyread rejects,
    raises ValueError (message containing the segyread error text, e.g.
    ``format not SEGY standard (1, 2, 3, 5, or 8)``).
    """
    from . import _native

    raw = Path(source).read_bytes() if isinstance(source, (str, os.PathLike)) else bytes(source)
    if len(raw) < 3600:
        raise ValueError("segyread: file shorter than 3600-byte text and binary headers")

    bh = np.frombuffer(raw, dtype=BHEADER_DTYPE, count=1, offset=3200)[0]
    fmt = int(bh["format"])
    if fmt not in (1, 2, 3, 5, 8):
        raise ValueError("format not SEGY standard (1, 2, 3, 5, or 8)")
    # The bhdr key table calls hns a signed short; the C bhed member is unsigned.
    ns = int(bh["hns"]) & 0xffff
    if not ns:
        raise ValueError("samples/trace not set in binary header")

    # segyread fills only TAPEBHDR's 27 fields in its zeroed native bhed.
    native_bh = np.zeros(1, dtype=BHEADER_DTYPE.newbyteorder("="))
    for key in BKEYS:
        native_bh[key] = bh[key]

    nextended = int.from_bytes(raw[3504:3506], "big", signed=True)
    start = 3600 + max(nextended, 0) * 3200
    sample_dtype = {1: ">u4", 2: ">i4", 3: ">i2", 5: ">u4", 8: "i1"}[fmt]
    record_size = 240 + ns * np.dtype(sample_dtype).itemsize
    ntr = max(0, (len(raw) - start + record_size - 1) // record_size)
    header_bytes = np.empty((ntr, 240), dtype=np.uint8)
    data = np.empty((ntr, ns), dtype=np.float32)
    previous = bytes(record_size)
    for i in range(ntr):
        pos = start + i * record_size
        record = raw[pos:pos + record_size]
        # C reads into the same initially zeroed tapetr on every iteration,
        # including a final short read, whose remaining bytes are left intact.
        if len(record) < record_size:
            record += previous[len(record):]
        previous = record
        header_bytes[i] = np.frombuffer(record, dtype=np.uint8, count=240)
        samples = record[240:]
        if fmt == 1:
            data[i] = np.frombuffer(_native.ibm_to_ieee(samples), dtype=np.float32)
        elif fmt == 5:
            # Copy words rather than perform float arithmetic: even signaling
            # NaNs and their payloads must survive the IEEE format path.
            data[i].view(np.uint32)[:] = np.frombuffer(samples, dtype=">u4")
        else:
            data[i] = np.frombuffer(samples, dtype=sample_dtype)
            trwf = int.from_bytes(record[168:170], "big", signed=True)
            if trwf:
                with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                    scale = np.float32(np.exp2(np.float64(-trwf)))
                    data[i] *= scale

    # Only TAPEHDR's first 180 bytes were decoded/swapped by segyread. The
    # optional 60 bytes were memcpy'd into native words before puttr's XDR.
    header_bytes[:, 180:208] = (
        header_bytes[:, 180:208].reshape(ntr, 7, 4)[:, :, ::-1].reshape(ntr, 28)
    )
    header_bytes[:, 208:240] = (
        header_bytes[:, 208:240].reshape(ntr, 16, 2)[:, :, ::-1].reshape(ntr, 32)
    )
    headers = header_bytes.view(HEADER_DTYPE).reshape(ntr)
    headers["ns"] = ns
    return SegyFile(raw[:3200], native_bh.tobytes(), Traces(data, headers))


def write_segy(dest, traces, *, text_header=None, binary_header=None, format=None):
    """Write ``traces`` to ``dest`` (path or binary file object) like ``segywrite``.

    ``text_header``: the ASCII header text (``str`` or bytes), i.e. the content of the
    ``hfile`` segywrite would read; converted exactly as segywrite does (80-byte ``fgets``
    lines, then ``zebc``). None -> 3200 ASCII spaces.
    ``binary_header``: either 400 raw bytes in the **big-endian (XDR)** layout that
    ``segyhdrs`` writes and segywrite reads as its ``bfile`` (note: segyread's bfile /
    ``SegyFile.binary_header`` is native-endian, pass ``binary_header_dict()`` instead),
    or a dict of bhdr.h keys (missing keys 0), or None -> the binary header ``segyhdrs``
    (default parameters) would write for these traces.
    ``format``: None -> keep the binary header's format code (segywrite without
    ``format=``); 1 or 5 -> set the code (segywrite ``format=``); anything else ->
    ValueError. As in the original, samples are ALWAYS written as IBM floats: ``format``
    changes only the binary-header code (a ``format=5`` file read back is not ``traces``).
    Fields segywrite itself adjusts (this build always sets ``ntrpr`` = 1; ``tracr`` of the
    trace headers is renumbered 1..ntr) are adjusted the same way; a
    binary header with ``hns`` 0 -> ValueError (``bh.hns not set in binary header``).
    Returns None.
    """
    from . import _native

    if format is not None and format not in (1, 5):
        raise ValueError("format must be 1 or 5")
    if not isinstance(traces, Traces):
        raise TypeError("traces must be a Traces instance")

    if binary_header is None:
        if not traces.ntr:
            raise ValueError("can't get first trace")
        binary_header = {
            "jobid": 1, "lino": 1, "reno": 1, "format": 1,
            "hdt": int(traces.headers["dt"][0]),
            "hns": int(traces.headers["ns"][0]),
        }
    if isinstance(binary_header, dict):
        unknown = binary_header.keys() - set(BKEYS)
        if unknown:
            raise KeyError("unknown binary header key %r" % next(iter(unknown)))
        bh = bytearray(400)
        for key, typecode, offset in BHEADER_KEYS:
            size = 4 if typecode == "i" else 2
            value = int(binary_header.get(key, 0)) & ((1 << (8 * size)) - 1)
            bh[offset:offset + size] = value.to_bytes(size, "big")
    else:
        bh = bytearray(binary_header)
        if len(bh) != 400:
            raise ValueError("binary_header must contain exactly 400 bytes")
    # XDR decode/re-encode preserves the supplied spare words as well.
    bh[12:14] = b"\x00\x01"
    if format is not None:
        bh[24:26] = int(format).to_bytes(2, "big")
    ns = int.from_bytes(bh[20:22], "big")
    if not ns:
        raise ValueError("bh.hns not set in binary header")

    text = b" " * 3200 if text_header is None else text_header
    if isinstance(text, str):
        text = text.encode("ascii")
    text = bytes(text)
    cards = bytearray(3201)
    pos = 0
    for offset in range(0, 3200, 80):
        if pos < len(text):
            end = min(pos + 80, len(text))
            newline = text.find(b"\n", pos, end)
            if newline >= 0:
                end = newline + 1
            chunk = text[pos:end]
            cards[offset:offset + len(chunk)] = chunk
            cards[offset + len(chunk)] = 0
            pos = end
        # fgets retains newlines and embedded NULs. C replaces only the byte
        # found by strlen, and an EOF leaves the previous buffer untouched.
        cards[cards.index(0, offset)] = 32
    tape_text = _native.zebc(bytes(cards[:3200]))

    def _write(stream):
        stream.write(tape_text)
        stream.write(bh)
        samples = np.zeros(ns, dtype=np.float32)
        count = min(ns, traces.ns)
        for i in range(traces.ntr):
            header = bytearray(traces.headers[i:i + 1].tobytes())
            header[4:8] = (i + 1).to_bytes(4, "big")
            samples[:count] = traces.data[i, :count]
            stream.write(header)
            stream.write(_native.ieee_to_ibm(samples.tobytes()))

    if hasattr(dest, "write"):
        _write(dest)
    else:
        with open(dest, "wb") as stream:
            _write(stream)
