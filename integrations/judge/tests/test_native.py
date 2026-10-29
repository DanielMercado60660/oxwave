# Tests for (py/src/lib.rs): the oxwave._native extension.
import glob
import os
import random
import zipfile

import numpy as np
import pytest

from conftest import FIX, ROOT, fixture_bytes

import oxwave
from oxwave.traces import Traces

SEG = FIX / "segy"


def sgy_split(raw, ns):
    """(text 3200, binary 400, [(trace header 240, data 4*ns)])"""
    body = raw[3600:]
    rec = 240 + 4 * ns
    assert len(body) % rec == 0
    return raw[:3200], raw[3200:3600], [(body[i:i + 240], body[i + 240:i + rec]) for i in range(0, len(body), rec)]


def test_version(native):
    assert native.version() == "0.1.0"


def test_header_keys_match_frozen_table(native):
    assert [tuple(x) for x in native.header_keys()] == list(oxwave.HEADER_KEYS)


def test_bheader_keys_match_frozen_table(native):
    assert [tuple(x) for x in native.bheader_keys()] == list(oxwave.BHEADER_KEYS)


@pytest.mark.parametrize("name", ["plane", "wild"])
def test_ibm_to_ieee_matches_segyread(native, name):
    raw = (SEG / ("%s.sgy" % name)).read_bytes()
    want = Traces.from_bytes((SEG / ("%s.sgy.su" % name)).read_bytes())
    _, _, trs = sgy_split(raw, want.ns)
    got = native.ibm_to_ieee(b"".join(d for _, d in trs))
    assert isinstance(got, bytes)
    assert got == want.data.astype("<f4").tobytes()


@pytest.mark.parametrize("name", ["plane", "wild"])
def test_ieee_to_ibm_matches_segywrite(native, name):
    src = Traces.from_bytes((SEG / ("%s.su" % name)).read_bytes() if name == "wild" else fixture_bytes("plane.su"))
    raw = (SEG / ("%s.sgy" % name)).read_bytes()
    _, _, trs = sgy_split(raw, src.ns)
    assert native.ieee_to_ibm(src.data.astype("<f4").tobytes()) == b"".join(d for _, d in trs)


def test_little_endian_words(native):
    src = Traces.from_bytes(fixture_bytes("gather.su")).data.astype("<f4")
    be = native.ieee_to_ibm(src.tobytes())
    le = native.ieee_to_ibm(src.tobytes(), big_endian=False)
    assert le == np.frombuffer(be, ">u4").astype("<u4").tobytes()
    assert native.ibm_to_ieee(le, False) == native.ibm_to_ieee(be, True)
    assert native.ibm_to_ieee(be) == native.ibm_to_ieee(be, big_endian=True)


def test_ebcdic_to_ascii_matches_segyread_hfile(native):
    for name in ("plane", "plane5", "wild"):
        raw = (SEG / ("%s.sgy" % name)).read_bytes()
        assert native.ebcdic_to_ascii(raw[:3200]) == (SEG / ("%s.sgy.hfile" % name)).read_bytes()


def test_zebc_matches_segywrite_text_header(native):
    hdr = (SEG / "header").read_bytes()
    assert native.zebc(hdr) == (SEG / "plane.sgy").read_bytes()[:3200]
    assert native.zebc(bytes(range(256)))[:32] == b"\x20" * 32


def test_dd_tables_are_inverse_on_text(native):
    txt = bytes(range(32, 127)) * 4
    assert native.ebcdic_to_ascii(native.ascii_to_ebcdic(txt)) == txt
    assert native.ascii_to_ebcdic(b"C 1 CLIENT") == bytes([0xC3, 0x40, 0xF1, 0x40, 0xC3, 0xD3, 0xC9, 0xC5, 0xD5, 0xE3])
    assert len(native.ascii_to_ebcdic(bytes(range(256)))) == 256


@pytest.mark.parametrize("name", ["plane.su", "spike.su", "gather.su", "one.su", "segy/wild.su"])
def test_su_decode_matches_reference(native, name):
    b = fixture_bytes(name)
    ref = Traces.from_bytes(b)
    h, s, ns, ntr = native.su_decode(b)
    assert (ns, ntr) == (ref.ns, ref.ntr)
    assert h == ref.headers.tobytes()
    assert s == ref.data.astype("=f4").tobytes()
    assert native.su_encode(h, s, ns) == b


def test_su_decode_errors_match_reference(native):
    b = fixture_bytes("plane.su")
    bad = bytearray(b)
    bad[496 + 114:496 + 116] = (63).to_bytes(2, "big")  # trace 1 (0-based) claims ns=63
    cases = [b"", b[:239], b[:240], b[:300], b[:len(b) - 10], bytes(bad)]
    for c in cases:
        with pytest.raises(ValueError) as ref:
            Traces.from_bytes(c)
        with pytest.raises(ValueError) as got:
            native.su_decode(c)
        assert str(got.value) == str(ref.value), (len(c), str(got.value))


def test_length_errors(native):
    for f in (native.ibm_to_ieee, native.ieee_to_ibm):
        with pytest.raises(ValueError):
            f(b"\x00\x01\x02")
    h, s, ns, ntr = native.su_decode(fixture_bytes("plane.su"))
    with pytest.raises(ValueError):
        native.su_encode(h[:-1], s, ns)
    with pytest.raises(ValueError):
        native.su_encode(h, s[:-4], ns)
    with pytest.raises(ValueError):
        native.su_encode(h, s, 0)


def test_garbage_never_panics(native):
    rnd = random.Random(7)
    for n in (240, 241, 496, 1000, 5000):
        for _ in range(20):
            buf = bytes(rnd.getrandbits(8) for _ in range(n))
            try:
                native.su_decode(buf)
            except ValueError:
                pass  # anything but a panic (PanicException is a BaseException)


def test_large_stream_roundtrip(native):
    rng = np.random.default_rng(3)
    tr = Traces(rng.standard_normal((2000, 501)).astype(np.float32), dt=0.002)
    tr["offset"] = np.arange(2000) * 12
    b = tr.to_bytes()
    h, s, ns, ntr = native.su_decode(b)
    assert (ns, ntr) == (501, 2000)
    assert s == tr.data.tobytes()
    assert native.su_encode(h, s, ns) == b


def test_wheel_contents():
    wheels = sorted(glob.glob(str(ROOT / "target" / "wheels" / "oxwave-*.whl")), key=os.path.getmtime)
    assert wheels, "no wheel: run tools/build_ext.py"
    with zipfile.ZipFile(wheels[-1]) as z:
        names = set(z.namelist())
        for need in ("oxwave/__init__.py", "oxwave/_programs.json", "oxwave/traces.py",
                     "oxwave/mcp_server/server.py", "oxwave/_native.pyi"):
            assert need in names, need
        assert any(n.startswith("oxwave/_native") and n.endswith(".so") for n in names)
        ep = [n for n in names if n.endswith("entry_points.txt")]
        assert ep and "oxwave-mcp=oxwave.mcp_server.server:main" in z.read(ep[0]).decode().replace(" ", "")
