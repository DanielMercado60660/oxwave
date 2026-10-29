# Tests for (py/python/oxwave/su.py). Pure NumPy; no binaries.
import io

import numpy as np
import pytest

from conftest import FIX, fixture_bytes

import oxwave
from oxwave import SUFormatError, Traces

SU_FIXTURES = ["plane.su", "spike.su", "gather.su", "one.su"]


def parse_geom(path, keys):
    rows = [l.split() for l in (FIX / path).read_text().splitlines() if l.strip()]
    return {k: [r[i] for r in rows] for i, k in enumerate(keys)}


@pytest.mark.parametrize("name", SU_FIXTURES)
def test_read_su_path_bytes_fileobj(name):
    b = fixture_bytes(name)
    ref = Traces.from_bytes(b)
    for src in (FIX / name, str(FIX / name), b, bytearray(b), io.BytesIO(b)):
        t = oxwave.read_su(src)
        assert isinstance(t, Traces)
        assert t.to_bytes() == b
        assert t.headers.dtype == oxwave.HEADER_DTYPE and t.data.dtype == np.float32
    assert ref.identical(oxwave.read_su(FIX / name))


def test_read_su_max_traces():
    b = fixture_bytes("gather.su")
    t = oxwave.read_su(FIX / "gather.su", max_traces=10)
    assert t.ntr == 10 and t.to_bytes() == b[:10 * (240 + 4 * 250)]
    assert oxwave.read_su(io.BytesIO(b), max_traces=500).ntr == 120
    assert oxwave.read_su(b, max_traces=0).ntr == 0


def test_read_su_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        oxwave.read_su(tmp_path / "nope.su")
    with pytest.raises(SUFormatError, match="bad first header"):
        oxwave.read_su(b"")
    b = fixture_bytes("plane.su")
    with pytest.raises(SUFormatError, match="incomplete"):
        oxwave.read_su(b[:-10])


@pytest.mark.parametrize("name", SU_FIXTURES)
def test_write_su_bit_exact(tmp_path, name):
    b = fixture_bytes(name)
    t = oxwave.read_su(b)
    oxwave.write_su(tmp_path / "x.su", t)
    assert (tmp_path / "x.su").read_bytes() == b
    oxwave.write_su(str(tmp_path / "y.su"), t)
    assert (tmp_path / "y.su").read_bytes() == b
    buf = io.BytesIO()
    buf.write(b"JUNK")
    oxwave.write_su(buf, t)
    assert buf.getvalue() == b"JUNK" + b


@pytest.mark.parametrize("chunk", [1, 7, 120, 1000])
def test_iter_su_chunks(chunk):
    b = fixture_bytes("gather.su")
    parts = list(oxwave.iter_su(FIX / "gather.su", chunk=chunk))
    assert all(isinstance(p, Traces) for p in parts)
    assert [p.ntr for p in parts] == [min(chunk, 120 - i) for i in range(0, 120, chunk)]
    assert b"".join(p.to_bytes() for p in parts) == b
    assert oxwave.concat(parts).to_bytes() == b
    with open(FIX / "gather.su", "rb") as f:
        assert b"".join(p.to_bytes() for p in oxwave.iter_su(f, chunk=chunk)) == b


def test_iter_su_error_after_good_traces():
    b = bytearray(fixture_bytes("plane.su"))
    k = 5
    b[k * 496 + 114:k * 496 + 116] = (63).to_bytes(2, "big")
    it = oxwave.iter_su(bytes(b), chunk=2)
    got = [next(it).ntr, next(it).ntr]
    assert got == [2, 2]
    with pytest.raises(SUFormatError) as e:
        for _ in it:
            pass
    assert str(e.value) == ("fgettr.c: on trace #5, number of samples in header (63) differs "
                            "from number for first trace (64)")
    it = oxwave.iter_su(b"")
    with pytest.raises(SUFormatError, match="bad first header"):
        next(it)


def test_concat_errors():
    a = oxwave.read_su(FIX / "plane.su")
    s = oxwave.read_su(FIX / "spike.su")
    with pytest.raises(ValueError):
        oxwave.concat([a, s])
    with pytest.raises(ValueError):
        oxwave.concat([])
    c = oxwave.concat([a[0:3], a[3:]])
    assert c.identical(a) and c is not a


def test_gethw_matches_sugethw_plane():
    t = oxwave.read_su(FIX / "plane.su")
    keys = ["tracl", "offset", "sx", "gx", "cdp", "d2", "dt", "ns"]
    g = oxwave.gethw(t, ",".join(keys))
    assert list(g) == keys
    want = parse_geom("plane.gethw.txt", keys)
    for k in keys:
        if k == "d2":
            assert np.array_equal(g[k], np.array([float(x) for x in want[k]], np.float32))
            assert g[k].dtype == np.float32
        else:
            assert g[k].tolist() == [int(x) for x in want[k]], k
    assert g["ns"].dtype == np.uint16 and g["offset"].dtype == np.int32


def test_gethw_matches_sugethw_gather():
    t = oxwave.read_su(FIX / "gather.su")
    keys = ["fldr", "tracf", "offset", "sx", "gx", "cdp", "trid"]
    g = oxwave.gethw(t, keys)
    want = parse_geom("gather.gethw.txt", keys)
    for k in keys:
        assert g[k].tolist() == [int(x) for x in want[k]], k
    assert g["trid"].dtype == np.int16
    with pytest.raises(KeyError):
        oxwave.gethw(t, "offset,bogus")


def test_sethw_matches_sushw_plane():
    t = oxwave.read_su(FIX / "plane.su")
    out = oxwave.sethw(t, "cdp,fldr,f2", a=[10, 1, 0.5], b=[2, 0, 0.25], c=[0, 1, 0], j=[3, 3, 3])
    assert isinstance(out, Traces) and out is not t
    assert out.to_bytes() == fixture_bytes("plane.sushw1.su")
    assert t.to_bytes() == fixture_bytes("plane.su"), "input must not be modified"


def test_sethw_matches_sushw_gather():
    t = oxwave.read_su(FIX / "gather.su")
    out = oxwave.sethw(t, ["ep", "offset", "d1"], a=[5, -100, 0.004], b=[1, 10, 0], c=[100, 0, 0],
                       d=[2, 0, 0], j=[10, 7, 0])
    assert out.to_bytes() == fixture_bytes("gather.sushw2.su")


def test_sethw_scalars_and_errors():
    t = oxwave.read_su(FIX / "plane.su")
    out = oxwave.sethw(t, "offset", a=-200, b=50)
    assert out["offset"].tolist() == list(range(-200, 200, 50))
    out = oxwave.sethw(t, "sx", a=6400, c=-100, j=3)
    assert out["sx"].tolist() == [6400, 6400, 6400, 6300, 6300, 6300, 6200, 6200]
    with pytest.raises(KeyError):
        oxwave.sethw(t, "nokey", a=1)
    with pytest.raises(ValueError):
        oxwave.sethw(t, "sx,gx", a=[1, 2, 3])


@pytest.mark.parametrize("name", SU_FIXTURES)
def test_surange_text_exact(name):
    t = oxwave.read_su(FIX / name)
    assert oxwave.surange_text(t) == (FIX / name.replace(".su", ".surange.txt")).read_text()


def test_surange_dict_plane():
    r = oxwave.surange(oxwave.read_su(FIX / "plane.su"))
    assert r["ntr"] == 8
    assert list(r["keys"]) == ["tracl", "tracr", "cdp", "offset", "sx", "gx", "ns", "dt", "d2"]
    assert r["keys"]["offset"] == (-200, 150)
    assert r["keys"]["ns"] == (64, 64)
    assert r["first"]["offset"] == -200 and r["last"]["offset"] == 150
    assert abs(r["keys"]["d2"][0] - 0.025) < 1e-9 and isinstance(r["keys"]["d2"][0], float)
    assert isinstance(r["keys"]["tracl"][0], int)
    assert r["coords"]["shot"]["East"] == (70, 0)
    assert r["coords"]["midpoint"]["East"] == (155, 0)


def test_surange_dict_one():
    r = oxwave.surange(oxwave.read_su(FIX / "one.su"))
    assert r["ntr"] == 1 and list(r["keys"]) == ["tracl", "offset", "ns", "dt"]
    assert r["coords"] == {}
