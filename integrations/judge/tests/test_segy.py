# Tests for (py/python/oxwave/segy.py).
# Goldens: original segywrite / segyread outputs (judge/fixtures/segy, see tools/make_fixtures.py).
import io

import pytest

from conftest import FIX, fixture_bytes, load_su

import oxwave
from oxwave import SegyFile, Traces

SEG = FIX / "segy"
FILES = ["plane.sgy", "plane5.sgy", "wild.sgy"]


@pytest.fixture(autouse=True)
def _need_native(native):
    return native


@pytest.mark.parametrize("name", FILES)
def test_read_segy_traces_match_segyread(name):
    s = oxwave.read_segy(SEG / name)
    assert isinstance(s, SegyFile) and isinstance(s.traces, Traces)
    assert s.traces.to_bytes() == (SEG / (name + ".su")).read_bytes()


@pytest.mark.parametrize("name", FILES)
def test_read_segy_headers_match_segyread_files(name):
    raw = (SEG / name).read_bytes()
    s = oxwave.read_segy(str(SEG / name))
    assert s.text_header == raw[:3200]
    assert s.text_header_ascii() == (SEG / (name + ".hfile")).read_bytes()
    assert s.binary_header == (SEG / (name + ".bfile")).read_bytes()


def test_binary_header_dict():
    d = oxwave.read_segy(SEG / "plane.sgy").binary_header_dict()
    assert list(d) == list(oxwave.BKEYS)
    assert d["format"] == 1 and d["hns"] == 64 and d["hdt"] == 4000 and d["ntrpr"] == 1
    assert oxwave.read_segy(SEG / "plane5.sgy").binary_header_dict()["format"] == 5


def test_read_segy_from_bytes():
    raw = (SEG / "plane.sgy").read_bytes()
    assert oxwave.read_segy(raw).traces.to_bytes() == (SEG / "plane.sgy.su").read_bytes()


@pytest.mark.parametrize("name,fmt", [("plane.sgy", None), ("plane5.sgy", 5)])
def test_write_segy_matches_segywrite(tmp_path, name, fmt):
    tr = load_su(FIX / "plane.su")
    out = tmp_path / "out.sgy"
    oxwave.write_segy(out, tr, text_header=(SEG / "header").read_bytes(),
                      binary_header=(SEG / "binary").read_bytes(), format=fmt)
    assert out.read_bytes() == (SEG / name).read_bytes()


def test_write_segy_wild_values(tmp_path):
    tr = Traces.from_bytes((SEG / "wild.su").read_bytes())
    out = tmp_path / "w.sgy"
    oxwave.write_segy(out, tr, text_header=(SEG / "header").read_text(),
                      binary_header=(SEG / "binary").read_bytes())
    assert out.read_bytes() == (SEG / "wild.sgy").read_bytes()


def test_write_segy_file_object_and_roundtrip(tmp_path):
    tr = load_su(FIX / "plane.su")
    buf = io.BytesIO()
    oxwave.write_segy(buf, tr, text_header=(SEG / "header").read_bytes())  # default binary = segyhdrs'
    assert buf.getvalue() == (SEG / "plane.sgy").read_bytes()
    d = oxwave.read_segy(SEG / "plane.sgy").binary_header_dict()
    buf2 = io.BytesIO()
    oxwave.write_segy(buf2, tr, text_header=(SEG / "header").read_bytes(), binary_header=d)
    assert buf2.getvalue() == (SEG / "plane.sgy").read_bytes()
    p = tmp_path / "r.sgy"
    oxwave.write_segy(p, tr, format=5)
    back = oxwave.read_segy(p)
    assert back.binary_header_dict()["format"] == 5
    assert back.text_header == oxwave._native.zebc(b" " * 3200)
    # segywrite format=5 only changes the binary-header code; samples stay IBM (as the original)
    assert back.traces.to_bytes() == (SEG / "plane5.sgy.su").read_bytes()


def test_errors(tmp_path):
    with pytest.raises(ValueError):
        oxwave.read_segy(b"\x40" * 3000)
    bad = bytearray((SEG / "plane.sgy").read_bytes())
    bad[3224:3226] = (4).to_bytes(2, "big")  # format code 4 is not supported by segyread
    with pytest.raises(ValueError, match="format not SEGY standard"):
        oxwave.read_segy(bytes(bad))
    tr = load_su(FIX / "plane.su")
    with pytest.raises(ValueError):
        oxwave.write_segy(tmp_path / "x.sgy", tr, format=3)
    with pytest.raises(TypeError):
        oxwave.read_segy(SEG / "plane.sgy", endian="little")  # no byte-order option (see segy.py)
    with pytest.raises(ValueError, match="hns"):
        oxwave.write_segy(tmp_path / "x.sgy", tr, binary_header={"format": 1})
    assert fixture_bytes("plane.su")  # fixture sanity
