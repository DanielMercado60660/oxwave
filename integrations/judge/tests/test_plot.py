# Tests for (py/python/oxwave/plot.py). No binaries.
import struct

import numpy as np
import pytest

from conftest import FIX, load_su

import oxwave
from oxwave import Traces
from oxwave import plot


def png_size(path):
    b = open(path, "rb").read(32)
    assert b[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    return struct.unpack(">II", b[16:24])  # (width, height)


def ref_raster(data, c):
    x = np.nan_to_num(data.astype(np.float64), nan=0.0, posinf=np.inf, neginf=-np.inf)
    v = np.clip(x / c, -1.0, 1.0)
    return np.round(127.5 - 127.5 * v).astype(np.uint8).T


def test_clip_value():
    tr = load_su(FIX / "plane.su")
    a = np.abs(tr.data)
    assert plot.clip_value(tr) == pytest.approx(float(a.max()), rel=0, abs=0)
    dense = Traces(np.random.default_rng(0).standard_normal((8, 64)).astype(np.float32))
    ad = np.abs(dense.data)
    assert plot.clip_value(dense, perc=90) == pytest.approx(float(np.percentile(ad, 90)), rel=1e-12)
    assert plot.clip_value(tr, perc=90) == 1.0  # plane.su is mostly zeros: 90th percentile 0 -> 1.0
    assert plot.clip_value(tr, clip=-3) == 3.0
    assert plot.clip_value(Traces(np.zeros((3, 10)))) == 1.0
    with pytest.raises(ValueError):
        plot.clip_value(tr, perc=101)
    t = tr.copy()
    t.data[0, :5] = np.nan
    t.data[1, :5] = np.inf
    fin = np.abs(t.data[np.isfinite(t.data)])
    assert plot.clip_value(t) == pytest.approx(float(fin.max()))


@pytest.mark.parametrize("name,perc,clip", [("plane.su", 100.0, None), ("spike.su", 100.0, None),
                                            ("gather.su", 95.0, None), ("plane.su", 100.0, 0.25),
                                            ("dense", 97.0, None)])
def test_raster_exact(name, perc, clip):
    if name == "dense":
        tr = Traces(np.random.default_rng(1).standard_normal((30, 50)).astype(np.float32))
    else:
        tr = load_su(FIX / name)
    r = plot.raster(tr, perc=perc, clip=clip)
    assert r.dtype == np.uint8 and r.shape == (tr.ns, tr.ntr)
    c = plot.clip_value(tr, perc, clip)
    assert np.array_equal(r, ref_raster(tr.data, c))


def test_raster_semantics():
    d = np.zeros((3, 4), np.float32)
    d[0, 1] = 2.0
    d[1, 2] = -5.0
    d[2, 3] = np.nan
    r = plot.raster(Traces(d), clip=1.0)
    assert r[1, 0] == 0 and r[2, 1] == 255 and r[0, 0] == 128 and r[3, 2] == 128
    with pytest.raises(ValueError):
        plot.raster(Traces(np.zeros((0, 4), np.float32)))


@pytest.mark.parametrize("kind", ["image", "wiggle"])
def test_plot_section_png(tmp_path, kind):
    tr = load_su(FIX / "gather.su")[:48]
    p = plot.plot_section(tr, tmp_path / ("s_%s.png" % kind), kind=kind, title="gather")
    assert p == tmp_path / ("s_%s.png" % kind) and p.exists()
    assert png_size(p) == (800, 600)
    import matplotlib.image as mpimg
    img = mpimg.imread(p)
    assert img.std() > 0.01  # not blank
    p2 = plot.plot_section(tr, str(tmp_path / ("again_%s.png" % kind)), kind=kind, title="gather")
    assert p2.read_bytes() == p.read_bytes()  # deterministic


def test_plot_kinds_differ_and_sizes(tmp_path):
    tr = load_su(FIX / "spike.su")
    a = plot.plot_section(tr, tmp_path / "a.png", kind="image")
    b = plot.plot_section(tr, tmp_path / "b.png", kind="wiggle")
    assert a.read_bytes() != b.read_bytes()
    c = plot.plot_section(tr, tmp_path / "c.png", width=4, height=3, dpi=50)
    assert png_size(c) == (200, 150)
    one = plot.plot_section(load_su(FIX / "one.su"), tmp_path / "one.png", kind="wiggle")
    assert png_size(one) == (800, 600)
    nodt = Traces(np.sin(np.arange(200, dtype=np.float32) / 7).reshape(4, 50))
    assert png_size(plot.plot_section(nodt, tmp_path / "nodt.png")) == (800, 600)


def test_plot_errors(tmp_path):
    tr = load_su(FIX / "plane.su")
    with pytest.raises(ValueError):
        plot.plot_section(tr, tmp_path / "x.png", kind="contour")
    with pytest.raises(ValueError):
        plot.plot_section(Traces(np.zeros((0, 4), np.float32)), tmp_path / "y.png")
