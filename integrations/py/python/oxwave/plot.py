# Part of OxWave, a Rust port of CWP/SU.
"""Seismic section plots (PNG) with matplotlib (Agg backend; never opens a window).

Contract (stable signatures/docstrings). ``matplotlib`` is imported lazily inside the
functions that need it, so ``import oxwave`` works without matplotlib.
"""
import numpy as np  # noqa: F401

__all__ = ["clip_value", "raster", "plot_section"]


def clip_value(traces, perc=100.0, clip=None):
    """The amplitude that maps to full scale.

    ``clip`` given -> ``abs(float(clip))``. Otherwise the ``perc`` percentile
    (``numpy.percentile``, default linear method) of ``abs(data)`` over finite samples.
    A result that is 0, NaN or infinite (or no finite samples) -> ``1.0``.
    ``perc`` outside [0, 100] -> ValueError.
    """
    if not 0 <= perc <= 100:
        raise ValueError("perc must be between 0 and 100")
    if clip is not None:
        value = abs(float(clip))
    else:
        data = traces.data
        finite = data[np.isfinite(data)]
        value = float(np.percentile(np.abs(finite), perc)) if finite.size else 0.0
    return value if np.isfinite(value) and value != 0 else 1.0


def raster(traces, *, perc=100.0, clip=None):
    """The section as an 8-bit grey image, no axes: ``uint8`` array of shape
    ``(ns, ntr)`` (row = time sample, column = trace).

    With ``c = clip_value(traces, perc, clip)`` and ``v = clip(x / c, -1, 1)`` (NaN
    samples count as 0): ``pixel = numpy.round(127.5 - 127.5 * v)`` (round half to
    even), so ``+c`` and above -> 0 (black), ``-c`` and below -> 255 (white), 0 -> 128.
    0 traces -> ValueError.
    """
    if traces.ntr == 0:
        raise ValueError("cannot plot zero traces")
    c = clip_value(traces, perc, clip)
    data = np.asarray(traces.data, dtype=np.float64)
    data = np.where(np.isnan(data), 0.0, data)
    v = np.clip(data / c, -1.0, 1.0)
    return np.round(127.5 - 127.5 * v).astype(np.uint8).T


def plot_section(traces, path, *, kind="image", perc=99.0, clip=None, title=None,
                 width=8.0, height=6.0, dpi=100):
    """Write a PNG of the section to ``path`` and return it as a ``pathlib.Path``.

    ``kind="image"``: the :func:`raster` grey image (``imshow``, ``cmap="gray"``,
    ``vmin=0, vmax=255``, ``aspect="auto"``), time increasing downwards.
    ``kind="wiggle"``: one wiggle per trace, positive lobes filled black, scaled by
    ``clip_value``. Axes: x = trace number (1-based), y = time in seconds from ``dt`` and
    the first trace's ``delrt`` (sample index if ``dt`` is unset). ``title`` on top.
    The PNG is exactly ``round(width*dpi)`` x ``round(height*dpi)`` pixels (no tight
    bbox), written with ``metadata={"Software": None}`` so identical inputs give
    identical bytes. Any other ``kind`` -> ValueError; 0 traces -> ValueError.
    """
    from pathlib import Path

    if kind not in ("image", "wiggle"):
        raise ValueError("unknown plot kind: %s" % kind)
    if traces.ntr == 0:
        raise ValueError("cannot plot zero traces")

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    c = clip_value(traces, perc, clip)
    dt = traces.dt
    step = dt if dt is not None else 1.0
    start = int(traces.headers["delrt"][0]) * 0.001 if dt is not None else 0.0
    times = start + np.arange(traces.ns) * step
    fig = Figure(figsize=(round(width * dpi) / dpi, round(height * dpi) / dpi), dpi=dpi)
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    if kind == "image":
        ax.imshow(raster(traces, perc=perc, clip=clip), cmap="gray", vmin=0, vmax=255,
                  aspect="auto", origin="upper",
                  extent=(0.5, traces.ntr + 0.5,
                          start + (traces.ns - 0.5) * step, start - 0.5 * step))
    else:
        for i, trace in enumerate(traces.data, 1):
            values = np.asarray(trace, dtype=np.float64)
            values = np.where(np.isnan(values), 0.0, values)
            values = np.clip(values / c, -1.0, 1.0)
            ax.plot(i + values, times, color="black", linewidth=0.5)
            ax.fill_betweenx(times, i, i + values, where=values > 0, color="black")
        ax.set_xlim(0.5, traces.ntr + 0.5)
        ax.set_ylim(start + max(traces.ns - 0.5, 0.5) * step, start - 0.5 * step)
    ax.set_xlabel("Trace number")
    ax.set_ylabel("Time (s)" if dt is not None else "Sample")
    if title is not None:
        ax.set_title(title)
    path = Path(path)
    fig.savefig(path, format="png", dpi=dpi, metadata={"Software": None})
    return path
