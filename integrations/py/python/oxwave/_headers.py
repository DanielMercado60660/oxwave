# Part of OxWave, a Rust port of CWP/SU.
"""SU trace header and SEG-Y binary header layouts as NumPy dtypes.

These mirror ``segy.h`` / ``hdr.h`` / ``bhdr.h`` exactly (and the frozen Rust table
``oxwave_core::header::HDR``). Field byte order is **big-endian** (XDR), i.e. the
on-disk layout of an SU stream in this build, so a header array round-trips bit-exactly.
NumPy converts to native values transparently when you compute with a field.
"""
import numpy as np

#: hdr.h key table in SU key-index order: (key, SU type code, byte offset).
#: Type codes: "i" int32, "h" int16, "u" uint16, "f" float32.
HEADER_KEYS = (
    ("tracl", "i", 0), ("tracr", "i", 4), ("fldr", "i", 8), ("tracf", "i", 12),
    ("ep", "i", 16), ("cdp", "i", 20), ("cdpt", "i", 24), ("trid", "h", 28),
    ("nvs", "h", 30), ("nhs", "h", 32), ("duse", "h", 34), ("offset", "i", 36),
    ("gelev", "i", 40), ("selev", "i", 44), ("sdepth", "i", 48), ("gdel", "i", 52),
    ("sdel", "i", 56), ("swdep", "i", 60), ("gwdep", "i", 64), ("scalel", "h", 68),
    ("scalco", "h", 70), ("sx", "i", 72), ("sy", "i", 76), ("gx", "i", 80),
    ("gy", "i", 84), ("counit", "h", 88), ("wevel", "h", 90), ("swevel", "h", 92),
    ("sut", "h", 94), ("gut", "h", 96), ("sstat", "h", 98), ("gstat", "h", 100),
    ("tstat", "h", 102), ("laga", "h", 104), ("lagb", "h", 106), ("delrt", "h", 108),
    ("muts", "h", 110), ("mute", "h", 112), ("ns", "u", 114), ("dt", "u", 116),
    ("gain", "h", 118), ("igc", "h", 120), ("igi", "h", 122), ("corr", "h", 124),
    ("sfs", "h", 126), ("sfe", "h", 128), ("slen", "h", 130), ("styp", "h", 132),
    ("stas", "h", 134), ("stae", "h", 136), ("tatyp", "h", 138), ("afilf", "h", 140),
    ("afils", "h", 142), ("nofilf", "h", 144), ("nofils", "h", 146), ("lcf", "h", 148),
    ("hcf", "h", 150), ("lcs", "h", 152), ("hcs", "h", 154), ("year", "h", 156),
    ("day", "h", 158), ("hour", "h", 160), ("minute", "h", 162), ("sec", "h", 164),
    ("timbas", "h", 166), ("trwf", "h", 168), ("grnors", "h", 170), ("grnofr", "h", 172),
    ("grnlof", "h", 174), ("gaps", "h", 176), ("otrav", "h", 178), ("d1", "f", 180),
    ("f1", "f", 184), ("d2", "f", 188), ("f2", "f", 192), ("ungpow", "f", 196),
    ("unscale", "f", 200), ("ntr", "i", 204), ("mark", "h", 208), ("shortpad", "h", 210),
)

#: The 80 SU key names, in key-index order (``unass`` is not a key).
KEYS = tuple(k for k, _, _ in HEADER_KEYS)
#: key -> SU type code.
KEY_TYPE = {k: t for k, t, _ in HEADER_KEYS}
#: key -> byte offset in the 240-byte header.
KEY_OFFSET = {k: o for k, _, o in HEADER_KEYS}

_BE = {"i": ">i4", "h": ">i2", "u": ">u2", "f": ">f4"}
_NATIVE = {"i": np.int32, "h": np.int16, "u": np.uint16, "f": np.float32}

#: 240-byte SU trace header, big-endian fields at the segy.h offsets, plus the
#: 14 unassigned shorts ``unass`` at bytes 212..239 (kept for bit-exact round trips).
HEADER_DTYPE = np.dtype({
    "names": list(KEYS) + ["unass"],
    "formats": [_BE[t] for _, t, _ in HEADER_KEYS] + [(">i2", (14,))],
    "offsets": [o for _, _, o in HEADER_KEYS] + [212],
    "itemsize": 240,
})

HDRBYTES = 240
SU_NFLTS = 65535

#: bhdr.h binary (reel) header keys: (key, type, offset). hdt/dto/hns/nso are
#: typed "h" in SU's table although segy.h declares them unsigned short.
BHEADER_KEYS = (
    ("jobid", "i", 0), ("lino", "i", 4), ("reno", "i", 8), ("ntrpr", "h", 12),
    ("nart", "h", 14), ("hdt", "h", 16), ("dto", "h", 18), ("hns", "h", 20),
    ("nso", "h", 22), ("format", "h", 24), ("fold", "h", 26), ("tsort", "h", 28),
    ("vscode", "h", 30), ("hsfs", "h", 32), ("hsfe", "h", 34), ("hslen", "h", 36),
    ("hstyp", "h", 38), ("schn", "h", 40), ("hstas", "h", 42), ("hstae", "h", 44),
    ("htatyp", "h", 46), ("hcorr", "h", 48), ("bgrcv", "h", 50), ("rcvm", "h", 52),
    ("mfeet", "h", 54), ("polyt", "h", 56), ("vpol", "h", 58),
)
BKEYS = tuple(k for k, _, _ in BHEADER_KEYS)

#: 400-byte SEG-Y binary header (big-endian, as on tape), plus ``hunass`` (170 shorts).
BHEADER_DTYPE = np.dtype({
    "names": list(BKEYS) + ["hunass"],
    "formats": [_BE[t] for _, t, _ in BHEADER_KEYS] + [(">i2", (170,))],
    "offsets": [o for _, _, o in BHEADER_KEYS] + [60],
    "itemsize": 400,
})


def key_type(key):
    """SU type code ("i", "h", "u", "f") of header key ``key``; KeyError if unknown."""
    try:
        return KEY_TYPE[key]
    except KeyError:
        raise KeyError("unknown SU header key %r" % (key,)) from None


def native_type(key):
    """NumPy scalar type (native byte order) of header key ``key``."""
    return _NATIVE[key_type(key)]
