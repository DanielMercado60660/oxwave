# Contract of the Rust extension `oxwave._native` (py/src/lib.rs).
# Built with pyo3 0.30.1 (abi3-py310) on the frozen crate ../crates/oxwave-core. Every function
# must release the GIL for its byte work and must not panic (map failures to ValueError).
from typing import List, Tuple

def version() -> str:
    """The crate version, "0.1.0"."""

def header_keys() -> List[Tuple[str, str, int]]:
    """oxwave_core::header::HDR as [(key, type code, byte offset)], 80 entries in key-index
    order; equals oxwave.HEADER_KEYS."""

def bheader_keys() -> List[Tuple[str, str, int]]:
    """oxwave_core::segy::BHDR as [(key, type code, byte offset)], 27 entries; equals
    oxwave.BHEADER_KEYS."""

def ibm_to_ieee(buf: bytes, big_endian: bool = True) -> bytes:
    """IBM System/360 float words -> IEEE float32, bit-exact with segyread.c (via
    oxwave_core::segy::ibm_to_float_word). Input words are stored big-endian (or
    little-endian if big_endian=False); output is native-endian float32 bytes, same length.
    len(buf) % 4 != 0 -> ValueError."""

def ieee_to_ibm(buf: bytes, big_endian: bool = True) -> bytes:
    """Native-endian IEEE float32 bytes -> IBM words, bit-exact with segywrite.c
    (oxwave_core::segy::float_to_ibm_word), stored big-endian (or little-endian if
    big_endian=False). len(buf) % 4 != 0 -> ValueError."""

def ebcdic_to_ascii(buf: bytes) -> bytes:
    """Byte-wise `dd conv=ascii` (oxwave_core::segy::dd_conv_ascii), same length."""

def ascii_to_ebcdic(buf: bytes) -> bytes:
    """Byte-wise `dd conv=ebcdic` (oxwave_core::segy::dd_conv_ebcdic), same length."""

def zebc(buf: bytes) -> bytes:
    """Byte-wise ASCII -> EBCDIC with cwp/lib/zebc.c's table (oxwave_core::segy::zebc), the
    conversion segywrite (SUXDR build) applies to the text header; bytes outside 32..127
    become 0x20. Same length."""

def su_decode(buf: bytes) -> Tuple[bytes, bytes, int, int]:
    """Decode an SU (XDR) stream: returns (headers, samples, ns, ntr) where headers is
    ntr*240 bytes in on-disk (big-endian) layout and samples is ntr*ns native float32 bytes.
    Same acceptance rules and ValueError messages as oxwave.Traces.from_bytes."""

def su_encode(headers: bytes, samples: bytes, ns: int) -> bytes:
    """Inverse of su_decode: interleave 240-byte headers (as is) with ns big-endian float32
    samples per trace. len(headers) % 240 != 0, len(samples) != ntr*ns*4 or ns == 0 ->
    ValueError. su_encode(*su_decode(b)[:3]) == b for every valid stream b."""
