// Copyright (c) 2008, Colorado School of Mines,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or
// without modification, are permitted provided that the following
// conditions are met:
//
//     *  Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     *  Redistributions in binary form must reproduce the above
//        copyright notice, this list of conditions and the following
//        disclaimer in the documentation and/or other materials provided
//        with the distribution.
//     *  Neither the name of the Colorado School of Mines nor the names of
//        its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
//
// Warranty Disclaimer:
// THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Rust port for OxWave, a Rust port of CWP/SU.

//! Python bindings for the OxWave core I/O routines.
use oxwave_core::{header, segy};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
#[pyfunction]
fn header_keys() -> Vec<(&'static str, &'static str, usize)> {
    header::HDR
        .iter()
        .map(|k| (k.key, k.type_, k.offs))
        .collect()
}
#[pyfunction]
fn bheader_keys() -> Vec<(&'static str, &'static str, usize)> {
    segy::BHDR
        .iter()
        .map(|k| (k.key, k.type_, k.offs))
        .collect()
}
#[pyfunction]
#[pyo3(signature = (buf, big_endian=true))]
fn ibm_to_ieee<'py>(
    py: Python<'py>,
    buf: &[u8],
    big_endian: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    if buf.len() % 4 != 0 {
        return Err(PyValueError::new_err(
            "buffer length must be a multiple of 4",
        ));
    }
    let out = py.detach(|| {
        let mut out = Vec::with_capacity(buf.len());
        for c in buf.chunks_exact(4) {
            let a = [c[0], c[1], c[2], c[3]];
            let word = if big_endian {
                u32::from_be_bytes(a)
            } else {
                u32::from_le_bytes(a)
            };
            out.extend_from_slice(&segy::ibm_to_float_word(word, 1).to_bits().to_ne_bytes());
        }
        out
    });
    Ok(PyBytes::new(py, &out))
}
#[pyfunction]
#[pyo3(signature = (buf, big_endian=true))]
fn ieee_to_ibm<'py>(
    py: Python<'py>,
    buf: &[u8],
    big_endian: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    if buf.len() % 4 != 0 {
        return Err(PyValueError::new_err(
            "buffer length must be a multiple of 4",
        ));
    }
    let out = py.detach(|| {
        let mut out = Vec::with_capacity(buf.len());
        for c in buf.chunks_exact(4) {
            let word = u32::from_ne_bytes([c[0], c[1], c[2], c[3]]);
            let converted = segy::float_to_ibm_word(word, 1);
            out.extend_from_slice(&if big_endian {
                converted.to_be_bytes()
            } else {
                converted.to_le_bytes()
            });
        }
        out
    });
    Ok(PyBytes::new(py, &out))
}
#[pyfunction]
fn ebcdic_to_ascii<'py>(py: Python<'py>, buf: &[u8]) -> Bound<'py, PyBytes> {
    let out = py.detach(|| {
        let mut out = buf.to_vec();
        segy::dd_conv_ascii(&mut out);
        out
    });
    PyBytes::new(py, &out)
}
#[pyfunction]
fn ascii_to_ebcdic<'py>(py: Python<'py>, buf: &[u8]) -> Bound<'py, PyBytes> {
    let out = py.detach(|| {
        let mut out = buf.to_vec();
        segy::dd_conv_ebcdic(&mut out);
        out
    });
    PyBytes::new(py, &out)
}
#[pyfunction]
fn zebc<'py>(py: Python<'py>, buf: &[u8]) -> Bound<'py, PyBytes> {
    let out = py.detach(|| {
        let mut out = vec![0; buf.len()];
        segy::zebc(buf, &mut out, buf.len());
        out
    });
    PyBytes::new(py, &out)
}
fn decode(buf: &[u8]) -> Result<(Vec<u8>, Vec<u8>, usize, usize), String> {
    const H: usize = header::HDRBYTES;
    let n = buf.len();
    if n < H {
        return Err("fgettr.c: bad first header".into());
    }
    let trid = i16::from_be_bytes([buf[28], buf[29]]);
    if trid == 201 || trid == 202 {
        return Err("packed SU traces (trid=201/202) are not supported".into());
    }
    let ns = u16::from_be_bytes([buf[114], buf[115]]) as usize;
    let nb = H + 4 * ns;
    if ns == 0 || n < nb {
        if ns > 0 && n > H {
            return Err(format!(
                "fgettr.c: first trace: tried to read {} bytes read {} bytes",
                4 * ns,
                n - H
            ));
        }
        return Err("fgettr.c: no data on first trace".into());
    }
    let ntr = n / nb;
    if n % nb >= H {
        return Err(format!("truncated SU stream: trace {} is incomplete", ntr));
    }
    for (k, rec) in buf[..ntr * nb].chunks_exact(nb).enumerate() {
        let actual = u16::from_be_bytes([rec[114], rec[115]]) as usize;
        if actual != ns {
            return Err(format!("fgettr.c: on trace #{}, number of samples in header ({}) differs from number for first trace ({})", k, actual, ns));
        }
    }
    let mut headers = Vec::with_capacity(ntr * H);
    let mut samples = Vec::with_capacity(ntr * ns * 4);
    for rec in buf[..ntr * nb].chunks_exact(nb) {
        headers.extend_from_slice(&rec[..H]);
        for c in rec[H..].chunks_exact(4) {
            samples.extend_from_slice(&u32::from_be_bytes([c[0], c[1], c[2], c[3]]).to_ne_bytes());
        }
    }
    Ok((headers, samples, ns, ntr))
}
#[pyfunction]
#[allow(clippy::type_complexity)]
fn su_decode<'py>(
    py: Python<'py>,
    buf: &[u8],
) -> PyResult<(Bound<'py, PyBytes>, Bound<'py, PyBytes>, usize, usize)> {
    let (h, s, ns, ntr) = py.detach(|| decode(buf)).map_err(PyValueError::new_err)?;
    Ok((PyBytes::new(py, &h), PyBytes::new(py, &s), ns, ntr))
}
#[pyfunction]
fn su_encode<'py>(
    py: Python<'py>,
    headers: &[u8],
    samples: &[u8],
    ns: usize,
) -> PyResult<Bound<'py, PyBytes>> {
    let result = py
        .detach(|| -> Result<Vec<u8>, &'static str> {
            const H: usize = header::HDRBYTES;
            if headers.len() % H != 0 {
                return Err("header length must be a multiple of 240");
            }
            if ns == 0 {
                return Err("cannot encode traces with ns=0");
            }
            let ntr = headers.len() / H;
            let trace_bytes = ns.checked_mul(4).ok_or("sample length overflow")?;
            let data_len = trace_bytes
                .checked_mul(ntr)
                .ok_or("sample length overflow")?;
            if samples.len() != data_len {
                return Err("sample length does not match headers and ns");
            }
            let capacity = headers
                .len()
                .checked_add(data_len)
                .ok_or("stream length overflow")?;
            let mut out = Vec::with_capacity(capacity);
            for (h, data) in headers
                .chunks_exact(H)
                .zip(samples.chunks_exact(trace_bytes))
            {
                out.extend_from_slice(h);
                for c in data.chunks_exact(4) {
                    out.extend_from_slice(
                        &u32::from_ne_bytes([c[0], c[1], c[2], c[3]]).to_be_bytes(),
                    );
                }
            }
            Ok(out)
        })
        .map_err(PyValueError::new_err)?;
    Ok(PyBytes::new(py, &result))
}
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(header_keys, m)?)?;
    m.add_function(wrap_pyfunction!(bheader_keys, m)?)?;
    m.add_function(wrap_pyfunction!(ibm_to_ieee, m)?)?;
    m.add_function(wrap_pyfunction!(ieee_to_ibm, m)?)?;
    m.add_function(wrap_pyfunction!(ebcdic_to_ascii, m)?)?;
    m.add_function(wrap_pyfunction!(ascii_to_ebcdic, m)?)?;
    m.add_function(wrap_pyfunction!(zebc, m)?)?;
    m.add_function(wrap_pyfunction!(su_decode, m)?)?;
    m.add_function(wrap_pyfunction!(su_encode, m)?)?;
    Ok(())
}
