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

//! Port of `su/lib/readkfile.c` (declarations: `su/include/readkfile.h`, revision 1.1).
//! Credits: Copyright (c) Colorado School of Mines, 2021. All rights reserved.
//! The C implementation supplies no named author credit.
//!
//! K-files contain ASCII names, printf formats, and one comma-separated K record.
//! Rust owns parsed strings; C's missing terminators and allocation overruns
//! are not reproduced. Partial outputs and error precedence follow C.

#![forbid(unsafe_code)]

use oxwave_core::cio::{fgetc, fputs, fputs_bytes, fseeko, CFile, EOF, SEEK_SET};

/// C: `void tparse(char *tbuf, char d, char **fields, int *numfields);`
///
/// Split at `d` and the first NUL, retaining leading/trailing empty fields
/// as `"null"`. The vector replaces `fields` and `numfields`.
/// Like C, this expects ASCII text with no spaces or tabs; `d` is ASCII.
pub fn tparse(tbuf: &str, d: u8) -> Vec<String> {
    assert!(d.is_ascii());
    let bytes = tbuf.as_bytes();
    let mut fields = Vec::new();
    let mut nbeg = 0;
    for n in 0..=bytes.len() {
        let c = bytes.get(n).copied().unwrap_or(0);
        if c == d || c == 0 {
            if n > nbeg {
                fields.push(tbuf[nbeg..n].to_owned());
            } else {
                fields.push("null".to_owned());
            }
            nbeg = n + 1;
        }
        if c == 0 {
            break;
        }
    }
    fields
}

/// C: `void getkvalscsv(char *textraw, char *textbeg, int maxtext,
/// char rdel, double *dfield, int *nspot, int numcases, int ncount,
/// int *comerr, int *morerr, int *numerr, int *nblank);`
///
/// `textraw.len()` replaces `maxtext`; `dfield.len()` replaces `numcases`.
/// `nspot` holds requested zero-based indices in increasing order and needs
/// at least `dfield.len()` elements. Scratch `textbeg` needs at least
/// `textraw.len()` bytes. Counters accumulate; `ncount` is unused by C.
/// The full input slice is scanned, including bytes after NUL: pass C's padded
/// buffer to retain short-record behavior. Empty output is a safe no-op
/// (C dereferences nspot[0] even when numcases is zero).
#[allow(clippy::too_many_arguments)]
pub fn getkvalscsv(
    textraw: &[u8],
    textbeg: &mut [u8],
    rdel: u8,
    dfield: &mut [f64],
    nspot: &[i32],
    _ncount: i32,
    comerr: &mut i32,
    morerr: &mut i32,
    numerr: &mut i32,
    nblank: &mut i32,
) {
    assert!(nspot.len() >= dfield.len());
    assert!(textbeg.len() >= textraw.len());
    if dfield.is_empty() {
        return;
    }
    let mut nbeg = 0;
    let mut nfield = 0i32;
    let mut ineed = 0;
    for n in 0..textraw.len() {
        let c = textraw[n];
        if c == rdel || c == 0 || c == b'\n' || c == b'\r' {
            if nfield == nspot[ineed] {
                let mut dval = 1.1e308;
                let mut nb = None;
                if n > nbeg {
                    let len = n - nbeg;
                    textbeg[..len].copy_from_slice(&textraw[nbeg..n]);
                    textbeg[len] = 0;
                    let mut ib = false;
                    for (m, &b) in textbeg[..len].iter().enumerate() {
                        if b != b' ' {
                            nb = Some(m);
                            if ib {
                                *morerr = morerr.wrapping_add(1);
                                nb = None;
                                break;
                            }
                        }
                        if b == b' ' && nb.is_some() {
                            ib = true;
                        }
                    }
                    if nb.is_some() {
                        if let Some(value) = scan_double(&textbeg[..len]) {
                            dval = value;
                        } else {
                            *numerr = numerr.wrapping_add(1);
                        }
                    }
                }
                if nb.is_none() {
                    *nblank = nblank.wrapping_add(1);
                    dval = 0.0;
                }
                dfield[ineed] = dval;
                ineed += 1;
                if ineed >= dfield.len() {
                    break;
                }
            }
            if (c == 0 || c == b'\n' || c == b'\r') && ineed < dfield.len() {
                *comerr = comerr.wrapping_add(1);
            }
            nbeg = n + 1;
            nfield = nfield.wrapping_add(1);
        }
    }
}

// fgets without the core wrapper's lossy UTF-8 conversion. Retain stale bytes
// after its terminating NUL, as C's fixed buffers do.
fn read_line(fp: &CFile, raw: &mut [u8; 10001]) -> bool {
    let mut n = 0;
    while n < 10000 {
        let c = fgetc(fp);
        if c == EOF {
            break;
        }
        raw[n] = c as u8;
        n += 1;
        if c == b'\n' as i32 {
            break;
        }
    }
    if n == 0 {
        return false;
    }
    raw[n] = 0;
    true
}

/// A Rust output representation or capacity error from [`readkfile`].
/// C's own parser errors continue to be reported through `errwarn`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadKFileError {
    /// Metadata cannot be represented by the API's UTF-8 strings.
    InvalidUtf8 {
        output: &'static str,
        source: std::str::Utf8Error,
    },
    /// An output slice is too small. Resize it to `required` and retry.
    InsufficientCapacity {
        output: &'static str,
        required: usize,
        available: usize,
    },
}

impl std::fmt::Display for ReadKFileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidUtf8 { output, source } => {
                write!(f, "readkfile: invalid UTF-8 in {output} metadata: {source}")
            }
            Self::InsufficientCapacity {
                output,
                required,
                available,
            } => write!(
                f,
                "readkfile: {output} output needs {required} elements, has {available}"
            ),
        }
    }
}

impl std::error::Error for ReadKFileError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidUtf8 { source, .. } => Some(source),
            Self::InsufficientCapacity { .. } => None,
        }
    }
}

fn check_capacity(
    output: &'static str,
    required: usize,
    available: usize,
) -> Result<(), ReadKFileError> {
    if available < required {
        return Err(ReadKFileError::InsufficientCapacity {
            output,
            required,
            available,
        });
    }
    Ok(())
}

fn parse_into(
    text: &[u8],
    fields: &mut [String],
    output: &'static str,
) -> Result<usize, ReadKFileError> {
    let text = std::str::from_utf8(text)
        .map_err(|source| ReadKFileError::InvalidUtf8 { output, source })?;
    let parsed = tparse(text, b',');
    let count = parsed.len();
    check_capacity(output, count, fields.len())?;
    for (dst, src) in fields.iter_mut().zip(parsed) {
        *dst = src;
    }
    Ok(count)
}

/// C: `void readkfile(FILE *fpR, cwp_String *names, cwp_String *forms,
/// double *dfield, int *numcasesout, int *errwarn);`
///
/// Read metadata and the first uppercase `K` record. Output slices must hold
/// all metadata fields, including C_SU_ID and null fields. The meaningful
/// prefix has `*numcasesout` elements. Errors 1..7 describe missing/duplicate
/// headers or unequal metadata counts; 9 duplicate names; 10 unreadable
/// numbers; 11 two space-separated tokens; 12 too few delimiters; -1 blank
/// fields assumed zero. Positive errors leave `numcasesout` unchanged and
/// may leave partial output changes. Code 8 is C's allocation failure;
/// Rust uses its normal allocation handling. Metadata must be UTF-8.
/// Initially indeterminate C scratch bytes are zeroed here.
///
/// Returns `Err` for metadata that cannot be represented as UTF-8 or for
/// insufficient output capacity, with the output name and required capacity.
/// These errors leave `numcasesout` unchanged, `errwarn` zero, and the stream
/// rewound, so a capacity error can be retried with larger slices. Earlier
/// output changes may remain; the failing metadata list is not copied and
/// numeric outputs are untouched. `Ok(())` includes C errors in `errwarn`.
pub fn readkfile(
    fp_r: &CFile,
    names: &mut [String],
    forms: &mut [String],
    dfield: &mut [f64],
    numcasesout: &mut i32,
    errwarn: &mut i32,
) -> Result<(), ReadKFileError> {
    *errwarn = 0;
    let mut textraw = [0u8; 10001];
    let mut textbeg = [0u8; 10001];
    let mut num_names = 0;
    let mut num_forms = 0;
    let mut num_c_su_names = 0;
    let mut num_c_su_forms = 0;
    let mut num_c_su_setid = 0;
    let mut read_names = 0;
    let mut read_forms = 0;
    while read_line(fp_r, &mut textraw) {
        // C consumes another line before noticing both lists are complete.
        if read_names == -1 && read_forms == -1 {
            break;
        }
        let mut tsize = 0;
        for &c in &textraw {
            if c == 0 || c == b'\n' || c == b'\r' {
                break;
            }
            if c != b' ' && c != b'\t' {
                textbeg[tsize] = c;
                tsize += 1;
            }
        }
        // No new NUL until a list is parsed: short lines see stale bytes.
        for b in &mut textbeg[..10001] {
            *b = b.to_ascii_lowercase();
        }
        if textbeg.starts_with(b"c_su_setid") {
            num_c_su_setid += 1;
        }
        if textbeg.starts_with(b"c_su_names") {
            num_c_su_names += 1;
        }
        if textbeg.starts_with(b"c_su_forms") {
            num_c_su_forms += 1;
        }
        if read_names > 0 {
            textbeg[tsize] = 0;
            num_names = parse_into(&textbeg[..tsize], names, "names").inspect_err(|_| {
                fseeko(fp_r, 0, SEEK_SET);
            })?;
            read_names = -1;
        }
        if textbeg.starts_with(b"c_su_names") {
            read_names = 1;
        }
        if read_forms > 0 {
            textbeg[tsize] = 0;
            num_forms = parse_into(&textbeg[..tsize], forms, "forms").inspect_err(|_| {
                fseeko(fp_r, 0, SEEK_SET);
            })?;
            read_forms = -1;
        }
        if textbeg.starts_with(b"c_su_forms") {
            read_forms = 1;
        }
    }
    fseeko(fp_r, 0, SEEK_SET);
    for (count, duplicate, missing) in [
        (num_c_su_names, 1, 2),
        (num_c_su_forms, 3, 4),
        (num_c_su_setid, 5, 6),
    ] {
        if count > 1 {
            *errwarn = duplicate;
            return Ok(());
        } else if count == 0 {
            *errwarn = missing;
            return Ok(());
        }
    }
    if num_forms != num_names {
        *errwarn = 7;
        return Ok(());
    }
    if names[0].starts_with("c_su_id") {
        names[0] = "null".to_owned();
    }
    for n in 0..num_names {
        if !names[n].starts_with("null") {
            for m in n + 1..num_names {
                if names[n] == names[m] {
                    *errwarn = 9;
                    return Ok(());
                }
            }
        }
    }
    let mut nspot = Vec::with_capacity(num_names);
    let mut numcases = 0;
    for i in 0..num_names {
        if !names[i].starts_with("null") {
            names[numcases] = names[i].clone();
            forms[numcases] = forms[i].clone();
            nspot.push(i as i32);
            numcases += 1;
        }
    }
    check_capacity("dfield", numcases, dfield.len())?;
    textraw.fill(0);
    textbeg.fill(0);
    let (mut ncount, mut comerr, mut morerr, mut numerr, mut nblank) = (0, 0, 0, 0, 0);
    let mut nextrow = false;
    while read_line(fp_r, &mut textraw) {
        ncount += 1;
        let mut textfront = [0u8; 10];
        for n in 0..10 {
            textfront[n] = textraw[n].to_ascii_lowercase();
        }
        if textfront.starts_with(b"c_su") || nextrow {
            nextrow = textfront.starts_with(b"c_su_names") || textfront.starts_with(b"c_su_forms");
        } else if textraw[0] == b'K' {
            getkvalscsv(
                &textraw,
                &mut textbeg,
                b',',
                &mut dfield[..numcases],
                &nspot,
                ncount,
                &mut comerr,
                &mut morerr,
                &mut numerr,
                &mut nblank,
            );
            break;
        }
    }
    for (count, error) in [(numerr, 10), (morerr, 11), (comerr, 12)] {
        if count > 0 {
            *errwarn = error;
            return Ok(());
        }
    }
    if nblank > 0 {
        *errwarn = -1;
    }
    *numcasesout = numcases as i32;
    Ok(())
}

/// C: `void writekfile(FILE *fpW, cwp_String *names, cwp_String *forms,
/// double *dfield, int numcasesout, int *errwarn);`
///
/// Write C_SU_SETID, forms, names, and a K record. `dfield.len()` replaces
/// `numcasesout`; names/forms must have at least that many elements.
/// Each form may consume at most one double, without `*` width/precision.
/// Values >= 1e308 and NaN become `*`; negative infinity is formatted normally.
/// As in C, `errwarn` becomes zero even when stream writes fail.
/// Rust grows buffers instead of overflowing C's 10001-byte arrays.
pub fn writekfile(
    fp_w: &CFile,
    names: &[String],
    forms: &[String],
    dfield: &[f64],
    errwarn: &mut i32,
) {
    *errwarn = 0;
    assert!(names.len() >= dfield.len() && forms.len() >= dfield.len());
    fputs("C_SU_SETID,K\n", fp_w);
    fputs("C_SU_FORMS\n", fp_w);
    for (i, fields) in [forms, names].into_iter().enumerate() {
        let mut textraw = String::from("C_SU_ID,");
        for field in &fields[..dfield.len()] {
            textraw.push_str(field.split('\0').next().unwrap());
            textraw.push(',');
        }
        textraw.pop();
        textraw.push('\n');
        fputs(&textraw, fp_w);
        if i == 0 {
            fputs("C_SU_NAMES\n", fp_w);
        }
    }
    let mut textraw2 = Vec::from(&b"K"[..]);
    for ineed in 0..dfield.len() {
        textraw2.push(b',');
        if dfield[ineed] < 1.0e308 {
            let textbeg = format_double(forms[ineed].split('\0').next().unwrap(), dfield[ineed]);
            textraw2.extend_from_slice(textbeg.split('\0').next().unwrap().as_bytes());
        } else {
            textraw2.push(b'*');
        }
    }
    textraw2.push(b'\n');
    fputs_bytes(&textraw2, fp_w);
}

// A single %lf conversion: whitespace and numeric prefixes follow strtod.
// Invalid input leaves the caller's sentinel unchanged, as sscanf does.
fn scan_double(text: &[u8]) -> Option<f64> {
    let end = text.iter().position(|&c| c == 0).unwrap_or(text.len());
    let text = String::from_utf8_lossy(&text[..end]);
    let unsigned = text.trim_start_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
    let negative = unsigned.starts_with('-');
    let unsigned = unsigned.strip_prefix(['+', '-']).unwrap_or(unsigned);
    if unsigned.starts_with("0x") || unsigned.starts_with("0X") {
        if let Some(value) = scan_hex_double(&unsigned.as_bytes()[2..]) {
            return Some(if negative { -value } else { value });
        }
    }
    let (value, consumed) = oxwave_core::par::strtod(&text);
    (consumed > 0).then_some(value)
}

// The core hex parser rounds a u64 mantissa before scaling it, which can
// round twice for subnormals. This %lf path retains guard and sticky bits
// until the final IEEE encoding instead.
fn scan_hex_double(text: &[u8]) -> Option<f64> {
    let mut digits = Vec::new();
    let mut before_dot = 0i32;
    let mut dot = false;
    let mut i = 0;
    while i < text.len() {
        if let Some(digit) = (text[i] as char).to_digit(16) {
            digits.push(digit as u8);
            if !dot {
                before_dot += 1;
            }
        } else if text[i] == b'.' && !dot {
            dot = true;
        } else {
            break;
        }
        i += 1;
    }
    if digits.is_empty() {
        return None;
    }
    let mut exponent = 0i32;
    if text.get(i).is_some_and(|&c| c == b'p' || c == b'P') {
        i += 1;
        let negative = text.get(i) == Some(&b'-');
        if text.get(i).is_some_and(|&c| c == b'+' || c == b'-') {
            i += 1;
        }
        while text.get(i).is_some_and(u8::is_ascii_digit) {
            exponent = (exponent * 10 + (text[i] - b'0') as i32).min(100_000);
            i += 1;
        }
        if negative {
            exponent = -exponent;
        }
    }
    let Some(first) = digits.iter().position(|&d| d != 0) else {
        return Some(0.0);
    };
    let top = 7 - digits[first].leading_zeros() as i32;
    exponent += 4 * (before_dot - first as i32 - 1) + top;
    if exponent > 1023 {
        return Some(f64::INFINITY);
    }
    if exponent < -1075 {
        return Some(0.0);
    }
    let keep = if exponent >= -1022 {
        53
    } else {
        exponent + 1075
    };
    let mut mant = 0u64;
    let mut count = 0;
    let mut guard = false;
    let mut sticky = false;
    for (j, &digit) in digits[first..].iter().enumerate() {
        let high = if j == 0 { top } else { 3 };
        for bit in (0..=high).rev() {
            let one = (digit >> bit) & 1 != 0;
            if count < keep {
                mant = (mant << 1) | u64::from(one);
            } else if count == keep {
                guard = one;
            } else {
                sticky |= one;
            }
            count += 1;
        }
    }
    if count < keep {
        mant <<= keep - count;
    }
    if guard && (sticky || mant & 1 != 0) {
        mant += 1;
    }
    if exponent < -1022 {
        return Some(f64::from_bits(mant));
    }
    if mant == 1u64 << 53 {
        mant >>= 1;
        exponent += 1;
    }
    let bits = ((exponent + 1023) as u64) << 52 | (mant & ((1u64 << 52) - 1));
    Some(f64::from_bits(bits))
}

fn format_double(form: &str, value: f64) -> String {
    let bytes = form.as_bytes();
    let mut out = String::new();
    let mut n = 0;
    let mut conversions = 0;
    let mut all_positional = true;
    while n < bytes.len() {
        let literal = n;
        while n < bytes.len() && bytes[n] != b'%' {
            n += 1;
        }
        out.push_str(&form[literal..n]);
        if n == bytes.len() {
            break;
        }
        n += 1;
        if bytes.get(n) == Some(&b'%') {
            out.push('%');
            n += 1;
            continue;
        }
        let start = n;
        while bytes.get(n).is_some_and(u8::is_ascii_digit) {
            n += 1;
        }
        let positional = bytes.get(n) == Some(&b'$');
        if positional {
            assert_eq!(&bytes[start..n], b"1", "K-file has only one argument");
            n += 1;
        } else {
            n = start;
        }
        assert!(
            conversions == 0 || (positional && all_positional),
            "K-file format consumes more than one double"
        );
        all_positional &= positional;
        let flags_start = n;
        while bytes.get(n).is_some_and(|b| b"-+ #0'".contains(b)) {
            n += 1;
        }
        let flags = &bytes[flags_start..n];
        let mut spec = String::from("%");
        // K-file output uses the C locale: grouping has no effect.
        for &flag in flags {
            if flag != b'\'' {
                spec.push(flag as char);
            }
        }
        let width_start = n;
        while bytes.get(n).is_some_and(u8::is_ascii_digit) {
            n += 1;
        }
        let width = form[width_start..n].parse::<usize>().unwrap_or(0);
        let mut precision = None;
        if bytes.get(n) == Some(&b'.') {
            n += 1;
            let precision_start = n;
            while bytes.get(n).is_some_and(u8::is_ascii_digit) {
                n += 1;
            }
            precision = Some(form[precision_start..n].parse::<usize>().unwrap_or(0));
        }
        spec.push_str(&form[width_start..n]);
        // printf's 'l' is a no-op for double conversions; 'L' would require
        // a long double, which the original call does not supply.
        if bytes.get(n) == Some(&b'l') {
            n += 1;
        }
        let conversion = bytes.get(n).copied().expect("incomplete K-file format");
        assert!(
            b"fFeEgGaA".contains(&conversion),
            "K-file format must accept a double"
        );
        conversions += 1;
        n += 1;
        if conversion == b'a' || conversion == b'A' {
            out.push_str(&hex_double(
                value,
                flags,
                width,
                precision,
                conversion == b'A',
            ));
        } else {
            spec.push(conversion as char);
            out.push_str(&oxwave_core::sprintf!(&spec, value));
        }
    }
    out
}

// BSD printf %a: normalize even subnormal numbers. Its rounding examines the
// first omitted nibble and (when it is 8) the parity of the following nibble,
// rather than that of the retained digit. Preserve this reference quirk.
fn hex_double(
    value: f64,
    flags: &[u8],
    width: usize,
    precision: Option<usize>,
    upper: bool,
) -> String {
    let sign = if value.is_nan() {
        ""
    } else if value.is_sign_negative() {
        "-"
    } else if flags.contains(&b'+') {
        "+"
    } else if flags.contains(&b' ') {
        " "
    } else {
        ""
    };
    let mut prefix = "";
    let mut body;
    if !value.is_finite() {
        body = if value.is_nan() { "nan" } else { "inf" }.to_owned();
    } else {
        prefix = if upper { "0X" } else { "0x" };
        let bits = value.abs().to_bits();
        let mut mant = bits & ((1u64 << 52) - 1);
        let mut exponent = ((bits >> 52) & 2047) as i32 - 1023;
        if bits == 0 {
            exponent = 0;
        } else if bits >> 52 == 0 {
            let shift = mant.leading_zeros() as i32 - 11;
            mant <<= shift;
            exponent = -1022 - shift;
        } else {
            mant |= 1u64 << 52;
        }
        let p = precision.unwrap_or(13);
        if p < 13 {
            let shift = 52 - 4 * p;
            let tail = mant & ((1u64 << shift) - 1);
            let mut retained = mant >> shift;
            let first = tail >> (shift - 4);
            let following_odd = shift >= 8 && ((tail >> (shift - 8)) & 1) != 0;
            let up = first > 8 || (first == 8 && following_odd);
            if up {
                retained += 1;
            }
            mant = retained << shift;
        }
        body = format!("{:x}", mant >> 52);
        let mut fraction = format!("{:013x}", mant & ((1u64 << 52) - 1));
        if precision.is_none() {
            while fraction.ends_with('0') {
                fraction.pop();
            }
        } else if p <= 13 {
            fraction.truncate(p);
        } else {
            fraction.extend(std::iter::repeat_n('0', p - 13));
        }
        if !fraction.is_empty() || flags.contains(&b'#') {
            body.push('.');
            body.push_str(&fraction);
        }
        body.push_str(&format!("p{exponent:+}"));
    }
    if upper {
        body.make_ascii_uppercase();
    }
    let pad = width.saturating_sub(sign.len() + prefix.len() + body.len());
    let spaces = " ".repeat(pad);
    if flags.contains(&b'-') {
        format!("{sign}{prefix}{body}{spaces}")
    } else if flags.contains(&b'0') && value.is_finite() {
        format!("{sign}{prefix}{}{body}", "0".repeat(pad))
    } else {
        format!("{spaces}{sign}{prefix}{body}")
    }
}
