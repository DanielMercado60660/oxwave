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

//! Port of comp/dwpt/1d/lib/wpc1coding.c. Credits: CWP: Tong Chen, 1995.

use crate::{Buffer, Compressed, Quant};

/// A prefix code (`entCODE` in C).
#[derive(Clone, Copy, Debug)]
pub struct EntCode {
    pub value: u32,
    pub bits: i32,
}

/// A variable length integer codebook (`codeBOOK` in C).
#[derive(Clone, Debug)]
pub struct CodeBook {
    pub maxrun: i32,
    pub bitmaxrun: i32,
    pub code: Vec<EntCode>,
    pub pref2bit: Vec<i32>,
}

/// State of C `wpcBITBUFF`; initialize `mask` to 0x80 and `rack` to zero.
pub struct BitBuffer<'a> {
    pub buff: &'a mut Buffer,
    pub mask: u8,
    pub rack: u8,
}

fn put_byte(buff: &mut Buffer, byte: u8) -> bool {
    if buff.pos == buff.code.len() {
        return false;
    }
    buff.code[buff.pos] = byte;
    buff.pos += 1;
    true
}

fn get_byte(buff: &mut Buffer, byte: &mut u8) -> bool {
    if buff.pos == buff.code.len() {
        return false;
    }
    *byte = buff.code[buff.pos];
    buff.pos += 1;
    true
}

fn output_bits(buff: &mut BitBuffer<'_>, code: u32, count: i32) {
    // Match arm64's masked unsigned shift, including the C macro's count=0.
    let mut lmask = 1u32.wrapping_shl((count - 1) as u32);
    while lmask != 0 {
        if lmask & code != 0 {
            buff.rack |= buff.mask;
        }
        buff.mask >>= 1;
        if buff.mask == 0 {
            if !put_byte(buff.buff, buff.rack) {
                break;
            }
            buff.rack = 0;
            buff.mask = 0x80;
        }
        lmask >>= 1;
    }
}

fn input_bits_partial(buff: &mut BitBuffer<'_>, count: i32) -> (u32, bool) {
    let mut lmask = 1u32.wrapping_shl((count - 1) as u32);
    let mut number = 0;
    while lmask != 0 {
        if buff.mask == 0x80 && !get_byte(buff.buff, &mut buff.rack) {
            return (number, false);
        }
        if buff.rack & buff.mask != 0 {
            number |= lmask;
        }
        lmask >>= 1;
        buff.mask >>= 1;
        if buff.mask == 0 {
            buff.mask = 0x80;
        }
    }
    (number, true)
}

fn input_bits(buff: &mut BitBuffer<'_>, count: i32) -> Option<u32> {
    let (number, success) = input_bits_partial(buff, count);
    if success {
        Some(number)
    } else {
        None
    }
}

fn output_run(buff: &mut BitBuffer<'_>, book: &CodeBook, mut runlen: i32) {
    while runlen > 0 {
        let code = book.code[0];
        output_bits(buff, code.value, code.bits);
        if runlen <= book.maxrun {
            output_bits(buff, (runlen - 1) as u32, book.bitmaxrun);
            runlen = 0;
        } else {
            output_bits(buff, (book.maxrun - 1) as u32, book.bitmaxrun);
            runlen -= book.maxrun;
        }
    }
}

/// C prototype: `void bookCoder(void *buff, void *book, void *q, int numsmp);`
/// Leaves a partial byte in `rack` so that subsequent blocks can share it.
pub fn bookCoder(buff: &mut BitBuffer<'_>, book: &CodeBook, q: &[i32]) {
    let mut runlen = 0;
    for &qvalue in q {
        let mut symbol = qvalue;
        if symbol == 0 {
            runlen += 1;
            continue;
        }
        if runlen != 0 {
            output_run(buff, book, runlen);
            runlen = 0;
        }
        let absq = symbol.wrapping_abs();
        let mut ibit = 0;
        while ibit < 16 {
            if absq < (1 << ibit) {
                break;
            }
            ibit += 1;
        }
        let code = book.code[ibit as usize];
        output_bits(buff, code.value, code.bits);
        if symbol < 0 {
            symbol += (1 << ibit) - 1;
        }
        output_bits(buff, symbol as u32, ibit);
    }
    if runlen != 0 {
        output_run(buff, book, runlen);
    }
    let code = book.code[book.code.len() - 1];
    output_bits(buff, code.value, code.bits);
}

/// C prototype: `int bookDecoder(void *bbuff, void *book, void *q, int numsmp);`
/// Returns zero on exhausted/inconsistent input. Malformed runs are bounds checked.
pub fn bookDecoder(buff: &mut BitBuffer<'_>, book: &CodeBook, q: &mut [i32]) -> i32 {
    let code0 = book.code[0];
    let eob = book.code[book.code.len() - 1];
    let fullpref = eob.bits == code0.bits;
    let mut i = 0;
    loop {
        let Some(mut pref) = input_bits(buff, code0.bits) else {
            return 0;
        };
        if pref == code0.value {
            let Some(mut runlen) = input_bits(buff, book.bitmaxrun).map(|x| x as i32) else {
                return 0;
            };
            while runlen >= 0 {
                if i == q.len() {
                    return 0;
                }
                q[i] = 0;
                i += 1;
                runlen -= 1;
            }
        } else {
            let Some(prefix_code) = book.code.get(pref as usize) else {
                return 0;
            };
            let ibit;
            let absq;
            if prefix_code.bits == code0.bits {
                if fullpref && prefix_code.value == eob.value {
                    break;
                }
                let Some(&bits) = book.pref2bit.get(pref as usize) else {
                    return 0;
                };
                ibit = bits;
                let Some(value) = input_bits(buff, ibit) else {
                    return 0;
                };
                absq = value;
                if !fullpref && (pref.wrapping_shl(ibit as u32) + absq == eob.value) {
                    break;
                }
            } else {
                let pbit = prefix_code.bits;
                let morebit = pbit - code0.bits;
                let Some(value) = input_bits(buff, morebit) else {
                    return 0;
                };
                pref = (pref << morebit) + value;
                if pbit == eob.bits && pref == eob.value {
                    break;
                }
                let Some(&bits) = book.pref2bit.get(pref as usize) else {
                    return 0;
                };
                ibit = bits;
                let Some(value) = input_bits(buff, ibit) else {
                    return 0;
                };
                absq = value;
                if pbit + ibit == eob.bits && (pref.wrapping_shl(ibit as u32) + absq == eob.value) {
                    break;
                }
            }
            if !(1..16).contains(&ibit) || i == q.len() {
                return 0;
            }
            let symbol = if absq & (1 << (ibit - 1)) != 0 {
                absq as i32
            } else {
                absq as i32 - (1 << ibit) + 1
            };
            q[i] = symbol;
            i += 1;
        }
    }
    (i == q.len()) as i32
}

/// C prototype: `void entCoder(void *qstate, int nblock, void *outbuff);`
pub fn entCoder(quant: &[Quant], buff: &mut Buffer) {
    let mut bitbuff = BitBuffer {
        buff,
        mask: 0x80,
        rack: 0,
    };
    for q in quant {
        output_bits(&mut bitbuff, q.nbit as u32, 4);
    }
    let dict = builddict();
    for q in quant {
        bookCoder(&mut bitbuff, &dict[q.nbit as usize], &q.qx[..q.numsmp]);
    }
    if bitbuff.mask != 0x80 {
        put_byte(bitbuff.buff, bitbuff.rack);
    }
}

/// C prototype: `int entDecoder(void *qstate, int nblock, void *inbuff);`
pub fn entDecoder(quant: &mut [Quant], buff: &mut Buffer) -> i32 {
    let mut bitbuff = BitBuffer {
        buff,
        mask: 0x80,
        rack: 0,
    };
    for q in quant.iter_mut() {
        // C ignores flag here; bookDecoder subsequently detects exhaustion.
        q.nbit = input_bits_partial(&mut bitbuff, 4).0 as i32;
    }
    let dict = builddict();
    let mut retval = 0;
    for q in quant {
        let Some(book) = dict.get(q.nbit as usize) else {
            return 0;
        };
        q.qx.resize(q.numsmp, 0);
        retval = bookDecoder(&mut bitbuff, book, &mut q.qx[..q.numsmp]);
        if retval == 0 {
            break;
        }
    }
    retval
}

fn nint(x: f32) -> i32 {
    if x > 0.0 {
        (x as f64 + 0.5) as i32
    } else {
        (x as f64 - 0.5) as i32
    }
}

/// C prototype: `void wpc1Encoder(void *qstate, int nblock, void *wpc1);`
/// `w.code.len()` is the output capacity on entry and encoded size on return.
pub fn wpc1Encoder(quant: &[Quant], w: &mut Compressed) {
    let mut buff = Buffer {
        code: std::mem::take(&mut w.code),
        pos: 0,
    };
    let mut maxstep = -1.0f32;
    let mut maxave = -1.0f32;
    for q in quant {
        if maxstep < q.step {
            maxstep = q.step;
        }
        let atmp = if q.ave < 0.0 { -q.ave } else { q.ave };
        if maxave < atmp {
            maxave = atmp;
        }
    }
    w.overallstep = maxstep;
    w.overallave = maxave;
    if maxstep > f32::MIN_POSITIVE {
        let dstep = (65535.0 / maxstep as f64) as f32;
        for q in quant {
            let nstep = nint(q.step * dstep);
            put_byte(&mut buff, nstep as u8);
            put_byte(&mut buff, (nstep >> 8) as u8);
        }
    }
    if maxave > f32::MIN_POSITIVE {
        let dave = (32767.0 / maxave as f64) as f32;
        for q in quant {
            let mut nave = nint(q.ave * dave);
            if nave < 0 {
                nave += 0x10000;
            }
            put_byte(&mut buff, nave as u8);
            put_byte(&mut buff, (nave >> 8) as u8);
        }
    }
    if maxstep > f32::MIN_POSITIVE {
        entCoder(quant, &mut buff);
    }
    buff.code.truncate(buff.pos);
    w.code = buff.code;
}

/// C prototype: `int wpc1Decoder(void *qstate, int nblock, void *wpc1);`
/// As in C, averages are left untouched when `overallave <= FLT_MIN`.
pub fn wpc1Decoder(quant: &mut [Quant], w: &Compressed) -> i32 {
    let mut buff = Buffer {
        code: w.code.clone(),
        pos: 0,
    };
    let mut c1 = 0;
    let mut c2 = 0;
    if w.overallstep <= f32::MIN_POSITIVE {
        if w.overallave > f32::MIN_POSITIVE {
            let dave = (w.overallave as f64 / 32767.0) as f32;
            for q in quant.iter_mut() {
                get_byte(&mut buff, &mut c1);
                get_byte(&mut buff, &mut c2);
                let mut nave = c1 as i32 + ((c2 as i32) << 8);
                if c2 & 0x80 != 0 {
                    nave -= 0x10000;
                }
                q.ave = nave as f32 * dave;
            }
        }
        for q in quant {
            q.step = f32::MIN_POSITIVE;
        }
        1
    } else {
        let dstep = (w.overallstep as f64 / 65535.0) as f32;
        for q in quant.iter_mut() {
            get_byte(&mut buff, &mut c1);
            get_byte(&mut buff, &mut c2);
            let nstep = c1 as i32 + ((c2 as i32) << 8);
            q.step = nstep as f32 * dstep;
        }
        if w.overallave > f32::MIN_POSITIVE {
            let dave = (w.overallave as f64 / 32767.0) as f32;
            for q in quant.iter_mut() {
                get_byte(&mut buff, &mut c1);
                get_byte(&mut buff, &mut c2);
                let mut nave = c1 as i32 + ((c2 as i32) << 8);
                if c2 & 0x80 != 0 {
                    nave -= 0x10000;
                }
                q.ave = nave as f32 * dave;
            }
        }
        entDecoder(quant, &mut buff)
    }
}

// Literal dictionary from builddict(); unused prefix entries are zero.
fn builddict() -> Vec<CodeBook> {
    vec![
        CodeBook {
            maxrun: 16,
            bitmaxrun: 4,
            code: vec![
                EntCode { value: 0, bits: 1 },
                EntCode { value: 2, bits: 2 },
                EntCode { value: 3, bits: 2 },
            ],
            pref2bit: vec![0, 0, 1, 0],
        },
        CodeBook {
            maxrun: 16,
            bitmaxrun: 4,
            code: vec![
                EntCode { value: 0, bits: 2 },
                EntCode { value: 1, bits: 2 },
                EntCode { value: 2, bits: 2 },
                EntCode { value: 3, bits: 2 },
            ],
            pref2bit: vec![0, 1, 2, 0],
        },
        CodeBook {
            maxrun: 16,
            bitmaxrun: 4,
            code: vec![
                EntCode { value: 0, bits: 2 },
                EntCode { value: 1, bits: 2 },
                EntCode { value: 2, bits: 2 },
                EntCode { value: 6, bits: 3 },
                EntCode { value: 7, bits: 3 },
            ],
            pref2bit: vec![0, 1, 2, 0, 0, 0, 3, 0],
        },
        CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                EntCode { value: 0, bits: 2 },
                EntCode { value: 1, bits: 2 },
                EntCode { value: 2, bits: 2 },
                EntCode { value: 6, bits: 3 },
                EntCode { value: 7, bits: 3 },
                EntCode {
                    value: 127,
                    bits: 7,
                },
            ],
            pref2bit: vec![0, 1, 2, 0, 0, 0, 3, 4],
        },
        CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                EntCode { value: 0, bits: 2 },
                EntCode { value: 1, bits: 2 },
                EntCode { value: 4, bits: 3 },
                EntCode { value: 5, bits: 3 },
                EntCode { value: 6, bits: 3 },
                EntCode { value: 7, bits: 3 },
                EntCode {
                    value: 255,
                    bits: 8,
                },
            ],
            pref2bit: vec![0, 1, 0, 0, 2, 3, 4, 5],
        },
        CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                EntCode { value: 0, bits: 2 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 3, bits: 3 },
                EntCode { value: 4, bits: 3 },
                EntCode { value: 5, bits: 3 },
                EntCode { value: 6, bits: 3 },
                EntCode { value: 7, bits: 3 },
                EntCode {
                    value: 511,
                    bits: 9,
                },
            ],
            pref2bit: vec![0, 0, 1, 2, 3, 4, 5, 6],
        },
        CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 3, bits: 3 },
                EntCode { value: 4, bits: 3 },
                EntCode { value: 5, bits: 3 },
                EntCode { value: 6, bits: 3 },
                EntCode { value: 7, bits: 3 },
                EntCode {
                    value: 1023,
                    bits: 10,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 6, 7],
        },
        CodeBook {
            maxrun: 8,
            bitmaxrun: 3,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 3, bits: 3 },
                EntCode { value: 4, bits: 3 },
                EntCode { value: 5, bits: 3 },
                EntCode { value: 6, bits: 3 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 4095,
                    bits: 12,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 7, 8],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 3, bits: 3 },
                EntCode { value: 4, bits: 3 },
                EntCode { value: 5, bits: 3 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 8191,
                    bits: 13,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 3, bits: 3 },
                EntCode { value: 4, bits: 3 },
                EntCode { value: 10, bits: 4 },
                EntCode { value: 11, bits: 4 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 16383,
                    bits: 14,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 5, 6, 7, 8, 9, 10],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 3, bits: 3 },
                EntCode { value: 8, bits: 4 },
                EntCode { value: 9, bits: 4 },
                EntCode { value: 10, bits: 4 },
                EntCode { value: 11, bits: 4 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 32767,
                    bits: 15,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 0, 0, 0, 0, 4, 5, 6, 7, 8, 9, 10, 11],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 2, bits: 3 },
                EntCode { value: 6, bits: 4 },
                EntCode { value: 7, bits: 4 },
                EntCode { value: 8, bits: 4 },
                EntCode { value: 9, bits: 4 },
                EntCode { value: 10, bits: 4 },
                EntCode { value: 11, bits: 4 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 65535,
                    bits: 16,
                },
            ],
            pref2bit: vec![0, 1, 2, 0, 0, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 1, bits: 3 },
                EntCode { value: 4, bits: 4 },
                EntCode { value: 5, bits: 4 },
                EntCode { value: 6, bits: 4 },
                EntCode { value: 7, bits: 4 },
                EntCode { value: 8, bits: 4 },
                EntCode { value: 9, bits: 4 },
                EntCode { value: 10, bits: 4 },
                EntCode { value: 11, bits: 4 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 131071,
                    bits: 17,
                },
            ],
            pref2bit: vec![0, 1, 0, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 3 },
                EntCode { value: 2, bits: 4 },
                EntCode { value: 3, bits: 4 },
                EntCode { value: 4, bits: 4 },
                EntCode { value: 5, bits: 4 },
                EntCode { value: 6, bits: 4 },
                EntCode { value: 7, bits: 4 },
                EntCode { value: 8, bits: 4 },
                EntCode { value: 9, bits: 4 },
                EntCode { value: 10, bits: 4 },
                EntCode { value: 11, bits: 4 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 262143,
                    bits: 18,
                },
            ],
            pref2bit: vec![0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        },
        CodeBook {
            maxrun: 4,
            bitmaxrun: 2,
            code: vec![
                EntCode { value: 0, bits: 4 },
                EntCode { value: 1, bits: 4 },
                EntCode { value: 2, bits: 4 },
                EntCode { value: 3, bits: 4 },
                EntCode { value: 4, bits: 4 },
                EntCode { value: 5, bits: 4 },
                EntCode { value: 6, bits: 4 },
                EntCode { value: 7, bits: 4 },
                EntCode { value: 8, bits: 4 },
                EntCode { value: 9, bits: 4 },
                EntCode { value: 10, bits: 4 },
                EntCode { value: 11, bits: 4 },
                EntCode { value: 12, bits: 4 },
                EntCode { value: 13, bits: 4 },
                EntCode { value: 14, bits: 4 },
                EntCode { value: 15, bits: 4 },
                EntCode {
                    value: 524287,
                    bits: 19,
                },
            ],
            pref2bit: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        },
    ]
}
