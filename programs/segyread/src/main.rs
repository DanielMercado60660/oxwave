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

//! SEGYREAD - read an SEG-Y tape.
//! Port of su/main/data_conversion/segyread.c (Revision 1.76).
//! Credits: SEP: Einar Kjartansson; CWP: Jack K. Cohen, Brian Sumner,
//! Chris Liner, John Stockwell (8mm tape and format 8).
//! Tony Kocurko, Department of Earth Sciences, Memorial University of
//! Newfoundland, St. John's, Newfoundland (conv and tape=-).
//! Remco Romijn and J.W. de Bruijn, Applied Geophysics, TU Delft (formats 2,3).
//! Matthias Imhof, Virginia Tech (header remapping).
//! IBM converter: SEP: Stewart A. Levin, c.1995; ibmflt.f translated by f2c
//! (version 1995/10/25). Short converter modified by Toralf Foerster,
//! Baltic Sea Research Institute, March 1997.

use oxwave_core::prelude::*;
use oxwave_core::segy::{
    gettapebhval, gettapehval, ibm_to_float, putbhval, swap_float_4, swap_int_4, swap_short_2,
    swapbhval, Bhed, BHED_NKEYS, BNYBYTES, EBCBYTES, TAPEHDR,
};
use oxwave_f_syscalls::eclose;
use std::fs::File;
use std::io::Read;
use std::os::fd::{AsFd, IntoRawFd};

/// Self-documentation, byte-identical to the C `sdoc[]` array.
#[rustfmt::skip]
const SDOC: &[&str] = &[
    "\t\t\t\t\t\t\t\t\t",
    " SEGYREAD - read an SEG-Y tape\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   segyread > stdout tape=\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   or\t\t\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "   SEG-Y data stream ... | segyread tape=-  > stdout\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Required parameter:\t\t\t\t\t\t\t",
    " tape=\t\tinput tape device or seg-y filename (see notes)\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Optional parameters:\t\t\t\t\t\t\t",
    " buff=1\tfor buffered device (9-track reel tape drive)\t\t",
    "\t\t=0 possibly useful for 8mm EXABYTE drives\t\t",
    " verbose=0\tsilent operation\t\t\t\t\t",
    "\t\t=1 ; echo every 'vblock' traces\t\t\t\t",
    " vblock=50\techo every 'vblock' traces under verbose option\t\t",
    " hfile=header\tfile to store ebcdic block (as ascii)\t\t\t",
    " bfile=binary\tfile to store binary block\t\t\t\t",
    " xfile=xhdrs\tfile to store extended text block\t\t\t",
    " over=0\tquit if bhed format not equal 1, 2, 3, 5, or 8\t\t",
    "\t\t= 1 ; override and attempt conversion\t\t\t",
    " format=bh.format\tif over=1 try to convert assuming format value  ",
    " conv=1\tconvert data to native format\t\t\t\t",
    "\t\t\t= 0 ; assume data is in native format\t\t",
    " ebcdic=1\tperform ebcdic to ascii conversion on 3200 byte textural",
    "               header. =0 do not perform conversion\t\t\t",
    " ns=bh.hns\tnumber of samples (use if bhed ns wrong)\t\t",
    " trcwt=1\tapply trace weighting factor (bytes 169-170)\t\t",
    "\t\t=0, do not apply.  (Default is 0 for formats 1 and 5)\t",
    " trmin=1\t\tfirst trace to read\t\t\t\t",
    " trmax=INT_MAX\tlast trace to read\t\t\t\t\t",
    " endian=(autodetected) =1 for big-endian,  =0 for little-endian byte order",
    " swapbhed=endian\tswap binary reel header?\t\t\t",
    " swaphdrs=endian\tswap trace headers?\t\t\t\t",
    " swapdata=endian\tswap data?\t\t\t\t\t",
    " nextended=(set from binary header) number of extended text headers\t",
    " errmax=0\tallowable number of consecutive tape IO errors\t\t",
    " remap=...,...\tremap key(s) \t\t\t\t\t\t",
    " byte=...,...\tformats to use for header remapping \t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Notes:\t\t\t\t\t\t\t\t",
    " Traditionally tape=/dev/rmt0.\t However, in the modern world tape device",
    " names are much less uniform.  The magic name can often be deduced by\t",
    " \"ls /dev\".  Likely man pages with the names of the tape devices are:",
    " \"mt\", \"sd\" \"st\".  Also try \"man -k scsi\", \" man mt\", etc.\t",
    " Sometimes \"mt status\" will tell the device name.\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " For a SEG-Y diskfile use tape=filename.\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " The xfile argument will only be used if the file contains extended\t",
    " text headers.\t\t\t\t\t\t\t\t",
    " \t\t\t\t\t\t\t\t\t",
    " Remark: a SEG-Y file is not the same as an su file. A SEG-Y file\t",
    " consists of three parts: an ebcdic header, a binary reel header, and\t",
    " the traces.  The traces are (usually) in 32 bit IBM floating point\t",
    " format.  An SU file consists only of the trace portion written in the ",
    " native binary floats.\t\t\t\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Formats supported:\t\t\t\t\t\t\t",
    " 1: IBM floating point, 4 byte (32 bits)\t\t\t\t",
    " 2: two's complement integer, 4 byte (32 bits)\t\t\t\t",
    " 3: two's complement integer, 2 byte (16 bits)\t\t\t\t",
    " 5: IEEE floating point, 4 byte (32 bits)\t\t\t\t",
    " 8: two's complement integer, 1 byte (8 bits)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " tape=-   read from standard input. Caveat, under Solaris, you will\t",
    " need to use the buff=1 option, as well.\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Header remap:\t\t\t\t\t\t\t\t",
    " The value of header word remap is mapped from the values of byte\t",
    "\t\t\t\t\t\t\t\t\t",
    " Map a float at location 221 to sample spacing d1:\t\t\t",
    "\tsegyread <data >outdata remap=d1 byte=221f\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Map a long at location 225 to source location sx:\t\t\t",
    "\tsegyread <data >outdata remap=sx byte=225l\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Map a short at location 229 to gain constant igc:\t\t\t",
    "\tsegyread <data >outdata remap=igc byte=229s\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    " Or all combined: \t\t\t\t\t\t\t",
    "\tsegyread <data >outdata remap=d1,sx,igc byte=221f,225l,229s\t",
    "\t\t\t\t\t\t\t\t\t",
    " Segy header words are accessed as Xt where X denotes the byte number\t",
    " starting at 1 in correspondance with the SEGY standard (1975)\t\t",
    " Known types include:\tf\tfloat (4 bytes)\t\t\t\t",
    " \t\t\tl\tlong int (4 bytes)\t\t\t",
    " \t\t\ts\tshort int (2 bytes)\t\t\t",
    " \t\t\tb\tbyte (1 bytes)\t\t\t\t",
    "\t\t\t\t\t\t\t\t\t",
    "\t  type:\t  sudoc segyread   for further information\t\t",
    "\t\t\t\t\t\t\t\t\t",
];

// The two SEGYREAD_TAPE_ERROR macros share this state. Descriptor reads
// consume a single tape record; stdio reads fill the requested byte count.
struct Tape {
    file: Option<File>,
    fp: Option<CFile>,
    verbose: i32,
    errmax: i32,
    errcount: i32,
}
impl Tape {
    fn read(&mut self, buf: &mut [u8], count: usize, which: &str, itr: Option<i32>) -> i32 {
        // A negative C record length becomes an invalid size_t read count.
        // Reject requests beyond the fixed tape buffer before slicing, using
        // the same error path as a failed tape read.
        let (n, failed) = if count > buf.len() {
            (-1, true)
        } else if let Some(fp) = &self.fp {
            let n = fread(&mut buf[..count], fp) as i32;
            (n, ferror(fp))
        } else {
            match self.file.as_mut().unwrap().read(&mut buf[..count]) {
                Ok(n) => (n as i32, false),
                Err(_) => (-1, true),
            }
        };
        if failed {
            if self.verbose != 0 {
                if let Some(itr) = itr {
                    warn!("tape read error on trace %d", itr);
                } else {
                    warn!("tape read error on %s", which);
                }
            }
            self.errcount = self.errcount.wrapping_add(1);
            if self.errcount > self.errmax {
                err!("exceeded maximum io errors");
            }
            if let Some(fp) = &self.fp {
                clearerr(fp);
            }
        } else {
            self.errcount = 0;
        }
        n
    }
}

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    let mut bh = Bhed::zeroed();
    let mut tapebh = [0u8; BNYBYTES];
    let mut tapetr = vec![0u8; HDRBYTES + 4 * SU_NFLTS];
    let mut ebcbuf = [0u8; EBCBYTES];
    initargs(&argv);
    requestdoc(0, SDOC); /* stdin not used */
    match filestat(STDOUT) {
        FileType::Tty => err!("stdout can't be tty"),
        FileType::Directory => err!("stdout must be a file, not a directory"),
        FileType::BadFileType => err!("stdout is illegal filetype"),
        _ => {}
    }
    /* Set filenames */
    let mut tape = String::new();
    mustgetparstring("tape", &mut tape);
    let mut hfile = String::new();
    let mut bfile = String::new();
    let mut xfile = String::new();
    if getparstring("hfile", &mut hfile) == 0 {
        hfile = "header".into();
    }
    if getparstring("bfile", &mut bfile) == 0 {
        bfile = "binary".into();
    }
    if getparstring("xfile", &mut xfile) == 0 {
        xfile = "xhdrs".into();
    }
    /* Set parameters in the original getpar order. */
    let mut trmin = 1i32;
    let mut trmax = i32::MAX;
    let mut verbose = 0i32;
    let mut vblock = 50i32;
    let mut endian = if cfg!(target_endian = "big") {
        1i32
    } else {
        0i32
    };
    getparint("trmin", &mut trmin);
    getparint("trmax", &mut trmax);
    getparint("verbose", &mut verbose);
    getparint("vblock", &mut vblock);
    getparint("endian", &mut endian);
    let mut swapdata = endian;
    let mut swapbhed = endian;
    let mut swaphdrs = endian;
    getparint("swapdata", &mut swapdata);
    getparint("swapbhed", &mut swapbhed);
    getparint("swaphdrs", &mut swaphdrs);
    let mut ebcdic = 1i32;
    let mut errmax = 0i32;
    let mut buff = 1i32;
    getparint("ebcdic", &mut ebcdic);
    getparint("errmax", &mut errmax);
    getparint("buff", &mut buff);
    let mut over = 0i32;
    let mut format = 0i32;
    getparint("over", &mut over);
    let mut format_set = getparint("format", &mut format) != 0;
    if over != 0 && format_set {
        bh.format = format as i16;
        if !matches!(format, 1 | 2 | 3 | 5 | 8) {
            warn!("Specified format=%d not supported", format);
            warn!("Assuming IBM floating point format, instead");
        }
    }
    let mut conv = 1i32;
    getparint("conv", &mut conv);
    let nkeys = countparval("remap") as usize;
    let mut key1 = Vec::<String>::new();
    let mut key2 = Vec::<String>::new();
    if nkeys != 0 {
        getparstringarray("remap", &mut key1);
    }
    let n = countparval("byte");
    if n != 0 {
        if n as usize != nkeys {
            err!("number of byte's and remap's must be equal!");
        }
        getparstringarray("byte", &mut key2);
    }
    let mut type1 = Vec::new();
    let mut index1 = Vec::new();
    for ikey in 0..nkeys {
        type1.push(hdtype(&key1[ikey]));
        index1.push(getindex(&key1[ikey]));
    }
    let mut ubyte = Vec::new();
    let mut type2 = Vec::new();
    for ikey in 0..nkeys {
        // sscanf("%d%c"): %c takes the next byte, including whitespace.
        // Missing byte= is undefined in C; use its incomplete-format error.
        let s = key2.get(ikey).map(String::as_str).unwrap_or("");
        let b = s.as_bytes();
        let mut start = 0;
        while start < b.len() && b[start].is_ascii_whitespace() {
            start += 1;
        }
        let mut end = start;
        if end < b.len() && matches!(b[end], b'+' | b'-') {
            end += 1;
        }
        let digits = end;
        while end < b.len() && b[end].is_ascii_digit() {
            end += 1;
        }
        if end == digits || end == b.len() {
            err!("user format XXXt");
        }
        ubyte.push(atoi(&s[start..end]));
        type2.push(b[end]);
    }
    /* Keep descriptor reads distinct from the fread path. */
    let mut input = Tape {
        file: None,
        fp: None,
        verbose,
        errmax,
        errcount: 0,
    };
    if buff != 0 {
        input.file = Some(if tape == "-" {
            // Core also owns stdin. Own a duplicate so dropping this File
            // cannot close the descriptor retained by its stream registry.
            match std::io::stdin().as_fd().try_clone_to_owned() {
                Ok(fd) => File::from(fd),
                Err(_) => syserr!("%s: eopen: open failed", "syscalls.c"),
            }
        } else {
            match File::open(&tape) {
                Ok(file) => file,
                Err(_) => syserr!("%s: eopen: open failed", "syscalls.c"),
            }
        });
    } else {
        input.fp = Some(if tape == "-" {
            stdin()
        } else {
            efopen(&tape, "r")
        });
    }
    let headerfp = efopen(&hfile, "w");
    if verbose != 0 {
        warn!("header file opened successfully");
    }
    let binaryfp = efopen(&bfile, "w");
    if verbose != 0 {
        warn!("binary file opened successfully");
    }
    input.read(&mut ebcbuf, EBCBYTES, "ebcdic header", None);
    let cmdbuf = if ebcdic == 1 {
        sprintf!("dd ibs=1 of=%s conv=ascii count=3200", hfile)
    } else {
        sprintf!("dd ibs=1 of=%s count=3200", hfile)
    };
    let pipefp = epopen(&cmdbuf, "w");
    efwrite(&[ebcbuf], &pipefp);
    input.read(&mut tapebh, BNYBYTES, "binary header", None);
    tapebhed_to_bhed(&tapebh, &mut bh);
    if swapbhed == 0 {
        for i in 0..BHED_NKEYS {
            swapbhval(&mut bh, i as i32);
        }
    }
    /* The source repeats the format getpars after reading the header. */
    if getparint("over", &mut over) == 0 {
        over = 0;
    }
    if getparint("format", &mut format) != 0 {
        format_set = true;
    }
    if over != 0 && format_set {
        bh.format = format as i16;
    }
    let mut trcwt = if bh.format == 1 || bh.format == 5 {
        0i32
    } else {
        1i32
    };
    getparint("trcwt", &mut trcwt);
    match bh.format {
        1 => {
            if verbose != 0 {
                warn!("assuming IBM floating point input");
            }
        }
        2 => {
            if verbose != 0 {
                warn!("assuming 4 byte integer input");
            }
        }
        3 => {
            if verbose != 0 {
                warn!("assuming 2 byte integer input");
            }
        }
        5 => {
            if verbose != 0 {
                warn!("assuming IEEE floating point input");
            }
        }
        8 => {
            if verbose != 0 {
                warn!("assuming 1 byte integer input");
            }
        }
        _ => {
            if over != 0 {
                warn!("ignoring bh.format ... continue");
            } else {
                err!("format not SEGY standard (1, 2, 3, 5, or 8)");
            }
        }
    }
    let mut ns = bh.hns as i32;
    getparint("ns", &mut ns);
    if ns == 0 {
        err!("samples/trace not set in binary header");
    }
    bh.hns = ns as u16;
    let nsegy = match bh.format {
        8 => ns.wrapping_add(HDRBYTES as i32) as usize,
        3 => ns.wrapping_mul(2).wrapping_add(HDRBYTES as i32) as usize,
        _ => ns.wrapping_mul(4).wrapping_add(HDRBYTES as i32) as usize,
    };
    efwrite(bh.as_bytes(), &binaryfp);
    efclose(&binaryfp);
    if verbose != 0 {
        warn!("binary file closed successfully");
    }
    efclose(&headerfp);
    epclose(&pipefp);
    if verbose != 0 {
        warn!("header file closed successfully");
    }
    let mut nextended = 0i16;
    // The unbraced if applies only to checkpars, not the assignment.
    if getparshort("nextended", &mut nextended) == 0 {
        checkpars();
    }
    nextended = i16::from_ne_bytes([tapebh[304], tapebh[305]]);
    if endian == 0 {
        swap_short_2(&mut nextended);
    }
    if verbose != 0 {
        warn!("Number of extended text headers: %d", nextended);
    }
    if nextended > 0 {
        let xheaderfp = efopen(&xfile, "w");
        if verbose != 0 {
            warn!("extended text header file opened successfully");
        }
        for _i in 0..nextended {
            input.read(&mut ebcbuf, EBCBYTES, "extended text header", None);
            efwrite(&ebcbuf, &xheaderfp);
        }
        efclose(&xheaderfp);
        if verbose != 0 {
            warn!("extended text header file closed successfully");
        }
    }
    let mut nsflag = false;
    let mut itr = 0i32;
    while itr < trmax {
        let nread = input.read(&mut tapetr, nsegy, "", Some(itr));
        if nread == 0 {
            break;
        }
        // Packed samples may fit the tape buffer while exceeding Trace's
        // float storage. C overruns that storage; fail before conversion.
        if ns > SU_NFLTS as i32 {
            err!("ns=%d exceeds trace capacity (%d)", ns, SU_NFLTS);
        }
        tapesegy_to_segy(&tapetr, &mut tr);
        if swaphdrs == 0 {
            for i in 0..TAPEHDR.len() {
                swaphval(&mut tr, i as i32);
            }
        }
        if !nsflag && ns != tr.ns as i32 {
            warn!(
                "discrepant tr.ns = %d with tape/user ns = %d\n\t... first noted on trace %d",
                tr.ns,
                ns,
                itr + 1
            );
            nsflag = true;
        }
        for ikey in 0..nkeys {
            let val1 = ugethval(
                type1[ikey],
                type2[ikey],
                ubyte[ikey].wrapping_sub(1),
                &tapetr,
                endian,
                conv,
                verbose,
            );
            puthval(&mut tr, index1[ikey], &val1);
        }
        itr = itr.wrapping_add(1);
        if itr >= trmin {
            if conv != 0 {
                match bh.format {
                    1 => {
                        // The C converter's signed loop is empty for ns <= 0.
                        if ns > 0 {
                            ibm_to_float(&mut tr.data, ns as usize, swapdata);
                        }
                    }
                    2 => int_to_float(&mut tr, ns, swapdata),
                    3 => short_to_float(&mut tr, ns, swapdata),
                    5 => {
                        if swapdata == 0 {
                            for i in 0..ns {
                                swap_float_4(&mut tr.data[i as usize]);
                            }
                        }
                    }
                    8 => integer1_to_float(&mut tr, ns),
                    _ => {}
                }
                if trcwt != 0 && tr.trwf != 0 {
                    let scale = 2.0f64.powf(-(tr.trwf as f64)) as f32;
                    for i in 0..ns {
                        tr.data[i as usize] *= scale;
                    }
                }
            } else {
                match bh.format {
                    1 | 5 => {
                        if swapdata == 0 {
                            for i in 0..ns {
                                swap_float_4(&mut tr.data[i as usize]);
                            }
                        }
                    }
                    2 => int_to_float(&mut tr, ns, endian),
                    3 => {
                        if swapdata == 0 {
                            // C advances by one float, swapping its first short.
                            for i in 0..ns {
                                let p = i as usize * 4;
                                let bytes = tr.data_bytes_mut();
                                let mut s = i16::from_ne_bytes([bytes[p], bytes[p + 1]]);
                                swap_short_2(&mut s);
                                bytes[p..p + 2].copy_from_slice(&s.to_ne_bytes());
                            }
                        }
                        tr.trid = SHORTPACK;
                    }
                    8 => integer1_to_float(&mut tr, ns),
                    _ => {}
                }
            }
            tr.ns = ns as u16;
            puttr(&mut tr);
            // arm64 computes remainder as itr - (itr / vblock) * vblock.
            // Division by zero yields zero, leaving itr as the remainder.
            if verbose != 0 && itr.checked_rem(vblock).unwrap_or(itr) == 0 {
                warn!(" %d traces from tape", itr);
            }
        }
    }
    if nsflag {
        warn!("discrepancy found in header and trace ns values\nthe value (%d) was used to extract traces", ns);
    }
    if let Some(fp) = &input.fp {
        efclose(fp);
    } else if let Some(file) = input.file.take() {
        eclose(file.into_raw_fd());
    }
    if verbose != 0 {
        warn!("tape closed successfully");
    }
    exit(CWP_Exit());
}

fn tapebhed_to_bhed(tapebh: &[u8], bh: &mut Bhed) {
    let mut val = Value::default();
    for i in 0..BHED_NKEYS {
        gettapebhval(tapebh, i as i32, &mut val);
        putbhval(bh, i as i32, &val);
    }
}
fn tapesegy_to_segy(tapetr: &[u8], tr: &mut Trace) {
    let mut val = Value::default();
    for i in 0..TAPEHDR.len() {
        gettapehval(tapetr, i as i32, &mut val);
        puthval(tr, i as i32, &val);
    }
    tr.as_bytes_mut()[180..240].copy_from_slice(&tapetr[180..240]);
    tr.data_bytes_mut()[..4 * SU_NFLTS].copy_from_slice(&tapetr[240..240 + 4 * SU_NFLTS]);
}
fn int_to_float(tr: &mut Trace, n: i32, endian: i32) {
    for i in 0..n {
        let mut from = tr.data[i as usize].to_bits() as i32;
        if endian == 0 {
            swap_int_4(&mut from);
        }
        tr.data[i as usize] = from as f32;
    }
}
fn short_to_float(tr: &mut Trace, n: i32, endian: i32) {
    for i in (0..n).rev() {
        let p = i as usize * 2;
        let bytes = tr.data_bytes();
        let mut from = i16::from_ne_bytes([bytes[p], bytes[p + 1]]);
        if endian == 0 {
            swap_short_2(&mut from);
        }
        tr.data[i as usize] = from as f32;
    }
}
fn integer1_to_float(tr: &mut Trace, n: i32) {
    for i in (0..n).rev() {
        tr.data[i as usize] = tr.data_bytes()[i as usize] as i8 as f32;
    }
}
fn ugethval(
    type1: &str,
    type2: u8,
    ubyte: i32,
    ptr2: &[u8],
    endian: i32,
    conv: i32,
    _verbose: i32,
) -> Value {
    let p = ubyte as usize;
    let width = match type2 {
        b'b' => 1,
        b's' => 2,
        b'l' | b'f' => 4,
        _ => 0, // Preserve the existing unsupported-type diagnostic below.
    };
    // Offsets refer to the fixed C tapesegy allocation, including its data.
    // C's unchecked pointer reads are undefined outside that allocation.
    if width != 0 && (p > ptr2.len() || width > ptr2.len() - p) {
        err!(
            "invalid remap byte offset %d for type %c",
            ubyte.wrapping_add(1),
            type2 as i32
        );
    }
    let dval1 = match type2 {
        b'b' => ptr2[p] as i8 as f64,
        b's' => {
            let mut s = i16::from_ne_bytes(ptr2[p..p + 2].try_into().unwrap());
            if endian == 0 {
                swap_short_2(&mut s);
            }
            s as f64
        }
        b'l' => {
            let mut l = i32::from_ne_bytes(ptr2[p..p + 4].try_into().unwrap());
            if endian == 0 {
                swap_int_4(&mut l);
            }
            l as f64
        }
        b'f' => {
            let mut f = f32::from_ne_bytes(ptr2[p..p + 4].try_into().unwrap());
            if conv != 0 {
                // C converts a copy in l, but never copies l back to f.
                let mut l = [f];
                ibm_to_float(&mut l, 1, endian);
            } else if endian == 0 {
                swap_float_4(&mut f);
            }
            f as f64
        }
        _ => {
            // C calls err("unknown type %s", type2): the promoted char is
            // passed where %s expects a pointer. On the reference this
            // faults after err's prologue, before printing the message.
            // Preserve the shell-visible exit status without invoking C or
            // reproducing its invalid pointer access.
            fflush(&stdout());
            fprintf!(stderr(), "\n%s: ", xargv()[0]);
            fflush(&stderr());
            exit(139);
        }
    };
    let mut val = Value::default();
    match type1.as_bytes()[0] {
        b's' => err!("can't change char header word"),
        b'h' => val.set_h(dval1 as i32 as i16),
        b'u' => val.set_u(dval1 as u32 as u16),
        b'l' => val.set_l(dval1 as i64),
        b'v' => val.set_v(dval1 as u64),
        b'i' => val.set_i(dval1 as i32),
        b'p' => val.set_p(dval1 as u32),
        b'f' => val.set_f(dval1 as f32),
        b'd' => val.set_d(dval1),
        _ => err!("unknown type %s", type1),
    }
    val
}
