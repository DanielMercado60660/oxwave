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

//! OxWave core: the frozen foundation every ported SU program builds on.
//!
//! Modules (each documents its C -> Rust mapping):
//! * [`header`]: `segy.h` trace header, `hdr.h` key table, `hdrpkge.c`
//! * [`value`]:  `Value` union and `valpkge.c`
//! * [`trace`]:  `fgettr.c` / `fputtr.c` (XDR, big-endian streams)
//! * [`segy`]:   binary header, text header tables, IBM floats, byte swaps
//! * [`par`]:    `getpars.c`, `atopkge.c`
//! * [`doc`]:    `docpkge.c` (selfdoc)
//! * [`err`]:    `errpkge.c` (`err!`, `warn!`, `syserr!`)
//! * [`cio`]:    C stdio streams + `subcalls.c` (efopen, efread, ...), `filestat.c`
//! * [`cfmt`]:   C printf formatting (`printf!`, `fprintf!`, `sprintf!`)
//!
//! A program does `use oxwave_core::prelude::*;`.

pub mod cfmt;
pub mod cio;
pub mod doc;
pub mod err;
pub mod header;
pub mod par;
pub mod segy;
pub mod trace;
pub mod value;

/// cwp/lib/exit.c `CWP_Exit()`: returns EXIT_SUCCESS (0).
#[allow(non_snake_case)]
pub fn CWP_Exit() -> i32 {
    cio::EXIT_SUCCESS
}

/// su.h `CHECK_NT(label, nt)`: err if nt > SU_NFLTS.
#[allow(non_snake_case)]
pub fn CHECK_NT(label: &str, nt: i32) {
    if nt > header::SU_NFLTS as i32 {
        err!("%s=%d must not exceed %d", label, nt, header::SU_NFLTS as i32);
    }
}

pub use cio::exit;

/// Everything a ported program normally needs.
pub mod prelude {
    pub use crate::cio::{
        clearerr, efclose, efflush, efopen, efread, efseek, efseeko, eftell, eftello, efwrite,
        epclose, epopen, eremove, erename, erewind, etmpfile, exit, fclose, feof, ferror, fflush,
        fgetc, fgets, fileno, filestat, fopen, fputc, fputs, fread, fscanf_f32, fscanf_f64,
        fscanf_i32, fseeko, ftello, fwrite, isatty, putchar, stderr, stdin, stdout, ungetc, CFile,
        FileType, EOF, EXIT_FAILURE, EXIT_SUCCESS, SEEK_CUR, SEEK_END, SEEK_SET, STDIN, STDOUT,
    };
    pub use crate::doc::{pagedoc, requestdoc};
    pub use crate::header::{
        getindex, getkey, gethdval, gethval, hdtype, printheader, puthdval, puthval, swaphval,
        SegyHeader, HDR, HDRBYTES, SU_NFLTS, SU_NKEYS,
    };
    pub use crate::par::{
        args, atof, atoi, checkpars, countnparval, countparname, countparval, eatod, eatof, eatoh,
        eatoi, eatol, eatop, eatou, eatov, getnpardouble, getnparfloat, getnparint, getnparlong,
        getnparshort, getnparstring, getnparstringarray, getnparuint, getnparulong,
        getnparushort, getpardouble, getparfloat, getparint, getparlong, getparshort,
        getparstring, getparstringarray, getparuint, getparulong, getparushort, initargs,
        mustgetpardouble, mustgetparfloat, mustgetparint, mustgetparstring, strtod, strtol,
        strtoul, xargc, xargv,
    };
    pub use crate::trace::{
        fgettr, fgettra, fputtr, fvgettr, fvputtr, gettr, gettra, puttr, vgettr, vputtr, Trace,
        CHARPACK, SHORTPACK,
    };
    pub use crate::value::{
        atoval, fprintfval, getparval, printfval, printftype, valcmp, valtoabs, vtod, vtof, vtoi,
        vtol, Value,
    };
    pub use crate::{err, fprintf, printf, sprintf, syserr, warn, CHECK_NT, CWP_Exit};
}
