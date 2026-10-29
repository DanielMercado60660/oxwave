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

//! `cio`: C stdio streams (`FILE *`) and SU's error-checking wrappers from
//! `par/lib/subcalls.c` (efopen, efread, ...), plus `filestat`.
//!
//! `CFile` is a cheap-to-clone handle to one C stream: clones are the *same*
//! `FILE *` (shared buffer, position and SU trace I/O state).  `stdin()`,
//! `stdout()`, `stderr()` return the process-wide streams on fds 0/1/2.
//! stdout is fully buffered (line buffered on a tty) and stderr unbuffered,
//! like C.  **Buffered data is only flushed by `fflush`/`efflush`/`efclose`,
//! by dropping the last handle of a file, or by [`exit`]** -- every program
//! must end with `oxwave_core::exit(CWP_Exit())` (or `exit(code)`), never by
//! returning from `main`.
//!
//! | C                                  | Rust                                    |
//! |------------------------------------|-----------------------------------------|
//! | `stdin`, `stdout`, `stderr`        | `stdin()`, `stdout()`, `stderr()`       |
//! | `efopen(name, mode)`               | `efopen(name, mode) -> CFile`           |
//! | `efclose(fp)`                      | `efclose(&fp)`                          |
//! | `etmpfile()`                       | `etmpfile()`                            |
//! | `epopen(cmd, type)` / `epclose(fp)`| `epopen(cmd, type)` / `epclose(&fp)`    |
//! | `efread(buf, size, n, fp)`         | `efread(&mut buf[..n], &fp) -> usize` (items, any `Pod` element) |
//! | `efwrite(buf, size, n, fp)`        | `efwrite(&buf[..n], &fp) -> usize`      |
//! | `fread`/`fwrite` (no check)        | `fread` / `fwrite`                      |
//! | `efseeko(fp, off, SEEK_SET)`       | `efseeko(&fp, off, SEEK_SET)`           |
//! | `efseek`, `eftell(o)`, `erewind`   | `efseek`, `eftell`/`eftello`, `erewind` |
//! | `efflush(fp)` / `fflush(fp)`       | `efflush(&fp)` / `fflush(&fp)`          |
//! | `fgets(buf, n, fp)`                | `fgets(n, &fp) -> Option<String>`       |
//! | `fputs(s, fp)`, `fputc`, `putchar` | `fputs(s, &fp)`, `fputc(c, &fp)`, `putchar(c)` |
//! | `getc`/`fgetc`, `ungetc`           | `fgetc(&fp) -> i32` (EOF = -1), `ungetc(c, &fp)` |
//! | `feof`, `ferror`, `clearerr`       | `feof(&fp)`, `ferror(&fp)`, `clearerr(&fp)` |
//! | `fscanf(fp, "%f"/"%lf"/"%d", &x)`  | `fscanf_f32/fscanf_f64/fscanf_i32(&fp) -> Option<_>` (None = EOF or no match) |
//! | `fileno(fp)`, `isatty(fd)`         | `fileno(&fp)`, `isatty(fd)`             |
//! | `filestat(fd)`                     | `filestat(fd) -> FileType`              |
//! | `eremove(f)`, `erename(a, b)`      | `eremove(f)`, `erename(a, b)`           |
//! | `exit(code)`                       | `exit(code)` (flushes all streams)      |

use std::cell::RefCell;
use std::fs::{File, OpenOptions};
use std::io::{IsTerminal, Read, Seek, SeekFrom, Write};
use std::os::fd::{AsRawFd, FromRawFd, IntoRawFd, OwnedFd};
use std::os::unix::fs::{FileTypeExt, MetadataExt};
use std::rc::{Rc, Weak};

pub const SEEK_SET: i32 = 0;
pub const SEEK_CUR: i32 = 1;
pub const SEEK_END: i32 = 2;
pub const EOF: i32 = -1;
pub const STDIN: i32 = 0;
pub const STDOUT: i32 = 1;
pub const STDERR: i32 = 2;
pub const EXIT_SUCCESS: i32 = 0;
pub const EXIT_FAILURE: i32 = 1;

/// cwp.h `FileType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum FileType {
    BadFileType = -1,
    Tty = 0,
    Disk = 1,
    Directory = 2,
    Tape = 3,
    Pipe = 4,
    Fifo = 5,
    Socket = 6,
    Symlink = 7,
}

/// Plain-old-data element types that `efread`/`efwrite` move as raw native
/// (little-endian) bytes, like C `fread`/`fwrite`.
///
/// # Safety
/// Implementors must be `Copy`, have no padding, and accept every bit pattern.
pub unsafe trait Pod: Copy + 'static {}
unsafe impl Pod for u8 {}
unsafe impl Pod for i8 {}
unsafe impl Pod for u16 {}
unsafe impl Pod for i16 {}
unsafe impl Pod for u32 {}
unsafe impl Pod for i32 {}
unsafe impl Pod for u64 {}
unsafe impl Pod for i64 {}
unsafe impl Pod for f32 {}
unsafe impl Pod for f64 {}
unsafe impl Pod for crate::header::SegyHeader {}
unsafe impl<T: Pod, const N: usize> Pod for [T; N] {}

/// View a slice of `Pod` values as bytes.
pub fn as_bytes<T: Pod>(v: &[T]) -> &[u8] {
    // SAFETY: Pod guarantees no padding.
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}
/// View a mutable slice of `Pod` values as bytes.
pub fn as_bytes_mut<T: Pod>(v: &mut [T]) -> &mut [u8] {
    // SAFETY: Pod guarantees every bit pattern is valid.
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, std::mem::size_of_val(v)) }
}

#[derive(PartialEq, Clone, Copy)]
enum Buffering {
    Full,
    Line,
    None,
}

/// Per-stream input state of fgettr (struct insegyinfo).
#[derive(Clone, Debug)]
pub(crate) struct TrIn {
    pub itr: i64,
    pub nsfirst: i32,
    pub bytesper: u16,
    pub nsegy: usize,
    pub ntr: i64,
    pub ftype: FileType,
}

/// Per-stream output state of fputtr (struct outsegyinfo).
#[derive(Clone, Debug)]
pub(crate) struct TrOut {
    pub itr: i64,
    pub nsfirst: u32,
    pub bytesper: u16,
}

pub(crate) struct Inner {
    file: Option<File>,
    fd: i32,
    owned: bool,
    rbuf: Vec<u8>,
    rpos: usize,
    wbuf: Vec<u8>,
    buffering: Buffering,
    eof: bool,
    error: bool,
    child: Option<std::process::Child>,
    pub(crate) trin: Option<TrIn>,
    pub(crate) trout: Option<TrOut>,
}

const BUFSZ: usize = 65536;

impl Inner {
    fn new(file: File, owned: bool, buffering: Buffering) -> Inner {
        let fd = file.as_raw_fd();
        Inner {
            file: Some(file),
            fd,
            owned,
            rbuf: Vec::new(),
            rpos: 0,
            wbuf: Vec::new(),
            buffering,
            eof: false,
            error: false,
            child: None,
            trin: None,
            trout: None,
        }
    }

    fn flush(&mut self) -> bool {
        if self.wbuf.is_empty() {
            return true;
        }
        let ok = match self.file.as_mut() {
            Some(f) => f.write_all(&self.wbuf).is_ok(),
            None => false,
        };
        self.wbuf.clear();
        if !ok {
            self.error = true;
        }
        ok
    }

    /// Drop read-ahead, repositioning the fd to the logical position.
    fn drop_readahead(&mut self) {
        let unread = self.rbuf.len() - self.rpos;
        if unread > 0 {
            if let Some(f) = self.file.as_mut() {
                let _ = f.seek(SeekFrom::Current(-(unread as i64)));
            }
        }
        self.rbuf.clear();
        self.rpos = 0;
    }

    fn read(&mut self, out: &mut [u8]) -> usize {
        if !self.wbuf.is_empty() {
            self.flush();
        }
        let mut n = 0;
        while n < out.len() {
            if self.rpos < self.rbuf.len() {
                let k = (self.rbuf.len() - self.rpos).min(out.len() - n);
                out[n..n + k].copy_from_slice(&self.rbuf[self.rpos..self.rpos + k]);
                self.rpos += k;
                n += k;
                continue;
            }
            let want = out.len() - n;
            let f = match self.file.as_mut() {
                Some(f) => f,
                None => {
                    self.error = true;
                    break;
                }
            };
            if want >= BUFSZ {
                match f.read(&mut out[n..]) {
                    Ok(0) => {
                        self.eof = true;
                        break;
                    }
                    Ok(k) => n += k,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                    Err(_) => {
                        self.error = true;
                        break;
                    }
                }
            } else {
                self.rbuf.resize(BUFSZ, 0);
                self.rpos = 0;
                match f.read(&mut self.rbuf[..]) {
                    Ok(0) => {
                        self.rbuf.clear();
                        self.eof = true;
                        break;
                    }
                    Ok(k) => self.rbuf.truncate(k),
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => self.rbuf.clear(),
                    Err(_) => {
                        self.rbuf.clear();
                        self.error = true;
                        break;
                    }
                }
            }
        }
        n
    }

    fn write(&mut self, data: &[u8]) -> usize {
        if self.rpos < self.rbuf.len() || !self.rbuf.is_empty() {
            self.drop_readahead();
        }
        if self.file.is_none() {
            self.error = true;
            return 0;
        }
        match self.buffering {
            Buffering::None => {
                let ok = self.file.as_mut().unwrap().write_all(data).is_ok();
                if !ok {
                    self.error = true;
                    return 0;
                }
            }
            Buffering::Line => {
                self.wbuf.extend_from_slice(data);
                if data.contains(&b'\n') || self.wbuf.len() >= BUFSZ {
                    if !self.flush() {
                        return 0;
                    }
                }
            }
            Buffering::Full => {
                self.wbuf.extend_from_slice(data);
                if self.wbuf.len() >= BUFSZ && !self.flush() {
                    return 0;
                }
            }
        }
        data.len()
    }

    fn seek(&mut self, off: i64, whence: i32) -> bool {
        if !self.flush() {
            return false;
        }
        let unread = (self.rbuf.len() - self.rpos) as i64;
        self.rbuf.clear();
        self.rpos = 0;
        let f = match self.file.as_mut() {
            Some(f) => f,
            None => return false,
        };
        let pos = match whence {
            SEEK_SET => {
                if off < 0 {
                    return false;
                }
                SeekFrom::Start(off as u64)
            }
            SEEK_CUR => SeekFrom::Current(off - unread),
            SEEK_END => SeekFrom::End(off),
            _ => return false,
        };
        // If seeking fails (pipe), restore nothing: C leaves the stream as is.
        if f.seek(pos).is_err() {
            return false;
        }
        self.eof = false;
        true
    }

    fn tell(&mut self) -> i64 {
        let unread = (self.rbuf.len() - self.rpos) as i64;
        let pend = self.wbuf.len() as i64;
        match self.file.as_mut().map(|f| f.stream_position()) {
            Some(Ok(p)) => p as i64 - unread + pend,
            _ => -1,
        }
    }

    fn close(&mut self) -> bool {
        let ok = self.flush();
        if let Some(f) = self.file.take() {
            if !self.owned {
                let _ = f.into_raw_fd();
            }
        }
        if let Some(mut c) = self.child.take() {
            let _ = c.wait();
        }
        ok
    }
}

impl Drop for Inner {
    fn drop(&mut self) {
        self.close();
    }
}

/// A C stream (`FILE *`). Clones refer to the same stream.
#[derive(Clone)]
pub struct CFile(pub(crate) Rc<RefCell<Inner>>);

impl PartialEq for CFile {
    fn eq(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.0, &other.0)
    }
}

impl std::fmt::Debug for CFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CFile(fd {})", self.0.borrow().fd)
    }
}

thread_local! {
    static REGISTRY: RefCell<Vec<Weak<RefCell<Inner>>>> = const { RefCell::new(Vec::new()) };
    static STD: (CFile, CFile, CFile) = {
        // SAFETY: fds 0/1/2 are owned by the process; `owned == false`
        // ensures we never close them.
        let mk = |fd: i32, b: Buffering| CFile::register(Inner::new(unsafe { File::from_raw_fd(fd) }, false, b));
        let out_buf = if std::io::stdout().is_terminal() { Buffering::Line } else { Buffering::Full };
        (mk(0, Buffering::Full), mk(1, out_buf), mk(2, Buffering::None))
    };
}

impl CFile {
    fn register(inner: Inner) -> CFile {
        let rc = Rc::new(RefCell::new(inner));
        REGISTRY.with(|r| {
            let mut r = r.borrow_mut();
            r.retain(|w| w.strong_count() > 0);
            r.push(Rc::downgrade(&rc));
        });
        CFile(rc)
    }

    /// Wrap an already-open `File` as a C stream (fully buffered, closed on drop).
    pub fn from_file(f: File) -> CFile {
        CFile::register(Inner::new(f, true, Buffering::Full))
    }

    pub(crate) fn with<R>(&self, f: impl FnOnce(&mut Inner) -> R) -> R {
        f(&mut self.0.borrow_mut())
    }
}

/// The process's standard input stream.
pub fn stdin() -> CFile {
    STD.with(|s| s.0.clone())
}
/// The process's standard output stream (buffered; see module docs).
pub fn stdout() -> CFile {
    STD.with(|s| s.1.clone())
}
/// The process's standard error stream (unbuffered).
pub fn stderr() -> CFile {
    STD.with(|s| s.2.clone())
}

/// Flush every open stream (C `fflush(NULL)`).
pub fn flush_all() {
    let list: Vec<Rc<RefCell<Inner>>> =
        REGISTRY.with(|r| r.borrow().iter().filter_map(|w| w.upgrade()).collect());
    for rc in list {
        if let Ok(mut i) = rc.try_borrow_mut() {
            i.flush();
        }
    }
}

/// C `exit(code)`: flush all stdio streams, then terminate.
pub fn exit(code: i32) -> ! {
    flush_all();
    std::process::exit(code)
}

/// Restore default SIGPIPE handling (C programs die silently when their
/// reader goes away; Rust ignores SIGPIPE by default).  Called by `initargs`.
pub fn restore_sigpipe() {
    extern "C" {
        fn signal(sig: i32, handler: usize) -> usize;
    }
    // SAFETY: SIGPIPE (13 on macOS/Linux) set to SIG_DFL (0).
    unsafe {
        signal(13, 0);
    }
}

fn mode_options(mode: &str) -> Option<(OpenOptions, Buffering)> {
    let m: String = mode.chars().filter(|&c| c != 'b').collect();
    let mut o = OpenOptions::new();
    match m.as_str() {
        "r" => o.read(true),
        "w" => o.write(true).create(true).truncate(true),
        "a" => o.append(true).create(true),
        "r+" => o.read(true).write(true),
        "w+" => o.read(true).write(true).create(true).truncate(true),
        "a+" => o.read(true).append(true).create(true),
        _ => return None,
    };
    Some((o, Buffering::Full))
}

/// fopen without error check (None = NULL).
pub fn fopen(file: &str, mode: &str) -> Option<CFile> {
    let (o, b) = mode_options(mode)?;
    let f = o.open(file).ok()?;
    Some(CFile::register(Inner::new(f, true, b)))
}

/// efopen: fopen with error check ("subcalls.c: efopen: fopen failed").
pub fn efopen(file: &str, mode: &str) -> CFile {
    match fopen(file, mode) {
        Some(f) => f,
        None => crate::err!("%s: efopen: fopen failed", "subcalls.c"),
    }
}

/// fclose; returns 0 or EOF.
pub fn fclose(fp: &CFile) -> i32 {
    if fp.with(|i| i.close()) {
        0
    } else {
        EOF
    }
}

/// efclose: fclose with error check.
pub fn efclose(fp: &CFile) -> i32 {
    if fclose(fp) == EOF {
        crate::err!("%s: efclose: fclose failed", "subcalls.c");
    }
    0
}

/// etmpfile: anonymous temporary file opened "w+" and already unlinked.
pub fn etmpfile() -> CFile {
    use std::sync::atomic::{AtomicU64, Ordering};
    static N: AtomicU64 = AtomicU64::new(0);
    let dir = std::env::temp_dir();
    for _ in 0..100 {
        let n = N.fetch_add(1, Ordering::Relaxed);
        let t = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let p = dir.join(format!("oxwave_tmp_{}_{}_{}", std::process::id(), n, t));
        if let Ok(f) = OpenOptions::new().read(true).write(true).create_new(true).open(&p) {
            let _ = std::fs::remove_file(&p);
            return CFile::register(Inner::new(f, true, Buffering::Full));
        }
    }
    crate::err!("%s: etmpfile: tmpfile failed", "subcalls.c")
}

/// epopen: popen via `/bin/sh -c` with error check; `type_` is "r" or "w".
pub fn epopen(command: &str, type_: &str) -> CFile {
    use std::process::{Command, Stdio};
    flush_all();
    let mut cmd = Command::new("/bin/sh");
    cmd.arg("-c").arg(command);
    let reading = type_.starts_with('r');
    if reading {
        cmd.stdout(Stdio::piped());
    } else {
        cmd.stdin(Stdio::piped());
    }
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        Err(_) => crate::err!("%s: epopen: popen failed", "subcalls.c"),
    };
    let fd: OwnedFd = if reading {
        child.stdout.take().unwrap().into()
    } else {
        child.stdin.take().unwrap().into()
    };
    let mut inner = Inner::new(File::from(fd), true, Buffering::Full);
    inner.child = Some(child);
    CFile::register(inner)
}

/// epclose: pclose with error check; returns the child's wait status.
pub fn epclose(fp: &CFile) -> i32 {
    let st = fp.with(|i| {
        i.flush();
        if let Some(f) = i.file.take() {
            drop(f);
        }
        i.child.take().map(|mut c| c.wait())
    });
    match st {
        Some(Ok(s)) => {
            use std::os::unix::process::ExitStatusExt;
            s.into_raw()
        }
        _ => crate::err!("%s: epclose: pclose failed", "subcalls.c"),
    }
}

/// fread without error check: reads `buf.len()` items, returns complete items read.
pub fn fread<T: Pod>(buf: &mut [T], fp: &CFile) -> usize {
    let sz = std::mem::size_of::<T>();
    let n = fp.with(|i| i.read(as_bytes_mut(buf)));
    n / sz
}

/// efread: fread with error check (error only on a stream error, not EOF).
pub fn efread<T: Pod>(buf: &mut [T], fp: &CFile) -> usize {
    let count = buf.len();
    let nread = fread(buf, fp);
    if nread != count && ferror(fp) {
        crate::err!("%s: efread: fread only %d items of %d", "subcalls.c", nread, count);
    }
    nread
}

/// fwrite without error check: returns complete items written.
pub fn fwrite<T: Pod>(buf: &[T], fp: &CFile) -> usize {
    let sz = std::mem::size_of::<T>();
    let n = fp.with(|i| i.write(as_bytes(buf)));
    n / sz
}

/// efwrite: fwrite with error check.
pub fn efwrite<T: Pod>(buf: &[T], fp: &CFile) -> usize {
    let count = buf.len();
    let nwrite = fwrite(buf, fp);
    if nwrite != count {
        crate::err!("%s: efwrite: fwrite only %d items of %d", "subcalls.c", nwrite, count);
    }
    nwrite
}

/// fseeko; returns 0 on success, -1 on failure.
pub fn fseeko(fp: &CFile, offset: i64, origin: i32) -> i32 {
    if fp.with(|i| i.seek(offset, origin)) {
        0
    } else {
        -1
    }
}

/// efseeko: fseeko with error check ("subcalls.c : efseeko: fseeko failed").
pub fn efseeko(fp: &CFile, offset: i64, origin: i32) -> i32 {
    if fseeko(fp, offset, origin) != 0 {
        crate::err!("%s : efseeko: fseeko failed", "subcalls.c");
    }
    0
}

/// efseek: fseek with error check ("subcalls.c: efseek: fseek failed").
pub fn efseek(fp: &CFile, offset: i64, origin: i32) -> i32 {
    if fseeko(fp, offset, origin) != 0 {
        crate::err!("%s: efseek: fseek failed", "subcalls.c");
    }
    0
}

/// ftello (-1 on a pipe).
pub fn ftello(fp: &CFile) -> i64 {
    fp.with(|i| i.tell())
}

/// eftello: ftello; on failure prints "sizeof(off_t)=8" to stderr like SU.
pub fn eftello(fp: &CFile) -> i64 {
    let p = ftello(fp);
    if p == -1 {
        fputs("sizeof(off_t)=8\n", &stderr());
    }
    p
}

/// eftell: ftell with error check.
pub fn eftell(fp: &CFile) -> i64 {
    let p = ftello(fp);
    if p == -1 {
        crate::syserr!("%s: eftell: ftell failed", "subcalls.c");
    }
    p
}

/// erewind / rewind: seek to 0 and clear the error and EOF indicators.
pub fn erewind(fp: &CFile) {
    fp.with(|i| {
        i.seek(0, SEEK_SET);
        i.error = false;
        i.eof = false;
    });
}
pub use erewind as rewind;

/// fflush; returns 0 or EOF.
pub fn fflush(fp: &CFile) -> i32 {
    if fp.with(|i| i.flush()) {
        0
    } else {
        EOF
    }
}

/// efflush: fflush with error check.
pub fn efflush(fp: &CFile) -> i32 {
    if fflush(fp) == EOF {
        crate::err!("%s: efflush: fflush failed", "subcalls.c");
    }
    0
}

pub fn feof(fp: &CFile) -> bool {
    fp.with(|i| i.eof)
}
pub fn ferror(fp: &CFile) -> bool {
    fp.with(|i| i.error)
}
pub fn clearerr(fp: &CFile) {
    fp.with(|i| {
        i.eof = false;
        i.error = false
    })
}
pub fn fileno(fp: &CFile) -> i32 {
    fp.with(|i| i.fd)
}

/// fputs (also used by printf!/fprintf!). Returns 0 or EOF.
pub fn fputs(s: &str, fp: &CFile) -> i32 {
    fputs_bytes(s.as_bytes(), fp)
}

/// fputs for raw bytes.
pub fn fputs_bytes(s: &[u8], fp: &CFile) -> i32 {
    if fp.with(|i| i.write(s)) == s.len() {
        0
    } else {
        EOF
    }
}

pub fn fputc(c: i32, fp: &CFile) -> i32 {
    if fp.with(|i| i.write(&[c as u8])) == 1 {
        c & 0xFF
    } else {
        EOF
    }
}
/// putchar(c) on stdout.
pub fn putchar(c: i32) -> i32 {
    fputc(c, &stdout())
}

/// fgetc / getc: next byte or EOF (-1).
pub fn fgetc(fp: &CFile) -> i32 {
    let mut b = [0u8; 1];
    if fp.with(|i| i.read(&mut b)) == 1 {
        b[0] as i32
    } else {
        EOF
    }
}

/// ungetc: push one byte back.
pub fn ungetc(c: i32, fp: &CFile) -> i32 {
    if c == EOF {
        return EOF;
    }
    fp.with(|i| {
        if i.rpos > 0 {
            i.rpos -= 1;
            i.rbuf[i.rpos] = c as u8;
        } else {
            i.rbuf.insert(0, c as u8);
        }
        i.eof = false;
    });
    c
}

/// fgets: read at most n-1 bytes, stopping after a newline. None at EOF
/// with nothing read (C returns NULL). Bytes are converted lossily to UTF-8.
pub fn fgets(n: usize, fp: &CFile) -> Option<String> {
    let mut v = Vec::new();
    while v.len() + 1 < n {
        let c = fgetc(fp);
        if c == EOF {
            break;
        }
        v.push(c as u8);
        if c == b'\n' as i32 {
            break;
        }
    }
    if v.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(&v).into_owned())
    }
}

fn scan_token(fp: &CFile, float: bool) -> Option<String> {
    // skip white space
    let mut c = fgetc(fp);
    while c != EOF && (c as u8).is_ascii_whitespace() {
        c = fgetc(fp);
    }
    if c == EOF {
        return None;
    }
    let mut t = String::new();
    let accept = |t: &str, ch: u8| -> bool {
        let last = t.as_bytes().last().copied();
        match ch {
            b'+' | b'-' => t.is_empty() || (float && matches!(last, Some(b'e') | Some(b'E'))),
            b'0'..=b'9' => true,
            b'.' => float && !t.contains('.') && !t.contains(['e', 'E']),
            b'e' | b'E' => float && !t.contains(['e', 'E']) && t.bytes().any(|b| b.is_ascii_digit()),
            _ => false,
        }
    };
    while c != EOF && accept(&t, c as u8) {
        t.push(c as u8 as char);
        c = fgetc(fp);
    }
    if c != EOF {
        ungetc(c, fp);
    }
    Some(t)
}

/// fscanf(fp, "%lf", &x): Some(x), or None on EOF / matching failure.
pub fn fscanf_f64(fp: &CFile) -> Option<f64> {
    let t = scan_token(fp, true)?;
    if !t.bytes().any(|b| b.is_ascii_digit()) {
        return None;
    }
    Some(crate::par::strtod(&t).0)
}
/// fscanf(fp, "%f", &x).
pub fn fscanf_f32(fp: &CFile) -> Option<f32> {
    fscanf_f64(fp).map(|x| x as f32)
}
/// fscanf(fp, "%d", &x).
pub fn fscanf_i32(fp: &CFile) -> Option<i32> {
    let t = scan_token(fp, false)?;
    if !t.bytes().any(|b| b.is_ascii_digit()) {
        return None;
    }
    Some(crate::par::strtol(&t).0 as i32)
}

/// isatty(fd).
pub fn isatty(fd: i32) -> bool {
    extern "C" {
        #[link_name = "isatty"]
        fn c_isatty(fd: i32) -> i32;
    }
    // SAFETY: isatty only inspects the descriptor.
    unsafe { c_isatty(fd) == 1 }
}

/// filestat: classify a file descriptor exactly like par/lib/filestat.c.
pub fn filestat(fd: i32) -> FileType {
    // SAFETY: borrowed fd, released with into_raw_fd below.
    let f = unsafe { File::from_raw_fd(fd) };
    let md = f.metadata();
    let _ = f.into_raw_fd();
    let md = match md {
        Ok(m) => m,
        Err(_) => return FileType::BadFileType,
    };
    if isatty(fd) {
        return FileType::Tty;
    }
    let ft = md.file_type();
    if ft.is_file() {
        return FileType::Disk;
    }
    if ft.is_dir() {
        return FileType::Directory;
    }
    if md.nlink() == 0 {
        return FileType::Pipe;
    }
    if ft.is_fifo() {
        return FileType::Fifo;
    }
    if ft.is_socket() {
        return FileType::Socket;
    }
    if ft.is_symlink() {
        return FileType::Symlink;
    }
    if let Ok(dn) = std::fs::metadata("/dev/null") {
        if md.rdev() == dn.rdev() {
            return FileType::Disk;
        }
    }
    FileType::BadFileType
}

/// eremove: remove a file with error check (syserr on failure).
pub fn eremove(file: &str) -> i32 {
    if std::fs::remove_file(file).is_err() && std::fs::remove_dir(file).is_err() {
        crate::syserr!("%s: eremove: remove failed", "subcalls.c");
    }
    0
}

/// erename: rename with error check.
pub fn erename(oldfile: &str, newfile: &str) -> i32 {
    if std::fs::rename(oldfile, newfile).is_err() {
        crate::syserr!("%s: erename: rename failed", "subcalls.c");
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tmpfile_rw_seek() {
        let fp = etmpfile();
        let x = [1.0f32, 2.0, 3.5];
        assert_eq!(efwrite(&x, &fp), 3);
        assert_eq!(ftello(&fp), 12);
        erewind(&fp);
        let mut y = [0f32; 4];
        assert_eq!(efread(&mut y, &fp), 3);
        assert_eq!(&y[..3], &x);
        assert!(feof(&fp));
        efseeko(&fp, 4, SEEK_SET);
        let mut z = [0f32; 1];
        efread(&mut z, &fp);
        assert_eq!(z[0], 2.0);
        assert_eq!(ftello(&fp), 8);
        efseeko(&fp, 0, SEEK_END);
        fputs("12 3.5e1 x\n", &fp);
        efseeko(&fp, 12, SEEK_SET);
        assert_eq!(fscanf_i32(&fp), Some(12));
        assert_eq!(fscanf_f64(&fp), Some(35.0));
        assert_eq!(fscanf_f64(&fp), None);
        assert_eq!(fgetc(&fp), b'x' as i32);
        assert_eq!(fgets(100, &fp).as_deref(), Some("\n"));
        assert_eq!(fgets(100, &fp), None);
        efclose(&fp);
    }
}
