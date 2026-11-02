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
//! like C. Buffers flush as they fill (and at newlines on a tty). Remaining
//! data is flushed by `fflush`/`efflush`/`efclose`, by dropping the last
//! handle of a file, or by [`exit`] -- every program
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
//! | `fgets` with arbitrary input bytes | `cio::fgets_bytes(n, &fp) -> Option<Vec<u8>>` |
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
    bufsize: usize,
    seek_optimizable: bool,
    eof: bool,
    error: bool,
    child: Option<std::process::Child>,
    pub(crate) trin: Option<TrIn>,
    pub(crate) trout: Option<TrOut>,
}

// Darwin's fallback BUFSIZ. Regular files normally use st_blksize (4096),
// and pipes use their own block size (16384 on Darwin).
const BUFSZ: usize = 1024;

impl Inner {
    fn new(file: File, owned: bool, buffering: Buffering) -> Inner {
        let fd = file.as_raw_fd();
        let metadata = file.metadata().ok();
        let bufsize = metadata
            .as_ref()
            .map(|m| m.blksize() as usize)
            .filter(|&n| n > 0)
            .unwrap_or(BUFSZ);
        extern "C" {
            fn fcntl(fd: i32, cmd: i32, ...) -> i32;
        }
        // SAFETY: F_GETFL (3 on Darwin/Linux) only inspects the descriptor.
        // libc optimizes seeks within read-only regular-file buffers, but
        // deliberately does not do so for update streams (r+/w+/a+).
        let flags = unsafe { fcntl(fd, 3) };
        let seek_optimizable = flags >= 0
            && flags & 3 == 0
            && metadata.as_ref().is_some_and(|m| m.is_file());
        let buffering = if buffering == Buffering::Full && file.is_terminal() {
            Buffering::Line
        } else {
            buffering
        };
        Inner {
            file: Some(file),
            fd,
            owned,
            rbuf: Vec::new(),
            rpos: 0,
            wbuf: Vec::new(),
            buffering,
            bufsize,
            seek_optimizable,
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
            if want >= self.bufsize {
                self.rbuf.clear();
                self.rpos = 0;
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
                self.rbuf.resize(self.bufsize, 0);
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
        if data.is_empty() {
            return 0;
        }
        if self.rpos < self.rbuf.len() || !self.rbuf.is_empty() {
            self.drop_readahead();
        }
        if self.file.is_none() {
            self.error = true;
            return 0;
        }
        // __sfvwrite: write whole blocks directly when the buffer is empty.
        // Otherwise fill it, flushing only when more data needs room. An
        // exact fill of an existing full buffer stays pending until next time.
        let mut n = 0;
        while n < data.len() {
            let remaining = &data[n..];
            if self.buffering == Buffering::None
                || (self.buffering == Buffering::Full
                    && self.wbuf.is_empty()
                    && remaining.len() >= self.bufsize)
            {
                let count = if self.buffering == Buffering::None {
                    remaining.len().min(BUFSZ)
                } else {
                    remaining.len() / self.bufsize * self.bufsize
                };
                match self.file.as_mut().unwrap().write(&remaining[..count]) {
                    Ok(k) if k > 0 => n += k,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    _ => {
                        self.error = true;
                        break;
                    }
                }
            } else {
                let room = self.bufsize - self.wbuf.len();
                let mut count = room.min(remaining.len());
                let mut flush = remaining.len() > room;
                if self.buffering == Buffering::Line {
                    if let Some(newline) = remaining[..count].iter().position(|&c| c == b'\n') {
                        count = newline + 1;
                        flush = true;
                    }
                    flush |= count == room;
                }
                self.wbuf.extend_from_slice(&remaining[..count]);
                // On a failed buffer flush libc does not count the bytes
                // just copied, even if part of the buffer reached the fd.
                if flush && !self.flush() {
                    break;
                }
                n += count;
            }
        }
        n
    }

    fn seek(&mut self, off: i64, whence: i32) -> bool {
        if !self.flush() {
            return false;
        }
        let unread = (self.rbuf.len() - self.rpos) as i64;
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
        if self.seek_optimizable {
            if let Ok(current) = f.stream_position() {
                let target = match whence {
                    SEEK_SET => Some(off),
                    SEEK_CUR => (current as i64 - unread).checked_add(off),
                    SEEK_END => f.metadata().ok().and_then(|m| (m.len() as i64).checked_add(off)),
                    _ => None,
                };
                let start = current as i64 - self.rbuf.len() as i64;
                if let Some(target) = target {
                    if !self.rbuf.is_empty() && (start..=current as i64).contains(&target) {
                        self.rpos = (target - start) as usize;
                        self.eof = false;
                        return true;
                    }
                    if target >= 0 {
                        // Outside the current buffer, Darwin seeks to the
                        // containing block and reads ahead to the target.
                        // This also matters if another fopen later truncates
                        // the file: the rest of this block remains readable.
                        let skip = target as usize % self.bufsize;
                        if f.seek(SeekFrom::Start(target as u64 - skip as u64)).is_ok() {
                            self.rbuf.clear();
                            self.rpos = 0;
                            self.eof = false;
                            if skip == 0 {
                                return true;
                            }
                            self.rbuf.resize(self.bufsize, 0);
                            match f.read(&mut self.rbuf) {
                                Ok(k) => {
                                    self.rbuf.truncate(k);
                                    if k >= skip {
                                        self.rpos = skip;
                                        return true;
                                    }
                                }
                                Err(_) => self.error = true,
                            }
                            // Seeking past EOF still succeeds, even when
                            // the read-ahead cannot reach the target.
                            self.rbuf.clear();
                        }
                        if f.seek(SeekFrom::Start(target as u64)).is_err() {
                            return false;
                        }
                        self.rbuf.clear();
                        self.rpos = 0;
                        self.eof = false;
                        return true;
                    }
                }
            }
        }
        // A failed lseek (notably rewind on a pipe) leaves read-ahead intact.
        if f.seek(pos).is_err() {
            return false;
        }
        self.rbuf.clear();
        self.rpos = 0;
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

    /// Wrap an open `File` as a C stream (line buffered on a tty, otherwise
    /// fully buffered; closed on drop).
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
    // libc's standard FILE slots precede the slots allocated by fopen/fdopen,
    // even when our first access to stdout occurs after an auxiliary open.
    STD.with(|s| {
        for fp in [&s.0, &s.1, &s.2] {
            if let Ok(mut i) = fp.0.try_borrow_mut() {
                i.flush();
            }
        }
    });
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
    // Darwin tmpfile does not interpret an explicitly empty TMPDIR as cwd.
    if std::env::var_os("TMPDIR").is_some_and(|s| s.is_empty()) {
        crate::err!("%s: etmpfile: tmpfile failed", "subcalls.c");
    }
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

/// fgets on raw bytes: read at most n-1 bytes, stopping after a newline.
/// The terminating C NUL is implicit; NULs in the input are preserved.
/// None if no bytes were read before EOF or a read error (C returns NULL).
/// Darwin returns a partial line even when the next refill fails.
pub fn fgets_bytes(n: usize, fp: &CFile) -> Option<Vec<u8>> {
    if n == 0 {
        return None;
    }
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
    if v.is_empty() && n > 1 {
        None
    } else {
        Some(v)
    }
}

/// Text compatibility interface for existing ported callers. A Rust String
/// cannot represent arbitrary C input bytes; callers needing byte fidelity
/// must use fgets_bytes and byte-oriented parsing/formatting throughout.
pub fn fgets(n: usize, fp: &CFile) -> Option<String> {
    fgets_bytes(n, fp).map(|v| String::from_utf8_lossy(&v).into_owned())
}

fn scan_token(fp: &CFile, float: bool) -> Option<String> {
    // C isspace includes vertical tab, unlike u8::is_ascii_whitespace.
    let mut c = fgetc(fp);
    while c != EOF && matches!(c as u8, b' ' | b'\t'..=b'\r') {
        c = fgetc(fp);
    }
    if c == EOF {
        return None;
    }
    #[derive(Clone, Copy)]
    enum State {
        Start,
        Mantissa { hex: bool, digits: bool, dot: bool },
        Exponent { digits: bool, sign: bool },
        Infinity(usize),
        Nan(usize),
        Payload,
        Done,
    }
    use State::*;

    let mut t = Vec::new();
    let mut state = Start;
    let mut valid = 0;
    while c != EOF {
        let ch = c as u8;
        let lower = ch.to_ascii_lowercase();
        let next = match state {
            Start if t.is_empty() && matches!(ch, b'+' | b'-') => Some((Start, false)),
            Start if ch.is_ascii_digit() => Some((
                Mantissa {
                    hex: false,
                    digits: true,
                    dot: false,
                },
                true,
            )),
            Start if float && ch == b'.' => Some((
                Mantissa {
                    hex: false,
                    digits: false,
                    dot: true,
                },
                false,
            )),
            Start if float && lower == b'i' => Some((Infinity(1), false)),
            Start if float && lower == b'n' => Some((Nan(1), false)),
            Mantissa { hex, dot, .. } if ch.is_ascii_digit() || (hex && ch.is_ascii_hexdigit()) => {
                Some((
                    Mantissa {
                        hex,
                        digits: true,
                        dot,
                    },
                    true,
                ))
            }
            Mantissa {
                hex: false,
                dot: false,
                ..
            } if float && lower == b'x' && (t == b"0" || t == b"+0" || t == b"-0") => Some((
                Mantissa {
                    hex: true,
                    digits: false,
                    dot: false,
                },
                false,
            )),
            Mantissa {
                hex,
                digits,
                dot: false,
            } if float && ch == b'.' => Some((
                Mantissa {
                    hex,
                    digits,
                    dot: true,
                },
                digits,
            )),
            Mantissa {
                hex, digits: true, ..
            } if float && lower == if hex { b'p' } else { b'e' } => Some((
                Exponent {
                    digits: false,
                    sign: false,
                },
                false,
            )),
            Exponent { .. } if ch.is_ascii_digit() => Some((
                Exponent {
                    digits: true,
                    sign: true,
                },
                true,
            )),
            Exponent {
                digits: false,
                sign: false,
            } if matches!(ch, b'+' | b'-') => Some((
                Exponent {
                    digits: false,
                    sign: true,
                },
                false,
            )),
            Infinity(n) if b"infinity".get(n) == Some(&lower) => {
                Some((Infinity(n + 1), n + 1 == 3 || n + 1 == 8))
            }
            Nan(n) if b"nan".get(n) == Some(&lower) => Some((Nan(n + 1), n + 1 == 3)),
            Nan(3) if ch == b'(' => Some((Payload, false)),
            Payload if ch.is_ascii_alphanumeric() || ch == b'_' => Some((Payload, false)),
            Payload if ch == b')' => Some((Done, true)),
            _ => None,
        };
        let Some((next, complete)) = next else { break };
        state = next;
        t.push(ch);
        if complete {
            valid = t.len();
        }
        c = fgetc(fp);
    }
    if c != EOF {
        ungetc(c, fp);
    }
    // Darwin scanf backs up to the last complete subject sequence: 1e+
    // converts 1 and leaves e+, and a lone sign is a matching failure with
    // the sign unread. ungetc also clears EOF when lookahead reached it.
    for &ch in t[valid..].iter().rev() {
        ungetc(ch as i32, fp);
    }
    t.truncate(valid);
    if t.is_empty() {
        None
    } else {
        Some(String::from_utf8(t).unwrap())
    }
}

/// fscanf(fp, "%lf", &x): Some(x), or None on EOF / matching failure.
pub fn fscanf_f64(fp: &CFile) -> Option<f64> {
    extern "C" {
        #[link_name = "strtod"]
        fn c_strtod(s: *const std::ffi::c_char, end: *mut *mut std::ffi::c_char) -> f64;
    }
    let mut t = scan_token(fp, true)?.into_bytes();
    t.push(0);
    // SAFETY: the scanner supplies a complete, NUL-terminated ASCII token.
    Some(unsafe { c_strtod(t.as_ptr().cast(), std::ptr::null_mut()) })
}
/// fscanf(fp, "%f", &x).
pub fn fscanf_f32(fp: &CFile) -> Option<f32> {
    extern "C" {
        #[link_name = "strtof"]
        fn c_strtof(s: *const std::ffi::c_char, end: *mut *mut std::ffi::c_char) -> f32;
    }
    let mut t = scan_token(fp, true)?.into_bytes();
    t.push(0);
    // Use the host libc, as scanf does, including hex floats and NaN payloads.
    // strtof rounds directly to float; strtod followed by a cast double-rounds
    // decimal inputs just above or below a binary32 midpoint.
    // SAFETY: the scanner supplies a complete, NUL-terminated ASCII token.
    Some(unsafe { c_strtof(t.as_ptr().cast(), std::ptr::null_mut()) })
}
/// fscanf(fp, "%d", &x).
pub fn fscanf_i32(fp: &CFile) -> Option<i32> {
    let t = scan_token(fp, false)?;
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
    fn fgets_preserves_raw_bytes_and_limits() {
        let fp = etmpfile();
        let input = b"\xff\xc3\xa9\0z\nnext";
        fwrite(input, &fp);
        erewind(&fp);
        assert_eq!(fgets_bytes(0, &fp), None);
        assert_eq!(fgets_bytes(1, &fp), Some(Vec::new()));
        assert_eq!(ftello(&fp), 0);
        assert_eq!(fgets_bytes(3, &fp).as_deref(), Some(&b"\xff\xc3"[..]));
        assert_eq!(fgets_bytes(100, &fp).as_deref(), Some(&b"\xa9\0z\n"[..]));
        assert_eq!(fgets_bytes(100, &fp).as_deref(), Some(&b"next"[..]));
        assert!(feof(&fp));
        assert_eq!(fgets_bytes(100, &fp), None);
        assert_eq!(fgets_bytes(1, &fp), Some(Vec::new()));
    }

    fn on_disk(fp: &CFile) -> u64 {
        fp.with(|i| i.file.as_ref().unwrap().metadata().unwrap().len())
    }

    #[test]
    fn full_buffer_flush_points() {
        let fp = etmpfile();
        let block = fp.with(|i| i.bufsize);
        #[cfg(target_os = "macos")]
        assert_eq!(block, 4096);
        // A block written to an empty buffer goes straight to the file.
        assert_eq!(fwrite(&vec![b'a'; block], &fp), block);
        assert_eq!(on_disk(&fp), block as u64);
        assert_eq!(fputc(b'b' as i32, &fp), b'b' as i32);
        assert_eq!(fwrite(&vec![b'b'; block - 1], &fp), block - 1);
        // Exactly filling an already occupied buffer does not flush it.
        assert_eq!(on_disk(&fp), block as u64);
        assert_eq!(ftello(&fp), (2 * block) as i64);
        fputc(b'c' as i32, &fp);
        assert_eq!(on_disk(&fp), (2 * block) as u64);
        assert_eq!(fflush(&fp), 0);
        assert_eq!(on_disk(&fp), (2 * block + 1) as u64);
    }

    #[test]
    fn fdopen_flushes_blocks_and_retains_tail() {
        let original = etmpfile();
        let file = original.with(|i| i.file.as_ref().unwrap().try_clone().unwrap());
        let fp = CFile::from_file(file);
        let block = fp.with(|i| i.bufsize);
        let data = vec![42u8; 3 * block + 17];
        assert_eq!(fwrite(&data, &fp), data.len());
        assert_eq!(on_disk(&fp), (3 * block) as u64);
        assert_eq!(ftello(&fp), data.len() as i64);
        assert_eq!(fclose(&fp), 0);
        erewind(&original);
        let mut result = vec![0u8; data.len()];
        assert_eq!(fread(&mut result, &original), data.len());
        assert_eq!(result, data);
    }

    #[test]
    fn line_buffer_keeps_suffix_after_last_newline() {
        let fp = etmpfile();
        fp.with(|i| i.buffering = Buffering::Line);
        fputs("first\nsecond\ntail", &fp);
        assert_eq!(on_disk(&fp), 13);
        assert_eq!(ftello(&fp), 17);
        fflush(&fp);
        assert_eq!(on_disk(&fp), 17);
    }

    #[test]
    #[cfg(target_os = "macos")]
    fn fdopen_terminal_is_line_buffered() {
        extern "C" {
            fn openpty(
                master: *mut i32,
                slave: *mut i32,
                name: *mut u8,
                termios: *const std::ffi::c_void,
                winsize: *const std::ffi::c_void,
            ) -> i32;
            fn fcntl(fd: i32, cmd: i32, ...) -> i32;
        }
        let (mut master, mut slave) = (-1, -1);
        // SAFETY: openpty initializes both descriptors with default terminal
        // settings. F_SETFL=4 and O_NONBLOCK=4 on Darwin keep reads bounded.
        let (mut reader, writer) = unsafe {
            assert_eq!(
                openpty(&mut master, &mut slave, std::ptr::null_mut(),
                    std::ptr::null(), std::ptr::null()),
                0
            );
            assert_eq!(fcntl(master, 4, 4), 0);
            (File::from_raw_fd(master), File::from_raw_fd(slave))
        };
        let fp = CFile::from_file(writer);
        fputs("prompt", &fp);
        let mut bytes = [0; 32];
        assert_eq!(
            reader.read(&mut bytes).unwrap_err().kind(),
            std::io::ErrorKind::WouldBlock
        );
        fputs("\ntail", &fp);
        let count = reader.read(&mut bytes).unwrap();
        assert_eq!(&bytes[..count], b"prompt\r\n");
        assert_eq!(
            reader.read(&mut bytes).unwrap_err().kind(),
            std::io::ErrorKind::WouldBlock
        );
    }

    #[test]
    fn failed_pipe_seek_preserves_input() {
        extern "C" {
            fn pipe(fds: *mut i32) -> i32;
        }
        let mut fds = [-1; 2];
        // SAFETY: pipe initializes both descriptors; each gets one owner.
        let (reader, mut writer) = unsafe {
            assert_eq!(pipe(fds.as_mut_ptr()), 0);
            (File::from_raw_fd(fds[0]), File::from_raw_fd(fds[1]))
        };
        writer.write_all(b"abcdef").unwrap();
        drop(writer);
        let fp = CFile::from_file(reader);
        assert_eq!(fgetc(&fp), b'a' as i32);
        assert_eq!(fseeko(&fp, 0, SEEK_SET), -1);
        assert_eq!(fgetc(&fp), b'b' as i32);
        fp.with(|i| i.error = true);
        erewind(&fp);
        assert!(!ferror(&fp));
        assert!(!feof(&fp));
        let mut rest = [0u8; 5];
        assert_eq!(fread(&mut rest, &fp), 4);
        assert_eq!(&rest[..4], b"cdef");
        assert!(feof(&fp));
    }

    #[test]
    fn readonly_seek_retains_input_across_truncation() {
        // Give the read-only stream a separate open file description, as in
        // sultt. The temporary pathname is removed before any assertions.
        for size in [2156, 10000] {
            let path = std::env::temp_dir().join(format!("oxwave_seek_{}", std::process::id()));
            let mut writer = OpenOptions::new().write(true).create_new(true).open(&path).unwrap();
            writer.write_all(&vec![b'x'; size]).unwrap();
            let reader = File::open(&path).unwrap();
            std::fs::remove_file(path).unwrap();
            let fp = CFile::from_file(reader);
            assert_eq!(fread(&mut [0u8; 308], &fp), 308);
            assert_eq!(fseeko(&fp, 0, SEEK_END), 0);
            assert_eq!(ftello(&fp), size as i64);
            assert_eq!(fseeko(&fp, 308, SEEK_SET), 0);
            writer.set_len(0).unwrap();
            let mut tail = vec![0u8; fp.with(|i| i.bufsize).min(size) - 308];
            assert_eq!(fread(&mut tail, &fp), tail.len());
            assert!(tail.iter().all(|&c| c == b'x'));
            assert_eq!(fgetc(&fp), EOF);
            assert_eq!(fseeko(&fp, 9000, SEEK_SET), 0);
            assert_eq!(ftello(&fp), 9000);
            assert!(!feof(&fp));
        }
    }

    fn child_test(name: &str) -> std::process::Command {
        let mut command = std::process::Command::new(std::env::current_exe().unwrap());
        command.args(["--exact", name, "--nocapture"]);
        command.env("OXWAVE_CIO_TEST_CHILD", name);
        command
    }

    fn is_child(name: &str) -> bool {
        std::env::var("OXWAVE_CIO_TEST_CHILD").as_deref() == Ok(name)
    }

    #[test]
    fn exit_flushes_stdout_before_earlier_opened_file() {
        const NAME: &str = "cio::tests::exit_flushes_stdout_before_earlier_opened_file";
        if is_child(NAME) {
            // dup works even where the sandbox denies fopen("/dev/stdout").
            use std::os::fd::AsFd;
            let file = File::from(std::io::stdout().as_fd().try_clone_to_owned().unwrap());
            let outpar = CFile::from_file(file);
            fputs("[outpar]", &outpar);
            fputs("[stdout]", &stdout());
            exit(0);
        }
        let result = child_test(NAME).output().unwrap();
        assert!(result.status.success(), "{:?}", result);
        assert!(result.stdout.ends_with(b"[stdout][outpar]"), "{:?}", result);
    }

    #[test]
    fn empty_tmpdir_fails() {
        const NAME: &str = "cio::tests::empty_tmpdir_fails";
        if is_child(NAME) {
            etmpfile();
            panic!("empty TMPDIR was accepted");
        }
        let result = child_test(NAME).env("TMPDIR", "").output().unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&result.stderr).contains("subcalls.c: etmpfile: tmpfile failed"));
    }

    #[test]
    #[cfg(target_os = "macos")]
    fn partial_writes_count_complete_items() {
        const NAME: &str = "cio::tests::partial_writes_count_complete_items";
        if is_child(NAME) {
            #[repr(C)]
            struct Rlimit {
                current: u64,
                maximum: u64,
            }
            extern "C" {
                fn setrlimit(resource: i32, limit: *const Rlimit) -> i32;
                fn signal(sig: i32, handler: usize) -> usize;
            }
            // SAFETY: Darwin RLIMIT_FSIZE=1 and SIGXFSZ=25. Isolate the
            // process-wide limit and signal disposition in this child.
            unsafe {
                signal(25, 1); // SIG_IGN
                assert_eq!(setrlimit(1, &Rlimit { current: 1024, maximum: 1024 }), 0);
            }
            for count in [4097, 32769] {
                let fp = etmpfile();
                assert_eq!(fwrite(&vec![0f32; count], &fp), 256);
                assert!(ferror(&fp));
                assert_eq!(on_disk(&fp), 1024);
            }
            let fp = etmpfile();
            // Header-sized items first fill the stdio buffer. A failed
            // flush of that buffer must not count the current partial item.
            for _ in 0..17 {
                assert_eq!(fwrite(&[[0u8; 240]], &fp), 1);
            }
            assert_eq!(fwrite(&[[0u8; 240]], &fp), 0);
            assert!(ferror(&fp));
            let fp = etmpfile();
            fp.with(|i| i.buffering = Buffering::None);
            assert_eq!(fwrite(&vec![0f32; 4097], &fp), 256);
            assert!(ferror(&fp));
            return;
        }
        let result = child_test(NAME).output().unwrap();
        assert!(result.status.success(), "{:?}", result);
    }

    fn scan_input(s: &str) -> CFile {
        let fp = etmpfile();
        fputs(s, &fp);
        erewind(&fp);
        fp
    }

    #[test]
    fn scanf_hex_floats() {
        for (s, want) in [
            ("0x1.8p2", 6.0),
            ("-0X.8P+2", -2.0),
            ("+0x1c2p-1", 225.0),
            ("0x1.e", 1.875),
            ("0x1p-149", f32::from_bits(1)),
            ("0x1.000001p0", 1.0),
            ("0x1.000001000001p0", f32::from_bits(0x3f800001)),
            ("0x1p128", f32::INFINITY),
            ("-0x0p0", -0.0),
        ] {
            let fp = scan_input(&format!("{s},7"));
            assert_eq!(fscanf_f32(&fp).unwrap().to_bits(), want.to_bits(), "{s}");
            assert_eq!(fgetc(&fp), b',' as i32);
            assert_eq!(fscanf_i32(&fp), Some(7));
        }
        for (s, want) in [
            ("0x1.fffffffffffffp1023", f64::MAX),
            ("0x1p-1074", f64::from_bits(1)),
            ("-0X1.8P2", -6.0),
        ] {
            assert_eq!(
                fscanf_f64(&scan_input(s)).unwrap().to_bits(),
                want.to_bits()
            );
        }
    }

    #[test]
    fn scanf_special_floats() {
        for s in ["inf", "+INF", "infinity", "InFiNiTy"] {
            assert_eq!(fscanf_f32(&scan_input(s)), Some(f32::INFINITY));
            assert_eq!(fscanf_f64(&scan_input(s)), Some(f64::INFINITY));
        }
        assert_eq!(
            fscanf_f32(&scan_input("-infinity")),
            Some(f32::NEG_INFINITY)
        );
        assert_eq!(fscanf_f64(&scan_input("-inf")), Some(f64::NEG_INFINITY));
        for s in ["nan", "NaN", "+nan()", "-NAN(payload_123)"] {
            let fp = scan_input(&format!("{s} 42"));
            assert!(fscanf_f32(&fp).unwrap().is_nan(), "{s}");
            assert_eq!(fscanf_i32(&fp), Some(42));
            assert!(fscanf_f64(&scan_input(s)).unwrap().is_nan(), "{s}");
        }
    }

    #[test]
    fn scanf_restores_incomplete_suffixes() {
        for (s, want, rest) in [
            ("1e+", 1.0, "e+"),
            ("25E- 3", 25.0, "E- 3"),
            ("1e 2", 1.0, "e 2"),
            ("0x1p+", 1.0, "p+"),
            ("0x1P-z", 1.0, "P-z"),
            ("0x.", 0.0, "x."),
            ("0xp1", 0.0, "xp1"),
            ("infinite", f64::INFINITY, "inite"),
        ] {
            for single in [false, true] {
                let fp = scan_input(s);
                let got = if single {
                    fscanf_f32(&fp).map(f64::from)
                } else {
                    fscanf_f64(&fp)
                };
                assert_eq!(got, Some(want), "{s}");
                assert!(!feof(&fp), "{s}");
                assert_eq!(fgets(100, &fp).as_deref(), Some(rest), "{s}");
            }
        }
        for s in ["nan(foo", "nan(a!)"] {
            let fp = scan_input(s);
            assert!(fscanf_f32(&fp).unwrap().is_nan());
            assert!(!feof(&fp));
            assert_eq!(fgets(100, &fp).unwrap(), &s[3..]);
        }
    }

    #[test]
    fn scanf_matching_failure_preserves_input_and_eof() {
        for s in ["+", "-", ".", "+.", "in", "na", "+x"] {
            let fp = scan_input(s);
            assert_eq!(fscanf_f32(&fp), None, "{s}");
            assert!(!feof(&fp), "{s}");
            assert_eq!(fscanf_f64(&fp), None, "{s}");
            assert!(!feof(&fp), "{s}");
            assert_eq!(fgets(100, &fp).as_deref(), Some(s), "{s}");
        }
        for s in ["+", "-", "+x"] {
            let fp = scan_input(s);
            assert_eq!(fscanf_i32(&fp), None);
            assert!(!feof(&fp));
            assert_eq!(fgets(100, &fp).as_deref(), Some(s));
        }
        for s in ["", " \t\n\r\x0b\x0c"] {
            let fp = scan_input(s);
            assert_eq!(fscanf_f32(&fp), None);
            assert!(feof(&fp));
        }
        let fp = scan_input("42");
        assert_eq!(fscanf_i32(&fp), Some(42));
        assert!(feof(&fp));
    }

    #[test]
    fn scanf_c_whitespace() {
        let fp = scan_input("\x0b-12\x0b0x1p2\x0c.5\r\n\t 6");
        assert_eq!(fscanf_i32(&fp), Some(-12));
        assert_eq!(fscanf_f32(&fp), Some(4.0));
        assert_eq!(fscanf_f64(&fp), Some(0.5));
        assert_eq!(fscanf_i32(&fp), Some(6));
    }

    #[test]
    fn scanf_adjacent_float_conversions() {
        // Whitespace is optional between conversions. A sign, special value,
        // or decimal following a closed NaN payload starts the next token.
        let fp = scan_input("0x1.8p2-3.5infinitynan(0x123)17e+4,9");
        assert_eq!(fscanf_f32(&fp), Some(6.0));
        assert_eq!(fscanf_f64(&fp), Some(-3.5));
        assert_eq!(fscanf_f32(&fp), Some(f32::INFINITY));
        assert!(fscanf_f64(&fp).unwrap().is_nan());
        assert_eq!(fscanf_f32(&fp), Some(170000.0));
        assert_eq!(fgetc(&fp), b',' as i32);
        assert_eq!(fscanf_i32(&fp), Some(9));
        assert!(feof(&fp));
    }

    #[test]
    fn scanf_decimal_integer_subject_sequences() {
        // %d remains decimal even though the same tokenizer accepts hex,
        // exponents and special values for floating-point conversions.
        for (s, want, rest) in [
            ("0x1p2", 0, "x1p2"),
            ("-0Xf", 0, "Xf"),
            ("010,", 10, ","),
            ("25e+2", 25, "e+2"),
            ("2147483647.", i32::MAX, "."),
            ("-2147483648!", i32::MIN, "!"),
        ] {
            let fp = scan_input(s);
            assert_eq!(fscanf_i32(&fp), Some(want), "{s}");
            assert_eq!(ftello(&fp), (s.len() - rest.len()) as i64, "{s}");
            assert!(!feof(&fp), "{s}");
            assert_eq!(fgets(100, &fp).as_deref(), Some(rest), "{s}");
        }
        for s in ["inf", "nan(0x123)", "-infinity", ".5"] {
            let fp = scan_input(s);
            assert_eq!(fscanf_i32(&fp), None, "{s}");
            assert_eq!(ftello(&fp), 0, "{s}");
            assert!(!feof(&fp), "{s}");
            assert_eq!(fgets(100, &fp).as_deref(), Some(s), "{s}");
        }
    }

    #[test]
    fn scanf_decimal_rounds_directly_to_f32() {
        for (s, bits) in [
            ("8.00000047683715820312500000001", 0x41000001),
            ("16.0000009536743164062500000000000001", 0x41800001),
            ("1.0000000596046447753906250000000001", 0x3f800001),
            ("1.0000001788139343261718749999999999", 0x3f800001),
            ("-8.00000047683715820312500000001", 0xc1000001),
        ] {
            assert_eq!(fscanf_f32(&scan_input(s)).unwrap().to_bits(), bits, "{s}");
        }
    }

    #[test]
    #[cfg(target_os = "macos")]
    fn scanf_matches_darwin_subject_sequences() {
        extern "C" {
            #[link_name = "sscanf"]
            fn c_sscanf(s: *const std::ffi::c_char, fmt: *const std::ffi::c_char, ...) -> i32;
        }
        // Compare both conversion and consumed bytes with the libc used by
        // the original binaries, including malformed prefixes and suffixes.
        for prefix in ["", "+", "-", " \x0b"] {
            for body in [
                "",
                ".",
                "0",
                "00",
                "1.",
                ".5",
                "0x",
                "0x.",
                "0x.8",
                "0x1.e",
                "0X1C2",
                "0x1p-149",
                "1e-45",
                "1e39",
                "1e309",
                "in",
                "inf",
                "infinity",
                "nan",
                "nan()",
                "nan(0x123)",
                "nan(payload_1)",
            ] {
                for suffix in ["", "e", "E+", "p-", "3", "x", ".5", "(foo", "(a!)", ",7"] {
                    let s = format!("{prefix}{body}{suffix}");
                    let cs = std::ffi::CString::new(s.as_str()).unwrap();
                    let mut expected32 = 0.0f32;
                    let mut expected64 = 0.0f64;
                    let mut n32 = -1i32;
                    let mut n64 = -1i32;
                    // SAFETY: NUL-terminated strings and correctly typed
                    // writable pointers for each scanf conversion.
                    let (rc32, rc64) = unsafe {
                        (
                            c_sscanf(cs.as_ptr(), c"%f%n".as_ptr(), &mut expected32, &mut n32),
                            c_sscanf(cs.as_ptr(), c"%lf%n".as_ptr(), &mut expected64, &mut n64),
                        )
                    };
                    let fp = scan_input(&s);
                    let got = fscanf_f32(&fp);
                    assert_eq!(got.is_some(), rc32 == 1, "f32: {s:?}");
                    if let Some(got) = got {
                        assert_eq!(got.to_bits(), expected32.to_bits(), "f32: {s:?}");
                        assert_eq!(ftello(&fp), n32 as i64, "f32: {s:?}");
                    }
                    let fp = scan_input(&s);
                    let got = fscanf_f64(&fp);
                    assert_eq!(got.is_some(), rc64 == 1, "f64: {s:?}");
                    if let Some(got) = got {
                        assert_eq!(got.to_bits(), expected64.to_bits(), "f64: {s:?}");
                        assert_eq!(ftello(&fp), n64 as i64, "f64: {s:?}");
                    }
                }
            }
        }
    }

    #[test]
    #[cfg(target_os = "macos")]
    fn scanf_matches_darwin_stream_state() {
        extern "C" {
            fn fmemopen(
                buf: *mut std::ffi::c_void,
                size: usize,
                mode: *const std::ffi::c_char,
            ) -> *mut std::ffi::c_void;
            #[link_name = "fscanf"]
            fn c_fscanf(fp: *mut std::ffi::c_void, fmt: *const std::ffi::c_char, ...) -> i32;
            #[link_name = "ftell"]
            fn c_ftell(fp: *mut std::ffi::c_void) -> std::ffi::c_long;
            #[link_name = "feof"]
            fn c_feof(fp: *mut std::ffi::c_void) -> i32;
            #[link_name = "fgetc"]
            fn c_fgetc(fp: *mut std::ffi::c_void) -> i32;
            #[link_name = "fclose"]
            fn c_fclose(fp: *mut std::ffi::c_void) -> i32;
        }
        // sscanf cannot expose EOF or the unread bytes on a matching failure.
        // Check real libc streams too, especially a lone sign at end of input:
        // callers use feof to distinguish that failure from exhausted input.
        for s in [
            " ",
            "\t\n\x0b\x0c\r",
            "+",
            "-",
            ".",
            "+.",
            "-x",
            "in",
            "na",
            "x",
            "42",
            "42 ",
            "-12\x0b7",
            "1e+",
            "25E- 3",
            "1e 2",
            "0x",
            "0x.",
            "0x1p+",
            "0x1P-z",
            "0x1.8p2",
            "infinite",
            "infinity",
            "nan",
            "nan()",
            "nan(foo",
            "nan(a!)",
            "nan(0x123)",
            "-nan(payload_1)",
            "8.00000047683715820312500000001",
            "\0",
            "1\0e2",
        ] {
            for format in [c"%f", c"%lf", c"%d"] {
                let mut input = s.as_bytes().to_vec();
                let mut x32 = 0.0f32;
                let mut x64 = 0.0f64;
                let mut xi = 0i32;
                // SAFETY: input remains live and writable until fclose. Each
                // conversion receives a pointer of the type its format needs.
                let (rc, expected, pos, eof, next) = unsafe {
                    let fp = fmemopen(input.as_mut_ptr().cast(), input.len(), c"r".as_ptr());
                    assert!(!fp.is_null());
                    let (rc, bits) = match format.to_bytes() {
                        b"%f" => {
                            let rc = c_fscanf(fp, format.as_ptr(), &mut x32);
                            (rc, x32.to_bits() as u64)
                        }
                        b"%lf" => {
                            let rc = c_fscanf(fp, format.as_ptr(), &mut x64);
                            (rc, x64.to_bits())
                        }
                        _ => {
                            let rc = c_fscanf(fp, format.as_ptr(), &mut xi);
                            (rc, xi as u64)
                        }
                    };
                    let state = (rc, bits, c_ftell(fp), c_feof(fp) != 0, c_fgetc(fp));
                    c_fclose(fp);
                    state
                };
                let fp = scan_input(s);
                let got = match format.to_bytes() {
                    b"%f" => fscanf_f32(&fp).map(|x| x.to_bits() as u64),
                    b"%lf" => fscanf_f64(&fp).map(f64::to_bits),
                    _ => fscanf_i32(&fp).map(|x| x as u64),
                };
                assert_eq!(got, (rc == 1).then_some(expected), "{format:?}: {s:?}");
                assert_eq!(ftello(&fp), pos as i64, "{format:?}: {s:?}");
                assert_eq!(feof(&fp), eof, "{format:?}: {s:?}");
                assert_eq!(fgetc(&fp), next, "{format:?}: {s:?}");
            }
        }
    }

    #[test]
    fn scanf_suffix_rollback_across_refill() {
        for padding in [BUFSZ - 2, BUFSZ - 1, BUFSZ, BUFSZ + 1] {
            for single in [false, true] {
                let fp = scan_input(&format!("{}1e+", " ".repeat(padding)));
                let got = if single {
                    fscanf_f32(&fp).map(f64::from)
                } else {
                    fscanf_f64(&fp)
                };
                assert_eq!(got, Some(1.0));
                assert!(!feof(&fp));
                assert_eq!(ftello(&fp), (padding + 1) as i64);
                assert_eq!(fscanf_f32(&fp), None);
                assert!(!feof(&fp));
                assert_eq!(fgetc(&fp), b'e' as i32);
                assert_eq!(fgetc(&fp), b'+' as i32);
                assert_eq!(fgetc(&fp), EOF);
                assert!(feof(&fp));
            }
        }
    }

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
