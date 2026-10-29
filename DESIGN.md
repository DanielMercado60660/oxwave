# OxWave design

OxWave is a Rust port of CWP/SU (Seismic Un*x). Every ported program must
behave like the original SU binary in `vendor/su-build/bin`: same stdout
bytes (floats within 1e-5 relative), same exit code, same error messages,
same selfdoc. Every program builds on the core crate (`oxwave-core`), so it
changes rarely; functionality it lacks goes into a separate library crate.

## Workspace layout

```
Cargo.toml                 workspace: members crates/* and programs/*; std only
crates/oxwave-core/        core (header, trace I/O, getpar, err, selfdoc, stdio, printf)
crates/oxwave-cwp/         port of cwp/lib: `complex`, `cwp` (cwp.h macros); one module per function group
crates/oxwave-cref/        TEST ONLY: links the original libsu.a/libpar.a/libcwp.a for differential tests
programs/<prog>/           one crate per SU program; binary name == SU program name
tools/make_goldens.py      record outputs of the original SU binaries
tools/check_program.py     compare an OxWave binary against those goldens
tools/sdoc2rs.py           convert a C `sdoc[]` array into a Rust `const SDOC`
tools/cwp_header.txt       license header every ported file starts with
goldens/<prog>/            cases.txt (per-program cases) + generated goldens
vendor/su-original/src     original C source (read only)
vendor/su-build            original build (CWPROOT): bin/, lib/, include/
```

Toolchain: Rust 1.88 or later. Build with
`cargo build --release --offline -p <prog>`. No external crates: the workspace
uses std only and builds offline.

## Per-program crate template

`programs/sufoo/Cargo.toml`:
```toml
[package]
name = "sufoo"
version.workspace = true
edition.workspace = true
license.workspace = true
description.workspace = true
publish.workspace = true

[[bin]]
name = "sufoo"
path = "src/main.rs"

[dependencies]
oxwave-core.workspace = true
oxwave-cwp.workspace = true      # only if needed

[dev-dependencies]
oxwave-cref.workspace = true     # only for differential unit tests
```

`programs/sufoo/src/main.rs` (see `programs/sugethw/src/main.rs`, the reference port):
```rust
<contents of tools/cwp_header.txt, ending with "// Rust port for OxWave, a Rust port of CWP/SU.">

//! SUFOO - one-line purpose. Port of su/main/<family>/sufoo.c (Revision x.y).
//! Credits (from the C source): ...

use oxwave_core::prelude::*;

/// Self-documentation, byte-identical to the C `sdoc[]` array.
const SDOC: &[&str] = &[ ... ];          // generate: tools/sdoc2rs.py <file.c>

fn main() {
    let argv = args();
    let mut tr = Trace::new();
    initargs(&argv);
    requestdoc(1, SDOC);                 // same flag as the C
    ...                                  // line-by-line translation
    checkpars();
    while gettr(&mut tr) != 0 { ...; puttr(&mut tr); }
    exit(CWP_Exit());                    // ALWAYS end with exit(): it flushes stdout
}
```
Keep the C structure, variable names, comments that matter, and statement
order. Helper functions of the C file become private Rust functions in
the same file (or `src/<name>.rs` modules if large). Globals become locals
passed explicitly or a `struct` owned by `main`.

## C to Rust mapping (oxwave-core; `use oxwave_core::prelude::*`)

| C | Rust |
|---|---|
| `int main(int argc, char **argv)` | `let argv = args();` (`Vec<String>`, argv[0] as invoked) |
| `initargs(argc, argv)` | `initargs(&argv)` |
| `requestdoc(flag)` / `pagedoc()` | `requestdoc(flag, SDOC)` / `pagedoc(SDOC)` |
| `xargc`, `xargv` | `xargc()`, `xargv()` |
| `getparint("n", &n)` (also uint/short/ushort/long/ulong/float/double) | `getparint("n", &mut n)` returns `i32` count; `if getparint(..) == 0 { n = default; }` |
| `getparfloat("x", x)` into an array | `getparfloat("x", &mut x)` with `x: Vec<f32>` (overwritten from 0, grows) or `&mut x[..]` |
| `getparstring("s", &s)` | `let mut s = String::new(); getparstring("s", &mut s)` |
| `getparstringarray("k", k)` | `getparstringarray("k", &mut k)` (`Vec<String>`) |
| `getnparint(i, "n", &n)` etc. | same name, same args |
| `countparval`, `countnparval`, `countparname` | same, return `i32` |
| `checkpars()` | `checkpars()` |
| `MUSTGETPARINT("n",&n)` (FLOAT, DOUBLE, STRING) | `mustgetparint("n", &mut n)` |
| `eatoi`, `eatof`, ... `atoi`, `atof`, `strtod(s,NULL)`, `strtol(s,NULL,10)` | same names; `strtod(s).0`, `strtol(s).0` |
| `err(fmt, ...)` / `warn(...)` / `syserr(...)` | `err!("fmt %d", x)` / `warn!(..)` / `syserr!(..)` with the C format string verbatim |
| `printf(fmt, ...)`, `fprintf(fp, ...)`, `sprintf(buf, ...)` | `printf!(..)`, `fprintf!(fp, ..)` (`fp: CFile`, e.g. `stderr()`), `let buf = sprintf!(..)` |
| `putchar(c)`, `fputs(s, fp)`, `fputc`, `fgetc/getc`, `ungetc`, `fgets(buf,n,fp)` | `putchar(c as i32)`, `fputs(s, &fp)`, ..., `fgets(n, &fp) -> Option<String>` |
| `fscanf(fp,"%f",&x)` / `"%lf"` / `"%d"` | `fscanf_f32(&fp)` / `fscanf_f64` / `fscanf_i32` -> `Option` (None = EOF/no match) |
| `segy tr;` | `let mut tr = Trace::new();` (`tr.data` has SU_NFLTS=65535 samples, zeroed) |
| `tr.ns`, `tr.offset`, `tr.data[i]` | `tr.ns`, `tr.offset` (Trace derefs to `SegyHeader`), `tr.data[i]` |
| `(char *) tr.data`, `(short *) tr.data` | `tr.data_bytes()`, `tr.data_bytes_mut()` |
| `gettr(&tr)` / `puttr(&tr)` | `gettr(&mut tr) != 0` / `puttr(&mut tr)` |
| `fgettr(fp,&tr)`, `fputtr(fp,&tr)`, `fvgettr`, `fvputtr`, `vgettr`, `vputtr` | same names, `(&fp, &mut tr)` |
| `fgettra(fp, &tr, itr)` / `gettra(&tr, itr)` | same names |
| `stdin`, `stdout`, `stderr` (FILE *) | `stdin()`, `stdout()`, `stderr()` (`CFile`; clones are the same stream) |
| `efopen(name, "r")`, `efclose(fp)` | `let fp = efopen(name, "r"); efclose(&fp);` |
| `etmpfile()`, `erewind(fp)`, `efseeko(fp,o,SEEK_SET)`, `eftello(fp)` | same names, `&fp` |
| `efread(buf, FSIZE, n, fp)` / `efwrite(buf, FSIZE, n, fp)` | `efread(&mut buf[..n], &fp)` / `efwrite(&buf[..n], &fp)`, any numeric element type, native byte order, returns items |
| `epopen(cmd,"r")`/`epclose(fp)`, `eremove`, `erename`, `filestat(fd)`, `isatty(fd)`, `fileno(fp)` | same names |
| `exit(code)` / `return EXIT_SUCCESS` from main | `exit(code)` (flushes every stream). `CWP_Exit()` returns 0 |
| `Value val; val.i = 3; val.f` | `let mut val = Value::default(); val.set_i(3); val.f()` (union semantics) |
| `gethval`, `puthval`, `gethdval`, `puthdval`, `hdtype`, `getindex`, `getkey`, `swaphval`, `printheader` | same names (`getkey` -> `Option<&str>`, index args `i32`) |
| `vtoi`, `vtol`, `vtof`, `vtod`, `valcmp`, `atoval`, `printfval`, `fprintfval`, `getparval`, `valtoabs` | same names, type arg is `&str` like `"i"` |
| `hdr[i].key/.type/.offs`, `SU_NKEYS`, `HDRBYTES`, `SU_NFLTS` | `HDR[i].key/.type_/.offs`, same constants |
| `bhed bh;`, `getbhval/putbhval/swapbhval`, `xdrbhdrsub` | `oxwave_core::segy::{Bhed, ..}`, `bh.to_xdr()`, `Bhed::from_xdr` |
| `tapesegy`/`tapebhed`, `gettapehval/puttapehval/gettapebhval/puttapebhval` | `oxwave_core::segy::{TAPEHDR, TAPEBHDR, gettapehval, ..}` on raw `[u8]` buffers |
| `ibm_to_float`, `float_to_ibm` (segyread/segywrite) | `oxwave_core::segy::{ibm_to_float, float_to_ibm}` (bit-exact, differential-tested) |
| `swap_short_2`, `swap_int_4`, `swap_float_4`, ... (cwp/lib/swapbyte.c) | `oxwave_core::segy::swap_*` (do not re-port) |
| `zebc` (segywrite), `dd conv=ascii` (segyread) | `oxwave_core::segy::{zebc, dd_conv_ascii, dd_conv_ebcdic}` |
| `ealloc1float(n)`, `alloc1float`, `ealloc1int`, `ealloc1complex` | `vec![0.0f32; n]`, `vec![0i32; n]`, `vec![Complex::default(); n]` |
| `ealloc2float(n1, n2)` -> `x[i2][i1]` | `vec![vec![0.0f32; n1]; n2]` (`x[i2][i1]`), or flat `vec![0.0; n1*n2]` with `x[i2*n1+i1]` when the C uses `x[0]` as one contiguous block (e.g. passes `x[0]` to a function, `memcpy`s it, or `fwrite`s `n1*n2`) |
| `free1float(x)` etc. | nothing (drop) |
| `memcpy`, `memset` | `copy_from_slice`, `fill(0.0)` |
| `STREQ(a,b)` | `a == b` |
| `complex`, `cmplx`, `cadd`, `cmul`, `cdiv`, `rcabs`, `crmul`, `conjg`, ... | `oxwave_cwp::complex::*` (same names; `Complex { r, i }`) |
| `PI`, `MAX`, `MIN`, `ABS`, `SGN`, `NINT`, `ISODD`, `CLOSETO` | `oxwave_cwp::cwp::{PI, max, min, abs, sgn, nint, isodd, closeto}` |

Library functions not listed (pfafft, intlin, ...) are ported in the
library crates (`oxwave-cwp` or a sibling crate) under their C
names; call them from there. Never re-implement a library routine inside a
program crate.

**Already in core, do not port again:** par/lib `errpkge.c`, `getpars.c`,
`atopkge.c`, `docpkge.c`, `subcalls.c`, `filestat.c`, `ealloc.c` (= `Vec`);
su/lib `fgettr.c`, `fputtr.c`, `hdrpkge.c` (incl. tape and binary header
accessors), `valpkge.c`, `xdrhdrsub.c`, `xdrbhdrsub.c`; cwp/lib
`complex.c`, `complexd.c` (clean-room), `swapbyte.c`, `zebc.c`, `exit.c`, and
the cwp.h macros. Deliberately not ported: `scanfval`, `isascii_txt` /
`isebcdic_txt` and every non-SUXDR branch (dead in this build). su/lib
`fgetgthr.c` / `fputgthr.c` (`get_gather` & co.) and the rest of su/lib,
par/lib and cwp/lib are ported in the library crates, on top of this core.

## Porting rules

1. **Line-by-line faithful translation.** Same control flow, same order of
   operations, same intermediate types. Do not "improve", vectorise,
   reassociate, use iterator sums, or replace a loop by a library call.
2. **float vs double exactly as C.** A C `float` is `f32`, `double` is `f64`.
   C promotes to double when *any* operand is double, and **unsuffixed
   literals are double**: `float y = x*0.5;` is `let y = ((x as f64)*0.5) as f32;`
   whereas `x*0.5f` stays `f32`. `PI` is a double. C90 math functions
   (`sqrt`, `sin`, `exp`, `pow`, `log`, `atan2`, `fabs`, `floor`...) take and
   return double: `sqrt(x)` with float x is `(x as f64).sqrt()`, then `as f32`
   if assigned to a float. Use `powf` for `pow` (never `powi`), `ln` for `log`,
   `%` on floats for `fmod`, `round_ties_even` for `rint`.
3. **Integer semantics.** `int` is `i32`, `long` is `i64`, `short` `i16`,
   `unsigned short` `u16`, `size_t` `usize`. Where C can overflow, use
   `wrapping_*`. `(int)x` from float is `x as i32` (truncates; saturates and
   maps NaN to 0, exactly what the arm64 instruction does). Integer division
   by zero returns 0 on arm64 in C but panics in Rust: if the C code can
   divide by zero, write `checked_div(..).unwrap_or(0)` and comment it.
   Integer `/` and `%` truncate toward zero in both languages.
4. **No fast-math, no FMA.** Never use `mul_add`. Rust does not contract.
5. **Out-of-bounds.** C may read past an array (garbage) where Rust panics.
   If a test case hits that, reproduce the C result deliberately (e.g. size
   the array like the C allocation) and comment it; never silently clamp.
6. **Messages and formats.** Copy every `err`/`warn`/`printf`/`fprintf`
   format string verbatim into `err!`/`warn!`/`printf!`/`fprintf!`. Library
   messages built from `__FILE__` keep the literal C file name.
7. **Selfdoc** is generated with `tools/sdoc2rs.py`; never retype it.
8. **Header comment.** Every ported file starts with `tools/cwp_header.txt`
   verbatim (it ends with `// Rust port for OxWave, a Rust port of CWP/SU.`),
   then a doc comment naming the C source and copying its credits. Never use
   CSM or contributor names to promote OxWave.
9. **Randomness/time.** If the C uses `rand`, `franuni`, time or pid, port
   the generator faithfully (seeded runs must match); unseeded/time cases are
   flagged UNSTABLE by the harness and skipped.
10. **Files the program writes** besides stdout (e.g. `outpar=`) are not
    compared automatically; add a case that `cat`s them to stdout.

## Testing a program

```sh
cargo build --release --offline -p sufoo
python3 tools/make_goldens.py sufoo      # original SU outputs -> goldens/sufoo/
python3 tools/check_program.py sufoo     # OxWave vs goldens; exit 0 = all pass
```
Cases = `@selfdoc` (no args on a pty, PAGER=cat) + the inventory's
`test_cases` + `goldens/sufoo/cases.txt` (one shell pipeline per
line, `#` comments). Each runs with `bash -o pipefail -c` in an empty temp
dir, stdin `/dev/null`, no controlling tty, with `vendor/su-build/bin`
first on PATH; for the check, a directory holding only a `sufoo` symlink to
the OxWave binary goes first, so every other stage is still the original. Always
invoke programs by bare name (argv[0] appears in error messages).

Add cases for: every parameter and its default, each output mode, error
exits (missing/invalid parameters, empty input), `verbose=1`, odd sizes
(nt=1, ntr=1, non-power-of-2), negative/zero values. Use SU generators for
input (`suplane`, `suspike`, `sunull`, `suaddnoise seed=..`, `sushw`/`suchw`
to set headers). Tags before a command: `[ws]` whitespace-normalised text,
`[exact]` exact bytes, `[rtol=1e-4]`/`[atol=1e-6]` tolerances (justify in a
comment), `[stderr]` compare whole stderr, `[f32]` raw float output.

Comparison: SU streams trace by trace, header bytes identical, samples
within `rtol*max|golden trace| + atol` (defaults 1e-5, 1e-30); text exact;
selfdoc exact; exit code always; first stderr line when the original fails.

Library crates have differential tests in their `tests/` using
`oxwave-cref` (declare the C prototype in an `extern "C"` block, compare
on thousands of `XorShift`-generated inputs; see
`crates/oxwave-cref/tests/diff_core.rs`). Functions `static` inside SU
programs can be copied verbatim into `crates/oxwave-cref/csrc/shim.c` with an
`oxcref_` prefix (cref is shared by every crate's tests).

## Gotchas (all verified against this build)

- **Byte order.** SU here is built with `-DSUXDR`: every trace stream is
  **big-endian** (each header field big-endian at its segy.h offset,
  samples big-endian IEEE). `trace` handles it; never write traces by hand.
  Native-endian input fails with `fgettr.c: no data on first trace`. Raw
  float outputs (`sustrip`, `sugethw output=binary`, `efwrite` of floats) are
  native little-endian. No SEG-Y line headers on SU streams.
- **First trace fixes ns.** gettr reads nsfirst samples per trace and errs
  `fgettr.c: on trace #<n>, number of samples in header (<ns>) differs from
  number for first trace (<nsfirst>)`; puttr always writes nsfirst samples.
  Use `vgettr`/`vputtr` only where the C does. Empty input:
  `fgettr.c: bad first header`. stdin a tty: `fgettr.c: segy input can't be tty`.
- **`tr.data` keeps stale samples** past ns from earlier traces, like C.
- **err/warn format**: stderr gets `"\n" + argv0 + ": " + msg + "\n"` (the
  first stderr line is empty); stdout is flushed first; err exits 1.
- **Selfdoc** goes through `sh -c "$PAGER 1>&2"` (default `more`), so it lands
  on **stderr**, exit status 1. `requestdoc(1)` only fires if argc==1 and
  stdin is a tty; `requestdoc(0)` needs stdin and stdout ttys; `requestdoc(n>1)`
  fires if argc<=n. Without a tty the program reads stdin instead.
- **checkpars()** prints `Unknown <argv0> argument <arg>` for unrequested
  `name=` args (does not exit), prints `Using Big Endian SU data format w/ XDR.`
  when `verbose=1`, and later getpar calls print
  `getpar() call after checkpars(): <name>`. Call it exactly where the C does.
- **getpar quirks**: last occurrence wins; `par=file` entries come before the
  command line; `n=` (empty) returns 0 for numbers but 1 (`""`) for strings;
  `1,,2` gives [1,0,2]; `abc` gives 0 silently; trailing garbage ignored;
  ints are base 10; out-of-range errs `atopkge.c: eatoi: overflow` (etc.);
  floats like `1e-310` err `atopkge.c: eatod: overflow` (macOS strtod
  underflow); `eatof` narrows to float before its range check.
- **printf**: `%f` of a float prints the value promoted to double; NaN prints
  `nan` (never signed), infinities `inf`/`-inf`. `printf!` matches macOS libc
  for all of d i u o x X c s f e g E G with flags, width, precision, `*`, h/l.
- **cwp.h macros**: `MAX/MIN` are `x > y ? x : y`, not `f32::max` (NaN
  differs); `NINT(x)` adds 0.5 in double; `ABS(-0.0)` is `-0.0`; `SGN` is a
  double.
- **Complex arithmetic** in oxwave-cwp is clean-room: `cadd/csub/cmul/cdiv`
  are bit-identical to libcwp; `rcabs`, `cwp_csqrt`, `cwp_cexp` differ from
  the C in the last bit in some cases (<= 1.3e-7 relative), well inside the
  harness tolerance.
- **sumax** (and any program with `outpar=/dev/tty` default) needs
  `outpar=/dev/stdout` in test cases (no tty in the harness).
- **SIGPIPE**: `initargs` restores the default, so `prog | head` dies
  silently like C.
- **stdout buffering**: anything not flushed is lost if you `return` from
  `main` or `std::process::exit`. Use `exit()` from the prelude.
