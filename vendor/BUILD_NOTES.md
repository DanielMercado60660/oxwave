# CWP/SU reference build notes

Build date: 2026-10-28.

- Host: Apple M5 Max (arm64), macOS (Darwin 27.0.0).
- Compiler: Apple clang 21.0.0 (clang-2100.3.34.2), invoked as `/usr/bin/gcc`.
- Make: GNU Make 3.81 (Xcode).
- Source: `vendor/su-original`, a clone of github.com/JohnWStockwellJr/SeisUnix at commit `4db4181`. The clone was not modified (`git status` is clean).
- Build tree: `vendor/su-build`, a copy of the clone without `.git`.
- `CWPROOT=vendor/su-build`
- Result: **success**. `make -k install` exited 0, and the full log has 0 `error:` lines. The log is in `su-build/build.log`.

## Read this first: SU was built with XDR, so .su files are big-endian

The Big Sur config sets `XDRFLAG = -DSUXDR` and `ENDIANFLAG = -DCWP_LITTLE_ENDIAN`. Both were left unchanged. As a result:

- **Every `.su` stream these programs read or write is big-endian (XDR) on this little-endian arm64 host.** This covers the 240-byte headers and the float samples. Checked on `suplane | sugain scale=2`: header `offset` reads as `00 00 01 90` (400) and `ns` as `00 40` (64), and trace 6's samples decode as `2.0` only when read as big-endian float32.
- In SUXDR mode `fgettr` reads **only** XDR/big-endian data. It does not auto-detect native-endian input. Tested: `suplane | suswapbytes > /tmp/le.su; surange < /tmp/le.su` prints `fgettr.c: no data on first trace` and reads nothing.
- OxWave golden-file comparisons must treat `.su` as big-endian, or convert with care.
- `sustrip` and `supaste` work on raw native floats. `sustrip` output was little-endian on this host: `2.0` appears as `00 00 00 40`.
- `LINEHDRFLAG` is empty, so there are no 3200/400-byte SEG-Y line headers on `.su` streams.
- `LARGE_FILE_FLAG` is empty. On 64-bit macOS `off_t` is 64-bit anyway.

## Exact steps

```sh
cd vendor
rsync -a --exclude=.git su-original/ su-build/
cd su-build/src
cp Makefile.config_MacOSX_Big_Sur Makefile.config
# edit OPTC in Makefile.config (diff below)
# skip X11 xplot inside su/graphics (diff below)

# license / mail-home handling (see section below)
touch LICENSE_45R00_ACCEPTED MAILHOME_45R00
cat > mailhome.sh <<'EOF'      # stub, belt-and-braces
#! /bin/sh
# OxWave build copy: mail-home step disabled (no email/network).
echo "mailhome.sh disabled in this build copy; no mail sent."
exit 0
EOF
echo "CWP/SU Release 45R00 30 April 2023" > cwp_su_version   # what license.sh would write; not compiled in

export CWPROOT=vendor/su-build
echo y | make -k install > ../build.log 2>&1      # "y" answers chkroot.sh's "Continue install?"
```

Only the `install` target was run, not `xtinstall`, `xminstall`, `mglinstall`, `finstall`, `sfinstall` or `utils`. Nothing was installed system-wide; sudo was not used and no shell rc file was touched.

## License and mail-home prompts

- The `LICENSE_44R25_ACCEPTED` and `LICENSE_44R18_ACCEPTED` files **do not** skip the prompt. The top-level `Makefile` `install` target depends on the marker targets `LICENSE_45R00_ACCEPTED` (runs `license.sh`) and `MAILHOME_45R00` (runs `mailhome.sh`), and the 44R* files have different names.
- Both 45R00 markers were created with `touch`, so neither script ran.
- `mailhome.sh` would mail the SU maintainers (`echo ... | mail <address>`) if it got a "y". In the build copy it was also replaced with a no-op stub. No mail or network request was sent.
- No other Makefile or script in the `install` path calls `mail`, `curl` or `wget`. The only hits are in `BHP_SU_V7.0/dev/*`, which `install` does not build.
- `chkroot.sh` still asks "Continue install? [y/n]". It got its "y" from `echo y |`.

## Makefile.config changes

These are against `src/Makefile.config_MacOSX_Big_Sur`:

```diff
--- su-original/src/Makefile.config_MacOSX_Big_Sur
+++ su-build/src/Makefile.config
@@ -91,7 +91,7 @@
 
 CC = /usr/bin/gcc
 #OPTC = -g -c90  -m64 -Wall -ansi  -Wno-long-long 
-OPTC = -g -std=c90  -m64 -Wall -ansi  -Wno-long-long 
+OPTC = -g -std=c90  -m64 -Wall -ansi  -Wno-long-long -ffp-contract=off
 CFLAGS = -I$I $(OPTC) $(CWP_FLAGS)
```

Other files changed in the build copy:

```diff
--- su-original/src/su/graphics/Makefile
+++ su-build/src/su/graphics/Makefile
 INSTALL:
-	cd xplot ; $(MAKE)
+	@echo "su/graphics/xplot skipped (X11, not built for OxWave judge)"
 	cd psplot ; $(MAKE)
 remake:
-	cd xplot ; $(MAKE) remake
+	@echo "su/graphics/xplot skipped"
 	cd psplot ; $(MAKE) remake
```

`src/mailhome.sh` was replaced by the stub. `src/cwp_su_version` was rewritten as `license.sh` would. `src/LICENSE_45R00_ACCEPTED` and `src/MAILHOME_45R00` were added.

## Final CC / CFLAGS

```
CC     = /usr/bin/gcc            (Apple clang 21.0.0, arm64-apple-darwin27.0.0)
OPTC   = -g -std=c90 -m64 -Wall -ansi -Wno-long-long -ffp-contract=off
CFLAGS = -I$CWPROOT/include $(OPTC) -DCWP_LITTLE_ENDIAN -DSUXDR
```

This is the actual compile line seen in `build.log`:

```
/usr/bin/gcc -I.../su-build/include -g -std=c90 -m64 -Wall -ansi -Wno-long-long -ffp-contract=off -DCWP_LITTLE_ENDIAN -DSUXDR  prog.c -L.../su-build/lib -lsu -lpar -lcwp -lm -o .../su-build/bin/prog
```

- **Optimization:** the Big Sur config has **no `-O` flag**, so everything is built at clang's default **-O0**. It was kept that way: the reference build uses whatever the config uses. `-ffp-contract=off` also guarantees that no FMA fusion happens if someone later adds `-O2`.
- **Language mode:** `-ansi` comes after `-std=c90`, so it wins and the C mode is C89/C90. Clang 21 compiled everything in that mode, so the Monterey config's switch to `-std=c99` was not needed.
- The tri, tetra, wpc and other sub-builds take `$(OPTC)` through `Rules/cflags.rules`, so they got `-ffp-contract=off` too.

## What was built

- **`su-build/bin`: 455 programs.**
  - 400 Mach-O arm64 executables.
  - 55 shell scripts (from `su/shell`, e.g. `suagc`, `sudoc`, `suhelp`, `sukeyword`).
  - bin also holds 399 `*.dSYM` debug-symbol directories. These come from `-g` with compile-and-link in one step, are harmless, and can be ignored or deleted.
- **`su-build/lib`:** `libcwp.a`, `libpar.a`, `libsu.a`, plus `libpsplot.a`, `libcomp.a`, `libwpc.a`, `libwpc1.a`, `libreflect.a`, `libtriang.a`, `libtrielas.a` and `libtetra.a`.
- **su/main:** 284 of 286 `.c` sources became binaries in `bin/`.
- **Other directories built:** `cwp`, `par`, `psplot` (PostScript plotting, no X), `su/{include,lib,main,shell,graphics/psplot}`, `tri`, `Trielas`, `tetra`, `comp`, `Refl`.

### su/main programs not built (one line each)

- `su2voxet`: `su/main/gocad/` has no Makefile and is not listed in `su/main/Makefile`, so it is not part of the standard install. It did not fail.
- `sudctcomp`: the source is in `su/main/data_compression/` but is not in that Makefile's `PROGS`. The equivalent `dctcomp` is built from `src/comp` and is in `bin/`. It did not fail.

### Skipped on purpose (X11/Motif/GL/Fortran/PVM)

- `su/graphics/xplot`: `suxwigb`, `suximage`, `suxgraph`, `suxmovie`, `suxpicker`, `suxcontour` and similar.
- The `Xtcwp`, `Xmcwp`, `xplot`, `xtri`, `Mesa`, `Fortran`, `PVM`, `tksu`, `Sfio`/`segdread` and `cwputils` targets were not run.
- X11 is not present on this Mac anyway: `/usr/X11` is a dangling symlink and `/opt/X11` does not exist.

No non-X part failed.

## Verification

```sh
export CWPROOT=vendor/su-build
export PATH=$CWPROOT/bin:$PATH
suplane | sugain scale=2 | surange
#   32 traces:
#   tracl    1 32 (1 - 32)
#   tracr    1 32 (1 - 32)
#   offset   400
#   ns       64
#   dt       4000
suplane > /tmp/x.su          # 15872 bytes = 32 traces x (240 + 64*4)
suplane | sumax outpar=/dev/stdout                 # 3.000000e+00 0.000000e+00
suplane | sugain scale=2 | sumax outpar=/dev/stdout  # 6.000000e+00 0.000000e+00
```

**Selfdoc:** when `sugain` is run with no arguments **on a terminal**, it prints its selfdoc ("SUGAIN - apply various types of gain ...") through a pager. Checked under a pty with `script`.

## Surprises and gotchas for using this as a test oracle

1. **XDR/big-endian `.su` files.** This matters most; details are in the first section.
2. **No optimization (-O0).** The config has no `-O` flag. Results are still deterministic, and `-ffp-contract=off` means no FMA in either case.
3. **Selfdoc only shows when stdin is a tty.** With stdin redirected or piped (for example `sugain </dev/null`, or in a harness), the program reads stdin instead. It fails with `fgettr.c: bad first header` and exits 1. To capture selfdoc in automation, use `sudoc sugain` (tested; it reads `src/` in su-build). Or run under a pty with a non-interactive pager: `PAGER=cat script -q /dev/null sugain`. The selfdoc code in `par/lib/docpkge.c` honours `$PAGER`, and the default pager waits for input.
4. **`sumax` defaults to `outpar=/dev/tty`.** Without a controlling terminal it fails with `efopen: fopen failed`. Pass `outpar=/dev/stdout` in harnesses. Other programs with an `outpar=` parameter may do the same.
5. **Compiler errors don't stop the build.** su/main Makefiles run `-$(CC)` (the leading `-` means errors are ignored). A "successful" make is therefore not proof that everything compiled. It was checked via `grep error: build.log` (0 hits) and by comparing sources against `bin/` (above).
6. **103 compiler warnings.** These may matter for reproducibility in the affected programs:
   - `-Wsometimes-uninitialized`:
     - `headcase.c` (`dval`)
     - `sutetraray.c` (`ifwf2dump`, `nwf2dump`)
   - `-Wliteral-range` (float comparison that is always false):
     - `upweik.c`
     - `suinvvxzco.c`
     - `susynvxzcs.c`
   - Constant int→float conversions:
     - `rayt2dan.c`
     - `suharlan.c`
     - `sucliphead.c`
     - `sushw.c`
   - The rest are unused-variable, misleading-indentation and deprecated-non-prototype warnings.
7. **Version note:** `cwp_su_version` in the clone said 44R24, while the Makefile and scripts are 45R00. The source tree is the 45R00 (30 April 2023) release plus later GitHub commits up to `4db4181`.
