# Part of OxWave, a Rust port of CWP/SU.
"""Run OxWave programs: binary lookup, argv encoding, run/pipe, selfdoc, catalogue.

Contract (stable signatures/docstrings).

Programs run as real subprocesses with argv lists (never a shell), stdin/stdout
piped, so results are exactly what the command line produces. The environment is
``os.environ`` plus ``LC_ALL=C`` plus the caller's ``env`` overrides.

Binary lookup (first hit wins): directories from :func:`set_bin_dir`; else
``$OXWAVE_BIN`` (``os.pathsep``-separated); else ``<package dir>/bin`` if it exists.
``$PATH`` is never searched. A program is *available* if an executable regular file
with its exact name exists in one of those directories.

The catalogue is the bundled ``oxwave/_programs.json`` (``{"programs": [entry...]}``;
entry keys: name, purpose, family, source, stdin, stdout, params [{name, type, default,
required, description}], selfdoc_first_line, related_programs, uses_randomness, status).
``stdout`` is one of ``su_traces``, ``ascii``, ``binary_floats``, ``other``, ``none``.

Argument encoding (``format_argv``): ``argv = [prog, *args, *encoded params]`` where
params are emitted in insertion order as ``key=value``:

* key must match ``[A-Za-z_][A-Za-z0-9_]*`` (else ValueError);
* ``None`` -> parameter omitted; ``bool`` -> TypeError (write 0/1 explicitly);
* int / numpy integer -> ``str(int(v))``; float / numpy floating -> ``repr(float(v))``
  (non-finite -> ValueError); ``str`` -> as is (NUL -> ValueError); ``os.PathLike`` ->
  ``os.fspath``;
* list / tuple / 1-D numpy array -> elements encoded as above, joined by ``,``
  (nested sequences -> TypeError);
* :class:`Repeat` ``([v1, v2, ...])`` -> the key repeated, ``key=v1 key=v2 ...``
  (each ``v`` encoded as above, so a ``v`` may itself be a list).

``args`` are extra positional argv words placed right after the program name (e.g.
file names for ``suop2``); each must be a ``str``/PathLike and must not start with ``-``.
"""
import copy
import errno
import json
import math
import os
import pty
import re
import select
import signal
import subprocess
import threading
import time
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from .traces import Traces

__all__ = [
    "Repeat", "RunResult", "OxWaveError", "ProgramNotFound", "ProgramError", "ProgramTimeout",
    "OutputLimitExceeded", "UnknownParameterError", "catalog", "program_info", "bin_dirs",
    "set_bin_dir", "which", "available_programs", "format_argv", "validate_params", "run",
    "pipe", "selfdoc", "help",
]


with Path(__file__).with_name("_programs.json").open(encoding="utf-8") as _catalog_file:
    _CATALOG = {entry["name"]: entry for entry in json.load(_catalog_file)["programs"]}
_BIN_DIRS = None


def _program_name(name):
    if (not isinstance(name, str) or not name or name in (".", "..")
            or any(c in name for c in ("/", "\\", "\0"))):
        raise ValueError(f"program must be a bare name, got {name!r}")
    return name


def _word(value):
    if not isinstance(value, (str, os.PathLike)):
        raise TypeError("argv words must be strings or PathLike objects")
    word = os.fsdecode(os.fspath(value))
    if "\0" in word:
        raise ValueError("argv words cannot contain NUL")
    return word


def _scalar(value):
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        raise TypeError("boolean parameters are ambiguous; write 0 or 1 explicitly")
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError("parameters must be finite")
        return repr(number)
    if isinstance(value, (str, os.PathLike)):
        return _word(value)
    raise TypeError(f"unsupported parameter value: {value!r}")


def _value(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        if isinstance(value, np.ndarray) and value.ndim != 1:
            raise TypeError("parameter arrays must be one-dimensional")
        parts = [_scalar(item) for item in value]
        return ",".join(part for part in parts if part is not None)
    return _scalar(value)


def _environment(overrides=None):
    env = dict(os.environ, LC_ALL="C")
    if overrides is not None:
        env.update(overrides)
    return env


def _file_path(path, cwd):
    path = Path(path)
    return Path(cwd) / path if cwd is not None and not path.is_absolute() else path


def _deadline(timeout, started):
    if timeout is None:
        return None
    timeout = float(timeout)
    if not math.isfinite(timeout):
        raise ValueError("timeout must be finite")
    return started + timeout


def _remaining(deadline):
    if deadline is None:
        return None
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise ProgramTimeout("OxWave run exceeded its timeout")
    return remaining


def _kill(processes):
    for process in processes:
        if process.poll() is None:
            try:
                process.kill()
            except ProcessLookupError:
                pass


class Repeat(list):
    """Marks a parameter whose key is repeated: ``ref=Repeat(["1:1,2;4,2", "1:0,3;4,3"])``
    -> ``ref=1:1,2;4,2 ref=1:0,3;4,3``."""


class OxWaveError(RuntimeError):
    """Base class of all errors raised by this module."""


class ProgramNotFound(OxWaveError, FileNotFoundError):
    """No binary for the program in the lookup directories. Message names the program,
    the directories searched, and mentions ``OXWAVE_BIN``."""


class ProgramError(OxWaveError):
    """A stage exited non-zero (with ``check=True``).

    Attributes: ``program`` (name of the failing stage), ``stage`` (0-based index),
    ``returncode`` (int; negative = killed by that signal), ``stderr`` (that stage's
    stderr, str), ``argv`` (that stage's argv list). ``str(e)`` contains the program name,
    the return code and the first non-empty stderr line.
    A non-last stage killed by SIGPIPE is not a failure (downstream stopped reading).
    """

    def __init__(self, message, *, program=None, stage=0, returncode=None, stderr="", argv=None):
        super().__init__(message)
        self.program = program
        self.stage = stage
        self.returncode = returncode
        self.stderr = stderr
        self.argv = argv or []


class ProgramTimeout(OxWaveError, TimeoutError):
    """The run exceeded ``timeout`` seconds; every stage was killed."""


class OutputLimitExceeded(OxWaveError):
    """Final stdout exceeded ``max_output_bytes``; every stage was killed and a partial
    output file (if ``output`` was a path) was deleted."""


class UnknownParameterError(OxWaveError, ValueError):
    """A parameter name is not in the program's catalogue entry. ``str(e)`` names the
    program and the parameter."""


class RunResult:
    """Outcome of :func:`run` / :func:`pipe` when the result is not returned as Traces.

    Attributes: ``argv`` (list of argv lists, one per stage), ``returncodes`` (list of
    int), ``returncode`` (the last stage's), ``stdout`` (bytes; ``b""`` when written to a
    file), ``stderrs`` (list of str, one per stage), ``stderr`` (all stage stderrs
    concatenated in stage order), ``output`` (Path or None), ``elapsed`` (seconds).
    """

    def __init__(self, argv, returncodes, stdout, stderrs, output=None, elapsed=0.0):
        self.argv = argv
        self.returncodes = list(returncodes)
        self.returncode = self.returncodes[-1] if self.returncodes else None
        self.stdout = stdout
        self.stderrs = list(stderrs)
        self.stderr = "".join(self.stderrs)
        self.output = output
        self.elapsed = elapsed

    def __repr__(self):
        return "RunResult(returncodes=%r, stdout=%d bytes)" % (self.returncodes, len(self.stdout or b""))


def catalog():
    """The bundled catalogue as ``{name: entry}`` (a fresh deep copy each call)."""
    return copy.deepcopy(_CATALOG)


def program_info(name):
    """Catalogue entry of ``name`` (deep copy); KeyError if not in the catalogue."""
    return copy.deepcopy(_CATALOG[name])


def bin_dirs():
    """The lookup directories currently in effect, as a list of ``Path`` (see module doc)."""
    if _BIN_DIRS is not None:
        return list(_BIN_DIRS)
    value = os.environ.get("OXWAVE_BIN")
    if value:
        return [Path(directory) for directory in value.split(os.pathsep) if directory]
    bundled = Path(__file__).parent / "bin"
    return [bundled] if bundled.is_dir() else []


def set_bin_dir(path):
    """Override the lookup: a directory, a sequence of directories, or None to go back
    to ``$OXWAVE_BIN`` / bundled lookup. Process-wide."""
    global _BIN_DIRS
    if path is None:
        _BIN_DIRS = None
    elif isinstance(path, (str, os.PathLike)):
        _BIN_DIRS = [Path(path)]
    else:
        _BIN_DIRS = [Path(directory) for directory in path]


def which(name):
    """Absolute ``Path`` of the binary for ``name``; ProgramNotFound if not available.
    ``name`` must be a bare program name (no path separators) -> else ValueError."""
    _program_name(name)
    directories = bin_dirs()
    for directory in directories:
        binary = directory / name
        if binary.is_file() and os.access(binary, os.X_OK):
            return binary.absolute()
    searched = ", ".join(map(str, directories)) or "(no directories)"
    raise ProgramNotFound(f"{name}: executable not found in {searched}; set OXWAVE_BIN")


def available_programs():
    """Sorted names of catalogue programs that are available (binary found)."""
    names = []
    for name in sorted(_CATALOG):
        try:
            which(name)
        except ProgramNotFound:
            continue
        names.append(name)
    return names


def format_argv(prog, params=None, args=()):
    """argv list ``[prog, *args, *key=value]`` per the module's encoding rules.
    ``params``: a mapping (insertion order kept) or None. Does not look up binaries."""
    argv = [_program_name(prog)]
    if isinstance(args, (str, bytes, os.PathLike)):
        raise TypeError("args must be a sequence of argv words")
    for arg in args:
        word = _word(arg)
        if word.startswith("-"):
            raise ValueError("positional arguments must not start with '-'")
        argv.append(word)
    if params is None:
        return argv
    if not isinstance(params, Mapping):
        raise TypeError("params must be a mapping or None")
    for key, value in params.items():
        if not isinstance(key, str) or re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key) is None:
            raise ValueError(f"invalid parameter name: {key!r}")
        values = value if isinstance(value, Repeat) else (value,)
        for item in values:
            encoded = _value(item)
            if encoded is not None:
                argv.append(f"{key}={encoded}")
    return argv


def validate_params(prog, params):
    """Raise UnknownParameterError for any key of ``params`` not in ``prog``'s catalogue
    ``params`` (``par`` is always allowed). Programs whose catalogue entry has an empty
    ``params`` list, or that are not in the catalogue, accept anything. Returns None."""
    entry = _CATALOG.get(prog)
    if entry is None or not entry["params"] or params is None:
        return
    allowed = {param["name"] for param in entry["params"]} | {"par"}
    for key in params:
        if key not in allowed:
            raise UnknownParameterError(f"{prog}: unknown parameter {key!r}")


def run(prog, input=None, *, output=None, args=(), params=None, timeout=None, check=True,
        validate=False, raw=False, cwd=None, env=None, max_output_bytes=None, **kwparams):
    """Run one program; same as ``pipe([(prog, merged_params, args)], ...)``.

    Parameters are ``params`` (mapping) followed by ``**kwparams`` (a key in both ->
    TypeError). Use ``params={"lambda": 1}`` for names that are Python keywords.
    See :func:`pipe` for ``input``, ``output``, return value and errors.
    """
    if params is not None and not isinstance(params, Mapping):
        raise TypeError("params must be a mapping or None")
    merged = dict(params) if params is not None else {}
    for key, value in kwparams.items():
        if key in merged:
            raise TypeError(f"parameter {key!r} supplied in both params and keyword arguments")
        merged[key] = value
    return pipe([(prog, merged, args)], input, output=output, timeout=timeout, check=check,
                validate=validate, raw=raw, cwd=cwd, env=env, max_output_bytes=max_output_bytes)


def pipe(stages, input=None, *, output=None, timeout=None, check=True, validate=False, raw=False,
         cwd=None, env=None, max_output_bytes=None):
    """Run a pipeline ``stage0 | stage1 | ...`` of OxWave programs.

    ``stages``: non-empty sequence; each stage is ``"prog"``, ``(prog,)``,
    ``(prog, params)`` or ``(prog, params, args)`` with ``params`` a mapping or None.
    All binaries are resolved (ProgramNotFound) and, if ``validate``, parameters
    checked (UnknownParameterError) before anything starts.

    ``input`` (stdin of stage 0): None -> no stdin (``/dev/null``); :class:`Traces` ->
    its SU bytes; bytes-like -> those bytes; ``str``/PathLike -> that file opened for
    reading (relative to ``cwd`` if given).
    ``output``: None -> stdout captured; ``str``/PathLike -> stdout streamed into that
    file (created/truncated; relative to ``cwd`` if given).
    ``cwd``: working directory of every stage. ``timeout``: seconds for the whole run
    (ProgramTimeout). ``max_output_bytes``: limit on final stdout (OutputLimitExceeded).

    Returns: if ``output`` is None, ``raw`` is False and the last program's catalogue
    ``stdout`` is ``su_traces``: ``Traces.from_bytes(stdout)`` (empty stdout -> a Traces
    with 0 traces and ns 0). Otherwise a :class:`RunResult`.
    Errors: with ``check=True`` the first failing stage (see ProgramError) raises
    ProgramError; with ``check=False`` a RunResult is returned instead of raising for
    exit codes (and is also returned instead of Traces when any stage failed).
    The stdout bytes are identical to ``bash -o pipefail -c "prog0 ... | prog1 ..."``
    run with the same binaries (stdout only: bash reports 141 for an upstream stage
    killed by SIGPIPE, which this function does not treat as a failure).
    """
    started = time.monotonic()
    deadline = _deadline(timeout, started)
    if max_output_bytes is not None and (
        not isinstance(max_output_bytes, (int, np.integer)) or max_output_bytes < 0
    ):
        raise ValueError("max_output_bytes must be a non-negative integer")
    if isinstance(stages, (str, bytes)):
        raise TypeError("stages must be a sequence of stages")
    argvs, binaries = [], []
    for stage in stages:
        if isinstance(stage, str):
            prog, params, args = stage, None, ()
        elif isinstance(stage, (tuple, list)) and 1 <= len(stage) <= 3:
            prog = stage[0]
            params = stage[1] if len(stage) >= 2 else None
            args = stage[2] if len(stage) == 3 else ()
        else:
            raise ValueError("a stage must be a program name or a tuple of 1 to 3 items")
        argv = format_argv(prog, params, args)
        binary = which(prog)
        if validate:
            validate_params(prog, params)
        argvs.append(argv)
        binaries.append(binary)
    if not argvs:
        raise ValueError("a pipeline must contain at least one stage")

    payload = None
    input_path = None
    if isinstance(input, Traces):
        payload = input.to_bytes()
    elif isinstance(input, (str, os.PathLike)):
        input_path = _file_path(input, cwd)
    elif input is not None:
        try:
            payload = memoryview(input).tobytes()
        except TypeError:
            raise TypeError("input must be Traces, bytes-like, a path, or None") from None
    output_path = None if output is None else _file_path(output, cwd)
    child_env = _environment(env)
    processes, threads = [], []
    stderr_chunks = [bytearray() for _ in argvs]
    stdout_chunks = []
    errors = []
    wake = threading.Event()
    input_file = output_file = None
    remove_output = False

    def worker(function, *args):
        try:
            function(*args)
        except Exception as error:  # noqa: BLE001 - re-raised in the calling thread
            errors.append(error)
        finally:
            wake.set()

    def start_worker(function, *args):
        thread = threading.Thread(target=worker, args=(function, *args), daemon=True)
        thread.start()
        threads.append(thread)

    def read_stderr(stream, chunks):
        with stream:
            while True:
                chunk = stream.read(65536)
                if not chunk:
                    return
                chunks.extend(chunk)

    def write_stdin(stream):
        with stream:
            view = memoryview(payload)
            offset = 0
            try:
                while offset < len(view):
                    offset += stream.write(view[offset:offset + 65536])
            except BrokenPipeError:
                pass  # A consumer may deliberately stop before all input is written.

    def read_stdout(stream):
        total = 0
        with stream:
            while True:
                chunk = stream.read(65536)
                if not chunk:
                    return
                total += len(chunk)
                if max_output_bytes is not None and total > max_output_bytes:
                    raise OutputLimitExceeded(f"final stdout exceeded {max_output_bytes} bytes")
                if output_file is None:
                    stdout_chunks.append(chunk)
                else:
                    output_file.write(chunk)

    try:
        _remaining(deadline)
        if input_path is not None:
            input_file = input_path.open("rb")
        if output_path is not None:
            output_file = output_path.open("wb")
        first_stdin = (input_file if input_file is not None else
                       subprocess.PIPE if payload is not None else subprocess.DEVNULL)
        for index, (argv, binary) in enumerate(zip(argvs, binaries)):
            _remaining(deadline)
            process = subprocess.Popen(
                argv, executable=str(binary),
                stdin=first_stdin if index == 0 else processes[-1].stdout,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd, env=child_env,
                close_fds=True, restore_signals=True, bufsize=0,
            )
            processes.append(process)
            # Only the downstream child may retain this read end. Otherwise the
            # upstream child cannot receive SIGPIPE when the downstream exits.
            if index:
                processes[index - 1].stdout.close()
            start_worker(read_stderr, process.stderr, stderr_chunks[index])
            if index == 0 and payload is not None:
                start_worker(write_stdin, process.stdin)
        start_worker(read_stdout, processes[-1].stdout)

        while True:
            finished = (all(process.poll() is not None for process in processes)
                        and all(not thread.is_alive() for thread in threads))
            # Check errors after checking thread completion, so a reader that
            # fails just as the last process exits cannot lose its exception.
            if errors:
                raise errors[0]
            if finished:
                break
            remaining = _remaining(deadline)
            wake.wait(0.01 if remaining is None else min(0.01, remaining))
            wake.clear()
    except OutputLimitExceeded:
        remove_output = True
        raise
    finally:
        # This also covers a launch failure midway through the pipeline, worker
        # exceptions, and interruptions in the caller.
        _kill(processes)
        for process in processes:
            process.wait()
        for thread in threads:
            thread.join()
        for process in processes:
            for stream in (process.stdin, process.stdout, process.stderr):
                if stream is not None:
                    stream.close()
        if input_file is not None:
            input_file.close()
        if output_file is not None:
            try:
                output_file.close()
            finally:
                if remove_output:
                    output_path.unlink(missing_ok=True)

    returncodes = [process.returncode for process in processes]
    stderrs = [bytes(chunks).decode("latin-1") for chunks in stderr_chunks]
    failing = [index for index, code in enumerate(returncodes)
               if code != 0 and not (index < len(processes) - 1 and code == -signal.SIGPIPE)]
    if check and failing:
        index = failing[0]
        first_line = next((line.strip() for line in stderrs[index].splitlines() if line.strip()), "")
        message = f"{argvs[index][0]} exited with return code {returncodes[index]}"
        if first_line:
            message += ": " + first_line
        raise ProgramError(message, program=argvs[index][0], stage=index,
                           returncode=returncodes[index], stderr=stderrs[index], argv=argvs[index])
    stdout = b"".join(stdout_chunks)
    if (output is None and not raw and not failing
            and _CATALOG.get(argvs[-1][0], {}).get("stdout") == "su_traces"):
        return Traces.from_bytes(stdout) if stdout else Traces(np.zeros((0, 0), np.float32))
    return RunResult(argvs, returncodes, stdout, stderrs, output=output_path,
                     elapsed=time.monotonic() - started)


def selfdoc(prog, timeout=30):
    """The program's self-documentation exactly as SU prints it: run the binary with no
    arguments, stdin and stdout on a pseudo-terminal, ``PAGER=cat``, and return what it
    writes to stderr (decoded as latin-1). ProgramNotFound if not available."""
    started = time.monotonic()
    deadline = _deadline(timeout, started)
    binary = which(prog)
    master, slave = pty.openpty()
    process = None
    stderr = bytearray()
    try:
        _remaining(deadline)
        process = subprocess.Popen(
            [prog], executable=str(binary), stdin=slave, stdout=slave, stderr=subprocess.PIPE,
            env=_environment({"PAGER": "cat"}), start_new_session=True,
            close_fds=True, restore_signals=True, bufsize=0,
        )
        os.close(slave)
        slave = None
        err_fd = process.stderr.fileno()
        descriptors = [master, err_fd]
        while descriptors:
            remaining = _remaining(deadline)
            ready, _, _ = select.select(
                descriptors, [], [], 0.1 if remaining is None else min(0.1, remaining),
            )
            for fd in ready:
                try:
                    chunk = os.read(fd, 65536)
                except OSError as error:
                    if fd != master or error.errno != errno.EIO:
                        raise
                    chunk = b""  # A closed pty reports EIO on some systems.
                if not chunk:
                    descriptors.remove(fd)
                elif fd == err_fd:
                    stderr.extend(chunk)
            if process.poll() is not None and err_fd not in descriptors:
                break
        try:
            process.wait(timeout=_remaining(deadline))
        except subprocess.TimeoutExpired:
            raise ProgramTimeout(f"{prog}: selfdoc exceeded its timeout") from None
        return bytes(stderr).decode("latin-1")
    finally:
        if process is not None:
            _kill([process])
            process.wait(timeout=1.0)
            process.stderr.close()
        if slave is not None:
            os.close(slave)
        os.close(master)


def help(prog):
    """Human help: :func:`selfdoc` when the binary is available; otherwise a text built
    from the catalogue (name, purpose, status, parameters with defaults). KeyError if the
    program is neither available nor in the catalogue."""
    try:
        return selfdoc(prog)
    except ProgramNotFound:
        entry = program_info(prog)
    lines = [f"{entry['name']} - {entry['purpose']}", "Status: " + entry["status"]]
    if entry["params"]:
        lines.append("Parameters:")
        for param in entry["params"]:
            required = " (required)" if param["required"] else ""
            lines.append(f"  {param['name']}={param['default']}{required}: {param['description']}")
    return "\n".join(lines) + "\n"
