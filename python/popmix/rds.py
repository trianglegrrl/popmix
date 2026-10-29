"""Read and write R's ``.rds`` files (R serialization, XDR format), for the
objects admixtools2 stores in f2 directories: numeric and integer vectors
and matrices, character vectors, lists, NULL and attributes.

Read objects come back as ``RObject`` (data and attributes), or as the bare
data when there are no attributes: a list of floats, ints, strings or
objects. Missing values are None (character and integer) or NaN.

Provenance: the serialization format as the "R Internals" manual describes
it ("Serialization Formats"), checked by reading the f2 directories
admixtools2 2.0.10 writes and by admixtools2 reading the ones popmix writes
(tests/golden/test_f2_dirs_admixtools2.py), and by round trips
(tests/python/test_rds.py). Written before the provenance rule (ec93f28);
not compared with R's source.
"""

from __future__ import annotations

import gzip
import struct
from dataclasses import dataclass, field
from typing import Any

NILSXP, SYMSXP, LISTSXP, CHARSXP, LGLSXP, INTSXP, REALSXP, STRSXP, VECSXP = 0, 1, 2, 9, 10, 13, 14, 16, 19
ALTREP_SXP, NILVALUE_SXP, REFSXP = 238, 254, 255
NA_INTEGER = -(2**31)
NA_REAL_BITS = 0x7FF00000000007A2  # R's NA_real_
IS_OBJECT, HAS_ATTR, HAS_TAG = 1 << 8, 1 << 9, 1 << 10
UTF8_MASK = 1 << 3


@dataclass
class RObject:
    data: Any
    attributes: dict[str, Any] = field(default_factory=dict)


class _Reader:
    def __init__(self, buf: bytes):
        self.buf, self.pos, self.refs = buf, 0, []

    def take(self, n: int) -> bytes:
        b = self.buf[self.pos:self.pos + n]
        if len(b) < n:
            raise ValueError("truncated .rds data")
        self.pos += n
        return b

    def int(self) -> int:
        return struct.unpack(">i", self.take(4))[0]

    def length(self) -> int:
        n = self.int()
        if n == -1:
            hi, lo = self.int(), self.int()
            n = (hi << 32) + lo
        return n

    def string(self, flags: int) -> str | None:
        n = self.int()
        # R writes UTF-8 or ASCII here; latin-1 strings are not expected
        return None if n == -1 else self.take(n).decode("utf-8")

    def item(self) -> Any:
        flags = self.int()
        t = flags & 0xFF
        if t == NILVALUE_SXP:
            return None
        if t == REFSXP:
            i = flags >> 8
            return self.refs[(i if i else self.int()) - 1]
        if t == SYMSXP:
            self.int()  # CHARSXP flags
            name = self.string(0)
            self.refs.append(name)
            return name
        if t == ALTREP_SXP:
            info, state, attr = self.item(), self.item(), self.item()
            cls = info[0][1] if isinstance(info, _Pairlist) else info
            state = [v for _, v in state] if isinstance(state, _Pairlist) else state
            if cls == "compact_intseq":
                n, start, step = (int(x) for x in _values(state))
                data: Any = [start + k * step for k in range(n)]
            elif cls == "compact_realseq":
                n, start, step = _values(state)
                data = [start + k * step for k in range(int(n))]
            elif cls in ("wrap_integer", "wrap_real", "wrap_string", "wrap_logical", "wrap_list"):
                data = state[0] if isinstance(state, list) else state
                return data
            else:
                raise ValueError(f"unsupported ALTREP class {cls!r} in .rds data")
            return _with_attributes(data, _pairlist_dict(attr))
        if t == LISTSXP:
            out = _Pairlist()
            while t == LISTSXP:
                if flags & HAS_ATTR:
                    self.item()
                tag = self.item() if flags & HAS_TAG else None
                out.append((tag, self.item()))
                flags = self.int()
                t = flags & 0xFF
            if t != NILVALUE_SXP:
                raise ValueError("malformed pairlist in .rds data")
            return out
        if t == CHARSXP:
            return self.string(flags)
        if t in (INTSXP, LGLSXP):
            n = self.length()
            vals = struct.unpack(f">{n}i", self.take(4 * n))
            data = [None if v == NA_INTEGER else (bool(v) if t == LGLSXP else v) for v in vals]
        elif t == REALSXP:
            n = self.length()
            data = list(struct.unpack(f">{n}d", self.take(8 * n)))
        elif t == STRSXP:
            n = self.length()
            data = []
            for _ in range(n):
                f = self.int()
                data.append(self.string(f))
        elif t == VECSXP:
            n = self.length()
            data = [self.item() for _ in range(n)]
        elif t == NILSXP:
            return None
        else:
            raise ValueError(f"unsupported R type {t} in .rds data")
        attrs = _pairlist_dict(self.item()) if flags & HAS_ATTR else {}
        return _with_attributes(data, attrs)


class _Pairlist(list):
    """A pairlist as (tag, value) pairs."""


def _values(x: Any) -> list:
    return x.data if isinstance(x, RObject) else x


def _pairlist_dict(x: Any) -> dict:
    return dict(x) if isinstance(x, _Pairlist) else {}


def _with_attributes(data: Any, attrs: dict) -> Any:
    return RObject(data, attrs) if attrs else data


def read_rds(path: str) -> Any:
    """Read an ``.rds`` file (gzip-compressed or not, XDR format)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    if raw[:2] != b"X\n":
        raise ValueError(f"{path}: not an XDR .rds file")
    r = _Reader(raw[2:])
    version = r.int()
    r.int()  # writer's R version
    r.int()  # minimal reader version
    if version == 3:
        r.take(r.int())  # native encoding
    elif version != 2:
        raise ValueError(f"{path}: unsupported serialization version {version}")
    return r.item()


class _Writer:
    def __init__(self):
        self.parts: list[bytes] = []
        self.symbols: dict[str, int] = {}

    def int(self, v: int):
        self.parts.append(struct.pack(">i", v))

    def length(self, n: int):
        if n < 2**31 - 1:
            self.int(n)
        else:
            self.int(-1)
            self.int(n >> 32)
            self.int(n & 0xFFFFFFFF)

    def charsxp(self, s: str | None):
        if s is None:
            self.int(CHARSXP)
            self.int(-1)
            return
        b = s.encode("utf-8")
        ascii_ = all(c < 128 for c in b)
        self.int(CHARSXP | ((1 << 6 if ascii_ else UTF8_MASK) << 12))
        self.int(len(b))
        self.parts.append(b)

    def symbol(self, name: str):
        if name in self.symbols:
            self.int(REFSXP | (self.symbols[name] << 8))
            return
        self.symbols[name] = len(self.symbols) + 1
        self.int(SYMSXP)
        self.charsxp(name)

    def item(self, x: Any):
        attrs: dict = {}
        if isinstance(x, RObject):
            x, attrs = x.data, x.attributes
        if x is None:
            self.int(NILVALUE_SXP)
            return
        if isinstance(x, (str, int, float, bool)):
            x = [x]
        x = list(x)
        flag = HAS_ATTR if attrs else 0
        if all(isinstance(v, str) or v is None for v in x) and any(isinstance(v, str) for v in x):
            self.int(STRSXP | flag)
            self.length(len(x))
            for v in x:
                self.charsxp(v)
        elif x and all(isinstance(v, bool) for v in x):
            self.int(LGLSXP | flag)
            self.length(len(x))
            self.parts.append(struct.pack(f">{len(x)}i", *(int(v) for v in x)))
        elif x and all(isinstance(v, int) and not isinstance(v, bool) or v is None for v in x) and \
                any(isinstance(v, int) for v in x):
            self.int(INTSXP | flag)
            self.length(len(x))
            self.parts.append(struct.pack(f">{len(x)}i", *(NA_INTEGER if v is None else v for v in x)))
        elif all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in x):
            self.int(REALSXP | flag)
            self.length(len(x))
            self.parts.append(struct.pack(f">{len(x)}d", *(float(v) for v in x)))
        else:
            self.int(VECSXP | flag)
            self.length(len(x))
            for v in x:
                self.item(v)
        if attrs:
            for name, value in attrs.items():
                self.int(LISTSXP | HAS_TAG)
                self.symbol(name)
                self.item(value)
            self.int(NILVALUE_SXP)


def write_rds(path: str, obj: Any, compress: bool = True):
    """Write ``obj`` as an ``.rds`` file R reads with ``readRDS``: lists of
    numbers become numeric vectors (integer if every value is an int),
    lists of strings character vectors, other lists R lists; ``RObject``
    adds attributes (``dim``, ``dimnames``, ``names``...)."""
    w = _Writer()
    w.parts.append(b"X\n")
    w.int(3)
    w.int(0x040502)
    w.int(0x030500)
    w.int(5)
    w.parts.append(b"UTF-8")
    w.item(obj)
    raw = b"".join(w.parts)
    with open(path, "wb") as fh:
        fh.write(gzip.compress(raw, mtime=0) if compress else raw)


def matrix(values_by_column: list[list[float]], colnames: list[str] | None = None) -> RObject:
    """A numeric R matrix from its columns."""
    nrow = len(values_by_column[0]) if values_by_column else 0
    attrs: dict[str, Any] = {"dim": [nrow, len(values_by_column)]}
    if colnames is not None:
        attrs["dimnames"] = [None, list(colnames)]
    return RObject([v for col in values_by_column for v in col], attrs)
