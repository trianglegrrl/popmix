"""convertf's parameter files (TASKS.md P2): ``popmix convertf -p file.par``.

The keys, from convertf's documentation (README.CONVERTF), each checked
against AdmixTools 8 convertf on the same file
(tests/golden/test_convertf_par.py):

- ``genotypename``, ``snpname``, ``indivname``: the input files.
- ``outputformat``: ANCESTRYMAP, EIGENSTRAT, PED, PACKEDPED,
  PACKEDANCESTRYMAP (the default) or TRANSPOSE_PACKED.
- ``genotypeoutname``, ``snpoutname``, ``indivoutname`` (or mergeit's
  spellings ``genooutfilename``, ``snpoutfilename``, ``indoutfilename``):
  the output files.
- ``poplistname``: a file of the populations to keep.
- ``badsnpname``: a file of SNPs to leave out (first column; SNPs the data
  does not have are ignored).
- ``chrom``, with ``lopos`` and ``hipos``: one chromosome, or a region of it.

A value may start with a name defined earlier in the file and a ``/``
(``DIR: /data`` then ``genotypename: DIR/x.geno``), as in EIGENSOFT's
parameter files; ``#`` starts a comment. Any other key is an error rather
than silently ignored.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

from popmix import _core
from popmix.errors import InputError

FORMATS = {"ANCESTRYMAP": "ancestrymap", "EIGENSTRAT": "eigenstrat", "PED": "ped", "PACKEDPED": "packedped",
           "PACKEDANCESTRYMAP": "packedancestrymap", "TRANSPOSE_PACKED": "tgeno"}
OUT_KEYS = {"genotypeoutname": "geno", "genooutfilename": "geno", "snpoutname": "snp", "snpoutfilename": "snp",
            "indivoutname": "ind", "indoutfilename": "ind"}
IN_KEYS = {"genotypename": "geno", "snpname": "snp", "indivname": "ind"}
OTHER_KEYS = {"outputformat", "poplistname", "badsnpname", "chrom", "lopos", "hipos"}
#: popmix's output files by extension, as the role convertf names them
ROLE = {".geno": "geno", ".ancestrymapgeno": "geno", ".bed": "geno", ".ped": "geno", ".snp": "snp", ".bim": "snp",
        ".map": "snp", ".pedsnp": "snp", ".ind": "ind", ".fam": "ind", ".pedind": "ind"}


def read_par(path: str | os.PathLike[str]) -> tuple[dict[str, str], set[str]]:
    """The keys of a convertf parameter file, names substituted, and the
    names that were substituted."""
    out: dict[str, str] = {}
    names: set[str] = set()
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if ":" not in line:
            raise InputError(f"{path}, line {n}: expected 'name: value', got {line!r}")
        key, value = (x.strip() for x in line.split(":", 1))
        head, sep, rest = value.partition("/")
        if sep and head in out:
            value = f"{out[head]}/{rest}"
            names.add(head)
        out[key] = value
    return out, names


def _ids(path: str) -> list[str]:
    return [line.split()[0] for line in Path(path).read_text().splitlines() if line.split()]


def _input_prefix(par: dict[str, str], work: Path) -> str:
    files = {role: par[key] for key, role in IN_KEYS.items() if key in par}
    if len(files) < 3:
        missing = [k for k in IN_KEYS if k not in par]
        raise InputError(f"the parameter file needs {', '.join(missing)}")
    plink = files["geno"].endswith(".bed")
    exts = {"geno": ".bed", "snp": ".bim", "ind": ".fam"} if plink else {"geno": ".geno", "snp": ".snp", "ind": ".ind"}
    stems = {files[r][: -len(e)] if files[r].endswith(e) else None for r, e in exts.items()}
    if len(stems) == 1 and None not in stems:
        return stems.pop()
    # files of different names: linked under one prefix
    for role, ext in exts.items():
        os.symlink(os.path.abspath(files[role]), work / f"in{ext}")
    return str(work / "in")


def convertf_par(path: str | os.PathLike[str]) -> dict:
    """Run a convertf parameter file.

    Parameters
    ----------
    path
        The parameter file (see the module documentation for its keys).

    Returns
    -------
    dict
        As ``convertf``'s: ``n_samples``, ``n_snps``, ``input_format``,
        ``output_format``.
    """
    par, names = read_par(path)
    unknown = sorted(k for k in par if k not in IN_KEYS and k not in OUT_KEYS and k not in OTHER_KEYS
                     and k not in names)
    if unknown:
        raise InputError(f"{path}: popmix does not support convertf's {', '.join(unknown)}")
    fmt = par.get("outputformat", "PACKEDANCESTRYMAP").upper()
    if fmt not in FORMATS:
        raise InputError(f"{path}: unknown outputformat {fmt}; expected one of {', '.join(FORMATS)}")
    targets = {role: par[key] for key, role in OUT_KEYS.items() if key in par}
    if len(targets) < 3:
        raise InputError(f"{path}: needs genotypeoutname, snpoutname and indivoutname")
    subset: dict = {}
    if "poplistname" in par:
        subset["populations"] = _ids(par["poplistname"])
    if "badsnpname" in par:
        # SNPs listed that the data does not have are nothing to leave out
        snp = par["snpname"]
        present = {line.split()[1 if snp.endswith(".bim") else 0] for line in Path(snp).read_text().splitlines()
                   if line.split()}
        subset["exclude_snps"] = [s for s in _ids(par["badsnpname"]) if s in present]
    if "chrom" in par:
        if "lopos" in par or "hipos" in par:
            subset["region"] = (par["chrom"], int(par.get("lopos", 0)), int(par.get("hipos", 2**32 - 1)))
        else:
            subset["chromosomes"] = [par["chrom"]]
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        result = _core.convertf(_input_prefix(par, work), str(work / "out"), FORMATS[fmt], **subset)
        for produced in work.glob("out.*"):
            role = ROLE.get(produced.suffix)
            if role is not None:
                shutil.move(str(produced), targets[role])
    return result
