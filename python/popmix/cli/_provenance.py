"""A provenance file beside every result a command writes (TASKS.md H3).

After a command succeeds, its main output (``--output``, ``--output-prefix``
or the f2 directory; else ``--output-graph``, ``--plot``, ``--heatmap`` or
``--plot-fit``) gets ``<output>.provenance.json`` (inside an f2 directory,
``provenance.json``): the popmix version, the command and its arguments,
every parameter with defaults filled in, and the size and modification time
of each input. Output on stdout gets none.
"""

from __future__ import annotations

import functools
import json
import os
import platform
import sys
from datetime import datetime, timezone
from typing import Any

import click

import popmix as pm

#: Output parameters, the main one first.
OUTPUTS = ("output", "output_prefix", "output_dir", "output_graph", "plot", "heatmap", "plot_fit")
#: Files that make up a genotype prefix.
PREFIX_FILES = (".geno", ".snp", ".ind", ".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam", ".vcf", ".vcf.gz")


def _jsonable(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _stat(path: str) -> dict[str, Any]:
    st = os.stat(path)
    return {"size": st.st_size,
            "modified": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds")}


def _inputs(params: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Files and directories the parameters name, outputs aside."""
    found: dict[str, dict[str, Any]] = {}
    for name, value in params.items():
        if name in OUTPUTS:
            continue
        for v in value if isinstance(value, (list, tuple)) else [value]:
            if not isinstance(v, str) or not v:
                continue
            paths = [v] if os.path.exists(v) else [v + ext for ext in PREFIX_FILES if os.path.exists(v + ext)]
            for p in paths:
                found[p] = _stat(p)
    return found


def _sidecar(name: str, path: str) -> str:
    return os.path.join(path, "provenance.json") if name == "output_dir" else f"{path}.provenance.json"


def record(command: click.Command) -> click.Command:
    """Make ``command`` write its provenance beside its main output."""
    callback = command.callback

    @functools.wraps(callback)
    def run(*args: Any, **params: Any) -> Any:
        result = callback(*args, **params)
        output = next(((n, params[n]) for n in OUTPUTS if params.get(n)), None)
        if output is not None:
            root = click.get_current_context().find_root()
            provenance = {
                "popmix_version": pm.__version__,
                "command": command.name,
                "argv": sys.argv[1:],
                "threads": (root.obj or {}).get("threads"),
                "parameters": {k: _jsonable(v) for k, v in params.items()},
                "inputs": _inputs(params),
                "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "python": platform.python_version(),
                "platform": platform.platform(),
            }
            with open(_sidecar(*output), "w") as f:
                json.dump(provenance, f, indent=2)
                f.write("\n")
        return result

    run.records_provenance = True  # type: ignore[attr-defined]
    command.callback = run
    return command

