#!/usr/bin/env python3
"""Hashes of the individuals and populations AADR v66.p1 removed.

v66.p1 withdrew the 161 present-day genomes of Jacobs et al. 2019
(publication JacobsCell2019 in the v62 annotation), and with them 27
populations made only of those genomes. Their data must not be used, so
tests/python/test_data_ethics.py fails if any of their sample IDs or
population labels appears in the repository or the reference cache. The
names themselves are not stored, only SHA-256 hashes of each ID and label,
with and without its data-type suffix (".DG"), so this file does not list
what it guards against.

Usage: python tests/data_ethics/generate_removed.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

V62 = Path("/references/aadr/v62.0_1240k_public.anno")
V66P1 = Path("/references/aadr/v66.p1/v66.p1_1240K.aadr.patch.PUB.ind")
OUT = Path(__file__).parent / "removed_v66p1.json"


def variants(name: str) -> set[str]:
    base = name.rsplit(".", 1)[0] if name.rsplit(".", 1)[-1] in ("DG", "SDG", "AG", "HO") else name
    return {name, base}


def main() -> None:
    header, *rows = [line.split("\t") for line in V62.read_text(encoding="utf-8").splitlines()]
    pub = next(i for i, h in enumerate(header) if h.startswith("Publication abbreviation"))
    group = next(i for i, h in enumerate(header) if h.startswith("Group ID"))
    kept = {line.split()[0] for line in V66P1.read_text(encoding="utf-8").splitlines()}
    jacobs = [r for r in rows if r[pub] == "JacobsCell2019"]
    removed = [r for r in jacobs if r[0] not in kept]
    assert len(removed) == len(jacobs) == 161, (len(removed), len(jacobs))
    pop_sizes: dict[str, int] = {}
    for r in rows:
        pop_sizes[r[group]] = pop_sizes.get(r[group], 0) + 1
    removed_pops = {r[group] for r in removed}
    whole = sorted(p for p in removed_pops if sum(r[group] == p for r in removed) == pop_sizes[p])
    names = set()
    for r in removed:
        names |= variants(r[0])
    for p in whole:
        names |= variants(p)
    OUT.write_text(json.dumps({
        "source": "AADR v62.0 annotation, JacobsCell2019, absent from v66.p1",
        "individuals": len(removed),
        "populations": len(whole),
        "sha256": sorted(hashlib.sha256(n.encode()).hexdigest() for n in names),
    }, indent=1) + "\n", encoding="utf-8")
    print(f"{len(removed)} individuals, {len(whole)} populations, {len(names)} hashes -> {OUT}")


if __name__ == "__main__":
    main()
