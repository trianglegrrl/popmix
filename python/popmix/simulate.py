"""msprime simulations of admixture graphs, as admixtools2's
``msprime_sim`` and ``msprime_genome``: each writes a Python script that
simulates genotypes under the graph with msprime and writes them in
EIGENSTRAT format, and can run it.

Nodes are dated as admixtools2 ``pseudo_dates`` does: one ``time`` step
per level of the graph drawn as a tree from the root (each node at its
shortest distance from the root), with the part of the graph below an
admixture node moved down until every admixture edge points down. Every
node has ``neff`` individuals, admixture proportions are ``admix_default``
for the parent that sorts first and 1 - ``admix_default`` for the other,
and the scripts are admixtools2's text.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from collections import deque
from collections.abc import Mapping, Sequence
from decimal import Decimal

Edge = tuple[str, str]


def _r(x: float | int) -> str:
    """A number as R's ``as.character`` writes it: up to 15 significant
    digits, in fixed or scientific notation, whichever is narrower."""
    x = float(x)
    if x == 0:
        return "0"
    if math.isinf(x) or math.isnan(x):
        return {math.inf: "Inf", -math.inf: "-Inf"}.get(x, "NaN")
    digits = f"{abs(x):.14e}"  # 15 significant digits
    mant, exp = digits.split("e")
    e = int(exp)
    sig = mant.replace(".", "").rstrip("0") or "0"
    sign = "-" if x < 0 else ""
    sci = sign + sig[0] + ("." + sig[1:] if len(sig) > 1 else "") + f"e{'-' if e < 0 else '+'}{abs(e):02d}"
    fixed = sign + format(Decimal(sig) * Decimal(10) ** (e - len(sig) + 1), "f")
    return fixed if len(fixed) <= len(sci) else sci


def _graph(edges: Sequence[Sequence[str]]):
    edges = [(a, b) for a, b, *_ in edges]
    nodes = list(dict.fromkeys(v for e in edges for v in e))
    children = {v: [] for v in nodes}
    parents = {v: [] for v in nodes}
    for a, b in edges:
        children[a].append(b)
        parents[b].append(a)
    return edges, nodes, children, parents


def pseudo_dates(edges: Sequence[Sequence[str]], time: float = 1000, fix_leaf: bool = False) -> dict[str, float]:
    """Each node's date in generations (admixtools2 ``pseudo_dates``)."""
    edges, nodes, children, parents = _graph(edges)
    root = next(v for v in nodes if not parents[v])
    depth, todo = {root: 0}, deque([root])
    while todo:
        v = todo.popleft()
        for c in children[v]:
            if c not in depth:
                depth[c] = depth[v] + 1
                todo.append(c)
    top = max(depth.values())
    y = {v: top - depth[v] for v in nodes}

    def below(v: str) -> set[str]:
        out, stack = {v}, [v]
        while stack:
            for c in children[stack.pop()]:
                if c not in out:
                    out.add(c)
                    stack.append(c)
        return out

    # fix_shiftdown: admixture edges in admixtools2 find_admixedges order
    for a in (v for v in nodes if len(parents[v]) == 2):
        for p in parents[a]:
            by = y[a] - y[p] + 1
            if by > 0:
                for v in below(a):
                    y[v] -= by
    low = min(y.values())
    out = {v: (y[v] - low) * time for v in nodes}
    if fix_leaf:
        for v in nodes:
            if not children[v]:
                out[v] = 0
    return out


def _model(edges, neff, ind_per_pop, time, fix_leaf, admix_default, ghost_lineages, shorten_admixed_leaves):
    edges, nodes, children, parents = _graph(edges)
    leaves = [v for v in nodes if not children[v]]
    admix = {v for v in nodes if len(parents[v]) > 1}
    if isinstance(time, Mapping):
        dates = dict(time)
    else:
        dates = pseudo_dates(edges, time, fix_leaf)
        if shorten_admixed_leaves:
            for a in admix:
                for c in children[a]:
                    if c in leaves:
                        dates[c] = 0.999 * dates[a]
    rows = sorted(((a, b, dates[a]) for a, b in edges), key=lambda r: r[2])
    weight = {}
    for a, b, _ in rows:
        if b in admix:
            weight[a, b] = admix_default if a == min(parents[b]) else 1 - admix_default
        else:
            weight[a, b] = 1.0
    if isinstance(neff, Mapping):
        if set(neff) < set(nodes):
            raise ValueError("'neff' has to be a single number or a named vector with a number for each node!")
        initsize = [neff[v] for v in nodes]
    else:
        popsize = {}
        for a, b, _ in rows:
            popsize.setdefault(b, 1.0 if b in admix else 1 / weight[a, b])
        for v in nodes:
            if not parents[v]:
                popsize[v] = 1.0
        initsize = [round(neff * popsize.get(v, neff * 100), 2) for v in nodes]
    if isinstance(ind_per_pop, Mapping):
        if sorted(ind_per_pop) != sorted(leaves):
            raise ValueError("'ind_per_pop' has to be a single number or a named vector with a number for each leaf node!")
        ipp = dict(ind_per_pop)
    else:
        ipp = dict.fromkeys(leaves, ind_per_pop)
    indnam = [f"{p}_{k}" for p, n in ipp.items() for k in range(1, n + 1)]

    lines = ['\nindnam = ["' + '", "'.join(indnam) + '"]']
    lines.append('\n\nsamples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( ('
                 + ", ".join(str(n) for n in ipp.values()) + '), ("' + '", "'.join(ipp) + '"), ('
                 + ", ".join(_r(dates[p]) for p in ipp) + ') )]')
    lines.append("\nnhap = int(len(indnam) * 2)\n")
    lines.append("\ndemo = msprime.Demography()\n")
    lines.append("".join(f'demo.add_population(initial_size = {_r(s)}, name = "{v}")\n' for v, s in zip(nodes, initsize)))  # noqa: B905
    # splits: by (ancestral, date), ancestral names in C order
    groups: dict[tuple[str, float], list[str]] = {}
    for a, b, d in rows:
        if b not in admix:
            groups.setdefault((a, d), []).append(b)
    for (a, d) in sorted(groups, key=lambda k: (k[0].encode(), k[1])):
        lines.append(f'demo.add_population_split(time={_r(d)}, ancestral="{a}", derived=["'
                     + '","'.join(groups[a, d]) + '"])\n')
    adm_rows = [(a, b, d) for a, b, d in rows if b in admix]
    if adm_rows:
        split_date = {a: d for a, b, d in rows if a in admix}
        first = {}
        for _, b, d in adm_rows:
            first[b] = min(first.get(b, d), d)
        agroups: dict[tuple[str, float], list[tuple[str, float]]] = {}
        for a, b, _ in adm_rows:
            d = split_date[b] if ghost_lineages else first[b]
            agroups.setdefault((b, d), []).append((a, weight[a, b]))
        for (b, d) in sorted(agroups, key=lambda k: (k[0].encode(), k[1])):
            src = agroups[b, d]
            lines.append(f'demo.add_admixture(time={_r(d)}, derived="{b}", ancestral=["'
                         + '","'.join(a for a, _ in src) + '"], proportions=[' + ",".join(_r(w) for _, w in src) + "])\n")
    lines.append("demo.sort_events()\n")
    return "".join(lines)


_IND_FILE = ("\n\nwith open('{o}.ind', 'w') as f:\n  for i in range(len(indnam)):\n    bytes = f.write(indnam[i] + "
             "'\\tU\\t' + indnam[i].rstrip(\"1234567890\").rstrip(\"_\") + '\\n')\n")


def _write(text: str, outpref: str, run: bool | str) -> str:
    path = f"{outpref}.py"
    with open(path, "w") as fh:
        fh.write(text)
    if run:
        subprocess.run([sys.executable if run is True else str(run), path], check=True)
    return os.path.abspath(path)


def msprime_sim(edges: Sequence[Sequence[str]], outpref: str | os.PathLike = "msprime_sim", nsnps: int = 1000,
                neff: float | Mapping[str, float] = 1000, ind_per_pop: int | Mapping[str, int] = 1,
                mutation_rate: float = 1e-3, time: float | Mapping[str, float] = 1000, fix_leaf: bool = False,
                admix_default: float = 0.5, run: bool | str = False, numcores: int | None = None,
                ghost_lineages: bool = False, shorten_admixed_leaves: bool = False) -> str:
    """Write a script that simulates ``nsnps`` unlinked SNPs under the graph
    (admixtools2 ``msprime_sim``), each from its own tree, and writes
    ``outpref``.geno/.snp/.ind; with ``run``, run it (with this Python, or
    the given interpreter).

    Parameters
    ----------
    edges
        The graph as (from, to) pairs.
    neff, ind_per_pop
        Population size of every node, and samples per leaf (or per-node /
        per-leaf mappings).
    time
        Generations per level of the graph, or each node's date.
    fix_leaf
        Sample every leaf at time 0.

    Returns
    -------
    str
        The script's path.
    """
    o = os.path.abspath(os.fspath(outpref))
    text = "import numpy\nimport msprime\nimport multiprocess\n"
    text += _model(edges, neff, ind_per_pop, time, fix_leaf, admix_default, ghost_lineages, shorten_admixed_leaves)
    text += f"\nnsnps = int({_r(nsnps)})"
    text += ("\n\ndef f(i):\n  ts = msprime.sim_ancestry(samples=samples, demography=demo)"
             f"\n  tree_sequence = msprime.sim_mutations(ts, rate={_r(mutation_rate)}, model='binary')"
             "\n  if tree_sequence.genotype_matrix().shape[0] > 0:"
             "\n    return (tree_sequence.genotype_matrix()[0,range(0,nhap,2)] + "
             "tree_sequence.genotype_matrix()[0,range(1,nhap,2)])"
             "\n  else:\n    return numpy.zeros(nhap//2, 'int')\n")
    text += f"\np = multiprocess.Pool(int(min(multiprocess.cpu_count(),{_r(100 if numcores is None else numcores)})))"
    text += "\ngt = numpy.array(p.map(f, range(nsnps)))"
    text += f"\n\nnumpy.savetxt('{o}.geno', gt, '%d', '')"
    text += (f"\n\nwith open('{o}.snp', 'w') as f:\n  for i in range(nsnps):\n    bytes = f.write('rs'+str(i+1)+"
             "'\\t1\\t' + str(i*50/float(nsnps-1)) + '\\t' + str(i*100) + '\\tA\\tC\\n')")
    text += _IND_FILE.format(o=o)
    return _write(text + "\n", o, run)


def msprime_genome(edges: Sequence[Sequence[str]], outpref: str | os.PathLike = "msprime_sim",
                   neff: float | Mapping[str, float] = 1000, ind_per_pop: int | Mapping[str, int] = 1,
                   mutation_rate: float = 1.25e-8, time: float | Mapping[str, float] = 1000, fix_leaf: bool = False,
                   nchr: int = 1, recomb_rate_chr: float = 2e-8, seq_length: float | Sequence[float] = 1e3,
                   admix_default: float = 0.5, run: bool | str = False, ghost_lineages: bool = False,
                   shorten_admixed_leaves: bool = False) -> str:
    """Write a script that simulates ``nchr`` recombining chromosomes of
    ``seq_length`` bases each under the graph (admixtools2
    ``msprime_genome``), unlinked from each other; otherwise as
    ``msprime_sim``."""
    o = os.path.abspath(os.fspath(outpref))
    text = "import numpy\nimport math\nimport msprime\nimport multiprocess\n"
    text += _model(edges, neff, ind_per_pop, time, fix_leaf, admix_default, ghost_lineages, shorten_admixed_leaves)
    text += f"\nr_chrom = {_r(recomb_rate_chr)}\nr_break = math.log(2)\n"
    if isinstance(seq_length, (int, float)):
        chr_pos = [i * seq_length for i in range(nchr + 1)]
    elif len(seq_length) == nchr:
        chr_pos = [0.0]
        for s in seq_length:
            chr_pos.append(chr_pos[-1] + s)
    else:
        raise ValueError("seq_length must be a number or one per chromosome")
    text += "chrom_positions = [" + ", ".join(_r(p) for p in chr_pos) + "]\n"
    map_pos = []
    for p in chr_pos:
        map_pos.append(p)
        if p != min(chr_pos) and p != max(chr_pos):
            map_pos.append(p + 1)
    text += "map_positions = [" + ", ".join(_r(p) for p in map_pos) + "]\n"
    rates = ["r_break" if i % 2 == 0 else "r_chrom" for i in range(1, len(map_pos))]
    text += "rates = [" + ", ".join(rates) + "]\n"
    text += "rate_map = msprime.RateMap(position=map_positions, rate=rates)\n"
    text += ("\nts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), "
             "msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)"
             f"\ntree_sequence = msprime.sim_mutations(ts, rate={_r(mutation_rate)}, model='binary')"
             "\nhap_gt = tree_sequence.genotype_matrix()"
             "\ngt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]"
             "\nnsnps = gt.shape[0]")
    text += ("\nts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, "
             "side='right') + 1")
    text += f"\n\nnumpy.savetxt('{o}.geno', gt, '%d', '')"
    text += (f"\n\nwith open('{o}.snp', 'w') as f:"
             "\n  for i, j in enumerate(tree_sequence.tables.sites.position):"
             "\n    chrm = int(ts_chroms[i])"
             "\n    if chrm == 1:"
             "\n      pos = int(j) "
             "\n    else:"
             "\n      pos = int(j - chrom_positions[chrm-1])"
             "\n    nsnps_chr = numpy.sum(ts_chroms == chrm)"
             "\n    pos_chr = numpy.sum(ts_chroms[:i] == chrm)"
             "\n    bytes = f.write('rs'+str(i+1)+'\\t' + str(chrm) + '\\t' + str(pos_chr*50/float(nsnps_chr-1)) + "
             "'\\t' + str(pos) + '\\tA\\tC\\n')")
    text += _IND_FILE.format(o=o)
    return _write(text + "\n", o, run)
