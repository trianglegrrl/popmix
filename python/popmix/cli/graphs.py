"""Admixture graph commands: qpgraph and find-graphs."""

from __future__ import annotations

import json

import click

import popmix as pm

from ._common import (
    f2_source_options,
    load_f2_blocks,
    parse_population_list,
)
from ._errors import EXPECTED_ERRORS, fail


def read_graph_file(path: str) -> list:
    """Read an admixture graph: one edge per line, "from to" (whitespace-
    separated). Blank lines and text after "#" are ignored."""
    edges = []
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            if len(fields) != 2:
                raise click.BadParameter(f"{path}:{n}: expected 'from to', got {line.strip()!r}",
                                         param_hint="--graph")
            edges.append((fields[0], fields[1]))
    if not edges:
        raise click.BadParameter(f"{path} has no edges", param_hint="--graph")
    return edges


def graph_leaves(edges: list) -> list:
    """Leaves (nodes without outgoing edges) in order of first appearance."""
    order, sources = [], set()
    for a, b in edges:
        sources.add(a)
        for v in (a, b):
            if v not in order:
                order.append(v)
    return [v for v in order if v not in sources]


@click.command()
@click.option("--graph", "graph_file", required=True, type=click.Path(exists=True, dir_okay=False),
              help='Graph file: one edge per line, "from to"; the leaves are the populations')
@click.option("--admix-weights", default=None,
              help="Fit at these admixture weights (comma-separated, one per admixture node) "
                   "instead of optimizing them")
@click.option("--numstart", type=click.IntRange(min=1), default=10,
              help="Random starts for the admixture weights (default: 10)")
@click.option("--seed", type=int, default=1, help="Seed for the random starts (default: 1)")
@click.option("--diag", type=float, default=1e-4, help="Edge-length regularization (default: 1e-4)")
@click.option("--diag-f3", type=float, default=1e-5, help="f3 covariance regularization (default: 1e-5)")
@click.option("--unconstrained", is_flag=True, default=False,
              help="Allow negative drift lengths and weights outside [0, 1]")
@click.option("--plot", default=None, help="Also draw the fitted graph to this image file (needs popmix[plot])")
@click.option("--plot-fit", default=None,
              help="Also plot each f3's residual z to this image file (needs popmix[plot])")
@click.option("--dot", default=None, help="Also write the fitted graph in Graphviz DOT format to this file")
@click.option("--out-of-sample/--no-out-of-sample", default=True,
              help="Also score the graph out of sample: fitted on alternate jackknife blocks, scored on the "
                   "others (default: on; two more fits)")
@click.option("--compat", type=click.Choice(["admixtools2-2.0.10"]), default=None,
              help="Reproduce admixtools2: stop the fit where L-BFGS-B stalls before a zero-bound kink, and "
                   "with --dot write its DOT text exactly (see DEVIATIONS.md)")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json", "tsv"]), default="text",
              help="Output format")
def qpgraph(graph_file: str, admix_weights: str | None, numstart: int, seed: int, diag: float,
            diag_f3: float, unconstrained: bool, plot: str | None, plot_fit: str | None, dot: str | None,
            out_of_sample: bool, compat: str | None, input_prefix: str | None, f2_dir: str | None, block_size: float, maxmiss: float,
            all_chromosomes: bool, output_format: str):
    """Fit an admixture graph, as admixtools2 qpgraph does.

    Fits drift lengths and admixture weights to f3(base; i, j) for all pairs
    of the graph's leaves; the first leaf in the file is the base. A node
    with two incoming edges is an admixture node, and --admix-weights gives
    the weight of its incoming edge whose source node appears first.
    """
    edges = read_graph_file(graph_file)
    pops = graph_leaves(edges)
    weights = None
    if admix_weights is not None:
        try:
            weights = [float(w) for w in admix_weights.split(",")]
        except ValueError:
            raise click.BadParameter(f"not a list of numbers: {admix_weights!r}",
                                     param_hint="--admix-weights") from None

    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, pops, block_size, maxmiss, all_chromosomes)
        result = pm.qpgraph(f2_blocks, edges, admix_weights=weights, numstart=numstart, seed=seed,
                            diag=diag, diag_f3=diag_f3, constrained=not unconstrained, compat=compat)
        oos = (pm.qpgraph_out_of_sample(f2_blocks, edges, admix_weights=weights, numstart=numstart, seed=seed,
                                        diag=diag, diag_f3=diag_f3, constrained=not unconstrained, compat=compat)
               if out_of_sample else None)
        if dot:
            pm.write_dot(result, path=dot, compat=compat)
        if plot or plot_fit:
            from popmix import plot as plots

            if plot:
                plots.graph(result, path=plot)
            if plot_fit:
                plots.qpgraph_fit(result, path=plot_fit)

        if output_format == "json":
            output = {
                "pops": result.pops,
                "score": result.score,
                "score_out_of_sample": oos,
                "n_f3": result.n_f3,
                "n_params": result.n_params,
                "dof": result.dof,
                "penalized_score": result.penalized_score,
                "edges": [{"from": e.from_node, "to": e.to_node, "type": e.kind, "weight": e.weight,
                           "low": e.low, "high": e.high} for e in result.edges],
                "f3": [{"pop1": f.pop1, "pop2": f.pop2, "pop3": f.pop3, "est": f.est, "se": f.se,
                        "fit": f.fit, "diff": f.diff, "z": f.z, "p": f.p} for f in result.f3],
                "opt": [{"start": o.start, "weights": o.weights, "score": o.score,
                         "evaluations": o.evaluations, "convergence": o.convergence} for o in result.opt],
            }
            click.echo(json.dumps(output, indent=2))
        elif output_format == "tsv":
            click.echo("from\tto\ttype\tweight\tlow\thigh")
            for e in result.edges:
                low = "" if e.low is None else f"{e.low:.6g}"
                high = "" if e.high is None else f"{e.high:.6g}"
                click.echo(f"{e.from_node}\t{e.to_node}\t{e.kind}\t{e.weight:.6g}\t{low}\t{high}")
        else:
            worst = max(result.f3, key=lambda f: abs(f.z))
            click.echo(f"qpGraph: {len(pops)} populations, {len(edges)} edges")
            click.echo(f"f3 base population: {result.pops[0]}")
            click.echo(f"Score: {result.score:.4f}")
            click.echo(f"Parameters: {result.n_params} for {result.n_f3} f3 statistics ({result.dof} dof); "
                       f"penalized score {result.penalized_score:.4f}")
            if oos is not None:
                click.echo(f"Out-of-sample score (alternate blocks): {oos:.4f}")
            click.echo(f"Worst f3 residual: f3({worst.pop1}; {worst.pop2}, {worst.pop3}) "
                       f"z = {worst.z:.2f}")
            click.echo()
            click.echo(f"  {'from':<16} {'to':<16} {'type':<6} weight")
            for e in result.edges:
                click.echo(f"  {e.from_node:<16} {e.to_node:<16} {e.kind:<6} {e.weight:.6f}")

    except EXPECTED_ERRORS as e:
        fail(e)


@click.command("find-graphs")
@click.option("--populations", required=True, callback=parse_population_list,
              help="The graph's leaves (comma-separated); the first is the f3 base of every fit")
@click.option("--nadmix", type=click.IntRange(min=0), default=0, help="Number of admixture nodes (default: 0)")
@click.option("--outgroup", default=None, help="Only graphs whose root has this population as a direct child")
@click.option("--method", type=click.Choice(["stepwise", "exhaustive"]), default="stepwise",
              help="stepwise (default): add populations, then admixture edges, one at a time, keeping the "
                   "--beam-width best graphs; exhaustive: every graph, as admixtools2 generate_all_graphs()")
@click.option("--plot", default=None, help="Also draw the best graph to this image file (needs popmix[plot])")
@click.option("--beam-width", type=click.IntRange(min=1), default=20,
              help="Graphs kept after each stepwise step (default: 20)")
@click.option("--top", type=click.IntRange(min=1), default=5, help="Graphs to report (default: 5)")
@click.option("--output-graph", type=click.Path(dir_okay=False, writable=True), default=None,
              help="Write the best graph here, in the format 'popmix qpgraph --graph' reads")
@click.option("--numstart", type=click.IntRange(min=1), default=10,
              help="Random starts for each graph's admixture weights (default: 10)")
@click.option("--seed", type=int, default=1, help="Seed for the random starts (default: 1)")
@click.option("--out-of-sample/--no-out-of-sample", default=True,
              help="Also score the graphs reported out of sample, fitted on alternate jackknife blocks "
                   "(default: on; two more fits per graph)")
@f2_source_options
@click.option("--format", "output_format", type=click.Choice(["text", "json"]), default="text",
              help="Output format")
def find_graphs(populations: list, nadmix: int, outgroup: str | None, method: str,
                plot: str | None, beam_width: int, top: int, output_graph: str | None,
                numstart: int, seed: int, out_of_sample: bool, input_prefix: str | None, f2_dir: str | None,
                block_size: float, maxmiss: float, all_chromosomes: bool, output_format: str):
    """Search for the admixture graphs that best fit the data.

    Each candidate graph has POPULATIONS as leaves and NADMIX admixture nodes
    and is fitted as 'popmix qpgraph' fits it.
    """
    try:
        f2_blocks = load_f2_blocks(input_prefix, f2_dir, populations, block_size, maxmiss, all_chromosomes)
        graphs, n_evaluated = pm.find_graphs(f2_blocks, populations, nadmix=nadmix, outgroup=outgroup,
                                             method=method, beam_width=beam_width, numstart=numstart, seed=seed)
        if plot:
            from popmix import plot as plots

            plots.graph(graphs[0].edges, path=plot, title=f"score {graphs[0].score:.4f}")
        if output_graph is not None:
            with open(output_graph, "w", encoding="utf-8") as f:
                f.write(f"# popmix find-graphs: score {graphs[0].score!r}\n")
                for a, b in graphs[0].edges:
                    f.write(f"{a}\t{b}\n")

        shown = graphs[:top]
        n_params = pm.qpgraph(f2_blocks, shown[0].edges, numstart=1).n_params if shown else 0
        oos = [pm.qpgraph_out_of_sample(f2_blocks, g.edges, numstart=numstart, seed=seed) if out_of_sample else None
               for g in shown]
        if output_format == "json":
            output = {
                "n_evaluated": n_evaluated,
                "graphs": [{"score": g.score, "score_out_of_sample": o, "n_params": n_params, "edges": g.edges}
                           for g, o in zip(shown, oos)],
            }
            click.echo(json.dumps(output, indent=2))
        else:
            click.echo(f"find-graphs: {len(populations)} populations, {nadmix} admixture, {method} search")
            click.echo(f"Graphs fitted: {n_evaluated}, each with {n_params} parameters. The best score of many "
                       "is optimistic: compare graphs out of sample or with compare_fits.")
            for rank, (g, o) in enumerate(zip(shown, oos), 1):
                click.echo()
                tail = "" if o is None else f", out of sample {o:.4f}"
                click.echo(f"#{rank}  score {g.score:.4f}{tail}")
                for a, b in g.edges:
                    click.echo(f"  {a}\t{b}")

    except EXPECTED_ERRORS as e:
        fail(e)
